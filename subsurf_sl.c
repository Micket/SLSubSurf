/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2005 Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/subsurf_sl.c
 *  \ingroup bke
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "MEM_guardedalloc.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BLI_utildefines.h"
#include "BLI_bitmap.h"
#include "BLI_blenlib.h"
#include "BLI_edgehash.h"
#include "BLI_math.h"
#include "BLI_memarena.h"
#include "BLI_pbvh.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_global.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_multires.h"
#include "BKE_paint.h"
#include "BKE_scene.h"
#include "BKE_subsurf.h"
#include "BKE_tessmesh.h"

#include "PIL_time.h"
#include "BLI_array.h"

#include "GL/glew.h"

#include "GPU_draw.h"
#include "GPU_extensions.h"
#include "GPU_material.h"

#include "SLSubSurf.h"

extern GLubyte stipple_quarttone[128]; /* glutil.c, bad level data */

typedef struct SLDerivedMesh {
	DerivedMesh dm;
	SLSubSurf *ss;
	struct PBVH *pbvh;
	DMFlagMat *faceFlags;
	int freeSS;
} SLDerivedMesh;

/*
 * This code copies over the data from the previous derived mesh on top of the existing (from scratch, nothing incremental)
 * Vertex coordinates can come from an external array "vertexCos" (used in for example sculpting in multisurf).
 */
static void slss_sync_from_derivedmesh(SLSubSurf *ss, DerivedMesh *dm, float (*vertexCos)[3])
{
	MVert *mv;
	MEdge *me;
	MLoop *ml, *mloop;
	MPoly *mp;
	int totvert = dm->getNumVerts(dm);
	int totedge = dm->getNumEdges(dm);
	int i, j;
	void **fVerts = NULL;
	BLI_array_declare(fVerts);

	mv = dm->getVertArray(dm);
	//index = (int *)dm->getVertDataArray(dm, CD_ORIGINDEX);
	for (i = 0; i < totvert; i++, mv++) {
		SL_addVert(ss, SET_INT_IN_POINTER(i), vertexCos ? vertexCos[i] : mv->co, 0);
		//((int *)ccgSubSurf_getVertUserData(ss, v))[1] = (index) ? *index++ : i;
	}

	me = dm->getEdgeArray(dm);
	//index = (int *)dm->getEdgeDataArray(dm, CD_ORIGINDEX);
	for (i = 0; i < totedge; i++, me++) {
		SL_addEdge(ss, SET_INT_IN_POINTER(i), SET_INT_IN_POINTER(me->v1), SET_INT_IN_POINTER(me->v2), me->crease);
		//((int *)ccgSubSurf_getEdgeUserData(ss, e))[1] = (index) ? *index++ : i;
	}

	mloop = dm->getLoopArray(dm);
	mp = dm->getPolyArray(dm);
	//index = DM_get_poly_data_layer(dm, CD_ORIGINDEX);
	for (i = 0; i < dm->numPolyData; i++, mp++) {
		BLI_array_empty(fVerts);
		BLI_array_grow_items(fVerts, mp->totloop);
		ml = mloop + mp->loopstart;
		for (j = 0; j < mp->totloop; j++, ml++) {
			fVerts[j] = SET_INT_IN_POINTER(ml->v);
		}
		SL_addFace(ss, SET_INT_IN_POINTER(i), mp->totloop, fVerts);
		//((int *)ccgSubSurf_getFaceUserData(ss, f))[1] = (index) ? *index++ : i;
	}
	BLI_array_free(fVerts);
	SL_processSync(ss);
	SL_renumberAll(ss);
	printf("Sync processed. Now time to copy it back...\n");
}

/*
 * This function only updates the existing vertex coordinates.
 */
static void slss_update_coords_from_derivedmesh(SLSubSurf *ss, DerivedMesh *dm, float (*vertexCos)[3])
{
	MVert *mv;
	MEdge *me;
	int totvert = dm->getNumVerts(dm);
	int totedge = dm->getNumEdges(dm);
	int i;

	mv = dm->getVertArray(dm);
	//index = (int *)dm->getVertDataArray(dm, CD_ORIGINDEX);
	for (i = 0; i < totvert; i++, mv++) {
		SL_updateVert(ss, SET_INT_IN_POINTER(i), vertexCos ? vertexCos[i] : mv->co, 0);
	}

	me = dm->getEdgeArray(dm);
	//index = (int *)dm->getEdgeDataArray(dm, CD_ORIGINDEX);
	for (i = 0; i < totedge; i++, me++) {
		SL_updateEdge(ss, SET_INT_IN_POINTER(i), me->crease);
	}

	SL_processSync(ss);
	printf("Sync processed. Now time to copy it back...\n");
}

/*
 * This function creates the underlying subsurf structure.
 */
static SLSubSurf *_getSLSubSurf(SLSubSurf *prevSS, int smoothing)
{
	if (prevSS) {
		prevSS->smoothing = smoothing;
		return prevSS;
	}
	return SL_SubSurf_new(smoothing);
}


static SLDerivedMesh *getSLDerivedMesh(SLSubSurf *ss,
									int drawInteriorEdges,
									int useSubsurfUv,
									DerivedMesh *dm);


/*
 * Entrypoint for modifier. Should return the output DM based on the input DM.
 */
struct DerivedMesh *sl_subsurf_make_derived_from_derived(
	struct DerivedMesh *dm,
	struct SubsurfModifierData *smd,
	float (*vertCos)[3],
	SubsurfFlags flags)
{
	SLSubSurf *ss;
	int useSimple = smd->subdivType == ME_SIMPLE_SL_SUBSURF;
	//SLFlags useAging = smd->flags & eSubsurfModifierFlag_DebugIncr ? CCG_USE_AGING : 0;
	int useSubsurfUv = smd->flags & eSubsurfModifierFlag_SubsurfUv;
	int drawInteriorEdges = !(smd->flags & eSubsurfModifierFlag_ControlEdges);
	SLDerivedMesh *result;
	if (flags & SUBSURF_FOR_EDIT_MODE) {
		int levels = (smd->modifier.scene) ? get_render_subsurf_level(&smd->modifier.scene->r, smd->levels) : smd->levels;

		ss = _getSLSubSurf(NULL, !useSimple);
		smd->emCache = NULL; // TODO: Deal with caches
		slss_sync_from_derivedmesh(ss, dm, vertCos);

		result = getSLDerivedMesh(ss, drawInteriorEdges, useSubsurfUv, dm);
		result->freeSS = 1; // TODO: Revert after mCache support is finalized
	}
	else if (flags & SUBSURF_USE_RENDER_PARAMS) {
		/* Do not use cache in render mode. */
		int levels = (smd->modifier.scene) ? get_render_subsurf_level(&smd->modifier.scene->r, smd->renderLevels) : smd->renderLevels;

		if (levels == 0)
			return dm;

		ss = _getSLSubSurf(NULL, !useSimple);

		slss_sync_from_derivedmesh(ss, dm, vertCos);

		result = getSLDerivedMesh(ss, drawInteriorEdges, useSubsurfUv, dm);
		result->freeSS = 1;
	}
	else {
		int useIncremental = (smd->flags & eSubsurfModifierFlag_Incremental);
		int levels = (smd->modifier.scene) ? get_render_subsurf_level(&smd->modifier.scene->r, smd->levels) : smd->levels;

		/* It is quite possible there is a much better place to do this. It
		 * depends a bit on how rigorously we expect this function to never
		 * be called in editmode. In semi-theory we could share a single
		 * cache, but the handles used inside and outside editmode are not
		 * the same so we would need some way of converting them. Its probably
		 * not worth the effort. But then why am I even writing this long
		 * comment that no one will read? Hmmm. - zr
		 *
		 * Addendum: we can't really ensure that this is never called in edit
		 * mode, so now we have a parameter to verify it. - brecht
		 */
		if (!(flags & SUBSURF_IN_EDIT_MODE) && smd->emCache) {
			smd->cacheFree(smd->emCache);
			smd->emCache = NULL;
		}

		if (useIncremental && (flags & SUBSURF_IS_FINAL_CALC)) {
			ss = _getSLSubSurf(NULL, !useSimple);
			//smd->mCache = ss = _getSLSubSurf(smd->mCache, !useSimple);
			//smd->cacheFree = SL_SubSurf_free;

			slss_sync_from_derivedmesh(ss, dm, vertCos);

			result = getSLDerivedMesh(ss, drawInteriorEdges, useSubsurfUv, dm);
			result->freeSS = 1; // TODO: Revert after mCache support is finalized
		}
		else {
			if (smd->mCache && (flags & SUBSURF_IS_FINAL_CALC)) {
				smd->cacheFree(smd->mCache);
				smd->mCache = NULL;
			}

			ss = _getSLSubSurf(NULL, !useSimple);
			slss_sync_from_derivedmesh(ss, dm, vertCos);

			result = getSLDerivedMesh(ss, drawInteriorEdges, useSubsurfUv, dm);

			if (flags & SUBSURF_IS_FINAL_CALC)
				;//smd->mCache = ss;
				//smd->cacheFree = SL_SubSurf_free;
			else
				result->freeSS = 1;
		}
	}
	/*
	{
		// Checking output;
		DerivedMesh *newdm = (DerivedMesh*)result;
		int i;
		MVert *mverts;
		printf("output verts = %d\n", newdm->getNumVerts(newdm));
		mverts = MEM_callocN(sizeof(MVert) *newdm->getNumVerts(newdm), "test verts");
		newdm->copyVertArray(newdm, mverts);
		for (i = 0; i < newdm->getNumVerts(newdm); i++) {
			printf("Vert %d; coords = [%e, %e, %e]\n", i, mverts[i].co[0],mverts[i].co[1],mverts[i].co[2]);
		}
		MEM_freeN(mverts);
	}
	*/
	return (DerivedMesh *)result;
}

// Functions for the DerivedMesh-class;
static void slDM_getMinMax(DerivedMesh *dm, float min_r[3], float max_r[3]) {
	SLDerivedMesh *sldm = (SLDerivedMesh*)dm;
	SL_getMinMax(sldm->ss, min_r, max_r);
}

static void slDM_recalcTessellation(DerivedMesh *UNUSED(dm)) { /* Nothing to do */ }
static void slDM_calcNormals(DerivedMesh *UNUSED(dm)) { /* Nothing to do */ }

// TODO: Can i do this or do i have to story this separately?
static int slDM_getNumVerts(DerivedMesh *dm) { return dm->numVertData; }
static int slDM_getNumEdges(DerivedMesh *dm) { return dm->numEdgeData; }
static int slDM_getNumLoops(DerivedMesh *dm) { return dm->numLoopData; }
static int slDM_getNumPolys(DerivedMesh *dm) { return dm->numPolyData; }
static int slDM_getNumTessFaces(DerivedMesh *dm) { return dm->numTessFaceData; }

// Copy stuff
static void slDM_copyFinalVertArray(DerivedMesh *dm, MVert *mvert) {
	SL_copyNewVerts(((SLDerivedMesh *)dm)->ss, mvert);
}
static void slDM_copyFinalEdgeArray(DerivedMesh *dm, MEdge *medge) {
	SLDerivedMesh *sldm = (SLDerivedMesh*)dm;
	SL_copyNewEdges(sldm->ss, medge);
}
static void slDM_copyFinalLoopArray(DerivedMesh *dm, MLoop *mloop) {
	SLDerivedMesh *sldm = (SLDerivedMesh*)dm;
	SL_copyNewLoops(sldm->ss, mloop);
}
static void slDM_copyFinalPolyArray(DerivedMesh *dm, MPoly *mpoly) {
	SLDerivedMesh *sldm = (SLDerivedMesh*)dm;
	DMFlagMat *faceFlags = sldm->faceFlags;
	SL_copyNewPolys(sldm->ss, faceFlags, mpoly);
}
static void slDM_copyFinalTessFaceArray(DerivedMesh *dm, MFace *mface) {
	SLDerivedMesh *sldm = (SLDerivedMesh*)dm;
	DMFlagMat *faceFlags = sldm->faceFlags;
	SL_copyNewTessFaces(sldm->ss, faceFlags, mface);
}

static void slDM_release(DerivedMesh *dm) {
	SLDerivedMesh *sldm = (SLDerivedMesh*)dm;
	if (DM_release(dm)) { // Just doing what CCG code does...
		if (sldm->freeSS) SL_SubSurf_free(sldm->ss);
		MEM_freeN(sldm->faceFlags);
		MEM_freeN(sldm);
	}
}

// Unsure about these three
static void *slDM_get_vert_data_layer(DerivedMesh *dm, int type) {
	if (type == CD_ORIGINDEX) {
		/* create origindex on demand to save memory */
		SLDerivedMesh *ssdm = (SLDerivedMesh *)dm;
		SLSubSurf *ss = ssdm->ss;
		int *origindex;
		int a, tot;

		/* Avoid re-creation if the layer exists already */
		origindex = DM_get_vert_data_layer(dm, CD_ORIGINDEX);
		if (origindex) {
			return origindex;
		}

		DM_add_vert_layer(dm, CD_ORIGINDEX, CD_CALLOC, NULL);
		origindex = DM_get_vert_data_layer(dm, CD_ORIGINDEX);

		// original vertices are at the beginning
		a = 0;
		for (BLI_ghashIterator_init(ss->it,ss->verts); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it), a++) {
			origindex[a] = GET_INT_FROM_POINTER(BLI_ghashIterator_getKey(ss->it));
		}
		// Then new verts, no original index
		tot = SL_giveTotalNumberOfSubVerts(ss);
		for (; a < tot; a++) {
			origindex[a] = ORIGINDEX_NONE;
		}

		return origindex;
	}

	return DM_get_vert_data_layer(dm, type);
}
static void *slDM_get_edge_data_layer(DerivedMesh *dm, int type) {
	if (type == CD_ORIGINDEX) {
		/* create origindex on demand to save memory */
		SLDerivedMesh *ssdm = (SLDerivedMesh *)dm;
		SLSubSurf *ss = ssdm->ss;
		int *origindex;
		int a, tot;

		/* Avoid re-creation if the layer exists already */
		origindex = DM_get_edge_data_layer(dm, CD_ORIGINDEX);
		if (origindex) {
			return origindex;
		}

		DM_add_edge_layer(dm, CD_ORIGINDEX, CD_CALLOC, NULL);
		origindex = DM_get_edge_data_layer(dm, CD_ORIGINDEX);

		a = 0;
		for (BLI_ghashIterator_init(ss->it, ss->edges); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it), a++) {
			int mapIndex = GET_INT_FROM_POINTER(BLI_ghashIterator_getKey(ss->it));
			origindex[a++] = mapIndex;
			origindex[a++] = mapIndex;
		}
		// Then new internal edges, no original index
		tot = SL_giveTotalNumberOfSubEdges(ss);
		for (; a < tot; a++) {
			origindex[a] = ORIGINDEX_NONE;
		}

		return origindex;
	}

	return DM_get_edge_data_layer(dm, type);
}
static void *slDM_get_tessface_data_layer(DerivedMesh *dm, int type) {
	if (type == CD_ORIGINDEX) {
		/* create origindex on demand to save memory */
		SLDerivedMesh *sldm = (CCGDerivedMesh *)dm;
		SLSubSurf *ss = sldm->ss;
		int *origindex;
		int a, i;

		/* Avoid re-creation if the layer exists already */
		origindex = DM_get_tessface_data_layer(dm, CD_ORIGINDEX);
		if (origindex) {
			return origindex;
		}

		DM_add_tessface_layer(dm, CD_ORIGINDEX, CD_CALLOC, NULL);
		origindex = DM_get_tessface_data_layer(dm, CD_ORIGINDEX);

		a = 0;
		for (BLI_ghashIterator_init(ss->it,ss->faces); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it), a++) {
			SLFace *f = BLI_ghashIterator_getValue(ss->it);
			int mapIndex = GET_INT_FROM_POINTER(BLI_ghashIterator_getKey(ss->it));
			for (i = 0; i < SL_giveNumberOfInternalFaces(f); i++, a++)
				origindex[a] = mapIndex;
		}

		return origindex;
	}

	return DM_get_tessface_data_layer(dm, type);
}

// I fail to see the point of these;
static void *slDM_get_vert_data(DerivedMesh *dm, int index, int type) {
	if (type == CD_ORIGINDEX)
		slDM_get_vert_data_layer(dm, type);
	return DM_get_vert_data(dm, index, type);
}
static void *slDM_get_edge_data(DerivedMesh *dm, int index, int type) {
	if (type == CD_ORIGINDEX)
		slDM_get_edge_data_layer(dm, type);
	return DM_get_edge_data(dm, index, type);
}

static void *slDM_get_tessface_data(DerivedMesh *dm, int index, int type) {
	if (type == CD_ORIGINDEX)
		slDM_get_tessface_data_layer(dm, type);
	return DM_get_tessface_data(dm, index, type);
}

static void slDM_getVertCos(DerivedMesh *dm, float (*cos)[3])
{
	// TODO: This is basically just like asking for the vertices. Why duplicate this functionality? Its not like it needs the performance.
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;
	int i, index;
	
	i = 0;
	for (index = 0; index < BLI_ghash_size(ss->verts); index++) {
		SLVert *v = BLI_ghash_lookup(ss->verts, SET_INT_IN_POINTER(index)); // Not sure if i could just use the ghash iterator instead...
		copy_v3_v3(cos[i++], v->sl_coords);
	}
	for (index = 0; index < BLI_ghash_size(ss->edges); index++) {
		SLEdge *e = BLI_ghash_lookup(ss->edges, SET_INT_IN_POINTER(index));
		copy_v3_v3(cos[i++], e->sl_coords);
	}
	for (index = 0; index < BLI_ghash_size(ss->faces); index++) {
		SLFace *f = BLI_ghash_lookup(ss->faces, SET_INT_IN_POINTER(index));
		if (f->numVerts > 3) {
			copy_v3_v3(cos[i++], f->centroid);
		}
	}
}

static struct PBVH *slDM_getPBVH(Object *ob, DerivedMesh *dm)
{
	// TODO: I have no idea about any of this code. The grids and all that crap just seems relevant to CCG.
	SLDerivedMesh *ssdm = (SLDerivedMesh *)dm;
	SLSubSurf *ss = ssdm->ss;
	
	if (!ob) {
		ssdm->pbvh = NULL;
		return NULL;
	}
	
	if (!ob->sculpt)
		return NULL;
	
	if (ob->sculpt->pbvh) {
		ssdm->pbvh = ob->sculpt->pbvh;
	}
	
	if (ssdm->pbvh)
		return ssdm->pbvh;

	if (ob->type == OB_MESH) {
		Mesh *me = ob->data;
		ob->sculpt->pbvh = ssdm->pbvh = BLI_pbvh_new();
		BLI_assert(!(me->mface == NULL && me->mpoly != NULL)); /* BMESH ONLY complain if mpoly is valid but not mface */
		BLI_pbvh_build_mesh(ssdm->pbvh, me->mface, me->mvert,
							me->totface, me->totvert, &me->vdata);
	}
	
	return ssdm->pbvh;
}

// OpenGL drawing stuff (self explanatory)
static void slDM_drawVerts(DerivedMesh *dm) {
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;

	glBegin(GL_POINTS);
	for (BLI_ghashIterator_init(ss->it, ss->verts); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it)) {
		SLVert *v = BLI_ghashIterator_getValue(ss->it);
		glVertex3fv(v->sl_coords);
	}
	for (BLI_ghashIterator_init(ss->it, ss->edges); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it)) {
		SLEdge *e = BLI_ghashIterator_getValue(ss->it);
		glVertex3fv(e->sl_coords);
	}
	for (BLI_ghashIterator_init(ss->it, ss->faces); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it)) {
		SLFace *f = BLI_ghashIterator_getValue(ss->it);
		if (f->numVerts > 3) {
			glVertex3fv(f->centroid);
		}
	}
	glEnd();
}

static void slDM_drawEdges(DerivedMesh *dm, int drawLooseEdges, int drawAllEdges) {
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;
	int i;

	for (BLI_ghashIterator_init(ss->it, ss->edges); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it)) {
		SLEdge *e = BLI_ghashIterator_getValue(ss->it);

		if (!drawLooseEdges && e->numFaces == 0)
			continue;

		// Not sure what this means.
		//if (!drawAllEdges && ssdm->edgeFlags && !(ssdm->edgeFlags[j] & ME_EDGEDRAW))
		//	continue;

		// Or what aging is
		/*if (useAging && !(G.f & G_BACKBUFSEL)) {
			int ageCol = 255 - ccgSubSurf_getEdgeAge(ss, e) * 4;
			glColor3ub(0, ageCol > 0 ? ageCol : 0, 0);
		}*/

		glBegin(GL_LINE_STRIP);
		glVertex3fv(e->v0->sl_coords);
		glVertex3fv(e->sl_coords);
		glVertex3fv(e->v1->sl_coords);
		glEnd();
	}

	/*if (useAging && !(G.f & G_BACKBUFSEL)) {
		glColor3ub(0, 0, 0);
	}*/

	if (/*ssdm->drawInteriorEdges*/ 1) { // TODO: Why wouldn't we?

		for (BLI_ghashIterator_init(ss->it, ss->faces); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it)) {
			SLFace *f = BLI_ghashIterator_getValue(ss->it);

			if (f->numVerts == 3) { // Triangles are split differently from quads and ngons
				glBegin(GL_LINE_LOOP);
				glVertex3fv(f->edges[0]->sl_coords);
				glVertex3fv(f->edges[1]->sl_coords);
				glVertex3fv(f->edges[2]->sl_coords);
				glEnd();
			}
			else {
				glBegin(GL_LINES);
				for (i = 0; i < f->numVerts; i++) {
					glVertex3fv(f->centroid);
					glVertex3fv(f->edges[i]->sl_coords);
				}
				glEnd();
			}
		}
	}
}

static void slDM_drawLooseEdges(DerivedMesh *dm) {
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;

	for (BLI_ghashIterator_init(ss->it, ss->edges); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it)) {
		SLEdge *e = BLI_ghashIterator_getValue(ss->it);
		if (e->numFaces == 0) {
			glBegin(GL_LINE_STRIP);
			glVertex3fv(e->v0->sl_coords);
			glVertex3fv(e->sl_coords);
			glVertex3fv(e->v1->sl_coords);
			glEnd();
		}
	}
}

static void ssDM_glNormalFastTri(float *a, float *b, float *c)
{
	float a_bX = b[0] - a[0], a_bY = b[1] - a[1], a_bZ = b[2] - a[2];
	float a_cX = c[0] - a[0], a_cY = c[1] - a[1], a_cZ = c[2] - a[2];
	float no[3];

	no[0] = a_cY * a_bZ - a_cZ * a_bY;
	no[1] = a_cZ * a_bX - a_cX * a_bZ;
	no[2] = a_cX * a_bY - a_cY * a_bX;
	/* don't normalize, GL_NORMALIZE is enabled */
	// TODO: is it? It looks wrong.
	glNormal3fv(no);
}
static void ssDM_glNormalFast(float *a, float *b, float *c, float *d)
{
	float a_cX = c[0] - a[0], a_cY = c[1] - a[1], a_cZ = c[2] - a[2];
	float b_dX = d[0] - b[0], b_dY = d[1] - b[1], b_dZ = d[2] - b[2];
	float no[3];

	no[0] = b_dY * a_cZ - b_dZ * a_cY;
	no[1] = b_dZ * a_cX - b_dX * a_cZ;
	no[2] = b_dX * a_cY - b_dY * a_cX;
	/* don't normalize, GL_NORMALIZE is enabled */
	glNormal3fv(no);
}
static void slDM_drawFacesSolid(DerivedMesh *dm, float (*partial_redraw_planes)[4], int fast, DMSetMaterial setMaterial) {
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;
	DMFlagMat *faceFlags = ssdm->faceFlags;
	int i;
	int drawcurrent = 0, matnr = -1, shademodel = -1;

	//CCG_key_top_level(&key, ss);
	//ssdm_pbvh_update(ssdm);

	for (BLI_ghashIterator_init(ss->it, ss->faces); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it)) {
		SLFace *f = BLI_ghashIterator_getValue(ss->it);

		int index = GET_INT_FROM_POINTER(BLI_ghashIterator_getKey(ss->it));
		int new_matnr, new_shademodel;

		if (faceFlags) {
			new_shademodel = (faceFlags[index].flag & ME_SMOOTH) ? GL_SMOOTH : GL_FLAT;
			new_matnr = faceFlags[index].mat_nr;
		}
		else {
			new_shademodel = GL_SMOOTH;
			new_matnr = 0;
		}

		if (shademodel != new_shademodel || matnr != new_matnr) {
			matnr = new_matnr;
			shademodel = new_shademodel;

			drawcurrent = setMaterial(matnr + 1, NULL);

			glShadeModel(shademodel);
		}

		if (!drawcurrent)
			continue;

		// TODO: Smooth version
		if (f->numVerts == 3) {
			// 1x3 strip + 1 single triangle, not worth it (probably?)
			glBegin(GL_TRIANGLES);

			ssDM_glNormalFastTri(f->verts[0]->sl_coords, f->edges[0]->sl_coords, f->edges[2]->sl_coords);
			glVertex3fv(f->verts[0]->sl_coords);
			glVertex3fv(f->edges[0]->sl_coords);
			glVertex3fv(f->edges[2]->sl_coords);

			ssDM_glNormalFastTri(f->verts[1]->sl_coords, f->edges[1]->sl_coords, f->edges[0]->sl_coords);
			glVertex3fv(f->verts[1]->sl_coords);
			glVertex3fv(f->edges[1]->sl_coords);
			glVertex3fv(f->edges[0]->sl_coords);

			ssDM_glNormalFastTri(f->verts[2]->sl_coords, f->edges[2]->sl_coords, f->edges[1]->sl_coords);
			glVertex3fv(f->verts[2]->sl_coords);
			glVertex3fv(f->edges[2]->sl_coords);
			glVertex3fv(f->edges[1]->sl_coords);

			ssDM_glNormalFastTri(f->verts[0]->sl_coords, f->edges[1]->sl_coords, f->edges[2]->sl_coords);
			glVertex3fv(f->edges[0]->sl_coords);
			glVertex3fv(f->edges[1]->sl_coords);
			glVertex3fv(f->edges[2]->sl_coords);
			glEnd();
		} else {
			glBegin(GL_QUADS);
			for (i = 0; i < f->numVerts; i++) {
				ssDM_glNormalFast(f->verts[i]->sl_coords, f->edges[i]->sl_coords, f->centroid, f->edges[(i - 1 + f->numVerts) % f->numVerts]->sl_coords);
				glVertex3fv(f->verts[i]->sl_coords);
				glVertex3fv(f->edges[i]->sl_coords);
				glVertex3fv(f->centroid);
				glVertex3fv(f->edges[(i + f->numVerts - 1) % f->numVerts]->sl_coords);
			}
			glEnd();
		}
	}
}

static void slDM_drawFacesTex(DerivedMesh *dm,
							  DMSetDrawOptionsTex setDrawOptions,
							  DMCompareDrawOptions compareDrawOptions,
							  void *userData) {
	// TODO
	printf("slDM_drawFacesTex\n");
}

static void slDM_drawFacesGLSL(DerivedMesh *dm, DMSetMaterial setMaterial) {
	// TODO
	printf("slDM_drawFacesGLSL\n");
}

static void slDM_drawMappedFacesGLSL(DerivedMesh *dm,
									 DMSetMaterial setMaterial,
									 DMSetDrawOptions setDrawOptions,
									 void *userData) {
	// TODO
	printf("slDM_drawMappedFacesGLSL\n");
}

static void slDM_drawMappedFaces(DerivedMesh *dm,
								 DMSetDrawOptions setDrawOptions,
								 DMSetMaterial setMaterial,
								 DMCompareDrawOptions compareDrawOptions,
								 void *userData, DMDrawFlag flag) {
	// TODO
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;
	MCol *mcol = NULL;
	DMFlagMat *faceFlags = ssdm->faceFlags;
	int useColors = flag & DM_DRAW_USE_COLORS;
	int i, drawSmooth;
	
	//CCG_key_top_level(&key, ss);
	
	/* currently unused -- each original face is handled separately */
	(void)compareDrawOptions;
	
	if (useColors) {
		mcol = dm->getTessFaceDataArray(dm, CD_PREVIEW_MCOL);
		if (!mcol)
			mcol = dm->getTessFaceDataArray(dm, CD_MCOL);
	}
	
	
	for (BLI_ghashIterator_init(ss->it, ss->faces); !BLI_ghashIterator_isDone(ss->it); BLI_ghashIterator_step(ss->it)) {
		SLFace *f = BLI_ghashIterator_getValue(ss->it);
		int origIndex = GET_INT_FROM_POINTER(BLI_ghashIterator_getKey(ss->it));
		unsigned char *cp = NULL;

		if (flag & DM_DRAW_ALWAYS_SMOOTH) drawSmooth = 1;
		else if (faceFlags) drawSmooth = (faceFlags[origIndex].flag & ME_SMOOTH);
		else drawSmooth = 1;
		
		if (mcol) {
			cp = (unsigned char *)mcol;
			mcol += f->numVerts * 4; // 4 is used for what?
		}
		
		{
			DMDrawOption draw_option = DM_DRAW_OPTION_NORMAL;
			
			// I don't even know what this mysterious "index" is supposed to be. There is only ONE index.
			//if (index == ORIGINDEX_NONE)
				draw_option = setMaterial(faceFlags ? faceFlags[origIndex].mat_nr + 1 : 1, NULL);  /* XXX, no faceFlags no material */
			//else if (setDrawOptions)
			//	draw_option = setDrawOptions(userData, index);
				
			if (draw_option != DM_DRAW_OPTION_SKIP) {
				if (draw_option == DM_DRAW_OPTION_STIPPLE) {
					glEnable(GL_POLYGON_STIPPLE);
					glPolygonStipple(stipple_quarttone);
				}
					
				/* no need to set shading mode to flat because
				 *  normals are already used to change shading */
				glShadeModel(GL_SMOOTH);
				if (f->numVerts == 3) {
					glBegin(GL_TRIANGLES);
					ssDM_glNormalFastTri(f->verts[0]->sl_coords, f->edges[0]->sl_coords, f->edges[2]->sl_coords);
					glVertex3fv(f->verts[0]->sl_coords);
					glVertex3fv(f->edges[0]->sl_coords);
					glVertex3fv(f->edges[2]->sl_coords);
					
					ssDM_glNormalFastTri(f->verts[1]->sl_coords, f->edges[1]->sl_coords, f->edges[0]->sl_coords);
					glVertex3fv(f->verts[1]->sl_coords);
					glVertex3fv(f->edges[1]->sl_coords);
					glVertex3fv(f->edges[0]->sl_coords);
					
					ssDM_glNormalFastTri(f->verts[2]->sl_coords, f->edges[2]->sl_coords, f->edges[1]->sl_coords);
					glVertex3fv(f->verts[2]->sl_coords);
					glVertex3fv(f->edges[2]->sl_coords);
					glVertex3fv(f->edges[1]->sl_coords);
					
					ssDM_glNormalFastTri(f->verts[0]->sl_coords, f->edges[1]->sl_coords, f->edges[2]->sl_coords);
					glVertex3fv(f->edges[0]->sl_coords);
					glVertex3fv(f->edges[1]->sl_coords);
					glVertex3fv(f->edges[2]->sl_coords);
					glEnd();
				} else {
					glBegin(GL_QUADS);
					for (i = 0; i < f->numVerts; i++) {
						ssDM_glNormalFast(f->verts[i]->sl_coords, f->edges[i]->sl_coords, f->centroid, f->edges[(i - 1 + f->numVerts) % f->numVerts]->sl_coords);
						glVertex3fv(f->verts[i]->sl_coords);
						glVertex3fv(f->edges[i]->sl_coords);
						glVertex3fv(f->centroid);
						glVertex3fv(f->edges[(i + f->numVerts - 1) % f->numVerts]->sl_coords);
					}
					glEnd();
				}
				if (draw_option == DM_DRAW_OPTION_STIPPLE)
					glDisable(GL_POLYGON_STIPPLE);
			}
		}
	}
}

static void slDM_drawMappedFacesTex(DerivedMesh *dm,
									DMSetDrawOptions setDrawOptions,
									DMCompareDrawOptions compareDrawOptions,
									void *userData) {
	// TODO
	printf("slDM_drawMappedFacesTex\n");
}

/* Only used by non-editmesh types */
static void slDM_drawMappedFacesMat(DerivedMesh *dm,
									void (*setMaterial)(void *userData, int, void *attribs),
									int (*setFace)(void *userData, int index), void *userData) {
	// TODO
	printf("slDM_drawMappedFacesMat\n");
}

static void slDM_drawUVEdges(DerivedMesh *dm) {
	// TODO
	printf("slDM_drawUVEdges\n");
}

static void slDM_drawMappedEdgesInterp(DerivedMesh *dm,
									   DMSetDrawOptions setDrawOptions,
									   DMSetDrawInterpOptions setDrawInterpOptions,
									   void *userData) {
	// TODO
	printf("slDM_drawMappedEdgesInterp\n");
}

static void slDM_drawMappedEdges(DerivedMesh *dm,
								 DMSetDrawOptions setDrawOptions,
								 void *userData) {
	// TODO
	printf("slDM_drawMappedEdges\n");
}

// End of drawing functions

static SLDerivedMesh *getSLDerivedMesh(SLSubSurf *ss,
									   int drawInteriorEdges,
									   int useSubsurfUv,
									   DerivedMesh *source_dm)
{
	SLDerivedMesh *ssdm = MEM_callocN(sizeof(*ssdm), "sldm");
	DerivedMesh *newdm;
	int totvert, totedge, totface;
	int totsubvert, totsubedge, totsubface, totsubloop;
	int numTex, numCol;
	int hasPCol, hasOrigSpace;
	int gridInternalEdges;
	int *polyidx;
	int i;

	ssdm->ss = ss;
	//ssdm->drawInteriorEdges = drawInteriorEdges;
	//ssdm->useSubsurfUv = useSubsurfUv;
	newdm = &(ssdm->dm);
	for (i = 0; i < sizeof(newdm); i++) {
		((char*)newdm)[i] = 0; // Nulling everything
	}

	//DM_init_funcs(newdm); // Sets some default functions we don't care/need to overload.

	totvert = ss->numVerts;
	totedge = ss->numEdges;
	totface = ss->numFaces;

	totsubvert = SL_giveTotalNumberOfSubVerts(ss);
	totsubedge = SL_giveTotalNumberOfSubEdges(ss);
	totsubface = SL_giveTotalNumberOfSubFaces(ss);
	totsubloop = SL_giveTotalNumberOfSubLoops(ss);

	DM_from_template(newdm, source_dm, DM_TYPE_CCGDM, // TODO Change type; (general subsurf-type perhaps?)
					totsubvert,
					totsubedge,
					totsubface, // faces
					totsubloop, // loops TODO: unsure
					totsubface); // polys TODO: unsure

	CustomData_free_layer_active(&newdm->polyData, CD_NORMAL,
								 newdm->numPolyData);

	numTex = CustomData_number_of_layers(&newdm->loopData, CD_MLOOPUV);
	numCol = CustomData_number_of_layers(&newdm->loopData, CD_MLOOPCOL);
	hasPCol = CustomData_has_layer(&newdm->loopData, CD_PREVIEW_MLOOPCOL);
	hasOrigSpace = CustomData_has_layer(&newdm->loopData, CD_ORIGSPACE_MLOOP);

	if (
		(numTex && CustomData_number_of_layers(&newdm->faceData, CD_MTFACE) != numTex)  ||
		(numCol && CustomData_number_of_layers(&newdm->faceData, CD_MCOL) != numCol)    ||
		(hasPCol && !CustomData_has_layer(&newdm->faceData, CD_PREVIEW_MCOL))           ||
		(hasOrigSpace && !CustomData_has_layer(&newdm->faceData, CD_ORIGSPACE)) )
	{
		CustomData_from_bmeshpoly(&newdm->faceData,
								  &newdm->polyData,
								  &newdm->loopData,
								  totsubface);
	}

	/* We absolutely need that layer, else it's no valid tessellated data! */
	polyidx = CustomData_add_layer(&newdm->faceData, CD_POLYINDEX, CD_CALLOC,
								   NULL, totsubface);
	for (i = 0; i < totsubface; i++) polyidx[i] = i; // TODO: Check this, i think its OK. (it should just be 1 face per poly)

	ssdm->faceFlags = MEM_callocN(sizeof(DMFlagMat) * totface, "faceFlags");

	// These functions are straight forward...
	newdm->getMinMax = slDM_getMinMax;
	newdm->getNumVerts = slDM_getNumVerts;
	newdm->getNumEdges = slDM_getNumEdges;
	newdm->getNumTessFaces = slDM_getNumTessFaces;
	newdm->getNumLoops = slDM_getNumLoops;
	newdm->getNumPolys = slDM_getNumPolys;

	// Individual item access. This would be as slow as just fetching the array itself, no good way around that.
	// TODO: Implement these later as they are not directly needed to try out the algorithm
	newdm->getVert = NULL;
	newdm->getEdge = NULL;
	newdm->getTessFace = NULL;
	newdm->getVertCo = NULL;
	newdm->getVertNo = NULL;

	newdm->copyVertArray = slDM_copyFinalVertArray;
	newdm->copyEdgeArray = slDM_copyFinalEdgeArray;
	newdm->copyLoopArray = slDM_copyFinalLoopArray;
	newdm->copyPolyArray = slDM_copyFinalPolyArray;
	newdm->copyTessFaceArray = slDM_copyFinalTessFaceArray;

	// Not sure how to deal with these;
	newdm->getVertData = slDM_get_vert_data;
	newdm->getEdgeData = slDM_get_edge_data;
	newdm->getTessFaceData = slDM_get_tessface_data;
	// or these
	newdm->getVertDataArray = slDM_get_vert_data_layer;
	newdm->getEdgeDataArray = slDM_get_edge_data_layer;
	newdm->getTessFaceDataArray = slDM_get_tessface_data_layer;

	// Grids aren't used here (this part is to CCG-centric, shouldn't be in derived mesh)
	newdm->getNumGrids = NULL;
	newdm->getGridSize = NULL;
	newdm->getGridData = NULL;
	newdm->getGridAdjacency = NULL;
	newdm->getGridOffset = NULL;
	newdm->getGridKey = NULL;
	newdm->getGridFlagMats = NULL;
	newdm->getGridHidden = NULL;
	
	newdm->getPolyMap = NULL; // TODO: No idea what this is for
	newdm->getPBVH = slDM_getPBVH; // TODO: No idea what this is for
	newdm->getVertCos = slDM_getVertCos; // TODO: Not sure about this one either.

	newdm->calcNormals = slDM_calcNormals;
	newdm->recalcTessellation = slDM_recalcTessellation;

	newdm->release = slDM_release;
	
	// More unknown...
	newdm->foreachMappedVert = NULL;
	newdm->foreachMappedEdge = NULL;
	newdm->foreachMappedFaceCenter = NULL;

	// Drawing stuff;
	newdm->drawVerts = slDM_drawVerts;
	newdm->drawEdges = slDM_drawEdges;
	newdm->drawLooseEdges = slDM_drawLooseEdges;
	newdm->drawFacesSolid = slDM_drawFacesSolid;
	newdm->drawFacesTex = slDM_drawFacesTex;
	newdm->drawFacesGLSL = slDM_drawFacesGLSL;
	newdm->drawMappedFaces = slDM_drawMappedFaces;
	newdm->drawMappedFacesTex = slDM_drawMappedFacesTex;
	newdm->drawMappedFacesGLSL = slDM_drawMappedFacesGLSL;
	newdm->drawMappedFacesMat = slDM_drawMappedFacesMat;
	newdm->drawUVEdges = slDM_drawUVEdges;
	newdm->drawMappedEdgesInterp = slDM_drawMappedEdgesInterp;
	newdm->drawMappedEdges = slDM_drawMappedEdges;

	// TODO: deal with uv's..
#if 0
	if (useSubsurfUv) {
		CustomData *ldata = &newdm->loopData;
		CustomData *dmldata = &source_dm->loopData;
		int numlayer = CustomData_number_of_layers(ldata, CD_MLOOPUV);
		int dmnumlayer = CustomData_number_of_layers(dmldata, CD_MLOOPUV);

		for (i = 0; i < numlayer && i < dmnumlayer; i++)
			set_subsurf_uv(ss, source_dm, newdm, i);
	}
#endif

	/* All tessellated CD layers were updated! */
	newdm->dirty &= ~DM_DIRTY_TESS_CDLAYERS;
	return ssdm;
}
