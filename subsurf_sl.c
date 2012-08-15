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
	DerivedMesh dm; // Output derived mesh.
	SLSubSurf *ss;
	int drawInteriorEdges;
} SLDerivedMesh;

/*
 * Entrypoint for modifier. Should return the output DM based on the input DM.
 */
struct DerivedMesh *sl_subsurf_make_derived_from_derived(
	struct DerivedMesh *input,
	struct SubsurfModifierData *smd,
	float (*vertCos)[3],
	SubsurfFlags flags)
{
	SLSubSurf *ss;
	int levels;
	int smoothing = smd->subdivType == ME_SL_SUBSURF;
	//SLFlags useAging = smd->flags & eSubsurfModifierFlag_DebugIncr ? CCG_USE_AGING : 0;
	int useSubsurfUv = smd->flags & eSubsurfModifierFlag_SubsurfUv;
	int drawInteriorEdges = !(smd->flags & eSubsurfModifierFlag_ControlEdges);
	DerivedMesh *result;
	
	// We don't use caches for this modifier
	if (smd->emCache) {
		smd->cacheFree(smd->emCache);
		smd->emCache = NULL;
	}
	if (smd->mCache) {
		smd->cacheFree(smd->mCache);
		smd->mCache = NULL;
	}

	if (flags & SUBSURF_USE_RENDER_PARAMS) {
		levels = (smd->modifier.scene) ? get_render_subsurf_level(&smd->modifier.scene->r, smd->renderLevels) : smd->renderLevels;
	}
	else /*if (flags & SUBSURF_FOR_EDIT_MODE)*/ {
		levels = (smd->modifier.scene) ? get_render_subsurf_level(&smd->modifier.scene->r, smd->levels) : smd->levels;
	}
	if (levels == 0)
		return input;
	
	ss = SL_SubSurf_new(smoothing, input, vertCos);
	result = SL_SubSurf_constructOutput(ss);
	SL_processSync(ss); // Actually computes coordinates and such.
	// TODO: NOW WE LEAK MEMORY! FIXME FIXME FIXME Need to use the SLDerivedMesh and overload the release function
	return result;
}

#if 0
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
/*static void slDM_copyFinalVertArray(DerivedMesh *dm, MVert *mvert) {
	int i;
	SLSubSurf *ss = ((SLDerivedMesh*)dm)->ss;
	MVert *overt = ss->o_vert;
	for (i = 0; i < SL_giveTotalNumberOfSubVerts(ss); i++) {
		mvert[i].bweight = overt[i].bweight;
		mvert[i].flag = overt[i].flag;
		copy_v3_v3(mvert[i].co = overt[i].co);
		copy_v3_v3_short(mvert[i].no = overt[i].no);
	}
}*/
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
	SL_copyNewPolys(sldm->ss, mpoly);
}
static void slDM_copyFinalTessFaceArray(DerivedMesh *dm, MFace *mface) {
	SLDerivedMesh *sldm = (SLDerivedMesh*)dm;
	SL_copyNewTessFaces(sldm->ss, mface);
}

// Unsure about these three
static void *slDM_get_vert_data_layer(DerivedMesh *dm, int type) {
	if (type == CD_ORIGINDEX) {
		/* create origindex on demand to save memory */
		SLDerivedMesh *ssdm = (SLDerivedMesh *)dm;
		SLSubSurf *ss = ssdm->ss;
		int *origindex;
		int a, i, tot;

		/* Avoid re-creation if the layer exists already */
		origindex = DM_get_vert_data_layer(dm, CD_ORIGINDEX);
		if (origindex) {
			return origindex;
		}

		DM_add_vert_layer(dm, CD_ORIGINDEX, CD_CALLOC, NULL);
		origindex = DM_get_vert_data_layer(dm, CD_ORIGINDEX);

		// original vertices are at the beginning
		a = 0;
		for (i = 0; i < ss->numVerts; i++, a++) {
			origindex[a] = i;
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
		int a, i, tot;

		/* Avoid re-creation if the layer exists already */
		origindex = DM_get_edge_data_layer(dm, CD_ORIGINDEX);
		if (origindex) {
			return origindex;
		}

		DM_add_edge_layer(dm, CD_ORIGINDEX, CD_CALLOC, NULL);
		origindex = DM_get_edge_data_layer(dm, CD_ORIGINDEX);

		a = 0;
		for (i = 0; i < ss->numEdges; i++, a++) {
			origindex[a++] = i;
			origindex[a++] = i;
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
		SLDerivedMesh *sldm = (SLDerivedMesh *)dm;
		SLSubSurf *ss = sldm->ss;
		int *origindex;
		int a, i, j;

		/* Avoid re-creation if the layer exists already */
		origindex = DM_get_tessface_data_layer(dm, CD_ORIGINDEX);
		if (origindex) {
			return origindex;
		}

		DM_add_tessface_layer(dm, CD_ORIGINDEX, CD_CALLOC, NULL);
		origindex = DM_get_tessface_data_layer(dm, CD_ORIGINDEX);

		a = 0;
		for (i = 0; i < ss->numFaces; i++) {
			MPoly *poly = &ss->mpoly[i];
			for (j = 0; j < SL_giveNumberOfInternalFaces(poly); j++, a++)
				origindex[a] = i;
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
	int i;
	for (i = 0; i < SL_giveTotalNumberOfSubVerts(ss); i++) {
		copy_v3_v3(cos[i], ss->o_vert[i].co);
	}
}

static struct PBVH *slDM_getPBVH(Object *ob, DerivedMesh *dm)
{
	// TODO: I have no idea about any of this code. The grids and all that crap just seems relevant to CCG.
	SLDerivedMesh *ssdm = (SLDerivedMesh *)dm;
	
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
	int i;
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;
	glBegin(GL_POINTS);
	for (i = 0; i < SL_giveTotalNumberOfSubVerts(ss); i++) {
		glVertex3fv(ss->o_vert[i].co);
	}
	glEnd();
}

static void slDM_drawEdges(DerivedMesh *dm, int drawLooseEdges, int drawAllEdges) {
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;
	int i,j;

	for (i = 0; i < ss->numEdges; i++) {
		MEdge *e;
		if (!drawLooseEdges && ss->edge2poly[i]->count == 0)
			continue;
		
		e = &ss->medge[i];
		
		if (!drawAllEdges && !(e->flag & ME_EDGEDRAW))
			continue;
		
		
		// Or what aging is
		/*if (useAging && !(G.f & G_BACKBUFSEL)) {
			i n*t ageCol = 255 - ccgSubSurf_getEdgeAge(ss, e) * 4;
			glColor3ub(0, ageCol > 0 ? ageCol : 0, 0);
		}*/

		glBegin(GL_LINE_STRIP);
		glVertex3fv(ss->o_vert[e->v1].co);
		glVertex3fv(ss->o_vert[ss->numVerts + i].co);
		glVertex3fv(ss->o_vert[e->v2].co);
		glEnd();

	}

	/*if (useAging && !(G.f & G_BACKBUFSEL)) {
		glColor3ub(0, 0, 0);
	}*/

	if (ssdm->drawInteriorEdges) {
		for (i = 0; i < ss->numFaces ; i++) {
			MPoly *poly = &ss->mpoly[i];
			MLoop *loop = &ss->mloop[poly->loopstart];
			
			// TODO:
			if (poly->totloop == 3) { // Triangles are split differently from quads and ngons
				glBegin(GL_LINE_LOOP);
				glVertex3fv(ss->o_vert[loop[0].e + ss->numVerts].co);
				glVertex3fv(ss->o_vert[loop[1].e + ss->numVerts].co);
				glVertex3fv(ss->o_vert[loop[1].e + ss->numVerts].co);
				glEnd();
			}
			else {
				float *centroid = ss->o_vert[ss->poly2vert[i]].co;
				glBegin(GL_LINES);
				for (j = 0; j < poly->totloop; j++) {
					glVertex3fv(centroid);
					glVertex3fv(ss->o_vert[loop[j].e + ss->numVerts].co);
				}
				glEnd();
			}
		}
	}
}

static void slDM_drawLooseEdges(DerivedMesh *dm) {
	int i;
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;

	for (i = 0; i < ss->numEdges; i++) {
		if (ss->edge2poly[i]->count == 0) {
			MEdge *e = &ss->medge[i];
			glBegin(GL_LINE_STRIP);
			glVertex3fv(ss->o_vert[e->v1].co);
			glVertex3fv(ss->o_vert[ss->numVerts + i].co);
			glVertex3fv(ss->o_vert[e->v2].co);
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

#if 0
static void slDM_drawFacesSolid(DerivedMesh *dm, float (*partial_redraw_planes)[4], int fast, DMSetMaterial setMaterial) {
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;
	int i, j;
	int drawcurrent = 0, matnr = -1, shademodel = -1;
	
	MVert *o_vert = ss->o_vert; // For convenience.

	for (i = 0; i < ss->numFaces; i++) {
		MPoly *poly = &ss->mpoly[i];
		MLoop *loop = &ss->mloop[i];
		int new_matnr, new_shademodel;

		new_shademodel = (poly->flag & ME_SMOOTH) ? GL_SMOOTH : GL_FLAT;
		new_matnr = poly->mat_nr;

		if (shademodel != new_shademodel || matnr != new_matnr) {
			matnr = new_matnr;
			shademodel = new_shademodel;

			drawcurrent = setMaterial(matnr + 1, NULL);

			glShadeModel(shademodel);
		}

		if (!drawcurrent)
			continue;

		// TODO: Smooth version
		if (shademodel == GL_SMOOTH) {
			if (poly->totloop == 3) {
				MVert *v1 = &o_vert[loop[0].v], *v2 = &o_vert[loop[1].v], *v3 = &o_vert[loop[2].v];
				MVert *e1 = &o_vert[loop[0].e + ss->numVerts], 
					  *e2 = &o_vert[loop[1].v + ss->numVerts],
					  *e3 = &o_vert[loop[2].v + ss->numVerts];

				// 1x3 strip + 1 single triangle, not worth it (probably?)
				glBegin(GL_TRIANGLES);
				glNormal3fv(v1->no);glVertex3fv(v1->co);
				glNormal3fv(e1->no);glVertex3fv(e1->co);
				glNormal3fv(e3->no);glVertex3fv(e3->co);
				
				glNormal3fv(v2->no);glVertex3fv(v2->co);
				glNormal3fv(e2->no);glVertex3fv(e2->co);
				glNormal3fv(e1->no);glVertex3fv(e1->co);
				
				glNormal3fv(v3->no);glVertex3fv(v3->co);
				glNormal3fv(e3->no);glVertex3fv(e3->co);
				glNormal3fv(e2->no);glVertex3fv(e2->co);
				
				glNormal3fv(e1->no);glVertex3fv(e1->co);
				glNormal3fv(e2->no);glVertex3fv(e2->co);
				glNormal3fv(e3->no);glVertex3fv(e3->co);
				glEnd();
			} else {
				MVert *p = &o_vert[ss->poly2vert[i]];
				glBegin(GL_QUADS);
				for (j = 0; j < poly->totloop; j++) {
					int prevJ = (j + poly->totloop - 1) % poly->totloop;
					MVert *vj = &o_vert[loop[j].v];
					MVert *ej = &o_vert[loop[j].e + ss->numVerts];
					MVert *eprevJ = &o_vert[loop[prevJ].e + ss->numVerts];
					
					glNormal3fv(vj->no);glVertex3fv(vj->co);
					glNormal3fv(ej->no);glVertex3fv(ej->co);
					glNormal3fv(p->no);glVertex3fv(p->co);
					glNormal3fv(eprevJ->no);glVertex3fv(eprevJ->co);
				}
				glEnd();
			}
		}
		else {
			if (poly->totloop == 3) {
				MVert *v1 = &o_vert[loop[0].v], *v2 = &o_vert[loop[1].v], *v3 = &o_vert[loop[2].v];
				MVert *e1 = &o_vert[loop[0].e + ss->numVerts], 
						*e2 = &o_vert[loop[1].v + ss->numVerts],
						*e3 = &o_vert[loop[2].v + ss->numVerts];
				
				// 1x3 strip + 1 single triangle, not worth it (probably?)
				glBegin(GL_TRIANGLES);

				ssDM_glNormalFastTri(v1->co, e1->co, e3->co);
				glVertex3fv(v1->co);
				glVertex3fv(e1->co);
				glVertex3fv(e3->co);

				ssDM_glNormalFastTri(v2->co,e2->co, e1->co);
				glVertex3fv(v2->co);
				glVertex3fv(e2->co);
				glVertex3fv(e1->co);

				ssDM_glNormalFastTri(v3->co, e3->co, e2->co);
				glVertex3fv(v3->co);
				glVertex3fv(e3->co);
				glVertex3fv(e2->co);

				ssDM_glNormalFastTri(e1->co, e2->co, e3->co);
				glVertex3fv(e1->co);
				glVertex3fv(e2->co);
				glVertex3fv(e3->co);
				glEnd();
			} else {
				MVert *p = &o_vert[ss->poly2vert[i]];
				glBegin(GL_QUADS);
				for (j = 0; j < poly->totloop; j++) {
					int prevJ = (j + poly->totloop - 1) % poly->totloop;
					MVert *vj = &o_vert[loop[j].v];
					MVert *ej = &o_vert[loop[j].e + ss->numVerts];
					MVert *eprevJ = &o_vert[loop[prevJ].e + ss->numVerts];

					ssDM_glNormalFast(vj->co, ej->co, p->co, eprevJ->co);
					glVertex3fv(vj->co);
					glVertex3fv(ej->co);
					glVertex3fv(p->co);
					glVertex3fv(eprevJ->co);
				}
				glEnd();
			}
		}
	}
}
#endif

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

#if 0
static void slDM_drawMappedFaces(DerivedMesh *dm,
								 DMSetDrawOptions setDrawOptions,
								 DMSetMaterial setMaterial,
								 DMCompareDrawOptions compareDrawOptions,
								 void *userData, DMDrawFlag flag) {
	// TODO
	SLDerivedMesh *ssdm = (SLDerivedMesh *) dm;
	SLSubSurf *ss = ssdm->ss;
	MCol *mcol = NULL;
	int useColors = flag & DM_DRAW_USE_COLORS;
	int i, j, drawSmooth;
	
	// currently unused -- each original face is handled separately
	(void)compareDrawOptions;
	
	if (useColors) {
		mcol = dm->getTessFaceDataArray(dm, CD_PREVIEW_MCOL);
		if (!mcol)
			mcol = dm->getTessFaceDataArray(dm, CD_MCOL);
	}
	
	
	for (i = 0; i < ss->numFaces; i++) {
		MPoly *poly = &ss->mpoly[i];
		unsigned char *cp = NULL;

		drawSmooth = (flag & DM_DRAW_ALWAYS_SMOOTH) || poly->flag & ME_SMOOTH;
		
		if (mcol) {
			cp = (unsigned char *)mcol;
			mcol += poly->totloop * 4; // 4 is used for what?
		}
		
		{
			DMDrawOption draw_option = DM_DRAW_OPTION_NORMAL;
			draw_option = setMaterial(poly->mat_nr, NULL);
				
			if (draw_option != DM_DRAW_OPTION_SKIP) {
				if (draw_option == DM_DRAW_OPTION_STIPPLE) {
					glEnable(GL_POLYGON_STIPPLE);
					glPolygonStipple(stipple_quarttone);
				}
					
				/* no need to set shading mode to flat because
				 *  normals are already used to change shading */
				glShadeModel(GL_SMOOTH);
				if (poly->totloop == 3) {
					int ev1, ev2, ev3;
					
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
					for (j = 0; j < poly->totloop; j++) {
						ssDM_glNormalFast(f->verts[j]->sl_coords, f->edges[j]->sl_coords, f->centroid, f->edges[(j - 1 + f->numVerts) % f->numVerts]->sl_coords);
						glVertex3fv(f->verts[j]->sl_coords);
						glVertex3fv(f->edges[j]->sl_coords);
						glVertex3fv(f->centroid);
						glVertex3fv(f->edges[(j + f->numVerts - 1) % f->numVerts]->sl_coords);
					}
					glEnd();
				}
				if (draw_option == DM_DRAW_OPTION_STIPPLE)
					glDisable(GL_POLYGON_STIPPLE);
			}
		}
	}
}
#endif

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
	int totsubvert, totsubedge, totsubface, totsubloop;
	int numTex, numCol;
	int hasPCol, hasOrigSpace;
	int *polyidx;
	int i;

	ssdm->ss = ss;
	ssdm->drawInteriorEdges = drawInteriorEdges;
	//ssdm->useSubsurfUv = useSubsurfUv;
	newdm = &(ssdm->dm);
	for (i = 0; i < sizeof(newdm); i++) {
		((char*)newdm)[i] = 0; // Nulling everything
	}

	//DM_init_funcs(newdm); // Sets some default functions we don't care/need to overload.

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

	//newdm->copyVertArray = ; // TODO
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
	//newdm->drawFacesSolid = slDM_drawFacesSolid;
	newdm->drawFacesTex = slDM_drawFacesTex;
	newdm->drawFacesGLSL = slDM_drawFacesGLSL;
	//newdm->drawMappedFaces = slDM_drawMappedFaces;
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
#endif