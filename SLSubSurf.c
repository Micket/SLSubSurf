/*
   DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
   Version 2, December 2004

   Copyright (C) 2004 Sam Hocevar <sam@hocevar.net>

   Everyone is permitted to copy and distribute verbatim or modified
   copies of this license document, and changing it is allowed as long
   as the name is changed.

   DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

   0. You just DO WHAT THE FUCK YOU WANT TO.
 */

// Contributor(s): Mikael Öhman

#include "SLSubSurf.h"
#include "stdlib.h"
#include "BLI_math_vector.h"
#include "BLI_linklist.h"
#include "BLI_ghash.h"
#include "BLI_memarena.h"
#include "MEM_guardedalloc.h"
#include "DNA_meshdata_types.h"
#include "BKE_DerivedMesh.h"
#include "BKE_cdderivedmesh.h"
#include "BKE_mesh.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////
// External helpers

int SL_giveTotalNumberOfSubVerts(SLSubSurf *ss) {
	// Simple things first, corners and interpolated edge verts;
	int i, totNodes = ss->numVerts + ss->numEdges;
	// Then faces, which varies;
	for (i = 0; i < ss->numFaces; i++) {
		if (ss->mpoly[i].totloop > 3) {
			totNodes++;
		}
	}
	return totNodes;
}

int SL_giveTotalNumberOfSubEdges(SLSubSurf *ss) {
	int i, totEdges = ss->numEdges * 2;
	for (i = 0; i < ss->numFaces; i++) {
		totEdges += ss->mpoly[i].totloop;
	}
	return totEdges;
}

int SL_giveTotalNumberOfSubFaces(SLSubSurf *ss) {
	int i, totFaces = 0;
	for (i = 0; i < ss->numFaces; i++) {
		totFaces += (ss->mpoly[i].totloop == 3) ? 4 : ss->mpoly[i].totloop;
	}
	return totFaces;
}

int SL_giveTotalNumberOfSubLoops(SLSubSurf *ss) {
	int i, totLoops = 0;
	for (i = 0; i < ss->numFaces; i++) {
		totLoops += ss->mpoly[i].totloop * 4;
	}
	return totLoops;
}

/////////////////////////////////////////////////////////////

void SL_copyNewLoops(SLSubSurf *ss, MLoop *mloops) {
	int i, j, k, n, prevJ, subEdgeNext, subEdgePrev, faceEdgeStartIdx, faceVertStartIdx;
	MPoly *poly;
	MLoop *loop;
	MEdge *eNext, *ePrev;

	k = 0;
	faceEdgeStartIdx = ss->numEdges*2; // Number the new internal edges as we go
	faceVertStartIdx = ss->numVerts + ss->numEdges; // Number the new internal verts as we go
	for (i = 0; i < ss->numFaces; i++) {
		poly = &ss->mpoly[i];
		loop = &ss->mloop[poly->loopstart];
		n = poly->totloop;
		if (n == 3) {
			// First the corner triangles
			for (j = 0; j < 3; j++) {
				prevJ = (j + 2) % 3;
				eNext = &ss->medge[loop[j].e];
				ePrev = &ss->medge[loop[prevJ].e];

				// Check ordering of edge (to determine which sub-edge should be used)
				subEdgeNext = eNext->v1 != loop[j].v;
				subEdgePrev = ePrev->v1 != loop[j].v;

				// Corner to next edge;
				mloops[k+0].v = loop[j].v;
				mloops[k+0].e = loop[j].e*2 + subEdgeNext;
				// next edge node to internal edge j;
				mloops[k+1].v = loop[j].e + ss->numVerts;
				mloops[k+1].e = faceEdgeStartIdx + j;
				// previous edge node to previous edge;
				mloops[k+2].v = loop[prevJ].e + ss->numVerts;
				mloops[k+2].e = loop[prevJ].e*2 + subEdgePrev;
				k += 3;
			}

			// Last poly is the center polygon, only internal edges and edge nodes
			mloops[k+0].v = loop[2].e + ss->numVerts;
			mloops[k+0].e = faceEdgeStartIdx;
			mloops[k+1].v = loop[0].e + ss->numVerts;
			mloops[k+1].e = faceEdgeStartIdx+1;
			mloops[k+2].v = loop[1].e + ss->numVerts;
			mloops[k+2].e = faceEdgeStartIdx+2;
			k += 3;
			
			faceEdgeStartIdx += 3;
			//faceVertStartIdx += 0;
		} else {
			for (j = 0; j < n; j++) { // Loop over each sub-quad;
				// Unsure if i negative values work, adding numVerts just in case.
				prevJ = (j + n - 1) % n;
				eNext = &ss->medge[loop[j].e];
				ePrev = &ss->medge[loop[prevJ].e];

				// Check ordering of edge (to determine which sub-edge should be used)
				subEdgeNext = eNext->v1 != loop[j].v;
				subEdgePrev = ePrev->v1 != loop[j].v;

				// Starting from the corner node
				mloops[k+0].v = loop[j].v;
				mloops[k+0].e = loop[j].e*2 + subEdgeNext;
				// go to next edge
				mloops[k+1].v = loop[j].e + ss->numVerts;
				mloops[k+1].e = faceEdgeStartIdx + j;
				// then to midpoint
				mloops[k+2].v = faceVertStartIdx;
				mloops[k+2].e = faceEdgeStartIdx + prevJ;
				// then to previous edge,
				mloops[k+3].v = loop[prevJ].e + ss->numVerts;
				mloops[k+3].e = loop[prevJ].e*2 + subEdgePrev;
				k += 4;
				//printf("mloops.v = %e, %e, %e, %e\n", mloops[k+0].v, mloops[k+1].v, mloops[k+2].v, mloops[k+3].v);
			}
			faceEdgeStartIdx += n;
			faceVertStartIdx += 1;
		}
	}
}

void SL_copyNewPolys(SLSubSurf *ss, MPoly *mpolys) {
	int i, j, k = 0, m = 0, n;
	for (i = 0; i < ss->numFaces; i++) {
		MPoly *poly = &ss->mpoly[i];
		n = poly->totloop;
		if (n == 3) {
			for (j = 0; j < 4; j++) { // note! 4 faces
				mpolys[m].loopstart = k;
				mpolys[m].totloop = 3;
				mpolys[m].mat_nr = poly->mat_nr;
				mpolys[m].flag = poly->flag;
				m += 1;
				k += 3;
			}
		} else {
			for (j = 0; j < n; j++) {
				mpolys[m].loopstart = k;
				mpolys[m].totloop = 4;
				mpolys[m].mat_nr = poly->mat_nr;
				mpolys[m].flag = poly->flag;
				m += 1;
				k += 4;
			}
		}
	}
}

void SL_copyNewTessFaces(SLSubSurf *ss, MFace *mfaces)
{
	int i, j, prevJ, k, n, faceVertStartIdx;
	MPoly *poly;
	MLoop *loop;

	k = 0;
	faceVertStartIdx = ss->numVerts + ss->numEdges; // Number the new internal verts as we go
	for (i = 0; i < ss->numFaces; i++) {
		poly = &ss->mpoly[i];
		n = poly->loopstart;
		loop = &ss->mloop[n];
		n = poly->totloop;
		if (n == 3) {
			// First the corner triangles
			for (j = 0; j < 3; j++) {
				prevJ = (j + 2) % 3;

				mfaces[k].v1 = loop[j].v;
				mfaces[k].v2 = loop[j].e + ss->numVerts;
				mfaces[k].v3 = loop[prevJ].e + ss->numVerts;
				mfaces[k].v4 = 0;
				mfaces[k].mat_nr = poly->mat_nr;
				mfaces[k].flag = poly->flag;
				k += 1;
			}

			mfaces[k].v1 = loop[2].e + ss->numVerts;
			mfaces[k].v2 = loop[0].e + ss->numVerts;
			mfaces[k].v3 = loop[1].e + ss->numVerts;
			mfaces[k].v4 = 0;
			mfaces[k].mat_nr = poly->mat_nr;
			mfaces[k].flag = poly->flag;

			k += 1;
		} else {
			for (j = 0; j < n; j++) { // Loop over each sub-quad;
				// Unsure if i negative values work, adding numVerts just in case.
				prevJ = (j + n - 1) % n;

				mfaces[k].v1 = loop[j].v;
				mfaces[k].v2 = loop[j].e + ss->numVerts;
				mfaces[k].v3 = faceVertStartIdx;
				mfaces[k].v4 = loop[prevJ].e + ss->numVerts;
				mfaces[k].mat_nr = poly->mat_nr;
				mfaces[k].flag = poly->flag;
				k += 1;
			}
			faceVertStartIdx++;
		}
	}
}

void SL_copyNewEdges(SLSubSurf *ss, MEdge *medges) {
	/* For triangles;
    (0)
	 |\
	 | \
	 |  \ e0
  e2 |_0_\
	 |\  |\
	 |2\ |1\
  (2)|__\|__\ (1)
       e1

	And ngons, the natural ordering is used
	*/
	int i, j, k, n, faceVertStartIdx;
	// First the original edges
	k = 0;
	for (i = 0; i < ss->numEdges; i++) {
		MEdge *edge = &ss->medge[i];
		medges[k+0].v1 = edge->v1;
		medges[k+0].v2 = ss->numVerts + i;
		medges[k+0].crease = edge->crease;
		medges[k+0].bweight = edge->bweight;
		medges[k+1].v1 = ss->numVerts + i;
		medges[k+1].v2 = edge->v2;
		medges[k+1].crease = edge->crease;
		medges[k+1].bweight = edge->bweight;
		//printf("medges[%d] = {%d, %d}\n", i, medges[i+0].v1, medges[i+0].v2);
		//printf("medges[%d] = {%d, %d}\n", i+1, medges[i+0].v1, medges[i+0].v2);
		k += 2;
	}
	// Then the faces
	k = 2*ss->numEdges;
	faceVertStartIdx = ss->numVerts + ss->numEdges; // Number the new internal verts as we go
	for (i = 0; i < ss->numFaces; i++) {
		MPoly *poly = &ss->mpoly[i];
		MLoop *loop = &ss->mloop[poly->loopstart];
		n = poly->totloop;
		if (n == 3) {
			medges[k+0].v1 = ss->numVerts + loop[0].e;
			medges[k+0].v2 = ss->numVerts + loop[2].e;
			medges[k+1].v1 = ss->numVerts + loop[0].e;
			medges[k+1].v2 = ss->numVerts + loop[1].e;
			medges[k+2].v1 = ss->numVerts + loop[1].e;
			medges[k+2].v2 = ss->numVerts + loop[2].e;
			//printf("medges[%d] = {%d, %d}\n", i, medges[k+0].v1, medges[k+0].v2);
			//printf("medges[%d] = {%d, %d}\n", i, medges[k+1].v1, medges[k+1].v2);
			//printf("medges[%d] = {%d, %d}\n", i, medges[k+2].v1, medges[k+2].v2);
			k += 3;
		} else {
			for (j = 0; j < n; j++) {
				medges[k].v1 = ss->numVerts + loop[j].e;
				medges[k].v2 = faceVertStartIdx;
				//printf("medges[%d] = {%d, %d}\n", i, medges[k].v1, medges[k].v2);
				k++;
			}
			faceVertStartIdx++;
		}
	}
}


MINLINE float  *_origCoord(SLSubSurf *ss, int vert) {
	return ss->vertexCos ? ss->vertexCos[vert] : ss->mvert[vert].co;
}

/////////////////////////////////////////////////////////////

#ifdef _OPENMP
static omp_lock_t alloc_lock;

static void OMP_lock_alloc_thread(void) {
	omp_set_lock(&alloc_lock);
}

static void OMP_unlock_alloc_thread(void) {
	omp_unset_lock(&alloc_lock);
}
#endif

void SL_constructMaps(SLSubSurf *ss) {
	int i, idx, numLoops;
	ss->poly2vert = MEM_callocN(sizeof(int)*ss->numFaces, "sl poly2vert");
	numLoops = ss->input->getNumLoops(ss->input);
	
	#ifdef _OPENMP
	omp_init_lock(&alloc_lock);
	MEM_set_lock_callback(OMP_lock_alloc_thread, OMP_unlock_alloc_thread);
	#endif
	
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			idx = ss->numVerts + ss->numEdges;
			for (i = 0; i < ss->numFaces; i++) {
				if (ss->mpoly[i].totloop > 3) { // Just ignoring triangles, they will not be accessed.
					ss->poly2vert[i] = idx;
					idx++;
				}
			}
		}
		// Necessary backwards mappings (could potentially be done without these, but that would make partial updates impossible (not that partial updates are supported yet))
		#pragma omp section
		create_vert_poly_map(&ss->vert2poly, &ss->vert2poly_mem,
							 ss->mpoly, ss->mloop, ss->numVerts, ss->numFaces, numLoops);
		#pragma omp section
		create_vert_edge_map(&ss->vert2edge, &ss->vert2edge_mem,
							 ss->medge, ss->numVerts, ss->numEdges);
		#pragma omp section
		create_edge_poly_map(&ss->edge2poly, &ss->edge2poly_mem,
							 ss->mpoly, ss->mloop, ss->numEdges, ss->numFaces, numLoops);
	}
	#ifdef _OPENMP
	omp_destroy_lock(&alloc_lock);
	#endif
}

void SL_freeMaps(SLSubSurf *ss) {
	MEM_freeN(ss->poly2vert);
	MEM_freeN(ss->vert2poly);
	MEM_freeN(ss->vert2poly_mem);
	MEM_freeN(ss->vert2edge);
	MEM_freeN(ss->vert2edge_mem);
	MEM_freeN(ss->edge2poly);
	MEM_freeN(ss->edge2poly_mem);
}

SLSubSurf* SL_SubSurf_new(int smoothing, DerivedMesh *input, float (*vertexCos)[3]) {
	SLSubSurf *ss = (SLSubSurf*)MEM_callocN(sizeof(SLSubSurf), "slsubsurf");
	ss->smoothing = smoothing;
	ss->numVerts = ss->numEdges = ss->numFaces = 0;
	
	// For convenient access.
	ss->input = input;
	ss->vertexCos = vertexCos;

	ss->mvert = input->getVertArray(input);
	ss->medge = input->getEdgeArray(input);
	ss->mpoly = input->getPolyArray(input);
	ss->mloop = input->getLoopArray(input);
	
	ss->numVerts = input->getNumVerts(input);
	ss->numEdges = input->getNumEdges(input);
	ss->numFaces = input->getNumPolys(input);

	SL_constructMaps(ss);
	return ss;
}

// HACK: Copied from cdderivedmesh.c. I want to use the goodies from CDDM, but i need to allocate memory to fit the SLSurfSurf into the structure as well.
// This is one of the many things that beg to be upgraded to C++ and proper OO
typedef struct {
	DerivedMesh dm;
	
	/* these point to data in the DerivedMesh custom data layers,
	 * they are only here for efficiency and convenience **/
	MVert *mvert;
	MEdge *medge;
	MFace *mface;
	MLoop *mloop;
	MPoly *mpoly;
	
	/* Cached */
	struct PBVH *pbvh;
	int pbvh_draw;
	
	/* Mesh connectivity */
	MeshElemMap *pmap;
	int *pmap_mem;
} CDDerivedMesh;

struct SLDerivedMesh {
	CDDerivedMesh cddm; // Output derived mesh.
	SLSubSurf *ss;
	int drawInteriorEdges; // is it worth the bother? I'm not going to repeat everything in CDDM for this.
};

void SL_SubSurf_free(SLSubSurf *ss) {
	SL_freeMaps(ss);
	MEM_freeN(ss);
}

// A few functions that should be replaced from the CDDM default functions;
static void slDM_release(DerivedMesh *dm) {
	SLDerivedMesh *sldm = (SLDerivedMesh*)dm;
	CDDerivedMesh *cddm = (CDDerivedMesh *)dm;
	
	if (DM_release(dm)) {
		if (cddm->pmap) MEM_freeN(cddm->pmap);
		if (cddm->pmap_mem) MEM_freeN(cddm->pmap_mem);
		
		SL_SubSurf_free(sldm->ss);
		MEM_freeN(sldm);
	}
}

//static void slDM_doNothing(DerivedMesh *UNUSED(dm)) {}

#if 1
void slDM_constructTessFaces(DerivedMesh *output) {
	SLSubSurf *ss = ((SLDerivedMesh*)output)->ss;
	int i, numFaces, *polyidx;
	MFace *faces = output->getTessFaceArray(output);
	SL_copyNewTessFaces(ss, faces);
	
	numFaces = output->getNumPolys(output); //SL_giveTotalNumberOfSubFaces(ss);
	
	if (!CustomData_get_layer(&output->faceData, CD_POLYINDEX))
		CustomData_add_layer(&output->faceData, CD_POLYINDEX, CD_CALLOC, NULL, numFaces);
	/* We absolutely need that layer, else it's no valid tessellated data! */
	polyidx = CustomData_get_layer(&output->faceData, CD_POLYINDEX);
	for (i = 0; i < numFaces; i++) polyidx[i] = i; // TODO: Check this, i think its OK. (it should just be 1 face per poly)
}

void slDM_calcNormals(DerivedMesh *output) {
	int i, numPolys, numVerts;
	float (*vertNo)[3];
	float (*polyNo)[3];
	MPoly *poly;
	MLoop *loop;
	MVert *vert;

	numPolys = output->getNumPolys(output);
	numVerts = output->getNumVerts(output);
	
	poly = output->getPolyArray(output);
	loop = output->getLoopArray(output);
	vert = output->getVertArray(output);
	
	vertNo = MEM_mallocN(sizeof(float)*numVerts*3, "temp vert normals");
	
	if (!CustomData_get_layer(&output->polyData, CD_NORMAL))
		CustomData_add_layer(&output->polyData, CD_NORMAL, CD_CALLOC, NULL, numPolys);
	polyNo = CustomData_get_layer(&output->polyData, CD_NORMAL);
	
	#pragma omp for
	for (i = 0; i < numPolys; i++) {
		MPoly *p = &poly[i];
		MLoop *l = &loop[p->loopstart];
		mesh_calc_poly_normal(p, l, vert, polyNo[i]);
	}
	
	CustomData_add_layer(&output->faceData, CD_NORMAL, CD_REFERENCE, polyNo, numPolys);

	#pragma omp for
	for (i = 0; i < numVerts; i++) {
		zero_v3(vertNo[i]);
	}

	for (i = 0; i < numPolys; i++) {
		int j;
		MPoly *p = &poly[i];
		MLoop *l = &loop[p->loopstart];
		
		for (j = 0; j < p->totloop; j++) {
			add_v3_v3(vertNo[l[j].v], polyNo[i]);
		}
	}
	
	#pragma omp for
	for (i = 0; i < numVerts; i++) {
		normalize_v3(vertNo[i]);
		normal_float_to_short_v3(vert[i].no, vertNo[i]);
	}
	
	MEM_freeN(vertNo);
}
#endif

DerivedMesh *SL_SubSurf_constructOutput(SLSubSurf *ss) {
	int numVerts, numEdges, numFaces, numLoops;
	int i, j, a, *input_index, *output_index;
	SLDerivedMesh *output;
	DerivedMesh *output_dm;
	
	numVerts = SL_giveTotalNumberOfSubVerts(ss);
	numEdges = SL_giveTotalNumberOfSubEdges(ss);
	numFaces = SL_giveTotalNumberOfSubFaces(ss);
	numLoops = SL_giveTotalNumberOfSubLoops(ss);
	// Outside is responsible for freeing this memory.
	//output = CDDM_from_template(ss->input, numVerts, numEdges, numFaces, numLoops, numFaces);
	// TODO: Had to duplicate entire template function just to change the cdDM_create command.
	output = MEM_callocN(sizeof(*output), "SLDM output");
	output->ss = ss;
	output_dm = &output->cddm.dm;
	CDDM_init_funcs(output_dm);
	
	// TODO: Start of what should be in CDDM_from_template
	//CDDM_from_template(ss->input, output, numVerts, numEdges, numFaces, numLoops, numFaces);
	/* ensure these are created if they are made on demand */
	ss->input->getVertDataArray(ss->input, CD_ORIGINDEX);
	ss->input->getEdgeDataArray(ss->input, CD_ORIGINDEX);
	ss->input->getTessFaceDataArray(ss->input, CD_ORIGINDEX); // TODO: Can i move this to constructTessFaces ?
	
	/* this does a copy of all non mvert/medge/mface layers */
	DM_from_template(output_dm, ss->input, DM_TYPE_CDDM, numVerts, numEdges, numFaces, numLoops, numFaces);
	
	/* now add mvert/medge/mface layers */
	CustomData_add_layer(&output_dm->vertData, CD_MVERT, CD_CALLOC, NULL, numVerts);
	CustomData_add_layer(&output_dm->edgeData, CD_MEDGE, CD_CALLOC, NULL, numEdges);
	CustomData_add_layer(&output_dm->faceData, CD_MFACE, CD_CALLOC, NULL, numFaces);
	CustomData_add_layer(&output_dm->loopData, CD_MLOOP, CD_CALLOC, NULL, numLoops);
	CustomData_add_layer(&output_dm->polyData, CD_MPOLY, CD_CALLOC, NULL, numFaces);
	
	if (!CustomData_get_layer(&output_dm->vertData, CD_ORIGINDEX))
		CustomData_add_layer(&output_dm->vertData, CD_ORIGINDEX, CD_CALLOC, NULL, numVerts);
	if (!CustomData_get_layer(&output_dm->edgeData, CD_ORIGINDEX))
		CustomData_add_layer(&output_dm->edgeData, CD_ORIGINDEX, CD_CALLOC, NULL, numEdges);
	if (!CustomData_get_layer(&output_dm->faceData, CD_ORIGINDEX))
		CustomData_add_layer(&output_dm->faceData, CD_ORIGINDEX, CD_CALLOC, NULL, numFaces);

	output->cddm.mvert = CustomData_get_layer(&output_dm->vertData, CD_MVERT);
	output->cddm.medge = CustomData_get_layer(&output_dm->edgeData, CD_MEDGE);
	output->cddm.mface = CustomData_get_layer(&output_dm->faceData, CD_MFACE);
	output->cddm.mloop = CustomData_get_layer(&output_dm->loopData, CD_MLOOP);
	output->cddm.mpoly = CustomData_get_layer(&output_dm->polyData, CD_MPOLY);
	// End of CDDM_from_template copy	
	
	// Replace some functions from CDDM
	//output_dm->recalcTessellation = slDM_constructTessFaces; // I could perhaps optimize some stuff here.
	output_dm->release = slDM_release; // Must be replaced since we have our own custom structure (!)
	// TODO: Possibly replace any other functions from default CDDM for performance reasons.
	// Just calculate them here, once (No good for multi-res applications, where normals aren't always necessary (but time consuming)
	//SL_calcNormals(output_dm);
	output_dm->calcNormals = slDM_calcNormals; // Less complex than what CDDM does, but since CC gets away with this, I'll add it.
	
	// This was done in in the CCG code, so i suppose I'll use that as well.
	CustomData_free_layer_active(&output_dm->polyData, CD_NORMAL, output_dm->numPolyData);

	// We can directly fill the structural part, which doesn't actually need the new vertex coordinates.
	// The time spent doing this is in some cases actually quite significant, so lets OMP it up.
	// Unfortunately, the loops in these functions themselves cannot easily be parallel, but we can at least use sections;
	#pragma omp parallel sections
	{
		#pragma omp section
		SL_copyNewPolys(ss, CDDM_get_polys(output_dm));
		#pragma omp section
		SL_copyNewLoops(ss, CDDM_get_loops(output_dm));
		#pragma omp section
		SL_copyNewEdges(ss, CDDM_get_edges(output_dm));

		#pragma omp section
		{
			// ORIGININDEX vertex;
			input_index = ss->input->getVertDataArray(ss->input, CD_ORIGINDEX);
			output_index = CustomData_get_layer(&output_dm->vertData, CD_ORIGINDEX);
			// original vertices are at the beginning
			for (i = 0; i < ss->numVerts; i++) {
				output_index[i] = input_index ? input_index[i] : i;
			}
			// Then new verts, no original index
			for (; i < numVerts; i++) { // Note, keep i from previous loop.
				output_index[i] = ORIGINDEX_NONE;
			}
			// ORIGININDEX edge;
			input_index = ss->input->getEdgeDataArray(ss->input, CD_ORIGINDEX);
			output_index = CustomData_get_layer(&output_dm->edgeData, CD_ORIGINDEX);
			// First split edges;
			a = 0;
			for (i = 0; i < ss->numEdges; i++) {
				output_index[a++] = input_index ? input_index[i] : i;
				output_index[a++] = input_index ? input_index[i] : i;
			}
			// Then new internal edges
			for (; a < numEdges; a++) { // Note, keep i from previous loop.
				output_index[i] = ORIGINDEX_NONE;
			}
			// ORIGININDEX face;
			input_index = ss->input->getTessFaceDataArray(ss->input, CD_ORIGINDEX);
			output_index = CustomData_get_layer(&output_dm->faceData, CD_ORIGINDEX);
			a = 0;
			for (i = 0; i < ss->numFaces; i++) {
				MPoly *poly = &ss->mpoly[i];
				int internalFaces = poly->totloop == 3 ? 4 : poly->totloop;
				for (j = 0; j < internalFaces; j++, a++)
					output_index[a] = input_index ? input_index[i] : i;
			}
			// Same originindex for polys as or faces, probably.
			//CustomData_add_layer(&output->polyData, CD_ORIGINDEX, CD_REFERENCE, NULL, numFaces);
		}
	}
	return output_dm;
}


/////////////////////////////////////////////////////////////
// Misc helpers

/*
static int _sharedEdge(MeshElemMap *v2e1, MeshElemMap *v2e2) {
	int i, j, e;
	for (i = 0; i < v2e1->count; i++) {
		e = v2e1->indices[i];
		for (j = 0; j < v2e2->count; j++) {
			if ( e == v2e1->indices[j]) {
				return e;
			}
		}
	}
	return -1;
}
*/

static int _edgeIsBoundary(SLSubSurf *ss, int e) {
	return ss->edge2poly[e].count < 2;
}

static int _vertIsBoundary(SLSubSurf *ss, int v) {
	int i;
	MeshElemMap *map = &ss->vert2edge[v];
	for (i = 0; i < map->count; i++) {
		if (_edgeIsBoundary(ss, map->indices[i])) {
			return 1;
		}
	}
	return 0;
}

// Some values taken from "Combining 4- and 3-direction Subdivision"
// For n < 3, these values will never be used.
//float alpha[7]    = { 0.f, 0.25f,  0.3f, 3.0f/8.0f, 0.5f, 3.0f/5.0f, 5.0f/8.0f };
float alphabar[7] = { 0.f, -0.5f, -0.4f, -0.25f, 0.0f, 0.2f, 0.25f }; // Suggest in article;
//float alphabar[7]   = { 0.f, -0.5f, -0.4f,     -0.4f, 0.0f,  0.16505f,     0.25f }; // Modified (arbitrarly) by me

static double _alphabar(int edges) {
	return alphabar[edges > 6 ? 6 : edges];
}

/////////////////////////////////////////////////////////////
// Actual smoothing stuff
// The code is written w.r.t optimization with the assumption that smoothing will be used.
void SL_syncVerts(SLSubSurf *ss, DerivedMesh *output) {
	MVert *o_vert = DM_get_vert_data_layer(output, CD_MVERT);
	int i;
	// Face centroids aren't moved during smoothing (by any algorithm i've seen), but edge centroids do. 
	// When computing the smoothing part, we need to keep track of both, so we can't use o_vert->co directly.
	float (*eco)[3] = MEM_callocN(sizeof(float[3])*ss->numEdges, "sl edge centroids");

	// Compute centroid, used for smoothing and other things;
	//printf("Computing face centroids\n");
	// Perhaps i need to split the function for faceVertIdx to deal with this.
	#pragma omp parallel for
	for (i = 0; i < ss->numFaces; i++) {
		MPoly *poly = &ss->mpoly[i];
		//if (!face->requiresUpdate) continue; // TODO partial updates..
		//face->requiresUpdate = 0;
		int n = poly->totloop;
		
		if (n > 3) {
			int j;
			MLoop *loop = &ss->mloop[poly->loopstart];
			// Might as well use the final vert coordinates. Face centers aren't smoothed.
			float *coord = o_vert[ss->poly2vert[i]].co;

			zero_v3(coord);
			for (j = 0; j < n; j++) {
				add_v3_v3(coord, _origCoord(ss, loop[j].v));
			}
			mul_v3_fl(coord, 1.0f / n );
			//printf("o_vert[%d].no = {%d, %d, %d}\n", faceVertIdx, o_vert[ss->poly2vert[i]].no[0],o_vert[ss->poly2vert[i]].no[1],o_vert[ss->poly2vert[i]].no[2]);
		}
	}
	
	// also for edges;
	#pragma omp parallel for
	for (i = 0; i < ss->numEdges; i++) {
		float *c1, *c2;
		MEdge *edge = &ss->medge[i];
		//if (!edge->requiresUpdate) continue;
		c1 = _origCoord(ss, edge->v1);
		c2 = _origCoord(ss, edge->v2);
		mean_v3_v3v3(eco[i], c1, c2);
	}

	//printf("Computing vert smoothing (smooth %s) for %d verts\n", ss->smoothing ? "on":"off", ss->numVerts);
	#pragma omp parallel for
	// Loop over vertices and smooth out the Stam-Loop subsurface coordinate;
	for (i = 0; i < ss->numVerts; i++) {
		float *origCoord = _origCoord(ss, i);
		//if (!vert->requiresUpdate) continue;
		//vert->requiresUpdate = 0;
		float *coord = o_vert[i].co;
		MeshElemMap *v2p = &ss->vert2poly[i];
		
		if (!ss->smoothing) {
			copy_v3_v3(coord, origCoord);
		} else {
			// Compute average sharpness and seam;
			MeshElemMap *v2e = &ss->vert2edge[i];
			int j;
			//int seamCount = 0;
			int sharpnessCount = 0;
			float avgSharpness = 0.0f;
			//int seam;
			for (j = 0; j < v2e->count; j++) {
				int e = v2e->indices[j];
				MEdge *edge = &ss->medge[e];
				//MeshElemMap *e2p = &ss->edge2poly[e];

				// This is just for smoothing the UV coordinates (?)
				//if ( (edge->flag & ME_SEAM) && e2p->count < 2)
				//	seamCount++;
				if (edge->crease != 0) {
					sharpnessCount++;
					avgSharpness += ss->medge[e].crease / 255.f;
				}
			}

			if (sharpnessCount) {
				avgSharpness /= sharpnessCount;
				if ( avgSharpness > 1.0f ) {
					avgSharpness = 1.0f;
				}
			}

			// TODO: Is this correct? I don't know how to deal with seams, or why they matter here.
			//seam = seamCount >= 2 && seamCount != v2e->count;

			// Prepare to fill the output coordinate
			zero_v3(coord);

			// Now do the smoothing;
			if (_vertIsBoundary(ss, i)) {
				// No articles cover open edges, so I'm copying what CCG does.
				// smooth_coord = 1/2 * orig_coord + 1/2 * edge_centroid_average (only boundary edges)
				int avgCount = 0;
				zero_v3(coord);
				for (j = 0; j < v2e->count; j++) {
					int e = v2e->indices[j];
					if (_edgeIsBoundary(ss, e)) {
						add_v3_v3(coord, eco[e]);
						avgCount++;
					}
				}
				mul_v3_fl(coord, 0.5f / avgCount);
				madd_v3_v3fl(coord, origCoord, 0.5f);
			}
			else {
				// Tried to copy the article by Stam & Loop, but it didn't look right. Copying Catmull-Clark behavior (see wikipedia)
				// This entire part is up for debate. Feel free to modify until a nice smoothing is obtained.
				int avgCount, edgeMult;

				avgCount = 0;

				// Weights for edges are multiple of shared faces;
				for (j = 0; j < v2e->count; j++) {
					int k;
					int e = v2e->indices[j];
					MeshElemMap *e2p = &ss->edge2poly[e];
					// Lose edges, not covered normally by descriptions for smoothing. Just making something up here
					if (e2p->count == 0) {
						edgeMult = 1;
					} else {
						edgeMult = 0;
						for (k = 0; k < e2p->count; k++) {
							int p = e2p->indices[k];
							MPoly *poly = &ss->mpoly[p];
							if (poly->totloop == 3)
								edgeMult++;
							else
								edgeMult+=2;
						}
					}
					madd_v3_v3fl(coord, eco[e], edgeMult);
					avgCount += edgeMult;
				}
				mul_v3_fl(coord, (1.0f - _alphabar(v2p->count)) / avgCount);
				madd_v3_v3fl(coord, origCoord, _alphabar(v2p->count));
				//printf("vert, avgCount = %d, alphabar = %f, coords = {%f, %f, %f}, ls_coords = {%f, %f, %f}\n", avgCount, alphabar[vert->numEdges], 
				//	   vert->coords[0], vert->coords[1], vert->coords[2], vert->sl_coords[0], vert->sl_coords[1], vert->sl_coords[2]);

				/*
				f or (j = 0; j < v2p->count; j++) {
					int p = v2p->indices[j];
					poly = &ss->mpoly[p];
					if (poly->totloop > 3) {
						add_v3_v3(coord, ???? face->centroid); // Need to generate index to new vert to pull this off
						avgCount++; // Note that the subdivided area is a quad for any ngon > 3
					}
				}*/
			}

			// Deal with sharpness and seams
			// Code snipped converted from CCG (undocumented mystery code)
			if ((sharpnessCount > 1 && v2p->count > 0) /*|| seam*/) {
				// TODO: Haven't checked this carefully yet.
				int x;
				float q[3];

				/*
				if (seam) {
					avgSharpness = 1.0f;
					sharpnessCount = seamCount;
				}*/

				zero_v3(q);
				for (j = 0; j < v2e->count; j++) {
					int e = v2e->indices[j];
					//MeshElemMap *e2p = &ss->edge2poly[e];
					MEdge *edge = &ss->medge[e];
					/*if ( edge->flag & ME_SEAM ) {
						if (e2p->count < 2)
							add_v3_v3(q, eco[e]);
					}
					else*/ if (edge->crease != 0) {
						add_v3_v3(q, eco[e]);
					}
				}

				mul_v3_fl(q, 1.0f / sharpnessCount);

				if (sharpnessCount > 2 || v2p->count == 1 /*|| seam*/) {
					/* q = q + (co - q) * avgSharpness */
					for (x = 0; x < 3; x++) q[x] += (origCoord[x] - q[x])*avgSharpness;
				}

				/* r = co * 0.5 + q * 0.5 */
				//mean_v3_v3v3(q, origCoord, q);

				/* nCo = nCo + (r - nCo) * avgSharpness */
				for (x = 0; x < 3; x++) coord[x] += (0.5f * (q[x] + origCoord[x]) - coord[x]) * avgSharpness;
			}
		}
	}

	//printf("Computing edge smoothing\n");
	// Loop over edges and smooth
	#pragma omp parallel for
	for (i = 0; i < ss->numEdges; i++) {
		int vertidx = ss->numVerts + i;
		MEdge *edge = &ss->medge[i];
		//if (!edge->requiresUpdate) continue; // TODO
		//edge->requiresUpdate = 0;
		MeshElemMap *e2p = &ss->edge2poly[i];		
		float *coord = o_vert[vertidx].co;
		// Create the interpolated coordinates
		if (!ss->smoothing || e2p->count < 2 || edge->crease >= 255) { // If its an edge, or maximum sharpness, then just average.
			copy_v3_v3(coord, eco[i]); // Already there.
		} else { // Otherwise smooth
			int j, avgCount, numTris = 0, numQuads = 0, numNeighbors;
			// This entire part is up for debate. Feel free to modify until a nice smoothing is obtained.
			
			zero_v3(coord);
			avgCount = 0;

			e2p = &ss->edge2poly[i];
			for (j = 0; j < e2p->count; j++) {
				int p = e2p->indices[j];
				MPoly *poly = &ss->mpoly[p];
				if (poly->totloop == 3) {
					int k;
					// Triangles are split differently from the rest;
					// There are connections to the center nodes of the two opposite edges
					MLoop *loop = &ss->mloop[poly->loopstart];
					for (k = 0; k < poly->totloop; k++) {
						int temp_e = loop[k].e;
						if ( temp_e != i ) { // Then opposite edge
							madd_v3_v3fl(coord, eco[temp_e], 2);
						}
					}
					numTris++;
					avgCount += 4; // 2 edges each
				} else {
					// Otherwise all ngons are split into quads, leaving one center node and edges
					madd_v3_v3fl(coord, o_vert[ss->poly2vert[p]].co, 4);
					// Now find the other edges that share a node;
					/*for (j = 0; j < poly->totloop; j++) {
						int temp_e = loop[k].e;
						if ( temp_e != j ) {
							MEdge *tempEdge = &ss->medge[temp_e];
							// Check for a shared node;
							if (tempEdge->v1 == edge->v1 || tempEdge->v1 == edge->v2 ||
								tempEdge->v2 == edge->v1 || tempEdge->v2 == edge->v2) {
								add_v3_v3(coord, eco[temp_e]);
							}
						}
					}*/
					numQuads++;
					avgCount += 4; // 1x2 from centroid + 2x1 from edges
				}
			}
			madd_v3_v3fl(coord, eco[i], 2*(numTris+numQuads*2));
			avgCount += 2*(numTris+numQuads*2);

			numNeighbors = numTris*2 + numQuads + 2;
			mul_v3_fl(coord, (1.0f - _alphabar(numNeighbors))/avgCount);
			madd_v3_v3fl(coord, eco[i], _alphabar(numNeighbors));
			//printf("edge, numTris = %d, numQuads = %d, avgCount = %d, alphabar = %f, coords = {%f, %f, %f}, ls_coords = {%f, %f, %f}\n", numTris, numQuads, avgCount, alphabar[numNeighbors],
			//	   edge->centroid[0], edge->centroid[1], edge->centroid[2], edge->sl_coords[0], edge->sl_coords[1], edge->sl_coords[2]);

			// And take into account sharpness
			if (edge->crease > 0 ) {
				int x;
				float crease = edge->crease / 255.f;
				for (x = 0; x < 3; x++) {
					coord[x] += crease * (eco[i][x] - coord[x]);
				}
			}
			//printf("Smoothed coordinate[%d] (edge) = {%e, %e, %e}\n", vertidx, o_vert[vertidx].co[0], o_vert[vertidx].co[1], o_vert[vertidx].co[2]);
		}
	}

	// Loop over faces and smooth
	/*
	for (i = 0; i < ss->numFaces; i++) {
		poly = &ss->mpoly[i];
		loop = &ss->mloop[poly->loopstart];
	}*/
	
	MEM_freeN(eco);
}

// Syncs the UV coordinates for layer n
void SL_syncUV(SLSubSurf *ss, DerivedMesh *output, int UNUSED(smoothing), int n) {
	// TODO: Support smoothing subsurf uvs
	int i,j;
	// No smoothing for now.
	DerivedMesh *input = ss->input;
	MLoopUV *loopuv = CustomData_get_layer_n(&input->loopData, CD_MLOOPUV, n);
	MTexPoly *tpoly = CustomData_get_layer_n(&input->polyData, CD_MTEXPOLY, n);
	MLoopUV *o_loopuv = CustomData_get_layer_n(&output->loopData, CD_MLOOPUV, n);
	MTexPoly *o_tpoly = CustomData_get_layer_n(&output->polyData, CD_MTEXPOLY, n);
	
	for (i = 0; i < ss->numFaces; i++) {
		MPoly *poly = &ss->mpoly[i];
		MTexPoly *tp = &tpoly[i];
		MLoopUV *luv = &loopuv[poly->loopstart];
		int n = poly->totloop;
		
		if (n == 3) {
			// Compute the edge coordinates;
			float edge_uv[3][2];
			for (j = 0; j < 3; j++) {
				mean_v2_v2v2(edge_uv[j], luv[j].uv, luv[(j+1) % 3].uv);
			}

			for (j = 0; j < 3; j++) {
				int k = (j + 2) % 3;
				copy_v2_v2(o_loopuv->uv, luv[j].uv); o_loopuv++;
				copy_v2_v2(o_loopuv->uv, edge_uv[j]); o_loopuv++;
				copy_v2_v2(o_loopuv->uv, edge_uv[k]); o_loopuv++;
				ME_MTEXFACE_CPY(o_tpoly, tp); o_tpoly++;
			}
			copy_v2_v2(o_loopuv->uv, edge_uv[2]); o_loopuv++;
			copy_v2_v2(o_loopuv->uv, edge_uv[0]); o_loopuv++;
			copy_v2_v2(o_loopuv->uv, edge_uv[1]); o_loopuv++;
			ME_MTEXFACE_CPY(o_tpoly, tp); o_tpoly++;
		} else {
			// Compute the edge and center coordinates;
			float edge_uv[n][2];
			float center_uv[2];
			zero_v3(center_uv);
			for (j = 0; j < n; j++) {
				mean_v2_v2v2(edge_uv[j], luv[j].uv, luv[(j+1) % n].uv);
				// And the center coordinate;
				add_v2_v2(center_uv, luv[j].uv);
			}
			mul_v2_fl(center_uv, 1.0f/n);

			// Construct the new loops
			for (j = 0; j < n; j++) {
				int k = (j + n - 1) % n;
				copy_v2_v2(o_loopuv->uv, luv[j].uv); o_loopuv++;
				copy_v2_v2(o_loopuv->uv, edge_uv[j]); o_loopuv++;
				copy_v2_v2(o_loopuv->uv, center_uv); o_loopuv++;
				copy_v2_v2(o_loopuv->uv, edge_uv[k]); o_loopuv++;
				ME_MTEXFACE_CPY(o_tpoly, tp); o_tpoly++;
			}
		}
	}
}

void SL_syncPaint(SLSubSurf *ss, DerivedMesh *output, int n) {
	// No documentation of mcol, and I can't seem to get it to do *anything*. Giving up on it.
	//MCol *input_paint = CustomData_get_layer_n(&input->vertData, CD_MCOL, n);
	MLoopCol *loopc = CustomData_get_layer_n(&ss->input->loopData, CD_MLOOPCOL, n);
	int numPolys, i, j;
	MPoly *poly;
	MLoopCol *loop;
	MLoopCol *o_loopc = CustomData_get_layer_n(&output->loopData, CD_MLOOPCOL, n);
	if (!o_loopc) {
		// Not sure how this works with "n" ? 
		o_loopc = CustomData_add_layer(&output->loopData, CD_MLOOPCOL, CD_CALLOC, NULL, output->numVertData);
	}
	
	numPolys = ss->numFaces;
	for (i = 0; i < numPolys; i++) {
		poly = &ss->mpoly[i];
		loop = &loopc[poly->loopstart];
		if (poly->totloop == 3) {
			// Start by computing the edge colors
			char a[3], r[3], g[3], b[3];
			for (j = 0; j < 3; j++) {
				int k = (j + 1) % 3;
				a[j] = (char)(((float)loop[k].a + (float)loop[j].a)*0.5f);
				r[j] = (char)(((float)loop[k].r + (float)loop[j].r)*0.5f);
				g[j] = (char)(((float)loop[k].g + (float)loop[j].g)*0.5f);
				b[j] = (char)(((float)loop[k].b + (float)loop[j].b)*0.5f);
			}
			
			// Outer triangles first;
			for (j = 0; j < 3; j++) {
				int k = (j + 2) % 3;
				o_loopc->a = loop[j].a; o_loopc->r = loop[j].r; o_loopc->g = loop[j].g; o_loopc->b = loop[j].b; o_loopc++;
				o_loopc->a = a[j];      o_loopc->r = r[j];      o_loopc->g = g[j];      o_loopc->b = b[j]; o_loopc++;
				o_loopc->a = a[k];      o_loopc->r = r[k];      o_loopc->g = g[k];      o_loopc->b = b[k]; o_loopc++;
			}
			// Center triangle
			o_loopc->a = a[2]; o_loopc->r = r[2]; o_loopc->g = g[2]; o_loopc->b = b[2]; o_loopc++;
			o_loopc->a = a[0]; o_loopc->r = r[0]; o_loopc->g = g[0]; o_loopc->b = b[0]; o_loopc++;
			o_loopc->a = a[1]; o_loopc->r = r[1]; o_loopc->g = g[1]; o_loopc->b = b[1]; o_loopc++;
		} else {
			int n = poly->totloop;
			// Center color;
			char ca, cr, cg, cb;
			float fca = 0, fcr = 0, fcg = 0, fcb = 0;
			// Edge color;
			char a[n], r[n], g[n], b[n];
			for (j = 0; j < n; j++) {
				int k = (j + 1) % n;
				fca += loop[j].a; fcr += loop[j].r; fcg += loop[j].g; fcb += loop[j].b;
				a[j] = (char)(((float)loop[k].a + (float)loop[j].a)*0.5f);
				r[j] = (char)(((float)loop[k].r + (float)loop[j].r)*0.5f);
				g[j] = (char)(((float)loop[k].g + (float)loop[j].g)*0.5f);
				b[j] = (char)(((float)loop[k].b + (float)loop[j].b)*0.5f);
			}
			ca = (char)(fca/n); cr = (char)(fcr/n); cg = (char)(fcg/n); cb = (char)(fcb/n);
			//printf("center color = %d, %d, %d, %d\n", ca, cr, cg, cb);
			// Now construct the loops, which all start from the corner node
			for (j = 0; j < n; j++) {
				int k = (j + n - 1) % n;
				o_loopc->a = loop[j].a; o_loopc->r = loop[j].r; o_loopc->g = loop[j].g; o_loopc->b = loop[j].b; o_loopc++;
				o_loopc->a = a[j];      o_loopc->r = r[j];      o_loopc->g = g[j];      o_loopc->b = b[j];      o_loopc++;
				o_loopc->a = ca;        o_loopc->r = cr;        o_loopc->g = cg;        o_loopc->b = cb;        o_loopc++;
				o_loopc->a = a[k];      o_loopc->r = r[k];      o_loopc->g = g[k];      o_loopc->b = b[k];      o_loopc++;
			}
		}
	}
}

/////////////////////////////////////////////////////////////
