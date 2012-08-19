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

/////////////////////////////////////////////////////////////
// Support functions for faces;

int SL_giveNumberOfInternalFaces(MPoly *poly) {
	// Undecided on how anything above quads should be subdivided. Default Catmull-Clark behavior seems reasonable.
	return (poly->totloop == 3) ? 4 : poly->totloop;
}

int SL_giveNumberOfInternalNodes(MPoly *poly) {
	// No center node for triangles.
	return (poly->totloop == 3) ? 0 : 1;
}

int SL_giveNumberOfInternalEdges(MPoly *poly) {
	return poly->totloop;
}

int SL_giveNumberOfInternalLoops(MPoly *poly) {
	// We could actually save some memory on loops, since its predictable and we know how they will overlap.
	// (note, this isn't done though, since the code is simpler is you ignore that)
	return poly->totloop*4; // (just happens to be correct for triangles as well)
}

/////////////////////////////////////////////////////////////
// External helpers

int SL_giveTotalNumberOfSubVerts(SLSubSurf *ss) {
	// Simple things first, corners and interpolated edge verts;
	int i, totNodes = ss->numVerts + ss->numEdges;
	// Then faces, which varies;
	for (i = 0; i < ss->numFaces; i++) {
		totNodes += SL_giveNumberOfInternalNodes(&ss->mpoly[i]);
	}
	return totNodes;
}

int SL_giveTotalNumberOfSubEdges(SLSubSurf *ss) {
	int i, totEdges = ss->numEdges * 2;
	for (i = 0; i < ss->numFaces; i++) {
		totEdges += SL_giveNumberOfInternalEdges(&ss->mpoly[i]);
	}
	return totEdges;
}

int SL_giveTotalNumberOfSubFaces(SLSubSurf *ss) {
	int i, totFaces = 0;
	for (i = 0; i < ss->numFaces; i++) {
		totFaces += SL_giveNumberOfInternalFaces(&ss->mpoly[i]);
	}
	return totFaces;
}

int SL_giveTotalNumberOfSubLoops(SLSubSurf *ss) {
	int i, totLoops = 0;
	for (i = 0; i < ss->numFaces; i++) {
		totLoops += SL_giveNumberOfInternalLoops(&ss->mpoly[i]);
	}
	return totLoops;
}

/////////////////////////////////////////////////////////////

void SL_copyNewLoops(SLSubSurf *ss, MLoop *mloops) {
	int i, j, k, prevJ, subEdgeNext, subEdgePrev, faceEdgeStartIdx, faceVertStartIdx;
	MPoly *poly;
	MLoop *loop;
	MEdge *eNext, *ePrev;

	k = 0;
	faceEdgeStartIdx = ss->numEdges*2; // Number the new internal edges as we go
	faceVertStartIdx = ss->numVerts + ss->numEdges; // Number the new internal verts as we go
	for (i = 0; i < ss->numFaces; i++) {
		poly = &ss->mpoly[i];
		loop = &ss->mloop[poly->loopstart];
		if (poly->totloop == 3) {
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
			for (j = 0; j < poly->totloop; j++) { // Loop over each sub-quad;
				// Unsure if i negative values work, adding numVerts just in case.
				prevJ = (j - 1 + poly->totloop) % poly->totloop;
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
			faceEdgeStartIdx += poly->totloop;
			faceVertStartIdx += 1;
		}
	}
}

void SL_copyNewPolys(SLSubSurf *ss, MPoly *mpolys) {
	int i, j, k = 0, m = 0;
	for (i = 0; i < ss->numFaces; i++) {
		MPoly *poly = &ss->mpoly[i];
		if (poly->totloop == 3) {
			for (j = 0; j < 4; j++) { // note! 4 faces
				mpolys[m].loopstart = k;
				mpolys[m].totloop = 3;
				mpolys[m].mat_nr = poly->mat_nr;
				mpolys[m].flag = poly->flag;
				m += 1;
				k += 3;
			}
		} else {
			for (j = 0; j < poly->totloop; j++) {
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
	int i, j, prevJ, k, faceVertStartIdx;
	MPoly *poly;
	MLoop *loop;

	k = 0;
	faceVertStartIdx = ss->numVerts + ss->numEdges; // Number the new internal verts as we go
	for (i = 0; i < ss->numFaces; i++) {
		poly = &ss->mpoly[i];
		loop = &ss->mloop[poly->loopstart];
		if (poly->totloop == 3) {
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
			for (j = 0; j < poly->totloop; j++) { // Loop over each sub-quad;
				// Unsure if i negative values work, adding numVerts just in case.
				prevJ = (j - 1 + poly->totloop) % poly->totloop;

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
	int i, j, k, faceVertStartIdx;
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
	faceVertStartIdx = ss->numVerts + ss->numEdges; // Number the new internal verts as we go
	for (i = 0; i < ss->numFaces; i++) {
		MPoly *poly = &ss->mpoly[i];
		MLoop *loop = &ss->mloop[poly->loopstart];
		if (poly->totloop == 3) {
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
			for (j = 0; j < poly->totloop; j++) {
				medges[k].v1 = ss->numVerts + loop[j].e;
				medges[k].v2 = faceVertStartIdx;
				//printf("medges[%d] = {%d, %d}\n", i, medges[k].v1, medges[k].v2);
				k++;
			}
			faceVertStartIdx++;
		}
	}
}


static float  *_origCoord(SLSubSurf *ss, int vert) {
	return ss->vertexCos ? ss->vertexCos[vert] : ss->mvert[vert].co;
}

/////////////////////////////////////////////////////////////

SLSubSurf* SL_SubSurf_new(int smoothing, DerivedMesh *input, float (*vertexCos)[3]) {
	int numLoops;
	SLSubSurf *ss = (SLSubSurf*)MEM_callocN(sizeof(SLSubSurf), "slsubsurf");
	ss->smoothing = smoothing;
	ss->numVerts = ss->numEdges = ss->numFaces = 0;
	
	ss->input = input;
	ss->vertexCos = vertexCos;
	ss->mvert = CustomData_get_layer(&input->vertData, CD_MVERT);
	ss->medge = CustomData_get_layer(&input->edgeData, CD_MEDGE);
	ss->mpoly = CustomData_get_layer(&input->polyData, CD_MPOLY);
	ss->mloop = CustomData_get_layer(&input->loopData, CD_MLOOP);
	
	ss->numVerts = input->getNumVerts(input);
	ss->numEdges = input->getNumEdges(input);
	ss->numFaces = input->getNumPolys(input);
	numLoops = input->getNumLoops(input);
	
	// Necessary backwards mapping (could potentially be done without these, but that would make partial updates impossible (not that partial updates are supported yet))
	ss->poly2vert = MEM_callocN(sizeof(int)*ss->numFaces, "sl poly2vert");
	create_vert_poly_map(&ss->vert2poly, &ss->vert2poly_mem,
		ss->mpoly, ss->mloop, ss->numVerts, ss->numFaces, numLoops);
	create_vert_edge_map(&ss->vert2edge, &ss->vert2edge_mem,
		ss->medge, ss->numVerts, ss->numEdges);
	create_edge_poly_map(&ss->edge2poly, &ss->edge2poly_mem,
		ss->mpoly, ss->mloop, ss->numEdges, ss->numFaces, numLoops);

	ss->o_vert = NULL;
	ss->o_edge = NULL;
	ss->o_poly = NULL;
	ss->o_loop = NULL;

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
	MEM_freeN(ss->poly2vert);
	MEM_freeN(ss->vert2poly);
	MEM_freeN(ss->vert2poly_mem);
	MEM_freeN(ss->vert2edge);
	MEM_freeN(ss->vert2edge_mem);
	MEM_freeN(ss->edge2poly);
	MEM_freeN(ss->edge2poly_mem);
	if (ss->eco) {
		MEM_freeN(ss->eco);
	}
	MEM_freeN(ss);
}

// A few fucntions that should be replaced from the CDDM default functions;
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

static void slDM_recalcTessellation(DerivedMesh *UNUSED(dm)) { /* Nothing to do */ }


// Not sure what the bounding box should account for, both original and refined mesh makes some sense.
// This can be done faster than the default implementation in CDDerivedMesh, but its just 1 subdivision level
/*
static void _minmax_v3_v3v3(const float vec[3], float min[3], float max[3]) {
	if (min[0] > vec[0]) min[0] = vec[0];
	if (min[1] > vec[1]) min[1] = vec[1];
	if (min[2] > vec[2]) min[2] = vec[2];
	if (max[0] < vec[0]) max[0] = vec[0];
	if (max[1] < vec[1]) max[1] = vec[1];
	if (max[2] < vec[2]) max[2] = vec[2];
}
static void slDM_getMinMax(DerivedMesh *dm, float min_r[3], float max_r[3]) {
	SLDerivedMesh *sldm = (SLDerivedMesh*)dm;
	SLSubSurf *ss = sldm->ss;
	int i;
	if (ss->numVerts == 0) {
		zero_v3(min_r);
		zero_v3(max_r);
		return;
	}
	for(i = 0; i < ss->numVerts; i++) { // TODO: Should i use the smoothed coordinates(?) (does anyone care?)
		_minmax_v3_v3v3(_origCoord(ss, i), min_r, max_r);
	}
}
*/

void SL_syncPaint(SLSubSurf *ss, DerivedMesh *output) {
	// No documentation of mcol, and I can't seem to get it to do *anything*. Giving up on it.
	//MCol *input_paint = CustomData_get_layer(&input->vertData, CD_MCOL);
	MLoopCol *loopc = CustomData_get_layer(&ss->input->loopData, CD_MLOOPCOL);
	if (loopc) {
		int numPolys, i, j;
		MPoly *poly;
		MLoopCol *loop;
		MLoopCol *o_loopc = CustomData_get_layer(&output->loopData, CD_MLOOPCOL);
		if (!o_loopc) {
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
				printf("center color = %d, %d, %d, %d\n", ca, cr, cg, cb);
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
		//printf("Vertex paint synced\n");
	}
}

DerivedMesh *SL_SubSurf_constructOutput(SLSubSurf *ss) {
	int numVerts, numEdges, numFaces, numLoops;
	int i, j, a, *polyidx, *input_index, *output_index;
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
	ss->input->getTessFaceDataArray(ss->input, CD_ORIGINDEX);
	
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
	if (!CustomData_get_layer(&output_dm->faceData, CD_POLYINDEX))
		CustomData_add_layer(&output_dm->faceData, CD_POLYINDEX, CD_CALLOC, NULL, numFaces);
	
	output->cddm.mvert = CustomData_get_layer(&output_dm->vertData, CD_MVERT);
	output->cddm.medge = CustomData_get_layer(&output_dm->edgeData, CD_MEDGE);
	output->cddm.mface = CustomData_get_layer(&output_dm->faceData, CD_MFACE);
	output->cddm.mloop = CustomData_get_layer(&output_dm->loopData, CD_MLOOP);
	output->cddm.mpoly = CustomData_get_layer(&output_dm->polyData, CD_MPOLY);
	// End of CDDM_from_template copy	
	
	// Replace some functions from CDDM
	output_dm->recalcTessellation = slDM_recalcTessellation;
	output_dm->release = slDM_release;
	// TODO: Replace any other functions from default CDDM, perhaps?
	
	// This was done in in the CCG code, so i suppose I'll use that as well.
	CustomData_free_layer_active(&output_dm->polyData, CD_NORMAL, output_dm->numPolyData);
	
	// Convenience. Directly access these.
	ss->o_vert = CDDM_get_verts(output_dm);
	ss->o_edge = CDDM_get_edges(output_dm);
	ss->o_loop = CDDM_get_loops(output_dm);
	ss->o_poly = CDDM_get_polys(output_dm);

	// Face centroids aren't moved during smoothing (by any algorithm i've seen), but edge centroids do. 
	// When computing the smoothing part, we need to keep track of both, so we can't use o_vert->co directly.
	ss->eco = MEM_callocN(sizeof(float[3])*numEdges, "sl edge centroids");
	
	// We can directly fill the structural part, which doesn't actually need the new vertex coordinates.
	SL_copyNewPolys(ss, ss->o_poly);
	SL_copyNewLoops(ss, ss->o_loop);
	SL_copyNewEdges(ss, ss->o_edge);
	SL_copyNewTessFaces(ss, CDDM_get_tessfaces(output_dm)); // We don't need the faces for anything.
	// TODO: Will probably need some other info copied over from here 
	
	// Other data we must fill;
	/* We absolutely need that layer, else it's no valid tessellated data! */
	polyidx = CustomData_get_layer(&output_dm->faceData, CD_POLYINDEX);
	for (i = 0; i < numFaces; i++) polyidx[i] = i; // TODO: Check this, i think its OK. (it should just be 1 face per poly)

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
		for (j = 0; j < SL_giveNumberOfInternalFaces(poly); j++, a++)
			output_index[a] = input_index ? input_index[i] : i;
	}
	// TODO: Not sure why this isn't needed (?)
	//CustomData_add_layer(&output->polyData, CD_ORIGINDEX, CD_CALLOC, NULL, numPolys);
	//CustomData_get_layer(&output->polyData, CD_ORIGINDEX);
	
	// Vertex paint
	SL_syncPaint(ss, output_dm);
	
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
void SL_processSync(SLSubSurf *ss) {
	MPoly *poly;
	MLoop *loop;
	MEdge *edge;
	MeshElemMap *v2e, *v2p, *e2p;
	float *coord;
	float avgSharpness;
	int i, j, k, seamCount, sharpnessCount, seam, faceVertIdx;

	// Compute centroid, used for smoothing and other things;
	printf("Computing face centroids\n");
	faceVertIdx = ss->numVerts + ss->numEdges;
	for (i = 0; i < ss->numFaces; i++) {
		poly = &ss->mpoly[i];
		//if (!face->requiresUpdate) continue; // TODO
		//face->requiresUpdate = 0;
		
		if (poly->totloop == 3) {
			// Nothing to do for triangles.
			ss->poly2vert[i] = -1;
		} else {
			loop = &ss->mloop[poly->loopstart];
			// Might as well use the final vert coordinates. Face centers aren't smoothed.
			coord = ss->o_vert[faceVertIdx].co;

			zero_v3(coord);
			for (j = 0; j < poly->totloop; j++) {
				add_v3_v3(coord, _origCoord(ss, loop[j].v));
			}
			mul_v3_fl(coord, 1.0f / poly->totloop );
			ss->poly2vert[i] = faceVertIdx;
			//printf("o_vert[%d].no = {%d, %d, %d}\n", faceVertIdx, ss->o_vert[faceVertIdx].no[0],ss->o_vert[faceVertIdx].no[1],ss->o_vert[faceVertIdx].no[2]);
			faceVertIdx++;
		}
	}
	
	// also for edges;
	for (i = 0; i < ss->numEdges; i++) {
		int x;
		float *c1, *c2;
		edge = &ss->medge[i];
		//if (!edge->requiresUpdate) continue;
		coord = ss->eco[i]; // Directly on to the final vertex
		c1 = _origCoord(ss, edge->v1);
		c2 = _origCoord(ss, edge->v2);
		for (x = 0; x < 3; x++)
			coord[x] = 0.5*c1[x] + 0.5*c2[x];
	}

	printf("Computing vert smoothing (smooth %s) for %d verts\n", ss->smoothing ? "on":"off", ss->numVerts);
	// Loop over vertices and smooth out the Stam-Loop subsurface coordinate;
	for (i = 0; i < ss->numVerts; i++) {
		float *origCoord = _origCoord(ss, i);
		//if (!vert->requiresUpdate) continue;
		//vert->requiresUpdate = 0;
		coord = ss->o_vert[i].co;
		v2p = &ss->vert2poly[i];
		
		if (!ss->smoothing) {
			copy_v3_v3(coord, origCoord);
		} else {
			// Compute average sharpness and seam;
			v2e = &ss->vert2edge[i];
			seamCount = 0;
			sharpnessCount = 0;
			avgSharpness = 0.0f;
			for (j = 0; j < v2e->count; j++) {
				int e = v2e->indices[j];
				edge = &ss->medge[e];
				e2p = &ss->edge2poly[e];

				if ( (edge->flag & ME_SEAM) && e2p->count < 2)
					seamCount++;

				if (ss->medge[e].crease != 0) {
					sharpnessCount++;
					avgSharpness += ss->medge[e].crease / 256.f;
				}
			}

			if (sharpnessCount) {
				avgSharpness /= sharpnessCount;
				if ( avgSharpness > 1.0f ) {
					avgSharpness = 1.0f;
				}
			}

			if (seamCount < 2 || seamCount != v2e->count)
				seam = 0;

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
						add_v3_v3(coord, ss->eco[e]);
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
					int e = v2e->indices[j];
					e2p = &ss->edge2poly[e];
					// Lose edges, not covered normally by descriptions for smoothing. Just making something up here
					edgeMult = e2p->count == 0 ? 1 : 0;
					for (k = 0; k < e2p->count; k++) {
						int p = e2p->indices[k];
						poly = &ss->mpoly[p];
						if (poly->totloop == 3)
							edgeMult++;
						else
							edgeMult+=2;
					}
					madd_v3_v3fl(coord, ss->eco[e], edgeMult);
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
				// Leftovers from various tests;
				//mul_v3_fl(vert->sl_coords, 1.0f - );
				// Original position weighted in
				//mul_v3_fl(vert->sl_coords, 0.5f / avgCount );
				//madd_v3_v3fl(vert->sl_coords, vert->coords, 0.5f);
				//madd_v3_v3fl(vert->sl_coords, vert->coords, avgCount/3.0f-3.0f);
				//mul_v3_fl(vert->sl_coords, 3.0f / avgCount);

				// As stated in the the article by Stam, vertex positions need to be corrected,
				// or some thing else needs to be done. This doesn't look very nice.
				// Unfortunately, the article doesn't cover all cases.
			}

			// Deal with sharpness and seams
			// Code snipped converted from CCG (undocumented mystery code)
			if ((sharpnessCount > 1 && v2p->count > 0) || seam) {
				// TODO: Haven't checked this carefully yet.
				int x;
				float q[3];

				if (seam) {
					avgSharpness = 1.0f;
					sharpnessCount = seamCount;
				}

				zero_v3(q);
				for (j = 0; j < v2e->count; j++) {
					int e = v2e->indices[j];
					e2p = &ss->edge2poly[e];
					edge = &ss->medge[e];
					if ( edge->flag & ME_SEAM ) {
						if (e2p->count < 2)
							add_v3_v3(q, ss->eco[e]);
					}
					else if (edge->crease != 0) {
						add_v3_v3(q, ss->eco[e]);
					}
				}

				mul_v3_fl(q, 1.0f / sharpnessCount);

				if (sharpnessCount != 2 || seam) {
					/* q = q + (co - q) * avgSharpness */
					for (x = 0; x < 3; x++) q[x] += (origCoord[x] - q[x])*avgSharpness;
				}

				/* r = co * 0.75 + q * 0.25 */
				for (x = 0; x < 3; x++) q[x] = origCoord[x]*0.5f + q[x]*0.5f;

				/* nCo = nCo + (r - nCo) * avgSharpness */
				for (x = 0; x < 3; x++) coord[x] += (q[x] - coord[x]) * avgSharpness;
			}
		}
	}

	printf("Computing edge smoothing\n");
	// Loop over edges and smooth
	for (i = 0; i < ss->numEdges; i++) {
		int vertidx = ss->numVerts + i;
		edge = &ss->medge[i];
		//if (!edge->requiresUpdate) continue; // TODO
		//edge->requiresUpdate = 0;
		e2p = &ss->edge2poly[i];		
		coord = ss->o_vert[vertidx].co;
		// Create the interpolated coordinates
		if (!ss->smoothing || e2p->count < 2 || edge->crease >= 256) { // If its an edge, or maximum sharpness, then just average.
			copy_v3_v3(coord, ss->eco[i]); // Already there.
		} else { // Otherwise smooth
			int avgCount, numTris = 0, numQuads = 0, numNeighbors;
			// This entire part is up for debate. Feel free to modify until a nice smoothing is obtained.
			
			// Copy over the centroid to a temporary variable.
			zero_v3(coord);
			avgCount = 0;

			e2p = &ss->edge2poly[i];
			for (j = 0; j < e2p->count; j++) {
				int p = e2p->indices[j];
				poly = &ss->mpoly[p];
				if (poly->totloop == 3) {
					// Triangles are split differently from the rest;
					// There are connections to the center nodes of the two opposite edges
					loop = &ss->mloop[poly->loopstart];
					for (k = 0; k < poly->totloop; k++) {
						int temp_e = loop[k].e;
						if ( temp_e != i ) { // Then opposite edge
							madd_v3_v3fl(coord, ss->eco[temp_e], 2);
						}
					}
					numTris++;
					avgCount += 4; // 2 edges each
				} else {
					// Otherwise all ngons are split into quads, leaving one center node and edges
					madd_v3_v3fl(coord, ss->o_vert[ss->poly2vert[p]].co, 4);
					// Now find the other edges that share a node;
					/*for (j = 0; j < poly->totloop; j++) {
						int temp_e = loop[k].e;
						if ( temp_e != j ) {
							MEdge *tempEdge = &ss->medge[temp_e];
							// Check for a shared node;
							if (tempEdge->v1 == edge->v1 || tempEdge->v1 == edge->v2 ||
								tempEdge->v2 == edge->v1 || tempEdge->v2 == edge->v2) {
								add_v3_v3(coord, ss->eco[temp_e]);
							}
						}
					}*/
					numQuads++;
					avgCount += 4; // 1x2 from centroid + 2x1 from edges
				}
			}
			madd_v3_v3fl(coord, ss->eco[i], 2*(numTris+numQuads*2));
			avgCount += 2*(numTris+numQuads*2);

			numNeighbors = numTris*2 + numQuads + 2;
			mul_v3_fl(coord, (1.0f - _alphabar(numNeighbors))/avgCount);
			madd_v3_v3fl(coord, ss->eco[i], _alphabar(numNeighbors));
			//printf("edge, numTris = %d, numQuads = %d, avgCount = %d, alphabar = %f, coords = {%f, %f, %f}, ls_coords = {%f, %f, %f}\n", numTris, numQuads, avgCount, alphabar[numNeighbors],
			//	   edge->centroid[0], edge->centroid[1], edge->centroid[2], edge->sl_coords[0], edge->sl_coords[1], edge->sl_coords[2]);
			
			//printf("edge avgCount = %d\n", avgCount);
			//printf("edge sl_coords (accumulated) = %e, %e, %e\n", edge->sl_coords[0], edge->sl_coords[1], edge->sl_coords[2]);
			//mul_v3_fl(edge->sl_coords, 3.0f / avgCount);

			// Trying this for now (should be equivalent to CCG if its all quads) (compare with wikipedia, avgCount = 3*n )
			//madd_v3_v3fl(edge->sl_coords, edge->centroid, avgCount/3.0f-3.0f);
			//mul_v3_fl(edge->sl_coords, 3.0f / avgCount);
			//mul_v3_fl(edge->sl_coords, 4 + 2*edge->numFaces);
			//avgCount = 4 + 2*edge->numFaces;

			// And take into account sharpness
			if (edge->crease > 0 ) {
				int x;
				for (x = 0; x < 3; x++) {
					coord[x] += (edge->crease / 256.f) * (ss->eco[i][x] - coord[x]);
				}
			}
			//printf("Smoothed coordinate[%d] (edge) = {%e, %e, %e}\n", vertidx, ss->o_vert[vertidx].co[0], ss->o_vert[vertidx].co[1], ss->o_vert[vertidx].co[2]);
		}
	}

	// Loop over faces and smooth
	/*
	for (i = 0; i < ss->numFaces; i++) {
		poly = &ss->mpoly[i];
		loop = &ss->mloop[poly->loopstart];
	}*/
}

/////////////////////////////////////////////////////////////
