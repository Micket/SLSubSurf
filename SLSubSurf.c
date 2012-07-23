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

// Convenient macro for looping through linked lists (which is done a lot)
#define FOR_LIST(it, list) for (it = list; it != NULL; it = it->next)
// and the same for hashmaps
#define FOR_HASH(it, list) for (BLI_ghashIterator_init(it,list); !BLI_ghashIterator_isDone(it); BLI_ghashIterator_step(it))

/////////////////////////////////////////////////////////////
// Support functions for faces;

int SL_giveNumberOfInternalFaces(SLFace *face) {
	// Undecided on how anything above quads should be subdivided. Default Catmull-Clark behavior seems reasonable.
	return (face->numVerts == 3) ? 4 : face->numVerts;
}

int SL_giveNumberOfInternalNodes(SLFace *face) {
	// No center node for triangles.
	return (face->numVerts == 3) ? 0 : 1;
}

int SL_giveNumberOfInternalEdges(SLFace *face) {
	return face->numVerts;
}

int SL_giveNumberOfInternalLoops(SLFace *face) {
	return (face->numVerts == 3) ? 12 : face->numVerts*4;
}

/////////////////////////////////////////////////////////////
// External helpers

int SL_giveTotalNumberOfSubVerts(SLSubSurf *ss) {
	// Simple things first, corners and interpolated edge verts;
	int totNodes = ss->numVerts + ss->numEdges; // One new node per edge
	// Then faces, which varies;
	FOR_HASH(ss->it, ss->faces) {
		totNodes += SL_giveNumberOfInternalNodes((SLFace*)BLI_ghashIterator_getValue(ss->it));
	}
	return totNodes;
}

int SL_giveTotalNumberOfSubEdges(SLSubSurf *ss) {
	int totEdges = ss->numEdges * 2;
	FOR_HASH(ss->it, ss->faces) {
		totEdges += ((SLFace*)BLI_ghashIterator_getValue(ss->it))->numVerts; // (Holds for both triangles and ngons)
	}
	return totEdges;
}

int SL_giveTotalNumberOfSubFaces(SLSubSurf *ss) {
	// Since we know that all subdivided elements have the same number of sub-faces;
	int totFaces = 0; // One new node per edge
	// Then faces, which varies;
	FOR_HASH(ss->it, ss->faces) {
		totFaces += SL_giveNumberOfInternalFaces((SLFace*)BLI_ghashIterator_getValue(ss->it));
	}
	return totFaces;
}

int SL_giveTotalNumberOfSubLoops(SLSubSurf *ss) {
	// Since we know that all subdivided elements have the same number of sub-faces;
	int totLoops = 0; // One new node per edge
	// Then faces, which varies;
	FOR_HASH(ss->it, ss->faces) {
		totLoops += SL_giveNumberOfInternalLoops((SLFace*)BLI_ghashIterator_getValue(ss->it));
	}
	return totLoops;
}

/////////////////////////////////////////////////////////////

// Sets all (new) indices necessary for copying stuff back into DerivedMesh
void SL_renumberAll(SLSubSurf *ss)
{
	int idxA, idxB = 0;

	FOR_HASH(ss->it, ss->verts) {
		SLVert *vert = (SLVert*)BLI_ghashIterator_getValue(ss->it);
		vert->newVertIdx = idxA++;
	}
	FOR_HASH(ss->it, ss->edges) {
		SLEdge *edge = (SLEdge*)BLI_ghashIterator_getValue(ss->it);
		edge->newMetaIdx = idxB++;
	}
	idxA += idxB;
	idxB *= 2;
	FOR_HASH(ss->it, ss->faces) {
		SLFace *face = (SLFace*)BLI_ghashIterator_getValue(ss->it);
		if (face->numVerts > 3)
			face->newVertIdx = idxA++;
		face->newEdgeStartIdx = idxB;
		idxB += face->numVerts;
	}
}

void SL_copyNewPolys(SLSubSurf *ss, DMFlagMat *faceFlags, MPoly *mpolys)
{
	int i = 0, j, k;
	FOR_HASH(ss->it, ss->faces) {
		SLFace *face = (SLFace*)BLI_ghashIterator_getValue(ss->it);
		int flag = (faceFlags) ? faceFlags[i].flag : ME_SMOOTH;
		int mat_nr = (faceFlags) ? faceFlags[i].mat_nr : 0;
		if (face->numVerts == 3) {
			for (j = 0; j < 3; j++) {
				mpolys[i].loopstart = k;
				mpolys[i].totloop = 3;
				mpolys[i].mat_nr = mat_nr;
				mpolys[i].flag = flag;
				i += 1;
				k += 3;
			}
		} else {
			for (j = 0; j < face->numVerts; j++) {
				mpolys[i].loopstart = k;
				mpolys[i].totloop = 4;
				mpolys[i].mat_nr = mat_nr;
				mpolys[i].flag = flag;
				i += 1;
				k += 4;
			}
		}
	}
}

void SL_copyNewLoops(SLSubSurf *ss, MLoop *mloops)
{
	int i = 0, j, prevJ, subEdgeNext, subEdgePrev;
	SLFace *face;
	SLVert *vert;
	SLEdge *eNext, *ePrev;

	FOR_HASH(ss->it, ss->faces) {
		face = (SLFace*)BLI_ghashIterator_getValue(ss->it);
		if (face->numVerts == 3) {
			// First the corner triangles
			for (j = 0; j < 3; j++) {
				prevJ = (j + 2) % 3;
				eNext = face->edges[j];
				ePrev = face->edges[prevJ];
				vert = face->verts[j];

				// Check ordering of edge (to determine which sub-edge should be used)
				subEdgeNext = eNext->v0 != vert;
				subEdgePrev = ePrev->v0 != vert;

				// Corner to next edge;
				mloops[i+0].v = vert->newVertIdx;
				mloops[i+0].e = eNext->newMetaIdx*2 + subEdgeNext;
				// next edge node to internal edge j;
				mloops[i+1].v = eNext->newMetaIdx + ss->numVerts;
				mloops[i+1].e = face->newEdgeStartIdx + j;
				// previous edge node to previous edge;
				mloops[i+2].v = ePrev->newMetaIdx + ss->numVerts;
				mloops[i+2].e = ePrev->newMetaIdx*2 + subEdgePrev;
				i += 3;
			}

			// Last poly is the center polygon, only internal edges and edge nodes
			mloops[i+0].v = face->edges[2]->newMetaIdx + ss->numVerts;
			mloops[i+0].e = face->newEdgeStartIdx;
			mloops[i+1].v = face->edges[0]->newMetaIdx + ss->numVerts;
			mloops[i+1].e = face->newEdgeStartIdx+1;
			mloops[i+2].v = face->edges[1]->newMetaIdx + ss->numVerts;
			mloops[i+2].e = face->newEdgeStartIdx+2;
			i += 3;
		} else {
			for (j = 0; j < face->numVerts; j++) { // Loop over each sub-quad;
				// Unsure if i negative values work, adding numVerts just in case.
				prevJ = (j - 1 + face->numVerts) % face->numVerts;
				eNext = face->edges[j];
				ePrev = face->edges[prevJ];
				vert = face->verts[j];

				// Check ordering of edge (to determine which sub-edge should be used)
				subEdgeNext = eNext->v0 != vert;
				subEdgePrev = ePrev->v0 != vert;

				// Starting from the corner node
				mloops[i+0].v = vert->newVertIdx;
				mloops[i+0].e = eNext->newMetaIdx*2 + subEdgeNext;
				// go to next edge
				mloops[i+1].v = eNext->newMetaIdx + ss->numVerts;
				mloops[i+1].e = face->newEdgeStartIdx + j;
				// then to midpoint
				mloops[i+2].v = vert->newVertIdx;
				mloops[i+2].e = face->newEdgeStartIdx + prevJ;
				// then to previous edge,
				mloops[i+3].v = ePrev->newMetaIdx + ss->numVerts;
				mloops[i+3].e = ePrev->newMetaIdx*2 + subEdgePrev;
				i += 4;
			}
		}
	}
}

void SL_copyNewEdges(SLSubSurf *ss, MEdge *medges)
{
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
	int i = 0;
	// First the original edges
	FOR_HASH(ss->it, ss->edges) {
		SLEdge *edge = (SLEdge*)BLI_ghashIterator_getValue(ss->it);
		medges[i+0].v1 = edge->v0->newVertIdx;
		medges[i+0].v2 = edge->newMetaIdx + ss->numVerts;
		medges[i+1].v1 = edge->newMetaIdx + ss->numVerts;
		medges[i+1].v2 = edge->v0->newVertIdx;
		i += 2;
	}
	// Then the faces
	FOR_HASH(ss->it, ss->faces) {
		SLFace *face = (SLFace*)BLI_ghashIterator_getValue(ss->it);
		if (face->numVerts == 3) {
			medges[i+0].v1 = ss->numVerts + face->edges[0]->newMetaIdx;
			medges[i+0].v2 = ss->numVerts + face->edges[2]->newMetaIdx;
			medges[i+1].v1 = ss->numVerts + face->edges[0]->newMetaIdx;
			medges[i+1].v2 = ss->numVerts + face->edges[1]->newMetaIdx;
			medges[i+2].v1 = ss->numVerts + face->edges[1]->newMetaIdx;
			medges[i+2].v2 = ss->numVerts + face->edges[2]->newMetaIdx;
			i += 3;
		} else {
			int j;
			for (j = 0; j < face->numVerts; j++) {
				medges[i].v1 = ss->numVerts + face->edges[j]->newMetaIdx;
				medges[i].v2 = face->newVertIdx;
				i++;
			}
		}
	}
}

void SL_copyNewVerts(SLSubSurf *ss, MVert *mverts)
{
	int i = 0;
	FOR_HASH(ss->it, ss->verts) {
		SLVert *vert = (SLVert*)BLI_ghashIterator_getValue(ss->it);
		copy_v3_v3(mverts[i].co, vert->sl_coords);
		//normal_float_to_short_v3(mverts[i].no, vert->normal);
		i++;
	}
	FOR_HASH(ss->it, ss->edges) {
		SLEdge *edge = (SLEdge*)BLI_ghashIterator_getValue(ss->it);
		copy_v3_v3(mverts[i].co, edge->sl_coords);
		//normal_float_to_short_v3(mverts[i].no, edge->normal);
		i++;
	}
	FOR_HASH(ss->it, ss->edges) {
		SLFace *face = (SLFace*)BLI_ghashIterator_getValue(ss->it);
		if (face->numVerts > 3) {
			copy_v3_v3(mverts[i].co, face->centroid);
			//normal_float_to_short_v3(mverts[i].no, face->normal);
			i++;
		}
	}
}

/////////////////////////////////////////////////////////////

void _nofreefp(void *x) {
	// Nothing to free, its just the pointer, or freed elsewhere
}

SLSubSurf* SL_SubSurf_new(int smoothing) {
	MemArena *ma = BLI_memarena_new((1<<16), "SL subsurf");
	SLSubSurf *ss = (SLSubSurf*)BLI_memarena_alloc(ma, sizeof(SLSubSurf));
	ss->verts = BLI_ghash_new(BLI_ghashutil_ptrhash, BLI_ghashutil_ptrcmp, "SL verts");
	ss->edges = BLI_ghash_new(BLI_ghashutil_ptrhash, BLI_ghashutil_ptrcmp, "SL edges");
	ss->faces = BLI_ghash_new(BLI_ghashutil_ptrhash, BLI_ghashutil_ptrcmp, "SL faces");

	ss->it = BLI_ghashIterator_new(ss->verts); // Note: Can be used for every hash anyway (after init)

	ss->numVerts = ss->numEdges = ss->numFaces = 0;
	ss->memArena = ma;
	ss->smoothing = smoothing;
	return ss;
}

void SL_SubSurf_free(SLSubSurf *ss) {
	BLI_ghashIterator_free(ss->it);

	BLI_ghash_free(ss->verts, _nofreefp, _nofreefp);
	BLI_ghash_free(ss->edges, _nofreefp, _nofreefp);
	BLI_ghash_free(ss->faces, _nofreefp, _nofreefp);

	BLI_memarena_free(ss->memArena);
}

/////////////////////////////////////////////////////////////
// Helpers for connecting & disconnecting edges-verts-faces
// Using linked lists here;

static void _vertAddFace(SLSubSurf *ss, SLVert *v, SLFace *f) {
	BLI_linklist_prepend_arena(&v->faces, f, ss->memArena);
	v->requiresUpdate = 1;
	v->numFaces++;
}

static void _vertAddEdge(SLSubSurf *ss, SLVert *v, SLEdge *e) {
	BLI_linklist_prepend_arena(&v->edges, e, ss->memArena);
	v->requiresUpdate = 1;
	v->numEdges++;
}

static void _edgeAddFace(SLSubSurf *ss, SLEdge *e, SLFace *f) {
	BLI_linklist_prepend_arena(&e->faces, f, ss->memArena);
	e->requiresUpdate = 1;
	e->numFaces++;
}

/////////////////////////////////////////////////////////////
// Misc helpers

static SLEdge *_sharedEdge(SLVert *v0, SLVert *v1) {
	LinkNode *edge0, *edge1;
	FOR_LIST(edge0, v0->edges) {
		FOR_LIST(edge1, v1->edges) {
			if ( edge0->link == edge1->link) {
				return (SLEdge*)edge0->link;
			}
		}
	}
	return NULL;
}

/////////////////////////////////////////////////////////////
// Note! Must be added as verts, then edges, then faces and removed in the opposite order

void SL_addVert(SLSubSurf *ss, void *hashkey, float coords[3], int seam) {
	SLVert *vert = BLI_memarena_alloc(ss->memArena, sizeof(SLVert));
	copy_v3_v3(vert->coords, coords);
	vert->edges = NULL;
	vert->faces = NULL;
	vert->numFaces = 0;
	vert->numEdges = 0;
	vert->seam = seam;
	// Add to hashmap
	BLI_ghash_insert(ss->verts, hashkey, vert);
	ss->numVerts++;
	vert->requiresUpdate = 1;
}

void SL_updateVert(SLSubSurf *ss, void *hashkey, float coords[3], int seam) {
	int i;
	LinkNode *it;
	SLVert *vert = BLI_ghash_lookup(ss->verts, hashkey);
	BLI_assert(vert != NULL);
	copy_v3_v3(vert->coords, coords);
	vert->seam = seam;
	// Connected edges and faces need to updated
	FOR_LIST(it, vert->faces) {
		SLFace *face = (SLFace*)it->link;
		face->requiresUpdate = 1;
		// Also effects all edges on the connected face (since the center point moves)
		for (i = 0; i < face->numVerts; i++) {
			face->edges[i]->requiresUpdate = 1;
			face->verts[i]->requiresUpdate = 1;
		}
	}
	vert->requiresUpdate = 1;
}

// Must be called after syncVert
void SL_addEdge(SLSubSurf *ss, void *hashkey, void *vertkey0, void *vertkey1, float sharpness) {
	SLEdge *edge = BLI_memarena_alloc(ss->memArena, sizeof(SLEdge));

	edge->v0 = BLI_ghash_lookup(ss->verts, vertkey0);
	edge->v1 = BLI_ghash_lookup(ss->verts, vertkey1);
	edge->faces = NULL;
	edge->numFaces = 0;
	edge->sharpness = sharpness;
	edge->requiresUpdate = 1;

	_vertAddEdge(ss, edge->v0, edge);
	_vertAddEdge(ss, edge->v1, edge);

	// Add to hashmap
	BLI_ghash_insert(ss->edges, hashkey, edge);
	ss->numEdges++;
}

void SL_updateEdge(SLSubSurf *ss, void *hashkey, float sharpness) {
	LinkNode *it;
	SLEdge *edge = BLI_ghash_lookup(ss->edges, hashkey);
	BLI_assert(edge != NULL);
	// This means that sharpness has changed (v0 and v1 aren't allowed to change, that would indicate a *new* edge)
	edge->sharpness = sharpness;
	// Affects the directly connected faces and verts
	edge->v0->requiresUpdate = 1;
	edge->v1->requiresUpdate = 1;
	FOR_LIST(it, edge->faces) {
		((SLFace*)it->link)->requiresUpdate = 1;
	}
}

// Must be called after syncEdge
void SL_addFace(SLSubSurf *ss, void *hashkey, int numVerts, void **vertkeys) {
	int i;
	SLEdge *edge;
	SLVert *vert, *nextVert;
	SLFace *face;

	// New face? Then;
	face = BLI_memarena_alloc(ss->memArena, sizeof(SLFace));
	// Static lists for faces (more convenient, predictable size)
	face->verts = BLI_memarena_alloc(ss->memArena, sizeof(SLVert*)*numVerts);
	face->edges = BLI_memarena_alloc(ss->memArena, sizeof(SLEdge*)*numVerts);

	face->numVerts = numVerts;
	face->requiresUpdate = 1;
	for (i = 0; i < numVerts; i++) {
		// Verts
		vert = (SLVert*)BLI_ghash_lookup(ss->verts, vertkeys[i]);
		face->verts[i] = vert;
		_vertAddFace(ss, vert, face);
		// Then edges
		nextVert = (SLVert*)BLI_ghash_lookup(ss->verts, vertkeys[(i+1) % numVerts]);
		edge = _sharedEdge(vert, nextVert);
		face->edges[i] = edge;
		_edgeAddFace(ss, edge, face);
	}

	// Add to hashmap
	BLI_ghash_insert(ss->faces, hashkey, face);
	ss->numFaces++;
}

/////////////////////////////////////////////////////////////
// Actual smoothing stuff

void SL_processSync(SLSubSurf *ss) {
	SLFace *face;
	SLEdge *edge;
	SLVert *vert;
	LinkNode *it;
	int i;
	float avgSharpness;
	int seamCount, sharpnessCount;
	int seam;

	// Compute centroid, used for smoothing and other things;
	printf("Computing face centroids\n");
	FOR_HASH(ss->it, ss->faces) {
		face = (SLFace*)BLI_ghashIterator_getValue(ss->it);
		if (face->requiresUpdate) continue;

		zero_v3(face->centroid);

		for (i = 0; i < face->numVerts; i++) {
			add_v3_v3(face->centroid, face->verts[i]->coords);
		}
		mul_v3_fl(face->centroid, 1.0f / face->numVerts );
	}
	// also for edges;
	FOR_HASH(ss->it, ss->edges) {
		int x;
		edge = (SLEdge*)BLI_ghashIterator_getValue(ss->it);
		if (!edge->requiresUpdate) continue;

		for (x = 0; x < 3; x++)
			edge->centroid[x] = 0.5*edge->v0->coords[x] + 0.5*edge->v1->coords[x];
	}

	printf("Computing vert smoothing\n");
	// Loop over vertices and smooth out the Stam-Loop subsurface coordinate;
	FOR_HASH(ss->it, ss->verts) {
		vert = (SLVert*)BLI_ghashIterator_getValue(ss->it);
		if (!vert->requiresUpdate) continue;
		if (!ss->smoothing) copy_v3_v3(vert->sl_coords, vert->coords);

		// Compute average sharpness and seam;
		seamCount = 0;
		sharpnessCount = 0;
		avgSharpness = 0.0f;
		seam = vert->seam;
		FOR_LIST(it, vert->edges) {
			edge = (SLEdge*)it->link;

			if (seam && edge->numFaces < 2)
				seamCount++;

			if (edge->sharpness != 0.0f) {
				sharpnessCount++;
				avgSharpness += edge->sharpness;
			}
		}

		if (sharpnessCount) {
			avgSharpness /= sharpnessCount;
			if ( avgSharpness > 1.0f ) {
				avgSharpness = 1.0f;
			}
		}

		if (seamCount < 2 || seamCount != vert->numEdges)
			seam = 0;


		// Now do the smoothing;
		{
			int avgCount, edgeMult;

			zero_v3(vert->sl_coords);

			// Original coordinate, weight 4 (is this correct?)
			madd_v3_v3fl(vert->sl_coords, vert->coords, 4);
			avgCount = 4;

			// Weights for edges are multiple of shared faces;
			FOR_LIST(it, vert->edges) {
				edge = (SLEdge*)it->link;
				edgeMult = edge->numFaces == 0 ? 1 : edge->numFaces;
				madd_v3_v3fl(vert->sl_coords, edge->centroid, edgeMult);
				avgCount += edgeMult;
			}

			FOR_LIST(it, vert->faces) {
				face = (SLFace*)it->link;
				if (face->numVerts > 3) {
					add_v3_v3(vert->sl_coords, face->centroid);
					avgCount++; // Note that the subdivided area is a quad for any ngon > 3
				}
			}

			mul_v3_fl(vert->sl_coords, 1.0f / avgCount );
		}

		// Deal with sharpness and seams
		// Code snipped converted from CCG (undocumented mystery code)
		if ((sharpnessCount > 1 && vert->numFaces) || seam) {
			int x;
			float q[3];

			if (seam) {
				avgSharpness = 1.0f;
				sharpnessCount = seamCount;
			}

			zero_v3(q);
			FOR_LIST(it, vert->edges) {
				edge = (SLEdge*)it->link;
				if (seam) {
					if (edge->numFaces < 2)
						add_v3_v3(q, edge->centroid);
				}
				else if (edge->sharpness != 0.0f) {
					add_v3_v3(q, edge->centroid);
				}
			}

			mul_v3_fl(q, 1.0f / sharpnessCount);

			if (sharpnessCount != 2 || seam) {
				/* q = q + (co - q) * avgSharpness */
				for (x = 0; x < 3; x++) q[x] += (vert->coords[x] - q[x])*avgSharpness;
			}

			/* r = co * 0.75 + q * 0.25 */
			for (x = 0; x < 3; x++) q[x] = vert->coords[x]*0.75f + q[x]*0.25f;

			/* nCo = nCo + (r - nCo) * avgSharpness */
			for (x = 0; x < 3; x++) vert->sl_coords[x] += (q[x] - vert->sl_coords[x]) * avgSharpness;
		}
	}

	printf("Computing edge smoothing\n");
	// Loop over edges and smooth
	FOR_HASH(ss->it, ss->edges) {
		edge = (SLEdge*)BLI_ghashIterator_getValue(ss->it);
		if (!edge->requiresUpdate) continue;

		// Create the interpolated coordinates
		if (!ss->smoothing || edge->numFaces < 2 || edge->sharpness >= 1.0f) { // If its an edge, or maximum sharpness, then just average.
			copy_v3_v3(edge->sl_coords, edge->centroid);
		} else { // Otherwise smooth
			int avgCount;
			// Now, this is a bit tricky to deal with ngons. Quads and tris only would be a lot simpler.
			// Problem is to take into account edges appropriately.

			copy_v3_v3(edge->sl_coords, edge->v0->coords);
			add_v3_v3(edge->sl_coords, edge->v1->coords);
			mul_v3_fl(edge->sl_coords, 2);
			avgCount = 4;

			FOR_LIST(it, edge->faces) {
				face = (SLFace*)it->link;
				if (face->numVerts == 3) {
					// Triangles are split differently from the rest;
					// There are connections to the center nodes of the two opposite edges
					for (i = 0; i < face->numVerts; i++) {
						SLEdge *tempEdge = face->edges[i];
						if ( tempEdge != edge ) { // Then opposite edge
							madd_v3_v3fl(edge->sl_coords, tempEdge->centroid, 2);
						}
					}
					avgCount += 4; // 2 edges each
				} else {
					// Otherwise all ngons are split into quads, leaving one center node and edges
					madd_v3_v3fl(edge->sl_coords, face->centroid, 2);
					avgCount += 2;
					// Now find the other edges that share a node;
					for (i = 0; i < face->numVerts; i++) {
						SLEdge *tempEdge = face->edges[i];
						if ( tempEdge != edge ) {
							// Check for a shared node;
							if ( tempEdge->v0 == edge->v0 ||
									tempEdge->v0 == edge->v1 ||
									tempEdge->v1 == edge->v0 ||
									tempEdge->v1 == edge->v1) {
								add_v3_v3(edge->sl_coords, tempEdge->centroid);
							}
						}
					}
					avgCount += 4; // 1x2 from centroid + 2x1 from edges
				}
				madd_v3_v3fl(edge->sl_coords, face->centroid, 2);
			}
			mul_v3_fl(edge->sl_coords, 1.0f / (2.0f + edge->numFaces));

			// And take into account sharpness
			if (edge->sharpness > 0.0f ) {
				int x;
				for (x = 0; x < 3; x++) {
					edge->sl_coords[x] += edge->sharpness * (edge->centroid[x] - edge->sl_coords[x]);
				}
			}
		}
	}

	// Loop over faces and smooth
	/*
	FOR_HASH(ss->it, ss->faces) {
		face = (SLFace*)BLI_ghashIterator_getValue(ss->it);
	}*/
}

/////////////////////////////////////////////////////////////
