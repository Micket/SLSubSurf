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
//#include "MEM_guardedalloc.h"
void MEM_freeN(void *ptr);

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

inline void Vec3Zero(float a[3]) {
	int x;
	for (x = 0; x < 3; x++) a[x] = 0.0;
}

inline void Vec3Mult(float a[3], float b) {
	int x;
	for (x = 0; x < 3; x++) a[x] *= b;
}

inline void Vec3Add(float a[3], float b[3]) {
	int x;
	for (x = 0; x < 3; x++) a[x] += b[x];
}

inline void Vec3AddMult(float a[3], float b[3], float mult) {
	int x;
	for (x = 0; x < 3; x++) a[x] += mult*b[x];
}

inline void Vec3Copy(float a[3], float b[3]) {
	int x;
	for (x = 0; x < 3; x++) a[x] = b[x];
}

/////////////////////////////////////////////////////////////

void _nofreefp(void *x) {
	// Nothing to free, its just the pointer, or freed elsewhere
}

// These free the memory allocated the the linknode itself, not the link, which is allocated/freed by the memory arena
static void _valfreeVert(void *val) {
	SLVert *vert = (SLVert*)val;
	BLI_linklist_free(vert->edges, NULL);
	BLI_linklist_free(vert->faces, NULL);
}
static void _valfreeEdge(void *val) {
	SLEdge *edge = (SLEdge*)val;
	BLI_linklist_free(edge->faces, NULL);
}
static void _valfreeFace(void *val) {
	SLFace *face = (SLFace*)val;
	BLI_linklist_free(face->verts, NULL);
	BLI_linklist_free(face->edges, NULL);
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

	BLI_ghash_free(ss->verts, _nofreefp, _valfreeVert);
	BLI_ghash_free(ss->edges, _nofreefp, _valfreeEdge);
	BLI_ghash_free(ss->faces, _nofreefp, _valfreeFace);

	BLI_memarena_free(ss->memArena);
}

/////////////////////////////////////////////////////////////
// Helpers for connecting & disconnecting edges-verts-faces
// Using linked lists here;

// Missing from BLI_linklist;
void BLI_linklist_remove(LinkNode **list, void *item) {
	if ((*list)->link == item) { // First entry should be removed
		// Called should free the item, we'll just stop the node;
		LinkNode *tmp = (*list)->next;
		MEM_freeN(*list);
		*list = tmp;
		return;
	} else {
		LinkNode *prev = *list;
		LinkNode *node = prev->next;
		while (node != NULL) {
			if (node->link == item) {
				prev->next = node->next; // Connect the previous node
				MEM_freeN(*list);
			}
			prev = node;
			node = node->next;
		}
	}
}


static void _vertAddFace(SLVert *v, SLFace *f, SLSubSurf *ss) {
	BLI_linklist_prepend(&v->faces, f);
	v->requiresUpdate = 1;
	v->numFaces++;
}
static void _vertRemoveFace(SLVert *v, SLFace *f, SLSubSurf *ss) {
	BLI_linklist_remove(&v->faces, f);
	v->requiresUpdate = 1;
	v->numFaces--;
}

/////

static void _vertAddEdge(SLVert *v, SLEdge *e, SLSubSurf *ss) {
	BLI_linklist_prepend(&v->edges, e);
	v->requiresUpdate = 1;
	v->numEdges++;
}

static void _vertRemoveEdge(SLVert *v, SLEdge *e, SLSubSurf *ss) {
	BLI_linklist_remove(&v->edges, e);
	v->requiresUpdate = 1;
	v->numEdges--;
}

/////

static void _edgeAddFace(SLEdge *e, SLFace *f, SLSubSurf *ss) {
	BLI_linklist_prepend(&e->faces, f);
	e->requiresUpdate = 1;
	e->numFaces++;
}
static void _edgeRemoveFace(SLEdge *e, SLFace *f, SLSubSurf *ss) {
	BLI_linklist_remove(&e->faces, f);
	e->requiresUpdate = 1;
	e->numFaces--;
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

void SL_SubSurf_syncVert(SLSubSurf *ss, void *hashkey, float coords[3], int seam) {
	SLVert *vert;

	if ( (vert = BLI_ghash_lookup(ss->verts, hashkey)) == NULL ) { // Then new vert
		SLVert *vert = BLI_memarena_alloc(ss->memArena, sizeof(SLVert));
		Vec3Copy(vert->coords, coords);
		vert->edges = NULL;
		vert->faces = NULL;
		vert->numFaces = 0;
		vert->numEdges = 0;
		vert->seam = seam;

		// Add to hashmap
		BLI_ghash_insert(ss->verts, hashkey, vert);
		ss->numVerts++;
	} else {
		// Then existing vert has moved (TODO: Should I split this function?)
		LinkNode *it;
		Vec3Copy(vert->coords, coords);
		vert->seam = seam;
		// Connected edges and faces need to updated
		FOR_LIST(it, vert->faces) {
			LinkNode *it2;
			SLFace *face = (SLFace*)it->link;
			face->requiresUpdate = 1;
			// Also effects all edges on the connected face (since the center point moves)
			FOR_LIST(it2, face->edges) {
				((SLEdge*)it->link)->requiresUpdate = 1;
			}
			FOR_LIST(it2, face->verts) {
				((SLVert*)it->link)->requiresUpdate = 1;
			}
		}

	}
	vert->requiresUpdate = 1;
}

// Must be called after syncVert
void SL_SubSurf_syncEdge(SLSubSurf *ss, void *hashkey, void *vertkey0, void *vertkey1, float sharpness) {
	SLEdge *edge;

	if ( (edge = BLI_ghash_lookup(ss->edges, hashkey)) == NULL ) {
		// Then new edge
		SLEdge *edge = BLI_memarena_alloc(ss->memArena, sizeof(SLEdge));

		edge->v0 = BLI_ghash_lookup(ss->verts, vertkey0);
		edge->v1 = BLI_ghash_lookup(ss->verts, vertkey1);
		edge->faces = NULL;
		edge->sharpness = sharpness;
		edge->requiresUpdate = 1;

		_vertAddEdge(edge->v0, edge, ss);
		_vertAddEdge(edge->v1, edge, ss);

		// Add to hashmap
		BLI_ghash_insert(ss->edges, hashkey, edge);
	} else {
		// This means that sharpness has changed (v0 and v1 aren't allowed to change, that would indicate a *new* edge)
		LinkNode *it;
		edge->sharpness = sharpness;
		// Affects the directly connected faces and verts
		edge->v0->requiresUpdate = 1;
		edge->v1->requiresUpdate = 1;
		FOR_LIST(it, edge->faces) {
			((SLFace*)it->link)->requiresUpdate = 1;
		}
	}
	ss->numEdges++;
}

// Must be called after syncEdge
void SL_SubSurf_syncFace(SLSubSurf *ss, void *hashkey, int numVerts, void **vertkeys) {
	int i;
	SLEdge *edge;
	SLVert *vert, *nextVert;
	SLFace *face;

	if ( BLI_ghash_lookup(ss->faces, hashkey) != NULL ) {
		return; // Nothing can change be updated for existing faces.
	}

	// New face? Then;

	face = BLI_memarena_alloc(ss->memArena, sizeof(SLFace));

	// Static lists for faces maybe?
	//face->verts = (SLVert**)MEM_allocN(sizeof(SLVert*)*numVerts);
	//face->edges = (SLEdge**)MEM_allocN(sizeof(SLEdge*)*numVerts);

	face->numVerts = numVerts;
	face->requiresUpdate = 1;

	for (i = 0; i < numVerts; i++) {
		// Verts
		vert = (SLVert*)BLI_ghash_lookup(ss->verts, vertkeys[i]);
		BLI_linklist_prepend(&face->verts, vert);
		_vertAddFace(vert, face, ss);
		// Then edges
		nextVert = (SLVert*)BLI_ghash_lookup(ss->verts, vertkeys[(i+1) % numVerts]);
		edge = _sharedEdge(vert, nextVert);
		BLI_linklist_prepend(&face->edges, edge);
		_edgeAddFace(edge, face, ss);
	}

	// Add to hashmap
	BLI_ghash_insert(ss->faces, hashkey, face);
	ss->numFaces++;
}

// And then deletion

// Must be called after syncEdgeDel
void SL_SubSurf_syncVertDel(SLSubSurf *ss, void *hashkey) {
	BLI_ghash_remove(ss->edges, hashkey, _nofreefp, _valfreeVert);
	ss->numVerts--;
}

// Must be called after syncFaceDel!
void SL_SubSurf_syncEdgeDel(SLSubSurf *ss, void *hashkey) {
	// Disconnect edge;
	SLEdge *edge = (SLEdge*)BLI_ghash_lookup(ss->edges, hashkey);
	_vertRemoveEdge(edge->v0, edge, ss);
	_vertRemoveEdge(edge->v1, edge, ss);
	// then delete it
	BLI_ghash_remove(ss->edges, hashkey, _nofreefp, _valfreeEdge);
	ss->numEdges--;
}

void SL_SubSurf_syncFaceDel(SLSubSurf *ss, void *hashkey) {
	// Disconnect face;
	LinkNode *itVert, *itEdge;
	SLFace *face = (SLFace*)BLI_ghash_lookup(ss->edges, hashkey);
	for (itVert = face->verts, itEdge = face->edges; itVert != NULL; itVert = itVert->next, itEdge = itEdge->next) {
		_vertRemoveFace((SLVert*)itVert->link, face, ss);
		_edgeRemoveFace((SLEdge*)itEdge->link, face, ss);
	}
	// then delete it
	BLI_ghash_remove(ss->edges, hashkey, _nofreefp, _valfreeFace);
	ss->numFaces--;
}

/////////////////////////////////////////////////////////////
// Actual smoothing stuff

void SL_SubSurf_processSync(SLSubSurf *ss) {
	SLFace *face;
	SLEdge *edge;
	SLVert *vert;
	LinkNode *it;
	float avgSharpness;
	int seamCount, sharpnessCount;
	int seam;

	// Compute centroid, used for smoothing and other things;
	FOR_HASH(ss->it, ss->faces) {
		face = (SLFace*)BLI_ghashIterator_getValue(ss->it);
		if (face->requiresUpdate) continue;

		Vec3Zero(face->centroid);

		FOR_LIST(it, face->verts) {
			Vec3Add(face->centroid, ((SLVert*)it->link)->coords);
		}
		Vec3Mult(face->centroid, 1.0f / face->numVerts );
	}
	// also for edges;
	FOR_HASH(ss->it, ss->edges) {
		int x;
		edge = (SLEdge*)BLI_ghashIterator_getValue(ss->it);
		if (!edge->requiresUpdate) continue;

		for (x = 0; x < 3; x++)
			edge->sl_coords[x] = 0.5*edge->v0->coords[x] + 0.5*edge->v1->coords[x];
	}

	// Loop over vertices and smooth out the Stam-Loop subsurface coordinate;
	FOR_HASH(ss->it, ss->verts) {
		vert = (SLVert*)BLI_ghashIterator_getValue(ss->it);
		if (!vert->requiresUpdate) continue;
		if (!ss->smoothing) Vec3Copy(vert->sl_coords, vert->coords);

		// Compute average sharpness and seam;
		seamCount = 0;
		sharpnessCount = 0;
		avgSharpness = 0.0f;
		seam = vert->seam;
		FOR_LIST(it, face->verts) {
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

			Vec3Zero(vert->sl_coords);

			// Original coordinate, weight 4 (is this correct?)
			Vec3AddMult(vert->sl_coords, vert->coords, 4);
			avgCount = 4;

			// Weights for edges are multiple of shared faces;
			FOR_LIST(it, vert->edges) {
				edge = (SLEdge*)it->link;
				edgeMult = edge->numFaces == 0 ? 1 : edge->numFaces;
				Vec3AddMult(vert->sl_coords, edge->centroid, edgeMult);
				avgCount += edgeMult;
			}

			FOR_LIST(it, vert->faces) {
				face = (SLFace*)it->link;
				if (face->numVerts > 3) {
					Vec3Add(vert->sl_coords, face->centroid);
					avgCount++; // Note that the subdivided area is a quad for any ngon > 3
				}
			}

			Vec3Mult(vert->sl_coords, 1.0f / ( 4 + avgCount) );
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

			Vec3Zero(q);
			FOR_LIST(it, vert->edges) {
				edge = (SLEdge*)it->link;
				if (seam) {
					if (edge->numFaces < 2)
						Vec3Add(q, edge->centroid);
				}
				else if (edge->sharpness != 0.0f) {
					Vec3Add(q, edge->centroid);
				}
			}

			Vec3Mult(q, 1.0f / sharpnessCount);

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

	// Loop over edges and smooth
	FOR_HASH(ss->it, ss->edges) {
		edge = (SLEdge*)BLI_ghashIterator_getValue(ss->it);
		if (!edge->requiresUpdate) continue;

		// Create the interpolated coordinates
		if (!ss->smoothing || edge->numFaces < 2 || edge->sharpness >= 1.0f) { // If its an edge, or maximum sharpness, then just average.
			Vec3Copy(edge->sl_coords, edge->centroid);
		} else { // Otherwise smooth
			int avgCount;
			// Now, this is a bit tricky to deal with ngons. Quads and tris only would be a lot simpler.
			// Problem is to take into account edges appropriately.

			Vec3Copy(edge->sl_coords, edge->v0->coords);
			Vec3Add(edge->sl_coords, edge->v1->coords);
			Vec3Mult(edge->sl_coords, 2);
			avgCount = 4;

			FOR_LIST(it, edge->faces) {
				face = (SLFace*)it->link;
				if (face->numVerts == 3) {
					LinkNode *it2;
					// Triangles are split differently from the rest;
					// There are connections to the center nodes of the two opposite edges
					FOR_LIST(it2, vert->edges) {
						SLEdge *tempEdge = (SLEdge*)it2->link;
						if ( tempEdge != edge ) { // Then opposite edge
							Vec3AddMult(edge->sl_coords, tempEdge->centroid, 2);
						}
					}
					avgCount += 4; // 2 edges each
				} else {
					LinkNode *it2;
					// Otherwise all ngons are split into quads, leaving one center node and edges
					Vec3AddMult(edge->sl_coords, face->centroid, 2);
					avgCount += 2;
					// Now find the other edges that share a node;
					FOR_LIST(it2, vert->edges) {
						SLEdge *tempEdge = (SLEdge*)it2->link;
						if ( tempEdge != edge ) {
							// Check for a shared node;
							if ( tempEdge->v0 == edge->v0 ||
									tempEdge->v0 == edge->v1 ||
									tempEdge->v1 == edge->v0 ||
									tempEdge->v1 == edge->v1) {
								Vec3Add(edge->sl_coords, tempEdge->centroid);
							}
						}
					}
					avgCount += 4; // 1x2 from centroid + 2x1 from edges
				}
				Vec3AddMult(edge->sl_coords, face->centroid, 2);
			}
			Vec3Mult(edge->sl_coords, 1.0f / (2.0f + edge->numFaces));

			// And take into account sharpness
			if (edge->sharpness > 0.0f ) {
				int x;
				for (x = 0; x < 3; x++) {
					edge->sl_coords[x] = edge->sl_coords[x] + edge->sharpness * (edge->centroid[x] - edge->sl_coords[x]);
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
