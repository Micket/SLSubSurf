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

/////////////////////////////////////////////////////////////
// External helpers 

int SL_giveTotalNumberOfSubVerts(SLSubSurf *ss) {
    // Simple things first, corners and interpolated edge verts;
    int totNodes = ss->numVerts + ss->numEdges; // One new node per edge
    // Then faces, which varies;
    BLI_ghashIterator_init(ss->faceIter, ss->faces);
    for (int i = 0; i < ss->numFaces; i++) {
        totNodes += SL_giveNumberOfInternalNodes((SLFace*)BLI_ghashIterator_getValue(ss->faceIter));
        BLI_ghashIterator_step(ss->faceIter);
    }
    return totNodes;
}

int SL_giveTotalNumberOfSubEdges(SLSubSurf *ss) {
    return ss->numEdges * 2;
}

int SL_giveTotalNumberOfSubFaces(SLSubSurf *ss) {
    // Since we know that all subdivided elements have the same number of sub-faces;
    int totFaces = 0; // One new node per edge
    // Then faces, which varies;
    BLI_ghashIterator_init(ss->faceIter, ss->faces);
    for (int i = 0; i < ss->numFaces; i++) {
        totFaces += SL_giveNumberOfInternalFaces((SLFace*)BLI_ghashIterator_getValue(ss->faceIter));
        BLI_ghashIterator_step(ss->faceIter);
    }

    return totFaces;
}

/////////////////////////////////////////////////////////////

inline void Vec3Zero(float a[3]) {
    for (int x = 0; x < 3; x++) a[x] = 0.0;
}

inline void Vec3Mult(float a[3], float b) {
    for (int x = 0; x < 3; x++) a[x] *= b;
}

inline void Vec3Add(float a[3], float b[3]) {
    for (int x = 0; x < 3; x++) a[x] += b[x];
}

inline void Vec3AddMult(float a[3], float b[3], float mult) {
    for (int x = 0; x < 3; x++) a[x] += mult*b[x];
}

inline void Vec3Copy(float a[3], float b[3]) {
    for (int x = 0; x < 3; x++) a[x] = b[x];
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
 

SLSubSurf* SL_SubSurf_new(int smoothing, MemArena *ma) {
    SLSubSurf *ss = (SLSubSurf*)BLI_memarena_alloc(ma, sizeof(SLSubSurf));
    ss->verts = BLI_ghash_new(BLI_ghashutil_ptrhash, BLI_ghashutil_ptrcmp, "SL verts");
    ss->edges = BLI_ghash_new(BLI_ghashutil_ptrhash, BLI_ghashutil_ptrcmp, "SL edges");
    ss->faces = BLI_ghash_new(BLI_ghashutil_ptrhash, BLI_ghashutil_ptrcmp, "SL faces");

    ss->vertIter = BLI_ghashIterator_new(ss->verts);
    ss->edgeIter = BLI_ghashIterator_new(ss->edges);
    ss->faceIter = BLI_ghashIterator_new(ss->faces);

    ss->numVerts = ss->numEdges = ss->numFaces = 0;
    ss->memArena = ma;
    ss->smoothing = smoothing;

    return ss;
}

void SL_SubSurf_free(SLSubSurf *ss) {
    BLI_ghashIterator_free(ss->vertIter);
    BLI_ghashIterator_free(ss->edgeIter);
    BLI_ghashIterator_free(ss->faceIter);

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
    }

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


static void _vertAddFace(SLVert *v, SLFace *f, SLSubSurf *ss) {
    BLI_linklist_prepend(&v->faces, f); // Am I using it correctly?!
    v->numFaces++;
}
static void _vertRemoveFace(SLVert *v, SLFace *f, SLSubSurf *ss) {
    BLI_linklist_remove(&v->faces, f);
    v->numFaces--;
}

///// 

static void _vertAddEdge(SLVert *v, SLEdge *e, SLSubSurf *ss) {
    BLI_linklist_prepend(&v->edges, e);
    v->numEdges++;
}

static void _vertRemoteEdge(SLVert *v, SLEdge *e, SLSubSurf *ss) {
    BLI_linklist_remove(&v->edges, e);
    v->numEdges--;
}

///// 

static void _edgeAddFace(SLEdge *e, SLFace *f, SLSubSurf *ss) {
    BLI_linklist_prepend(&e->faces, f);
    e->numFaces++;
}
static void _edgeRemoveFace(SLEdge *e, SLFace *f, SLSubSurf *ss) {
    BLI_linklist_remove(&e->faces, f);
    e->numFaces--;
}

/////////////////////////////////////////////////////////////
// Misc helpers

static SLEdge *_sharedEdge(SLVert *v0, SLVert *v1) {
    LinkNode *edge0 = v0->edges;
    for (int i = 0; i < v0->numEdges; i++) {
        LinkNode *edge1 = v1->edges;
        for (int j = 0; j < v1->numEdges; j++) {
            if ( edge0->link == edge1->link) {
                return (SLEdge*)edge0->link;
            }
            edge1 = edge1->next;
        }
        edge0 = edge0->next;
    }
    return NULL; 
}

/////////////////////////////////////////////////////////////
// Note! Must be added as verts, then edges, then faces. 

void SL_SubSurf_syncVert(SLSubSurf *ss, void *hashkey, float coords[3], int seam) {
    // Naive code for now, should use a hashmap of some sort.
    SLVert *vert = BLI_memarena_alloc(ss->memArena, sizeof(SLVert));
    Vec3Copy(vert->coords, coords);
    vert->numFaces = 0;
    vert->numEdges = 0;
    vert->requiresUpdate = 1;
    vert->seam = 1;
    
    // Add to hashmap
    BLI_ghash_insert(ss->verts, hashkey, vert);
}

void SL_SubSurf_syncEdge(SLSubSurf *ss, void *hashkey, SLVert *v0, SLVert *v1, float crease) {
    // Naive code for now, should use a hashmap of some sort.
    SLEdge *edge = BLI_memarena_alloc(ss->memArena, sizeof(SLEdge));

    edge->v0 = v0;
    edge->v1 = v1;
    edge->crease = crease;
    edge->requiresUpdate = 1;

    _vertAddEdge(v0, edge, ss);
    _vertAddEdge(v1, edge, ss);

    // Add to hashmap
    BLI_ghash_insert(ss->edges, hashkey, edge);
}

void SL_SubSurf_syncFace(SLSubSurf *ss, void *hashkey, int numVerts, SLVert **vs) {
    // Naive code for now, should use a hashmap of some sort, not just a simple linked list.

    // New face? Then;
    SLEdge *edge;
    SLFace *face = BLI_memarena_alloc(ss->memArena, sizeof(SLFace));

    // Static lists for faces maybe?
    //face->verts = (SLVert**)malloc(sizeof(SLVert*)*numVerts);
    //face->edges = (SLEdge**)malloc(sizeof(SLEdge*)*numVerts);
    // Allocate memory;
    //face->edges = (SLEdge**)malloc(sizeof(SLEdge*)*numVerts);

    face->numVerts = numVerts;
    face->requiresUpdate = 1;

    for (int i = 0; i < numVerts; i++) {
        // Verts
        BLI_linklist_prepend(&face->verts, vs[i]);
        _vertAddFace(vs[i], face, ss);
        // Then edges
        edge = _sharedEdge(vs[i], vs[(i+1) % numVerts]);
        BLI_linklist_prepend(&face->edges, edge);
        _edgeAddFace(edge, face, ss);
    }
    
    // Add to hashmap
    BLI_ghash_insert(ss->faces, hashkey, face);
}

/////////////////////////////////////////////////////////////
// Actual smoothing stuff

void SL_SubSurf_subdivideAll(SLSubSurf *ss) {
    SLFace *face;
    SLEdge *edge;
    SLVert *vert;
    LinkNode *temp;
    float avgCrease;
    int seamCount, creaseCount;
    int seam;

    // Compute centroid, used for smoothing and other things;
    BLI_ghashIterator_init(ss->faceIter, ss->faces);
    for (; !BLI_ghashIterator_isDone(ss->faceIter); BLI_ghashIterator_step(ss->faceIter))  {
        face = (SLFace*)BLI_ghashIterator_getValue(ss->faceIter);
        if (face->requiresUpdate) continue;

        Vec3Zero(face->centroid);
        for (int j = 0; j < face->numVerts; j++) {
            temp = face->verts;
            Vec3Add(face->centroid, ((SLVert*)temp->link)->coords);
            temp = temp->next;
        }
        Vec3Mult(face->centroid, 1.0f / face->numVerts );
    }
    // also for edges;
    BLI_ghashIterator_init(ss->edgeIter, ss->edges);
    for (; !BLI_ghashIterator_isDone(ss->edgeIter); BLI_ghashIterator_step(ss->edgeIter))  {
        edge = (SLEdge*)BLI_ghashIterator_getValue(ss->edgeIter);
        if (!edge->requiresUpdate) continue;

        for (int x = 0; x < 3; x++)
           edge->sl_coords[x] = 0.5*edge->v0->coords[x] + 0.5*edge->v1->coords[x];
    }

    // Loop over vertices and smooth out the Stam-Loop subsurface coordinate;
    BLI_ghashIterator_init(ss->vertIter, ss->verts);
    for (; !BLI_ghashIterator_isDone(ss->vertIter); BLI_ghashIterator_step(ss->vertIter))  {
        vert = (SLVert*)BLI_ghashIterator_getValue(ss->vertIter);
        if (!vert->requiresUpdate) continue;

        // Compute average crease and seam;
        seamCount = 0;
        creaseCount = 0;
        avgCrease = 0.0f;
        seam = vert->seam;
        temp = vert->edges;
		for (int j = 0; j < vert->numEdges; j++) {
			edge = (SLEdge*)temp->link;

			if (seam && edge->numFaces < 2)
				seamCount++;

			if (edge->crease != 0.0f) {
				creaseCount++;
				avgCrease += edge->crease;
			}
            temp = temp->next;
		}

		if (creaseCount) {
			avgCrease /= creaseCount;
			if ( avgCrease > 1.0f ) {
				avgCrease = 1.0f;
			}
		}

		if (seamCount < 2 || seamCount != vert->numEdges)
			seam = 0;


        // Now do the smoothing;
		{
			int numQuads = 0, numEdges = 0, edgeMult;

		    Vec3Zero(vert->sl_coords);
            temp = vert->faces;

            // Original coordinate, weight 4 (is this correct?)
			Vec3AddMult(vert->sl_coords, vert->coords, 4);

            // Weights for edges are multiple of shared faces;
            temp = vert->edges;
			for (int j = 0; j < vert->numEdges; j++) {
				edge = (SLEdge*)temp->link;
                edgeMult = edge->numFaces == 0 ? 1 : edge->numFaces;
				Vec3AddMult(vert->sl_coords, edge->centroid, edgeMult);
                numEdges += edgeMult;
                temp = temp->next;
			}

			for (int j = 0; j < vert->numFaces; j++) {
				face = (SLFace*)temp->link;
                if (face->numVerts > 3) {
    				Vec3Add(vert->sl_coords, face->centroid);
		    		numQuads++; // Note that the subdivided area is a quad for any ngon > 3
                }
                temp = temp->next;
			}

			Vec3Mult(vert->sl_coords, 1.0f / ( 4 + numQuads + numEdges) );
		}

        for (int x = 0; x < 3; x++) 
            vert->sl_coords[x] = vert->coords[x];
        BLI_ghashIterator_step(ss->vertIter);
    }

    // Loop over edges and smooth
    BLI_ghashIterator_init(ss->edgeIter, ss->edges);
    for (; !BLI_ghashIterator_isDone(ss->edgeIter); BLI_ghashIterator_step(ss->edgeIter))  {
        edge = (SLEdge*)BLI_ghashIterator_getValue(ss->edgeIter);
        if (!edge->requiresUpdate) continue;

        // Create the interpolated coordinates
		if (edge->numFaces < 2 || edge->crease >= 1.0f) { // If its an edge, or maximum crease, then just average.
            for (int x = 0; x < 3; x++)
                edge->sl_coords[x] = 0.5*edge->v0->coords[x] + 0.5*edge->v1->coords[x];
		} else { // Otherwise smooth
            int avgCount;
            // Now, this is a bit tricky to deal with ngons. Quads and tris only would be a lot simpler.
            // Problem is to take into account edges appropriately.

			Vec3Copy(edge->sl_coords, edge->v0->coords);
			Vec3Add(edge->sl_coords, edge->v1->coords);
            Vec3Mult(edge->sl_coords, 2);
            avgCount = 4;

            temp = edge->faces;
			for (int j = 0; j < edge->numFaces; j++) {
				face = (SLFace*)temp->link;
                if (face->numVerts == 3) {
                    // Triangles are split differently from the rest;
                    // There are connections to the center nodes of the two opposite edges
                    LinkNode *tempEdgeIt = face->edges;
                    for (int k = 0; k < 3; k++) {
                        SLEdge *tempEdge = (SLEdge*)tempEdgeIt->link;
                        if ( tempEdge != edge ) { // Then opposite edge
        				    Vec3AddMult(edge->sl_coords, tempEdge->centroid, 2);
                        }
                        tempEdgeIt = tempEdgeIt->next;
                    }
                    avgCount += 4; // 2 edges each
                } else {
                    // Otherwise all ngons are split into quads, leaving one center node and edges
        			Vec3AddMult(edge->sl_coords, face->centroid, 2);
                    avgCount += 2;
                    // Now find the other edges that share a node;
                    LinkNode *tempEdgeIt = face->edges;
                    for (int k = 0; k < face->numVerts; k++) {
                        SLEdge *tempEdge = (SLEdge*)tempEdgeIt->link;
                        if ( tempEdge != edge ) {
                            // Check for a shared node;
                            if ( tempEdge->v0 == edge->v0 || 
                                 tempEdge->v0 == edge->v1 ||
                                 tempEdge->v1 == edge->v0 ||
                                 tempEdge->v1 == edge->v1) {
            				    Vec3Add(edge->sl_coords, tempEdge->centroid);
                            }
                        }
                        tempEdgeIt = tempEdgeIt->next;
                    }
                    avgCount += 4; // 1x2 from centroid + 2x1 from edges
                }
				Vec3AddMult(edge->sl_coords, face->centroid, 2);
                temp = temp->next;
			}
			Vec3Mult(edge->sl_coords, 1.0f / (2.0f + edge->numFaces));

            // And take into account crease
            if (edge->crease > 0.0f ) {
                for (int x = 0; x < 3; x++) {
                    edge->sl_coords[x] = edge->sl_coords[x] + edge->crease * 
                        ((edge->v0->coords[x] + edge->v1->coords[x])*0.5 - edge->sl_coords[x]);
                }
            }
		}

        BLI_ghashIterator_step(ss->edgeIter);
    }

    // Loop over faces and smooth
    /*BLI_ghashIterator_init(ss->faceIter, ss->faces);
    for (int i = 0; i < ss->numFaces; i++) {
        face = (SLFace*)BLI_ghashIterator_getValue(ss->faceIter);
        BLI_ghashIterator_step(ss->faceIter);
    }*/
}

/////////////////////////////////////////////////////////////
