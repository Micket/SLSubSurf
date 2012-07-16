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

SLSubSurf* SL_SubSurf_new(int smoothing) {
    SLSubSurf *ss = (SLSubSurf*)malloc(sizeof(SLSubSurf));
    ss->smoothing = smoothing;
    return ss;
}
 
void SL_SubSurf_free(SLSubSurf *ss) {
    //TODO
}

/////////////////////////////////////////////////////////////
// Helpers for connecting & disconnecting edges-verts-faces
// Using linked lists here;

static void _vertAddFace(SLVert *v, SLFace *f, SLSubSurf *ss)
{
    BLI_linklist_prepend(&v->faces, f); // Am I using it correctly?!
    v->numFaces++;
}
static void _vertRemoveFace(SLVert *v, SLFace *f, SLSubSurf *ss)
{
    //BLI_linklist_remove(&v->faces, f);
    v->numFaces--;
}

///// 

static void _vertAddEdge(SLVert *v, SLEdge *e, SLSubSurf *ss)
{
    BLI_linklist_prepend(&v->edges, e);
    v->numEdges++;
}
static void _vertRemoteEdge(SLVert *v, SLEdge *e, SLSubSurf *ss)
{
    //BLI_linklist_remove(&v->edges, e);
    v->numEdges--;
}

///// 

static void _edgeAddFace(SLEdge *e, SLFace *f, SLSubSurf *ss)
{
    BLI_linklist_prepend(&e->faces, f);
    e->numFaces++;
}
static void _edgeRemoveFace(SLEdge *e, SLFace *f, SLSubSurf *ss)
{
    //BLI_linklist_remove(&e->faces, f);
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
    SLVert *vert = malloc(sizeof(SLVert));
    vert->coords[0] = coords[0];
    vert->coords[1] = coords[1];
    vert->coords[2] = coords[2];
    vert->numFaces = 0;
    vert->numEdges = 0;
    vert->requiresUpdate = 1;
    vert->seam = 1;
    
    // Add to hashmap
    BLI_ghash_insert(ss->verts, hashkey, vert);
}

void SL_SubSurf_syncEdge(SLSubSurf *ss, void *hashkey, SLVert *v0, SLVert *v1, float crease) {
    // Naive code for now, should use a hashmap of some sort.
    SLEdge *edge = malloc(sizeof(SLEdge));

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
    SLFace *face = malloc(sizeof(SLFace));

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
// Actual smoothing stuff

void SL_SubSurf_subdivideAll(SLSubSurf *ss) {
    SLFace *face;
    SLEdge *edge;
    SLVert *vert;
    LinkNode *temp;
    float weight;
    float avgCrease;
    int seamCount, creaseCount;
    int seam;

    // Compute centroid, used for smoothing and other things;
    BLI_ghashIterator_init(ss->faceIter, ss->faces);
    for (int i = 0; i < ss->numFaces; i++) {
        face = (SLFace*)BLI_ghashIterator_getValue(ss->faceIter);
        weight = 1.0 / face->numVerts;
        face->centroid[0] = face->centroid[1] = face->centroid[2] = 0.0;
        for (int j = 0; j < face->numVerts; j++) {
            temp = face->verts;
            for (int x = 0; x < 3; x++) {
                face->centroid[x] += weight*((SLVert*)temp->link)->coords[x];
                temp = temp->next;
            }
        }
        BLI_ghashIterator_step(ss->faceIter);
    }

    // Loop over vertices and smooth out the Stam-Loop subsurface coordinate;
    BLI_ghashIterator_init(ss->vertIter, ss->verts);
    for (int i = 0; i < ss->numEdges; i++) {
        vert = (SLVert*)BLI_ghashIterator_getValue(ss->vertIter);
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

            // Original coordinate, weight ???
			Vec3AddMult(vert->coords, edge->sl_coords, 4);

            // Weights for edges are multiple of shared faces;
            temp = vert->edges;
			for (i = 0; i < vert->numEdges; i++) {
				edge = (SLEdge*)temp->link;
                edgeMult = edge->numFaces == 0 ? 1 : edge->numFaces;
				Vec3AddMult(vert->sl_coords, edge->sl_coords, edgeMult);
                numEdges += edgeMult;
                temp = temp->next;
			}

			for (i = 0; i < vert->numFaces; i++) {
				face = (SLFace*)temp->link;
                if (face->numVerts > 3) {
    				Vec3Add(vert->sl_coords, face->centroid);
		    		numQuads++; // Note that the subdivided area is a quad for any ngon > 3
                }
                temp = temp->next;
			}

			Vec3Mult(vert->sl_coords, 1.0f / ( 4 + numQuads + numEdges) );
		}


        // TODO: Actually do the smoothing part...
        for (int x = 0; x < 3; x++) 
            vert->sl_coords[x] = vert->coords[x];
        BLI_ghashIterator_step(ss->vertIter);
    }

    // Loop over edges and smooth
    BLI_ghashIterator_init(ss->edgeIter, ss->edges);
    for (int i = 0; i < ss->numEdges; i++) {
        edge = (SLEdge*)BLI_ghashIterator_getValue(ss->edgeIter);
        // Loop and create the interpolated coordinates

		if (edge->numFaces < 2 || edge->crease >= 1.0f) { // If its an edge, or maximum crease, then just average.
            for (int x = 0; x < 3; x++)
                edge->sl_coords[x] = 0.5*edge->v0->coords[x] + 0.5*edge->v1->coords[x];
		} else { // Otherwise smooth

            // TODO: Change this to Catmull-Clark to the proposed smoothing scheme by Stam-Loop
			Vec3Copy(edge->sl_coords, edge->v0->coords);
			Vec3Add(edge->sl_coords, edge->v1->coords);
            temp = edge->faces;
			for (int j = 0; j < edge->numFaces; j++) {
				SLFace *f = (SLFace*)temp->link;
				Vec3Add(edge->sl_coords, f->centroid);
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
