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
// Support functions

int SL_giveNumberOfSubFaces(int subdivLevel, SLFace *face) {
    return (1 << (subdivLevel << 1)); // 4^s = 2^(2*s) (correct for quads and tris, not sure how to deal with larger polygons)
}

int SL_giveNumberOfInternalFaceNodes(int subdivLevel, SLFace *face) {
    if (face->numVerts == 3) { // Triangles are special
        int e = 1 << subdivLevel;
        return ( (e - 2) * (e - 1) ) >> 1;
    } else { // For anything higher, replicate Catmull-Clark
        // Only correct for subdivlevel > 0;
        int en = SL_giveNumberOfInternalEdgeNodes(subdivLevel);
        // One subgrid for each vertex;
        return 1 + en * face->numVerts;
    }
}

int SL_giveNumberSubEdges(int subdivLevel) {
    return (1 << subdivLevel); // 2^s
}

int SL_giveNumberInternalEdgeNodes(int subdivLevel) {
    return (1 << subdivLevel) - 1; // 2^s
}


/////////////////////////////////////////////////////////////

int SL_giveTotalNumberOfSubVerts(SLSubSurf *ss) {
    // Simple things first, corners and interpolated edge verts;
    int totNodes = ss->numVerts + ss->numEdges * ( SL_giveNumberOfSubEdges(ss->subdivLevel) - 1 );
    // Then faces, which varies;
    BLI_ghashIterator_init(ss->faceIter, ss->faces);
    for (int i = 0; i < ss->numFaces; i++) {
        totNodes += SL_giveNumberOfInternalFaceNodes(ss->subdivLevel, (SLFace*)BLI_ghashIterator_getValue(ss->faceIter));
        BLI_ghashIterator_step(ss->faceIter);
    }
    return totNodes;
}

int SL_giveTotalNumberOfSubEdges(SLSubSurf *ss) {
    return ss->numEdges * SL_giveNumberOfSubEdges(ss->subdivLevel);
}

int SL_giveTotalNumberOfSubFaces(SLSubSurf *ss) {
    // Since we know that all subdivided elements have the same number of sub-faces;
    return ss->numFaces * SL_giveNumberOfSubFaces(ss->subdivLevel, NULL);
}

/////////////////////////////////////////////////////////////

SLSubSurf* SL_SubSurf_new(int subdivisionLevels, int smoothing) {
    SLSubSurf *ss = (SLSubSurf*)malloc(sizeof(SLSubSurf));
    ss->subdivLevel = subdivisionLevels;
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
    // TODO: Just loop and compare, for any reasonable mesh the number of edges should be very few.
    return NULL; 
}

/////////////////////////////////////////////////////////////
// Note! Must be added as verts, then edges, then faces. 

void SL_SubSurf_syncVert(SLSubSurf *ss, void *hashkey, double coords[3], int seam) {
    // Naive code for now, should use a hashmap of some sort.
    SLVert *vert = malloc(sizeof(SLVert));
    vert->coords[0] = coords[0];
    vert->coords[1] = coords[1];
    vert->coords[2] = coords[2];
    vert->numFaces = 0;
    vert->numEdges = 0;
    vert->requiresUpdate = 1;
    
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
// Actual smoothing stuff

void SL_SubSurf_subdivideAll(SLSubSurf *ss) {
    SLFace *face;
    SLEdge *edge;
    SLVert *vert;
    LinkNode *temp;
    double weight;
    double s0, s1, s2; // Interpolated coordinates;
    int intEdgeNodes;

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

    // Loop over vertices
    BLI_ghashIterator_init(ss->vertIter, ss->verts);
    for (int i = 0; i < ss->numEdges; i++) {
        vert = (SLVert*)BLI_ghashIterator_getValue(ss->vertIter);
        // TODO: Actually do the smoothing part...
        for (int x = 0; x < 3; x++) 
            vert->ls_coords[x] = vert->coords[x];
        BLI_ghashIterator_step(ss->vertIter);
    }

    // Loop over edges;
    BLI_ghashIterator_init(ss->edgeIter, ss->edges);
    intEdgeNodes = SL_giveNumberInternalEdgeNodes(ss->subdivLevel);
    for (int i = 0; i < ss->numEdges; i++) {
        edge = (SLEdge*)BLI_ghashIterator_getValue(ss->edgeIter);
        // Loop and create the interpolated coordinates
        // TODO: Actually do the smoothing part...
        for (int j = 1; j <= intEdgeNodes; j++) {
            s0 = ((double)j)/intEdgeNodes;
            for (int x = 0; x < 3; x++)
                edge->ls_coords[i][x] = (1-s0)*edge->v0->coords[x] + s0*edge->v1->coords[x];
        }
        BLI_ghashIterator_step(ss->edgeIter);
    }

    BLI_ghashIterator_init(ss->faceIter, ss->faces);
    for (int i = 0; i < ss->numFaces; i++) {
        face = (SLFace*)BLI_ghashIterator_getValue(ss->faceIter);
        SLVert *v0, *v1, *v2;
        temp = face->verts;
        v0 = (SLVert*)temp->link;
        temp = temp->next;
        v1 = (SLVert*)temp->link;
        temp = temp->next;
        v2 = (SLVert*)temp->link;
        // Loop and create the interpolated coordinates
        // TODO: Deal with more than just triangles...
        // TODO: Actually do the smoothing part...
        int k = 0;
        for (int j0 = 1; j0 < intEdgeNodes; j0++) {
            s0 = ((double)j0)/intEdgeNodes;
            for (int j1 = 1; j1 < intEdgeNodes-j0; j1++) {
                s1 = ((double)j1)/intEdgeNodes;
                s2 = 1.0 - s1 - s2; // Trenary coordinates

                for (int x = 0; x < 3; x++) 
                    face->ls_coords[k][x] = s0*v0->coords[x] + 
                                            s1*v1->coords[x] + 
                                            s2*v2->coords[x];
                k++;
            }
        }
        BLI_ghashIterator_step(ss->faceIter);
    }
}

/////////////////////////////////////////////////////////////
