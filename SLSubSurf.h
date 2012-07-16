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

#include "BLI_linklist.h"
#include "BLI_ghash.h"

typedef struct SLSubSurf SLSubSurf;
typedef struct SLFace SLFace;
typedef struct SLEdge SLEdge;
typedef struct SLVert SLVert;

/**
 * \file LSSurf.h
   This code is meant to implement a basic Stam/Loop subdivision surface
 * \note The code is intentionally written to avoid recursion. Complicates/limits some things.
 */


struct SLVert {
    int hashkey;
    float coords[3]; // Initial coordinate

    LinkNode *edges;
    LinkNode *faces;
    unsigned short numEdges, numFaces;

    unsigned short requiresUpdate, seam;
    float crease; // Support node creases as well (why not?)

    // Smoothed position;
    float sl_coords[3];
};

struct SLEdge {
    int hashkey;

    SLVert *v0, *v1;
    LinkNode *faces;
    unsigned short numFaces;
    unsigned short requiresUpdate;

    float crease;

    float centroid[3];

    // Smoothing center node position
    float sl_coords[3];
};

struct SLFace {
    int hashkey;

    LinkNode *verts;
    LinkNode *edges;
    unsigned short numVerts; // note: numVerts same as numEdges
    unsigned short requiresUpdate;

    float centroid[3];
};

struct SLSubSurf {
    int smoothing; // Boolean, nonzero for smoothing.

 	GHash *verts, *edges, *faces;
 	GHashIterator *vertIter, *edgeIter, *faceIter;

    int numVerts;
    int numEdges;
    int numFaces;
};

int SL_giveNumberOfInternalFaces(SLFace *face);
int SL_giveNumberOfInternalNodes(SLFace *face);
int SL_giveNumberOfInternalEdges(SLFace *face);

int SL_giveTotalNumberOfSubVerts(SLSubSurf *ss);
int SL_giveTotalNumberOfSubEdges(SLSubSurf *ss);
int SL_giveTotalNumberOfSubFaces(SLSubSurf *ss);


SLSubSurf* SL_SubSurf_new(int smoothing); // Allocators here? 
void SL_SubSurf_free(SLSubSurf *ss);


// This code basically adds verts i suppose.
void SL_SubSurf_syncVert(SLSubSurf *ss, void* hashkey, float coords[3], int seam);
void SL_SubSurf_syncEdge(SLSubSurf *ss, void* hashkey, SLVert *v0, SLVert *v1, float crease);
void SL_SubSurf_syncFace(SLSubSurf *ss, void* hashkey, int numVerts, SLVert **vs);


