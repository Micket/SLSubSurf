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
#include "BLI_memarena.h"

typedef struct MLoop MLoop;
typedef struct MEdge MEdge;
typedef struct MPoly MPoly;

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
	int newVertIdx; // Not sure if necessary (numbering is predictable)

    float coords[3]; // Initial coordinate

    LinkNode *edges;
    LinkNode *faces;
    unsigned short numEdges, numFaces;

    unsigned short requiresUpdate, seam;

    // Smoothed position;
    float sl_coords[3];
};

struct SLEdge {
	// Meta-index. Starts from 0 for the first edge, which contains the first 2 subedges.
	int newMetaIdx; // Not sure if necessary (numbering is predictable)

    SLVert *v0, *v1;
    LinkNode *faces;
    unsigned short numFaces;
    unsigned short requiresUpdate;

    float sharpness;

    float centroid[3];

    // Smoothing center node position
    float sl_coords[3];
};

struct SLFace {
	// New indices are given for original verts, edge nodes, then face nodes (some faces have no new node)
    int newVertIdx; // (unused for triangles)
	int newEdgeStartIdx;

    SLVert **verts;
    SLEdge **edges;
    unsigned short numVerts; // note: numVerts same as numEdges
    unsigned short requiresUpdate;

    float centroid[3];
};

struct SLSubSurf {
    MemArena *memArena;

    int smoothing; // Boolean, nonzero for smoothing.

    GHash *verts, *edges, *faces;
 	GHashIterator *it;

    int numVerts;
    int numEdges;
    int numFaces;
};

int SL_giveNumberOfInternalFaces(SLFace *face);
int SL_giveNumberOfInternalNodes(SLFace *face);
int SL_giveNumberOfInternalEdges(SLFace *face);
int SL_giveNumberOfInternalLoops(SLFace *face);

int SL_giveTotalNumberOfSubVerts(SLSubSurf *ss);
int SL_giveTotalNumberOfSubEdges(SLSubSurf *ss);
int SL_giveTotalNumberOfSubFaces(SLSubSurf *ss);
int SL_giveTotalNumberOfSubLoops(SLSubSurf *ss);

// Methods for obtaining the loops and edges internal to the face. Uses the newIdx variables for numbering.
// Sufficient memory should be allocated by caller. Returns number of loops for sub face (will always be 3 or 4)
void SL_giveSubLoopInFace(SLSubSurf *ss, SLFace *face, int *loopCount, int *polyCount, MLoop *mloops, MPoly *mpolys);
void SL_giveSubEdgeInFace(SLSubSurf *ss, SLFace *face,  MEdge *medges);

SLSubSurf* SL_SubSurf_new(int smoothing);
void SL_SubSurf_free(SLSubSurf *ss);

void SL_SubSurf_syncVert(SLSubSurf *ss, void* hashkey, float coords[3], int seam);
void SL_SubSurf_syncEdge(SLSubSurf *ss, void* hashkey, void *vertkey0, void *vertkey1, float sharpness);
void SL_SubSurf_syncFace(SLSubSurf *ss, void* hashkey, int numVerts, void **vertkeys);

void SL_SubSurf_processSync(SLSubSurf *ss);
