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

typedef struct MVert MVert;
typedef struct MEdge MEdge;
typedef struct MLoop MLoop;
typedef struct MPoly MPoly;
typedef struct MFace MFace;
typedef struct DMFlagMat DMFlagMat;
typedef struct MemArena MemArena;
typedef struct GHash GHash;
typedef struct GHashIterator GHashIterator;
typedef struct LinkNode LinkNode;

typedef struct SLSubSurf SLSubSurf;
typedef struct SLFace SLFace;
typedef struct SLEdge SLEdge;
typedef struct SLVert SLVert;

/**
 * \file LSSurf.h
 * This code implements a basic Stam/Loop subdivision surface
 * Only subdivides 1 level!
 */


struct SLVert {
	int newVertIdx; // Not sure if necessary (numbering is predictable)

    float coords[3]; // Initial coordinate

    LinkNode *edges;
    LinkNode *faces;
    unsigned short numEdges, numFaces;

    unsigned short requiresUpdate, seam;

	float sl_coords[3]; // Smoothed position
	float normal[3];
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
	
	float sl_coords[3]; // Smoothing center node position
	float normal[3];
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
	
	float normal[3];
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
void SL_copyNewPolys(SLSubSurf *ss, DMFlagMat *faceFlags, MPoly *mpolys);
void SL_copyNewLoops(SLSubSurf *ss, MLoop *mloops);
void SL_copyNewEdges(SLSubSurf *ss, MEdge *medges);
void SL_copyNewVerts(SLSubSurf *ss, MVert *mverts);
void SL_copyNewTessFaces(SLSubSurf *ss, DMFlagMat *faceFlags, MFace *mfaces);

SLSubSurf* SL_SubSurf_new(int smoothing);
void SL_SubSurf_free(SLSubSurf *ss);

void SL_addVert(SLSubSurf *ss, void* hashkey, float coords[3], int seam);
void SL_addEdge(SLSubSurf *ss, void* hashkey, void *vertkey0, void *vertkey1, float sharpness);
void SL_addFace(SLSubSurf *ss, void* hashkey, int numVerts, void **vertkeys);
// Updates coordinate or sharpness in existing vert/edge
void SL_updateVert(SLSubSurf *ss, void* hashkey, float coords[3], int seam);
void SL_updateEdge(SLSubSurf *ss, void* hashkey, float sharpness);

void SL_getMinMax(SLSubSurf *ss, float min_r[3], float max_r[3]);

void SL_processSync(SLSubSurf *ss);
void SL_renumberAll(SLSubSurf *ss);