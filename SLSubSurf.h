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
typedef struct MLoopCol MLoopCol;
typedef struct MeshElemMap MeshElemMap;
typedef struct DerivedMesh DerivedMesh;
typedef struct SLSubSurf SLSubSurf;

/// The output DM for the Stam/Loop subdivision modifier.
typedef struct SLDerivedMesh SLDerivedMesh;

/**
 * \file LSSurf.h
 * This code implements a basic subdivision surface with triangle -> 4 * triangle, ngons -> n * quads.
 * Named after Stam/Loop, but many possible smoothing schemes exist.
 * Only subdivides 1 level!
 */

struct SLSubSurf {
    int smoothing; // Boolean, nonzero for smoothing.

    int numVerts;
    int numEdges;
    int numFaces;

	// If dm changes (fundamentally) the entire SLSubSurf should be removed and recreated.
	DerivedMesh *input;
	float (*vertexCos)[3]; // Vertex coordinates.
	// For convenience
	MVert *mvert;
	MEdge *medge;
	MLoop *mloop;
	MPoly *mpoly;

	// Maps necessary for quick access to the smoothing part
	// (TODO: it is possible that this could be performed without these, but its much simpler with them).
	int *poly2vert; // Maps old faces to new vert indices.
	MeshElemMap *edge2poly;
	int *edge2poly_mem;
	MeshElemMap *vert2poly;
	int *vert2poly_mem;
	MeshElemMap *vert2edge;
	int *vert2edge_mem;
};

int SL_giveTotalNumberOfSubVerts(SLSubSurf *ss);
int SL_giveTotalNumberOfSubEdges(SLSubSurf *ss);
int SL_giveTotalNumberOfSubFaces(SLSubSurf *ss);
int SL_giveTotalNumberOfSubLoops(SLSubSurf *ss);

// Methods for obtaining the loops and edges internal to the face. Uses the newIdx variables for numbering.
// Sufficient memory should be allocated by caller. Returns number of loops for sub face (will always be 3 or 4)
void SL_copyNewPolys(SLSubSurf *ss, MPoly *mpolys);
void SL_copyNewLoops(SLSubSurf *ss, MLoop *mloops);
void SL_copyNewEdges(SLSubSurf *ss, MEdge *medges);
void SL_copyNewTessFaces(SLSubSurf *ss, MFace *mfaces);

void SL_getMinMax(SLSubSurf *ss, float min_r[3], float max_r[3]);

SLSubSurf* SL_SubSurf_new(int smoothing, DerivedMesh *input, float (*vertexCos)[3]);
DerivedMesh *SL_SubSurf_constructOutput(SLSubSurf *ss);
void SL_SubSurf_free(SLSubSurf *ss);

void SL_syncVerts(SLSubSurf *ss, DerivedMesh *output);
void SL_syncUV(SLSubSurf *ss, DerivedMesh *output, int useSubsurfUv, int n);
void SL_syncPaint(SLSubSurf *ss, DerivedMesh *output, int n);
