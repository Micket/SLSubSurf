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
    SLVert *next; // Linked list
    double coords[3]; // Initial coordinate

    SLEdge **edges;
    SLFace **faces;

    unsigned short numEdges, numFaces;
    unsigned short requiresUpdate;
    
    float crease; // Support node creases as well (why not?)

    // Smoothed position;
    double ls_coords[3];
};

struct SLEdge {
    SLEdge *next; // Linked list
    SLVert *v0, *v1;
    SLFace **faces;
    unsigned short numFaces;
    unsigned short requiresUpdate;

    float crease;

    // The subdivided data below (node positions)
     // Size is implicitly defined from the subdivlevel.
    double *ls_coords[3];
};

struct SLFace {
    SLFace *next; // Linked list
    SLVert **verts;
    SLEdge **edges;
    unsigned short numVerts; // note: numVerts same as numEdges
    unsigned short requiresUpdate;


    double centroid[3];
 
    // The subdivided data (internal to face, excluding initial edges)
     // Size is implicitly defined from the subdivlevel.
    double *ls_coords[3];
};

struct SLSubSurf {
    int subdivLevel;
    int smoothing; // Boolean, nonzero for smoothing.

 	SLVert *lastVert;
	SLEdge *lastEdge;
	SLFace *lastFace;

    int numVerts;
    int numEdges;
    int numFaces;
};

int SL_giveNumberOfSubFaces(int subdivLevel, SLFace *face);
int SL_giveNumberOfInternalFaceNodes(int subdivLevel, SLFace *face);
int SL_giveNumberOfSubEdges(int subdivLevel);
int SL_giveNumberOfInternalEdgeNodes(int subdivLevel);

int SL_giveTotalNumberOfSubVerts(SLSubSurf *ss);
int SL_giveTotalNumberOfSubEdges(SLSubSurf *ss);
int SL_giveTotalNumberOfSubFaces(SLSubSurf *ss);


SLSubSurf* SL_SubSurf_new(int subdivisionLevels, int smoothing); // Allocators here? 
void SL_SubSurf_free(SLSubSurf *ss);


// This code basically adds verts i suppose.
void SL_SubSurf_syncVert(SLSubSurf *ss, double coords[3], int seam);
void SL_SubSurf_syncEdge(SLSubSurf *ss, SLVert *v0, SLVert *v1, float crease);
void SL_SubSurf_syncFace(SLSubSurf *ss, int numVerts, SLVert **vs);


