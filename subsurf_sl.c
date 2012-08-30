/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2005 Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

// Contributor(s): Mikael Öhman

/** \file blender/blenkernel/intern/subsurf_sl.c
 *  \ingroup bke
 */

#include <stdlib.h>
#include <stdio.h>

#include "DNA_mesh_types.h"
#include "DNA_scene_types.h"

#include "BKE_modifier.h"
#include "BKE_scene.h"
#include "BKE_subsurf.h"

#include "SLSubSurf.h"

#include "BLI_math_vector.h"

/**
 * Copies data from loop/poly to face (assumes all triangles or quads, as generated by most subdivisions)
 */
static void loops_to_customdata_corners(DerivedMesh *output)
{
	int i, numPolys, numTex, numCol, hasPCol;

	MPoly *polys;
	
	//if (!CustomData_get_layer_n(&output->loopData, CD_MLOOPCOL, n)) {
		// Not sure how this works with "n" ? 
	//	CustomData_add_layer(&output->loopData, CD_MLOOPCOL, CD_CALLOC, NULL, output->numVertData);
	//}

	numCol = CustomData_number_of_layers(&output->loopData, CD_MLOOPCOL);
	numTex = CustomData_number_of_layers(&output->loopData, CD_MLOOPUV);
	hasPCol = CustomData_has_layer(&output->loopData, CD_PREVIEW_MLOOPCOL);
	
	polys = output->getPolyArray(output);
	numPolys = output->getNumPolys(output);

	for (i = 0; i < numTex; i++) {
		int k;
		MTexPoly *texpoly = CustomData_get_layer_n(&output->polyData, CD_MTEXPOLY, i);
		MLoopUV *mloopuv = CustomData_get_layer_n(&output->loopData, CD_MLOOPUV, i);
		MTFace *texface = CustomData_get_layer_n(&output->faceData, CD_MTFACE, i);
		if (!texface) {
			texface = CustomData_add_layer(&output->faceData, CD_MTFACE, CD_CALLOC, NULL, numPolys);
		}

		for (k = 0; k < numPolys; k++) {
			int j;
			MPoly *poly = &polys[k];
			MLoopUV *ml = &mloopuv[poly->loopstart];
			
			ME_MTEXFACE_CPY(&texface[k], &texpoly[k]);
			for (j = 0; j < poly->totloop; j++) {
				copy_v2_v2(texface[k].uv[j], ml[j].uv);
			}
		}
	}

	for (i = 0; i < numCol; i++) {
		int k;
		MLoopCol *mloopcol = CustomData_get_layer_n(&output->loopData, CD_MLOOPCOL, i);
		MCol *mcol = CustomData_get_layer_n(&output->faceData, CD_MCOL, i);
		if (!mcol) {
			mcol = CustomData_add_layer(&output->faceData, CD_MCOL, CD_CALLOC, NULL, numPolys * 4);
		}
		for (k = 0; k < numPolys; k++) {
			int j;
			MPoly *poly = &polys[k];
			MLoopCol *ml = &mloopcol[poly->loopstart];
			for (j = 0; j < poly->totloop; j++) {
				MESH_MLOOPCOL_TO_MCOL(&ml[j], &mcol[j]);
			}
			mcol += 4;
		}
	}
	
	if (hasPCol) {
		int k;
		MLoopCol *mloopcol = CustomData_get_layer(&output->loopData, CD_PREVIEW_MLOOPCOL);
		MCol *mcol = CustomData_get_layer(&output->faceData, CD_PREVIEW_MCOL);
		if (!mcol) {
			mcol = CustomData_add_layer(&output->faceData, CD_MCOL, CD_CALLOC, NULL, numPolys * 4);
		}
		for (k = 0; k < numPolys; k++) {
			int j;
			MPoly *poly = &polys[k];
			MLoopCol *ml = &mloopcol[poly->loopstart];
			for (j = 0; j < poly->totloop; j++) {
				MESH_MLOOPCOL_TO_MCOL(&ml[j], &mcol[j]);
			}
			mcol += 4;
		}
	}
}


/*
 * Entrypoint for modifier. Should return the output DM based on the input DM.
 */
struct DerivedMesh *sl_subsurf_make_derived_from_derived(
	struct DerivedMesh *input,
	struct SubsurfModifierData *smd,
	float (*vertCos)[3],
	SubsurfFlags flags)
{
	SLSubSurf *ss;
	int i,j;
	int numCol, numTex;
	int levels;
	int smoothing = smd->subdivType == ME_SL_SUBSURF;
	int useSubsurfUv = smd->flags & eSubsurfModifierFlag_SubsurfUv;
	// TODO drawInterior
	//int drawInteriorEdges = !(smd->flags & eSubsurfModifierFlag_ControlEdges);
	DerivedMesh *result, *prevResult;
	
	// We don't use caches for this modifier
	if (smd->emCache) {
		smd->cacheFree(smd->emCache);
		smd->emCache = NULL;
	}
	if (smd->mCache) {
		smd->cacheFree(smd->mCache);
		smd->mCache = NULL;
	}

	if (flags & SUBSURF_USE_RENDER_PARAMS) {
		levels = (smd->modifier.scene) ? get_render_subsurf_level(&smd->modifier.scene->r, smd->renderLevels) : smd->renderLevels;
	}
	else /*if (flags & SUBSURF_FOR_EDIT_MODE)*/ {
		levels = (smd->modifier.scene) ? get_render_subsurf_level(&smd->modifier.scene->r, smd->levels) : smd->levels;
	}
	// Not sure if this is acceptable, but it needs to be dealt with:
	if (levels == 0) return input;

	// Now create the structures recursively, freeing all the middle steps (this gets a bit messy...)
	numCol = CustomData_number_of_layers(&input->loopData, CD_MLOOPCOL);
	numTex = CustomData_number_of_layers(&input->loopData, CD_MLOOPUV);
	
	prevResult = input;
	for (i = 0; i < levels; i++) {
		if (i == 0) {
			prevResult = input;
		} else {
			prevResult = result;
		}
		ss = SL_SubSurf_new(smoothing, prevResult, vertCos);
		result = SL_SubSurf_constructOutput(ss);
		SL_syncVertData(ss, result); // Have to before syncing verts themselves.
		SL_syncVerts(ss, result);
		for (j = 0; j < numCol; j++) {
			SL_syncPaint(ss, result, j);
		}
		for (j = 0; j < numTex; j++) {
			SL_syncUV(ss, result, useSubsurfUv, j);
		}
		if (prevResult != input && i < levels - 1) { // We are not allowed to remove the input DM, and not the last DM either (we need to tesselate)
			prevResult->release(prevResult);
		}
	}
	// We know the what the output DM contains, so we can optimize the tess face generation;
	//slDM_constructTessFaces(result);
	// TODO: All poly custom data should be copied over to tess faces as CD_REFERENCE, which would save time and memory.
	//result->recalcTessellation(result);
	slDM_constructTessFaces(result);
	if (prevResult != input) { // We are not allowed to remove the input DM, and not the last DM either (we need to tesselate)
		prevResult->release(prevResult);
	}
	
	// The preview color (for weight painting) needs to be copied over manually;
	
	
	loops_to_customdata_corners(result);
	
	
	return result;
}
