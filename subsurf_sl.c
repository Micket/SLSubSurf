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

	prevResult = input;
	for (i = 0; i < levels; i++) {
		ss = SL_SubSurf_new(smoothing, prevResult, vertCos);
		result = SL_SubSurf_constructOutput(ss);
		SL_syncVerts(ss, result);
		numCol = CustomData_number_of_layers(&input->loopData, CD_MLOOPCOL);
		for (j = 0; j < numCol; j++) {
			SL_syncPaint(ss, result, j);
		}
		numTex = CustomData_number_of_layers(&input->loopData, CD_MLOOPUV);
		for (j = 0; j < numTex; j++) {
			SL_syncUV(ss, result, useSubsurfUv, j);
		}
		if (prevResult != input) {
			prevResult->release(prevResult);
		}
		prevResult = result;
	}
	// We know the what the output DM contains, so we can optimize the tess face generation;
	SL_constructTessFaces(ss, result);
	
	// Not clear where this should be done. Its inconsistent if its in object mode vs edit mode.
	result->calcNormals(result);
	
	// TODO: All poly custom data should be copied over to tess faces as CD_REFERENCE, which would save time and memory.
	//result->recalcTessellation(result);
	
	return result;
}
