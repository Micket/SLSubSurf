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
	int i;
	int numCol, numTex;
	int levels;
	int smoothing = smd->subdivType == ME_SL_SUBSURF;
	int useSubsurfUv = smd->flags & eSubsurfModifierFlag_SubsurfUv;
	// TODO drawInterior
	//int drawInteriorEdges = !(smd->flags & eSubsurfModifierFlag_ControlEdges);
	DerivedMesh *result;
	
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
	// Not sure if this is acceptable;
	if (levels == 0) return input;
	
	ss = SL_SubSurf_new(smoothing, input, vertCos);
	result = SL_SubSurf_constructOutput(ss);
	SL_syncVerts(ss, result);
	numCol = CustomData_number_of_layers(&input->loopData, CD_MLOOPCOL);
	for (i = 0; i < numCol; i++) {
		SL_syncPaint(ss, result, i);
	}
	numTex = CustomData_number_of_layers(&input->loopData, CD_MLOOPUV);
	for (i = 0; i < numTex; i++) {
		SL_syncUV(ss, result, useSubsurfUv, i);
	}
	result->calcNormals(result);
	result->recalcTessellation(result);
	
	return result;
}

#if 0

// End of drawing functions

static SLDerivedMesh *getSLDerivedMesh(SLSubSurf *ss,
									   int drawInteriorEdges,
									   int useSubsurfUv,
									   DerivedMesh *source_dm)
{
	SLDerivedMesh *ssdm = MEM_callocN(sizeof(*ssdm), "sldm");
	DerivedMesh *newdm;
	int totsubvert, totsubedge, totsubface, totsubloop;
	int numTex, numCol;
	int hasPCol, hasOrigSpace;
	int *polyidx;
	int i;

	ssdm->ss = ss;
	ssdm->drawInteriorEdges = drawInteriorEdges;
	//ssdm->useSubsurfUv = useSubsurfUv;
	newdm = &(ssdm->dm);
	for (i = 0; i < sizeof(newdm); i++) {
		((char*)newdm)[i] = 0; // Nulling everything
	}

	//DM_init_funcs(newdm); // Sets some default functions we don't care/need to overload.

	totsubvert = SL_giveTotalNumberOfSubVerts(ss);
	totsubedge = SL_giveTotalNumberOfSubEdges(ss);
	totsubface = SL_giveTotalNumberOfSubFaces(ss);
	totsubloop = SL_giveTotalNumberOfSubLoops(ss);

	DM_from_template(newdm, source_dm, DM_TYPE_CCGDM, // TODO Change type; (general subsurf-type perhaps?)
					totsubvert,
					totsubedge,
					totsubface, // faces
					totsubloop, // loops TODO: unsure
					totsubface); // polys TODO: unsure

	CustomData_free_layer_active(&newdm->polyData, CD_NORMAL,
								 newdm->numPolyData);

	numTex = CustomData_number_of_layers(&newdm->loopData, CD_MLOOPUV);
	numCol = CustomData_number_of_layers(&newdm->loopData, CD_MLOOPCOL);
	hasPCol = CustomData_has_layer(&newdm->loopData, CD_PREVIEW_MLOOPCOL);
	hasOrigSpace = CustomData_has_layer(&newdm->loopData, CD_ORIGSPACE_MLOOP);

	if (
		(numTex && CustomData_number_of_layers(&newdm->faceData, CD_MTFACE) != numTex)  ||
		(numCol && CustomData_number_of_layers(&newdm->faceData, CD_MCOL) != numCol)    ||
		(hasPCol && !CustomData_has_layer(&newdm->faceData, CD_PREVIEW_MCOL))           ||
		(hasOrigSpace && !CustomData_has_layer(&newdm->faceData, CD_ORIGSPACE)) )
	{
		CustomData_from_bmeshpoly(&newdm->faceData,
								  &newdm->polyData,
								  &newdm->loopData,
								  totsubface);
	}

	// TODO: deal with uv's..
#if 0
	if (useSubsurfUv) {
		CustomData *ldata = &newdm->loopData;
		CustomData *dmldata = &source_dm->loopData;
		int numlayer = CustomData_number_of_layers(ldata, CD_MLOOPUV);
		int dmnumlayer = CustomData_number_of_layers(dmldata, CD_MLOOPUV);

		for (i = 0; i < numlayer && i < dmnumlayer; i++)
			set_subsurf_uv(ss, source_dm, newdm, i);
	}
#endif

	/* All tessellated CD layers were updated! */
	newdm->dirty &= ~DM_DIRTY_TESS_CDLAYERS;
	return ssdm;
}
#endif