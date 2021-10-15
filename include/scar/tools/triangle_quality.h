// Copyright (c) 2021 ASGA and Universite de Lorraine. All Rights Reserved.

// This program (SCAR) was developed in the frame of the RING project managed by
// ASGA and Universite de Lorraine.

// It is distributed under a dual licensing scheme:

// 1. RING Consortium License
// Members of the RING-COCAD Consortium may only use this file in
// accordance with the terms of described in the GOCAD Advancement Agreement,
// without the prior written authorization of the ASGA.
// Licencee agrees to attach or embed this Notice on all copies
// of the program, including partial copies or modified versions thereof.
// Please use: contact at ring dash team dot org, for more information.

// 2. GNU General Public License Usage
// Alternatively, this file may be used under the terms of the
// GNU General Public license version 3. The licenses are as published by
// the Free Software Foundation and appearing in the file
// https://www.gnu.org/licenses/gpl-3.0.html

#pragma once

#include <scar/basic/common.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModel);
    FORWARD_DECLARATION_DIMENSION_CLASS(SurfaceMesh);
    ALIAS_3D(SurfaceMesh);
}

namespace SCAR
{

    enum struct TriangleQualityMode
    {
        ALL,
        MIN_ANGLE,
        MAX_ANGLE,
        ANGLE_DIFFERENCE,
        RADIUS_RATIO,
        EDGE_RATIO,
        EDGE_TO_CIRCUMRADIUS,
        EDGE_TO_INRADIUS,
        MIN_HEIGHT
    };

    /*!
     * @brief Computes and stores triangular mesh quality in the GeoModel.
     *
     * @param[in] mesh_qual_mode mesh quality to compute.
     * @param[in,out] geomodel GeoModel in which the mesh quality is performed.
     * The quality is stored on the polygons of each Surface.
     *
     * @warning The GeoModel must have at least one surface. All the surfaces
     * must be meshed by simplexes (triangles).
     */
    template <index_t DIMENSION>
    void scar_api compute_prop_triangle_mesh_quality(
        TriangleQualityMode mesh_qual_mode,
        const RINGMesh::GeoModel<DIMENSION> &geomodel);

    /*!
     * @brief Fill the /p output_mesh with cells of quality below \p min_quality
     * @param[in] mesh_qual_mode mesh quality of cells.
     * @param[in] min_quality Value of quality below which a cell is add in the
     * mesh.
     * @param[in] geomodel GeoModel in which the mesh quality is read.
     * @param[out] output_mesh VolumeMesh to fill with low quality cells.
     * @returns The minimum value of cell quality
     *
     * @warning The GeoModel must have at least one region. All the regions
     * must be meshed by simplexes (tetrahedra).
     */
    template <index_t DIMENSION>
    void scar_api output_quality_values(
        const std::string &filename,
        TriangleQualityMode mesh_qual_mode,
        const RINGMesh::GeoModel<DIMENSION> &geomodel);
}
