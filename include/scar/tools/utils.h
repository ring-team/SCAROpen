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

#include <ringmesh/geomodel/core/entity_type.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModel);
    FORWARD_DECLARATION_DIMENSION_CLASS(SurfaceMesh);
}

namespace SCAR
{

    /*!
     * \name Queries on containers
     * @{
     */

    template <typename T>
    bool set_contains(
        const T element,
        const std::set<T> &container)
    {
        return container.find(element) != container.end();
    }

    bool scar_api is_element_in(
        const index_t element,
        const std::vector<index_t> &container);

    bool scar_api is_element_in(
        const index_t element,
        const std::set<index_t> &container);

    bool scar_api add_element_if_not_in(
        const index_t element,
        std::vector<index_t> &container);

    void scar_api add_elements_if_not_in(
        const std::vector<index_t> &elements,
        std::vector<index_t> &container);

    bool scar_api are_all_elements_in(
        const std::vector<index_t> &elements,
        const std::vector<index_t> &container);

    bool scar_api are_all_elements_not_in(
        const std::vector<index_t> &elements,
        const std::vector<index_t> &container);

    bool scar_api are_vectors_egal_up_to_permutations(
        const std::vector<index_t> &first,
        const std::vector<index_t> &second);

    /*! @}
     */

    /*!
     * \name Queries on GeoModels
     * @{
     */

    template <index_t DIMENSION>
    bool scar_api is_in_mesh_entity_boundaries(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        const RINGMesh::gmme_id &mesh_entity_id,
        const index_t query);

    template <index_t DIMENSION>
    std::set<index_t> scar_api surface_corners(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        index_t surface_id);

    template <index_t DIMENSION>
    bool scar_api is_point_on_border(
        const RINGMesh::SurfaceMesh<DIMENSION> &surface_mesh,
        index_t vertex_id,
        index_t first_polygon = NO_ID);

    template <index_t DIMENSION>
    std::set<index_t> scar_api get_one_ring_vertices(
        const RINGMesh::SurfaceMesh<DIMENSION> &surface_mesh,
        index_t vertex_id,
        index_t first_polygon = NO_ID);

    /*! @}
     */

    inline std::pair<index_t, index_t> scar_api double_integer_bounds(
        double value)
    {
        return std::pair<index_t, index_t>(
            static_cast<index_t>(std::floor(value)),
            static_cast<index_t>(std::ceil(value)));
    }
}
