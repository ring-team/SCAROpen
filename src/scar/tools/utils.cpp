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

#include <scar/tools/utils.h>

#include <algorithm>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/mesh_index.h>

namespace SCAR
{
    bool is_element_in(
        const index_t element,
        const std::vector<index_t> &container)
    {
        if (container.empty())
        {
            return false;
        }
        return std::find(container.begin(), container.end(), element) != container.end() ? true : false;
    }

    bool is_element_in(
        const index_t element,
        const std::set<index_t> &container)
    {
        if (container.empty())
        {
            return false;
        }
        return std::find(container.begin(), container.end(), element) != container.end() ? true : false;
    }

    bool add_element_if_not_in(
        const index_t element,
        std::vector<index_t> &container)
    {
        if (!is_element_in(element, container))
        {
            container.push_back(element);
            return true;
        }
        return false;
    }

    void add_elements_if_not_in(
        const std::vector<index_t> &elements,
        std::vector<index_t> &container)
    {
        for (index_t element = 0; element < elements.size(); element++)
        {
            add_element_if_not_in(elements[element], container);
        }
    }

    bool are_all_elements_in(
        const std::vector<index_t> &elements,
        const std::vector<index_t> &container)
    {
        /// @todo Test and profile with std::includes
        for (index_t element = 0; element < elements.size(); element++)
        {
            if (!is_element_in(elements[element], container))
            {
                return false;
            }
        }
        return true;
    }

    bool are_all_elements_not_in(
        const std::vector<index_t> &elements,
        const std::vector<index_t> &container)
    {
        /// @todo Test and profile with std::includes
        for (index_t element = 0; element < elements.size(); element++)
        {
            if (is_element_in(elements[element], container))
            {
                return false;
            }
        }
        return true;
    }

    bool are_vectors_egal_up_to_permutations(
        const std::vector<index_t> &first,
        const std::vector<index_t> &second)
    {
        /// Check first vector cardinalities
        if (first.size() != second.size())
        {
            return false;
        }
        /// Check if A is included in B
        if (!are_all_elements_in(first, second))
        {
            return false;
        }
        /// Check if B is included in A
        if (!are_all_elements_in(second, first))
        {
            return false;
        }
        /// They are equals up to a permutation of their elements
        return true;
    }

    template <index_t DIMENSION>
    bool is_in_mesh_entity_boundaries(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        const RINGMesh::gmme_id &mesh_entity_id,
        const index_t query)
    {
        const RINGMesh::GeoModelMeshEntity<DIMENSION> &entity =
            geomodel.mesh_entity(mesh_entity_id);
        for (index_t b = 0; b < entity.nb_boundaries(); ++b)
        {
            if (entity.boundary_gmme(b).index() == query)
            {
                return true;
            }
        }
        return false;
    }

    template <index_t DIMENSION>
    std::set<index_t> surface_corners(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        index_t surface_id)
    {
        std::set<index_t> results;
        const auto &surface = geomodel.surface(surface_id);
        for (auto b_i : RINGMesh::range(surface.nb_boundaries()))
        {
            const auto &line = surface.boundary(b_i);
            for (auto b_b_i : RINGMesh::range(line.nb_boundaries()))
            {
                results.insert(line.boundary(b_b_i).index());
            }
        }
        return results;
    }

    template <index_t DIMENSION>
    bool is_point_on_border(
        const RINGMesh::SurfaceMesh<DIMENSION> &surface_mesh,
        index_t vertex_id,
        index_t first_polygon)
    {
        std::vector<index_t> polygons_around = surface_mesh.polygons_around_vertex(
            vertex_id, true, first_polygon);
        return !polygons_around.empty();
    }

    template <index_t DIMENSION>
    std::set<index_t> get_one_ring_vertices(
        const RINGMesh::SurfaceMesh<DIMENSION> &surface_mesh,
        index_t vertex_id,
        index_t first_polygon)
    {
        std::vector<index_t> polygons_around = surface_mesh.polygons_around_vertex(
            vertex_id, false, first_polygon);
        std::set<index_t> one_ring_vertices;
        for (auto p : polygons_around)
        {
            for (auto v : RINGMesh::range(surface_mesh.nb_polygon_vertices(p)))
            {
                auto cur_neigh_v_id = surface_mesh.polygon_vertex({p, v});
                if (cur_neigh_v_id == vertex_id)
                {
                    continue;
                }
                one_ring_vertices.insert(cur_neigh_v_id);
            }
        }
        return one_ring_vertices;
    }

    template bool scar_api is_in_mesh_entity_boundaries(
        const RINGMesh::GeoModel2D &,
        const RINGMesh::gmme_id &,
        const index_t);
    template bool scar_api is_point_on_border(
        const RINGMesh::SurfaceMesh2D &,
        index_t,
        index_t);
    template std::set<index_t> scar_api surface_corners(
        const RINGMesh::GeoModel2D &,
        index_t);
    template std::set<index_t> scar_api get_one_ring_vertices(
        const RINGMesh::SurfaceMesh2D &,
        index_t,
        index_t);

    template bool scar_api is_in_mesh_entity_boundaries(
        const RINGMesh::GeoModel3D &,
        const RINGMesh::gmme_id &,
        const index_t);
    template bool scar_api is_point_on_border(
        const RINGMesh::SurfaceMesh3D &,
        index_t,
        index_t);
    template std::set<index_t> scar_api surface_corners(
        const RINGMesh::GeoModel3D &,
        index_t);
    template std::set<index_t> scar_api get_one_ring_vertices(
        const RINGMesh::SurfaceMesh3D &,
        index_t,
        index_t);

}
