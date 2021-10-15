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

#include <scar/tools/distance.h>

#include <queue>

#include <geogram/basic/attributes.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <scar/tools/convex_shape.h>
#include <scar/tools/geometry.h>
#include <scar/tools/utils.h>

namespace
{
    using namespace SCAR;

    template <index_t DIMENSION>
    double compute_line_line_one_sided_hausdorff_distance(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        const index_t line_from_id,
        const index_t line_to_id)
    {
        double max_distance = 0;
        const RINGMesh::Line<DIMENSION> &line_from = geomodel.line(line_from_id);
        const RINGMesh::Line<DIMENSION> &line_to = geomodel.line(line_to_id);
        for (index_t v = 0; v < line_from.nb_vertices(); ++v)
        {
            double cur_min_dist;
            std::tie(std::ignore, std::ignore, cur_min_dist) =
                line_to.edge_aabb().closest_edge(line_from.vertex(v));
            max_distance = std::max(max_distance, cur_min_dist);
        }
        return max_distance;
    }

    double compute_surface_surface_one_sided_hausdorff_distance(
        const RINGMesh::GeoModel3D &geomodel,
        const index_t surface_from_id,
        const index_t surface_to_id)
    {
        double max_distance = 0;
        const RINGMesh::Surface3D &surface_from = geomodel.surface(
            surface_from_id);
        const RINGMesh::Surface3D &surface_to = geomodel.surface(surface_to_id);
        for (index_t v = 0; v < surface_from.nb_vertices(); ++v)
        {
            double cur_min_dist;
            std::tie(std::ignore, std::ignore, cur_min_dist) =
                surface_to.polygon_aabb().closest_triangle(
                    surface_from.vertex(v));
            max_distance = std::max(max_distance, cur_min_dist);
        }
        return max_distance;
    }
}

namespace SCAR
{

    template <index_t DIMENSION>
    double compute_line_line_symmetric_hausdorff_distance(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        const index_t line1_id,
        const index_t line2_id)
    {
        return std::max(
            compute_line_line_one_sided_hausdorff_distance(geomodel, line1_id,
                                                           line2_id),
            compute_line_line_one_sided_hausdorff_distance(geomodel, line2_id,
                                                           line1_id));
    }

    double compute_surface_surface_symmetric_hausdorff_distance(
        const RINGMesh::GeoModel3D &geomodel,
        const index_t surface1_id,
        const index_t surface2_id)
    {
        return std::max(
            compute_surface_surface_one_sided_hausdorff_distance(geomodel,
                                                                 surface1_id, surface2_id),
            compute_surface_surface_one_sided_hausdorff_distance(geomodel,
                                                                 surface2_id, surface1_id));
    }

    template <index_t DIMENSION>
    double compute_curvilinear_distance(
        const RINGMesh::Line<DIMENSION> &line,
        const index_t line_vertex,
        const index_t boundary_corner_id)
    {
        scar_assert(boundary_corner_id == 0 || boundary_corner_id == 1);
        double curv_dist = 0;
        for (auto e_id : RINGMesh::range(line_vertex))
        {
            curv_dist += line.mesh_element_size(e_id);
        }
        if (boundary_corner_id == 0)
        {
            return curv_dist;
        }
        return line.size() - curv_dist;
    }

    void compute_geodesic_to_borders_djikstra_approximation(
        const RINGMesh::GeoModel3D &geomodel,
        const index_t surface_id)
    {
        const auto &surf = geomodel.surface(surface_id);
        GEO::AttributesManager &attr_manager = surf.vertex_attribute_manager();
        GEO::Attribute<double> dijk(attr_manager, geodesic_att_name);
        dijk.fill(max_double());

        std::queue<index_t> border_vertices;
        std::set<index_t> vertices_to_determine;

        for (auto e : RINGMesh::range(surf.nb_mesh_elements()))
        {
            for (auto v : RINGMesh::range(surf.nb_mesh_element_vertices(e)))
            {
                auto v_id = surf.mesh_element_vertex_index({e, v});
                if (is_point_on_border(surf.mesh(), v_id, e))
                {
                    dijk[v_id] = 0;
                    border_vertices.push(v_id);
                }
                else
                {
                    vertices_to_determine.insert(v_id);
                }
            }
        }

        while (!border_vertices.empty())
        {
            auto b_id = border_vertices.front();
            border_vertices.pop();
            auto cur_v = surf.vertex(b_id);
            auto neighbors = get_one_ring_vertices(surf.mesh(), b_id);
            for (auto n_id : neighbors)
            {
                if (is_element_in(n_id, vertices_to_determine))
                {
                    auto v_n = surf.vertex(n_id);
                    // dijk[b_id] = 0
                    dijk[n_id] = std::min(dijk[n_id], (v_n - cur_v).length());
                }
            }
        }

        while (!vertices_to_determine.empty())
        {
            // Find min value
            auto min_value = max_double();
            auto v_id = NO_ID;
            for (auto i : vertices_to_determine)
            {
                if (dijk[i] < min_value)
                {
                    v_id = i;
                    min_value = dijk[i];
                }
            }
            vertices_to_determine.erase(v_id);
            auto cur_v = surf.vertex(v_id);
            auto neighbors = get_one_ring_vertices(surf.mesh(), v_id);
            for (auto n_id : neighbors)
            {
                if (is_element_in(n_id, vertices_to_determine))
                {
                    auto v_n = surf.vertex(n_id);
                    dijk[n_id] = std::min(dijk[n_id],
                                          dijk[v_id] + (v_n - cur_v).length());
                }
            }
        }
    }

    double compute_line_vertices_surface_distance(
        const RINGMesh::GeoModel3D &geomodel,
        const index_t line_id,
        const index_t surface_id)
    {
        double max_distance = 0;
        const RINGMesh::Line3D &line = geomodel.line(line_id);
        const RINGMesh::Surface3D &surface = geomodel.surface(surface_id);
        for (index_t lv = 0; lv < line.nb_vertices(); ++lv)
        {
            double cur_distance;
            std::tie(std::ignore, std::ignore, cur_distance) =
                surface.polygon_aabb().closest_triangle(line.vertex(lv));
            max_distance = std::max(max_distance, cur_distance);
        }
        return max_distance;
    }

    bool line_line_minimum_distance(
        const vec3 &l0_pt0,
        const vec3 &l0_pt1,
        const vec3 &l1_pt0,
        const vec3 &l1_pt1,
        double &distance,
        vec3 &midpoint)
    {
        bool parallel_lines = false;
        // http://www.geometrictools.com/Source/Distance3D.html
        vec3 diff = l0_pt0 - l1_pt0;
        vec3 l0_dir = normalize(l0_pt1 - l0_pt0);
        vec3 l1_dir = normalize(l1_pt1 - l1_pt0);
        double a01 = -dot(l0_dir, l1_dir);
        double b0 = dot(diff, l0_dir);
        double s0, s1;

        if (std::abs(a01) < 1.)
        {
            // Lines not parallel
            double det = 1. - a01 * a01;
            double b1 = -dot(diff, l1_dir);
            s0 = (a01 * b1 - b0) / det;
            s1 = (a01 * b0 - b1) / det;
        }
        else
        {
            // Lines are parallel, select any pair of closest points.
            s0 = 0;
            s1 = 0;
            parallel_lines = true;
        }
        vec3 close_point_on_l0 = l0_pt0 + s0 * l0_dir;
        vec3 close_point_on_l1 = l1_pt0 + s1 * l1_dir;

        midpoint = 0.5 * (close_point_on_l0 + close_point_on_l1);
        distance = (close_point_on_l0 - close_point_on_l1).length();

        return parallel_lines;
    }

    bool segment_segment_distance(
        const vec3 &s0_pt0,
        const vec3 &s0_pt1,
        const vec3 &s1_pt0,
        const vec3 &s1_pt1,
        double &distance,
        vec3 &nearest_s0_pt,
        vec3 &nearest_s1_pt,
        bool &are_parallel)
    {
        nearest_s0_pt = NO_POINT_3D;
        nearest_s1_pt = NO_POINT_3D;

        vec3 line_midpoint;
        are_parallel = line_line_minimum_distance(s0_pt0, s0_pt1, s1_pt0, s1_pt1,
                                                  distance, line_midpoint);

        bool is_projection_possible_s1_on_s0;
        std::tie(is_projection_possible_s1_on_s0, nearest_s0_pt) =
            RINGMesh::point_segment_projection(line_midpoint, s0_pt0, s0_pt1);
        bool is_projection_possible_s0_on_s1;
        std::tie(is_projection_possible_s0_on_s1, nearest_s1_pt) =
            RINGMesh::point_segment_projection(line_midpoint, s1_pt0, s1_pt1);

        // Direct projection between segments
        if (is_projection_possible_s1_on_s0 && is_projection_possible_s0_on_s1)
        {
            scar_assert(
                nearest_s0_pt != NO_POINT_3D && nearest_s0_pt != NO_POINT_3D);
            distance = (nearest_s0_pt - nearest_s1_pt).length();
            return true;
        }

        // At least one segment could not be projected on the other
        double sq_dist_pt_00_to_pt_10 = (s0_pt0 - s1_pt0).length2();
        double sq_dist_pt_00_to_pt_11 = (s0_pt0 - s1_pt1).length2();
        double sq_dist_pt_01_to_pt_10 = (s0_pt1 - s1_pt0).length2();
        double sq_dist_pt_01_to_pt_11 = (s0_pt1 - s1_pt1).length2();

        if (!is_projection_possible_s1_on_s0)
        {
            double sq_dist_pt_00_to_seg1 = std::min(sq_dist_pt_00_to_pt_10,
                                                    sq_dist_pt_00_to_pt_11);
            double sq_dist_pt_01_to_seg1 = std::min(sq_dist_pt_01_to_pt_10,
                                                    sq_dist_pt_01_to_pt_11);
            if (sq_dist_pt_00_to_seg1 < sq_dist_pt_01_to_seg1)
            {
                nearest_s0_pt = s0_pt0;
            }
            else
            {
                nearest_s0_pt = s0_pt1;
            }
        }

        if (!is_projection_possible_s0_on_s1)
        {
            double sq_dist_pt_10_to_seg0 = std::min(sq_dist_pt_00_to_pt_10,
                                                    sq_dist_pt_01_to_pt_10);
            double sq_dist_pt_11_to_seg0 = std::min(sq_dist_pt_00_to_pt_11,
                                                    sq_dist_pt_01_to_pt_11);
            if (sq_dist_pt_10_to_seg0 < sq_dist_pt_11_to_seg0)
            {
                nearest_s1_pt = s1_pt0;
            }
            else
            {
                nearest_s1_pt = s1_pt1;
            }
        }

        scar_assert(nearest_s0_pt != NO_POINT_3D && nearest_s0_pt != NO_POINT_3D);
        distance = (nearest_s0_pt - nearest_s1_pt).length();
        return false;
    }

    template <index_t DIMENSION>
    NearbyLineChunksFromPoint<DIMENSION>::NearbyLineChunksFromPoint(
        const ConvexShape<DIMENSION> &tolerance_shape,
        const RINGMesh::LineMesh<DIMENSION> &line_mesh)
    {
        LineChunk cur_line_chunk;
        bool chunk_in_progress{false};
        // Initialization
        bool v0_inside = tolerance_shape.is_point_inside(line_mesh.vertex(0));
        if (v0_inside)
        {
            cur_line_chunk.begin = {0, 0.};
            chunk_in_progress = true;
        }
        for (index_t v = 1; v < line_mesh.nb_vertices(); ++v)
        {
            RINGMesh::Geometry::Segment<DIMENSION> previous_edge(
                line_mesh.vertex(v - 1), line_mesh.vertex(v));
            bool does_segment_intersect_tolerance;
            std::vector<vecn<DIMENSION>> intersections;

            std::tie(does_segment_intersect_tolerance, intersections) =
                tolerance_shape.segment_intersections(previous_edge);
            if (does_segment_intersect_tolerance)
            {
                std::vector<PointOnLineCoordinates<DIMENSION>> intersections_coord;
                intersections_coord.reserve(intersections.size());
                for (auto i : intersections)
                {
                    intersections_coord.emplace_back(i,
                                                     LineRelativeCoordinates{
                                                         v - 1, segment_point_normalized_coordinates(
                                                                    previous_edge, i)});
                }
                std::sort(intersections_coord.begin(), intersections_coord.end());
                for (auto i : intersections_coord)
                {
                    if (i.relative_coord == cur_line_chunk.begin)
                    {
                        continue;
                    }
                    if (chunk_in_progress)
                    {
                        cur_line_chunk.end = i.relative_coord;
                        chunks_.push_back(cur_line_chunk);
                        chunk_in_progress = false;
                    }
                    else
                    {
                        cur_line_chunk.begin = i.relative_coord;
                        chunk_in_progress = true;
                    }
                }
            }
        }

        // End of the line
        bool v_last_inside = tolerance_shape.is_point_inside(
            line_mesh.vertex(line_mesh.nb_vertices() - 1));
        if (v_last_inside)
        {
            cur_line_chunk.end = {line_mesh.nb_edges() - 1, 1.};
            chunks_.push_back(cur_line_chunk);
            chunk_in_progress = false;
        }
    }

    void show_chunks(const LineChunks &chunks)
    {
        for (auto chunk : chunks)
        {
            Logger::out("Chunk", "New chunk : ", "\n", "from : ",
                        chunk.begin.edge + chunk.begin.norm, "\n", "to : ",
                        chunk.end.edge + chunk.end.norm);
        }
    }

    template <index_t DIMENSION>
    double segment_point_normalized_coordinates(
        const RINGMesh::Geometry::Segment<DIMENSION> &oriented_segment,
        const vecn<DIMENSION> &point)
    {
        if (oriented_segment.length() < global_epsilon && (point - oriented_segment.p0).length() < global_epsilon)
        {
            return 0;
        }
        return std::sqrt(
            (point - oriented_segment.p0).length2() / (oriented_segment.p1 - oriented_segment.p0).length2());
    }

    template class scar_api NearbyLineChunksFromPoint<2>;
    template double scar_api compute_line_line_symmetric_hausdorff_distance(
        const RINGMesh::GeoModel2D &geomodel,
        const index_t,
        const index_t);
    template double scar_api segment_point_normalized_coordinates(
        const RINGMesh::Geometry::Segment2D &,
        const vec2 &);
    template double scar_api compute_curvilinear_distance(
        const RINGMesh::Line2D &line,
        const index_t line_vertex,
        const index_t boundary_corner_id);

    template class scar_api NearbyLineChunksFromPoint<3>;
    template double compute_line_line_symmetric_hausdorff_distance(
        const RINGMesh::GeoModel3D &,
        const index_t,
        const index_t);
    template double scar_api segment_point_normalized_coordinates(
        const RINGMesh::Geometry::Segment3D &,
        const vec3 &);
    template double scar_api compute_curvilinear_distance(
        const RINGMesh::Line3D &line,
        const index_t line_vertex,
        const index_t boundary_corner_id);
}
