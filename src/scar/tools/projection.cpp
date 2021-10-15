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

#include <scar/tools/projection.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/mesh_index.h>

#include <scar/tools/distance.h>
#include <scar/tools/geometry.h>
#include <scar/tools/utils.h>

namespace
{
    using namespace SCAR;

    index_t is_point_on_polygon_vertex(
        const RINGMesh::Surface3D &surface,
        const index_t polygon,
        const vec3 &point)
    {
        for (index_t fv = 0; fv < surface.nb_mesh_element_vertices(polygon);
             ++fv)
        {
            if ((point - surface.mesh_element_vertex({polygon, fv})).length2() < surface.geomodel().epsilon2())
            {
                return fv;
            }
        }
        return NO_ID;
    }

    /*!
     * @warn Only on triangles
     */
    index_t is_point_on_polygon_edge(
        const RINGMesh::Surface3D &surface,
        const index_t polygon,
        const vec3 &point)
    {
        std::array<double, 3> bary_coord;
        bool success;
        std::tie(success, bary_coord) = RINGMesh::triangle_barycentric_coordinates(
            point, surface.mesh_element_vertex({polygon, 0}),
            surface.mesh_element_vertex({polygon, 1}),
            surface.mesh_element_vertex({polygon, 2}));
        scar_assert(success);
        for (index_t fe = 0; fe < surface.nb_mesh_element_vertices(polygon);
             ++fe)
        {
            if (bary_coord[fe] < surface.geomodel().epsilon())
            {
                return (fe + 1) % surface.nb_mesh_element_vertices(polygon);
            }
        }
        return NO_ID;
    }
}

namespace SCAR
{

    template <index_t DIMENSION>
    PointOnLineProjection<DIMENSION> PointOnLineProjectionTool<DIMENSION>::find_projection(
        const vecn<DIMENSION> &point,
        const index_t line_id)
    {
        PointOnLineProjection<DIMENSION> projection;
        std::tie(std::ignore, projection.coord, projection.distance) =
            geomodel_.line(line_id).edge_aabb().closest_edge(point);
        projection.rcoord =
            PointOnLineCoordinates<DIMENSION>::relative_coordinates(
                projection.coord, geomodel_.line(line_id).mesh());
        return projection;
    }

    void LineOnSurfaceProjectionTool::find_projection(
        const index_t line_id,
        const index_t surface_id,
        LineOnSurfaceProjection &line_projection,
        bool filter_line_projection_vertices)
    {
        line_projection.vertices.clear();
        std::vector<vec3> line_vertices_projection;
        double line_vertex_distance = find_line_vertices_projection(line_id,
                                                                    surface_id, line_vertices_projection);
        std::vector<vec3> complete_line_projection;
        std::vector<bool> line_projection_vertices_inside_polygon;
        double line_segment_distance = stick_line_segments_to_surface_triangulation(
            line_vertices_projection, surface_id, complete_line_projection,
            line_projection_vertices_inside_polygon);
        line_projection.max_distance = std::max(line_vertex_distance,
                                                line_segment_distance);
        // Filter line projection to remove vertices with straight angles (180 degrees)
        if (filter_line_projection_vertices)
        {
            filter_vertices(complete_line_projection,
                            line_projection_vertices_inside_polygon, line_projection.vertices);
        }
        else
        {
            line_projection.vertices = complete_line_projection;
        }
    }

    double LineOnSurfaceProjectionTool::find_line_vertices_projection(
        const index_t line_id,
        const index_t surface_id,
        std::vector<vec3> &line_vertices_projection)
    {
        double max_distance = 0;
        const RINGMesh::Line3D &line = geomodel_.line(line_id);
        const RINGMesh::Surface3D &surface = geomodel_.surface(surface_id);
        line_vertices_projection.clear();
        line_vertices_projection.resize(line.nb_vertices(), NO_POINT_3D);

        for (index_t lv = 0; lv < line.nb_vertices(); ++lv)
        {
            double cur_distance;
            std::tie(std::ignore, line_vertices_projection[lv], cur_distance) =
                surface.polygon_aabb().closest_triangle(line.vertex(lv));
            max_distance = std::max(max_distance, cur_distance);
        }
        return max_distance;
    }

    double LineOnSurfaceProjectionTool::stick_line_segments_to_surface_triangulation(
        const std::vector<vec3> &line_vertices_projection,
        const index_t surface_id,
        std::vector<vec3> &cohesive_line_vertices,
        std::vector<bool> &line_vertices_inside_polygon)
    {
        cohesive_line_vertices.clear();
        line_vertices_inside_polygon.clear();
        double max_segment_distance = 0;
        for (index_t s = 0; s < line_vertices_projection.size() - 1; ++s)
        {
            vec3 proj_relocation;
            line_vertices_inside_polygon.push_back(
                is_point_inside_surface_polygon(surface_id,
                                                line_vertices_projection[s], proj_relocation));
            cohesive_line_vertices.push_back(proj_relocation);

            double cur_segment_distance =
                stick_line_segment_to_surface_triangulation(
                    line_vertices_projection[s], line_vertices_projection[s + 1],
                    surface_id, cohesive_line_vertices,
                    line_vertices_inside_polygon);
            max_segment_distance = std::max(max_segment_distance,
                                            cur_segment_distance);
        }
        vec3 proj_relocation;
        line_vertices_inside_polygon.push_back(
            is_point_inside_surface_polygon(surface_id,
                                            line_vertices_projection.back(), proj_relocation));
        cohesive_line_vertices.push_back(proj_relocation);

        return max_segment_distance;
    }

    double LineOnSurfaceProjectionTool::stick_line_segment_to_surface_triangulation(
        const vec3 &segment_v0,
        const vec3 &segment_v1,
        const index_t surface_id,
        std::vector<vec3> &cohesive_line_vertices,
        std::vector<bool> &line_vertices_inside_polygon)
    {
        std::vector<index_t> v0_containing_polygons;
        get_containing_polygons(surface_id, segment_v0, v0_containing_polygons);
        std::vector<index_t> v1_containing_polygons;
        get_containing_polygons(surface_id, segment_v1, v1_containing_polygons);
        bool different_polygons = are_all_elements_not_in(v1_containing_polygons,
                                                          v0_containing_polygons);
        scar_assert(
            !v0_containing_polygons.empty() && !v1_containing_polygons.empty());
        if (!different_polygons)
        {
            // This segment is included in a surface polygon, the distance between
            // the segment and the surface is null.
            return 0;
        }
        else
        {
            const RINGMesh::Surface3D &surface = geomodel_.surface(surface_id);
            for (index_t f = 0; f < v0_containing_polygons.size(); ++f)
            {
                // For each polygon containing segment first point
                index_t cur_v0_polygon = v0_containing_polygons[f];
                index_t nb_polygon_edges = surface.nb_mesh_element_vertices(
                    cur_v0_polygon);
                for (index_t polygon_edge = 0; polygon_edge < nb_polygon_edges;
                     ++polygon_edge)
                {
                    double cur_distance_line_segment_to_surface_edge;
                    vec3 proj_on_line_segment, proj_on_polygon_edge;
                    bool parallel;
                    bool orthogonal_distance =
                        segment_segment_distance(segment_v0, segment_v1,
                                                 surface.mesh_element_vertex({cur_v0_polygon,
                                                                              polygon_edge}),
                                                 surface.mesh_element_vertex(
                                                     {cur_v0_polygon, (polygon_edge + 1) % nb_polygon_edges}),
                                                 cur_distance_line_segment_to_surface_edge,
                                                 proj_on_line_segment, proj_on_polygon_edge, parallel);
                    // Checking projection on line segment
                    if (!orthogonal_distance)
                    {
                        continue;
                    }
                    double normalized_distance_on_line_segment =
                        (proj_on_line_segment - segment_v0).length() / (segment_v1 - segment_v0).length();
                    std::vector<index_t> debug_vector;
                    get_containing_polygons(surface_id, proj_on_polygon_edge,
                                            debug_vector);
                    if (normalized_distance_on_line_segment > global_epsilon && (1 - normalized_distance_on_line_segment) > global_epsilon)
                    {
                        cohesive_line_vertices.push_back(proj_on_polygon_edge);
                        line_vertices_inside_polygon.push_back(false);
                        // Recursive call on the segment from the found projection on
                        // polygon edge to the end of the segment.
                        double max_distance_line_segment_to_surface =
                            stick_line_segment_to_surface_triangulation(
                                proj_on_polygon_edge, segment_v1, surface_id,
                                cohesive_line_vertices,
                                line_vertices_inside_polygon);
                        return std::max(cur_distance_line_segment_to_surface_edge,
                                        max_distance_line_segment_to_surface);
                    }
                }
            }
            // If following is reached, all v0 containing polygons have been processed
            // and none projection have been found, which is unexpected.
            scar_assert_not_reached;
            return -1;
        }
    }

    bool LineOnSurfaceProjectionTool::is_point_inside_surface_polygon(
        const index_t surface_id,
        const vec3 &point,
        vec3 &relocation)
    {
        const RINGMesh::Surface3D &surface = geomodel_.surface(surface_id);
        vec3 proj;
        double distance;
        index_t containing_polygon_id;
        std::tie(containing_polygon_id, proj, distance) =
            surface.polygon_aabb().closest_triangle(point);
        scar_assert(distance < geomodel_.epsilon());
        index_t close_vertex_id = is_point_on_polygon_vertex(surface,
                                                             containing_polygon_id, point);
        if (close_vertex_id != NO_ID)
        {
            relocation = surface.mesh_element_vertex({containing_polygon_id,
                                                      close_vertex_id});
            scar_assert((point - relocation).length2() < geomodel_.epsilon2());
            return false;
        }
        index_t close_edge = is_point_on_polygon_edge(surface,
                                                      containing_polygon_id, point);
        if (close_edge != NO_ID)
        {
            std::array<double, 3> bary_coord;
            bool success;
            std::tie(success, bary_coord) =
                RINGMesh::triangle_barycentric_coordinates(point,
                                                           surface.mesh_element_vertex({containing_polygon_id, 0}),
                                                           surface.mesh_element_vertex({containing_polygon_id, 1}),
                                                           surface.mesh_element_vertex({containing_polygon_id, 2}));
            scar_assert(success);
            index_t null_bary_coord_vertex = (close_edge + 2) % surface.nb_mesh_element_vertices(containing_polygon_id);
            double bary_coord_close_edge = bary_coord[null_bary_coord_vertex];
            for (index_t e = 0; e < 3; ++e)
            {
                if (e != null_bary_coord_vertex)
                {
                    bary_coord[e] += 0.5 * bary_coord_close_edge;
                }
                else
                {
                    bary_coord[e] = 0;
                }
            }
            relocation = bary_coord[0] * surface.mesh_element_vertex({containing_polygon_id, 0}) + bary_coord[1] * surface.mesh_element_vertex({containing_polygon_id, 1}) + bary_coord[2] * surface.mesh_element_vertex({containing_polygon_id, 2});
            //@todo Relocation could be greater then epsilon. What to do?
            return false;
        }
        relocation = point;
        return true;
    }

    void LineOnSurfaceProjectionTool::get_containing_polygons(
        const index_t surface_id,
        const vec3 &point,
        std::vector<index_t> &containing_polygons)
    {
        containing_polygons.clear();
        const RINGMesh::Surface3D &surface = geomodel_.surface(surface_id);
        const auto &surface_mesh = surface.mesh();
        vec3 proj;
        double distance;
        index_t first_containing_polygon;
        std::tie(first_containing_polygon, proj, distance) =
            surface.polygon_aabb().closest_triangle(point);
        scar_assert(distance < geomodel_.epsilon());
        index_t local_vertex_id = is_point_on_polygon_vertex(surface,
                                                             first_containing_polygon, point);
        if (local_vertex_id != NO_ID)
        {
            containing_polygons = surface_mesh.polygons_around_vertex(
                surface.mesh_element_vertex_index({first_containing_polygon,
                                                   local_vertex_id}),
                false,
                first_containing_polygon);
            return;
        }
        containing_polygons.push_back(first_containing_polygon);
        index_t polygon_edge_id = is_point_on_polygon_edge(surface,
                                                           first_containing_polygon, point);
        if (polygon_edge_id != NO_ID)
        {
            containing_polygons.push_back(surface_mesh.polygon_adjacent({first_containing_polygon, polygon_edge_id}));
        }
    }

    void LineOnSurfaceProjectionTool::filter_vertices(
        const std::vector<vec3> &line_vertices,
        const std::vector<bool> &line_vertex_inside_surface_polygon,
        std::vector<vec3> &filtered_line_vertices)
    {
        scar_assert(
            line_vertices.size() == line_vertex_inside_surface_polygon.size());

        std::vector<bool> vertex_to_delete(line_vertices.size(), false);
        for (index_t v = 1; v < line_vertices.size() - 1; ++v)
        {
            if (!line_vertex_inside_surface_polygon[v])
            {
                continue;
            }
            double diff_lengthes =
                (line_vertices[v - 1] - line_vertices[v + 1]).length() - ((line_vertices[v] - line_vertices[v - 1]).length() + (line_vertices[v] - line_vertices[v + 1]).length());
            double triar = triangle_aspect_ratio(
                RINGMesh::Geometry::Triangle3D(line_vertices[v - 1],
                                               line_vertices[v], line_vertices[v - 1]));
            if (fabs(diff_lengthes) < geomodel_.epsilon() || triar > 1e100)
            {
                vertex_to_delete[v] = true;
            }
        }
        for (index_t v = 0; v < line_vertices.size(); ++v)
        {
            if (!vertex_to_delete[v])
            {
                filtered_line_vertices.push_back(line_vertices[v]);
            }
        }
    }

    template class scar_api PointOnLineProjectionTool<2>;

    template class scar_api PointOnLineProjectionTool<3>;
}
