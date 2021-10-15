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

#include <scar/tools/geometry.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/mesh/mesh_index.h>

namespace SCAR
{

    double get_signed_distance_point_to_plane(
        const vec3 &p,
        const vec3 &plane_norm,
        const vec3 &plane_origin)
    {
        vec3 vec_OP{p - plane_origin};
        return dot(vec_OP, normalize(plane_norm));
    }

    template <index_t DIMENSION>
    vecn<DIMENSION> get_point_line_projection(
        const vecn<DIMENSION> &p,
        const RINGMesh::Geometry::Line<DIMENSION> line)
    {
        vecn<DIMENSION> diff{p - line.origin};
        double d{dot(normalize(line.direction), diff)};
        return (line.origin + d * normalize(line.direction));
    }

    template <index_t DIMENSION>
    std::tuple<index_t, index_t, double> find_closest_point_pair(
        const std::vector<vecn<DIMENSION>> &points)
    {
        double shortest_dist{std::numeric_limits<double>::max()};
        index_t first_pt{NO_ID};
        index_t second_pt{NO_ID};
        for (index_t pt0 = 0; pt0 < points.size() - 1; pt0++)
        {
            for (index_t pt1 = pt0 + 1; pt1 < points.size(); pt1++)
            {
                double cur_dist{(points[pt0] - points[pt1]).length()};
                if (cur_dist < shortest_dist)
                {
                    shortest_dist = cur_dist;
                    first_pt = pt0;
                    second_pt = pt1;
                }
            }
        }
        return std::make_tuple(first_pt, second_pt, shortest_dist);
    }

    std::tuple<index_t, index_t, double> find_closest_point_pair(
        const std::vector<vec2> &points)
    {
        return find_closest_point_pair<2>(points);
    }

    std::tuple<index_t, index_t, double> find_closest_point_pair(
        const std::vector<vec3> &points)
    {
        return find_closest_point_pair<3>(points);
    }

    template <index_t DIMENSION>
    std::tuple<index_t, std::vector<vecn<DIMENSION>>> intersections_sphere_line(
        const RINGMesh::Line<DIMENSION> &line,
        const RINGMesh::Geometry::Sphere<DIMENSION> &sphere)
    {
        std::vector<vecn<DIMENSION>> intersections;
        //@todo use an aabb tree
        for (auto cur_edge : RINGMesh::range(line.nb_mesh_elements()))
        {
            std::vector<vecn<DIMENSION>> edge_intersections;
            RINGMesh::Geometry::Segment<DIMENSION> edge{line.mesh_element_vertex(
                                                            {cur_edge, 0}),
                                                        line.mesh_element_vertex(
                                                            {cur_edge, 1})};
            std::tie(std::ignore, edge_intersections) =
                RINGMesh::Intersection::segment_sphere(edge, sphere);
            intersections.insert(intersections.end(), edge_intersections.begin(),
                                 edge_intersections.end());
        }
        return std::make_tuple(static_cast<index_t>(intersections.size()),
                               intersections);
    }

    template <index_t DIMENSION>
    double triangle_aspect_ratio(
        const RINGMesh::Geometry::Triangle<DIMENSION> &triangle)
    {
        double l01{(triangle.p2 - triangle.p0).length()};
        double l12{(triangle.p2 - triangle.p1).length()};
        double l20{(triangle.p0 - triangle.p2).length()};

        double num{l01 * l12 * l20};
        double denom{(l01 + l12 - l20) * (l01 - l12 + l20) * (-l01 + l12 + l20)};

        if (denom < global_epsilon)
        {
            return std::numeric_limits<double>::max();
        }

        return num / denom;
    }

    template <>
    std::tuple<bool, std::vector<vec2>> intersection_line_line<2>(
        const RINGMesh::Line2D &line1,
        const RINGMesh::Line2D &line2)
    {
        std::vector<vec2> intersection_points;
        for (auto e1 : RINGMesh::range(line1.nb_mesh_elements()))
        {
            RINGMesh::Box2D e1_box;
            auto e1_v0 = line1.mesh_element_vertex({e1, 0});
            auto e1_v1 = line1.mesh_element_vertex({e1, 1});
            e1_box.add_point(e1_v0);
            e1_box.add_point(e1_v1);
            for (auto e2 : RINGMesh::range(line2.nb_mesh_elements()))
            {
                RINGMesh::Box2D e2_box;
                auto e2_v0 = line2.mesh_element_vertex({e2, 0});
                auto e2_v1 = line2.mesh_element_vertex({e2, 1});
                e2_box.add_point(e2_v0);
                e2_box.add_point(e2_v1);
                if (!e1_box.bboxes_overlap(e2_box))
                {
                    continue;
                }
                bool intersect;
                vec2 pt_intersection;
                std::tie(intersect, pt_intersection) =
                    RINGMesh::Intersection::segment_segment(
                        RINGMesh::Geometry::Segment2D(e1_v0, e1_v1),
                        RINGMesh::Geometry::Segment2D(e2_v0, e2_v1));
                if (intersect)
                {
                    intersection_points.push_back(pt_intersection);
                }
            }
        }

        return std::make_tuple(!intersection_points.empty(), intersection_points);
    }

    template <>
    std::tuple<bool, std::vector<vec3>> intersection_line_line<3>(
        const RINGMesh::Line3D &line1,
        const RINGMesh::Line3D &line2)
    {
        scar_unused(line1);
        scar_unused(line2);
        std::vector<vec3> intersection_points;
        throw SCARException("3D", "Line-Line intersection KO in 3D");
        return std::make_tuple(false, intersection_points);
    }

    //    template< index_t DIMENSION >
    //    std::tuple< bool, std::vector< vecn< DIMENSION > > > line_line_intersection(
    //        const RINGMesh::Line< DIMENSION >& line1,
    //        const RINGMesh::Line< DIMENSION >& line2 );

    template vec2 scar_api get_point_line_projection(
        const vec2 &,
        const RINGMesh::Geometry::Line2D);
    template double scar_api triangle_aspect_ratio(
        const RINGMesh::Geometry::Triangle2D &);

    template std::tuple<index_t, std::vector<vecn<2>>> intersections_sphere_line(
        const RINGMesh::Line<2> &,
        const RINGMesh::Geometry::Sphere<2> &);

    template vec3 scar_api get_point_line_projection(
        const vec3 &,
        const RINGMesh::Geometry::Line3D);
    template double scar_api triangle_aspect_ratio(
        const RINGMesh::Geometry::Triangle3D &);

    template std::tuple<index_t, std::vector<vecn<3>>> intersections_sphere_line(
        const RINGMesh::Line<3> &,
        const RINGMesh::Geometry::Sphere<3> &);
}
