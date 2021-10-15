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

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <scar/tools/distance.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS(Line);
    ALIAS_2D_AND_3D(Line);
}

namespace SCAR
{

    double get_signed_distance_point_to_plane(
        const vec3 &p,
        const vec3 &plane_norm,
        const vec3 &plane_origin);

    /*!
     * Computes the orthogonal projection of a point on a line
     * @param[in] p the point to project
     * @param[in] line_origin a point on a line
     * @param[in] line_dir Parallel vector of the line
     * @return the projected point
     */
    template <index_t DIMENSION>
    vecn<DIMENSION> get_point_line_projection(
        const vecn<DIMENSION> &p,
        const RINGMesh::Geometry::Line<DIMENSION> line);

    std::tuple<index_t, index_t, double> find_closest_point_pair(
        const std::vector<vec2> &points);

    std::tuple<index_t, index_t, double> find_closest_point_pair(
        const std::vector<vec3> &points);

    template <index_t DIMENSION>
    std::tuple<index_t, std::vector<vecn<DIMENSION>>> intersections_sphere_line(
        const RINGMesh::Line<DIMENSION> &line,
        const RINGMesh::Geometry::Sphere<DIMENSION> &sphere);

    template <index_t DIMENSION>
    double triangle_aspect_ratio(
        const RINGMesh::Geometry::Triangle<DIMENSION> &triangle);

    template <index_t DIMENSION>
    std::tuple<bool, std::vector<vecn<DIMENSION>>> intersection_line_line(
        const RINGMesh::Line<DIMENSION> &line1,
        const RINGMesh::Line<DIMENSION> &line2);

    // Copy of geogram

    /**
     * \brief Computes the cosine of the angle between two 3d vectors.
     * \param[in] v0 first vector
     * \param[in] v1 second vector
     * \return the cosine of the angle between \p v0 and \p v1
     */
    template <index_t DIMENSION>
    double cos_angle(const vecn<DIMENSION> &v0, const vecn<DIMENSION> &v1)
    {
        double nv0_2{length2(v0)};
        double nv1_2{length2(v1)};
        scar_assert(length2(v0) > global_epsilon);
        scar_assert(length2(v1) > global_epsilon);
        double result{dot(v0, v1) / std::sqrt(nv0_2 * nv1_2)};
        if (result > 1 && result - 1 < global_epsilon)
        {
            return 1;
        }
        else if (result < -1 && -1 - result < global_epsilon)
        {
            return -1;
        }
        else
        {
            return result;
        }
    }

    /**
     * \brief Computes the angle between two 3d vectors.
     * \param[in] v0 first vector
     * \param[in] v1 second vector
     * \return the angle between \p v0 and \p v1 in radians, in
     *  the interval \f$ [ 0 \ldots \pi ] \f$.
     */
    template <index_t DIMENSION>
    double angle(const vecn<DIMENSION> &v0, const vecn<DIMENSION> &v1)
    {
        return acos(cos_angle(v0, v1));
    }

    // In degrees
    template <index_t DIMENSION>
    double compute_angle_between_two_adjacent_lines(
        const RINGMesh::Corner<DIMENSION> &corner,
        const RINGMesh::Line<DIMENSION> &line1,
        const RINGMesh::Line<DIMENSION> &line2)
    {
        vecn<DIMENSION> l1_vertex =
            line1.boundary(0).gmme() == corner.gmme() ? line1.vertex(1) : line1.vertex(line1.nb_vertices() - 2);
        vecn<DIMENSION> l2_vertex =
            line2.boundary(0).gmme() == corner.gmme() ? line2.vertex(1) : line2.vertex(line2.nb_vertices() - 2);
        vecn<DIMENSION> corner_vertex = corner.vertex(0);
        double angle_in_radians = angle(l1_vertex - corner_vertex,
                                        l2_vertex - corner_vertex);
        return angle_in_radians * 180. / M_PI;
    }

}
