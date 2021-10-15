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
#include <scar/tools/distance.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModel);
    ALIAS_2D_AND_3D(GeoModel);
}

namespace SCAR
{

    template <index_t DIMENSION>
    struct scar_api PointOnLineProjection
    {
        vecn<DIMENSION> coord{};
        double distance{-1};
        LineRelativeCoordinates rcoord{};
    };

    struct scar_api LineOnSurfaceProjection
    {
        LineOnSurfaceProjection()
            : vertices(), max_distance(-1)
        {
        }
        std::vector<vec3> vertices{};
        double max_distance{};
    };

    template <index_t DIMENSION>
    class scar_api PointOnLineProjectionTool
    {
    public:
        PointOnLineProjectionTool(const RINGMesh::GeoModel<DIMENSION> &geomodel)
            : geomodel_(geomodel)
        {
        }

        PointOnLineProjection<DIMENSION> find_projection(
            const vecn<DIMENSION> &point,
            const index_t line_id);

    private:
        const RINGMesh::GeoModel<DIMENSION> &geomodel_;
    };

    class scar_api LineOnSurfaceProjectionTool
    {
    public:
        LineOnSurfaceProjectionTool(const RINGMesh::GeoModel3D &geomodel)
            : geomodel_(geomodel)
        {
        }

        /*!
         * @brief Computes the projection of a line on a surface
         * @details Orthogonal (to the surface) projection
         * @param[in] line_id Index of the line which is projected
         * @param[in] surface_id Index of the surface on which is
         * projected \p line_id
         * @param[out] projection_line_vertices Vertices of the line
         * projection on the surface
         */
        void find_projection(
            const index_t line_id,
            const index_t surface_id,
            LineOnSurfaceProjection &line_projection,
            bool filter_line_projection_vertices = true);

    private:
        double find_line_vertices_projection(
            const index_t line_id,
            const index_t surface_id,
            std::vector<vec3> &line_vertices_projection);

        /*!
         * @brief Splits given line segments to glue the line discretization
         * to be cohesive with surface triangulation.
         * @param[in] line_vertices_projection Not cohesize line projection vertices
         * @param[in] surface_id Index of the surface on which is
         * projected the line
         * @param[out] cohesive_line_vertices Vertices of the cohesive line
         * projection on the surface
         * @param[out] line_vertices_inside_polygon True if the i-th line projection
         * vertex is inside a surface polygon, false if it is on an edge or a surface
         * vertex
         * @return the distance
         */
        double stick_line_segments_to_surface_triangulation(
            const std::vector<vec3> &line_vertices_projection,
            const index_t surface_id,
            std::vector<vec3> &cohesive_line_vertices,
            std::vector<bool> &line_vertices_inside_polygon);

        double stick_line_segment_to_surface_triangulation(
            const vec3 &segment_v0,
            const vec3 &segment_v1,
            const index_t surface_id,
            std::vector<vec3> &cohesive_line_vertices,
            std::vector<bool> &line_vertices_inside_polygon);

        bool is_point_inside_surface_polygon(
            const index_t surface_id,
            const vec3 &point,
            vec3 &relocation);

        void get_containing_polygons(
            const index_t surface_id,
            const vec3 &point,
            std::vector<index_t> &containing_polygons);

        void filter_vertices(
            const std::vector<vec3> &line_vertices,
            const std::vector<bool> &line_vertex_on_surface_edge,
            std::vector<vec3> &filtered_line_vertices);

    private:
        const RINGMesh::GeoModel3D &geomodel_;
    };

}
