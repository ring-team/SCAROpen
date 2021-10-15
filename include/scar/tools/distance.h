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

#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/basic/geometry.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModel);
    ALIAS_2D_AND_3D(GeoModel);
    FORWARD_DECLARATION_DIMENSION_CLASS(Line);
}

namespace SCAR
{
    FORWARD_DECLARATION_DIMENSION_CLASS(ConvexShape);
}

namespace SCAR
{

    template <index_t DIMENSION>
    class scar_api DistancePointToLine
    {
    public:
        DistancePointToLine(
            const index_t line_index,
            const vecn<DIMENSION> &proj,
            const double dist)
            : line_id(line_index), projection(proj), distance(dist)
        {
        }

        bool operator<(const DistancePointToLine &rhs) const
        {
            return this->distance < rhs.distance;
        }

        index_t line_id{};
        vecn<DIMENSION> projection{};
        double distance{};
    };

    struct scar_api LineRelativeCoordinates
    {
        LineRelativeCoordinates() = default;
        LineRelativeCoordinates(index_t edge, double edge_norm)
            : edge(edge), norm(edge_norm)
        {
            scar_assert(edge_norm >= 0 && edge_norm <= 1);
        }

        LineRelativeCoordinates(double position)
        {
            double int_part;
            norm = std::modf(position, &int_part);
            edge = static_cast<index_t>(int_part);
            scar_assert(norm >= 0 && norm <= 1);
        }

        operator double() const
        {
            return norm + edge;
        }

        index_t edge{};
        double norm{};
    };

    template <index_t DIMENSION>
    struct scar_api PointOnLineCoordinates
    {
        PointOnLineCoordinates(
            const vecn<DIMENSION> &point,
            const LineRelativeCoordinates &relative_coord)
            : absolute_coord(point), relative_coord(relative_coord)
        {
        }

        bool operator<(const PointOnLineCoordinates &rhs) const
        {
            return this->relative_coord < rhs.relative_coord;
        }

        static vecn<DIMENSION> absolute_coordinates(
            const LineRelativeCoordinates &relative_coord,
            const RINGMesh::LineMesh<DIMENSION> &line_mesh)
        {
            if (relative_coord.norm == 0)
            {
                return line_mesh.vertex(relative_coord.edge);
            }
            else
            {
                return (1 - relative_coord.norm) * line_mesh.vertex(relative_coord.edge) + relative_coord.norm * line_mesh.vertex(relative_coord.edge + 1);
            }
        }

        static LineRelativeCoordinates relative_coordinates(
            const vecn<DIMENSION> &absolute_coord,
            const RINGMesh::LineMesh<DIMENSION> &line_mesh,
            index_t edge_id = NO_ID)
        {
            if (edge_id == NO_ID)
            {
                double distance;
                std::tie(edge_id, std::ignore, distance) =
                    line_mesh.edge_aabb().closest_edge(absolute_coord);
                if (distance > RINGMesh::global_epsilon)
                {
                    throw SCARException("POLCoord",
                                        "Query point is not on the line");
                }
            }
            scar_assert(edge_id != NO_ID);
            auto edge_v0 = line_mesh.vertex(
                line_mesh.edge_vertex({edge_id, 0}));
            auto edge_v1 = line_mesh.vertex(
                line_mesh.edge_vertex({edge_id, 1}));
            double distance_p_v0{(edge_v0 - absolute_coord).length()};
            double distance_p_v1{(edge_v1 - absolute_coord).length()};
            double edge_length{(edge_v1 - edge_v0).length()};
            if (std::fabs(distance_p_v0 + distance_p_v1 - edge_length) > RINGMesh::global_epsilon)
            {
                throw SCARException("POLCoord", "Query point in not on line edge");
            }
            scar_assert(edge_length > RINGMesh::global_epsilon);
            double norm{distance_p_v0 / edge_length};
            // Small hack to prevent norm to be greater than 1.
            if (norm > 1 && norm - 1 < RINGMesh::global_epsilon)
            {
                norm = 1;
            }

            return LineRelativeCoordinates(edge_id, norm);
        }

        vecn<DIMENSION> absolute_coord{};
        LineRelativeCoordinates relative_coord{};
    };
    ALIAS_2D_AND_3D(PointOnLineCoordinates);

    struct scar_api LineChunk
    {
        LineChunk() = default;
        LineChunk(
            const LineRelativeCoordinates &from,
            const LineRelativeCoordinates &to)
            : begin(from), end(to)
        {
        }

        std::tuple<bool, LineChunk> intersect(const LineChunk &rhs) const
        {
            LineChunk result = rhs;
            result.begin = std::max(begin, rhs.begin);
            result.end = std::min(end, rhs.end);
            if (result.end > result.begin)
            {
                return std::make_tuple(true, result);
            }
            return std::make_tuple(false, LineChunk());
        }

        RINGMesh::range get_range() const
        {
            index_t first_vertex = begin.edge;
            if (begin.norm == 1.)
            {
                ++first_vertex;
            }
            index_t last_vertex = end.edge;
            if (end.norm == 1.)
            {
                ++last_vertex;
            }
            // +1 is because end is the last vertex we want to process
            return RINGMesh::range(first_vertex, last_vertex + 1);
        }

        RINGMesh::range get_internal_range() const
        {
            index_t first_vertex = static_cast<index_t>(std::ceil(begin));
            index_t last_vertex = static_cast<index_t>(std::floor(end));
            // +1 is because end is the last vertex we want to process
            return RINGMesh::range(first_vertex, last_vertex + 1);
        }

        bool fully_include(const LineChunk &rhs) const
        {
            if (rhs.begin < begin)
            {
                return false;
            }
            if (rhs.end > end)
            {
                return false;
            }
            return true;
        }

        bool is_almost_null(double epsilon) const
        {
            return std::fabs(end - begin) < epsilon;
        }

        bool is_null() const
        {
            return begin == end;
        }

        bool is_inside(const LineRelativeCoordinates &query) const
        {
            return query >= begin && query <= end;
        }

        bool is_before(const LineRelativeCoordinates &query) const
        {
            return query < begin;
        }

        bool is_after(const LineRelativeCoordinates &query) const
        {
            return query > end;
        }

        bool operator!=(const LineChunk &rhs) const
        {
            if (rhs.begin != this->begin)
            {
                return true;
            }
            if (rhs.end != this->end)
            {
                return true;
            }
            return false;
        }

        bool operator==(const LineChunk &rhs) const
        {
            return !operator!=(rhs);
        }

        LineRelativeCoordinates begin{};
        LineRelativeCoordinates end{};
    };

    inline std::ostream &operator<<(std::ostream &out, const LineChunk &chunk)
    {
        out << "[" << chunk.begin << " - ";
        out << chunk.end << "]";
        return out;
    }

    template <index_t DIMENSION>
    struct scar_api Segment
    {
        Segment(const vecn<DIMENSION> &point0, const vecn<DIMENSION> &point1)
            : p0(point0), p1(point1)
        {
        }
        vecn<DIMENSION> p0{};
        vecn<DIMENSION> p1{};
    };

    ALIAS_2D_AND_3D(Segment);

    using LineChunks = std::vector<LineChunk>;
    void show_chunks(const LineChunks &);

    template <index_t DIMENSION>
    class scar_api NearbyLineChunksFromPoint
    {
    public:
        NearbyLineChunksFromPoint(
            const ConvexShape<DIMENSION> &tolerance_shape,
            const RINGMesh::LineMesh<DIMENSION> &line_mesh);

        index_t nb_chunks() const
        {
            return static_cast<index_t>(chunks_.size());
        }
        LineChunks chunks() const
        {
            return chunks_;
        }

    private:
        LineChunks chunks_{};
    };

    ALIAS_2D_AND_3D(NearbyLineChunksFromPoint);

    class DistanceLineToLine
    {
    public:
        DistanceLineToLine(const index_t line_index, const double dist)
            : line_id(line_index), distance(dist)
        {
        }

        DistanceLineToLine(const index_t line_index)
            : line_id(line_index), distance(-1)
        {
        }

        bool operator<(const DistanceLineToLine &rhs) const
        {
            return this->distance < rhs.distance;
        }

        index_t line_id{};
        double distance{};
    };

    class MinDistanceLineToLine : public DistanceLineToLine
    {
    public:
        MinDistanceLineToLine(
            const index_t line_index,
            const double dist,
            const vec3 &nearest_pt)
            : DistanceLineToLine(line_index, dist), nearest_point(nearest_pt)
        {
        }

        MinDistanceLineToLine(const index_t line_index)
            : DistanceLineToLine(line_index), nearest_point(NO_POINT_3D)
        {
        }

        vec3 nearest_point{};
    };

    template <index_t DIMENSION>
    double segment_point_normalized_coordinates(
        const RINGMesh::Geometry::Segment<DIMENSION> &oriented_segment,
        const vecn<DIMENSION> &point);

    template <index_t DIMENSION>
    double scar_api compute_line_line_symmetric_hausdorff_distance(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        const index_t line1_id,
        const index_t line2_id);

    double scar_api compute_surface_surface_symmetric_hausdorff_distance(
        const RINGMesh::GeoModel3D &geomodel,
        const index_t surface1_id,
        const index_t surface2_id);

    /*
     * boundary_corner_id should be 0 or 1
     */
    template <index_t DIMENSION>
    double scar_api compute_curvilinear_distance(
        const RINGMesh::Line<DIMENSION> &line,
        const index_t line_vertex,
        const index_t boundary_corner_id);

    double scar_api compute_line_vertices_surface_distance(
        const RINGMesh::GeoModel3D &geomodel,
        const index_t line_id,
        const index_t surface_id);

    void scar_api compute_geodesic_to_borders_djikstra_approximation(
        const RINGMesh::GeoModel3D &geomodel,
        const index_t surface_id);

    /*!
     * Calculates the minimal distance between 2 lines
     * @param[in] l0_pt0 the first vertex of the first line
     * @param[in] l0_pt1 the second vertex of the first line
     * @param[in] l1_pt0 the first vertex of the second line
     * @param[in] l1_pt1 the second vertex of the second line
     * @param[out] midpoint the center of the shortest segment
     * between the two lines
     * @return distance between two lines
     * @warning round operation errors
     * @todo Needs to be tested
     */
    bool scar_api line_line_minimum_distance(
        const vec3 &l0_pt0,
        const vec3 &l0_pt1,
        const vec3 &l1_pt0,
        const vec3 &l1_pt1,
        double &distance,
        vec3 &midpoint);

    /*!
     * Calculates the minimal distance between two segments
     * @param[in] s0_pt0 the first vertex of the first segment
     * @param[in] s0_pt1 the second vertex of the first segment
     * @param[in] s1_pt0 the first vertex of the second segment
     * @param[in] s1_pt1 the second vertex of the second segment
     * @param[out] distance Minimal distance between segments
     * @param[out] nearest_s0_pt Point on first segment which is the closest
     * to the second segment
     * @param[out] nearest_s1_pt Point on second segment which is the closest
     * to the first segment
     * @param[out] are parallel True if the segments are parallel
     * @return True if the minimum distance between segments is the same than the
     * minimum distance between the infinite lines defined by the segments
     * @warning round operation errors
     * @todo Needs to be tested
     */
    bool scar_api segment_segment_distance(
        const vec3 &s0_pt0,
        const vec3 &s0_pt1,
        const vec3 &s1_pt0,
        const vec3 &s1_pt1,
        double &distance,
        vec3 &nearest_s0_pt,
        vec3 &nearest_s1_pt,
        bool &are_parallel);
}
