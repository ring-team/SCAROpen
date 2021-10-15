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

#include <scar/repair/topology_recovery.h>

#include <stack>

#include <geogram/basic/progress.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/io/io.h>
#include <ringmesh/mesh/mesh_aabb.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>

#include <scar/tools/geometry.h>
#include <scar/tools/projection.h>
#include <scar/tools/utils.h>

#include <ringpcl/correlation_map.h>

namespace
{

    using namespace SCAR;

    index_t nb_active_nodes(const std::vector<NodeInformation> &nodes)
    {
        index_t counter = 0;
        for (const auto &node : nodes)
        {
            if (node.is_active())
            {
                counter++;
            }
        }
        return counter;
    }

    bool corner_is_free_horizon_border(const RINGMesh::Corner2D &corner)
    {
        if (corner.nb_incident_entities() != 1)
        {
            return false;
        }
        if (!corner.incident_entity(0).has_parent())
        {
            return false;
        }
        const auto &parent_interface = corner.incident_entity(0).parent(
            RINGMesh::Interface2D::type_name_static());
        if (parent_interface.geological_feature() == RINGMesh::GeoModelGeologicalEntity2D::GEOL_FEATURE::STRATI)
        {
            return true;
        }
        return false;
    }

    template <index_t DIMENSION>
    double convex_distance(
        const std::unique_ptr<ConvexShape<DIMENSION>> &convex,
        const vecn<DIMENSION> &point)
    {
        auto vector_from_center = point - convex->center();
        auto direction = normalize(vector_from_center);
        auto convex_width = convex->width(direction);
        auto euclidean_distance = vector_from_center.length2();
        scar_assert(convex_width > global_epsilon);
        return euclidean_distance / convex_width;
    }

    template <index_t DIMENSION>
    vecn<DIMENSION> convex_distance_coord(
        const std::unique_ptr<ConvexShape<DIMENSION>> &convex,
        const vecn<DIMENSION> &point)
    {
        auto vector_from_center = point - convex->center();
        auto direction = normalize(vector_from_center);
        auto convex_width = convex->width(direction);
        scar_assert(convex_width > global_epsilon);
        return vector_from_center / convex_width;
    }

    template <index_t DIMENSION>
    vecn<DIMENSION> inverse_convex_distance_coord(
        const std::unique_ptr<ConvexShape<DIMENSION>> &convex,
        const vecn<DIMENSION> &point)
    {
        auto direction = normalize(point);
        auto convex_width = convex->width(direction);
        return convex->center() + point * convex_width;
    }

    template <index_t DIMENSION>
    class ConvexDistanceFunction
    {
    public:
        ConvexDistanceFunction(
            const RINGMesh::LineMesh<DIMENSION> &mesh,
            const std::unique_ptr<ConvexShape<DIMENSION>> &convex)
            : mesh_(mesh), convex_(convex)
        {
        }

        std::tuple<double, vecn<DIMENSION>> operator()(
            const vecn<DIMENSION> &query,
            index_t cur_box) const
        {
            const auto &v0 = convex_distance_coord(convex_,
                                                   mesh_.vertex(mesh_.edge_vertex({cur_box, 0})));
            const auto &v1 = convex_distance_coord(convex_,
                                                   mesh_.vertex(mesh_.edge_vertex({cur_box, 1})));
            double convex_distance;
            vecn<DIMENSION> nearest_point_on_segment_in_convex_space;
            std::tie(convex_distance, nearest_point_on_segment_in_convex_space) =
                RINGMesh::Distance::point_to_segment(
                    convex_distance_coord(convex_, query),
                    RINGMesh::Geometry::Segment<DIMENSION>{v0, v1});
            auto nearest_point_on_segment = inverse_convex_distance_coord(convex_,
                                                                          nearest_point_on_segment_in_convex_space);
            return std::make_tuple(convex_distance, nearest_point_on_segment);
        }

        double operator()(
            const vecn<DIMENSION> &pt1,
            const vecn<DIMENSION> &pt2) const
        {
            const auto &convex_coord_pt1 = convex_distance_coord(convex_, pt1);
            const auto &convex_coord_pt2 = convex_distance_coord(convex_, pt2);
            return length2(convex_coord_pt1 - convex_coord_pt2);
        }

    private:
        const RINGMesh::LineMesh<DIMENSION> &mesh_;
        const std::unique_ptr<ConvexShape<DIMENSION>> &convex_;
    };

    bool corner_is_free_horizon_border(const RINGMesh::Corner3D &corner)
    {
        scar_unused(corner);
        return false;
    }

    template <index_t DIMENSION>
    bool horizon_lines(
        const RINGMesh::Line<DIMENSION> &line1,
        const RINGMesh::Line<DIMENSION> &line2)
    {
        if (line1.parent(RINGMesh::Interface<DIMENSION>::type_name_static()).geological_feature() != RINGMesh::GeoModelGeologicalEntity<DIMENSION>::GEOL_FEATURE::STRATI)
        {
            return false;
        }
        if (line2.parent(RINGMesh::Interface<DIMENSION>::type_name_static()).geological_feature() != RINGMesh::GeoModelGeologicalEntity<DIMENSION>::GEOL_FEATURE::STRATI)
        {
            return false;
        }
        return true;
    }

    vec2 find_oriented_projection_on_line(
        const RINGMesh::Corner2D &corner,
        const RINGMesh::Line2D &line)
    {
        vec2 pt_corner = corner.mesh().vertex(0);
        const auto &incident_line = corner.incident_entity(0);
        scar_assert(line.index() != incident_line.index());
        vec2 pt_line;
        if (incident_line.mesh().vertex(0) == pt_corner)
        {
            pt_line = incident_line.mesh().vertex(1);
        }
        else
        {
            pt_line = incident_line.mesh().vertex(incident_line.nb_vertices() - 2);
        }

        std::vector<vec2> intersections;
        for (auto edge_id : RINGMesh::range(line.nb_mesh_elements()))
        {
            bool intersect;
            vec2 where;
            std::tie(intersect, where) = RINGMesh::Intersection::segment_line(
                RINGMesh::Geometry::Segment2D(line.vertex(edge_id),
                                              line.vertex(edge_id + 1)),
                RINGMesh::Geometry::Line2D(
                    RINGMesh::Geometry::Segment2D(pt_corner, pt_line)));
            if (intersect)
            {
                intersections.push_back(where);
            }
        }
        if (intersections.size() == 1)
        {
            return intersections[0];
        }
        else if (intersections.size() > 1)
        {
            std::sort(intersections.begin(), intersections.end(),
                      [pt_corner](const vec2 &pt1, const vec2 &pt2)
                      { return (pt_corner - pt1).length2() < (pt_corner - pt2).length2(); });
            return intersections[0];
        }
        else
        {
            throw SCARException("Projection", "Not implemented yet");
            return vec2();
        }
    }

    vec3 find_oriented_projection_on_line(
        const RINGMesh::Corner3D &corner,
        const RINGMesh::Line3D &line)
    {
        scar_unused(corner);
        scar_unused(line);
        scar_assert_not_reached;
        return vec3();
    }

    template <index_t DIMENSION>
    vecn<DIMENSION> find_oriented_projection_on_line(
        const RINGMesh::Corner<DIMENSION> &corner,
        const RINGMesh::Line<DIMENSION> &line)
    {
        if (DIMENSION == 2)
        {
            return find_oriented_projection_on_line<2>(corner, line);
        };
        return vecn<DIMENSION>();
    }

}

namespace SCAR
{
    template <index_t DIMENSION>
    GeoModelTopologyRecovererBase<DIMENSION>::GeoModelTopologyRecovererBase(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        ConvexShapeType form_type,
        bool verbose)
        : invalidity_graph(8 * (geomodel.nb_corners() + geomodel.nb_lines()),
                           NO_ID),
          connectivity_graph(8 * (geomodel.nb_corners() + geomodel.nb_lines()),
                             false),
          geomodel_(geomodel),
          shape_type_(form_type),
          verbose_(verbose)
    {
        node_information.reserve(
            10 * (geomodel.nb_corners() + geomodel.nb_lines()));
        edge_information.reserve(
            (geomodel.nb_corners() + geomodel.nb_lines()) * (geomodel.nb_corners() + geomodel.nb_lines()));
    }

    template <index_t DIMENSION>
    GeoModelTopologyRecoverer<DIMENSION>::GeoModelTopologyRecoverer(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        ConvexShapeType form_type,
        bool verbose)
        : GeoModelTopologyRecovererBase<DIMENSION>(geomodel, form_type, verbose)
    {
    }

    GeoModelTopologyRecoverer<2>::GeoModelTopologyRecoverer(
        RINGMesh::GeoModel2D &geomodel,
        ConvexShapeType form_type,
        bool verbose)
        : GeoModelTopologyRecovererBase<2>(geomodel, form_type, verbose)
    {
    }

    GeoModelTopologyRecoverer<3>::GeoModelTopologyRecoverer(
        RINGMesh::GeoModel3D &geomodel,
        ConvexShapeType form_type,
        bool verbose)
        : GeoModelTopologyRecovererBase<3>(geomodel, form_type, verbose)
    {
    }

    template <index_t DIMENSION>
    GeoModelTopologyRecovererBase<DIMENSION>::~GeoModelTopologyRecovererBase()
    {
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::initialize_graph_nodes()
    {
        for (const auto &corner : geomodel_.corners())
        {
            NodeInformation corner_info(RINGMesh::corner_type_name_static());
            corner_info.parts.emplace_back(corner.gmme(), LineChunk(0, 0));
            add_node_into_invalidity_graph(corner_info);
        }
        for (const auto &line : geomodel_.lines())
        {
            NodeInformation line_info(RINGMesh::line_type_name_static());
            line_info.parts.emplace_back(line.gmme(),
                                         LineChunk(0, line.nb_mesh_elements()));
            line_info.ordered_boundary_nodes.push_back(line.boundary(0).index());
            line_info.ordered_boundary_nodes.push_back(line.boundary(1).index());
            add_node_into_invalidity_graph(line_info);
        }
    }

    void GeoModelTopologyRecoverer<3>::initialize_graph_nodes()
    {
        GeoModelTopologyRecovererBase<3>::initialize_graph_nodes();
        for (const auto &surface : geomodel_.surfaces())
        {
            NodeInformation surf_info(RINGMesh::surface_type_name_static());
            surf_info.parts.emplace_back(surface.gmme(),
                                         LineChunk(0, surface.nb_mesh_elements()));
            for (auto b_id : RINGMesh::range(surface.nb_boundaries()))
            {
                surf_info.ordered_boundary_nodes.push_back(
                    surface.boundary(b_id).index());
            }
            add_node_into_invalidity_graph(surf_info);
        }
    }

    template <index_t DIMENSION>
    index_t GeoModelTopologyRecovererBase<DIMENSION>::add_edge_into_invalidity_graph(
        const EdgeInformation &new_edge)
    {
        index_t new_edge_id = static_cast<index_t>(edge_information.size());
        // Check on edges
        scar_assert(new_edge.node1_id != NO_ID);
        scar_assert(new_edge.node2_id != NO_ID);
        scar_assert(new_edge.is_active());
        scar_assert(new_edge.node1.id.is_defined());
        scar_assert(!new_edge.node1.parts.empty());
        for (const auto &node1_part : new_edge.node1.parts)
        {
            scar_unused(node1_part);
            scar_assert(node1_part.begin <= node1_part.end);
        }
        scar_assert(new_edge.node2.id.is_defined());
        scar_assert(!new_edge.node2.parts.empty());
        for (const auto &node2_part : new_edge.node2.parts)
        {
            scar_unused(node2_part);
            scar_assert(node2_part.begin <= node2_part.end);
        }
        // If ok, add it and return position in vector
        edge_information.push_back(new_edge);
        invalidity_graph(new_edge.node1_id, new_edge.node2_id) = new_edge_id;
        return new_edge_id;
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::remove_edge_from_invalidity_graph(
        EdgeInformation &edge)
    {
        edge.set_inactive();
        invalidity_graph(edge.node1_id, edge.node2_id) = NO_ID;
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::remove_edge_from_invalidity_graph(
        index_t edge_id)
    {
        remove_edge_from_invalidity_graph(edge_information[edge_id]);
    }

    template <index_t DIMENSION>
    index_t GeoModelTopologyRecovererBase<DIMENSION>::add_node_into_invalidity_graph(
        const NodeInformation &new_node)
    {
        index_t new_node_id = static_cast<index_t>(node_information.size());
        // Check on edges
        scar_assert(new_node.is_active());
        scar_assert(!new_node.parts.empty());
        for (const auto &part : new_node.parts)
        {
            scar_assert(part.id.is_defined());
            if (part.part.begin > part.part.end)
            {
                scar_assert(part.part.begin <= part.part.end);
            }
        }
        if (new_node.type == RINGMesh::corner_type_name_static())
        {
            // No checks
        }
        else if (new_node.type == RINGMesh::line_type_name_static())
        {
            scar_assert(new_node.ordered_boundary_nodes.size() == 2);
        }
        else if (new_node.type == RINGMesh::surface_type_name_static())
        {
            scar_assert(!new_node.ordered_boundary_nodes.empty());
        }
        else
        {
            scar_assert_not_reached
        }
        // If ok, add it and return position in vector
        node_information.push_back(new_node);

        // Check capacity of graphs' correlation maps
        if (new_node_id >= invalidity_graph.nb_elem())
        {
            CorrelationMap<index_t> new_invalidity_graph(2 * new_node_id, NO_ID);
            CorrelationMap<bool> new_connectivity_graph(2 * new_node_id, false);
            for (auto i : RINGMesh::range(invalidity_graph.nb_elem()))
            {
                for (auto j : RINGMesh::range(invalidity_graph.nb_elem()))
                {
                    if (j <= i)
                    {
                        continue;
                    }
                    new_invalidity_graph(i, j) = invalidity_graph(i, j);
                    new_connectivity_graph(i, j) = connectivity_graph(i, j);
                }
            }
            invalidity_graph.reset_and_resize(2 * new_node_id, NO_ID);
            connectivity_graph.reset_and_resize(2 * new_node_id, false);
            for (auto i : RINGMesh::range(invalidity_graph.nb_elem()))
            {
                for (auto j : RINGMesh::range(invalidity_graph.nb_elem()))
                {
                    if (j <= i)
                    {
                        continue;
                    }
                    invalidity_graph(i, j) = new_invalidity_graph(i, j);
                    connectivity_graph(i, j) = new_connectivity_graph(i, j);
                }
            }
        }
        return new_node_id;
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::analyze_corner_corner()
    {
        index_t nb_fhc{0};
        for (const auto &corner1 : geomodel_.corners())
        {
            bool free_horizon_border = corner_is_free_horizon_border(corner1);
            auto corner1_exclusion = get_exclusion_shape(shape_type_,
                                                         corner1.mesh(), 0);
            if (free_horizon_border)
            {
                nb_fhc++;
                find_nearest_entity_and_add_graph_edge(corner1);
            }
            else
            {
                for (const auto &corner2 : geomodel_.corners())
                {
                    if (corner2.index() <= corner1.index())
                    {
                        continue;
                    }
                    if (invalidity_graph(corner1.index(), corner2.index()) != NO_ID)
                    {
                        // Already an edge between them (free_horizon border)
                        continue;
                    }
                    auto corner2_exclusion = get_exclusion_shape(shape_type_,
                                                                 corner2.mesh(), 0);
                    bool intersection = are_convex_intersecting(*corner1_exclusion,
                                                                *corner2_exclusion);
                    if (intersection)
                    {
                        // Add an edge in the graph between corner1 and corner2
                        EntityMultipleParts node_corner1;
                        node_corner1.id = corner1.gmme();
                        node_corner1.parts.emplace_back(0, 0);
                        EntityMultipleParts node_corner2;
                        node_corner2.id = corner2.gmme();
                        node_corner2.parts.emplace_back(0, 0);

                        EdgeInformation new_corner_corner_edge(node_corner1,
                                                               node_corner2);
                        new_corner_corner_edge.node1_id = corner1.index();
                        new_corner_corner_edge.node2_id = corner2.index();
                        add_edge_into_invalidity_graph(new_corner_corner_edge);

                        invalid_vertices.insert(
                            geomodel_.mesh.vertices.geomodel_vertex_id(
                                corner1.gmme(), 0));
                        invalid_vertices.insert(
                            geomodel_.mesh.vertices.geomodel_vertex_id(
                                corner2.gmme(), 0));
                    }
                }
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::analyze_corner_line()
    {
        for (const auto &corner : geomodel_.corners())
        {
            auto corner_exclusion = get_exclusion_shape(shape_type_, corner.mesh(),
                                                        0);
            auto connected_corners = CorrelationMapAPI::directly_linked(
                invalidity_graph, corner.index(), NO_ID, true);
            for (const auto &line : geomodel_.lines())
            {
                if (invalidity_graph(corner.index(),
                                     geomodel_.nb_corners() + line.index()) != NO_ID)
                {
                    // Already an edge between them (free horizon border)
                    continue;
                }
                if (!respect_geological_rules(corner, line))
                {
                    continue;
                }
                // Computes line intersections with corner exclusion zone
                NearbyLineChunksFromPoint<DIMENSION> line_chunks_inside_tolerance(
                    *corner_exclusion, line.mesh());
                // Remove LineChunk having a line corner directly linked to current corner
                if (line_chunks_inside_tolerance.nb_chunks() > 0)
                {
                    LineChunks edge_line_chunks;
                    edge_line_chunks.reserve(
                        line_chunks_inside_tolerance.nb_chunks());
                    for (const auto &chunk : line_chunks_inside_tolerance.chunks())
                    {
                        if (chunk.begin == 0)
                        {
                            if (RINGMesh::contains(connected_corners,
                                                   line.boundary(0).index()))
                            {
                                continue;
                            }
                        }
                        if (chunk.end == line.nb_mesh_elements())
                        {
                            if (RINGMesh::contains(connected_corners,
                                                   line.boundary(1).index()))
                            {
                                continue;
                            }
                        }
                        edge_line_chunks.push_back(std::move(chunk));
                    }
                    if (!edge_line_chunks.empty())
                    {
                        // Add an edge in the graph between the corner and the line
                        EntityMultipleParts node_corner;
                        node_corner.id = corner.gmme();
                        node_corner.parts.emplace_back(0, 0);
                        EntityMultipleParts node_line;
                        node_line.id = line.gmme();
                        node_line.parts = std::move(edge_line_chunks);

                        EdgeInformation new_corner_line_edge(node_corner,
                                                             node_line);
                        new_corner_line_edge.node1_id = corner.index();
                        new_corner_line_edge.node2_id = geomodel_.nb_corners() + line.index();
                        ;
                        add_edge_into_invalidity_graph(new_corner_line_edge);

                        invalid_vertices.insert(
                            geomodel_.mesh.vertices.geomodel_vertex_id(
                                corner.gmme(), 0));
                        for (const auto &line_part : new_corner_line_edge.node2.parts)
                        {
                            for (auto i : line_part.get_range())
                            {
                                invalid_vertices.insert(
                                    geomodel_.mesh.vertices.geomodel_vertex_id(
                                        line.gmme(), i));
                            }
                        }

                        for (const auto &line_part : new_corner_line_edge.node2.parts)
                        {
                            index_t v_begin = static_cast<index_t>(std::floor(
                                line_part.begin));
                            index_t v_end = static_cast<index_t>(std::ceil(
                                line_part.end));
                            for (auto i : RINGMesh::range(v_begin, v_end))
                            {
                                invalid_segments.insert(
                                    geomodel_.mesh.edges.edge(line.index(), i));
                            }
                        }
                    }
                }
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::analyze_line_line()
    {
        GEO::ProgressTask progress_bar("Line-Line analysis", geomodel_.nb_lines());
        if (!enum_contains(geological_rules,
                           GeologicalRule::NO_LINE_LINE_CHECK))
        {
            // @todo Find a better way to optimize (and to not check all possibilities)
            for (auto &line1 : geomodel_.lines())
            {
                progress_bar.next();
                RINGMesh::Box<DIMENSION> line1_bbox;
                for (const auto &v : RINGMesh::range(line1.nb_vertices()))
                {
                    line1_bbox.add_point(line1.vertex(v));
                }
                for (auto &line2 : geomodel_.lines())
                {
                    if (line2.index() <= line1.index())
                    {
                        continue;
                    }
                    RINGMesh::Box<DIMENSION> line2_bbox;
                    for (const auto &v : RINGMesh::range(line2.nb_vertices()))
                    {
                        line2_bbox.add_point(line2.vertex(v));
                    }
                    if (!line2_bbox.bboxes_overlap(line1_bbox))
                    {
                        continue;
                    }

                    if (enum_contains(geological_rules,
                                      GeologicalRule::CHECK_INTERSECTIONS))
                    {
                        bool segment_intersection;
                        std::vector<vecn<DIMENSION>> pts;
                        std::tie(segment_intersection, pts) =
                            intersection_line_line(line1, line2);
                        if (segment_intersection)
                        {
                            EdgeInformation line_line_intersections;
                            EntityMultipleParts node1_parts;
                            node1_parts.id = line1.gmme();
                            EntityMultipleParts node2_parts;
                            node2_parts.id = line2.gmme();
                            for (auto &pt : pts)
                            {
                                auto rel_coord1 =
                                    PointOnLineCoordinates<DIMENSION>::relative_coordinates(
                                        pt, line1.mesh());
                                node1_parts.parts.emplace_back(rel_coord1,
                                                               rel_coord1);
                                auto rel_coord2 =
                                    PointOnLineCoordinates<DIMENSION>::relative_coordinates(
                                        pt, line2.mesh());
                                node2_parts.parts.emplace_back(rel_coord2,
                                                               rel_coord2);
                            }

                            line_line_intersections.node1 = node1_parts;
                            line_line_intersections.node1_id = geomodel_.nb_corners() + line1.index();
                            line_line_intersections.node2 = node2_parts;
                            line_line_intersections.node2_id = geomodel_.nb_corners() + line2.index();
                            add_edge_into_invalidity_graph(
                                line_line_intersections);

                            for (const auto &line_part : line_line_intersections.node1.parts)
                            {
                                for (auto i : line_part.get_internal_range())
                                {
                                    invalid_vertices.insert(
                                        geomodel_.mesh.vertices.geomodel_vertex_id(
                                            line1.gmme(), i));
                                }
                            }
                            for (const auto &line_part : line_line_intersections.node2.parts)
                            {
                                for (auto i : line_part.get_internal_range())
                                {
                                    invalid_vertices.insert(
                                        geomodel_.mesh.vertices.geomodel_vertex_id(
                                            line2.gmme(), i));
                                }
                            }

                            for (const auto &line_part : line_line_intersections.node1.parts)
                            {
                                index_t v_begin = static_cast<index_t>(std::floor(
                                    line_part.begin));
                                index_t v_end = static_cast<index_t>(std::ceil(
                                    line_part.end));
                                for (auto i : RINGMesh::range(v_begin, v_end))
                                {
                                    invalid_segments.insert(
                                        geomodel_.mesh.edges.edge(line1.index(),
                                                                  i));
                                }
                            }

                            for (const auto &line_part : line_line_intersections.node2.parts)
                            {
                                index_t v_begin = static_cast<index_t>(std::floor(
                                    line_part.begin));
                                index_t v_end = static_cast<index_t>(std::ceil(
                                    line_part.end));
                                for (auto i : RINGMesh::range(v_begin, v_end))
                                {
                                    invalid_segments.insert(
                                        geomodel_.mesh.edges.edge(line2.index(),
                                                                  i));
                                }
                            }
                        }
                    }

                    bool build_edge;
                    EdgeInformation line_line_overlaps;
                    std::tie(build_edge, line_line_overlaps) =
                        analyze_line_pair_exclusion_zones(line1, line2);
                    if (build_edge)
                    {
                        line_line_overlaps.node1_id = geomodel_.nb_corners() + line1.index();
                        line_line_overlaps.node2_id = geomodel_.nb_corners() + line2.index();
                        add_edge_into_invalidity_graph(line_line_overlaps);

                        // For invalidity output
                        for (const auto &line_part : line_line_overlaps.node1.parts)
                        {
                            for (auto i : line_part.get_internal_range())
                            {
                                invalid_vertices.insert(
                                    geomodel_.mesh.vertices.geomodel_vertex_id(
                                        line1.gmme(), i));
                            }
                        }
                        for (const auto &line_part : line_line_overlaps.node2.parts)
                        {
                            for (auto i : line_part.get_internal_range())
                            {
                                invalid_vertices.insert(
                                    geomodel_.mesh.vertices.geomodel_vertex_id(
                                        line2.gmme(), i));
                            }
                        }

                        for (const auto &line_part : line_line_overlaps.node1.parts)
                        {
                            index_t v_begin = static_cast<index_t>(std::floor(
                                line_part.begin));
                            index_t v_end = static_cast<index_t>(std::ceil(
                                line_part.end));
                            for (auto i : RINGMesh::range(v_begin, v_end))
                            {
                                invalid_segments.insert(
                                    geomodel_.mesh.edges.edge(line1.index(), i));
                            }
                        }

                        for (const auto &line_part : line_line_overlaps.node2.parts)
                        {
                            index_t v_begin = static_cast<index_t>(std::floor(
                                line_part.begin));
                            index_t v_end = static_cast<index_t>(std::ceil(
                                line_part.end));
                            for (auto i : RINGMesh::range(v_begin, v_end))
                            {
                                invalid_segments.insert(
                                    geomodel_.mesh.edges.edge(line2.index(), i));
                            }
                        }
                    }
                }
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::build_graph_edges_from_geomodel_analysis()
    {
        if (verbose_)
        {
            Logger::out("Analysis", "Analyze Corner - Corner");
        }
        analyze_corner_corner();
        if (verbose_)
        {
            Logger::out("Analysis", "Analyze Corner - Line");
        }
        analyze_corner_line();
        if (verbose_)
        {
            Logger::out("Analysis", "Analyze Line - Line");
        }
        analyze_line_line();
    }

    void GeoModelTopologyRecoverer<3>::build_graph_edges_from_geomodel_analysis()
    {
        GeoModelTopologyRecovererBase<3>::build_graph_edges_from_geomodel_analysis();
        if (verbose_)
        {
            Logger::out("Analysis", "Analyze Corner - Surface");
        }
        analyze_corner_surface();
        if (verbose_)
        {
            Logger::out("Analysis", "Analyze Line - Surface");
        }
        analyze_line_surface();
        if (verbose_)
        {
            Logger::out("Analysis", "Analyze Surface - Surface");
        }
        analyze_surface_surface();
    }

    void GeoModelTopologyRecoverer<3>::analyze_corner_surface()
    {
        GEO::ProgressTask progress_bar("Corner-Surf analysis",
                                       geomodel_.nb_corners());
        for (const auto &corner : geomodel_.corners())
        {
            progress_bar.next();
            auto corner_exclusion = get_exclusion_shape(shape_type_, corner.mesh(),
                                                        0);
            RINGMesh::Box3D corner_v_bbox;
            vec3 direction;
            direction.x = 1.;
            corner_v_bbox.add_point(
                corner_exclusion->center() + 3 * corner_exclusion->width(direction) * direction);
            corner_v_bbox.add_point(
                corner_exclusion->center() - 3 * corner_exclusion->width(direction) * direction);
            direction.x = 0.;
            direction.y = 1.;
            corner_v_bbox.add_point(
                corner_exclusion->center() + 3 * corner_exclusion->width(direction) * direction);
            corner_v_bbox.add_point(
                corner_exclusion->center() - 3 * corner_exclusion->width(direction) * direction);
            direction.z = 1.;
            direction.y = 0.;
            corner_v_bbox.add_point(
                corner_exclusion->center() + 3 * corner_exclusion->width(direction) * direction);
            corner_v_bbox.add_point(
                corner_exclusion->center() - 3 * corner_exclusion->width(direction) * direction);

            auto connected_corners = CorrelationMapAPI::directly_linked(
                invalidity_graph, corner.index(), NO_ID, true);
            for (const auto &surface : geomodel_.surfaces())
            {
                if (invalidity_graph(corner.index(),
                                     geomodel_.nb_corners() + geomodel_.nb_lines() + surface.index()) != NO_ID)
                {
                    // Already an edge between them (free horizon border)
                    continue;
                }

                auto surf_corners = surface_corners(geomodel_, surface.index());
                if (is_element_in(corner.index(), surf_corners))
                {
                    continue;
                }

                for (auto e_id : RINGMesh::range(surface.nb_mesh_elements()))
                {
                    auto v0 = surface.mesh_element_vertex({e_id, 0});
                    auto v1 = surface.mesh_element_vertex({e_id, 1});
                    auto v2 = surface.mesh_element_vertex({e_id, 2});
                    RINGMesh::Box3D cur_bbox;
                    cur_bbox.add_point(v0);
                    cur_bbox.add_point(v1);
                    cur_bbox.add_point(v2);
                    if (!corner_v_bbox.bboxes_overlap(cur_bbox))
                    {
                        continue;
                    }

                    RINGMesh::Geometry::Triangle3D cur_elmt(v0, v1, v2);
                    double distance;
                    vec3 nearest_pt;
                    std::tie(distance, nearest_pt) =
                        RINGMesh::Distance::point_to_triangle(corner.vertex(0),
                                                              cur_elmt);

                    bool success;
                    std::array<double, 3> bary_coords;
                    std::tie(success, bary_coords) =
                        RINGMesh::triangle_barycentric_coordinates(nearest_pt, v0,
                                                                   v1, v2);

                    if (!success)
                    {
                        continue;
                    }
                    auto surf_exclusion = get_exclusion_shape(shape_type_,
                                                              surface.mesh(), e_id, bary_coords);

                    bool projector_intersect = are_convex_intersecting(
                        *corner_exclusion, *surf_exclusion);
                    if (projector_intersect)
                    {
                        // Output for invalidity
                        invalid_vertices.insert(
                            geomodel_.mesh.vertices.geomodel_vertex_id(
                                corner.gmme(), 0));
                        invalid_triangles.insert({surface.index(), e_id});
                    }
                }
            }
        }
    }

    bool line_is_boundary_of_surface(
        const RINGMesh::Line3D &line,
        const RINGMesh::Surface3D &surface)
    {
        for (auto i : RINGMesh::range(line.nb_incident_entities()))
        {
            if (line.incident_entity(i).index() == surface.index())
            {
                return true;
            }
        }
        return false;
    }

    void GeoModelTopologyRecoverer<3>::analyze_line_surface()
    {
        std::vector<RINGMesh::Box3D> surface_boxes(geomodel_.nb_surfaces());

        for (const auto &surface : geomodel_.surfaces())
        {
            index_t s_id = surface.index();
            for (const auto &v : RINGMesh::range(surface.nb_vertices()))
            {
                surface_boxes[s_id].add_point(surface.vertex(v));
            }
        }

        GEO::ProgressTask progress_bar("Line-Surf analysis", geomodel_.nb_lines());
        for (const auto &line : geomodel_.lines())
        {
            progress_bar.next();
            for (auto v : RINGMesh::range(line.nb_vertices()))
            {
                auto line_v_exclusion = get_exclusion_shape(shape_type_,
                                                            line.mesh(), v);
                RINGMesh::Box3D line_v_bbox;
                vec3 direction;
                direction.x = 1.;
                line_v_bbox.add_point(
                    line_v_exclusion->center() + 3 * line_v_exclusion->width(direction) * direction);
                line_v_bbox.add_point(
                    line_v_exclusion->center() - 3 * line_v_exclusion->width(direction) * direction);
                direction.x = 0.;
                direction.y = 1.;
                line_v_bbox.add_point(
                    line_v_exclusion->center() + 3 * line_v_exclusion->width(direction) * direction);
                line_v_bbox.add_point(
                    line_v_exclusion->center() - 3 * line_v_exclusion->width(direction) * direction);
                direction.z = 1.;
                direction.y = 0.;
                line_v_bbox.add_point(
                    line_v_exclusion->center() + 3 * line_v_exclusion->width(direction) * direction);
                line_v_bbox.add_point(
                    line_v_exclusion->center() - 3 * line_v_exclusion->width(direction) * direction);

                for (const auto &surface : geomodel_.surfaces())
                {
                    if (line_is_boundary_of_surface(line, surface))
                    {
                        continue;
                    }

                    if (!line_v_bbox.bboxes_overlap(
                            surface_boxes[surface.index()]))
                    {
                        continue;
                    }

                    for (auto e_id : RINGMesh::range(surface.nb_mesh_elements()))
                    {
                        auto v0 = surface.mesh_element_vertex({e_id, 0});
                        auto v1 = surface.mesh_element_vertex({e_id, 1});
                        auto v2 = surface.mesh_element_vertex({e_id, 2});
                        RINGMesh::Box3D cur_bbox;
                        cur_bbox.add_point(v0);
                        cur_bbox.add_point(v1);
                        cur_bbox.add_point(v2);
                        if (!line_v_bbox.bboxes_overlap(cur_bbox))
                        {
                            continue;
                        }

                        RINGMesh::Geometry::Triangle3D cur_elmt(v0, v1, v2);
                        double distance;
                        vec3 nearest_pt;
                        std::tie(distance, nearest_pt) =
                            RINGMesh::Distance::point_to_triangle(line.vertex(v),
                                                                  cur_elmt);

                        bool success;
                        std::array<double, 3> bary_coords;
                        std::tie(success, bary_coords) =
                            RINGMesh::triangle_barycentric_coordinates(nearest_pt,
                                                                       v0, v1, v2);

                        if (!success)
                        {
                            continue;
                        }
                        auto surf_exclusion = get_exclusion_shape(shape_type_,
                                                                  surface.mesh(), e_id, bary_coords);

                        bool projector_intersect = are_convex_intersecting(
                            *line_v_exclusion, *surf_exclusion);
                        if (projector_intersect)
                        {
                            // Output for invalidity
                            invalid_triangles.insert({surface.index(), e_id});
                        }
                    }
                }
            }
        }
    }
    void GeoModelTopologyRecoverer<3>::analyze_surface_surface()
    {
        std::vector<RINGMesh::Box3D> surface_boxes(geomodel_.nb_surfaces());

        for (const auto &surface : geomodel_.surfaces())
        {
            index_t s_id = surface.index();
            for (const auto &v : RINGMesh::range(surface.nb_vertices()))
            {
                surface_boxes[s_id].add_point(surface.vertex(v));
            }
        }

        GEO::ProgressTask progress_bar("Surf-Surf analysis",
                                       geomodel_.nb_surfaces());
        for (const auto &surface : geomodel_.surfaces())
        {
            progress_bar.next();
            for (auto v : RINGMesh::range(surface.nb_vertices()))
            {
                auto surface_v_exclusion = get_exclusion_shape(shape_type_,
                                                               surface.mesh(), v);
                RINGMesh::Box3D surface_v_bbox;
                vec3 direction;
                direction.x = 1.;
                surface_v_bbox.add_point(
                    surface_v_exclusion->center() + 2 * surface_v_exclusion->width(direction) * direction);
                surface_v_bbox.add_point(
                    surface_v_exclusion->center() - 2 * surface_v_exclusion->width(direction) * direction);
                direction.x = 0.;
                direction.y = 1.;
                surface_v_bbox.add_point(
                    surface_v_exclusion->center() + 2 * surface_v_exclusion->width(direction) * direction);
                surface_v_bbox.add_point(
                    surface_v_exclusion->center() - 2 * surface_v_exclusion->width(direction) * direction);
                direction.z = 1.;
                direction.y = 0.;
                surface_v_bbox.add_point(
                    surface_v_exclusion->center() + 2 * surface_v_exclusion->width(direction) * direction);
                surface_v_bbox.add_point(
                    surface_v_exclusion->center() - 2 * surface_v_exclusion->width(direction) * direction);

                for (const auto &surface2 : geomodel_.surfaces())
                {

                    if (surface2.index() == surface.index())
                    {
                        continue;
                    }

                    if (!surface_v_bbox.bboxes_overlap(
                            surface_boxes[surface2.index()]))
                    {
                        continue;
                    }

                    for (auto e_id : RINGMesh::range(surface2.nb_mesh_elements()))
                    {
                        auto v0 = surface2.mesh_element_vertex({e_id, 0});
                        auto v1 = surface2.mesh_element_vertex({e_id, 1});
                        auto v2 = surface2.mesh_element_vertex({e_id, 2});
                        RINGMesh::Box3D cur_bbox;
                        cur_bbox.add_point(v0);
                        cur_bbox.add_point(v1);
                        cur_bbox.add_point(v2);
                        if (!surface_v_bbox.bboxes_overlap(cur_bbox))
                        {
                            continue;
                        }

                        RINGMesh::Geometry::Triangle3D cur_elmt(v0, v1, v2);
                        double distance;
                        vec3 nearest_pt;
                        std::tie(distance, nearest_pt) =
                            RINGMesh::Distance::point_to_triangle(
                                surface.vertex(v), cur_elmt);

                        bool success;
                        std::array<double, 3> bary_coords;
                        std::tie(success, bary_coords) =
                            RINGMesh::triangle_barycentric_coordinates(nearest_pt,
                                                                       v0, v1, v2);

                        if (!success)
                        {
                            continue;
                        }
                        auto surf2_exclusion = get_exclusion_shape(shape_type_,
                                                                   surface2.mesh(), e_id, bary_coords);

                        bool projector_intersect = are_convex_intersecting(
                            *surface_v_exclusion, *surf2_exclusion);
                        if (projector_intersect)
                        {
                            // Output for invalidity
                            invalid_triangles.insert({surface2.index(), e_id});
                        }
                    }
                }
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::build_graph_edges_from_geomodel_connectivity()
    {
        auto nb_corners = geomodel_.nb_corners();
        for (auto &line : geomodel_.lines())
        {
            for (auto i : RINGMesh::range(line.nb_boundaries()))
            {
                auto boundary_corner_id = line.boundary(i).index();
                connectivity_graph(boundary_corner_id, nb_corners + line.index()) =
                    true;
            }
        }
    }

    void GeoModelTopologyRecoverer<3>::build_graph_edges_from_geomodel_connectivity()
    {
        GeoModelTopologyRecovererBase<3>::build_graph_edges_from_geomodel_connectivity();
        auto nb_corners = geomodel_.nb_corners();
        auto nb_lines = geomodel_.nb_lines();
        for (auto &surf : geomodel_.surfaces())
        {
            for (auto i : RINGMesh::range(surf.nb_boundaries()))
            {
                auto boundary_line_id = surf.boundary(i).index();
                connectivity_graph(nb_corners + boundary_line_id,
                                   nb_corners + nb_lines + surf.index()) = true;
            }
        }
    }

    template <index_t DIMENSION>
    bool GeoModelTopologyRecovererBase<DIMENSION>::nodes_are_all_corners(
        std::vector<index_t> connected_component) const
    {
        for (const auto &node_id : connected_component)
        {
            scar_assert(node_information[node_id].is_active());
            if (node_information[node_id].type != RINGMesh::corner_type_name_static())
            {
                return false;
            }
        }
        return true;
    }

    template <index_t DIMENSION>
    bool GeoModelTopologyRecovererBase<DIMENSION>::nodes_are_all_lines(
        std::vector<index_t> connected_component) const
    {
        for (const auto &node_id : connected_component)
        {
            scar_assert(node_information[node_id].is_active());
            if (node_information[node_id].type != RINGMesh::line_type_name_static())
            {
                return false;
            }
        }
        return true;
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::print_graph_edges(
        bool verbose) const
    {
        for (const auto &edge : edge_information)
        {
            if (!edge.is_active())
            {
                continue;
            }
            std::cout << "Edge: " << edge.node1_id << " - " << edge.node2_id << " ("
                      << node_information[edge.node1_id].type << "-"
                      << node_information[edge.node2_id].type << ")" << std::endl;
            if (verbose)
            {
                std::cout << "(" << edge.node1.parts[0] << "-" << edge.node2.parts[0]
                          << ")" << std::endl;
            }
            scar_assert(node_information[edge.node1_id].is_active());
            scar_assert(node_information[edge.node2_id].is_active());
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::print_graph_nodes(
        bool verbose) const
    {
        for (auto i : RINGMesh::range(node_information.size()))
        {
            const auto &node = node_information[i];
            if (!node.is_active())
            {
                continue;
            }
            std::cout << "Node: " << i << " (" << node.type << ")" << std::endl;
            for (const auto &part : node.parts)
            {
                std::cout << "     " << part.id << std::endl;
                if (verbose)
                {
                    std::cout << "(" << part.part.begin << "-" << part.part.end
                              << ")" << std::endl;
                }
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::print_connectivity_graph_edges() const
    {
        for (const auto &node1 : RINGMesh::range(connectivity_graph.nb_elem() - 1))
        {
            for (const auto &node2 : RINGMesh::range(node1 + 1,
                                                     connectivity_graph.nb_elem()))
            {
                if (connectivity_graph(node1, node2))
                {
                    std::cout << "Edge: " << node1 << " - " << node2 << std::endl;
                }
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::output_invalid_areas_base(
        RINGMesh::GeoModel<DIMENSION> &invalidity_model) const
    {
        RINGMesh::GeoModelBuilder<DIMENSION> builder(invalidity_model);

        for (auto p_id : invalid_vertices)
        {
            builder.topology.find_or_create_corner(
                geomodel_.mesh.vertices.vertex(p_id));
        }

        auto cur_gmme = builder.topology.create_mesh_entity(
            RINGMesh::line_type_name_static());
        std::vector<vecn<DIMENSION>> vertex_list;
        vertex_list.reserve(2 * invalid_segments.size());
        for (auto e_id : invalid_segments)
        {
            auto next_vertex = geomodel_.mesh.vertices.vertex(
                geomodel_.mesh.edges.vertex({e_id, 0}));
            if (next_vertex != vertex_list.back())
            {
                builder.geometry.set_line(cur_gmme.index(), vertex_list);
                vertex_list.clear();
                cur_gmme = builder.topology.create_mesh_entity(
                    RINGMesh::line_type_name_static());
                vertex_list.push_back(next_vertex);
            }
            vertex_list.push_back(
                geomodel_.mesh.vertices.vertex(
                    geomodel_.mesh.edges.vertex({e_id, 1})));
        }
        builder.geometry.set_line(cur_gmme.index(), vertex_list);
    }
    void GeoModelTopologyRecoverer<3>::output_invalid_areas() const
    {
        RINGMesh::GeoModel3D invalidity_model;
        GeoModelTopologyRecovererBase<3>::output_invalid_areas_base(
            invalidity_model);

        RINGMesh::GeoModelBuilder3D builder(invalidity_model);
        auto cur_gmme = builder.topology.create_mesh_entity(
            RINGMesh::surface_type_name_static());
        std::vector<vec3> vertex_list;
        vertex_list.reserve(2 * invalid_triangles.size());
        std::vector<index_t> polygon_ptr;
        polygon_ptr.push_back(0);
        for (auto se_id : invalid_triangles)
        {
            if (geomodel_.surface(se_id.first).is_on_voi())
            {
                continue;
            }

            auto v0 = geomodel_.surface(se_id.first).mesh_element_vertex({se_id.second, 0});
            auto v1 = geomodel_.surface(se_id.first).mesh_element_vertex({se_id.second, 1});
            auto v2 = geomodel_.surface(se_id.first).mesh_element_vertex({se_id.second, 2});
            vertex_list.push_back(v0);
            vertex_list.push_back(v1);
            vertex_list.push_back(v2);
            polygon_ptr.push_back(static_cast<index_t>(vertex_list.size()));
        }
        std::vector<index_t> polygons(vertex_list.size());
        std::iota(polygons.begin(), polygons.end(), 0);
        builder.geometry.set_surface_geometry(cur_gmme.index(), vertex_list,
                                              polygons, polygon_ptr);

        RINGMesh::geomodel_save(invalidity_model, "/tmp/inv_model.gm");
    }
    void GeoModelTopologyRecoverer<2>::output_invalid_areas() const
    {
        RINGMesh::GeoModel2D invalidity_model;
        GeoModelTopologyRecovererBase<2>::output_invalid_areas_base(
            invalidity_model);
        RINGMesh::geomodel_save(invalidity_model, "/tmp/inv_model.gm");
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::output_graphs() const
    {
        RINGMesh::Box<DIMENSION> bbox;
        for (auto v : RINGMesh::range(geomodel_.mesh.vertices.nb()))
        {
            bbox.add_point(geomodel_.mesh.vertices.vertex(v));
        }
        std::ofstream outfile;
        outfile.open("/tmp/graph.dot");

        // Header and viz parameters
        outfile << "graph G { \n";

        outfile << "graph [mindist=0.5] \n";
        outfile << "node [shape = circle, style = filled] \n";
        outfile << "\n";

        // Output nodes
        outfile << "node[fillcolor = white] \n";
        for (const auto &node_id : RINGMesh::range(node_information.size()))
        {
            const auto &node = node_information[node_id];
            if (!node.is_active())
            {
                continue;
            }
            outfile << node_id;

            if (node.type == RINGMesh::corner_type_name_static())
            {
                outfile << " [ label = \"C" << node.parts.front().id.index()
                        << "\" , shape = circle];\n";
            }
            else if (node.type == RINGMesh::line_type_name_static())
            {
                outfile << " [ label = \"L" << node.parts.front().id.index()
                        << "\" , shape = square];\n";
            }
        }

        // Output connectivity edges
        outfile << "edge [color = blue, penwidth=2, len = 3, constraint=false ] \n";
        for (const auto &node1 : RINGMesh::range(connectivity_graph.nb_elem() - 1))
        {
            for (const auto &node2 : RINGMesh::range(node1 + 1,
                                                     connectivity_graph.nb_elem()))
            {
                if (connectivity_graph(node1, node2))
                {
                    outfile << node1 << "--" << node2 << "; \n";
                }
            }
        }

        outfile
            << "edge [color = red, penwidth = 5 , len = 1.5, minlen = 0.5, constraint=true ] \n";
        for (const auto &edge : edge_information)
        {
            if (!edge.is_active())
            {
                continue;
            }

            outfile << edge.node1_id << "--" << edge.node2_id << "; \n";
        }

        outfile << "} \n";
        outfile.close();
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::build_connectivity_graph()
    {
        // Graph initialization and building
        build_graph_edges_from_geomodel_connectivity();
        if (verbose_)
        {
            print_connectivity_graph_edges();
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::build_invalidity_graph()
    {
        // Graph initialization and building
        build_graph_edges_from_geomodel_analysis();

        index_t nb_triangles{0};
        for (auto t : RINGMesh::range(geomodel_.mesh.polygons.nb()))
        {
            for (auto i : RINGMesh::range(3))
            {
                if (invalid_vertices.find(
                        geomodel_.mesh.polygons.vertex({t, i})) != invalid_vertices.end())
                {
                    ++nb_triangles;
                    break;
                }
            }
        }

        if (verbose_)
        {
            print_graph_nodes(verbose_);
            print_graph_edges(verbose_);
        }
    }

    template <index_t DIMENSION>
    bool GeoModelTopologyRecovererBase<DIMENSION>::respect_geological_rules(
        const RINGMesh::Corner<DIMENSION> &free_horizon_corner,
        const RINGMesh::Corner<DIMENSION> &other_corner) const
    {
        if (corner_is_free_horizon_border(free_horizon_corner))
        {
            if (geological_rules == GeologicalRule::EMPTY)
            {
                return true;
            }
            else
            {
                if (enum_contains(geological_rules,
                                  GeologicalRule::CONFORMABLE_HORIZONS))
                {
                    for (auto incident_line_i : RINGMesh::range(
                             other_corner.nb_incident_entities()))
                    {
                        const auto &incident_line = other_corner.incident_entity(
                            incident_line_i);
                        if (incident_line.parent(
                                             RINGMesh::Interface<DIMENSION>::type_name_static())
                                .geological_feature() != RINGMesh::GeoModelGeologicalEntity<DIMENSION>::GEOL_FEATURE::STRATI)
                        {
                            return true;
                        }
                    }
                    return false;
                }
                return true;
            }
        }
        else
        {
            return true;
        }
    }

    template <index_t DIMENSION>
    bool GeoModelTopologyRecovererBase<DIMENSION>::respect_geological_rules(
        const RINGMesh::Corner<DIMENSION> &free_horizon_corner,
        const RINGMesh::Line<DIMENSION> &line) const
    {
        if (line.boundary(0).index() == free_horizon_corner.index())
        {
            return false;
        }
        if (line.boundary(1).index() == free_horizon_corner.index())
        {
            return false;
        }
        if (corner_is_free_horizon_border(free_horizon_corner))
        {
            if (geological_rules == GeologicalRule::EMPTY)
            {
                return true;
            }
            else
            {
                if (enum_contains(geological_rules,
                                  GeologicalRule::CONFORMABLE_HORIZONS))
                {
                    if (line.parent(
                                RINGMesh::Interface<DIMENSION>::type_name_static())
                            .geological_feature() != RINGMesh::GeoModelGeologicalEntity<DIMENSION>::GEOL_FEATURE::STRATI)
                    {
                        return true;
                    }
                    return false;
                }
                return true;
            }
        }
        else
        {
            if (enum_contains(geological_rules,
                              GeologicalRule::MERGE_WHOLE_LINES))
            {
                // If geological rules force merging only whole lines,
                // Corner - Line edges are useless
                return false;
            } // Else no rules, so return true;
            return true;
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::find_nearest_entity_and_add_graph_edge(
        const RINGMesh::Corner<DIMENSION> &free_horizon_corner)
    {
        double nearest_convex_distance = RINGMesh::max_float64();
        RINGMesh::gmme_id nearest_entity;
        auto fh_corner_exclusion = get_exclusion_shape(shape_type_,
                                                       free_horizon_corner.mesh(), 0);

        // Check for other corners
        for (const auto &corner : geomodel_.corners())
        {
            if (corner.index() == free_horizon_corner.index())
            {
                continue;
            }
            if (!respect_geological_rules(free_horizon_corner, corner))
            {
                continue;
            }
            auto convex_dist = convex_distance(fh_corner_exclusion,
                                               corner.vertex(0));
            if (convex_dist < nearest_convex_distance)
            {
                nearest_convex_distance = convex_dist;
                nearest_entity = corner.gmme();
            }
        }

        vecn<DIMENSION> nearest_to_corner_pt;
        LineRelativeCoordinates nearest_line_rel_coord;
        // Check for lines
        for (const auto &line : geomodel_.lines())
        {
            if (!respect_geological_rules(free_horizon_corner, line))
            {
                continue;
            }
            ConvexDistanceFunction<DIMENSION> metric(line.mesh(),
                                                     fh_corner_exclusion);
            index_t nearest_edge;
            vecn<DIMENSION> nearest_point;
            double convex_dist;

            std::tie(nearest_edge, nearest_point, convex_dist) =
                line.edge_aabb().closest_edge(fh_corner_exclusion->center(),
                                              metric);
            if (convex_dist < nearest_convex_distance)
            {
                nearest_convex_distance = convex_dist;
                nearest_entity = line.gmme();
                nearest_to_corner_pt = nearest_point;
            }
        }
        // Add a graph edge
        if (nearest_entity.type() == RINGMesh::Corner<DIMENSION>::type_name_static())
        {
            // Add an edge in the graph between the two corners
            EntityMultipleParts node_fh_corner;
            node_fh_corner.id = free_horizon_corner.gmme();
            node_fh_corner.parts.emplace_back(0, 0);
            EntityMultipleParts node_other_corner;
            node_other_corner.id = nearest_entity;
            node_other_corner.parts.emplace_back(0, 0);

            if (free_horizon_corner.index() < nearest_entity.index())
            {
                EdgeInformation new_edge(node_fh_corner, node_other_corner);
                new_edge.node1_id = free_horizon_corner.index();
                new_edge.node2_id = nearest_entity.index();
                add_edge_into_invalidity_graph(new_edge);
            }
            else
            {
                EdgeInformation new_edge(node_other_corner, node_fh_corner);
                new_edge.node1_id = nearest_entity.index();
                new_edge.node2_id = free_horizon_corner.index();
                add_edge_into_invalidity_graph(new_edge);
            }
        }
        else if (nearest_entity.type() == RINGMesh::Line<DIMENSION>::type_name_static())
        {

            vecn<DIMENSION> projection_along_horizon_line =
                find_oriented_projection_on_line(free_horizon_corner,
                                                 geomodel_.line(nearest_entity.index()));

            LineRelativeCoordinates nearest_line_rel_coord = PointOnLineCoordinates<
                DIMENSION>::relative_coordinates(projection_along_horizon_line,
                                                 geomodel_.line(nearest_entity.index()).mesh());

            // Add an edge in the graph between the corner and the line
            EntityMultipleParts node_fh_corner;
            node_fh_corner.id = free_horizon_corner.gmme();
            node_fh_corner.parts.emplace_back(0, 0);
            EntityMultipleParts node_line;
            node_line.id = nearest_entity;
            node_line.parts.emplace_back(nearest_line_rel_coord /*- 0.01*/,
                                         nearest_line_rel_coord /*+ 0.01*/);

            EdgeInformation new_corner_line_edge(node_fh_corner, node_line);
            new_corner_line_edge.node1_id = free_horizon_corner.index();
            new_corner_line_edge.node2_id = geomodel_.nb_corners() + nearest_entity.index();
            add_edge_into_invalidity_graph(new_corner_line_edge);
        }
    }

    template <index_t DIMENSION>
    std::vector<std::vector<index_t>> GeoModelTopologyRecovererBase<DIMENSION>::get_connected_components_with_edges()
    {
        std::vector<std::vector<index_t>> connected_components;
        CorrelationMapAPI::get_connected_components(invalidity_graph,
                                                    connected_components, NO_ID);
        std::vector<std::vector<index_t>> connected_components_with_edges;
        connected_components_with_edges.reserve(connected_components.size());
        for (const auto &cc : connected_components)
        {
            if (cc.size() > 1)
            {
                connected_components_with_edges.push_back(std::move(cc));
            }
        }
        return connected_components_with_edges;
    }

    template <index_t DIMENSION>
    EdgeType GeoModelTopologyRecovererBase<DIMENSION>::get_edge_type(
        const EdgeInformation &edge) const
    {
        const auto &type_node1 = node_information[edge.node1_id].type;
        const auto &type_node2 = node_information[edge.node2_id].type;
        if (type_node1 == RINGMesh::corner_type_name_static())
        {
            if (type_node2 == RINGMesh::corner_type_name_static())
            {
                return EdgeType::CornerCorner;
            }
            if (type_node2 == RINGMesh::line_type_name_static())
            {
                return EdgeType::CornerLine;
            }
            if (type_node2 == RINGMesh::surface_type_name_static())
            {
                return EdgeType::CornerSurface;
            }
            scar_assert_not_reached;
        }
        if (type_node1 == RINGMesh::line_type_name_static())
        {
            if (type_node2 == RINGMesh::line_type_name_static())
            {
                return EdgeType::LineLine;
            }
            if (type_node2 == RINGMesh::surface_type_name_static())
            {
                return EdgeType::LineSurface;
            }
            scar_assert_not_reached;
        }
        if (type_node1 == RINGMesh::surface_type_name_static())
        {
            if (type_node2 == RINGMesh::surface_type_name_static())
            {
                return EdgeType::SurfaceSurface;
            }
            scar_assert_not_reached;
        }
        scar_assert_not_reached;
        return EdgeType::CornerCorner;
    }

    // TODO : reimplemente this function
    template <index_t DIMENSION>
    bool GeoModelTopologyRecovererBase<DIMENSION>::are_line_line_reverse(
        const EdgeInformation &edge) const
    {
        scar_assert(get_edge_type(edge) == EdgeType::LineLine);
        scar_assert(edge.node1.parts.size() == edge.node2.parts.size());
        const auto &line_node1 = geomodel_.line(edge.node1.id.index());
        const auto &line_node2 = geomodel_.line(edge.node2.id.index());

        // Not safe request based on geometry
        const auto &corner10_position = line_node1.boundary(0).vertex(0);
        const auto &corner11_position = line_node1.boundary(1).vertex(0);
        const auto &corner20_position = line_node2.boundary(0).vertex(0);
        const auto &corner21_position = line_node2.boundary(1).vertex(0);

        if (corner10_position == corner20_position && corner11_position != corner21_position)
        {
            return false;
        }
        if (corner11_position == corner21_position && corner10_position != corner20_position)
        {
            return false;
        }
        if (corner11_position == corner20_position && corner10_position != corner21_position)
        {
            return true;
        }
        if (corner10_position == corner21_position && corner11_position != corner20_position)
        {
            return true;
        }

        if ((corner10_position - corner20_position).length2() < (corner10_position - corner21_position).length2())
        {
            if ((corner11_position - corner21_position).length2() < (corner11_position - corner20_position).length2())
            {
                return false;
            }
            scar_assert_not_reached;
        }
        return true;
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::remove_edges_from_invalidity_graph()
    {
        // Graph modification to remove edges
        remove_corner_line_edges();

        bool to_reiterate = true;
        while (to_reiterate)
        {
            to_reiterate = false;
            // Analysis of the graph
            std::vector<std::vector<index_t>> connected_components_with_edges =
                get_connected_components_with_edges();
            for (const auto &cc : connected_components_with_edges)
            {
                if (nodes_are_all_corners(cc))
                {
                    EntityParts new_node;
                    for (const auto &cc_node : cc)
                    {
                        for (const auto &entity_part : node_information[cc_node].parts)
                        {
                            new_node.push_back(entity_part);
                        }
                    }
                    NodeInformation new_corner_node(
                        RINGMesh::corner_type_name_static(), new_node);
                    index_t new_node_id = add_node_into_invalidity_graph(
                        new_corner_node);
                    for (const auto &cc_node : cc)
                    {
                        std::vector<index_t> incident_line_nodes =
                            CorrelationMapAPI::directly_linked(connectivity_graph,
                                                               cc_node, false, false);
                        for (auto line_node_id : incident_line_nodes)
                        {
                            connectivity_graph(cc_node, line_node_id) = false;
                            connectivity_graph(new_node_id, line_node_id) = true;
                            auto &line_node_info = node_information[line_node_id];
                            scar_assert(line_node_info.is_active());
                            scar_assert(
                                line_node_info.ordered_boundary_nodes.size() == 2);
                            if (line_node_info.ordered_boundary_nodes[0] == cc_node)
                            {
                                line_node_info.ordered_boundary_nodes[0] =
                                    new_node_id;
                            }
                            if (line_node_info.ordered_boundary_nodes[1] == cc_node)
                            {
                                line_node_info.ordered_boundary_nodes[1] =
                                    new_node_id;
                            }
                            if (line_node_info.ordered_boundary_nodes[0] == line_node_info.ordered_boundary_nodes[1])
                            {
                                // The line has to be removed since its boundaries are merged
                                line_node_info.set_inactive();
                                connectivity_graph(
                                    line_node_info.ordered_boundary_nodes[0],
                                    line_node_id) = false;
                                // Deactivate all connected edges
                                std::vector<index_t> line_linked_nodes =
                                    CorrelationMapAPI::directly_linked(
                                        invalidity_graph, line_node_id, NO_ID,
                                        false);
                                for (auto other_node : line_linked_nodes)
                                {
                                    auto edge_to_deactivate = invalidity_graph(
                                        line_node_id, other_node);
                                    remove_edge_from_invalidity_graph(
                                        edge_to_deactivate);
                                }
                                to_reiterate = true;
                            }
                        }
                    }
                    // Remove all edges and set as inactive node and edges
                    for (const auto &cc_node1 : cc)
                    {
                        node_information[cc_node1].set_inactive();
                        for (const auto &cc_node2 : cc)
                        {
                            if (cc_node1 == cc_node2)
                            {
                                continue;
                            }
                            if (invalidity_graph(cc_node1, cc_node2) != NO_ID)
                            {
                                remove_edge_from_invalidity_graph(
                                    invalidity_graph(cc_node1, cc_node2));
                            }
                        }
                    }
                    continue;
                }
                if (to_reiterate)
                {
                    break;
                }
            }
        }
        if (verbose_)
        {
            print_graph_nodes();
            print_graph_edges();
        }
        split_line_line_edges_and_nodes();

        if (verbose_)
        {
            print_graph_nodes();
            print_graph_edges();
        }
        perform_edge_contractions();

        if (verbose_)
        {
            print_graph_nodes();
            print_graph_edges();
        }

        scar_assert(CorrelationMapAPI::size(invalidity_graph, NO_ID) == 0);

        for (auto node_id : RINGMesh::range(node_information.size()))
        {
            const auto &node = node_information[node_id];
            if (!node.is_active())
            {
                continue;
            }
            if (node.type != RINGMesh::line_type_name_static())
            {
                continue;
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::split_line_line_edges_and_nodes()
    {
        for (auto n : node_information)
        {
            if (!n.is_active())
            {
                continue;
            }
            if (n.type != RINGMesh::line_type_name_static())
            {
                continue;
            }
        }
        // Check if each line has only one part to merge
        // Hold by information on edges
        bool be_modified{false};
        do
        {
            be_modified = false;
            for (auto &edge_info : edge_information)
            {
                if (!edge_info.is_active())
                {
                    continue;
                }
                scar_assert(
                    node_information[edge_info.node1_id].type == node_information[edge_info.node2_id].type);
                if (node_information[edge_info.node1_id].type != RINGMesh::line_type_name_static())
                {
                    continue;
                }

                bool edge_ok = check_line_line_edge_compared_to_nodes(edge_info);
                if (!edge_ok)
                {
                    perform_line_line_edge_split(edge_info);
                    be_modified = true;
                }
            }
        } while (be_modified);
    }

    template <index_t DIMENSION>
    std::vector<LineRelativeCoordinates> GeoModelTopologyRecovererBase<DIMENSION>::get_node_line_hinge_points(
        const NodeInformation &line_node,
        const EntityMultipleParts &edge_info_for_node)
    {
        std::vector<LineRelativeCoordinates> ordered_hinge_points;
        ordered_hinge_points.push_back(line_node.parts[0].part.begin);
        for (const auto &part : edge_info_for_node.parts)
        {
            if (part.begin != ordered_hinge_points.back())
            {
                ordered_hinge_points.push_back(part.begin);
            }
            if (part.end != ordered_hinge_points.back())
            {
                ordered_hinge_points.push_back(part.end);
            }
        }
        if (line_node.parts[0].part.end != ordered_hinge_points.back())
        {
            ordered_hinge_points.push_back(line_node.parts[0].part.end);
        }
        return ordered_hinge_points;
    }

    template <index_t DIMENSION>
    std::tuple<std::vector<index_t>, std::vector<index_t>,
               std::vector<index_t>>
    GeoModelTopologyRecovererBase<DIMENSION>::split_line_node_at_hinge_points(
        const EntityMultipleParts &edge_info_for_node,
        const NodeInformation &node_info,
        const std::vector<LineRelativeCoordinates> &hinge_points)
    {
        // Cut the line node in several new nodes
        // Do not forget to add node (Corner-type) at
        // hinge points (to update connectivity graph)
        scar_assert(
            node_info.ordered_boundary_nodes[0] != node_info.ordered_boundary_nodes[1]);
        index_t first_boundary = node_info.ordered_boundary_nodes[0];
        index_t last_boundary = node_info.ordered_boundary_nodes[1];
        std::vector<index_t> nodes_to_link_node(edge_info_for_node.parts.size(),
                                                NO_ID);
        std::vector<index_t> hinge_nodes_to_link;
        std::vector<index_t> new_nodes;
        for (auto i : RINGMesh::range(hinge_points.size() - 1))
        {
            LineChunk cur_line_chunk(hinge_points[i], hinge_points[i + 1]);
            if (cur_line_chunk.is_null())
            {
                continue; // Two identical successive hinge points
            }
            for (auto part_id : RINGMesh::range(edge_info_for_node.parts.size()))
            {
                if (cur_line_chunk.begin == edge_info_for_node.parts[part_id].begin && cur_line_chunk.end == edge_info_for_node.parts[part_id].end)
                {
                    nodes_to_link_node[part_id] =
                        static_cast<index_t>(node_information.size());
                }
            }
            EntityParts cur_parts;
            cur_parts.push_back(
                EntityPart(edge_info_for_node.id, cur_line_chunk));

            NodeInformation new_line_node(RINGMesh::line_type_name_static(),
                                          cur_parts);
            new_line_node.ordered_boundary_nodes.clear();
            new_line_node.ordered_boundary_nodes.push_back(first_boundary);
            new_line_node.ordered_boundary_nodes.push_back(last_boundary);
            index_t new_line_node_id = add_node_into_invalidity_graph(
                new_line_node);
            new_nodes.push_back(new_line_node_id);

            index_t second_boundary = last_boundary;
            if (i + 2 < hinge_points.size())
            {
                if (hinge_points[i + 1] != hinge_points.back())
                {
                    // No create of hinge point if it corresponds to the last boundary
                    EntityParts corner_parts;
                    corner_parts.push_back(
                        EntityPart(edge_info_for_node.id,
                                   LineChunk(hinge_points[i + 1],
                                             hinge_points[i + 1])));
                    NodeInformation hinge_node(RINGMesh::corner_type_name_static(),
                                               corner_parts);
                    second_boundary = add_node_into_invalidity_graph(hinge_node);
                    node_information[new_line_node_id].ordered_boundary_nodes[1] =
                        second_boundary;
                }
            }

            connectivity_graph(new_line_node_id, first_boundary) = true;
            connectivity_graph(new_line_node_id, second_boundary) = true;
            first_boundary = second_boundary;
        }
        // If no line-line true intersections
        if (std::count(nodes_to_link_node.begin(), nodes_to_link_node.end(), NO_ID) == 0)
        {
            for (auto i : nodes_to_link_node)
            {
                hinge_nodes_to_link.push_back(
                    node_information[i].ordered_boundary_nodes[0]);
                hinge_nodes_to_link.push_back(
                    node_information[i].ordered_boundary_nodes[1]);
            }
        }
        else
        {
            // Only intersections of lines
            scar_assert(
                std::count(nodes_to_link_node.begin(), nodes_to_link_node.end(),
                           NO_ID) == (int)nodes_to_link_node.size());
            nodes_to_link_node.clear();
        }
        return std::make_tuple(nodes_to_link_node, hinge_nodes_to_link, new_nodes);
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::perform_line_line_edge_split(
        EdgeInformation &edge_info)
    {
        const auto &node1_edge_info = edge_info.node1;
        auto &node1_node_info = node_information[edge_info.node1_id];
        scar_assert(node1_node_info.is_active());

        const auto &node2_edge_info = edge_info.node2;
        auto &node2_node_info = node_information[edge_info.node2_id];
        scar_assert(node2_node_info.is_active());

        scar_assert(node1_edge_info.parts.size() == node2_edge_info.parts.size());

        if (node1_node_info.parts.size() > 1 || node2_node_info.parts.size() > 1)
        {
            scar_assert_not_reached;
        }

        auto ordered_hinge_points_node1 = get_node_line_hinge_points(
            node1_node_info, node1_edge_info);
        auto ordered_hinge_points_node2 = get_node_line_hinge_points(
            node2_node_info, node2_edge_info);

        if (ordered_hinge_points_node1.size() == 2 || ordered_hinge_points_node2.size() == 2)
        {
            remove_edge_from_invalidity_graph(edge_info);
            return;
        }

        // Cut the line node in several new nodes
        std::vector<index_t> nodes_to_link_node1;
        std::vector<index_t> hinge_node_to_link_node1;
        std::vector<index_t> new_nodes1;
        std::tie(nodes_to_link_node1, hinge_node_to_link_node1, new_nodes1) =
            split_line_node_at_hinge_points(node1_edge_info, node1_node_info,
                                            ordered_hinge_points_node1);
        std::vector<index_t> nodes_to_link_node2;
        std::vector<index_t> hinge_node_to_link_node2;
        std::vector<index_t> new_nodes2;
        std::tie(nodes_to_link_node2, hinge_node_to_link_node2, new_nodes2) =
            split_line_node_at_hinge_points(node2_edge_info, node2_node_info,
                                            ordered_hinge_points_node2);

        // Reallocate this edge
        scar_assert(nodes_to_link_node1.size() == nodes_to_link_node2.size());
        for (auto i : RINGMesh::range(nodes_to_link_node1.size()))
        {
            EdgeInformation new_edge;
            new_edge.node1_id = nodes_to_link_node1[i];
            new_edge.node1.id = edge_info.node1.id;
            new_edge.node1.parts.emplace_back(
                node_information[nodes_to_link_node1[i]].parts[0].part.begin,
                node_information[nodes_to_link_node1[i]].parts[0].part.end);
            new_edge.node2_id = nodes_to_link_node2[i];
            new_edge.node2.id = edge_info.node2.id;
            new_edge.node2.parts.emplace_back(
                node_information[nodes_to_link_node2[i]].parts[0].part.begin,
                node_information[nodes_to_link_node2[i]].parts[0].part.end);
            add_edge_into_invalidity_graph(new_edge);
        }

        // Add edges between hinge points
        scar_assert(
            hinge_node_to_link_node1.size() == hinge_node_to_link_node2.size());

        for (auto i_node1 : RINGMesh::range(hinge_node_to_link_node1.size()))
        {

            // Check if lines are oriented in the same direction
            auto i_node2 =
                are_line_line_reverse(edge_info) ? hinge_node_to_link_node2.size() - 1 - i_node1 : i_node1;
            if (hinge_node_to_link_node1[i_node1] == hinge_node_to_link_node2[i_node2])
            {
                continue;
            }
            if (invalidity_graph(hinge_node_to_link_node1[i_node1],
                                 hinge_node_to_link_node2[i_node2]) == NO_ID)
            {
                EntityMultipleParts hinge_parts1;
                hinge_parts1.id = edge_info.node1.id;
                hinge_parts1.parts.push_back(
                    node_information[hinge_node_to_link_node1[i_node1]].parts[0].part);
                EntityMultipleParts hinge_parts2;
                hinge_parts2.id = edge_info.node2.id;
                hinge_parts2.parts.push_back(
                    node_information[hinge_node_to_link_node2[i_node2]].parts[0].part);
                EdgeInformation hinge_point_edge(hinge_parts1, hinge_parts2);
                hinge_point_edge.node1_id = hinge_node_to_link_node1[i_node1];
                hinge_point_edge.node2_id = hinge_node_to_link_node2[i_node2];
                add_edge_into_invalidity_graph(hinge_point_edge);
            }
        }
        // Remove this edge
        remove_edge_from_invalidity_graph(edge_info);

        // Reallocate other incident edges
        reallocate_other_incident_edges_after_split(edge_info, new_nodes1,
                                                    new_nodes2);

        // Remove old nodes
        node1_node_info.set_inactive();
        node2_node_info.set_inactive();
        connectivity_graph(edge_info.node1_id,
                           node1_node_info.ordered_boundary_nodes[0]) = false;
        connectivity_graph(edge_info.node1_id,
                           node1_node_info.ordered_boundary_nodes[1]) = false;
        connectivity_graph(edge_info.node2_id,
                           node2_node_info.ordered_boundary_nodes[0]) = false;
        connectivity_graph(edge_info.node2_id,
                           node2_node_info.ordered_boundary_nodes[1]) = false;
    }

    template <index_t DIMENSION>
    bool GeoModelTopologyRecovererBase<DIMENSION>::check_line_line_edge_compared_to_nodes(
        EdgeInformation &edge_info)
    {
        // Check parts concerned by the node and parts concerned by the edge
        // They must be the same to perform edge contraction
        auto node1_id = edge_info.node1.id;
        if (edge_info.node1.parts.size() > 1)
        {
            return false;
        }
        if (node_information[edge_info.node1_id].parts.size() > 1)
        {
            return false;
        }
        if (node_information[edge_info.node1_id].parts[0].id == node1_id)
        {
            if (node_information[edge_info.node1_id].parts[0].part != edge_info.node1.parts[0])
            {
                return false;
            }
        }
        auto node2_id = edge_info.node2.id;
        if (edge_info.node2.parts.size() > 1)
        {
            return false;
        }
        if (node_information[edge_info.node2_id].parts.size() > 1)
        {
            return false;
        }
        if (node_information[edge_info.node2_id].parts[0].id == node2_id)
        {
            if (node_information[edge_info.node2_id].parts[0].part != edge_info.node2.parts[0])
            {
                return false;
            }
        }
        return true;
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::reallocate_linked_edge_to_node_after_split(
        const index_t node_id,
        const std::vector<index_t> &node_linked_nodes,
        const std::vector<index_t> &new_nodes_after_split)
    {
        for (auto linked_node : node_linked_nodes)
        {
            index_t cur_edge_id = invalidity_graph(node_id, linked_node);
            if (cur_edge_id == NO_ID)
            {
                continue; // Current edge already removed
            }
            EdgeInformation cur_edge_info = edge_information[cur_edge_id];

            // Find split node and the other
            index_t split_node, other_node;
            EntityMultipleParts split_node_edge_info, other_node_edge_info;
            if (cur_edge_info.node1_id == node_id)
            {
                split_node = cur_edge_info.node1_id;
                split_node_edge_info = cur_edge_info.node1;
                other_node = cur_edge_info.node2_id;
                other_node_edge_info = cur_edge_info.node2;
            }
            else
            {
                other_node = cur_edge_info.node1_id;
                other_node_edge_info = cur_edge_info.node1;
                split_node = cur_edge_info.node2_id;
                split_node_edge_info = cur_edge_info.node2;
            }

            std::vector<EdgeInformation> reallocated_edges(
                new_nodes_after_split.size());
            std::vector<LineChunk> new_nodes_chunk(new_nodes_after_split.size());
            std::vector<LineChunk> other_nodes_corresponding_chunk(
                new_nodes_after_split.size());
            // Find begin/end corresponding coord
            for (auto i : RINGMesh::range(new_nodes_chunk.size()))
            {
                index_t split_line_id =
                    node_information[split_node].parts[0].id.index();
                index_t other_line_id =
                    node_information[other_node].parts[0].id.index();

                new_nodes_chunk[i] =
                    node_information[new_nodes_after_split[i]].parts[0].part;
                vecn<DIMENSION> input_begin =
                    PointOnLineCoordinates<DIMENSION>::absolute_coordinates(
                        new_nodes_chunk[i].begin,
                        geomodel_.line(split_line_id).mesh());
                vecn<DIMENSION> corresponding_begin_coord;
                index_t holding_edge;
                std::tie(holding_edge, corresponding_begin_coord, std::ignore) =
                    geomodel_.line(other_line_id).edge_aabb().closest_edge(input_begin);
                LineRelativeCoordinates corresponding_begin{PointOnLineCoordinates<
                    DIMENSION>::relative_coordinates(corresponding_begin_coord,
                                                     geomodel_.line(other_line_id).mesh(), holding_edge)};
                other_nodes_corresponding_chunk[i].begin = corresponding_begin;

                vecn<DIMENSION> input_end =
                    PointOnLineCoordinates<DIMENSION>::absolute_coordinates(
                        new_nodes_chunk[i].end,
                        geomodel_.line(split_line_id).mesh());
                vecn<DIMENSION> corresponding_end_coord;
                std::tie(holding_edge, corresponding_end_coord, std::ignore) =
                    geomodel_.line(other_line_id).edge_aabb().closest_edge(input_end);
                LineRelativeCoordinates corresponding_end{PointOnLineCoordinates<
                    DIMENSION>::relative_coordinates(corresponding_end_coord,
                                                     geomodel_.line(other_line_id).mesh(), holding_edge)};
                other_nodes_corresponding_chunk[i].end = corresponding_end;
            }

            for (auto part_id : RINGMesh::range(split_node_edge_info.parts.size()))
            {
                scar_assert(
                    other_node_edge_info.parts.size() == split_node_edge_info.parts.size());
                const auto &part_split_side = split_node_edge_info.parts[part_id];
                const auto &part_other_side = other_node_edge_info.parts[part_id];

                // Locate in which new node(s) begin and end belong
                index_t node_bearing_begin = NO_ID;
                index_t node_bearing_end = NO_ID;
                for (auto i : RINGMesh::range(new_nodes_chunk.size()))
                {
                    if (new_nodes_chunk[i].is_inside(part_split_side.begin))
                    {
                        node_bearing_begin = i;
                    }
                    if (new_nodes_chunk[i].is_inside(part_split_side.end))
                    {
                        node_bearing_end = i;
                    }
                    if (node_bearing_begin != NO_ID && node_bearing_end != NO_ID)
                    {
                        break;
                    }
                }
                scar_assert(node_bearing_begin != NO_ID);
                scar_assert(node_bearing_end != NO_ID);

                if (node_bearing_begin == node_bearing_end)
                {
                    // No need to split
                    reallocated_edges[node_bearing_begin].node1.parts.push_back(
                        part_split_side);
                    reallocated_edges[node_bearing_begin].node2.parts.push_back(
                        part_other_side);
                }
                else
                {
                    // Need to split the edge into several reallocated edges
                    { // For first node
                        LineChunk cur_reallocated_edge_split_side;
                        LineChunk cur_reallocated_edge_other_side;
                        cur_reallocated_edge_split_side.begin =
                            part_split_side.begin;
                        cur_reallocated_edge_other_side.begin =
                            part_other_side.begin;
                        cur_reallocated_edge_split_side.end =
                            new_nodes_chunk[node_bearing_begin].end;
                        cur_reallocated_edge_other_side.end =
                            other_nodes_corresponding_chunk[node_bearing_begin].end;
                        // Always check that the part is not null
                        if (cur_reallocated_edge_split_side.is_null() && cur_reallocated_edge_other_side.is_null())
                        {
                            continue;
                        }
                        if (cur_reallocated_edge_split_side.is_almost_null(0.5) && cur_reallocated_edge_other_side.is_almost_null(
                                                                                       0.5))
                        {
                            Logger::out("Graph edit",
                                        "Very small edges (due to split) had been removed... ");
                            continue;
                        }
                        reallocated_edges[node_bearing_begin].node1.parts.push_back(
                            cur_reallocated_edge_split_side);
                        reallocated_edges[node_bearing_begin].node2.parts.push_back(
                            cur_reallocated_edge_other_side);
                    }
                    // Nodes inbetween
                    for (auto p : RINGMesh::range(node_bearing_begin + 1,
                                                  node_bearing_end))
                    {
                        LineChunk cur_reallocated_edge_split_side =
                            new_nodes_chunk[p];
                        LineChunk cur_reallocated_edge_other_side =
                            other_nodes_corresponding_chunk[p];
                        if (cur_reallocated_edge_split_side.is_null() && cur_reallocated_edge_other_side.is_null())
                        {
                            continue;
                        }
                        if (cur_reallocated_edge_split_side.is_almost_null(0.5) && cur_reallocated_edge_other_side.is_almost_null(
                                                                                       0.5))
                        {
                            Logger::out("Graph edit",
                                        "Very small edges (due to split) had been removed... ");
                            continue;
                        }
                        reallocated_edges[p].node1.parts.push_back(
                            cur_reallocated_edge_split_side);
                        reallocated_edges[p].node2.parts.push_back(
                            cur_reallocated_edge_other_side);
                    }
                    { // For last node
                        LineChunk cur_reallocated_edge_split_side;
                        LineChunk cur_reallocated_edge_other_side;
                        cur_reallocated_edge_split_side.begin =
                            new_nodes_chunk[node_bearing_end].begin;
                        cur_reallocated_edge_other_side.begin =
                            other_nodes_corresponding_chunk[node_bearing_end].begin;
                        cur_reallocated_edge_split_side.end = part_split_side.end;
                        cur_reallocated_edge_other_side.end = part_other_side.end;
                        // Always check that the part is not null
                        if (cur_reallocated_edge_split_side.is_null() && cur_reallocated_edge_other_side.is_null())
                        {
                            continue;
                        }
                        if (cur_reallocated_edge_split_side.is_almost_null(0.5) && cur_reallocated_edge_other_side.is_almost_null(
                                                                                       0.5))
                        {
                            Logger::out("Graph edit",
                                        "Very small edges (due to split) had been removed... ");
                            continue;
                        }
                        reallocated_edges[node_bearing_end].node1.parts.push_back(
                            cur_reallocated_edge_split_side);
                        reallocated_edges[node_bearing_end].node2.parts.push_back(
                            cur_reallocated_edge_other_side);
                    }
                }
            }

            // Add reallocated edges after adding all needed info
            for (auto i : RINGMesh::range(reallocated_edges.size()))
            {
                if (reallocated_edges[i].node2.parts.empty())
                {
                    continue;
                }

                reallocated_edges[i].node1_id = new_nodes_after_split[i];
                reallocated_edges[i].node2_id = linked_node;
                reallocated_edges[i].node1.id =
                    node_information[split_node].parts[0].id;
                reallocated_edges[i].node2.id =
                    node_information[other_node].parts[0].id;

                add_edge_into_invalidity_graph(reallocated_edges[i]);
            }
            remove_edge_from_invalidity_graph(cur_edge_id);
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::reallocate_other_incident_edges_after_split(
        EdgeInformation &edge_info,
        const std::vector<index_t> &new_nodes1,
        const std::vector<index_t> &new_nodes2)
    {
        auto node1_linked_nodes = CorrelationMapAPI::directly_linked(
            invalidity_graph, edge_info.node1_id, NO_ID, false);
        if (!node1_linked_nodes.empty())
        {
            reallocate_linked_edge_to_node_after_split(edge_info.node1_id,
                                                       node1_linked_nodes, new_nodes1);
        }

        auto node2_linked_nodes = CorrelationMapAPI::directly_linked(
            invalidity_graph, edge_info.node2_id, NO_ID, false);
        if (!node2_linked_nodes.empty())
        {
            reallocate_linked_edge_to_node_after_split(edge_info.node2_id,
                                                       node2_linked_nodes, new_nodes2);
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::perform_edge_contractions()
    {
        bool modified = true;
        while (modified)
        {
            modified = false;
            for (auto &edge_info : edge_information)
            {
                if (!edge_info.is_active())
                {
                    continue;
                }
                scar_assert(
                    node_information[edge_info.node1_id].type == node_information[edge_info.node2_id].type);
                if (node_information[edge_info.node1_id].type == RINGMesh::corner_type_name_static())
                {
                    continue;
                }
                else if (node_information[edge_info.node1_id].type == RINGMesh::line_type_name_static())
                {
                    do_line_line_edge_contraction(edge_info);
                    modified = true;
                }
                else
                {
                    scar_assert_not_reached;
                }
            }
        }
        print_graph_nodes(true);
        print_graph_edges(false);
        std::vector<std::vector<index_t>> connected_components_with_edges =
            get_connected_components_with_edges();

        for (const auto &cc : connected_components_with_edges)
        {
            if (nodes_are_all_corners(cc))
            {
                EntityParts new_node;
                for (const auto &cc_node : cc)
                {
                    for (const auto &entity_part : node_information[cc_node].parts)
                    {
                        new_node.push_back(entity_part);
                    }
                }
                NodeInformation new_corner_node(RINGMesh::corner_type_name_static(),
                                                new_node);
                index_t new_node_id = add_node_into_invalidity_graph(
                    new_corner_node);
                for (const auto &cc_node : cc)
                {
                    std::vector<index_t> incident_line_nodes =
                        CorrelationMapAPI::directly_linked(connectivity_graph,
                                                           cc_node, false, false);
                    for (auto line_node_id : incident_line_nodes)
                    {
                        connectivity_graph(cc_node, line_node_id) = false;
                        connectivity_graph(new_node_id, line_node_id) = true;
                        auto &line_node_info = node_information[line_node_id];
                        scar_assert(line_node_info.is_active());
                        scar_assert(
                            line_node_info.ordered_boundary_nodes.size() == 2);
                        if (line_node_info.ordered_boundary_nodes[0] == cc_node)
                        {
                            line_node_info.ordered_boundary_nodes[0] = new_node_id;
                        }
                        if (line_node_info.ordered_boundary_nodes[1] == cc_node)
                        {
                            line_node_info.ordered_boundary_nodes[1] = new_node_id;
                        }
                        if (line_node_info.ordered_boundary_nodes[0] == line_node_info.ordered_boundary_nodes[1])
                        {
                            // The line has to be removed since its boundaries are merged
                            line_node_info.set_inactive();
                            connectivity_graph(
                                line_node_info.ordered_boundary_nodes[0],
                                line_node_id) = false;
                        }
                    }
                }

                // Remove all edges and set as inactive node and edges
                for (const auto &cc_node1 : cc)
                {
                    node_information[cc_node1].set_inactive();
                    for (const auto &cc_node2 : cc)
                    {
                        if (cc_node1 == cc_node2)
                        {
                            continue;
                        }
                        if (invalidity_graph(cc_node1, cc_node2) != NO_ID)
                        {
                            remove_edge_from_invalidity_graph(
                                invalidity_graph(cc_node1, cc_node2));
                        }
                    }
                }
                continue;
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::do_line_line_edge_contraction(
        EdgeInformation &edge_info)
    {
        auto &node1_info = node_information[edge_info.node1_id];
        auto &node2_info = node_information[edge_info.node2_id];
        scar_assert(node1_info.type == node2_info.type);
        scar_assert(node1_info.ordered_boundary_nodes.size() == 2);
        scar_assert(node2_info.ordered_boundary_nodes.size() == 2);
        NodeInformation new_node(node1_info.type);
        for (const auto &part : node1_info.parts)
        {
            new_node.parts.push_back(part);
        }
        // Remove node 1
        node1_info.set_inactive();
        connectivity_graph(edge_info.node1_id,
                           node1_info.ordered_boundary_nodes[0]) = false;
        connectivity_graph(edge_info.node1_id,
                           node1_info.ordered_boundary_nodes[1]) = false;

        for (const auto &part : node2_info.parts)
        {
            new_node.parts.push_back(part);
        }
        // Remove node 2
        node2_info.set_inactive();
        connectivity_graph(edge_info.node2_id,
                           node2_info.ordered_boundary_nodes[0]) = false;
        connectivity_graph(edge_info.node2_id,
                           node2_info.ordered_boundary_nodes[1]) = false;

        // Add connectivity for the new node
        new_node.ordered_boundary_nodes.push_back(
            node1_info.ordered_boundary_nodes[0]);
        new_node.ordered_boundary_nodes.push_back(
            node1_info.ordered_boundary_nodes[1]);

        // Add node into node information table
        index_t new_node_id = add_node_into_invalidity_graph(new_node);

        // Update connectivity graph
        connectivity_graph(new_node_id, new_node.ordered_boundary_nodes[0]) = true;
        connectivity_graph(new_node_id, new_node.ordered_boundary_nodes[1]) = true;

        // Remove edge
        remove_edge_from_invalidity_graph(edge_info);

        // Reallocate edges
        auto node1_linked_nodes = CorrelationMapAPI::directly_linked(
            invalidity_graph, edge_info.node1_id, NO_ID, false);
        for (auto linked_node_id : node1_linked_nodes)
        {
            auto &old_edge_info = edge_information[invalidity_graph(
                edge_info.node1_id, linked_node_id)];
            // Create new edge
            EdgeInformation new_edge_info = old_edge_info;
            if (old_edge_info.node1_id == edge_info.node1_id)
            {
                new_edge_info.node1_id = new_node_id;
            }
            else
            {
                new_edge_info.node2_id = new_node_id;
            }
            // Remove old edge
            remove_edge_from_invalidity_graph(old_edge_info);

            // Add new edge
            add_edge_into_invalidity_graph(new_edge_info);
        }

        auto node2_linked_nodes = CorrelationMapAPI::directly_linked(
            invalidity_graph, edge_info.node2_id, NO_ID, false);
        for (auto linked_node_id : node2_linked_nodes)
        {
            auto &old_edge_info = edge_information[invalidity_graph(
                edge_info.node2_id, linked_node_id)];
            // Create new edge
            EdgeInformation new_edge_info = old_edge_info;
            if (old_edge_info.node1_id == edge_info.node2_id)
            {
                new_edge_info.node1_id = new_node_id;
            }
            else
            {
                new_edge_info.node2_id = new_node_id;
            }
            // Remove old edge
            remove_edge_from_invalidity_graph(edge_info);

            // Add new edge (only if the edge does not already exist)
            if (invalidity_graph(new_edge_info.node1_id, new_edge_info.node2_id) == NO_ID)
            {
                add_edge_into_invalidity_graph(new_edge_info);
            }
        }
    }

    template <index_t DIMENSION>
    LineChunk GeoModelTopologyRecovererBase<DIMENSION>::get_smallest_common_line_part_with_linked_nodes(
        const EdgeInformation &edge)
    {
        LineChunk common_line_part = edge.node2.parts[0];
        auto linked_nodes = CorrelationMapAPI::directly_linked(invalidity_graph,
                                                               edge.node2_id, NO_ID, false);
        std::vector<index_t> linked_nodes_with_common_line_part;
        for (const auto &linked_node : linked_nodes)
        {
            if (linked_node == edge.node1_id)
            {
                continue;
            }
            auto &other_edge_info = edge_information[invalidity_graph(edge.node2_id,
                                                                      linked_node)];
            scar_assert(other_edge_info.is_active());
            if (other_edge_info.node1.id.type() == RINGMesh::corner_type_name_static())
            {
                auto &other_edge_line_parts = other_edge_info.node2.parts;
                for (index_t i = 0; i < other_edge_line_parts.size(); i++)
                {
                    auto &other_edge_part = other_edge_line_parts[i];
                    bool intersect;
                    LineChunk intersection;
                    std::tie(intersect, intersection) = other_edge_part.intersect(
                        common_line_part);
                    if (intersect)
                    {
                        add_element_if_not_in(linked_node,
                                              linked_nodes_with_common_line_part);
                        common_line_part = intersection;
                        // Remove the edge_part
                        other_edge_line_parts.erase(
                            other_edge_line_parts.begin() + i);
                        i--;
                    }
                }
            }
        }
        // If there is linked nodes in common, there must be an edge between these corners
        for (auto sharing_line_part_node : linked_nodes_with_common_line_part)
        {
            scar_unused(sharing_line_part_node);
            scar_assert(
                invalidity_graph(edge.node2_id, sharing_line_part_node) != NO_ID);
        }
        return common_line_part;
    }

    template <index_t DIMENSION>
    std::tuple<index_t, index_t> GeoModelTopologyRecovererBase<DIMENSION>::split_line_node_into_two_nodes(
        NodeInformation &line_node_info,
        const LineRelativeCoordinates &begin,
        const LineRelativeCoordinates &end,
        const LineChunk &corner_line_intersection,
        const index_t corner_node_id)
    {
        EntityParts line_part_before;
        scar_assert(end > corner_line_intersection.end);
        LineChunk part_before(begin, corner_line_intersection.end);
        line_part_before.emplace_back(line_node_info.parts[0].id, part_before);
        NodeInformation line_node_before(RINGMesh::line_type_name_static(),
                                         line_part_before);
        line_node_before.ordered_boundary_nodes.push_back(
            line_node_info.ordered_boundary_nodes[0]);
        line_node_before.ordered_boundary_nodes.push_back(corner_node_id);
        index_t new_line_node_before = add_node_into_invalidity_graph(
            line_node_before);

        EntityParts line_part_after;
        scar_assert(begin < corner_line_intersection.begin);
        LineChunk part_after(corner_line_intersection.begin, end);
        line_part_after.emplace_back(line_node_info.parts[0].id, part_after);
        NodeInformation line_node_after(RINGMesh::line_type_name_static(),
                                        line_part_after);
        line_node_after.ordered_boundary_nodes.push_back(corner_node_id);
        line_node_after.ordered_boundary_nodes.push_back(
            line_node_info.ordered_boundary_nodes[1]);
        index_t new_line_node_after = add_node_into_invalidity_graph(
            line_node_after);

        line_node_info.set_inactive();
        return std::make_tuple(new_line_node_before, new_line_node_after);
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::corner_line_edge_reallocation(
        EdgeInformation &reallocated_edge_info,
        index_t new_line_node_before,
        index_t new_line_node_after)
    {
        // C-L edge
        EdgeInformation edge_before_info;
        edge_before_info.node1_id = reallocated_edge_info.node1_id;
        edge_before_info.node1 = reallocated_edge_info.node1;
        edge_before_info.node2_id = new_line_node_before;
        EntityMultipleParts edge_before_line_node;
        edge_before_line_node.id = reallocated_edge_info.node2.id;
        EdgeInformation edge_after_info;
        edge_after_info.node1_id = reallocated_edge_info.node1_id;
        edge_after_info.node1 = reallocated_edge_info.node1;
        edge_after_info.node2_id = new_line_node_after;
        EntityMultipleParts edge_after_line_node;
        edge_after_line_node.id = reallocated_edge_info.node2.id;
        for (const auto &other_part : reallocated_edge_info.node2.parts)
        {
            const LineChunk &part_before =
                node_information[new_line_node_before].parts[0].part;
            if (other_part.end < part_before.end)
            {
                // Line before
                edge_before_line_node.parts.push_back(other_part);
            }
            else
            {
                // Line after
                const LineChunk &part_after =
                    node_information[new_line_node_after].parts[0].part;
                scar_unused(part_after);
                scar_assert(other_part.begin > part_after.begin);
                edge_after_line_node.parts.push_back(other_part);
            }
        }
        edge_before_info.node2 = std::move(edge_before_line_node);
        edge_after_info.node2 = std::move(edge_after_line_node);

        remove_edge_from_invalidity_graph(reallocated_edge_info);
        if (!edge_before_info.node2.parts.empty())
        {
            add_edge_into_invalidity_graph(edge_before_info);
        }
        if (!edge_after_info.node2.parts.empty())
        {
            add_edge_into_invalidity_graph(edge_after_info);
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::line_line_edge_reallocation(
        EdgeInformation &reallocated_edge_info,
        index_t previous_line_node,
        const EdgeInformation &split_edge,
        index_t new_line_node_before,
        index_t new_line_node_after)
    {
        EntityMultipleParts previous_line_parts;
        EntityMultipleParts other_line_parts;
        if (reallocated_edge_info.node1_id == previous_line_node)
        {
            previous_line_parts = reallocated_edge_info.node1;
            other_line_parts = reallocated_edge_info.node2;
        }
        else if (reallocated_edge_info.node2_id == previous_line_node)
        {
            previous_line_parts = reallocated_edge_info.node2;
            other_line_parts = reallocated_edge_info.node1;
        }
        //@todo check if there are 1 or 2 boundary corner linked with
        const auto &other_line = geomodel_.line(other_line_parts.id.index());
        const auto &other_line_boundary0_id = other_line.boundary(0).gmme();
        const auto &other_line_boundary1_id = other_line.boundary(1).gmme();
        bool other_line_boundary0_linked{false};
        bool other_line_boundary1_linked{false};
        for (const auto &cur_edge : edge_information)
        {
            if (!cur_edge.is_active() || cur_edge.node1.id.type() != RINGMesh::corner_type_name_static())
            {
                // Cur_edge must be active and be of type Corner-XX (?)
                continue;
            }
            if (cur_edge.node1.id.index() != split_edge.node1.id.index() && cur_edge.node2.id.index() != split_edge.node1.id.index())
            {
                // Cur_edge must concern the corner of the split CL edge
                continue;
            }
            if (cur_edge.node1.id == other_line_boundary0_id || cur_edge.node2.id == other_line_boundary0_id)
            {
                other_line_boundary0_linked = true;
            }
            if (cur_edge.node1.id == other_line_boundary1_id || cur_edge.node2.id == other_line_boundary1_id)
            {
                other_line_boundary1_linked = true;
            }
        }

        EdgeInformation edge_before_info;
        // New nodes corresponding to split lines have index greater then the other line
        // @todo Not sure of the following lines
        edge_before_info.node1_id =
            reallocated_edge_info.node1_id == previous_line_node ? reallocated_edge_info.node2_id : reallocated_edge_info.node1_id;
        edge_before_info.node2_id = new_line_node_before;

        EntityMultipleParts edge_before_new_line_node;
        edge_before_new_line_node.id = split_edge.node2.id;
        EntityMultipleParts edge_before_other_line_node;
        edge_before_other_line_node.id = other_line.gmme();

        EdgeInformation edge_after_info;
        edge_after_info.node1_id = edge_before_info.node1_id;
        edge_after_info.node2_id = new_line_node_after;

        EntityMultipleParts edge_after_new_line_node;
        edge_after_new_line_node.id = split_edge.node2.id;
        EntityMultipleParts edge_after_other_line_node;
        edge_after_other_line_node.id = other_line.gmme();

        scar_assert(
            previous_line_parts.parts.size() == other_line_parts.parts.size());
        for (index_t i = 0; i < previous_line_parts.parts.size(); i++)
        {
            const auto &prev_line_part = previous_line_parts.parts[i];
            const auto &other_line_part = other_line_parts.parts[i];
            const LineChunk &part_before =
                node_information[new_line_node_before].parts[0].part;
            const LineChunk &part_after =
                node_information[new_line_node_after].parts[0].part;
            if (prev_line_part.begin == 0. && other_line_boundary0_linked)
            {
                if (prev_line_part.end < part_before.end)
                {
                    edge_before_new_line_node.parts.push_back(prev_line_part);
                    edge_before_other_line_node.parts.push_back(other_line_part);
                }
                else
                {
                    edge_after_new_line_node.parts.push_back(prev_line_part);
                    edge_after_other_line_node.parts.push_back(other_line_part);
                }
            }
            else if (prev_line_part.end == other_line.nb_mesh_elements() && other_line_boundary1_linked)
            {
                if (prev_line_part.begin > part_after.begin)
                {
                    edge_after_new_line_node.parts.push_back(prev_line_part);
                    edge_after_other_line_node.parts.push_back(other_line_part);
                }
                else
                {
                    edge_before_new_line_node.parts.push_back(prev_line_part);
                    edge_before_other_line_node.parts.push_back(other_line_part);
                }
            }
            else
            {
                if (prev_line_part.end < part_before.end)
                {
                    edge_before_new_line_node.parts.push_back(prev_line_part);
                    edge_before_other_line_node.parts.push_back(other_line_part);
                }
                else if (prev_line_part.begin > part_after.begin)
                {
                    edge_after_new_line_node.parts.push_back(prev_line_part);
                    edge_after_other_line_node.parts.push_back(other_line_part);
                }
                else
                {
                    // Must split other_line_part and reallocation
                    // the part before and after
                    auto part_before_end_point =
                        PointOnLineCoordinates<DIMENSION>::absolute_coordinates(
                            part_before.end,
                            geomodel_.line(previous_line_parts.id.index()).mesh());
                    auto part_after_begin_point =
                        PointOnLineCoordinates<DIMENSION>::absolute_coordinates(
                            part_after.begin,
                            geomodel_.line(previous_line_parts.id.index()).mesh());
                    // @todo This part does not handled anisotropic primitive
                    vecn<DIMENSION> part_before_end_point_equivalent;
                    index_t be_point_edge;
                    std::tie(be_point_edge, part_before_end_point_equivalent,
                             std::ignore) = other_line.edge_aabb().closest_edge(part_before_end_point);
                    vecn<DIMENSION> part_after_begin_point_equivalent;
                    index_t ab_point_edge;
                    std::tie(ab_point_edge, part_after_begin_point_equivalent,
                             std::ignore) = other_line.edge_aabb().closest_edge(part_after_begin_point);
                    LineRelativeCoordinates other_part_before_end{
                        PointOnLineCoordinates<DIMENSION>::relative_coordinates(
                            part_before_end_point_equivalent, other_line.mesh(),
                            be_point_edge)};
                    LineRelativeCoordinates other_part_after_begin{
                        PointOnLineCoordinates<DIMENSION>::relative_coordinates(
                            part_after_begin_point_equivalent, other_line.mesh(),
                            ab_point_edge)};

                    if (other_part_after_begin > other_part_before_end)
                    {
                        DEBUG(
                            "WARNING: the other line seems to be oriented in the other direction");
                    }
                    //@todo maybe a mistake if the line is oriented in the other direction
                    bool intersect_before;
                    LineChunk intersection_before;
                    std::tie(intersect_before, intersection_before) =
                        part_before.intersect(prev_line_part);
                    if (intersect_before)
                    {
                        edge_before_new_line_node.parts.push_back(
                            intersection_before);
                        edge_before_other_line_node.parts.push_back(
                            LineChunk(other_line_part.begin,
                                      other_part_before_end));
                    }
                    bool intersect_after;
                    LineChunk intersection_after;
                    std::tie(intersect_after, intersection_after) =
                        part_after.intersect(prev_line_part);
                    if (intersect_after)
                    {
                        edge_after_new_line_node.parts.push_back(
                            intersection_after);
                        edge_after_other_line_node.parts.push_back(
                            LineChunk(other_part_after_begin,
                                      other_line_part.end));
                    }
                }
            }
        }
        edge_before_info.node1 = std::move(edge_before_other_line_node);
        edge_before_info.node2 = std::move(edge_before_new_line_node);
        edge_after_info.node1 = std::move(edge_after_other_line_node);
        edge_after_info.node2 = std::move(edge_after_new_line_node);

        // Remove old edges
        remove_edge_from_invalidity_graph(reallocated_edge_info);
        if (!edge_before_info.node2.parts.empty())
        {
            add_edge_into_invalidity_graph(edge_before_info);
        }
        if (!edge_after_info.node2.parts.empty())
        {
            add_edge_into_invalidity_graph(edge_after_info);
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::incident_edges_reallocation(
        index_t node,
        const EdgeInformation &edge,
        index_t new_line_node_before,
        index_t new_line_node_after)
    {
        auto linked_nodes = CorrelationMapAPI::directly_linked(invalidity_graph,
                                                               node, NO_ID, false);
        // @todo Use a std::deque for new edges to avoid the modification of
        // edge_information vector (enclosing iterative loop)
        for (const auto &linked_node : linked_nodes)
        {
            auto &other_edge_info = edge_information[invalidity_graph(node,
                                                                      linked_node)];
            scar_assert(other_edge_info.is_active());
            if (other_edge_info.node1.id.type() == RINGMesh::corner_type_name_static())
            {
                corner_line_edge_reallocation(other_edge_info, new_line_node_before,
                                              new_line_node_after);
            }
            else
            {
                scar_assert(
                    other_edge_info.node1.id.type() == RINGMesh::line_type_name_static() && other_edge_info.node2.id.type() == RINGMesh::line_type_name_static());
                // L-L edge
                line_line_edge_reallocation(other_edge_info, node, edge,
                                            new_line_node_before, new_line_node_after);
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::remove_corner_line_edges()
    {
        for (index_t i = 0; i < edge_information.size(); i++)
        {
            auto &edge = edge_information[i];
            if (!edge.is_active())
            {
                continue;
            }
            if (edge.node1.id.type() == RINGMesh::corner_type_name_static() && edge.node2.id.type() == RINGMesh::line_type_name_static())
            {
                // Split the line in several parts (so in several nodes)

                auto line_node_id = edge.node2_id;
                auto &line_node_info = node_information[line_node_id];
                scar_assert(line_node_info.is_active())
                    scar_assert(line_node_info.parts.size() == 1)
                        scar_assert(
                            line_node_info.parts[0].id.type() == RINGMesh::line_type_name_static());
                auto corner_node_id = edge.node1_id;
                auto &corner_node_info = node_information[corner_node_id];
                scar_assert(corner_node_info.is_active())
                    scar_assert(
                        corner_node_info.parts[0].id.type() == RINGMesh::corner_type_name_static());

                auto linked_nodes = CorrelationMapAPI::directly_linked(
                    invalidity_graph, line_node_id, NO_ID, false);

                auto begin = line_node_info.parts[0].part.begin;
                auto end = line_node_info.parts[0].part.end;
                // Each part contains the line_portion included in
                // the corner exclusion zone
                const auto &first_edge_part = edge.node2.parts[0];

                // Check if there is another edge sharing this edge_part
                LineChunk common_line_part =
                    get_smallest_common_line_part_with_linked_nodes(edge);

                // Cut the line into two parts
                index_t new_line_node_before;
                index_t new_line_node_after;
                std::tie(new_line_node_before, new_line_node_after) =
                    split_line_node_into_two_nodes(line_node_info, begin, end,
                                                   common_line_part, corner_node_id);

                // Add the line part information on only one corner (all corners will be merged after)
                corner_node_info.parts.emplace_back(edge.node2.id,
                                                    first_edge_part);

                // Current edge (C-L)
                // Erase the processed line_part
                edge.node2.parts.erase(edge.node2.parts.begin());

                // Reallocation of the edges
                incident_edges_reallocation(line_node_id, edge,
                                            new_line_node_before, new_line_node_after);

                // Remove the old line node
                line_node_info.set_inactive();

                // Reallocation of connectivity graph edge
                scar_assert(line_node_info.ordered_boundary_nodes.size() == 2);
                connectivity_graph(line_node_id,
                                   line_node_info.ordered_boundary_nodes[0]) = false;
                connectivity_graph(line_node_id,
                                   line_node_info.ordered_boundary_nodes[1]) = false;

                connectivity_graph(corner_node_id, new_line_node_before) = true;
                connectivity_graph(line_node_info.ordered_boundary_nodes[0],
                                   new_line_node_before) = true;

                connectivity_graph(corner_node_id, new_line_node_after) = true;
                connectivity_graph(line_node_info.ordered_boundary_nodes[1],
                                   new_line_node_after) = true;
            }
        }
    }

    template <index_t DIMENSION>
    std::tuple<bool, EdgeInformation> GeoModelTopologyRecovererBase<DIMENSION>::analyze_line_pair_exclusion_zones(
        const RINGMesh::Line<DIMENSION> &line1,
        const RINGMesh::Line<DIMENSION> &line2)
    {
        PointOnLineProjectionTool<DIMENSION> projector(geomodel_);

        // If only whole lines can me merged else their boundaries should be equal
        // or linked by invalidity edges
        if (enum_contains(geological_rules, GeologicalRule::MERGE_WHOLE_LINES))
        {
            index_t l1b0 = line1.boundary(0).index();
            index_t l1b1 = line1.boundary(1).index();
            index_t l2b0 = line2.boundary(0).index();
            index_t l2b1 = line2.boundary(1).index();
            if (l1b0 != l2b0 && l1b0 != l2b1 && invalidity_graph(l1b0, l2b0) == NO_ID && invalidity_graph(l1b0, l2b1) == NO_ID)
            {
                return std::make_tuple(false, EdgeInformation());
            }
            if (l1b1 != l2b0 && l1b1 != l2b1 && invalidity_graph(l1b1, l2b0) == NO_ID && invalidity_graph(l1b1, l2b1) == NO_ID)
            {
                return std::make_tuple(false, EdgeInformation());
            }
        }

        std::map<double, double> collision_l1_to_l2;

        for (auto l1_v_id : RINGMesh::range(line1.nb_vertices()))
        {
            auto l1_v_id_exclusion = get_exclusion_shape(shape_type_, line1.mesh(),
                                                         l1_v_id);

            // Exclusion shapes should be NSpheres
            scar_assert(shape_type_ == ConvexShapeType::NSphere);
            auto projection_info = projector.find_projection(
                line1.vertex(l1_v_id), line2.index());
            // Quit if projection is the corner which is common
            if (projection_info.rcoord == 0. || projection_info.rcoord == line2.nb_vertices())
            {
                if (projector.find_projection(projection_info.coord, line1.index()).distance == 0.)
                {
                    continue;
                }
            }

            auto l2_exclusion = get_exclusion_shape(shape_type_, line2.mesh(),
                                                    projection_info.rcoord);
            bool projector_intersect = are_convex_intersecting(*l1_v_id_exclusion,
                                                               *l2_exclusion);
            if (projector_intersect)
            {
                collision_l1_to_l2[l1_v_id] = projection_info.rcoord;
            }
        }

        for (auto l2_v_id : RINGMesh::range(line2.nb_vertices()))
        {
            auto l2_v_id_exclusion = get_exclusion_shape(shape_type_, line2.mesh(),
                                                         l2_v_id);

            // Exclusion shapes should be NSpheres
            scar_assert(shape_type_ == ConvexShapeType::NSphere);
            auto projection_info = projector.find_projection(
                line2.vertex(l2_v_id), line1.index());

            // Quit if projection if corner which is common
            if (projection_info.rcoord == 0. || projection_info.rcoord == line1.nb_vertices())
            {
                if (projector.find_projection(projection_info.coord, line2.index()).distance == 0.)
                {
                    continue;
                }
            }

            auto l1_exclusion = get_exclusion_shape(shape_type_, line1.mesh(),
                                                    projection_info.rcoord);
            bool projector_intersect = are_convex_intersecting(*l2_v_id_exclusion,
                                                               *l1_exclusion);
            if (projector_intersect)
            {
                collision_l1_to_l2[projection_info.rcoord] = (double)l2_v_id;
            }
        }

        // Process for angle if lines are adjacent
        std::set<double> l1_in_collisions;
        std::set<double> l2_in_collisions;
        for (const auto &p : collision_l1_to_l2)
        {
            l1_in_collisions.insert(p.first);
            l2_in_collisions.insert(p.second);
        }

        if (line1.boundary(0).index() == line2.boundary(0).index())
        {
            collision_l1_to_l2[0] = 0;
        }
        else
        {
            if (invalidity_graph(line1.boundary(0).index(),
                                 line2.boundary(0).index()) != NO_ID)
            {
                collision_l1_to_l2[0] = 0;
            }
        }

        if (line1.boundary(0).index() == line2.boundary(1).index())
        {
            collision_l1_to_l2[0] = line2.nb_vertices() - 1;
        }
        else
        {
            if (invalidity_graph(line1.boundary(0).index(),
                                 line2.boundary(1).index()) != NO_ID)
            {
                collision_l1_to_l2[0] = line2.nb_vertices() - 1;
            }
        }

        if (line1.boundary(1).index() == line2.boundary(0).index())
        {
            collision_l1_to_l2[line1.nb_vertices() - 1] = 0;
        }
        else
        {
            if (invalidity_graph(line1.boundary(1).index(),
                                 line2.boundary(0).index()) != NO_ID)
            {
                collision_l1_to_l2[line1.nb_vertices() - 1] = 0;
            }
        }

        if (line1.boundary(1).index() == line2.boundary(1).index())
        {
            collision_l1_to_l2[line1.nb_vertices() - 1] = line2.nb_vertices() - 1;
        }
        else
        {
            if (invalidity_graph(line1.boundary(1).index(),
                                 line2.boundary(1).index()) != NO_ID)
            {
                collision_l1_to_l2[line1.nb_vertices() - 1] = line2.nb_vertices() - 1;
            }
        }

        LineChunks l1_chunks;
        LineChunks l2_chunks;
        LineChunk cur_l1_chunk{};
        LineChunk cur_l2_chunk{};
        bool init_chunks = true;

        // l1 points are ordered but not necessarily l2_points
        for (const auto &p : collision_l1_to_l2)
        {
            LineRelativeCoordinates cur_l1_pt(p.first);
            LineRelativeCoordinates cur_l2_pt(p.second);
            if (init_chunks)
            {
                cur_l1_chunk.begin = cur_l1_pt;
                cur_l1_chunk.end = cur_l1_pt;
                cur_l2_chunk.begin = cur_l2_pt;
                cur_l2_chunk.end = cur_l2_pt;
                init_chunks = false;
                continue;
            }

            if (cur_l1_pt - cur_l1_chunk.end <= 1. && (cur_l1_pt == std::floor(cur_l1_pt) || std::floor(cur_l1_pt) == std::floor(cur_l1_chunk.end)))
            {
                cur_l1_chunk.end = cur_l1_pt;
                cur_l2_chunk.begin = std::min(cur_l2_chunk.begin, cur_l2_pt);
                cur_l2_chunk.end = std::max(cur_l2_chunk.end, cur_l2_pt);
            }
            else
            {
                // Check lengthes of line_chunks, if not null add them, else delete them
                init_chunks = true;
                if (cur_l1_chunk.is_almost_null(global_epsilon) && (cur_l1_chunk.begin == 0 || cur_l1_chunk.end == line1.nb_vertices() - 1))
                {
                    continue;
                }
                if (cur_l2_chunk.is_almost_null(global_epsilon) && (cur_l2_chunk.begin == 0 || cur_l2_chunk.end == line2.nb_vertices() - 1))
                {
                    continue;
                }

                l1_chunks.push_back(cur_l1_chunk);
                l2_chunks.push_back(cur_l2_chunk);
            }
        }
        // Last chunk
        if (!(cur_l1_chunk.is_almost_null(global_epsilon) && (cur_l1_chunk.begin == 0 || cur_l1_chunk.end == line1.nb_vertices() - 1)) && !(cur_l2_chunk.is_almost_null(global_epsilon) && (cur_l2_chunk.begin == 0 || cur_l2_chunk.end == line2.nb_vertices() - 1)))
        {
            l1_chunks.push_back(cur_l1_chunk);
            l2_chunks.push_back(cur_l2_chunk);
        }

        EntityMultipleParts l1_edge_info = {line1.gmme(), l1_chunks};
        EntityMultipleParts l2_edge_info = {line2.gmme(), l2_chunks};

        scar_assert(l1_chunks.size() == l2_chunks.size());

        if (l1_chunks.empty())
        {
            return std::make_tuple(false, EdgeInformation());
        }

        if (l1_chunks.size() == 1)
        {
            if (l1_chunks[0].begin == 0)
            {
                auto linked_nodes = CorrelationMapAPI::directly_linked(
                    invalidity_graph, line1.boundary(0).index(), NO_ID, false);
                if (RINGMesh::contains(linked_nodes,
                                       geomodel_.nb_corners() + line1.index()))
                {
                    return std::make_tuple(false, EdgeInformation());
                }
            }
            if (l1_chunks[0].end == line1.nb_vertices() - 1)
            {
                auto linked_nodes = CorrelationMapAPI::directly_linked(
                    invalidity_graph, line1.boundary(1).index(), NO_ID, false);
                if (RINGMesh::contains(linked_nodes,
                                       geomodel_.nb_corners() + line1.index()))
                {
                    return std::make_tuple(false, EdgeInformation());
                }
            }
        }

        // If only whole lines can be merged, we check that the whole
        // lines have intersection exclusion zones
        if (enum_contains(geological_rules, GeologicalRule::MERGE_WHOLE_LINES))
        {
            if (l1_chunks.size() > 1)
            {
                return std::make_tuple(false, EdgeInformation());
            }
            if (l1_chunks[0].begin != 0 || l1_chunks[0].end != line1.nb_mesh_elements())
            {
                return std::make_tuple(false, EdgeInformation());
            }
            if (l2_chunks[0].begin != 0 || l2_chunks[0].end != line2.nb_mesh_elements())
            {
                return std::make_tuple(false, EdgeInformation());
            }
        }
        return std::make_tuple(true, EdgeInformation(l1_edge_info, l2_edge_info));
    }

    template <index_t DIMENSION>
    void GeoModelTopologyRecovererBase<DIMENSION>::perform_graph_based_repair()
    {
        initialize_graph_nodes();
        build_connectivity_graph();
        build_invalidity_graph();

        output_invalid_areas();
        output_graphs();

        if (DIMENSION == 3)
        {
            Logger::err("3D",
                        "Graph based repair and simplification not implemented in 3D.");
            return;
        }

        if (!enum_contains(geological_rules, GeologicalRule::CONSTANT_TOPOLOGY))
        {
            remove_edges_from_invalidity_graph();
        }
        if (verbose_)
        {
            print_graph_nodes(true);
            print_connectivity_graph_edges();
        }
    }

    template class scar_api GeoModelTopologyRecoverer<2>;
    template class scar_api GeoModelTopologyRecovererBase<2>;

    template class scar_api GeoModelTopologyRecoverer<3>;
    template class scar_api GeoModelTopologyRecovererBase<3>;
}
