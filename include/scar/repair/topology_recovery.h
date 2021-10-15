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

#include <ringmesh/basic/box.h>
#include <ringmesh/geomodel/core/geomodel.h>

#include <scar/repair/additional_rules.h>
#include <scar/tools/convex_shape.h>
#include <scar/tools/distance.h>

#include <ringpcl/correlation_map.h>

namespace SCAR
{
    struct EntityPart
    {
        EntityPart(const RINGMesh::gmme_id &gmme_id, LineChunk entity_part)
            : id(gmme_id), part(entity_part)
        {
        }
        RINGMesh::gmme_id id;
        LineChunk part;
    };

    using EntityParts = std::vector<EntityPart>;

    struct EntityMultipleParts
    {
        RINGMesh::gmme_id id;
        LineChunks parts;
    };

    struct NodeInformation
    {
        NodeInformation(RINGMesh::MeshEntityType node_type)
            : type(node_type)
        {
        }
        NodeInformation(
            RINGMesh::MeshEntityType node_type,
            EntityParts entity_parts)
            : type(node_type), parts(std::move(entity_parts))
        {
        }

        NodeInformation(const NodeInformation &rhs)
            : type(rhs.type),
              parts(rhs.parts),
              ordered_boundary_nodes(rhs.ordered_boundary_nodes)
        {
        }
        bool is_active() const
        {
            return active;
        }

        void set_inactive()
        {
            active = false;
        }
        RINGMesh::MeshEntityType type;
        EntityParts parts;
        std::vector<index_t> ordered_boundary_nodes;

    private:
        bool active{true};
    };

    enum struct EdgeType
    {
        CornerCorner,
        CornerLine,
        CornerSurface,
        LineLine,
        LineSurface,
        SurfaceSurface
    };

    struct EdgeInformation
    {
        EdgeInformation()
        {
        }
        EdgeInformation(
            EntityMultipleParts entity_part1,
            EntityMultipleParts entity_part2)
            : node1(std::move(entity_part1)), node2(std::move(entity_part2))
        {
        }

        bool is_active() const
        {
            return active;
        }
        void set_inactive()
        {
            active = false;
        }

        EntityMultipleParts node1;
        EntityMultipleParts node2;
        index_t node1_id{NO_ID};
        index_t node2_id{NO_ID};

    private:
        bool active{true};
    };

    /*!
     * @brief This class aims to analysis and recover the topology of a simplified or
     * repaired version of a given model according to tolerance and angle information.
     */
    template <index_t DIMENSION>
    class scar_api GeoModelTopologyRecovererBase
    {
        scar_disable_copy(GeoModelTopologyRecovererBase);
        scar_template_assert_2d_or_3d(DIMENSION);

    public:
        virtual ~GeoModelTopologyRecovererBase();

        void add_geological_rule(GeologicalRule new_rule)
        {
            geological_rules = geological_rules | new_rule;
        }

        void perform_graph_based_repair();

        const RINGMesh::GeoModel<DIMENSION> &geomodel() const
        {
            return geomodel_;
        }

        void build_connectivity_graph();
        virtual void build_graph_edges_from_geomodel_connectivity();

        void build_invalidity_graph();

        index_t add_edge_into_invalidity_graph(const EdgeInformation &new_edge);
        void remove_edge_from_invalidity_graph(index_t edge_id);
        void remove_edge_from_invalidity_graph(EdgeInformation &edge);
        index_t add_node_into_invalidity_graph(const NodeInformation &new_node);

        virtual void initialize_graph_nodes();
        virtual void build_graph_edges_from_geomodel_analysis();

        void output_invalid_areas_base(
            RINGMesh::GeoModel<DIMENSION> &invalidity_model) const;

        EdgeType get_edge_type(const EdgeInformation &edge) const;

    private:
        bool respect_geological_rules(
            const RINGMesh::Corner<DIMENSION> &free_horizon_corner,
            const RINGMesh::Corner<DIMENSION> &other_corner) const;
        bool respect_geological_rules(
            const RINGMesh::Corner<DIMENSION> &free_horizon_corner,
            const RINGMesh::Line<DIMENSION> &line) const;

        void find_nearest_entity_and_add_graph_edge(
            const RINGMesh::Corner<DIMENSION> &corner);

        void remove_edges_from_invalidity_graph();
        void remove_corner_line_edges();

        std::vector<std::vector<index_t>> get_connected_components_with_edges();

        LineChunk get_smallest_common_line_part_with_linked_nodes(
            const EdgeInformation &edge);

        std::tuple<index_t, index_t> split_line_node_into_two_nodes(
            NodeInformation &line_node_info,
            const LineRelativeCoordinates &begin,
            const LineRelativeCoordinates &end,
            const LineChunk &corner_line_intersection,
            const index_t corner_node_id);

        void incident_edges_reallocation(
            index_t node,
            const EdgeInformation &edge,
            index_t new_line_node_before,
            index_t new_line_node_after);
        void corner_line_edge_reallocation(
            EdgeInformation &reallocated_edge_info,
            index_t new_line_node_before,
            index_t new_line_node_after);
        void line_line_edge_reallocation(
            EdgeInformation &reallocated_edge_info,
            index_t previous_line_node,
            const EdgeInformation &split_edge,
            index_t new_line_node_before,
            index_t new_line_node_after);

        void perform_line_line_edge_split(EdgeInformation &edge_info);
        bool check_line_line_edge_compared_to_nodes(EdgeInformation &edge_info);
        void reallocate_other_incident_edges_after_split(
            EdgeInformation &edge_info,
            const std::vector<index_t> &new_nodes1,
            const std::vector<index_t> &new_nodes2);
        void reallocate_linked_edge_to_node_after_split(
            const index_t node_id,
            const std::vector<index_t> &node_linked_nodes,
            const std::vector<index_t> &new_nodes_after_split);
        std::vector<LineRelativeCoordinates> get_node_line_hinge_points(
            const NodeInformation &line_node,
            const EntityMultipleParts &edge_info_for_node);
        std::tuple<std::vector<index_t>, std::vector<index_t>,
                   std::vector<index_t>>
        split_line_node_at_hinge_points(
            const EntityMultipleParts &edge_info_for_node,
            const NodeInformation &node_info,
            const std::vector<LineRelativeCoordinates> &hinge_points);

        void split_line_line_edges_and_nodes();

        void perform_edge_contractions();
        void do_line_line_edge_contraction(EdgeInformation &edge_info);

        void analyze_corner_corner();
        void analyze_corner_line();
        void analyze_line_line();

        std::tuple<bool, EdgeInformation> analyze_line_pair_exclusion_zones(
            const RINGMesh::Line<DIMENSION> &line1,
            const RINGMesh::Line<DIMENSION> &line2);

        bool nodes_are_all_corners(
            std::vector<index_t> connected_component) const;
        bool nodes_are_all_lines(std::vector<index_t> connected_component) const;

        bool are_line_line_reverse(const EdgeInformation &edge) const;

        void print_graph_edges(bool verbose = false) const;
        void print_graph_nodes(bool verbose = false) const;

        void print_connectivity_graph_edges() const;

        virtual void output_invalid_areas() const = 0;
        void output_graphs() const;

    protected:
        GeoModelTopologyRecovererBase(
            RINGMesh::GeoModel<DIMENSION> &geomodel,
            ConvexShapeType form_type,
            bool verbose);

    public:
        CorrelationMap<index_t> invalidity_graph;
        std::vector<NodeInformation> node_information;
        std::vector<EdgeInformation> edge_information;

        CorrelationMap<bool> connectivity_graph;
        GeologicalRule geological_rules{GeologicalRule::EMPTY};

    protected:
        const RINGMesh::GeoModel<DIMENSION> &geomodel_;
        ConvexShapeType shape_type_;
        bool verbose_;
        std::set<index_t> invalid_vertices;
        std::set<index_t> invalid_segments;
    };

    ALIAS_2D_AND_3D(GeoModelTopologyRecovererBase);

    template <index_t DIMENSION>
    class scar_api GeoModelTopologyRecoverer final : public GeoModelTopologyRecovererBase<
                                                         DIMENSION>
    {
    public:
        GeoModelTopologyRecoverer(
            RINGMesh::GeoModel<DIMENSION> &geomodel,
            ConvexShapeType form_type,
            bool verbose);
    };

    template <>
    class scar_api GeoModelTopologyRecoverer<3> final : public GeoModelTopologyRecovererBase<
                                                            3>
    {
    public:
        GeoModelTopologyRecoverer(
            RINGMesh::GeoModel3D &geomodel,
            ConvexShapeType form_type,
            bool verbose);

    private:
        void initialize_graph_nodes() override;
        void build_graph_edges_from_geomodel_connectivity() override;
        void build_graph_edges_from_geomodel_analysis() override;

        void analyze_corner_surface();
        void analyze_line_surface();
        void analyze_surface_surface();

        void output_invalid_areas() const override;

    private:
        std::set<std::pair<index_t, index_t>> invalid_triangles;
    };

    template <>
    class scar_api GeoModelTopologyRecoverer<2> final : public GeoModelTopologyRecovererBase<
                                                            2>
    {
    public:
        GeoModelTopologyRecoverer(
            RINGMesh::GeoModel2D &geomodel,
            ConvexShapeType form_type,
            bool verbose);

    private:
        void output_invalid_areas() const override;
    };
    ALIAS_2D_AND_3D(GeoModelTopologyRecoverer);
}
