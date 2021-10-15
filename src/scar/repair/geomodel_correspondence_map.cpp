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

#include <scar/repair/geomodel_correspondence_map.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <ringpcl/correlation_map.h>

#include <scar/repair/topology_recovery.h>
#include <scar/tools/utils.h>

namespace
{

    using namespace SCAR;

    template <index_t DIMENSION>
    void build_mesh_entities(
        const RINGMesh::MeshEntityType &type,
        const GeoModelTopologyRecoverer<DIMENSION> &topology_recoverer,
        std::vector<index_t> &nodes_to_mesh_entity_ids,
        RINGMesh::GeoModel<DIMENSION> &output_geomodel)
    {
        scar_assert(output_geomodel.nb_mesh_entities(type) == 0);
        index_t nb_mesh_entities{0};
        for (auto node_id : RINGMesh::range(
                 topology_recoverer.node_information.size()))
        {
            const auto &node = topology_recoverer.node_information[node_id];
            if (!node.is_active())
            {
                continue;
            }
            if (node.type != type)
            {
                continue;
            }
            nodes_to_mesh_entity_ids[node_id] = nb_mesh_entities++;
        }
        RINGMesh::GeoModelBuilder<DIMENSION> builder(output_geomodel);
        builder.topology.create_mesh_entities(type, nb_mesh_entities);
    }

    template <index_t DIMENSION>
    void build_line_boundaries(
        const GeoModelTopologyRecoverer<DIMENSION> &topology_recoverer,
        const std::vector<index_t> &nodes_to_corner_ids,
        const std::vector<index_t> &nodes_to_line_ids,
        RINGMesh::GeoModel<DIMENSION> &output_geomodel)
    {
        RINGMesh::GeoModelBuilder<DIMENSION> builder(output_geomodel);
        for (const auto node_id : RINGMesh::range(nodes_to_line_ids.size()))
        {
            index_t line_id{nodes_to_line_ids[node_id]};
            if (line_id == NO_ID)
            {
                continue;
            }
            const auto &boundary_nodes =
                topology_recoverer.node_information[node_id].ordered_boundary_nodes;
            scar_assert(boundary_nodes.size() == 2);
            for (auto boundary_node : boundary_nodes)
            {
                builder.topology.add_line_corner_boundary_relation(line_id,
                                                                   nodes_to_corner_ids[boundary_node]);
            }
            std::cout << "Line " << line_id << " : "
                      << nodes_to_corner_ids[boundary_nodes[0]] << " - "
                      << nodes_to_corner_ids[boundary_nodes[1]] << std::endl;
        }
    }
}

namespace SCAR
{

    template <index_t DIMENSION>
    GeoModelTopologyMaker<DIMENSION>::GeoModelTopologyMaker(
        const GeoModelTopologyRecoverer<DIMENSION> &topology_recoverer,
        RINGMesh::GeoModel<DIMENSION> &output_geomodel,
        bool verbose)
        : topology_recoverer_(topology_recoverer),
          nodes_to_corner_ids_(topology_recoverer.node_information.size(),
                               NO_ID),
          nodes_to_line_ids_(topology_recoverer.node_information.size(), NO_ID),
          output_geomodel_(output_geomodel),
          verbose_(verbose)
    {
    }

    template <index_t DIMENSION>
    void GeoModelTopologyMaker<DIMENSION>::build_geomodel_topology()
    {
        build_mesh_entities(RINGMesh::corner_type_name_static(),
                            topology_recoverer_, nodes_to_corner_ids_, output_geomodel_);
        if (verbose_)
        {
            DEBUG(output_geomodel_.nb_corners());
        }

        build_mesh_entities(RINGMesh::line_type_name_static(), topology_recoverer_,
                            nodes_to_line_ids_, output_geomodel_);
        if (verbose_)
        {
            DEBUG(output_geomodel_.nb_lines());
        }

        // Corner-Line Connectivity
        build_line_boundaries(topology_recoverer_, nodes_to_corner_ids_,
                              nodes_to_line_ids_, output_geomodel_);

        for (auto i : RINGMesh::range(output_geomodel_.nb_corners()))
        {
            std::cout << i << " : " << output_geomodel_.corner(i).nb_incident_entities() << std::endl;
        }
    }

    template class scar_api GeoModelTopologyMaker<2>;

    template class scar_api GeoModelTopologyMaker<3>;
}
