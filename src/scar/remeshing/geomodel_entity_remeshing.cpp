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

#include <scar/remeshing/geomodel_entity_remeshing.h>

#include <deque>

#ifdef SCAR_WITH_MMG
extern "C"
{
#include <mmg/libmmg.h>
}
#endif

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/geomodel/core/entity_type.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/mesh/mesh_aabb.h>
#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/mesh_index.h>

#include <scar/tools/geometry.h>
#include <scar/repair/topology_recovery.h>
#include <scar/repair/geomodel_correspondence_map.h>
#include <scar/tools/distance.h>
#include <scar/tools/utils.h>

namespace SCAR
{

    template <index_t DIMENSION>
    std::vector<double> determine_line_resolutions(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        const std::vector<double> &surface_resolutions)
    {
        std::vector<double> line_resolutions(geomodel.nb_lines());
        for (auto &line : geomodel.lines())
        {
            double min_mesh_size = max_double();
            for (auto incident_surf : RINGMesh::range(line.nb_incident_entities()))
            {
                double cur_mesh_size = surface_resolutions[line.incident_entity(
                                                                   incident_surf)
                                                               .index()];
                if (cur_mesh_size < min_mesh_size)
                {
                    min_mesh_size = cur_mesh_size;
                }
            }
            scar_assert(min_mesh_size != max_double());
            line_resolutions[line.index()] = min_mesh_size;
        }
        return line_resolutions;
    }

    static const int MMG_START_ID = 1;

#ifdef SCAR_WITH_MMG
    void mesh_surfaces_using_mmg2d(
        RINGMesh::GeoModel2D &geomodel,
        const std::vector<double> &surface_resolutions)
    {

        auto line_resolutions = determine_line_resolutions(geomodel,
                                                           surface_resolutions);

        for (index_t s_id : RINGMesh::range(geomodel.nb_surfaces()))
        {

            const auto &surface = geomodel.surface(s_id);
            std::set<index_t> boundary_line_ids;
            for (auto b : RINGMesh::range(surface.nb_boundaries()))
            {
                boundary_line_ids.insert(surface.boundary(b).index());
            }
            std::vector<index_t> boundary_corners;
            std::vector<index_t> corner_vertex_ids;

            std::vector<vec2> vertices;
            std::vector<double> vertex_local_mesh_sizes;
            std::vector<index_t> edges;

            for (auto l : boundary_line_ids)
            {
                const auto &line = geomodel.line(l);
                double cur_line_resolution = line_resolutions[l];

                // First corner
                index_t found_first_corner = RINGMesh::find(boundary_corners,
                                                            geomodel.line(l).boundary(0).index());
                if (found_first_corner != NO_ID)
                {
                    edges.push_back(corner_vertex_ids[found_first_corner]);
                    vertex_local_mesh_sizes[corner_vertex_ids[found_first_corner]] =
                        std::min(
                            vertex_local_mesh_sizes[corner_vertex_ids[found_first_corner]],
                            cur_line_resolution);
                }
                else
                {
                    boundary_corners.push_back(
                        geomodel.line(l).boundary(0).index());
                    corner_vertex_ids.push_back(
                        static_cast<index_t>(vertices.size()));

                    edges.push_back(static_cast<index_t>(vertices.size()));
                    vertices.push_back(line.vertex(0));
                    vertex_local_mesh_sizes.push_back(cur_line_resolution);
                }

                for (auto ei : RINGMesh::range(1, line.nb_mesh_elements()))
                {
                    edges.push_back(static_cast<index_t>(vertices.size()));
                    edges.push_back(static_cast<index_t>(vertices.size()));
                    vertices.push_back(line.mesh_element_vertex({ei, 0}));
                    vertex_local_mesh_sizes.push_back(cur_line_resolution);
                }
                // Last corner
                index_t found_last_corner = RINGMesh::find(boundary_corners,
                                                           geomodel.line(l).boundary(1).index());
                if (found_last_corner != NO_ID)
                {
                    edges.push_back(corner_vertex_ids[found_last_corner]);
                }
                else
                {
                    boundary_corners.push_back(
                        geomodel.line(l).boundary(1).index());
                    corner_vertex_ids.push_back(
                        static_cast<index_t>(vertices.size()));
                    edges.push_back(static_cast<index_t>(vertices.size()));
                    vertices.push_back(line.vertex(line.nb_vertices() - 1));
                    vertex_local_mesh_sizes.push_back(cur_line_resolution);
                }
            }
            bool test_next_mmg_subdomain = true;
            int mmg_subdomain = 1;
            while (test_next_mmg_subdomain)
            {
                MMG5_pMesh mmg_mesh_ptr = nullptr;
                MMG5_pSol mmg_metric_ptr = nullptr;

                MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmg_mesh_ptr,
                                MMG5_ARG_ppMet, &mmg_metric_ptr, MMG5_ARG_end);

                // Keep only first subdomain
                MMG2D_Set_iparameter(mmg_mesh_ptr, mmg_metric_ptr,
                                     MMG2D_IPARAM_numsubdomain, mmg_subdomain);

                // TODO : using boundary_line_edges
                MMG2D_Set_meshSize(mmg_mesh_ptr,
                                   static_cast<int>(vertices.size()), 0,
                                   static_cast<int>(edges.size() / 2));

                /// Set the points
                for (index_t p : RINGMesh::range(vertices.size()))
                {
                    ///@todo Find what is the meaning of ref (set to 0)
                    MMG2D_Set_vertex(mmg_mesh_ptr, vertices[p].x, vertices[p].y, 0,
                                     static_cast<int>(MMG_START_ID + p));
                    MMG2D_Set_requiredVertex(mmg_mesh_ptr,
                                             static_cast<int>(p));
                }

                for (auto e : RINGMesh::range(edges.size() / 2))
                {
                    ///@todo Find what is the meaning of ref (set to 0)
                    MMG2D_Set_edge(mmg_mesh_ptr,
                                   static_cast<int>(MMG_START_ID + edges[2 * e]),
                                   static_cast<int>(MMG_START_ID + edges[2 * e + 1]), 0,
                                   static_cast<int>(MMG_START_ID + e));
                    MMG2D_Set_requiredEdge(mmg_mesh_ptr,
                                           MMG_START_ID + static_cast<int>(e));
                }

                /// Set the metric
                MMG2D_Set_solSize(mmg_mesh_ptr, mmg_metric_ptr, MMG5_Vertex,
                                  static_cast<int>(vertices.size()), MMG5_Scalar);

                for (index_t p = 0; p < vertices.size(); p++)
                {
                    MMG2D_Set_scalarSol(mmg_metric_ptr, vertex_local_mesh_sizes[p],
                                        static_cast<int>(MMG_START_ID + p));
                }
                ///////////////////////////////////////////////
                // Generate triangular mesh
                MMG2D_mmg2dmesh(mmg_mesh_ptr, mmg_metric_ptr);
                ////////////////////////////////////////////

                RINGMesh::GeoModelBuilder2D geomodel_builder(geomodel);
                int nbv, nbt, nbb;
                MMG2D_Get_meshSize(mmg_mesh_ptr, &nbv, &nbt, &nbb);
                std::vector<vec2> surface_vertices(
                    static_cast<index_t>(nbv));
                for (auto &vertex : surface_vertices)
                {
                    double x, y;
                    int ref, is_corner, is_required;
                    MMG2D_Get_vertex(mmg_mesh_ptr, &x, &y, &ref, &is_corner,
                                     &is_required);
                    vertex = vec2(x, y);
                }
                std::vector<index_t> surface_polygons(
                    static_cast<index_t>(3 * nbt));
                std::vector<index_t> surface_polygon_ptr;
                surface_polygon_ptr.reserve(static_cast<index_t>(nbt + 1));
                for (auto triangle : RINGMesh::range(nbt))
                {
                    int v0, v1, v2, ref, is_required;
                    MMG2D_Get_triangle(mmg_mesh_ptr, &v0, &v1, &v2, &ref,
                                       &is_required);
                    surface_polygon_ptr.push_back(3 * triangle);
                    surface_polygons[3 * triangle] = static_cast<index_t>(v0 - MMG_START_ID);
                    surface_polygons[3 * triangle + 1] = static_cast<index_t>(v1 - MMG_START_ID);
                    surface_polygons[3 * triangle + 2] = static_cast<index_t>(v2 - MMG_START_ID);
                }

                // Test if all the vertices of the boundaries
                // are in the triangular mesh
                std::vector<bool> vertex_is_in_mesh(surface_vertices.size(),
                                                    false);
                for (auto v_id : surface_polygons)
                {
                    vertex_is_in_mesh[v_id] = true;
                }
                test_next_mmg_subdomain = (std::find(vertex_is_in_mesh.begin(),
                                                     vertex_is_in_mesh.end(), false) != vertex_is_in_mesh.end());

                if (test_next_mmg_subdomain)
                {
                    mmg_subdomain++;
                }
                else
                {

                    surface_polygon_ptr.push_back(
                        static_cast<index_t>(surface_polygons.size()));
                    geomodel_builder.geometry.set_surface_geometry(s_id,
                                                                   surface_vertices, surface_polygons, surface_polygon_ptr);
                }
            }
        }
    }

#endif

    void set_surface_new_line_boundary_relation(
        RINGMesh::GeoModelBuilder2D &builder,
        const RINGMesh::Line2D &old_line,
        index_t new_line_id)
    {
        for (auto i : RINGMesh::range(old_line.nb_incident_entities()))
        {
            builder.topology.add_surface_line_boundary_relation(
                old_line.incident_entity(i).index(), new_line_id, true);
        }
    }

    void set_surface_new_line_boundary_relation(
        RINGMesh::GeoModelBuilder3D &builder,
        const RINGMesh::Line3D &old_line,
        index_t new_line_id)
    {
        for (auto i : RINGMesh::range(old_line.nb_incident_entities()))
        {
            builder.topology.add_surface_line_boundary_relation(
                old_line.incident_entity(i).index(), new_line_id);
        }
    }

    //@author Thomas Saglio (2017)
    template <index_t DIMENSION>
    void merge_lines(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        const RINGMesh::Corner<DIMENSION> &corner)
    {
        // Prepare lines
        std::vector<RINGMesh::gmme_id> lines;
        for (auto e : RINGMesh::range(corner.nb_incident_entities()))
        {
            if (corner.incident_entity(e).type_name() == RINGMesh::line_type_name_static())
            {
                lines.push_back(corner.incident_entity_gmme(e));
            }
        }
        scar_assert(lines.size() == 2);

        const RINGMesh::Line<DIMENSION> &l1 = geomodel.line(lines[0].index());
        const RINGMesh::Line<DIMENSION> &l2 = geomodel.line(lines[1].index());

        // Prepare the list of vertices corresponding to the result of
        // merging l1 and l2. The vertices are ordered appropriately.
        std::vector<RINGMesh::vecn<DIMENSION>> vertices;

        if (l1.vertex(l1.nb_vertices() - 1) == corner.vertex(0))
        {
            for (RINGMesh::index_t v{0}; v < l1.nb_vertices() - 1; ++v) // Discard the last vertex (common corner).
                vertices.push_back(l1.vertex(v));
        }
        else
        {
            for (RINGMesh::index_t v{l1.nb_vertices() - 1}; v > 0; --v) // Discard the first vertex (common corner), reverse iteration.
                vertices.push_back(l1.vertex(v));
        }
        vertices.push_back(corner.vertex(0));
        if (l2.vertex(0) == corner.vertex(0))
        {
            for (RINGMesh::index_t v{1}; v < l2.nb_vertices(); ++v) // Discard the first vertex (common corner).
                vertices.push_back(l2.vertex(v));
        }
        else
        {
            for (RINGMesh::index_t v{l2.nb_vertices() - 1}; v > 0; --v) // Discard the last vertex (common corner), reverse iteration.
                vertices.push_back(l2.vertex(v - 1));
        }

        // Make changes to the geomodel.
        RINGMesh::GeoModelBuilder<DIMENSION> builder(geomodel);

        // Create the new line which is the result of merging l1 and l2.
        RINGMesh::gmme_id merged_line = builder.topology.create_mesh_entity(
            RINGMesh::line_type_name_static());
        builder.geometry.set_line(merged_line.index(), vertices);

        // Bind corners to merged_line.
        RINGMesh::gmme_id c1{builder.topology.find_or_create_corner(
            vertices.front())};
        RINGMesh::gmme_id c2{builder.topology.find_or_create_corner(
            vertices.back())};
        builder.topology.add_line_corner_boundary_relation(merged_line.index(),
                                                           c1.index());
        builder.topology.add_line_corner_boundary_relation(merged_line.index(),
                                                           c2.index());
        scar_assert(l1.nb_incident_entities() == l2.nb_incident_entities());
        set_surface_new_line_boundary_relation(builder, l1, merged_line.index());

        builder.remove.remove_mesh_entities(std::set<RINGMesh::gmme_id>{
            l1.gmme(), l2.gmme()});
    }

    template <index_t DIMENSION>
    void remove_corners_with_valence_2(RINGMesh::GeoModel<DIMENSION> &geomodel)
    {
        RINGMesh::GeoModelBuilder<DIMENSION> builder(geomodel);

        bool restart = true;
        while (restart)
        {
            restart = false;
            for (const auto &corner : geomodel.corners())
            {
                if (corner.nb_incident_entities() == 2)
                {
                    // Conserve corners ( angle < 120 degrees)
                    double angle = compute_angle_between_two_adjacent_lines(corner,
                                                                            corner.incident_entity(0), corner.incident_entity(1));
                    if (angle < 120)
                    {
                        continue;
                    }
                    merge_lines(geomodel, corner);
                    std::set<RINGMesh::gmme_id> to_delete{corner.gmme()};
                    builder.remove.remove_mesh_entities(to_delete);
                    restart = true;
                    break;
                }
                if (corner.nb_incident_entities() == 0)
                {
                    std::set<RINGMesh::gmme_id> to_delete{corner.gmme()};
                    builder.remove.remove_mesh_entities(to_delete);
                    restart = true;
                    break;
                }
            }
        }
    }

    template void scar_api remove_corners_with_valence_2(RINGMesh::GeoModel2D &);

    template void scar_api remove_corners_with_valence_2(RINGMesh::GeoModel3D &);
}

namespace
{

    using namespace SCAR;

    template <index_t DIMENSION>
    vecn<DIMENSION> compute_new_corner_position(
        const RINGMesh::GeoModel<DIMENSION> &input_geomodel,
        const std::vector<index_t> &mapped_input_corners)
    {
        if (mapped_input_corners.empty())
        {
            Logger::warn("Corner remeshing",
                         "Input list of corner for remeshing is empty");
            return vecn<DIMENSION>();
        }
        vecn<DIMENSION> new_position;
        for (index_t input_corner_id : mapped_input_corners)
        {
            new_position += input_geomodel.corner(input_corner_id).vertex(0);
        }

        return (new_position / mapped_input_corners.size());
    }

    template <index_t DIMENSION>
    vecn<DIMENSION> point_set_barycenter(
        const std::vector<vecn<DIMENSION>> &points)
    {
        vecn<DIMENSION> result;
        for (auto &point : points)
        {
            result += point;
        }
        return (result / points.size());
    }

    template <index_t DIMENSION>
    void remesh_corners(
        RINGMesh::GeoModel<DIMENSION> &remeshed_geomodel,
        const RINGMesh::GeoModel<DIMENSION> &input_geomodel,
        const GeoModelTopologyMaker<DIMENSION> &geomodel_topology_maker)
    {
        RINGMesh::GeoModelBuilder<DIMENSION> builder(remeshed_geomodel);
        for (auto node_id : RINGMesh::range(
                 geomodel_topology_maker.nodes_to_corner_ids_.size()))
        {
            const auto corner_id =
                geomodel_topology_maker.nodes_to_corner_ids_[node_id];
            if (corner_id == NO_ID)
            {
                continue;
            }
            const auto &node =
                geomodel_topology_maker.topology_recoverer_.node_information[node_id];

            std::vector<index_t> input_corners;
            std::vector<EntityPart> input_line_parts;
            for (const auto &part : node.parts)
            {
                if (part.id.type() == RINGMesh::corner_type_name_static())
                {
                    input_corners.push_back(part.id.index());
                }
                else if (part.id.type() == RINGMesh::line_type_name_static())
                {
                    input_line_parts.push_back(part);
                }
                else
                {
                    scar_assert_not_reached;
                }
            }
            bool determined = false;

            if (input_corners.empty() && !input_line_parts.empty())
            {
                std::vector<vecn<DIMENSION>> input_positions;
                for (const auto &line_part : input_line_parts)
                {
                    if (enum_contains(
                            geomodel_topology_maker.topology_recoverer_.geological_rules,
                            GeologicalRule::KEEP_BOUNDARIES_UNMOVED))
                    {
                        if (!input_geomodel.line(line_part.id.index()).is_on_voi())
                        {
                            continue;
                        }
                    }
                    scar_assert(line_part.part.begin == line_part.part.end);
                    index_t v_id = static_cast<index_t>(line_part.part.begin);
                    scar_assert(
                        (double)v_id - line_part.part.begin < global_epsilon);
                    if (static_cast<double>(v_id) == line_part.part.begin)
                    {
                        input_positions.push_back(
                            input_geomodel.line(line_part.id.index()).vertex(v_id));
                    }
                    else
                    {
                        vecn<DIMENSION> corner_position = (1 - line_part.part.begin.norm) * input_geomodel.line(line_part.id.index()).vertex(v_id) + (line_part.part.begin.norm) * input_geomodel.line(line_part.id.index()).vertex(v_id + 1);
                        input_positions.push_back(corner_position);
                    }
                }
                if (!input_positions.empty())
                {
                    auto barycenter = point_set_barycenter(input_positions);
                    builder.geometry.set_corner(corner_id, barycenter);
                    determined = true;
                }
                else
                {
                    for (const auto &line_part : input_line_parts)
                    {
                        scar_assert(line_part.part.begin == line_part.part.end);
                        index_t v_id = static_cast<index_t>(line_part.part.begin);
                        scar_assert(
                            (double)v_id - line_part.part.begin < global_epsilon);
                        input_positions.push_back(
                            input_geomodel.line(line_part.id.index()).vertex(v_id));
                    }
                    auto barycenter = point_set_barycenter(input_positions);
                    builder.geometry.set_corner(corner_id, barycenter);
                    determined = true;
                }
            }
            if (determined)
            {
                continue;
            }
            if (enum_contains(
                    geomodel_topology_maker.topology_recoverer_.geological_rules,
                    GeologicalRule::KEEP_BOUNDARIES_UNMOVED))
            {
                std::vector<bool> input_corners_on_boundaries(
                    input_corners.size(), false);
                for (auto c : RINGMesh::range(input_corners.size()))
                {
                    const auto &corner = input_geomodel.corner(input_corners[c]);
                    for (auto il : RINGMesh::range(corner.nb_incident_entities()))
                    {
                        const auto &incident_line = corner.incident_entity(il);
                        if (incident_line.is_on_voi())
                        {
                            input_corners_on_boundaries[c] = true;
                            break;
                        }
                    }
                }
                index_t nb = static_cast<index_t>(std::count(
                    input_corners_on_boundaries.begin(),
                    input_corners_on_boundaries.end(), true));
                if (nb == 1)
                {
                    index_t voi_corner =
                        static_cast<index_t>(input_corners_on_boundaries.begin() - std::find(input_corners_on_boundaries.begin(),
                                                                                             input_corners_on_boundaries.end(), true));
                    builder.geometry.set_corner(corner_id,
                                                input_geomodel.corner(input_corners[voi_corner]).vertex(0));
                    determined = true;
                }
                else if (nb > 1)
                {
                    vecn<DIMENSION> new_location;
                    for (auto c : RINGMesh::range(input_corners.size()))
                    {
                        if (input_corners_on_boundaries[c])
                        {
                            new_location +=
                                input_geomodel.corner(input_corners[c]).vertex(0);
                        }
                    }
                    new_location /= nb;
                    builder.geometry.set_corner(corner_id, new_location);
                    determined = true;
                }
            }
            if (determined)
            {
                continue;
            }

            // DFN
            vecn<DIMENSION> barycenter;
            if (enum_contains(
                    geomodel_topology_maker.topology_recoverer_.geological_rules,
                    GeologicalRule::REMESH_DFN))
            {

                // Compute powers
                std::vector<index_t> input_corner_powers(input_corners.size(),
                                                         0);
                for (auto c : RINGMesh::range(input_corners.size()))
                {
                    input_corner_powers[c] = input_geomodel.corner(
                                                               input_corners[c])
                                                 .nb_incident_entities();
                }
                for (const auto &line : input_geomodel.lines())
                {
                    auto it_b0 = RINGMesh::find(input_corners,
                                                line.boundary(0).index());
                    if (it_b0 == NO_ID)
                    {
                        continue;
                    }
                    auto it_b1 = RINGMesh::find(input_corners,
                                                line.boundary(1).index());
                    if (it_b1 == NO_ID)
                    {
                        continue;
                    }
                    // the line has its two corners to merge
                    // so do not count this line in power
                    --input_corner_powers[it_b0];
                    --input_corner_powers[it_b1];
                }

                // Set barycenter
                index_t max_power = 0;
                index_t nb_max_powered_corners = 0;
                for (auto i : RINGMesh::range(input_corners.size()))
                {
                    if (input_corner_powers[i] < max_power)
                    {
                        continue;
                    }
                    if (input_corner_powers[i] == max_power)
                    {
                        nb_max_powered_corners++;
                        barycenter +=
                            input_geomodel.corner(input_corners[i]).vertex(0);
                        continue;
                    }
                    // New max -> reset
                    max_power = input_corner_powers[i];
                    nb_max_powered_corners = 1;
                    barycenter = input_geomodel.corner(input_corners[i]).vertex(0);
                }
                scar_assert(max_power > 0);
                scar_assert(nb_max_powered_corners > 0);
                barycenter /= nb_max_powered_corners;
            }
            else
            {
                barycenter = compute_new_corner_position(input_geomodel,
                                                         input_corners);
            }
            if (input_line_parts.empty())
            {
                // Only input corners
                builder.geometry.set_corner(corner_id, barycenter);
            }
            else
            {
                // Corners to project on lines
                std::vector<vecn<DIMENSION>> barycenter_projections;
                barycenter_projections.reserve(input_line_parts.size());
                for (const auto &line_part : input_line_parts)
                {
                    const RINGMesh::Line<DIMENSION> &input_line =
                        input_geomodel.line(line_part.id.index());
                    vecn<DIMENSION> barycenter_projection;
                    index_t edge_id;
                    std::tie(edge_id, barycenter_projection, std::ignore) =
                        input_line.edge_aabb().closest_edge(barycenter);
                    RINGMesh::Geometry::Segment<DIMENSION> edge(
                        input_line.vertex(edge_id),
                        input_line.vertex(edge_id + 1));
                    LineRelativeCoordinates projection_coord(edge_id,
                                                             segment_point_normalized_coordinates(edge,
                                                                                                  barycenter_projection));
                    if (line_part.part.is_inside(projection_coord))
                    {
                        barycenter_projections.push_back(barycenter_projection);
                    }
                    else if (line_part.part.is_before(projection_coord))
                    {
                        vecn<DIMENSION> line_chunk_begin = PointOnLineCoordinates<
                            DIMENSION>::absolute_coordinates(line_part.part.begin,
                                                             input_line.mesh());
                        barycenter_projections.push_back(line_chunk_begin);
                    }
                    else if (line_part.part.is_after(projection_coord))
                    {
                        vecn<DIMENSION> line_chunk_end = PointOnLineCoordinates<
                            DIMENSION>::absolute_coordinates(line_part.part.end,
                                                             input_line.mesh());
                        barycenter_projections.push_back(line_chunk_end);
                    }
                    else
                    {
                        scar_assert_not_reached;
                    }
                }
                auto new_corner_location = point_set_barycenter(
                    barycenter_projections);
                builder.geometry.set_corner(corner_id, new_corner_location);
            }
        }
    }

    enum struct LineScanDirection
    {
        NORMAL,
        INVERSE
    };

    template <index_t DIMENSION>
    double compute_coordinates_along_line(
        const RINGMesh::Line<DIMENSION> &line,
        const vecn<DIMENSION> &query_point,
        const LineScanDirection scan_direction)
    {
        /// Point must be on the line
        index_t closest_edge;
        vecn<DIMENSION> nearest_point;
        double distance;
        std::tie(closest_edge, nearest_point, distance) =
            line.edge_aabb().closest_edge(query_point);
        scar_assert(distance < line.geomodel().epsilon());

        // Get total curvilinear distance between
        // query point and line first boundary
        double curv_distance_until_query_point = 0.0;
        for (index_t e = 0; e < closest_edge; ++e)
        {
            curv_distance_until_query_point += line.mesh_element_size(e);
        }
        double distance_on_closest_edge =
            (query_point - line.mesh().vertex(line.mesh().edge_vertex({closest_edge, 0}))).length();
        curv_distance_until_query_point += distance_on_closest_edge;
        scar_assert(
            curv_distance_until_query_point / line.size() <= 1 && curv_distance_until_query_point / line.size() >= 0);
        if (scan_direction == LineScanDirection::NORMAL)
        {
            return curv_distance_until_query_point / line.size();
        }
        else
        {
            return 1 - (curv_distance_until_query_point / line.size());
        }
    }

    template <index_t DIMENSION>
    void perform_line_screening(
        const RINGMesh::GeoModel<DIMENSION> &input_geomodel,
        vecn<DIMENSION> &cur_point,
        const vecn<DIMENSION> &last_point,
        const std::vector<EntityPart> &input_line_parts,
        const std::vector<LineScanDirection> &input_line_scan_direction,
        std::vector<index_t> &input_line_scan_next_vertex,
        vecn<DIMENSION> &scan_direction,
        std::vector<vecn<DIMENSION>> &line_vertices);

    template <>
    void perform_line_screening<2>(
        const RINGMesh::GeoModel2D &input_geomodel,
        vec2 &cur_point,
        const vec2 &last_point,
        const std::vector<EntityPart> &input_line_parts,
        const std::vector<LineScanDirection> &input_line_scan_direction,
        std::vector<index_t> &input_line_scan_next_vertex,
        vec2 &scan_direction,
        std::vector<vec2> &line_vertices)
    {
        // scan_direction is the normal to the screening line
        while (true)
        {
            vec2 next_point = NO_POINT_2D;
            double smallest_dot = RINGMesh::max_float64();
            vec2 next_scan_direction = NO_POINT_2D;
            index_t progressing_line_part_id = NO_ID;

            //Find next_point
            for (auto line_part_id : RINGMesh::range(input_line_parts.size()))
            {
                while (true)
                {
                    const auto &line_part = input_line_parts[line_part_id];
                    const auto &next_cur_line_vertex_rel_coord =
                        LineRelativeCoordinates(
                            input_line_scan_next_vertex[line_part_id], 0);
                    if (!line_part.part.is_inside(
                            next_cur_line_vertex_rel_coord))
                    {
                        // continue;
                        break;
                    }

                    const auto &next_cur_line_vertex_abs_coord = input_geomodel.line(
                                                                                   line_part.id.index())
                                                                     .vertex(
                                                                         input_line_scan_next_vertex[line_part_id]);
                    // Compute dot product and compare it
                    auto cur_next_direction = next_cur_line_vertex_abs_coord - cur_point;
                    auto cur_dot_product = dot(scan_direction, cur_next_direction);
                    if (cur_dot_product <= 0)
                    {
                        // Do another step in the while
                        if (input_line_scan_direction[line_part_id] == LineScanDirection::NORMAL)
                        {
                            input_line_scan_next_vertex[line_part_id]++;
                        }
                        else
                        {
                            input_line_scan_next_vertex[line_part_id]--;
                        }
                        continue;
                    }
                    if (cur_dot_product < smallest_dot)
                    {
                        smallest_dot = cur_dot_product;
                        next_point = next_cur_line_vertex_abs_coord;
                        next_scan_direction = cur_next_direction;
                        progressing_line_part_id = line_part_id;
                    }
                    break;
                }
            }

            // End of the loop
            if (next_point == NO_POINT_2D)
            {
                line_vertices.push_back(last_point);
                return;
            }

            // Compare also with last_point
            auto last_point_dot_product = dot(scan_direction,
                                              last_point - cur_point);
            scar_assert(last_point_dot_product != 0.);
            double relative_diff = (last_point_dot_product - smallest_dot) / last_point_dot_product;
            if (last_point_dot_product >= 0 && relative_diff < 0.01)
            {
                line_vertices.push_back(last_point);
                return;
            }
            scar_assert(next_point != NO_POINT_2D);

            // Compute barycenter between line part intersections and scan screen
            if (input_line_parts.size() > 1)
            {
                RINGMesh::Geometry::Line2D scan_screen(
                    RINGMesh::normalized_perp(scan_direction), next_point);
                std::vector<vec2> scan_screen_intersections;
                scan_screen_intersections.reserve(input_line_parts.size());
                scan_screen_intersections.push_back(next_point);
                for (auto line_part_id : RINGMesh::range(input_line_parts.size()))
                {
                    if (line_part_id == progressing_line_part_id)
                    {
                        continue;
                    }
                    const auto &line_part = input_line_parts[line_part_id];
                    const auto &cur_line = input_geomodel.line(
                        line_part.id.index());
                    index_t edge_id = 0;
                    if (input_line_scan_next_vertex[line_part_id] != 0)
                    {
                        edge_id = input_line_scan_next_vertex[line_part_id] - 1;
                    }
                    RINGMesh::Geometry::Segment2D line_edge(
                        cur_line.mesh_element_vertex(RINGMesh::ElementLocalVertex{
                            edge_id, 0}),
                        cur_line.mesh_element_vertex(RINGMesh::ElementLocalVertex{
                            edge_id, 1}));
                    vec2 cur_screen_intersection;
                    bool is_intersecting;
                    std::tie(is_intersecting, cur_screen_intersection) =
                        RINGMesh::Intersection::segment_line(line_edge,
                                                             scan_screen);
                    if (is_intersecting)
                    {
                        scan_screen_intersections.push_back(
                            cur_screen_intersection);
                    }
                    else
                    {
                        // Numerical precision
                        RINGMesh::Geometry::Line2D line_edge_infinite(
                            normalize(line_edge.p1 - line_edge.p0), line_edge.p1);
                        std::tie(is_intersecting, cur_screen_intersection) =
                            RINGMesh::Intersection::line_line(line_edge_infinite,
                                                              scan_screen);
                        double distance;
                        std::tie(distance, std::ignore) =
                            RINGMesh::Distance::point_to_segment(
                                cur_screen_intersection, line_edge);
                        if (distance < RINGMesh::global_epsilon)
                        {
                            scan_screen_intersections.push_back(
                                cur_screen_intersection);
                        }
                    }
                }
                scar_assert(!scan_screen_intersections.empty());
                vec2 barycenter = point_set_barycenter(scan_screen_intersections);
                line_vertices.push_back(barycenter);
                scan_direction = barycenter - cur_point;
                cur_point = barycenter;
            }
            else
            {
                line_vertices.push_back(next_point);
                scan_direction = next_scan_direction;
                cur_point = next_point;
            }

            // For all the line part, if next line vertex is very close to the
            // screening line, we make progressing this line vertex.
            for (auto line_part_id : RINGMesh::range(input_line_parts.size()))
            {
                const auto &line_part = input_line_parts[line_part_id];
                const auto &next_cur_line_vertex_abs_coord = input_geomodel.line(
                                                                               line_part.id.index())
                                                                 .vertex(
                                                                     input_line_scan_next_vertex[line_part_id]);
                auto cur_next_direction = next_cur_line_vertex_abs_coord - cur_point;
                auto cur_dot_product = dot(scan_direction, cur_next_direction);
                if (std::fabs(cur_dot_product) < 10 - 4)
                {
                    if (input_line_scan_direction[line_part_id] == LineScanDirection::NORMAL)
                    {
                        input_line_scan_next_vertex[line_part_id]++;
                    }
                    else
                    {
                        input_line_scan_next_vertex[line_part_id]--;
                    }
                }
            }
        }
    }

    template <>
    void perform_line_screening<3>(
        const RINGMesh::GeoModel3D &input_geomodel,
        vec3 &cur_point,
        const vec3 &last_point,
        const std::vector<EntityPart> &input_line_parts,
        const std::vector<LineScanDirection> &input_line_scan_direction,
        std::vector<index_t> &input_line_scan_next_vertex,
        vec3 &scan_direction,
        std::vector<vec3> &line_vertices)
    {
        scar_unused(input_geomodel);
        scar_unused(cur_point);
        scar_unused(last_point);
        scar_unused(input_line_parts);
        scar_unused(input_line_scan_direction);
        scar_unused(input_line_scan_next_vertex);
        scar_unused(scan_direction);
        scar_unused(line_vertices);
        throw SCARException("3D",
                            "Perform line screening in 3D is not yet implemented");
    }

    template <index_t DIMENSION>
    std::vector<vecn<DIMENSION>> compute_new_line_vertices(
        const RINGMesh::Line<DIMENSION> &remeshed_line,
        const RINGMesh::GeoModel<DIMENSION> &input_geomodel,
        const std::vector<EntityPart> &input_line_parts)
    {
        // Initialize line vertices
        std::vector<vecn<DIMENSION>> line_vertices;
        // Add first boundary line vertex
        scar_assert(remeshed_line.nb_boundaries() == 2);
        vecn<DIMENSION> cur_point = remeshed_line.boundary(0).vertex(0);
        vecn<DIMENSION> last_point = remeshed_line.boundary(1).vertex(0);
        scar_assert(cur_point != last_point);
        line_vertices.push_back(cur_point);

        // Determine scanning direction for each line part
        // If the input line is scanned in the sense b0 -> b1, NORMAL
        // else (sense b1 -> b0), INVERSE
        std::vector<LineScanDirection> input_line_scan_direction(
            input_line_parts.size(), LineScanDirection::NORMAL);
        std::vector<index_t> input_line_scan_next_vertex(input_line_parts.size(),
                                                         NO_ID);
        vecn<DIMENSION> scan_direction;
        index_t cur_line_part = 0;
        for (auto line_part : input_line_parts)
        {
            scar_assert(line_part.part.begin < line_part.part.end);
            vecn<DIMENSION> cur_point_projection;
            std::tie(std::ignore, cur_point_projection, std::ignore) =
                input_geomodel.line(line_part.id.index()).edge_aabb().closest_edge(cur_point);
            double cur_point_rel_coordinates =
                PointOnLineCoordinates<DIMENSION>::relative_coordinates(
                    cur_point_projection,
                    input_geomodel.line(line_part.id.index()).mesh());
            double distance_to_line_part_begin = std::fabs(
                cur_point_rel_coordinates - line_part.part.begin);
            double distance_to_line_part_end = std::fabs(
                cur_point_rel_coordinates - line_part.part.end);
            if (distance_to_line_part_begin < distance_to_line_part_end)
            {
                input_line_scan_direction[cur_line_part] = LineScanDirection::NORMAL;
                input_line_scan_next_vertex[cur_line_part] =
                    static_cast<index_t>(std::ceil(cur_point_rel_coordinates));
                if (std::ceil(cur_point_rel_coordinates) == cur_point_rel_coordinates)
                {
                    input_line_scan_next_vertex[cur_line_part]++;
                }
                scan_direction +=
                    (input_geomodel.line(line_part.id.index()).vertex(input_line_scan_next_vertex[cur_line_part]) - cur_point);
            }
            else
            {
                input_line_scan_direction[cur_line_part] =
                    LineScanDirection::INVERSE;
                input_line_scan_next_vertex[cur_line_part] =
                    std::max(static_cast<index_t>(0),
                             static_cast<index_t>(std::floor(
                                 cur_point_rel_coordinates)));
                if (input_geomodel.line(line_part.id.index()).vertex(input_line_scan_next_vertex[cur_line_part]) == cur_point)
                {
                    input_line_scan_next_vertex[cur_line_part]--;
                }
                scan_direction +=
                    (input_geomodel.line(line_part.id.index()).vertex(input_line_scan_next_vertex[cur_line_part]) - cur_point);
            }
            ++cur_line_part;
        }
        scar_assert(cur_line_part > 0);
        scan_direction /= static_cast<double>(cur_line_part);

        perform_line_screening<DIMENSION>(input_geomodel, cur_point, last_point,
                                          input_line_parts, input_line_scan_direction, input_line_scan_next_vertex,
                                          scan_direction, line_vertices);

        return line_vertices;
    }

    template <index_t DIMENSION>
    void remesh_lines(
        RINGMesh::GeoModel<DIMENSION> &remeshed_geomodel,
        const RINGMesh::GeoModel<DIMENSION> &initial_geomodel,
        const GeoModelTopologyMaker<DIMENSION> &topology_maker)
    {
        RINGMesh::GeoModelBuilder<DIMENSION> builder(remeshed_geomodel);
        for (auto node_id : RINGMesh::range(
                 topology_maker.nodes_to_line_ids_.size()))
        {
            const auto line_id = topology_maker.nodes_to_line_ids_[node_id];
            if (line_id == NO_ID)
            {
                continue;
            }
            const auto &node =
                topology_maker.topology_recoverer_.node_information[node_id];

            std::vector<EntityPart> input_line_parts;
            for (const auto &part : node.parts)
            {
                if (part.id.type() == RINGMesh::line_type_name_static())
                {
                    input_line_parts.push_back(part);
                }
                else
                {
                    scar_assert_not_reached;
                }
            }

            if (enum_contains(topology_maker.topology_recoverer_.geological_rules,
                              GeologicalRule::KEEP_BOUNDARIES_UNMOVED))
            {
                std::vector<EntityPart> updated_input_line_parts;
                for (auto line_part : input_line_parts)
                {
                    if (initial_geomodel.line(line_part.id.index()).is_on_voi())
                    {
                        updated_input_line_parts.push_back(line_part);
                    }
                }
                if (!updated_input_line_parts.empty())
                {
                    input_line_parts = updated_input_line_parts;
                }
            }

            std::vector<vecn<DIMENSION>> new_line_vertices =
                compute_new_line_vertices(remeshed_geomodel.line(line_id),
                                          initial_geomodel, input_line_parts);

            // Check angle not so right
            if (enum_contains(topology_maker.topology_recoverer_.geological_rules,
                              GeologicalRule::REMESH_DFN))
            {
                std::vector<vecn<DIMENSION>> new_line_vertices_filtered;
                new_line_vertices_filtered.reserve(new_line_vertices.size());
                for (auto i : RINGMesh::range(0, new_line_vertices.size()))
                {
                    if (i == 0 || i == new_line_vertices.size() - 1)
                    {
                        new_line_vertices_filtered.push_back(new_line_vertices[i]);
                        continue;
                    }
                    auto cur_v = new_line_vertices[i];
                    auto prev_v = new_line_vertices[i - 1];
                    auto next_v = new_line_vertices[i + 1];
                    auto dot_prod = dot(normalize(prev_v - cur_v),
                                        normalize(cur_v - next_v));
                    double threshold = 0.258; // Correspond to 75 degrees
                    if (dot_prod > threshold)
                    {
                        new_line_vertices_filtered.push_back(new_line_vertices[i]);
                    }
                }
                new_line_vertices = new_line_vertices_filtered;
            }

            builder.geometry.set_line(line_id, new_line_vertices);
        }
    }
}

namespace SCAR
{

    template <index_t DIMENSION>
    GeoModelRemesherBase<DIMENSION>::GeoModelRemesherBase(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        const RINGMesh::GeoModel<DIMENSION> &init_geomodel,
        const GeoModelTopologyMaker<DIMENSION> &topology_maker,
        bool verbose)
        : geomodel_(geomodel),
          init_geomodel_(init_geomodel),
          topology_maker_(topology_maker),
          verbose_(verbose)
    {
    }

    template <index_t DIMENSION>
    GeoModelRemesher<DIMENSION>::GeoModelRemesher(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        const RINGMesh::GeoModel<DIMENSION> &init_geomodel,
        const GeoModelTopologyMaker<DIMENSION> &topology_maker,
        bool verbose)
        : GeoModelRemesherBase<DIMENSION>(geomodel, init_geomodel,
                                          topology_maker, verbose)
    {
    }

    GeoModelRemesher<2>::GeoModelRemesher(
        RINGMesh::GeoModel2D &geomodel,
        const RINGMesh::GeoModel2D &init_geomodel,
        const GeoModelTopologyMaker2D &topology_maker,
        bool verbose)
        : GeoModelRemesherBase<2>(geomodel, init_geomodel, topology_maker,
                                  verbose)
    {
    }

    GeoModelRemesher<3>::GeoModelRemesher(
        RINGMesh::GeoModel3D &geomodel,
        const RINGMesh::GeoModel3D &init_geomodel,
        const GeoModelTopologyMaker3D &topology_maker,
        bool verbose)
        : GeoModelRemesherBase<3>(geomodel, init_geomodel, topology_maker,
                                  verbose)
    {
    }

    template <index_t DIMENSION>
    void GeoModelRemesherBase<DIMENSION>::remesh_geomodel_entities()
    {
        remesh_corners(geomodel_, init_geomodel_, topology_maker_);
        remesh_lines(geomodel_, init_geomodel_, topology_maker_);
    }

    void GeoModelRemesher<2>::remesh_geomodel_entities()
    {
        GeoModelRemesherBase<2>::remesh_geomodel_entities();
    }

    void GeoModelRemesher<3>::remesh_geomodel_entities()
    {
        GeoModelRemesherBase<3>::remesh_geomodel_entities();
    }

    template class scar_api GeoModelRemesher<2>;
    template class scar_api GeoModelRemesher<3>;
}
