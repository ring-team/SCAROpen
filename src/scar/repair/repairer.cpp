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

#include <scar/repair/repairer.h>

#include <geogram/basic/stopwatch.h>
#include <geogram/basic/attributes.h>

#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

#include <scar/repair/geomodel_correspondence_map.h>
#include <scar/repair/topology_recovery.h>
#include <scar/remeshing/geomodel_entity_remeshing.h>
#include <scar/remeshing/geomodel_vertex_shifter.h>
#include <scar/tools/distance.h>
#include <scar/tools/entity_analysis.h>

namespace SCAR
{

    template <index_t DIMENSION>
    void assign_constant_tolerance_to_input_geomodel(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        const double &tolerance)
    {
        for (index_t t = 0;
             t < geomodel.entity_type_manager().mesh_entity_manager.nb_mesh_entity_types();
             t++)
        {
            const RINGMesh::MeshEntityType &type =
                geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types()[t];
            for (index_t e = 0; e < geomodel.nb_mesh_entities(type); e++)
            {
                const RINGMesh::GeoModelMeshEntity<DIMENSION> &E =
                    geomodel.mesh_entity(type, e);
                GEO::Attribute<double> tolerance_att(E.vertex_attribute_manager(),
                                                     tolerance_att_name);
                tolerance_att.fill(tolerance);
            }
        }
    }

    template <index_t DIMENSION>
    void assign_exclusion_convex_to_corner(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        const index_t corner_id,
        const ConvexShapeType &tolerance_shape_type,
        const double tolerance)
    {
        const auto &corner = geomodel.corner(corner_id);
        if (tolerance_shape_type == ConvexShapeType::NSphere)
        {
            GEO::Attribute<double> tolerance_att(
                corner.vertex_attribute_manager(), tolerance_att_name);
            tolerance_att.fill(tolerance / 2.);
        }
        else
        {
            Logger::out("Convex", "Not implemented yet");
            scar_assert_not_reached;
        }
    }

    // To improve with incidence relationships
    template <index_t DIMENSION>
    void assign_exclusion_convex_to_line(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        const index_t line_id,
        const ConvexShapeType &tolerance_shape_type,
        const double tolerance,
        const double angle)
    {
        const auto &line = geomodel.line(line_id);
        if (tolerance_shape_type == ConvexShapeType::NSphere)
        {
            GEO::Attribute<double> tolerance_att(line.vertex_attribute_manager(),
                                                 tolerance_att_name);
            for (auto v_id : RINGMesh::range(line.nb_vertices()))
            {
                // Compute angle terms
                double curv_dist0 = compute_curvilinear_distance(line, v_id, 0);
                double curv_dist1 = compute_curvilinear_distance(line, v_id, 1);
                double diameter_angle0 = std::tan(angle / 2.) * curv_dist0;
                double diameter_angle1 = std::tan(angle / 2.) * curv_dist1;
                // Compute size term from
                double max_size = tolerance / 2.;
                tolerance_att[v_id] = std::min(
                    std::min(diameter_angle0, diameter_angle1), max_size);
            }
        }
        else
        {
            Logger::out("Convex", "Not implemented yet..");
            scar_assert_not_reached;
        }
    }

    // To improve with incidence relationships
    void assign_exclusion_convex_to_surface(
        RINGMesh::GeoModel3D &geomodel,
        const index_t surface_id,
        const ConvexShapeType &tolerance_shape_type,
        const double tolerance,
        const double angle)
    {
        const auto &surface = geomodel.surface(surface_id);
        if (tolerance_shape_type == ConvexShapeType::NSphere)
        {
            compute_geodesic_to_borders_djikstra_approximation(geomodel,
                                                               surface_id);

            GEO::Attribute<double> tolerance_att(
                surface.vertex_attribute_manager(), tolerance_att_name);
            GEO::Attribute<double> geodesic_att(
                surface.vertex_attribute_manager(), geodesic_att_name);
            for (auto v_id : RINGMesh::range(surface.nb_vertices()))
            {
                // Compute angle terms
                double geodesic_dist = geodesic_att[v_id];
                double diameter_angle = std::tan(angle / 2.) * geodesic_dist;
                // Compute size term from
                double max_size = tolerance / 2.;
                tolerance_att[v_id] = std::min(diameter_angle, max_size);
            }
        }
        else
        {
            Logger::out("Convex", "Not implemented yet..");
            scar_assert_not_reached;
        }
    }

    template <index_t DIMENSION>
    void assign_exclusion_convex_to_entities_base(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        const ConvexShapeType &tolerance_shape_type,
        const double tolerance,
        const double angle)
    {
        for (auto corner_id : RINGMesh::range(geomodel.nb_corners()))
        {
            assign_exclusion_convex_to_corner(geomodel, corner_id,
                                              tolerance_shape_type, tolerance);
        }

        for (auto line_id : RINGMesh::range(geomodel.nb_lines()))
        {
            assign_exclusion_convex_to_line(geomodel, line_id, tolerance_shape_type,
                                            tolerance, angle);
        }
    }

    void assign_exclusion_convex_to_entities(
        RINGMesh::GeoModel2D &geomodel,
        const ConvexShapeType &tolerance_shape_type,
        const double tolerance,
        const double angle)
    {
        assign_exclusion_convex_to_entities_base<2>(geomodel,
                                                    tolerance_shape_type, tolerance, angle);
    }

    void assign_exclusion_convex_to_entities(
        RINGMesh::GeoModel3D &geomodel,
        const ConvexShapeType &tolerance_shape_type,
        const double tolerance,
        const double angle)
    {
        assign_exclusion_convex_to_entities_base<3>(geomodel,
                                                    tolerance_shape_type, tolerance, angle);

        for (auto surface_id : RINGMesh::range(geomodel.nb_surfaces()))
        {
            assign_exclusion_convex_to_surface(geomodel, surface_id,
                                               tolerance_shape_type, tolerance, angle);
        }
    }

    void retrieve_volumetric_entities(RINGMesh::GeoModel2D &output_geomodel)
    {
        RINGMesh::GeoModelBuilder2D builder(output_geomodel);
        builder.build_surfaces_from_corners_and_lines();
    }

    void retrieve_volumetric_entities(RINGMesh::GeoModel3D &output_geomodel)
    {
        RINGMesh::GeoModelBuilder3D builder(output_geomodel);
        builder.build_regions_from_lines_and_surfaces();
    }

    template <index_t DIMENSION>
    void repair_geomodel(
        RINGMesh::GeoModel<DIMENSION> &input_geomodel,
        const ConvexShapeType &tolerance_shape_type,
        const double &tolerance,
        const double &angle,
        const GeologicalRule &geological_rules,
        RINGMesh::GeoModel<DIMENSION> &output_geomodel,
        bool retrieve_volumes,
        bool verbose)
    {
        GEO::Stopwatch timer("Elapsed Time");
        Logger::out("Step1", "Repaired geomodel topology recovery...");
        assign_exclusion_convex_to_entities(input_geomodel, tolerance_shape_type,
                                            tolerance, angle);
        GeoModelTopologyRecoverer<DIMENSION> topology_recoverer(input_geomodel,
                                                                tolerance_shape_type, verbose);
        topology_recoverer.add_geological_rule(geological_rules);
        topology_recoverer.perform_graph_based_repair();

        if (DIMENSION == 3)
        {
            Logger::err("3D",
                        "Graph based repair and simplification not implemented in 3D.");
            return;
        }

        if (enum_contains(geological_rules, GeologicalRule::CONSTANT_TOPOLOGY))
        {
            Logger::out("Step2", "Constant topology - "
                                 "Copy geomodel topology and geometry...");
            RINGMesh::copy_geomodel(input_geomodel, output_geomodel);

            Logger::out("Step3", "Constant topology - "
                                 "Geometrical edition...");
            GeoModelVertexShifter<DIMENSION> geomodel_remesher(output_geomodel,
                                                               topology_recoverer, tolerance_shape_type, angle, verbose);
            geomodel_remesher.shift_model_vertices();

            timer.elapsed_time();
        }
        else
        {
            Logger::out("Step2", "Non constant topology - "
                                 "Creation and validation of repair geomodel topology...");
            GeoModelTopologyMaker<DIMENSION> topology_maker(topology_recoverer,
                                                            output_geomodel, verbose);
            topology_maker.build_geomodel_topology();

            Logger::out("Step3", "Non constant topology - "
                                 "Repaired geomodel entity remeshing...");
            GeoModelRemesher<DIMENSION> geomodel_remesher(output_geomodel,
                                                          input_geomodel, topology_maker, verbose);
            geomodel_remesher.remesh_geomodel_entities();

            timer.elapsed_time();
            Logger::out("Step4", "Post process: finding surface and cleaning...");
            clean_geomodel(output_geomodel, retrieve_volumes);
        }
    }

    template <index_t DIMENSION>
    void delete_duplicated_corners(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        index_t kept_corner_id,
        const std::vector<RINGMesh::GMEVertex> &duplicated_corner_vertices)
    {
        scar_assert(duplicated_corner_vertices.size() > 1);
        RINGMesh::GeoModelBuilder<DIMENSION> builder(geomodel);
        std::set<RINGMesh::gmme_id> to_delete;
        for (const auto &dupl_corner : duplicated_corner_vertices)
        {
            if (dupl_corner.gmme.index() == kept_corner_id)
            {
                continue;
            }

            for (auto i : RINGMesh::range(
                     geomodel.corner(dupl_corner.gmme.index()).nb_incident_entities()))
            {
                const auto &incident_line = geomodel.corner(
                                                        dupl_corner.gmme.index())
                                                .incident_entity(i);
                scar_assert(incident_line.nb_boundaries() == 2);
                builder.topology.remove_mesh_entity_boundary_relation(
                    incident_line.gmme(), dupl_corner.gmme);
                builder.topology.add_line_corner_boundary_relation(
                    incident_line.index(), kept_corner_id);
            }

            to_delete.insert(dupl_corner.gmme);
            builder.remove.remove_mesh_entities(to_delete);
        }
    }

    template <index_t DIMENSION>
    void scar_api clean_geomodel(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        bool retrieve_volumes)
    {
        RINGMesh::GeoModelBuilder<DIMENSION> builder(geomodel);
        Logger::out("Step0", "Clean colocated line vertices...");
        for (auto v_id : RINGMesh::range(geomodel.mesh.vertices.nb()))
        {
            auto corresponding_line_vertices =
                geomodel.mesh.vertices.gme_type_vertices(
                    RINGMesh::line_type_name_static(), v_id);
            if (corresponding_line_vertices.size() >= 2)
            {
                auto corresponding_corner_vertices =
                    geomodel.mesh.vertices.gme_type_vertices(
                        RINGMesh::corner_type_name_static(), v_id);
                if (corresponding_corner_vertices.empty())
                {
                    for (const auto &line_vertex : corresponding_line_vertices)
                    {
                        scar_assert(line_vertex.v_index > 0);
                        auto previous_vertex_location = geomodel.line(
                                                                    line_vertex.gmme.index())
                                                            .vertex(
                                                                line_vertex.v_index - 1);
                        builder.geometry.set_mesh_entity_vertex(line_vertex.gmme,
                                                                line_vertex.v_index, previous_vertex_location, false);
                    }
                }
            }
        }
        geomodel.mesh.vertices.clear();

        Logger::out("Step0bis", "Remove duplicated line vertices...");
        for (const auto &line : geomodel.lines())
        {
            std::vector<vecn<DIMENSION>> cleaned_line_vertices;
            cleaned_line_vertices.push_back(line.vertex(0));
            for (auto v : RINGMesh::range(1, line.nb_vertices()))
            {
                if (line.vertex(v) != cleaned_line_vertices.back())
                {
                    cleaned_line_vertices.push_back(line.vertex(v));
                }
            }
            if (cleaned_line_vertices.size() < line.nb_vertices())
            {
                builder.geometry.set_line(line.index(), cleaned_line_vertices);
            }
        }
        std::set<RINGMesh::gmme_id> lines_to_remove;
        for (const auto &line : geomodel.lines())
        {
            if (line.nb_vertices() == 1)
            {
                lines_to_remove.insert(line.gmme());
            }
        }
        builder.remove.remove_mesh_entities(lines_to_remove);
        geomodel.mesh.vertices.clear();

        Logger::out("Step1", "Check Corners/Lines connectivity and geometry...");
        for (auto c_id : RINGMesh::range(geomodel.nb_corners()))
        {
            const auto &corner = geomodel.corner(c_id);
            bool modified = true;
            while (modified)
            {
                modified = false;
                auto geomodel_v_id = geomodel.mesh.vertices.geomodel_vertex_id(
                    corner.gmme());
                auto corresponding_line_vertices =
                    geomodel.mesh.vertices.gme_type_vertices(
                        RINGMesh::line_type_name_static(), geomodel_v_id);
                for (const auto &line_v_id : corresponding_line_vertices)
                {
                    if (line_v_id.v_index == 0)
                    {
                        continue;
                    }
                    if (line_v_id.v_index == geomodel.line(line_v_id.gmme.index()).nb_vertices() - 1)
                    {
                        continue;
                    }
                    split_line_in_two_lines(geomodel, line_v_id.gmme.index(),
                                            corner.vertex(0));
                    modified = true;
                    break;
                }
            }
        }
        geomodel.mesh.vertices.clear();

        Logger::out("Step2", "Remove duplicated entities...");
        index_t corner_id = 0;
        index_t nb_corners = geomodel.nb_corners();
        index_t nb_corners_init = geomodel.nb_corners();

        while (corner_id < nb_corners)
        {
            const auto &corner = geomodel.corner(corner_id);
            ++corner_id;
            auto geomodel_v_id = geomodel.mesh.vertices.geomodel_vertex_id(
                corner.gmme());
            auto corresponding_corner_vertices =
                geomodel.mesh.vertices.gme_type_vertices(
                    RINGMesh::corner_type_name_static(), geomodel_v_id);
            if (corresponding_corner_vertices.size() == 1)
            {
                if (corresponding_corner_vertices[0].gmme == corner.gmme())
                {
                    continue;
                }
            }
            delete_duplicated_corners(geomodel, corner.index(),
                                      corresponding_corner_vertices);
            nb_corners = geomodel.nb_corners();
        }
        Logger::out("Step2", nb_corners_init - nb_corners, " Corner(s) removed.");

        std::set<RINGMesh::gmme_id> to_delete;
        for (const auto &line_id : RINGMesh::range(geomodel.nb_lines() - 1))
        {
            for (auto other_line_id : RINGMesh::range(line_id + 1,
                                                      geomodel.nb_lines()))
            {
                if (are_lines_identical(geomodel.line(line_id),
                                        geomodel.line(other_line_id)))
                {
                    to_delete.insert({RINGMesh::line_type_name_static(),
                                      other_line_id});
                }
            }
        }
        builder.remove.remove_mesh_entities(to_delete);
        geomodel.mesh.vertices.clear();

        Logger::out("Step3", "Remove useless corners");
        // Valence 2 or 0
        remove_corners_with_valence_2(geomodel);
        geomodel.mesh.vertices.clear();

        to_delete.clear();

        if (retrieve_volumes)
        {
            try
            {
                Logger::out("Step4", "Retrieve volumetric entities...");
                for (const auto &surface : geomodel.surfaces())
                {
                    to_delete.insert(surface.gmme());
                }
                builder.remove.remove_mesh_entities(to_delete);
                retrieve_volumetric_entities(geomodel);
                geomodel.mesh.vertices.clear();
            }
            catch (const SCARException &e)
            {
                Logger::err(e.category(), e.what());
                Logger::warn("Exception", "Aborting");
                geomodel.mesh.vertices.clear();
            }
        }
    }

    template void scar_api repair_geomodel(
        RINGMesh::GeoModel2D &,
        const ConvexShapeType &,
        const double &,
        const double &,
        const GeologicalRule &,
        RINGMesh::GeoModel2D &,
        bool,
        bool);

    template void scar_api clean_geomodel(RINGMesh::GeoModel2D &, bool);

    template void scar_api repair_geomodel(
        RINGMesh::GeoModel3D &,
        const ConvexShapeType &,
        const double &,
        const double &,
        const GeologicalRule &,
        RINGMesh::GeoModel3D &,
        bool,
        bool);

    template void scar_api clean_geomodel(RINGMesh::GeoModel3D &, bool);
}
