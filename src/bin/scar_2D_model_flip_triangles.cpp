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

#include <ringmesh/basic/box.h>
#include <iostream>

#include <geogram/basic/command_line.h>

#include <ringmesh/io/io.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/mesh/mesh_index.h>

#include <scar/basic/command_line.h>
#include <scar/remeshing/geomodel_entity_remeshing.h>

using namespace SCAR;

void hello()
{
    Logger::div("SCAR");
    Logger::out("", "Welcome to SCAR !");
    Logger::out("", "Sealed model Creation and Automatic Repair");
}

void import_arg_groups()
{
    CmdLine::import_arg_group("global");
    CmdLine::import_arg_group("in");
    CmdLine::import_arg_group("out");
}

template <index_t DIMENSION>
void save_output_model(const RINGMesh::GeoModel<DIMENSION> &output_geomodel)
{
    std::string geomodel_out_name{GEO::CmdLine::get_arg("out:geomodel")};
    if (geomodel_out_name.empty())
    {
        throw SCARException("I/O", "No filename in out:geomodel");
    }
    else
    {
        Logger::out("Step1", "Saving geomodel...");
        geomodel_save(output_geomodel, geomodel_out_name);
    }
}

void check_geomodel_normals(RINGMesh::GeoModel2D &geomodel)
{
    const auto &triangles = geomodel.mesh.polygons;
    const auto &points = geomodel.mesh.vertices;

    // Init
    auto v0_to_v1 = GEO::normalize(
        points.vertex(triangles.vertex({0, 1})) - points.vertex(triangles.vertex({0, 0})));
    auto v0_to_v2 = GEO::normalize(
        points.vertex(triangles.vertex({0, 2})) - points.vertex(triangles.vertex({0, 0})));
    vec3 v01{v0_to_v1.x, v0_to_v1.y, 0};
    vec3 v02{v0_to_v2.x, v0_to_v2.y, 0};
    bool positive = cross(v01, v02).z > 0;

    for (auto i : RINGMesh::range(triangles.nb_triangle()))
    {
        auto v0_to_v1 = GEO::normalize(
            points.vertex(triangles.vertex({i, 1})) - points.vertex(triangles.vertex({i, 0})));
        auto v0_to_v2 = GEO::normalize(
            points.vertex(triangles.vertex({i, 2})) - points.vertex(triangles.vertex({i, 0})));
        vec3 v01{v0_to_v1.x, v0_to_v1.y, 0};
        vec3 v02{v0_to_v2.x, v0_to_v2.y, 0};

        if ((cross(v01, v02).z > 0) != positive)
        {
            Logger::err("Normals",
                        "All the triangles are not oriented the same way. Abort!");
        }
    }
    if (positive)
    {
        Logger::out("Normals", "From DIRECT into INDIRECT");
    }
    else
    {
        Logger::out("Normals", "From INDIRECT into DIRECT");
    }
}

void flip_geomodel(RINGMesh::GeoModel2D &geomodel)
{
    RINGMesh::GeoModelBuilder2D gm_builder(geomodel);

    //    const auto& triangles = geomodel.mesh.polygons;
    //    const auto& points = geomodel.mesh.vertices;

    for (const auto &surface : geomodel.surfaces())
    {
        for (auto triangle : RINGMesh::range(surface.nb_mesh_elements()))
        {
            std::vector<index_t> reordered_corners(3);
            reordered_corners[0] = surface.mesh_element_vertex_index(
                {triangle, 1});
            reordered_corners[1] = surface.mesh_element_vertex_index(
                {triangle, 0});
            reordered_corners[2] = surface.mesh_element_vertex_index(
                {triangle, 2});

            gm_builder.geometry.set_surface_element_geometry(surface.index(),
                                                             triangle, reordered_corners);
        }
    }
}

void load_and_flip_geomodel2d(const std::string &filename)
{
    RINGMesh::GeoModel2D geomodel_in;
    geomodel_load(geomodel_in, filename);

    check_geomodel_normals(geomodel_in);
    flip_geomodel(geomodel_in);

    save_output_model(geomodel_in);
}

int main(int argc, char **argv)
{
    try
    {

        hello();
        import_arg_groups();

        if (argc == 1)
        {
            GEO::CmdLine::show_usage();
            return 0;
        }

        std::vector<std::string> filenames;
        if (!GEO::CmdLine::parse(argc, argv, filenames))
        {
            return 1;
        }

        std::string geomodel_in_name{GEO::CmdLine::get_arg("in:geomodel")};
        if (geomodel_in_name.empty())
        {
            throw SCARException("I/O", "Give at least a filename in in:geomodel");
        }

        index_t dimension = RINGMesh::find_geomodel_dimension(geomodel_in_name);

        if (dimension == 2)
        {
            load_and_flip_geomodel2d(geomodel_in_name);
        }
        else
        {
            throw SCARException("I/O", "Can not deal with GeoModel", dimension,
                                "D.");
        }
    }
    catch (const SCARException &e)
    {
        Logger::err(e.category(), e.what());
        return 1;
    }
    catch (const std::exception &e)
    {
        Logger::err("Exception", e.what());
        return 1;
    }
    return 0;
}
