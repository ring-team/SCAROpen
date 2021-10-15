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
#include <ringmesh/geomodel/builder/geomodel_builder.h>

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
    GEO::CmdLine::declare_arg_group("mult", "Multipliers for each coordinates "
                                            "(. = 1: no modification, . > 1: dilation, "
                                            "0 > . > 1: shrinkage, . < -1: rotation)");
    GEO::CmdLine::declare_arg("mult:x", 1, "Multiplier for x coordinates");
    GEO::CmdLine::declare_arg("mult:y", 1, "Multiplier for y coordinates");
    GEO::CmdLine::declare_arg("mult:z", 1, "Multiplier for z coordinates");
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
        Logger::out("Step1", "Saving repaired geomodel...");
        geomodel_save(output_geomodel, geomodel_out_name);
    }
}

void dilate_geomodel(
    RINGMesh::GeoModel2D &geomodel,
    const std::vector<double> &multipliers)
{
    scar_assert(multipliers.size() == 2);
    RINGMesh::GeoModelBuilder2D gm_builder(geomodel);
    for (auto v_id : RINGMesh::range(geomodel.mesh.vertices.nb()))
    {
        auto old_point = geomodel.mesh.vertices.vertex(v_id);
        vec2 new_point(multipliers[0] * old_point.x, multipliers[1] * old_point.y);
        gm_builder.geometry.set_mesh_entity_vertex(v_id, new_point);
    }
}

void dilate_geomodel(
    RINGMesh::GeoModel3D &geomodel,
    const std::vector<double> &multipliers)
{
    scar_assert(multipliers.size() == 3);
    RINGMesh::GeoModelBuilder3D gm_builder(geomodel);
    for (auto v_id : RINGMesh::range(geomodel.mesh.vertices.nb()))
    {
        auto old_point = geomodel.mesh.vertices.vertex(v_id);
        vec3 new_point(multipliers[0] * old_point.x, multipliers[1] * old_point.y,
                       multipliers[2] * old_point.z);
        gm_builder.geometry.set_mesh_entity_vertex(v_id, new_point);
    }
}

template <index_t DIMENSION>
void load_and_dilate_geomodel(
    const std::string &filename,
    const std::vector<double> &multipliers)
{
    RINGMesh::GeoModel<DIMENSION> geomodel_in;
    geomodel_load(geomodel_in, filename);

    dilate_geomodel(geomodel_in, multipliers);

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

        double mult_x = GEO::String::to_double(GEO::CmdLine::get_arg("mult:x"));
        double mult_y = GEO::String::to_double(GEO::CmdLine::get_arg("mult:y"));
        std::vector<double> multipliers;
        multipliers.push_back(mult_x);
        multipliers.push_back(mult_y);

        if (dimension == 2)
        {
            load_and_dilate_geomodel<2>(geomodel_in_name, multipliers);
        }
        else if (dimension == 3)
        {
            double mult_z = GEO::String::to_double(
                GEO::CmdLine::get_arg("mult:z"));
            multipliers.push_back(mult_z);
            load_and_dilate_geomodel<3>(geomodel_in_name, multipliers);
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
