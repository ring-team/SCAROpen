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

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/io/io.h>

#include <scar/basic/command_line.h>
#include <scar/repair/repairer.h>

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
    GEO::CmdLine::declare_arg("build_volumes", true,
                              "Set false to desactivate automatic build of volumes (or surfaces in 2D)");
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

template <index_t DIMENSION>
void load_and_clean_geomodel(
    const std::string &filename)
{
    RINGMesh::GeoModel<DIMENSION> geomodel_in;
    geomodel_load(geomodel_in, filename);

    bool retrieve_volumes{GEO::CmdLine::get_arg_bool("build_volumes")};

    clean_geomodel(geomodel_in, retrieve_volumes);

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
            load_and_clean_geomodel<2>(geomodel_in_name);
        }
        else if (dimension == 3)
        {
            load_and_clean_geomodel<3>(geomodel_in_name);
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
