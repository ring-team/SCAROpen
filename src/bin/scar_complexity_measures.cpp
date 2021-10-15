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

#include <scar/basic/common.h>

#include <geogram/basic/command_line.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/io/io.h>

#include <scar/basic/command_line.h>
#include <scar/tools/complexity_measures.h>

using namespace SCAR;

namespace
{

    void hello()
    {
        Logger::div("SCAR");
        Logger::out("", "Welcome to SCAR !");
        Logger::out("", "Sealed model Creation and Automatic Repair");
    }

    void import_arg_group_quality()
    {
        GEO::CmdLine::declare_arg_group("complexity", "Model complexity measures");
        GEO::CmdLine::declare_arg("complexity:mode", 0,
                                  "Complexity measures mode (1 for global measures, 2 for local measures "
                                  "(using RINGGraph), 0 for ALL)");
        GEO::CmdLine::declare_arg("complexity:min_size", 1.,
                                  "Minimum value of the feature sizes to compare with (GLOBAL MEASURES)");
        GEO::CmdLine::declare_arg("complexity:min_angle", 15.,
                                  "Minimum value of the feature angles (in degrees) to "
                                  "compare with (GLOBAL MEASURES)");
        GEO::CmdLine::declare_arg("complexity:grid_resolution", 5.,
                                  "Size of the grid cells (LOCAL MEASURES)");
        GEO::CmdLine::declare_arg("complexity:output", "",
                                  "Output path for info on complexity measures");
    }

    void import_arg_groups()
    {
        CmdLine::import_arg_group("in");
        import_arg_group_quality();
        CmdLine::import_arg_group("out");
    }

    template <index_t DIMENSION>
    void run(const std::string &geomodel_in_name)
    {
        auto geomodel_out_name = GEO::CmdLine::get_arg("out:geomodel");
        if (geomodel_out_name.empty())
        {
            throw SCARException("I/O", "Give at least a filename in out:geomodel");
        }

        RINGMesh::GeoModel<DIMENSION> geomodel;
        geomodel_load(geomodel, geomodel_in_name);
        auto filepath = GEO::CmdLine::get_arg("complexity:output");
        if (filepath.empty())
        {
            throw SCARException("I/O",
                                "Give the filepath to output complexity information");
        }

        auto complexity_mode = GEO::CmdLine::get_arg_uint("complexity:mode");
        if (complexity_mode == 1 || complexity_mode == 0)
        {
            auto min_size = GEO::CmdLine::get_arg_double("complexity:min_size");
            auto min_angle = GEO::CmdLine::get_arg_double("complexity:min_angle");
            std::string filename = filepath + "/global_measures.txt";
            compute_model_complexity_global_measures(geomodel, min_size, min_angle,
                                                     filename);
        }
        if (complexity_mode == 2 || complexity_mode == 0)
        {
            auto grid_resolution = GEO::CmdLine::get_arg_double(
                "complexity:grid_resolution");
            std::string filename = filepath + "/local_measures.txt";
            compute_model_complexity_local_measures(geomodel, grid_resolution,
                                                    filename);
        }
        geomodel_save(geomodel, geomodel_out_name);
    }
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
        auto geomodel_in_name = GEO::CmdLine::get_arg("in:geomodel");
        if (geomodel_in_name.empty())
        {
            throw SCARException("I/O", "Give at least a filename in in:geomodel");
        }
        index_t dimension = RINGMesh::find_geomodel_dimension(geomodel_in_name);

        if (dimension == 2)
        {
            run<2>(geomodel_in_name);
        }
        else if (dimension == 3)
        {
            run<3>(geomodel_in_name);
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
