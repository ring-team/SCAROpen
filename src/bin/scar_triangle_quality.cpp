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
#include <ringmesh/mesh/surface_mesh.h>

#include <scar/basic/command_line.h>
#include <scar/tools/triangle_quality.h>

using namespace SCAR;

namespace
{

    void hello()
    {
        Logger::div("SCAR");
        Logger::out("", "Welcome to SCAR !");
        Logger::out("", "Sealed model Creation and Automatic Repair");
        Logger::out("", "People working on the project in RING");
        Logger::out("", "Pierre Anquez <pierre.anquez@univ-lorraine.fr> ");
    }

    void import_arg_group_quality()
    {
        GEO::CmdLine::declare_arg_group("quality", "Triangular mesh quality");
        GEO::CmdLine::declare_arg("quality:mode", 0,
                                  "Mesh quality mode (0 computes ALL)");
        GEO::CmdLine::declare_arg("quality:min_value", 0.05,
                                  "Cell quality is defined as low if below this minimum value");
        GEO::CmdLine::declare_arg("quality:output", "",
                                  "Output filename for a mesh containing low quality triangles");
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

        auto quality_mode = GEO::CmdLine::get_arg_uint("quality:mode");
        compute_prop_triangle_mesh_quality(
            static_cast<TriangleQualityMode>(quality_mode), geomodel);

        auto histo_out_name = GEO::CmdLine::get_arg("quality:output");
        if (!histo_out_name.empty())
        {
            output_quality_values(histo_out_name,
                                  static_cast<TriangleQualityMode>(quality_mode), geomodel);
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
