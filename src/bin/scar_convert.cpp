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

#include <iostream>

#include <geogram/basic/command_line.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/box.h>
#include <ringmesh/geomodel/builder/geomodel_builder_2d_from_3d.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/io/io.h>

#include <scar/basic/command_line.h>

int main(int argc, char **argv)
{
    using namespace SCAR;

    try
    {

        Logger::div("SCAR");
        Logger::out("", "Sealed model Creation and Automatic Repair");
        Logger::div("SCAR Convert");
        Logger::out("", "Welcome to SCAR Convert !");

        CmdLine::import_arg_group("global");
        CmdLine::import_arg_group("in");
        CmdLine::import_arg_group("out");
        GEO::CmdLine::declare_arg("dim", "", "Input GeoModel dimension");

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

        std::string geomodel_in_name = GEO::CmdLine::get_arg("in:geomodel");
        if (geomodel_in_name.empty())
        {
            throw RINGMesh::RINGMeshException("I/O",
                                              "Give at least a filename in in:geomodel");
        }

        std::string geomodel_out_name = GEO::CmdLine::get_arg("out:geomodel");
        if (geomodel_out_name.empty())
        {
            throw RINGMesh::RINGMeshException("I/O",
                                              "Give at least a filename in out:geomodel");
        }

        std::string dimension_str = GEO::CmdLine::get_arg("dim");
        index_t dimension{3};
        if (!dimension_str.empty())
        {
            dimension = GEO::String::to_uint(dimension_str);
        }

        if (dimension == 2)
        {
            if (RINGMesh::find_geomodel_dimension(geomodel_in_name) == 2)
            {
                RINGMesh::GeoModel2D geomodel;
                // Convert is loading and then saving in new format
                geomodel_load(geomodel, geomodel_in_name);
                geomodel_save(geomodel, geomodel_out_name);
            }
            else
            {
                RINGMesh::GeoModel3D geomodel3d;
                // Convert is loading and then saving in new format
                geomodel_load(geomodel3d, geomodel_in_name);
                RINGMesh::Geometry::Plane plane;

                RINGMesh::Box3D geomodel_box;
                geomodel_box.add_point(geomodel3d.mesh.vertices.vertex(0));
                double a = 0;
                double b = 0;
                double c = 0;
                for (auto v_id : RINGMesh::range(1, geomodel3d.mesh.vertices.nb()))
                {
                    auto v = geomodel3d.mesh.vertices.vertex(v_id);
                    geomodel_box.add_point(v);
                    auto prev_v = geomodel3d.mesh.vertices.vertex(v_id - 1);
                    a += (prev_v.y - v.y) * (v.z - prev_v.z);
                    b += (prev_v.z - v.z) * (v.x - prev_v.x);
                    c += (prev_v.x - v.x) * (v.y - prev_v.y);
                }
                plane.origin = geomodel_box.center();
                // The normal is deduced of an approximation of the
                // mean plane of GeoModel vertices
                if (c < 0)
                {
                    plane.normal = normalize(vec3(-a, -b, -c));
                }
                else
                {
                    plane.normal = normalize(vec3(a, b, c));
                }

                RINGMesh::GeoModel2D geomodel;
                RINGMesh::GeoModelBuilder2DProjection builder2d(geomodel,
                                                                geomodel3d, plane);
                builder2d.build_geomodel();

                geomodel_save(geomodel, geomodel_out_name);
            }
        }
        else if (dimension == 3)
        {
            RINGMesh::GeoModel3D geomodel;
            // Convert is loading and then saving in new format
            geomodel_load(geomodel, geomodel_in_name);
            geomodel_save(geomodel, geomodel_out_name);
        }
        else
        {
            throw SCARException("I/O", "Can not deal with GeoModel", dimension,
                                "D.");
        }
    }
    catch (const RINGMesh::RINGMeshException &e)
    {
        GEO::Logger::err(e.category()) << e.what() << std::endl;
        return 1;
    }
    catch (const std::exception &e)
    {
        GEO::Logger::err("Exception") << e.what() << std::endl;
        return 1;
    }
    return 0;
}
