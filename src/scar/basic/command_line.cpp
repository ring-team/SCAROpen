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

#include <scar/basic/command_line.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <ringmesh/basic/command_line.h>

namespace SCAR
{

    namespace CmdLine
    {

        void import_arg_group_in()
        {

            GEO::CmdLine::declare_arg_group("in", "Input data");
            GEO::CmdLine::declare_arg("in:geomodel", "",
                                      "Input geological model");
            GEO::CmdLine::declare_arg("in:mesh", "",
                                      "Input surface mesh to remesh.");
        }

        void import_arg_group_out()
        {

            GEO::CmdLine::declare_arg_group("out", "Output data");
            GEO::CmdLine::declare_arg("out:geomodel", "",
                                      "Saves the repaired geological model");
            GEO::CmdLine::declare_arg("out:mesh", "",
                                      "Output remeshed surface mesh.");
        }

        void import_arg_group_tol()
        {
            GEO::CmdLine::declare_arg_group("tol", "Tolerance");
            GEO::CmdLine::declare_arg("tol:percent", 0.01,
                                      "Repairing tolerance express as a percentage of the model "
                                      "diagonal length");
            GEO::CmdLine::declare_arg("tol:distance", "",
                                      "Repairing tolerance expressed in meters");
        }

        void import_arg_group_res()
        {
            GEO::CmdLine::declare_arg_group("res", "Resolution");
            GEO::CmdLine::declare_arg("res:distance", "",
                                      "Mesh element resolution expressed in meters");
        }

        void import_arg_group_angle()
        {
            GEO::CmdLine::declare_arg_group("angle", "Minimum angle");
            GEO::CmdLine::declare_arg("angle:rad", 0.349,
                                      "Minimum angle between repaired model lines express in radians");
            GEO::CmdLine::declare_arg("angle:deg", "",
                                      "Minimum angle between repaired model lines express in degrees");
        }

        void import_arg_group_ellipsoid()
        {
            GEO::CmdLine::declare_arg_group("ellipsoid", "Ellipsoid info");
            GEO::CmdLine::declare_arg("in:ellipsoid", "",
                                      "File containing information about ellipsoid (anisotropic rescaling)");
        }

        bool import_arg_group(const std::string &name)
        {
            if (name == "in")
            {
                import_arg_group_in();
            }
            else if (name == "out")
            {
                import_arg_group_out();
            }
            else if (name == "tol")
            {
                import_arg_group_tol();
            }
            else if (name == "res")
            {
                import_arg_group_res();
            }
            else if (name == "angle")
            {
                import_arg_group_angle();
            }
            else if (name == "ellipsoid")
            {
                import_arg_group_ellipsoid();
            }
            else
            {
                return RINGMesh::CmdLine::import_arg_group(name);
            }
            return true;
        }

    }

}
