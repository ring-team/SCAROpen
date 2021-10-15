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

#pragma once

#include <scar/basic/common.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModel);
    ALIAS_2D(GeoModel);
}

namespace SCAR
{
    /*!
     * @brief Save a GeoModel to GMSH file format (.msh)
     */
    class scar_api GeoModelBuilderMSH
    {
    public:
        static void save_geomodel(
            const RINGMesh::GeoModel2D &geomodel,
            const std::string &filename);

        //        static void save_geomodel(
        //            const RINGMesh::GeoModel3D& geomodel,
        //            const std::string& filename );

    private:
        GeoModelBuilderMSH() = default;
    };

}
