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

#include <scar/tools/convex_shape.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModel);
    ALIAS_2D_AND_3D(GeoModel);
}

namespace SCAR
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModelTopologyRecoverer);
}

namespace SCAR
{
    template <index_t DIMENSION>
    class scar_api GeoModelVertexShifter
    {
        scar_disable_copy(GeoModelVertexShifter);

    public:
        GeoModelVertexShifter(
            RINGMesh::GeoModel<DIMENSION> &output_geomodel,
            const GeoModelTopologyRecoverer<DIMENSION> &topology_recoverer,
            ConvexShapeType convex_type,
            const double &angle,
            bool verbose);

        void shift_model_vertices();

    private:
        void edit_corners();
        void repel_corners_from_lines();
        void edit_lines();

    private:
        RINGMesh::GeoModel<DIMENSION> &geomodel_;
        const GeoModelTopologyRecoverer<DIMENSION> &topology_recoverer_;
        ConvexShapeType shape_type_;
        double angle_;
        bool verbose_;
    };

}
