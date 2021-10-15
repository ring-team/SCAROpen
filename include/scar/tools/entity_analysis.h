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

#include <ringmesh/geomodel/core/entity_type.h>

namespace SCAR
{
    template <index_t DIMENSION>
    bool same_number_vertices(
        const RINGMesh::Line<DIMENSION> &line1,
        const RINGMesh::Line<DIMENSION> &line2)
    {
        return line1.nb_vertices() == line2.nb_vertices();
    }

    template <index_t DIMENSION>
    bool same_boundaries(
        const RINGMesh::Line<DIMENSION> &line1,
        const RINGMesh::Line<DIMENSION> &line2)
    {
        scar_assert(line1.nb_boundaries() == 2);
        scar_assert(line2.nb_boundaries() == 2);
        if (line1.boundary_gmme(0).index() != line2.boundary_gmme(0).index() && line1.boundary_gmme(0).index() != line2.boundary_gmme(1).index())
        {
            return false;
        }
        if (line1.boundary_gmme(1).index() != line2.boundary_gmme(0).index() && line1.boundary_gmme(1).index() != line2.boundary_gmme(1).index())
        {
            return false;
        }
        return true;
    }

    template <index_t DIMENSION>
    bool same_internal_vertices(
        const RINGMesh::Line<DIMENSION> &line1,
        const RINGMesh::Line<DIMENSION> &line2)
    {
        bool same_order = true;
        if (line1.boundary_gmme(1).index() == line2.boundary_gmme(0).index() && line1.boundary_gmme(0).index() == line2.boundary_gmme(1).index())
        {
            same_order = false;
        }
        index_t nb_line_vertices = line1.nb_vertices();
        for (index_t v = 1; v < nb_line_vertices - 1; ++v)
        {
            // If same order, compare vertex of same id. Else, we must compare,
            // first of line1 with last (of index nb_line_vertices - 1) of line2,
            // and so on...
            if ((same_order && line1.vertex(v) != line2.vertex(v)) || (!same_order && line1.vertex(v) != line2.vertex((nb_line_vertices - 1) - v)))
            {
                return false;
            }
        }
        return true;
    }

    template <index_t DIMENSION>
    bool are_lines_identical(
        const RINGMesh::Line<DIMENSION> &line1,
        const RINGMesh::Line<DIMENSION> &line2)
    {
        // To be identical 2 lines must have same number of points, same boundaries
        // (corners) and same other points.
        return same_number_vertices(line1, line2) && same_boundaries(line1, line2) && same_internal_vertices(line1, line2);
    }

}
