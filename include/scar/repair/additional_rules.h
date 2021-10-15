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

namespace SCAR
{

    /// Option to select what are the global geological rules for repair
    enum struct GeologicalRule
    {
        EMPTY = 0,
        CONFORMABLE_HORIZONS = 1,
        MERGE_WHOLE_LINES = 1 << 1,       // 2
        NO_LINE_LINE_CHECK = 1 << 2,      // 4
        KEEP_BOUNDARIES_UNMOVED = 1 << 3, // 8
        CONSTANT_TOPOLOGY = 1 << 4,       // 16
        CHECK_INTERSECTIONS = 1 << 5,     // 32
        REMESH_DFN = 1 << 6               // 64
        // OTHER_RULE = 1 << x, etc.
    };
    ENABLE_BITMASK_OPERATORS(GeologicalRule);
}
