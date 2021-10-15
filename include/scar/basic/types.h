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

#include <ringmesh/basic/types.h>

/*!
 * @file Re-defintions of basic types similar to those of Geogram
 */

namespace SCAR
{
    /* If you need integer of 8bits of any other one
     * it is sufficient to write using GEO::Numeric::uint8 in your file.
     *
     * Dummy variables were removed, the pollute the namespace and
     * it is quite easy to do without them.
     */

    // Basic types used in SCAR
    inline double max_double()
    {
        return std::numeric_limits<double>::max();
    }

    // Using definitions of RINGMesh
    using RINGMesh::global_epsilon;
    using RINGMesh::global_epsilon_3;
    using RINGMesh::global_epsilon_sq;

    // This is an unsigned int
    using RINGMesh::index_t;
    // This is an int
    using RINGMesh::signed_index_t;

    // This is an array template of doubles
    template <index_t DIMENSION>
    using vecn = RINGMesh::vecn<DIMENSION>;
    using vec2 = vecn<2>;
    using vec3 = vecn<3>;

    // This is the value used in SCAR for a invalid index and invalid point
    using RINGMesh::NO_ID;

    template <index_t DIMENSION>
    inline const vecn<DIMENSION> NO_POINT()
    {
        vecn<DIMENSION> no_point;
        for (auto c : RINGMesh::range(DIMENSION))
        {
            no_point[c] = -999999.9;
        }
        return no_point;
    }
    static const vec2 NO_POINT_2D = NO_POINT<2>();
    static const vec3 NO_POINT_3D = NO_POINT<3>();

    template <index_t DIMENSION>
    bool operator==(const vecn<DIMENSION> &u, const vecn<DIMENSION> &v)
    {
        for (auto i : RINGMesh::range(DIMENSION))
        {
            if (u[i] != v[i])
                return false;
        }
        return true;
    }

    template <index_t DIMENSION>
    bool operator!=(const vecn<DIMENSION> &u, const vecn<DIMENSION> &v)
    {
        return !(u == v);
    }

    const std::string tolerance_att_name = "tolerance";
    const std::string ellipse_major_ppl_att_name = "ell_major_ppl";
    const std::string ellipse_minor_ppl_att_name = "ell_minor_ppl";
    const std::string geodesic_att_name = "geodesic";

    // RINGMesh bitwise operators
    template <typename Enum>
    auto to_underlying_type(Enum e) ->
        typename std::underlying_type<Enum>::type
    {
        return static_cast<typename std::underlying_type<Enum>::type>(e);
    }

    template <typename Enum>
    struct EnableBitMaskOperators
    {
        static const bool enable = false;
    };

    template <typename Enum>
    typename std::enable_if<EnableBitMaskOperators<Enum>::enable,
                            Enum>::type
    operator|(Enum lhs, Enum rhs)
    {
        using underlying = typename std::underlying_type<Enum>::type;
        return static_cast<Enum>(static_cast<underlying>(lhs) | static_cast<underlying>(rhs));
    }

    template <typename Enum>
    typename std::enable_if<EnableBitMaskOperators<Enum>::enable,
                            Enum>::type
    operator&(Enum lhs, Enum rhs)
    {
        using underlying = typename std::underlying_type<Enum>::type;
        return static_cast<Enum>(static_cast<underlying>(lhs) & static_cast<underlying>(rhs));
    }

    template <typename Enum>
    typename std::enable_if<EnableBitMaskOperators<Enum>::enable,
                            Enum>::type
    operator^(Enum lhs, Enum rhs)
    {
        using underlying = typename std::underlying_type<Enum>::type;
        return static_cast<Enum>(static_cast<underlying>(lhs) ^ static_cast<underlying>(rhs));
    }

    template <typename Enum>
    typename std::enable_if<EnableBitMaskOperators<Enum>::enable,
                            Enum>::type
    operator~(Enum lhs)
    {
        using underlying = typename std::underlying_type<Enum>::type;
        return static_cast<Enum>(~static_cast<underlying>(lhs));
    }

    template <typename Enum>
    bool enum_contains(Enum lhs, Enum rhs)
    {
        return (lhs & rhs) != Enum::EMPTY;
    }

#define ENABLE_BITMASK_OPERATORS(Enum)   \
    template <>                          \
    struct EnableBitMaskOperators<Enum>  \
    {                                    \
        static const bool enable = true; \
    }
}

namespace RINGpcl
{
    // This is an unsigned int
    using RINGMesh::index_t;
}
