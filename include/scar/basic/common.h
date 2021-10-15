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

/* Include the configuration file generated by CMAKE
 from cmake config options and the
 include/scar/config_scar.h.in file.
 File is in SCAR_BIN/scar directory.
 */
#include <SCAR_config.h>
#include <scar/export.h>

#include <ringmesh/basic/common.h>
#include <memory>

#if defined(_WIN32)
#ifndef WIN32
#define WIN32
#endif
#endif

#ifndef NDEBUG
#define SCAR_DEBUG
#else
#undef SCAR_DEBUG
#endif

#define scar_disable_copy(Class) \
    ringmesh_disable_copy(Class)

#define scar_template_assert_2d_or_3d(type) \
    ringmesh_template_assert_2d_or_3d(type)

#ifdef WIN32
#pragma warning(disable : 4267)
#endif

#ifdef SCAR_DEBUG
#define scar_assert(x)                                               \
    if (!(x))                                                        \
    {                                                                \
        RINGMesh::ringmesh_assertion_failed(#x, __FILE__, __LINE__); \
    }
#define scar_assert_not_reached \
    RINGMesh::ringmesh_should_not_have_reached(__FILE__, __LINE__);
#else
#define scar_assert(x)
#define scar_assert_not_reached
#endif

using SCARException = RINGMesh::RINGMeshException;

#define scar_unused(x) \
    ringmesh_unused(x)

#include <scar/basic/types.h>
#include <omp.h>

#include <ringmesh/basic/logger.h>

using Logger = RINGMesh::Logger;

template <typename Container>
inline void vector_DEBUG(const Container &c)
{
    for (auto i : c)
    {
        DEBUG(i);
    }
}

template <typename Map>
inline void map_DEBUG(const Map &c)
{
    for (auto i : c)
    {
        Logger::out("MapDebug", i.first, " ", i.second);
    }
}

#define VECDEBUG(c)                          \
    Logger::out("Debug", #c, " elements: "); \
    vector_DEBUG(c);

#define MAPDEBUG(c)                             \
    Logger::out("MapDebug", #c, " elements: "); \
    map_DEBUG(c);

namespace SCAR
{
    void scar_api configure_scar();
}

#ifdef SCAR_DEBUG
#define ringpcl_assert(x)                                            \
    if (!(x))                                                        \
    {                                                                \
        RINGMesh::ringmesh_assertion_failed(#x, __FILE__, __LINE__); \
    }
#define ringpcl_assert_not_reached \
    RINGMesh::ringmesh_should_not_have_reached(__FILE__, __LINE__);
#else
#define ringpcl_assert(x)
#define ringpcl_assert_not_reached
#endif
