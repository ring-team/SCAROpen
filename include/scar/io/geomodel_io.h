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

#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

#include <scar/io/geomodel_builder_ts.h>
#include <scar/io/geomodel_builder_pl.h>
#include <scar/io/geomodel_builder_msh.h>

namespace SCAR
{
    /* These classes are extensions of RINGMesh IO factories 
     * used for developping SCAR 
     */

    class TSIOHandler3D final : public RINGMesh::GeoModelInputHandler3D,
                                public RINGMesh::GeoModelOutputHandler3D
    {
    public:
        /*! Load a .ts (Gocad file)
         * @pre Filename is valid
         */
        void load(const std::string &filename, RINGMesh::GeoModel3D &geomodel) final
        {
            std::ifstream input(filename.c_str());
            if (!input)
            {
                throw SCARException("I/O", "Failed to open file " + filename);
            }

            GeoModelBuilderTS builder(geomodel, filename);
            builder.build_geomodel();
        }

        virtual void save(
            const RINGMesh::GeoModel3D &geomodel,
            const std::string &filename)
        {
            scar_unused(geomodel);
            scar_unused(filename);
            GeoModelBuilderTS::save_geomodel(geomodel, filename);
        }
    };

    class TSIOHandler2D final : public RINGMesh::GeoModelOutputHandler2D
    {
    public:
        virtual void save(
            const RINGMesh::GeoModel2D &geomodel,
            const std::string &filename)
        {
            GeoModelBuilderTS::save_geomodel(geomodel, filename);
        }
    };

#define ringmesh_register_geomodel_io_creator(type, name) \
    geo_register_creator(RINGMesh::GeoModelIOHandlerFactory, type, name)

    class PLIOHandler3D final : public RINGMesh::GeoModelInputHandler3D,
                                public RINGMesh::GeoModelOutputHandler3D
    {
    public:
        /*! Load a .pl (Gocad file)
         * @pre Filename is valid
         */
        void load(const std::string &filename, RINGMesh::GeoModel3D &geomodel) final
        {
            std::ifstream input(filename.c_str());
            if (!input)
            {
                throw SCARException("I/O", "Failed to open file " + filename);
            }

            GeoModelBuilderPL builder(geomodel, filename);
            builder.build_geomodel();
        }

        virtual void save(
            const RINGMesh::GeoModel3D &geomodel,
            const std::string &filename)
        {
            GeoModelBuilderPL::save_geomodel(geomodel, filename);
        }
    };

    class PLIOHandler2D final : public RINGMesh::GeoModelOutputHandler2D
    {
    public:
        virtual void save(
            const RINGMesh::GeoModel2D &geomodel,
            const std::string &filename)
        {
            GeoModelBuilderPL::save_geomodel(geomodel, filename);
        }
    };

    class MSHIOHandler2D final : public RINGMesh::GeoModelOutputHandler2D
    {
    public:
        virtual void save(
            const RINGMesh::GeoModel2D &geomodel,
            const std::string &filename)
        {
            GeoModelBuilderMSH::save_geomodel(geomodel, filename);
        }
    };

    void register_builder_scar_input()
    {
        // 3D GeoModels
        RINGMesh::GeoModelInputHandlerFactory3D::register_creator<TSIOHandler3D>(
            "ts");
        RINGMesh::GeoModelInputHandlerFactory3D::register_creator<PLIOHandler3D>(
            "pl");
    }
    void register_builder_scar_output()
    {
        // 2D GeoModels
        RINGMesh::GeoModelOutputHandlerFactory2D::register_creator<PLIOHandler2D>(
            "pl");
        RINGMesh::GeoModelOutputHandlerFactory2D::register_creator<TSIOHandler2D>(
            "ts");
        RINGMesh::GeoModelOutputHandlerFactory2D::register_creator<MSHIOHandler2D>(
            "msh");

        // 3D GeoModels
        RINGMesh::GeoModelOutputHandlerFactory3D::register_creator<PLIOHandler3D>(
            "pl");
        RINGMesh::GeoModelOutputHandlerFactory3D::register_creator<TSIOHandler3D>(
            "ts");
    }
}
