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

#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/io/geomodel_builder_gocad.h>

namespace SCAR
{
    class GeoModelBuilderTS;
}

namespace SCAR
{
    void scar_api tsurf_import_factory_initialize();

    struct TSurfLoadingStorage : public RINGMesh::GocadLoadingStorage
    {
        TSurfLoadingStorage();

        // Gocad index of the first vertex of current surface
        index_t cur_surf_first_vertex_index_;
    };

    class scar_api TSurfLineParser : public RINGMesh::GocadBaseParser
    {
        scar_disable_copy(TSurfLineParser);

    public:
        virtual ~TSurfLineParser() = default;
        virtual void execute(
            GEO::LineInput &line,
            TSurfLoadingStorage &load_storage) = 0;

    protected:
        TSurfLineParser(
            GeoModelBuilderTS &gm_builder,
            RINGMesh::GeoModel3D &geomodel);
    };

    using TSurfLineParserFactory = RINGMesh::Factory<std::string, TSurfLineParser, GeoModelBuilderTS &, RINGMesh::GeoModel3D &>;

    /*!
     * @brief Build a GeoModel from a Gocad set of surfaces (surfaces.ts)
     */
    class scar_api GeoModelBuilderTS : public RINGMesh::GeoModelBuilderGocad
    {
    public:
        GeoModelBuilderTS(
            RINGMesh::GeoModel3D &geomodel,
            const std::string &filename)
            : RINGMesh::GeoModelBuilderGocad(geomodel, std::move(filename))
        {
        }

        virtual ~GeoModelBuilderTS() = default;

        static void save_geomodel(
            const RINGMesh::GeoModel2D &geomodel,
            const std::string &filename);

        static void save_geomodel(
            const RINGMesh::GeoModel3D &geomodel,
            const std::string &filename);

    private:
        void load_file();

        /*!
         * @brief Reads the first word of the current line (keyword)
         * and executes the good action with the information of the line
         * @details Uses the MLLineParser factory
         */
        virtual void read_line();

        void remove_duplicate_lines();
        void find_duplicate_lines(
            std::set<RINGMesh::gmme_id> &duplicate_lines_to_remove);
        bool is_line_duplicated(const index_t line_id);
        bool are_lines_identical(
            const RINGMesh::Line3D &line1,
            const RINGMesh::Line3D &line2);
        bool same_number_vertices(
            const RINGMesh::Line3D &line1,
            const RINGMesh::Line3D &line2);
        bool same_boundaries(
            const RINGMesh::Line3D &line1,
            const RINGMesh::Line3D &line2);
        bool same_internal_vertices(
            const RINGMesh::Line3D &line1,
            const RINGMesh::Line3D &line2);

        void retrieve_geomodel_topology();

    private:
        TSurfLoadingStorage ts_load_storage_;
    };

}
