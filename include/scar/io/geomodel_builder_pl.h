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
    class GeoModelBuilderPL;
}

namespace SCAR
{

    void scar_api pline_import_factory_initialize();

    struct PLineLoadingStorage : public RINGMesh::GocadLoadingStorage
    {
        PLineLoadingStorage();

        index_t cur_line_first_vertex_index_;
        index_t cur_contact_;
        index_t cur_line_;
    };

    class scar_api PLineLineParser : public RINGMesh::GocadBaseParser
    {
        scar_disable_copy(PLineLineParser);

    public:
        virtual ~PLineLineParser() = default;
        virtual void execute(
            GEO::LineInput &line,
            PLineLoadingStorage &load_storage) = 0;

    protected:
        PLineLineParser(
            GeoModelBuilderPL &gm_builder,
            RINGMesh::GeoModel3D &geomodel);
    };

    using PLineLineParserFactory = RINGMesh::Factory<std::string, PLineLineParser, GeoModelBuilderPL &, RINGMesh::GeoModel3D &>;

    /*!
     * @brief Build a GeoModel from a Gocad set of lines (line.pl)
     */
    class scar_api GeoModelBuilderPL : public RINGMesh::GeoModelBuilderGocad
    {
    public:
        GeoModelBuilderPL(
            RINGMesh::GeoModel3D &geomodel,
            const std::string &filename)
            : RINGMesh::GeoModelBuilderGocad(geomodel, std::move(filename))
        {
        }

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

        void compute_corners();

    private:
        PLineLoadingStorage pl_load_storage_;
    };

}
