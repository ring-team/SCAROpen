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

#include <scar/io/geomodel_builder_pl.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/builder/geomodel_builder_topology.h>

#include <scar/tools/entity_analysis.h>

namespace
{

    using namespace SCAR;

    // Indices begin to 1 in Gocad
    index_t GOCAD_OFFSET = 1;

    void set_cur_line_mesh_from_load_storage(
        RINGMesh::GeoModelBuilder3D &builder,
        PLineLoadingStorage &load_storage)
    {
        builder.geometry.set_line(load_storage.cur_line_, load_storage.vertices_);
        load_storage.cur_line_first_vertex_index_ +=
            static_cast<index_t>(load_storage.vertices_.size());
        load_storage.vertices_.clear();
        load_storage.cur_surf_polygon_corners_gocad_id_.clear();
        load_storage.cur_surf_polygon_ptr_.clear();
        load_storage.cur_surf_polygon_ptr_.push_back(0);
    }

    class LoadPLineNewGObj : public PLineLineParser
    {
    public:
        LoadPLineNewGObj(
            GeoModelBuilderPL &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : PLineLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            PLineLoadingStorage &load_storage)
        {
            if (!line.field_matches(1, "PLine"))
            {
                std::string gobj_type(line.field(1));
                throw SCARException("I/O",
                                    "While loading PLines, other object in file : " + gobj_type);
            }

            // Creation of a new contact ( <-> PLine )
            RINGMesh::gmge_id cur_contact =
                builder().geology.create_geological_entity(
                    RINGMesh::Contact3D::type_name_static());
            load_storage.cur_contact_ = cur_contact.index();
        }
    };

    class LoadPLineName : public PLineLineParser
    {
    public:
        LoadPLineName(
            GeoModelBuilderPL &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : PLineLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            PLineLoadingStorage &load_storage)
        {
            // Set the current contact name
            builder().info.set_geological_entity_name(
                RINGMesh::gmge_id(RINGMesh::Contact3D::type_name_static(),
                                  load_storage.cur_contact_),
                line.field(1));
        }
    };

    class LoadPLineLine : public PLineLineParser
    {
    public:
        LoadPLineLine(
            GeoModelBuilderPL &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : PLineLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            PLineLoadingStorage &load_storage)
        {
            scar_unused(line);
            if (!load_storage.vertices_.empty())
            {
                set_cur_line_mesh_from_load_storage(builder(), load_storage);
            }
            // Creation of new Line ( <-> ILine )
            RINGMesh::gmme_id cur_line = builder().topology.create_mesh_entity(
                RINGMesh::line_type_name_static());
            load_storage.cur_line_ = cur_line.index();

            // Set relationship between line and contact
            RINGMesh::gmge_id parent(RINGMesh::Contact3D::type_name_static(),
                                     load_storage.cur_contact_);
            builder().geology.add_parent_children_relation(parent, cur_line);
        }
    };

    class LoadPLineEndContact : public PLineLineParser
    {
    public:
        LoadPLineEndContact(
            GeoModelBuilderPL &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : PLineLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            PLineLoadingStorage &load_storage)
        {
            scar_unused(line);
            if (!load_storage.vertices_.empty())
            {
                set_cur_line_mesh_from_load_storage(builder(), load_storage);
            }
        }
    };
}

namespace SCAR
{
    PLineLineParser::PLineLineParser(
        GeoModelBuilderPL &gm_builder,
        RINGMesh::GeoModel3D &geomodel)
        : RINGMesh::GocadBaseParser(gm_builder, geomodel)
    {
    }

    PLineLoadingStorage::PLineLoadingStorage()
    {
        cur_line_first_vertex_index_ = GOCAD_OFFSET;
        cur_contact_ = 0;
        cur_line_ = 0;
    }

    void GeoModelBuilderPL::load_file()
    {
        // Read file and import information contained herein
        read_file();
        compute_corners();

        remove_duplicate_lines();
    }

    void GeoModelBuilderPL::read_line()
    {
        std::string keyword = file_line().field(0);
        std::unique_ptr<PLineLineParser> pline_parser(
            PLineLineParserFactory::create(keyword, *this, geomodel_));
        if (pline_parser)
        {
            pline_parser->execute(file_line(), pl_load_storage_);
        }
        else
        {
            std::unique_ptr<RINGMesh::GocadLineParser> gocad_parser =
                RINGMesh::GocadLineFactory::create(keyword, *this, geomodel_);
            if (gocad_parser)
            {
                gocad_parser->execute(file_line(), pl_load_storage_);
            }
        }
    }

    void GeoModelBuilderPL::compute_corners()
    {
        for (index_t l = 0; l < geomodel_.nb_lines(); ++l)
        {
            const RINGMesh::Line3D &line = geomodel_.line(l);
            // First vertex
            vec3 first_vertex = line.vertex(0);
            RINGMesh::gmme_id first_corner_id = topology.find_or_create_corner(
                first_vertex);
            // Last vertex
            vec3 last_vertex = line.vertex(line.nb_vertices() - 1);
            RINGMesh::gmme_id last_corner_id = topology.find_or_create_corner(
                last_vertex);
            // Update connectivity
            topology.add_line_corner_boundary_relation(line.index(),
                                                       first_corner_id.index());
            topology.add_line_corner_boundary_relation(line.index(),
                                                       last_corner_id.index());
        }
    }

    void GeoModelBuilderPL::save_geomodel(
        const RINGMesh::GeoModel2D &geomodel,
        const std::string &filename)
    {
        std::ofstream out(filename.c_str());
        out.precision(16);

        // Print Model3d headers
        out << "GOCAD PLine 1" << std::endl
            << "HEADER {" << std::endl
            << "name:"
            << geomodel.name() << std::endl
            << "}" << std::endl;

        out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl
            << "NAME Default"
            << std::endl
            << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
            << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl
            << "ZPOSITIVE Elevation"
            << std::endl
            << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl;

        index_t vertex_count{0};
        index_t edge_count{0};

        for (index_t line_id = 0; line_id < geomodel.nb_lines(); ++line_id)
        {
            ++edge_count;
            out << "ILINE" << std::endl;
            const RINGMesh::Line2D &cur_line = geomodel.line(line_id);
            for (index_t v = 0; v < cur_line.nb_vertices(); ++v)
            {
                out << "VRTX " << ++vertex_count << " " << cur_line.vertex(v).x
                    << " " << cur_line.vertex(v).y << " " << 0
                    << std::endl;
            }
            for (index_t v = 0; v < cur_line.nb_mesh_elements(); ++v)
            {
                out << "SEG " << edge_count;
                edge_count++;
                out << " " << edge_count << std::endl;
            }
        }
        out << "END" << std::endl;
    }

    void GeoModelBuilderPL::save_geomodel(
        const RINGMesh::GeoModel3D &geomodel,
        const std::string &filename)
    {
        std::ofstream out(filename.c_str());
        out.precision(16);

        // Print Model3d headers
        out << "GOCAD PLine 1" << std::endl
            << "HEADER {" << std::endl
            << "name:"
            << geomodel.name() << std::endl
            << "}" << std::endl;

        out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl
            << "NAME Default"
            << std::endl
            << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
            << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl
            << "ZPOSITIVE Elevation"
            << std::endl
            << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl;

        index_t vertex_count{0};
        index_t edge_count{0};

        for (index_t line_id = 0; line_id < geomodel.nb_lines(); ++line_id)
        {
            ++edge_count;
            out << "ILINE" << std::endl;
            const RINGMesh::Line3D &cur_line = geomodel.line(line_id);
            for (index_t v = 0; v < cur_line.nb_vertices(); ++v)
            {
                out << "VRTX " << ++vertex_count << " " << cur_line.vertex(v).x
                    << " " << cur_line.vertex(v).y << " " << cur_line.vertex(v).z
                    << std::endl;
            }
            for (index_t v = 0; v < cur_line.nb_mesh_elements(); ++v)
            {
                out << "SEG " << edge_count;
                edge_count++;
                out << " " << edge_count << std::endl;
            }
        }
        out << "END" << std::endl;
    }

    void GeoModelBuilderPL::remove_duplicate_lines()
    {
        std::set<RINGMesh::gmme_id> duplicate_lines_to_remove;
        find_duplicate_lines(duplicate_lines_to_remove);
        remove.remove_mesh_entities(duplicate_lines_to_remove);
    }

    void GeoModelBuilderPL::find_duplicate_lines(
        std::set<RINGMesh::gmme_id> &duplicate_lines_to_remove)
    {
        for (index_t line_id = 0; line_id < geomodel_.nb_lines(); ++line_id)
        {
            bool is_line_to_remove = is_line_duplicated(line_id);
            if (is_line_to_remove)
            {
                duplicate_lines_to_remove.insert(
                    RINGMesh::gmme_id(RINGMesh::Line3D::type_name_static(), line_id));
            }
        }
    }

    bool GeoModelBuilderPL::is_line_duplicated(const index_t line_id)
    {
        const RINGMesh::Line3D &query_line = geomodel_.line(line_id);
        for (index_t cur_line_id = 0; cur_line_id < line_id; ++cur_line_id)
        {
            const RINGMesh::Line3D &cur_line = geomodel_.line(cur_line_id);
            if (are_lines_identical(cur_line, query_line))
            {
                return true;
            }
        }
        return false;
    }

    void pline_import_factory_initialize()
    {
        PLineLineParserFactory::register_creator<LoadPLineNewGObj>("GOCAD");
        PLineLineParserFactory::register_creator<LoadPLineName>("name");
        PLineLineParserFactory::register_creator<LoadPLineLine>("ILINE");
        PLineLineParserFactory::register_creator<LoadPLineEndContact>("END");
    }

}
