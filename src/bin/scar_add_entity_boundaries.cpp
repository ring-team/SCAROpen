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

#include <iostream>

#include <geogram/basic/command_line.h>
#include <geogram/basic/line_stream.h>

#include <ringmesh/basic/box.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

#include <scar/basic/command_line.h>

using namespace SCAR;

void import_arg_group_add_info()
{
    GEO::CmdLine::declare_arg_group("mesh_entity",
                                    "Add a mesh entity with its boundaries");
    GEO::CmdLine::declare_arg("mesh_entity:type", "Surface",
                              "The kind of mesh entity to add");
    GEO::CmdLine::declare_arg("mesh_entity:file", "",
                              "The path to the file giving all boundaries for each new mesh entity");
    GEO::CmdLine::declare_arg("clear", false,
                              "Clear all entities of mesh_entity:type before add new ones");
}

void add_boundaries(
    RINGMesh::GeoModel2D &geomodel,
    const RINGMesh::gmme_id &me_id,
    GEO::LineInput &file_reader)
{
    RINGMesh::GeoModelBuilder2D geomodel_builder(geomodel);
    for (auto i : RINGMesh::range(file_reader.nb_fields()))
    {
        if (me_id.type() == RINGMesh::line_type_name_static())
        {
            geomodel_builder.topology.add_line_corner_boundary_relation(
                me_id.index(), file_reader.field_as_uint(i));
        }
        else if (me_id.type() == RINGMesh::surface_type_name_static())
        {
            geomodel_builder.topology.add_surface_line_boundary_relation(
                me_id.index(), file_reader.field_as_uint(i), true);
        }
        else
        {
            scar_assert_not_reached;
        }
    }
}

void add_boundaries(
    RINGMesh::GeoModel3D &geomodel,
    const RINGMesh::gmme_id &me_id,
    GEO::LineInput &file_reader)
{
    RINGMesh::GeoModelBuilder3D geomodel_builder(geomodel);
    for (auto i : RINGMesh::range(file_reader.nb_fields()))
    {
        if (me_id.type() == RINGMesh::line_type_name_static())
        {
            geomodel_builder.topology.add_line_corner_boundary_relation(
                me_id.index(), file_reader.field_as_uint(i));
        }
        else if (me_id.type() == RINGMesh::surface_type_name_static())
        {
            geomodel_builder.topology.add_surface_line_boundary_relation(
                me_id.index(), file_reader.field_as_uint(i));
        }
        else if (me_id.type() == RINGMesh::region_type_name_static())
        {
            geomodel_builder.topology.add_region_surface_boundary_relation(
                me_id.index(), file_reader.field_as_uint(i), true);
        }
        else
        {
            scar_assert_not_reached;
        }
    }
}

template <index_t DIMENSION>
void add_mesh_entity_boundaries(RINGMesh::GeoModel<DIMENSION> &geomodel)

{
    RINGMesh::GeoModelBuilder<DIMENSION> geomodel_builder(geomodel);
    auto mesh_entity_type = GEO::CmdLine::get_arg("mesh_entity:type");
    if (mesh_entity_type.empty())
    {
        throw SCARException("I/O",
                            "Give type of mesh entities to add (mesh_entity:type)");
    }

    std::string filepath_boundaries = GEO::CmdLine::get_arg("mesh_entity:file");
    if (filepath_boundaries.empty())
    {
        throw SCARException("I/O",
                            "Give info for new entities boundaries (mesh_entity:file)");
    }
    Logger::out("MeshEntity", "Using file for boundaries: ", filepath_boundaries);
    std::ifstream input(filepath_boundaries.c_str());
    if (!input)
    {
        throw SCARException("I/O", "Failed to open file " + filepath_boundaries);
    }
    GEO::LineInput file_reader(filepath_boundaries);
    scar_assert(file_reader.OK());

    bool do_clear = GEO::CmdLine::get_arg_bool("clear");
    if (do_clear)
    {

        std::set<RINGMesh::gmme_id> to_delete;
        for (auto i : RINGMesh::range(
                 geomodel.nb_mesh_entities(
                     RINGMesh::MeshEntityType(mesh_entity_type))))
        {
            to_delete.emplace(RINGMesh::MeshEntityType(mesh_entity_type), i);
        }
        geomodel_builder.remove.remove_mesh_entities(to_delete);
        geomodel.mesh.vertices.clear();
    }

    while (!file_reader.eof() && file_reader.get_line())
    {
        // Create a new entity and add its boundaries
        file_reader.get_fields();

        auto me_id = geomodel_builder.topology.create_mesh_entity(
            RINGMesh::MeshEntityType(mesh_entity_type));
        add_boundaries(geomodel, me_id, file_reader);
    }
    geomodel.mesh.vertices.clear();
}

template <index_t DIMENSION>
void save_output_model(const RINGMesh::GeoModel<DIMENSION> &output_geomodel)
{
    std::string geomodel_out_name{GEO::CmdLine::get_arg("out:geomodel")};
    if (geomodel_out_name.empty())
    {
        throw SCARException("I/O", "No filename in out:geomodel");
    }
    else
    {
        Logger::out("Step1", "Saving updated geomodel...");
        geomodel_save(output_geomodel, geomodel_out_name);
    }
}

template <index_t DIMENSION>
void load_geomodel_and_add_mesh_entity_boundaries(const std::string &filename)
{
    RINGMesh::GeoModel<DIMENSION> geomodel;
    geomodel_load(geomodel, filename);

    add_mesh_entity_boundaries(geomodel);
    RINGMesh::is_geomodel_valid(geomodel);
    save_output_model(geomodel);
}

int main(int argc, char **argv)
{
    try
    {
        import_arg_group_add_info();

        Logger::div("SCAR");
        Logger::out("", "Sealed model Creation and Automatic Repair");
        Logger::div("SCAR - Add information to the GeoModel");
        Logger::out("", "Welcome! In this executable, ",
                    "you can information to your geological model ");

        CmdLine::import_arg_group("in");
        CmdLine::import_arg_group("out");

        if (argc == 1)
        {
            GEO::CmdLine::show_usage();
            return 0;
        }

        std::vector<std::string> filenames;
        if (!GEO::CmdLine::parse(argc, argv, filenames))
        {
            return 1;
        }

        std::string geomodel_in_name = GEO::CmdLine::get_arg("in:geomodel");
        if (geomodel_in_name.empty())
        {
            throw SCARException("I/O", "Give at least a filename in in:geomodel");
        }

        std::string geomodel_out_name = GEO::CmdLine::get_arg("out:geomodel");
        if (geomodel_out_name.empty())
        {
            throw SCARException("I/O", "Give at least a filename in out:geomodel");
        }

        index_t dimension = RINGMesh::find_geomodel_dimension(geomodel_in_name);

        if (dimension == 2)
        {
            load_geomodel_and_add_mesh_entity_boundaries<2>(geomodel_in_name);
        }
        else if (dimension == 3)
        {
            load_geomodel_and_add_mesh_entity_boundaries<3>(geomodel_in_name);
        }
        else
        {
            throw SCARException("I/O", "Can not deal with GeoModel", dimension,
                                "D.");
        }
    }
    catch (const RINGMesh::RINGMeshException &e)
    {
        GEO::Logger::err(e.category()) << e.what() << std::endl;
        return 1;
    }
    catch (const std::exception &e)
    {
        GEO::Logger::err("Exception") << e.what() << std::endl;
        return 1;
    }
    return 0;
}
