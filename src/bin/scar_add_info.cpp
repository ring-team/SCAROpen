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
#include <ringmesh/geomodel/builder/geomodel_builder_geology.h>
#include <ringmesh/geomodel/builder/geomodel_builder_remove.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/io/io.h>

#include <scar/basic/command_line.h>

using namespace SCAR;

void import_arg_group_add_info()
{
    GEO::CmdLine::declare_arg_group("geol_type",
                                    "Add info about geological type of entities");
    GEO::CmdLine::declare_arg("geol_type:file", "",
                              "The path to the file giving all elements for each new geological entity");
    GEO::CmdLine::declare_arg("geol_type:type", "none",
                              "The geological type to add");
    GEO::CmdLine::declare_arg("geol_type:name", "NewGeolEntity",
                              "The name of added geological entity");
    GEO::CmdLine::declare_arg("geol_type:entity_type",
                              RINGMesh::line_type_name_static().string(),
                              "Type of the entity(ies) on which add a geol_type");
    GEO::CmdLine::declare_arg("geol_type:entity_id", double(NO_ID),
                              "Index of the entity on which add a geol_type");
    GEO::CmdLine::declare_arg("geol_type:entity_ids_begin", double(NO_ID),
                              "First index of range of indices of the entities on which add a geol_type");
    GEO::CmdLine::declare_arg("geol_type:entity_ids_end", double(NO_ID),
                              "Last index of range of indices of the entities on which add a geol_type");
    GEO::CmdLine::declare_arg("geol_type:entity_ids", "",
                              "Unordered set of indices of the entities on which add a geol_type");
}

template <index_t DIMENSION>
RINGMesh::GeologicalEntityType determine_parent_type(
    const RINGMesh::MeshEntityType &mesh_entity_type)
{
    if (RINGMesh::Contact<DIMENSION>::child_type_name_static() == mesh_entity_type)
    {
        return RINGMesh::Contact<DIMENSION>::type_name_static();
    }
    else if (RINGMesh::Interface<DIMENSION>::child_type_name_static() == mesh_entity_type)
    {
        return RINGMesh::Interface<DIMENSION>::type_name_static();
    }
    else if (RINGMesh::Layer<DIMENSION>::child_type_name_static() == mesh_entity_type)
    {
        return RINGMesh::Layer<DIMENSION>::type_name_static();
    }
    else
    {
        scar_assert_not_reached;
        return RINGMesh::GeologicalEntityType();
    }
}

RINGMesh::MeshEntityType determine_mesh_entity_type(
    const std::string &mesh_entity_type_string)
{
    if (RINGMesh::corner_type_name_static().string() == mesh_entity_type_string)
    {
        return RINGMesh::corner_type_name_static();
    }
    else if (RINGMesh::line_type_name_static().string() == mesh_entity_type_string)
    {
        return RINGMesh::line_type_name_static();
    }
    else if (RINGMesh::surface_type_name_static().string() == mesh_entity_type_string)
    {
        return RINGMesh::surface_type_name_static();
    }
    else if (RINGMesh::region_type_name_static().string() == mesh_entity_type_string)
    {
        return RINGMesh::region_type_name_static();
    }
    else
    {
        scar_assert_not_reached;
        return RINGMesh::corner_type_name_static();
    }
}

template <index_t DIMENSION>
void add_info(
    RINGMesh::GeoModel<DIMENSION> &geomodel,
    const std::string &geol_type,
    const std::string &geol_name,
    const RINGMesh::MeshEntityType &mesh_entity_type,
    const std::vector<index_t> &mesh_entity_ids)
{
    RINGMesh::GeoModelBuilder<DIMENSION> geomodel_builder(geomodel);
    auto geol_feature =
        RINGMesh::GeoModelGeologicalEntity<DIMENSION>::determine_geological_type(
            geol_type);

    RINGMesh::gmge_id new_geol_entity_id =
        geomodel_builder.geology.create_geological_entity(
            determine_parent_type<DIMENSION>(mesh_entity_type));
    geomodel_builder.info.set_geological_entity_name(new_geol_entity_id,
                                                     geol_name);
    geomodel_builder.geology.set_geological_entity_geol_feature(new_geol_entity_id,
                                                                geol_feature);
    for (auto id : mesh_entity_ids)
    {
        const auto &cur_mesh_entity = geomodel.mesh_entity(mesh_entity_type, id);
        if (!cur_mesh_entity.has_parent(new_geol_entity_id.type()))
        {
            geomodel_builder.geology.add_parent_children_relation(
                new_geol_entity_id, cur_mesh_entity.gmme());
        }
        else
        {
            geomodel_builder.geology.set_mesh_entity_parent(cur_mesh_entity.gmme(),
                                                            0, new_geol_entity_id, true);
        }
    }

    std::set<RINGMesh::gmge_id> empty_interfaces;
    for (auto &geol_entity : geomodel.geol_entities(
             RINGMesh::Interface<DIMENSION>::type_name_static()))
    {
        if (geol_entity.nb_children() == 0)
        {
            empty_interfaces.insert(geol_entity.gmge());
        }
    }
    geomodel_builder.remove.remove_geological_entities(empty_interfaces);
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

std::vector<index_t> determine_mesh_entity_ids()
{
    std::vector<index_t> results;
    index_t unique_id = GEO::String::to_uint(
        GEO::CmdLine::get_arg("geol_type:entity_id"));
    if (unique_id != NO_ID)
    {
        results.push_back(unique_id);
        return results;
    }

    index_t first_id = GEO::String::to_uint(
        GEO::CmdLine::get_arg("geol_type:entity_ids_begin"));
    if (first_id != NO_ID)
    {
        index_t last_id = GEO::String::to_uint(
            GEO::CmdLine::get_arg("geol_type:entity_ids_end"));
        if (first_id == NO_ID)
        {
            throw SCARException("I/O",
                                "Give at least one index in geol_type:entity_ids_begin AND "
                                "one index in geol_type:entity_ids_end");
        }
        for (auto i : RINGMesh::range(first_id, last_id + 1))
        {
            results.push_back(i);
        }
        return results;
    }

    std::string set_ids = GEO::CmdLine::get_arg("geol_type:entity_ids");
    if (set_ids.empty())
    {
        throw SCARException("I/O",
                            "Give at least an index in geol_type:entity_id or "
                            "geol_type:entity_ids or by using range");
    }
    std::vector<std::string> split_ids;
    GEO::String::split_string(set_ids, ' ', split_ids, true);
    for (auto split_id_string : split_ids)
    {
        auto id = GEO::String::to_uint(split_id_string);
        results.push_back(id);
    }
    return results;
}

template <index_t DIMENSION>
void load_geomodel_and_add_info(
    const std::string &filename,
    const std::string &geol_type,
    const std::string &geol_name,
    const RINGMesh::MeshEntityType &mesh_entity_type,
    const std::vector<index_t> &mesh_entity_ids)
{
    RINGMesh::GeoModel<DIMENSION> geomodel;
    geomodel_load(geomodel, filename);

    add_info(geomodel, geol_type, geol_name, mesh_entity_type, mesh_entity_ids);

    save_output_model(geomodel);
}

template <index_t DIMENSION>
void load_geomodel_and_add_infos_from_file(
    const std::string &filename_geomodel,
    const std::string &filename_entity)
{
    RINGMesh::GeoModel<DIMENSION> geomodel;
    geomodel_load(geomodel, filename_geomodel);

    GEO::LineInput file_reader(filename_entity);
    scar_assert(file_reader.OK());
    while (!file_reader.eof() && file_reader.get_line())
    {
        // Read each line of the file to obtain the info to had
        // (geol_type, geol_name, mesh_entity_type, mesh_entity_ids)
        file_reader.get_fields();
        std::string geol_type = file_reader.field(0);
        std::string geol_name = file_reader.field(1);
        auto mesh_entity_type = determine_mesh_entity_type(
            file_reader.field(2));
        std::vector<index_t> mesh_entity_ids;
        for (auto i : RINGMesh::range(3, file_reader.nb_fields()))
        {
            mesh_entity_ids.push_back(file_reader.field_as_uint(i));
        }
        add_info(geomodel, geol_type, geol_name, mesh_entity_type, mesh_entity_ids);
    }

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

        std::string entity_file = GEO::CmdLine::get_arg("geol_type:file");
        std::string geol_type = GEO::CmdLine::get_arg("geol_type:type");
        std::string geol_name = GEO::CmdLine::get_arg("geol_type:name");
        std::string me_type_string = GEO::CmdLine::get_arg(
            "geol_type:entity_type");
        if (entity_file.empty())
        {
            if (geol_type.empty())
            {
                throw SCARException("I/O", "Give an info file or at least a geol entity type in geol_type:type");
            }
            if (geol_name.empty())
            {
                throw SCARException("I/O",
                                    "Give an info file or at least an entity name in geol_type:name");
            }
            if (me_type_string.empty())
            {
                throw SCARException("I/O",
                                    "Give an info file or at least a mesh entity type in geol_type:entity_type");
            }

            auto mesh_entity_type = determine_mesh_entity_type(
                me_type_string);

            auto mesh_entity_ids = determine_mesh_entity_ids();

            index_t dimension = RINGMesh::find_geomodel_dimension(geomodel_in_name);

            if (dimension == 2)
            {
                load_geomodel_and_add_info<2>(geomodel_in_name, geol_type, geol_name,
                                              mesh_entity_type, mesh_entity_ids);
            }
            else if (dimension == 3)
            {
                load_geomodel_and_add_info<3>(geomodel_in_name, geol_type, geol_name,
                                              mesh_entity_type, mesh_entity_ids);
            }
            else
            {
                throw SCARException("I/O", "Can not deal with GeoModel", dimension,
                                    "D.");
            }
        }
        else
        {

            const index_t dimension = RINGMesh::find_geomodel_dimension(geomodel_in_name);
            if (dimension == 2)
            {
                load_geomodel_and_add_infos_from_file<2>(geomodel_in_name, entity_file);
            }
            else if (dimension == 3)
            {
                load_geomodel_and_add_infos_from_file<3>(geomodel_in_name, entity_file);
            }
            else
            {
                throw SCARException("I/O", "Can not deal with GeoModel", dimension,
                                    "D.");
            }
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
