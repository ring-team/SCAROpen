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

#include <ringmesh/basic/box.h>
#include <iostream>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/line_stream.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>

#include <ringmesh/io/io.h>

#include <scar/basic/command_line.h>

using namespace SCAR;

void hello()
{
    Logger::div("SCAR");
    Logger::out("", "Welcome to SCAR !");
    Logger::out("", "Sealed model Creation and Automatic Repair");
}

void import_arg_groups()
{
    CmdLine::import_arg_group("global");
    CmdLine::import_arg_group("in");
    CmdLine::import_arg_group("out");
    GEO::CmdLine::declare_arg("points", "",
                              "Path to file where polygon points are given");
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
        Logger::out("Step1", "Saving subtracted geomodel...");
        geomodel_save(output_geomodel, geomodel_out_name);
    }
}

template <index_t DIMENSION>
vecn<DIMENSION> points_barycenter(const std::vector<vecn<DIMENSION>> &points)
{
    vecn<DIMENSION> result;
    for (auto pt : points)
    {
        result += pt;
    }
    return result / points.size();
}

std::vector<index_t> cut_line_by_segment(
    RINGMesh::GeoModel2D &geomodel,
    index_t line_id,
    const RINGMesh::Geometry::Segment2D &segment)
{
    std::vector<index_t> new_lines;
    RINGMesh::GeoModelBuilder2D builder(geomodel);
    const RINGMesh::Line2D &line = geomodel.line(line_id);
    std::vector<vec2> cur_vertices;
    RINGMesh::gmme_id begin_corner = line.boundary_gmme(0);
    bool split_line = false;
    for (auto edge_i : RINGMesh::range(line.nb_mesh_elements()))
    {
        RINGMesh::Geometry::Segment2D edge(line.vertex(edge_i),
                                           line.vertex(edge_i + 1));
        bool intersection;
        vec2 point;
        std::tie(intersection, point) = RINGMesh::Intersection::segment_segment(
            segment, edge);
        cur_vertices.push_back(line.vertex(edge_i));
        if (!intersection)
        {
            continue;
        }
        else
        {
            split_line = true;
            cur_vertices.push_back(point);
            RINGMesh::gmme_id end_corner = builder.topology.find_or_create_corner(
                point);
            auto new_line_id = builder.topology.create_mesh_entity(
                RINGMesh::line_type_name_static());
            new_lines.push_back(new_line_id.index());
            builder.geometry.set_line(new_line_id.index(), cur_vertices);
            cur_vertices.clear();
            cur_vertices.push_back(point);
            builder.topology.add_line_corner_boundary_relation(new_line_id.index(),
                                                               begin_corner.index());
            builder.topology.add_line_corner_boundary_relation(new_line_id.index(),
                                                               end_corner.index());
            begin_corner = end_corner;
            continue;
        }
    }
    if (split_line)
    {
        cur_vertices.push_back(line.vertex(line.nb_vertices() - 1));
        RINGMesh::gmme_id end_corner = builder.topology.find_or_create_corner(
            line.vertex(line.nb_vertices() - 1));
        auto new_line_id = builder.topology.create_mesh_entity(
            RINGMesh::line_type_name_static());
        new_lines.push_back(new_line_id.index());
        builder.geometry.set_line(new_line_id.index(), cur_vertices);
        builder.topology.add_line_corner_boundary_relation(new_line_id.index(),
                                                           begin_corner.index());
        builder.topology.add_line_corner_boundary_relation(new_line_id.index(),
                                                           end_corner.index());
    }

    return new_lines;
}

void flag_outside_lines(
    const RINGMesh::GeoModel2D &geomodel,
    RINGMesh::Sign inside_sign,
    const RINGMesh::Geometry::Segment2D &segment,
    const std::vector<index_t> &new_line_ids,
    std::vector<index_t> &outside_lines)
{
    for (auto new_line_id : new_line_ids)
    {
        const auto &new_line = geomodel.line(new_line_id);
        RINGMesh::Sign sign0 = RINGMesh::Position::point_side_to_segment(
            new_line.vertex(0), segment);
        index_t vertex_id = 1;
        while (sign0 == RINGMesh::Sign::ZERO)
        {
            sign0 = RINGMesh::Position::point_side_to_segment(
                new_line.vertex(vertex_id++), segment);
        }
        if (sign0 == inside_sign)
        {
            continue;
        }
        else
        {
            outside_lines.push_back(new_line_id);
            continue;
        }
    }
}

void subtract_2D_geomodel(
    RINGMesh::GeoModel2D &geomodel,
    const std::vector<vec2> &polygon_points)
{
    std::vector<index_t> old_lines;
    auto nb_lines = geomodel.nb_lines();
    for (index_t i : RINGMesh::range(polygon_points.size()))
    {
        const vec2 &pt1 = polygon_points[i];
        const vec2 &pt2 = polygon_points[(i + 1) % polygon_points.size()];
        RINGMesh::Geometry::Segment2D cur_segment(pt1, pt2);
        for (auto line_id : RINGMesh::range(nb_lines))
        {
            std::vector<index_t> new_line_ids = cut_line_by_segment(geomodel,
                                                                    line_id, cur_segment);
            if (!new_line_ids.empty())
            {
                old_lines.push_back(line_id);
            }
        }
    }
    RINGMesh::GeoModelBuilder2D builder(geomodel);
    std::set<RINGMesh::gmme_id> lines_to_remove;
    for (auto i : old_lines)
    {
        lines_to_remove.insert(
            RINGMesh::gmme_id(RINGMesh::line_type_name_static(), i));
    }
    builder.remove.remove_mesh_entities(lines_to_remove);

    // @todo May be KO for non convex polygons
    vec2 polygon_barycenter = points_barycenter(polygon_points);
    std::set<RINGMesh::gmme_id> outside_lines_to_remove;
    for (auto line_id : RINGMesh::range(geomodel.nb_lines()))
    {
        bool inside = true;
        for (index_t i : RINGMesh::range(polygon_points.size()))
        {
            const vec2 &pt1 = polygon_points[i];
            const vec2 &pt2 = polygon_points[(i + 1) % polygon_points.size()];
            RINGMesh::Geometry::Segment2D cur_segment(pt1, pt2);
            RINGMesh::Sign inside_sign = RINGMesh::Position::point_side_to_segment(
                polygon_barycenter, cur_segment);

            const auto &line = geomodel.line(line_id);
            RINGMesh::Sign sign0 = RINGMesh::Position::point_side_to_segment(
                line.vertex(0), cur_segment);

            index_t vertex_id = 0;
            while (sign0 == RINGMesh::Sign::ZERO || RINGMesh::Position::point_inside_segment(
                                                        line.vertex(vertex_id), cur_segment))
            {
                vertex_id++;
                if (vertex_id >= line.nb_vertices())
                {
                    break;
                }
                sign0 = RINGMesh::Position::point_side_to_segment(
                    line.vertex(vertex_id), cur_segment);
            }
            if (sign0 == inside_sign)
            {
                continue;
            }
            else
            {
                inside = false;
                break;
            }
        }
        if (!inside)
        {
            outside_lines_to_remove.insert(
                RINGMesh::gmme_id(RINGMesh::line_type_name_static(), line_id));
        }
    }
    builder.remove.remove_mesh_entities(outside_lines_to_remove);

    std::set<RINGMesh::gmme_id> corners_to_remove;
    for (const auto &corner : geomodel.corners())
    {
        if (corner.nb_incident_entities() == 0)
        {
            corners_to_remove.insert(corner.gmme());
        }
    }
    builder.remove.remove_mesh_entities(corners_to_remove);
}

void load_and_subtract_2D_geomodel(
    const std::string &geomodel_in_name,
    const std::vector<vec2> &polygon_points)
{
    RINGMesh::GeoModel2D geomodel;
    geomodel_load(geomodel, geomodel_in_name);

    subtract_2D_geomodel(geomodel, polygon_points);
    save_output_model(geomodel);
}

std::vector<vec2> read_point_file(const std::string &filename)
{
    std::vector<vec2> points;
    GEO::LineInput line_input(filename);
    while (!line_input.eof() && line_input.get_line())
    {
        line_input.get_fields();
        if (line_input.nb_fields() == 0)
        {
            continue;
        }
        vec2 point;
        point.x = line_input.field_as_double(0);
        point.y = line_input.field_as_double(1);
        points.push_back(point);
    }

    return points;
}

int main(int argc, char **argv)
{
    try
    {

        hello();
        import_arg_groups();

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

        std::string geomodel_in_name{GEO::CmdLine::get_arg("in:geomodel")};
        if (geomodel_in_name.empty())
        {
            throw SCARException("I/O", "Give at least a filename in in:geomodel");
        }

        std::string point_file{GEO::CmdLine::get_arg("points")};
        if (point_file.empty())
        {
            throw SCARException("I/O", "Give at least a filename in points");
        }
        std::vector<vec2> polygon_points = read_point_file(point_file);

        load_and_subtract_2D_geomodel(geomodel_in_name, polygon_points);
    }
    catch (const SCARException &e)
    {
        Logger::err(e.category(), e.what());
        return 1;
    }
    catch (const std::exception &e)
    {
        Logger::err("Exception", e.what());
        return 1;
    }
    return 0;
}
