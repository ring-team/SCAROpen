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

#include <scar/tools/complexity_measures.h>

#include <geogram/basic/attributes.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/mesh_aabb.h>
#include <ringmesh/mesh/surface_mesh.h>

#ifdef SCAR_WITH_RINGGrid
#include <ringmesh/mesh/cartesian_grid.h>
#include <ringgrid/algorithm/rasterizer.h>
#include <ringgrid/basic/voxel_actions.h>
#endif

#include <scar/tools/utils.h>

namespace
{
    using namespace SCAR;

    const char SPACE = ' ';
    const char ENDL = '\n';
    const std::string BLANK_LINE = "\n\n";

    void compute_model_region_sizes(
        const RINGMesh::GeoModel2D &geomodel,
        double min_size,
        std::ofstream &out)
    {
        scar_unused(geomodel);
        scar_unused(min_size);
        scar_unused(out);
        return; // No region in 2D models
    }

    void compute_model_region_sizes(
        const RINGMesh::GeoModel3D &geomodel,
        double min_size,
        std::ofstream &out)
    {
        index_t nb_short_regions{0};
        double global_complexity{0};
        for (auto &region : geomodel.regions())
        {
            auto region_size = region.size();
            out << region.gmme() << SPACE << region_size << SPACE;
            if (region_size < min_size * min_size * min_size)
            {
                ++nb_short_regions;
                out << "COMPLEX";
            }
            out << ENDL;
        }
        out << "COMPLEX REGIONS: " << nb_short_regions << " ON "
            << geomodel.nb_regions() << ENDL;
        global_complexity = double(nb_short_regions) / double(geomodel.nb_regions());
        out << "GLOBAL COMPLEXITY: " << global_complexity << BLANK_LINE;
        return;
    }

    template <index_t DIMENSION>
    void compute_model_sizes(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        double min_size,
        std::ofstream &out)
    {
        out << "==== SIZES" << ENDL;
        out << "== LINES" << ENDL;
        index_t nb_short_lines{0};
        double global_complexity{0};
        for (auto &line : geomodel.lines())
        {
            auto line_size = line.size();
            out << line.gmme() << SPACE << line_size << SPACE;
            if (line_size < min_size)
            {
                ++nb_short_lines;
                out << "COMPLEX";
            }
            out << ENDL;
        }
        out << "COMPLEX LINES: " << nb_short_lines << " ON " << geomodel.nb_lines()
            << ENDL;
        global_complexity = double(nb_short_lines) / double(geomodel.nb_lines());
        out << "GLOBAL COMPLEXITY: " << global_complexity << BLANK_LINE;

        index_t nb_short_surfaces{0};
        global_complexity = 0;
        for (auto &surface : geomodel.surfaces())
        {
            auto surface_size = surface.size();
            out << surface.gmme() << SPACE << surface_size << SPACE;
            if (surface_size < min_size * min_size)
            {
                ++nb_short_surfaces;
                out << "COMPLEX";
            }
            out << ENDL;
        }
        out << "COMPLEX SURFACES: " << nb_short_surfaces << " ON "
            << geomodel.nb_surfaces() << ENDL;
        global_complexity = double(nb_short_surfaces) / double(geomodel.nb_surfaces());
        out << "GLOBAL COMPLEXITY: " << global_complexity << BLANK_LINE;

        compute_model_region_sizes(geomodel, min_size, out);
        return;
    }

    void compute_model_region_thicknesses(
        const RINGMesh::GeoModel2D &geomodel,
        double min_size,
        std::ofstream &out)
    {
        scar_unused(geomodel);
        scar_unused(min_size);
        scar_unused(out);
        return; // No region in 2D models
    }

    void compute_model_region_thicknesses(
        const RINGMesh::GeoModel3D &geomodel,
        double min_size,
        std::ofstream &out)
    {
        scar_unused(geomodel);
        scar_unused(min_size);
        scar_unused(out);
        return; // Not implemented in 3D
    }

    template <index_t DIMENSION>
    std::tuple<double, double> analyze_surface_thickness(
        const RINGMesh::Surface<DIMENSION> &surface,
        double min_size)
    {
        // Approximation on the number of points tested
        index_t nb_points_total{0};
        index_t nb_points_complex{0};
        double min_thickness = max_double();
        for (auto i : RINGMesh::range(surface.nb_boundaries()))
        {
            const auto &cur_line = surface.boundary(i);
            for (auto edge : RINGMesh::range(cur_line.nb_mesh_elements()))
            {
                auto point = cur_line.mesh_element_barycenter(edge);
                double min_distance = max_double();
                for (auto j : RINGMesh::range(surface.nb_boundaries()))
                {
                    if (i == j)
                    {
                        continue;
                    }
                    const auto &other_line = surface.boundary(j);
                    double distance;
                    std::tie(std::ignore, std::ignore, distance) =
                        other_line.edge_aabb().closest_edge(point);
                    min_distance = std::min(min_distance, distance);
                }
                // Analyze
                ++nb_points_total;
                if (min_distance < min_size)
                {
                    ++nb_points_complex;
                }
                min_thickness = std::min(min_thickness, min_distance);
            }
        }
        return std::make_tuple(
            double(nb_points_complex) / double(nb_points_total), min_thickness);
    }

    template <index_t DIMENSION>
    void compute_model_thicknesses(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        double min_size,
        std::ofstream &out)
    {
        out << "==== THICKNESSES" << ENDL;
        out << "== SURFACES" << ENDL;
        index_t nb_thin_surfaces{0};
        double global_complexity{0};
        for (auto &surface : geomodel.surfaces())
        {
            double thin_boundary_percent, min_thickness;
            std::tie(thin_boundary_percent, min_thickness) =
                analyze_surface_thickness(surface, min_size);
            out << surface.gmme() << SPACE << thin_boundary_percent << SPACE
                << " (min: " << min_thickness << ")" << SPACE;
            if (thin_boundary_percent > 0.5)
            {
                ++nb_thin_surfaces;
                out << "COMPLEX";
            }
            out << ENDL;
            global_complexity += thin_boundary_percent;
        }
        out << "COMPLEX SURFACES: " << nb_thin_surfaces << " ON "
            << geomodel.nb_surfaces() << ENDL;
        out << "GLOBAL COMPLEXITY: " << global_complexity << BLANK_LINE;

        compute_model_region_thicknesses(geomodel, min_size, out);
        return;
    }
    template <index_t DIMENSION>
    double get_surface_angle_at_corner(
        const RINGMesh::Surface<DIMENSION> &surface,
        index_t corner_id,
        index_t &line1,
        index_t &line2)
    {
        vecn<DIMENSION> corner;
        vecn<DIMENSION> point_on_line1;
        vecn<DIMENSION> point_on_line2;
        for (auto b_i : RINGMesh::range(surface.nb_boundaries()))
        {
            const auto &line = surface.boundary(b_i);
            for (auto b_b_i : RINGMesh::range(line.nb_boundaries()))
            {
                if (line.boundary(b_b_i).index() == corner_id)
                {
                    if (line1 == NO_ID)
                    {
                        line1 = line.index();
                        if (b_b_i == 0)
                        {
                            corner = line.vertex(0);
                            point_on_line1 = line.vertex(1);
                        }
                        else
                        {
                            auto nb_v = line.nb_vertices();
                            corner = line.vertex(nb_v - 1);
                            point_on_line1 = line.vertex(nb_v - 2);
                        }
                    }
                    else
                    {
                        line2 = line.index();
                        if (b_b_i == 0)
                        {
                            point_on_line2 = line.vertex(1);
                        }
                        else
                        {
                            auto nb_v = line.nb_vertices();
                            point_on_line2 = line.vertex(nb_v - 2);
                        }
                        break;
                    }
                }
            }
            if (line2 != NO_ID)
            {
                break;
            }
        }
        scar_assert(line1 != NO_ID && line2 != NO_ID);
        vecn<DIMENSION> e1 = normalize(point_on_line1 - corner);
        vecn<DIMENSION> e2 = normalize(point_on_line2 - corner);
        double angle = std::acos(dot(e1, e2)) * 180. / M_PI;
        return angle;
    }

    std::tuple<double, double> compute_model_angles_along_lines(
        const RINGMesh::GeoModel2D &geomodel,
        double min_size,
        std::ofstream &out)
    {
        scar_unused(geomodel);
        scar_unused(min_size);
        scar_unused(out);
        return std::make_tuple(0, 0); // No region in 2D models
    }

    std::tuple<double, double> compute_model_angles_along_lines(
        const RINGMesh::GeoModel3D &geomodel,
        double min_size,
        std::ofstream &out)
    {
        // TODO (iterate on line_segments)
        DEBUG("Compute angles along lines between surfaces, not implemented yet!");
        scar_unused(geomodel);
        scar_unused(min_size);
        scar_unused(out);
        return std::make_tuple(0, 0);
    }

    template <index_t DIMENSION>
    void compute_model_angles(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        double min_angle,
        std::ofstream &out)
    {
        out << "==== ANGLES" << ENDL;
        out << "== SURFACE BOUNDARIES at CORNERS" << ENDL;
        index_t nb_small_angles{0};
        double global_complexity{0};
        for (auto &surface : geomodel.surfaces())
        {
            index_t nb_small_angles_cur_surface{0};
            std::set<index_t> boundary_corners = surface_corners(geomodel,
                                                                 surface.index());
            out << surface.gmme() << ":" << ENDL;
            for (auto corner_id : boundary_corners)
            {
                index_t line1{NO_ID}, line2{NO_ID};
                double angle = get_surface_angle_at_corner(surface, corner_id,
                                                           line1, line2);
                out << "AT Corner " << corner_id << SPACE << "(Line " << line1
                    << " / Line " << line2 << ")" << SPACE << angle;
                if (angle < min_angle)
                {
                    ++nb_small_angles;
                    ++nb_small_angles_cur_surface;
                    out << SPACE << "COMPLEX";
                }
                out << ENDL;
            }
            global_complexity += double(nb_small_angles_cur_surface) / double(boundary_corners.size());
        }
        out << "COMPLEX SURFACE CORNERS: " << nb_small_angles << ENDL;
        out << "GLOBAL COMPLEXITY: " << global_complexity << BLANK_LINE;

        double percent_complex, minimum_angle;
        std::tie(percent_complex, minimum_angle) =
            compute_model_angles_along_lines(geomodel, min_angle, out);
    }

#ifdef SCAR_WITH_RINGGrid
    // Two versions of this same code due to RINGGrid code (non generic).
    void paint_complexity_values(
        const RINGMesh::GeoModel2D &geomodel,
        const RINGMesh::CartesianGrid2D &cartesian_grid,
        std::vector<index_t> &complexity_values)
    {
        RINGGrid::IncrementVoxelAction incrementer(complexity_values);
        RINGGrid::Rasterizer<2, RINGGrid::IncrementVoxelAction> model_rasterizer(
            cartesian_grid, incrementer);
        for (const auto &corner : geomodel.corners())
        {
            model_rasterizer.paint(corner.mesh());
        }
        for (const auto &line : geomodel.lines())
        {
            std::vector<bool> painted_cells(complexity_values.size(), false);
            RINGGrid::ValueVoxelAction<bool> act(painted_cells, true);
            RINGGrid::Rasterizer<2, RINGGrid::ValueVoxelAction<bool>> model_current_rasterizer(
                cartesian_grid, act);
            model_current_rasterizer.paint_precise(line.mesh());
            for (auto c : RINGMesh::range(complexity_values.size()))
            {
                if (painted_cells[c])
                {
                    ++complexity_values[c];
                }
            }
        }
        for (const auto &surface : geomodel.surfaces())
        {
            std::vector<bool> painted_cells(complexity_values.size(), false);
            RINGGrid::ValueVoxelAction<bool> act(painted_cells, true);
            RINGGrid::Rasterizer<2, RINGGrid::ValueVoxelAction<bool>> model_current_rasterizer(
                cartesian_grid, act);
            model_current_rasterizer.paint(surface.mesh());
            for (auto c : RINGMesh::range(complexity_values.size()))
            {
                if (painted_cells[c])
                {
                    ++complexity_values[c];
                }
            }
        }
    }

    void paint_complexity_values(
        const RINGMesh::GeoModel3D &geomodel,
        const RINGMesh::CartesianGrid3D &cartesian_grid,
        std::vector<index_t> &complexity_values)
    {
        RINGGrid::IncrementVoxelAction incrementer(complexity_values);
        RINGGrid::Rasterizer<3, RINGGrid::IncrementVoxelAction> model_rasterizer(
            cartesian_grid, incrementer);
        for (const auto &corner : geomodel.corners())
        {
            std::vector<RINGMesh::Geometry::Point3D> corner_point(1);
            corner_point[0] = RINGMesh::Geometry::Point3D(corner.vertex(0));
            model_rasterizer.paint(corner_point);
        }
        for (const auto &line : geomodel.lines())
        {
            std::vector<bool> painted_cells(complexity_values.size(), false);
            RINGGrid::ValueVoxelAction<bool> act(painted_cells, true);
            RINGGrid::Rasterizer<3, RINGGrid::ValueVoxelAction<bool>> model_current_rasterizer(
                cartesian_grid, act);
            std::vector<RINGMesh::Geometry::Segment3D> line_segments(
                line.nb_mesh_elements());
            for (auto edge_id : RINGMesh::range(line.nb_mesh_elements()))
            {
                line_segments[edge_id] = RINGMesh::Geometry::Segment3D(
                    line.mesh_element_vertex({edge_id, 0}),
                    line.mesh_element_vertex({edge_id, 1}));
            }
            model_current_rasterizer.paint(line_segments);
            for (auto c : RINGMesh::range(complexity_values.size()))
            {
                if (painted_cells[c])
                {
                    ++complexity_values[c];
                }
            }
        }
    }

    template <index_t DIMENSION>
    void compute_model_neighborhood_measure(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        double grid_resolution,
        std::ofstream &out)
    {
        // Compute Bounding Box of the model to compute
        // the number of Cartesian grid columns
        RINGMesh::Box<DIMENSION> bbox;
        for (const auto &line : geomodel.lines())
        {
            for (auto v : RINGMesh::range(line.nb_vertices()))
            {
                bbox.add_point(line.vertex(v));
            }
        }
        for (const auto &surface : geomodel.surfaces())
        {
            for (auto v : RINGMesh::range(surface.nb_vertices()))
            {
                bbox.add_point(surface.vertex(v));
            }
        }
        vecn<DIMENSION> bbox_dimension = bbox.diagonal();
        RINGMesh::ivecn<DIMENSION> grid_dimension;
        for (auto d : RINGMesh::range(DIMENSION))
        {
            grid_dimension[d] = static_cast<index_t>(std::ceil(
                                                         bbox_dimension[d] / grid_resolution) +
                                                     20);
        }
        // Frame aligned along X, Y ( and Z) axes.
        RINGMesh::ReferenceFrame<DIMENSION> frame;
        for (auto d : RINGMesh::range(DIMENSION))
        {
            vecn<DIMENSION> cur_dir;
            cur_dir[d] = grid_resolution;
            frame[d] = cur_dir;
        }
        vecn<DIMENSION> origin = bbox.center();
        for (auto d : RINGMesh::range(DIMENSION))
        {
            origin[d] -= grid_dimension[d] * 0.5 * frame[d][d];
        }
        frame.origin() = origin;
        // Create the Cartesian grid
        RINGMesh::CartesianGrid<DIMENSION> cartesian_grid(grid_dimension, frame);
        std::vector<index_t> complexity_values(cartesian_grid.nb_cells(), 0);
        paint_complexity_values(geomodel, cartesian_grid, complexity_values);

        // Output
        out << "GRID DIMENSIONS: ";
        for (auto d : RINGMesh::range(DIMENSION))
        {
            out << grid_dimension[d] << SPACE;
        }
        out << ENDL;
        for (auto v : complexity_values)
        {
            out << v << ENDL;
        }
    }
#endif
}

namespace SCAR
{

    template <index_t DIMENSION>
    void compute_model_complexity_global_measures(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        double min_size,
        double min_angle,
        const std::string &filename)
    {
        std::ofstream out;
        out.open(filename);
        compute_model_sizes(geomodel, min_size, out);
        compute_model_thicknesses(geomodel, min_size, out);
        compute_model_angles(geomodel, min_angle, out);
        //TODO: see if we need to compute model shapes (with projection) [Pellerin2015]
        out << std::flush;
        out.close();
    }

    template <index_t DIMENSION>
    void scar_api compute_model_complexity_local_measures(
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        double grid_resolution,
        const std::string &filename)
    {
#ifdef SCAR_WITH_RINGGrid
        std::ofstream out;
        out.open(filename);
        compute_model_neighborhood_measure(geomodel, grid_resolution, out);
        out << std::flush;
        out.close();
#else
        scar_unused(geomodel);
        scar_unused(grid_resolution);
        scar_unused(filename);
        throw SCARException("RINGGrid",
                            "You need RINGGrid for running local complexity measures.");
#endif
    }

    //Explicite

    template void scar_api compute_model_complexity_global_measures(
        const RINGMesh::GeoModel2D &,
        double,
        double,
        const std::string &);

    template void scar_api compute_model_complexity_local_measures(
        const RINGMesh::GeoModel2D &geomodel,
        double grid_resolution,
        const std::string &filename);

    template void scar_api compute_model_complexity_global_measures(
        const RINGMesh::GeoModel3D &,
        double,
        double,
        const std::string &);

    template void scar_api compute_model_complexity_local_measures(
        const RINGMesh::GeoModel3D &geomodel,
        double grid_resolution,
        const std::string &filename);
}
