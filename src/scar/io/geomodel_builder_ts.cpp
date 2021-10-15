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

#include <scar/io/geomodel_builder_ts.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/builder/geomodel_builder_topology.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/mesh_index.h>

namespace
{

    using namespace SCAR;

    // Indices begin to 1 in Gocad
    index_t GOCAD_OFFSET = 1;

    void set_cur_surface_geometry_from_load_storage(
        RINGMesh::GeoModelBuilder3D &builder,
        TSurfLoadingStorage &load_storage)
    {
        builder.geometry.set_surface_geometry(load_storage.cur_surface_,
                                              load_storage.vertices_, load_storage.cur_surf_polygon_corners_gocad_id_,
                                              load_storage.cur_surf_polygon_ptr_);
        load_storage.cur_surf_first_vertex_index_ +=
            static_cast<index_t>(load_storage.vertices_.size());
        load_storage.vertices_.clear();
        load_storage.cur_surf_polygon_corners_gocad_id_.clear();
        load_storage.cur_surf_polygon_ptr_.clear();
        load_storage.cur_surf_polygon_ptr_.push_back(0);
    }

    void get_gme_vertex_in_interface(
        const RINGMesh::GeoModel3D &geomodel,
        const index_t index_in_interface,
        const TSurfLoadingStorage &load_storage,
        index_t &surface_id,
        index_t &id_in_surface)
    {
        index_t prev_nb_interface_vertices = 0;
        const RINGMesh::GeoModelGeologicalEntity3D &cur_interface =
            geomodel.geological_entity(RINGMesh::Interface3D::type_name_static(),
                                       load_storage.cur_interface_);

        for (index_t s : RINGMesh::range(cur_interface.nb_children()))
        {
            index_t cur_nb_vertices = cur_interface.child(s).nb_vertices();
            if (prev_nb_interface_vertices + cur_nb_vertices > index_in_interface)
            {
                surface_id = cur_interface.child_gmme(s).index();
                id_in_surface = index_in_interface - prev_nb_interface_vertices;
                return;
            }
            else
            {
                prev_nb_interface_vertices += cur_nb_vertices;
            }
        }
        scar_assert_not_reached;
        return;
    }

    /*!
     * @brief Check if a point is a corner of the geomodel
     * @todo use an aabb tree (optimization)
     */
    bool is_corner(const RINGMesh::GeoModel3D &geomodel, const vec3 &point)
    {
        for (index_t i = 0; i < geomodel.nb_corners(); ++i)
        {
            if (geomodel.corner(i).vertex(0) == point)
            {
                return true;
            }
        }
        return false;
    }
    class LoadTSurfNewGObj : public TSurfLineParser
    {
    public:
        LoadTSurfNewGObj(
            GeoModelBuilderTS &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : TSurfLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            TSurfLoadingStorage &load_storage)
        {
            if (!line.field_matches(1, "TSurf"))
            {
                std::string gobj_type(line.field(1));
                throw SCARException("I/O",
                                    "While loading TSurfs, other object in file : " + gobj_type);
            }

            // Creation of a new interface ( <-> TSurf )
            RINGMesh::gmge_id cur_interface =
                builder().geology.create_geological_entity(
                    RINGMesh::Interface3D::type_name_static());
            load_storage.cur_interface_ = cur_interface.index();
        }
    };

    class LoadTSurfName : public TSurfLineParser
    {
    public:
        LoadTSurfName(
            GeoModelBuilderTS &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : TSurfLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            TSurfLoadingStorage &load_storage)
        {
            // Set the current interface name
            builder().info.set_geological_entity_name(
                RINGMesh::gmge_id(RINGMesh::Interface3D::type_name_static(),
                                  load_storage.cur_interface_),
                line.field(1));
        }
    };

    class LoadTSurfAtom final : public TSurfLineParser
    {
    public:
        LoadTSurfAtom(
            GeoModelBuilderTS &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : TSurfLineParser(gm_builder, geomodel)
        {
        }

    private:
        void execute(GEO::LineInput &line, TSurfLoadingStorage &load_storage) final
        {
            index_t global_vertex_id = line.field_as_uint(2) - GOCAD_OFFSET;
            // Find surface and local index in surface
            index_t surface_id{NO_ID};
            index_t vertex_id_in_surface{NO_ID};
            index_t vertex_count{0};
            for (auto surf_id : RINGMesh::range(load_storage.cur_surface_))
            {
                index_t vertex_count_save{vertex_count};
                vertex_count += geomodel().surface(surf_id).nb_vertices();
                if (vertex_count > global_vertex_id)
                {
                    surface_id = surf_id;
                    vertex_id_in_surface = global_vertex_id - vertex_count_save;
                    break;
                }
            }
            scar_assert(surface_id != NO_ID && vertex_id_in_surface != NO_ID);
            const vec3 &vertex = geomodel().surface(surface_id).vertex(vertex_id_in_surface);
            load_storage.vertices_.push_back(vertex);
        }
    };

    class LoadTSurfTriangle : public TSurfLineParser
    {
    public:
        LoadTSurfTriangle(
            GeoModelBuilderTS &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : TSurfLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            TSurfLoadingStorage &load_storage)
        {
            load_storage.cur_surf_polygon_corners_gocad_id_.push_back(
                line.field_as_uint(1) - load_storage.cur_surf_first_vertex_index_);
            load_storage.cur_surf_polygon_corners_gocad_id_.push_back(
                line.field_as_uint(2) - load_storage.cur_surf_first_vertex_index_);
            load_storage.cur_surf_polygon_corners_gocad_id_.push_back(
                line.field_as_uint(3) - load_storage.cur_surf_first_vertex_index_);
            load_storage.end_polygon();
        }
    };

    class LoadTSurfSurface : public TSurfLineParser
    {
    public:
        LoadTSurfSurface(
            GeoModelBuilderTS &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : TSurfLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            TSurfLoadingStorage &load_storage)
        {
            scar_unused(line);
            if (!load_storage.vertices_.empty())
            {
                set_cur_surface_geometry_from_load_storage(builder(),
                                                           load_storage);
            }

            // Creation of new Surface ( <-> TFace )
            RINGMesh::gmme_id cur_surface = builder().topology.create_mesh_entity(
                RINGMesh::Surface3D::type_name_static());
            load_storage.cur_surface_ = cur_surface.index();

            // Set relationship between surface and interface
            RINGMesh::gmge_id parent(RINGMesh::Interface3D::type_name_static(),
                                     load_storage.cur_interface_);
            builder().geology.add_parent_children_relation(parent, cur_surface);
        }
    };

    class LoadTSurfEndInterface : public TSurfLineParser
    {
    public:
        LoadTSurfEndInterface(
            GeoModelBuilderTS &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : TSurfLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            TSurfLoadingStorage &load_storage)
        {
            scar_unused(line);
            if (!load_storage.vertices_.empty())
            {
                set_cur_surface_geometry_from_load_storage(builder(),
                                                           load_storage);
                load_storage.cur_surf_first_vertex_index_ = GOCAD_OFFSET;
            }
        }
    };

    class LoadTSurfCorner : public TSurfLineParser
    {
    public:
        LoadTSurfCorner(
            GeoModelBuilderTS &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : TSurfLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            TSurfLoadingStorage &load_storage)
        {
            if (!load_storage.vertices_.empty())
            {
                set_cur_surface_geometry_from_load_storage(builder(),
                                                           load_storage);
                load_storage.cur_surf_first_vertex_index_ = GOCAD_OFFSET;
            }

            index_t corner_id_in_interface = line.field_as_uint(1) - GOCAD_OFFSET;
            vec3 position = get_corner_position(corner_id_in_interface,
                                                load_storage);
            if (!is_corner(geomodel(), position))
            {
                RINGMesh::gmme_id new_corner = builder().topology.create_mesh_entity(
                    RINGMesh::Corner3D::type_name_static());

                builder().geometry.set_corner(new_corner.index(), position);
            }
        }

        vec3 get_corner_position(
            const index_t index_in_interface,
            const TSurfLoadingStorage &load_storage)
        {
            index_t surface_id;
            index_t id_in_surface;
            get_gme_vertex_in_interface(geomodel(), index_in_interface,
                                        load_storage, surface_id, id_in_surface);

            vec3 corner_position = NO_POINT_3D;
            corner_position = geomodel().surface(surface_id).vertex(id_in_surface);

            scar_assert(corner_position != NO_POINT_3D);
            return corner_position;
        }
    };

    class LoadTSurfLine : public TSurfLineParser
    {
    public:
        LoadTSurfLine(
            GeoModelBuilderTS &gm_builder,
            RINGMesh::GeoModel3D &geomodel)
            : TSurfLineParser(gm_builder, geomodel)
        {
        }

    private:
        virtual void execute(
            GEO::LineInput &line,
            TSurfLoadingStorage &load_storage)
        {
            RINGMesh::gmme_id new_line = builder().topology.create_mesh_entity(
                RINGMesh::Line3D::type_name_static());

            index_t line_first_vertex = line.field_as_uint(2) - GOCAD_OFFSET;
            index_t line_second_vertex = line.field_as_uint(3) - GOCAD_OFFSET;

            // Find first corner (corresponding to line_first_vertex)
            // + find involved surface
            index_t incident_surface;
            index_t line_first_vertex_surf_index;
            get_gme_vertex_in_interface(geomodel(), line_first_vertex, load_storage,
                                        incident_surface, line_first_vertex_surf_index);

            // Propagation toward direction given by second vertex
            index_t incident_surface_compare;
            index_t line_second_vertex_surf_index;
            get_gme_vertex_in_interface(geomodel(), line_second_vertex,
                                        load_storage, incident_surface_compare,
                                        line_second_vertex_surf_index);
            scar_assert(incident_surface == incident_surface_compare);

            // Propagation to find line vertices
            std::vector<vec3> line_vertices;
            find_line_vertices_by_propagation(geomodel(), incident_surface,
                                              line_first_vertex_surf_index, line_second_vertex_surf_index,
                                              line_vertices);

            // Set line geometry
            builder().geometry.set_line(new_line.index(), line_vertices);

            // Topology information
            /// Line - Surface relationship
            RINGMesh::gmme_id incident_surface_id(
                RINGMesh::Surface3D::type_name_static(), incident_surface);
            builder().topology.add_surface_line_boundary_relation(
                incident_surface_id.index(), new_line.index());
            /// Line - Corner relationship
            index_t first_corner_id = find_corner(geomodel(),
                                                  line_vertices.front());
            index_t last_corner_id = find_corner(geomodel(), line_vertices.back());
            builder().topology.add_line_corner_boundary_relation(new_line.index(),
                                                                 first_corner_id);
            builder().topology.add_line_corner_boundary_relation(new_line.index(),
                                                                 last_corner_id);
        }

        /*!
         * @brief Return true if the propagation direction is toward next vertex,
         * false if is toward previous vertex
         * param[in] surface   Surface to which vertices belong
         * param[in] first_vertex_id Index in the surface of the first point of
         * the line
         * param[in] second_vertex_id Index in the surface of the second point of
         * the line
         */
        bool determine_propagation_direction(
            const RINGMesh::Surface3D &surface,
            const index_t &first_vertex_id,
            const index_t &second_vertex_id)
        {
            const RINGMesh::SurfaceMesh3D &surface_mesh = surface.mesh();
            index_t polygon_id = surface_mesh.polygon_from_vertex_ids(
                first_vertex_id, second_vertex_id);
            index_t first_vertex_local_index = surface_mesh.vertex_index_in_polygon(
                polygon_id, first_vertex_id);
            RINGMesh::ElementLocalVertex first_vertex_polygon_id(polygon_id,
                                                                 first_vertex_local_index);

            RINGMesh::ElementLocalVertex next_polygon_vertex_id =
                surface_mesh.next_polygon_vertex(first_vertex_polygon_id);
            RINGMesh::ElementLocalVertex prev_polygon_vertex_id =
                surface_mesh.prev_polygon_vertex(first_vertex_polygon_id);
            if (surface_mesh.polygon_vertex(next_polygon_vertex_id) == second_vertex_id)
            {
                return true;
            }
            else if (surface_mesh.polygon_vertex(prev_polygon_vertex_id) == second_vertex_id)
            {
                return false;
            }
            else
            {
                // Problem
                scar_assert_not_reached;
                return false;
            }
        }

        /*!
     * @brief Computes the line vertices
     * @details By propagation on the \p incident_surface border, this function
     * finds the vertices until it reach a corner of the /p geomodel
     * param[in] geomodel           Involved GeoModel
     * param[in] incident_surface   Surface on which line vertices are
     * param[in] first_vertex_id    Index in the surface of the first point of
     * the line
     * param[in] second_vertex_id   Index in the surface of the second point of
     * the line
     * param[out] vertices          Vertices composing the line
     */
        void find_line_vertices_by_propagation(
            const RINGMesh::GeoModel3D &geomodel,
            const index_t incident_surface,
            const index_t first_vertex_id,
            const index_t second_vertex_id,
            std::vector<vec3> &vertices)
        {
            const RINGMesh::Surface3D &surface = geomodel.surface(incident_surface);
            vertices.push_back(surface.vertex(first_vertex_id));
            vertices.push_back(surface.vertex(second_vertex_id));

            bool propagation_toward_next = determine_propagation_direction(surface,
                                                                           first_vertex_id, second_vertex_id);

            const RINGMesh::SurfaceMesh3D &surface_mesh = surface.mesh();
            index_t cur_polygon_id = surface_mesh.polygon_from_vertex_ids(
                first_vertex_id, second_vertex_id);

            index_t cur_edge_id;
            if (propagation_toward_next)
            {
                cur_edge_id = surface_mesh.vertex_index_in_polygon(cur_polygon_id,
                                                                   first_vertex_id);
            }
            else
            {
                cur_edge_id = surface_mesh.vertex_index_in_polygon(cur_polygon_id,
                                                                   second_vertex_id);
            }

            bool line_last_vertex_is_a_corner = is_corner(geomodel,
                                                          vertices.back());

            while (!line_last_vertex_is_a_corner)
            {
                RINGMesh::PolygonLocalEdge next_polygon_edge;
                if (propagation_toward_next)
                {
                    next_polygon_edge = surface_mesh.next_on_border({cur_polygon_id, cur_edge_id});
                    index_t next_vertex_id = surface_mesh.polygon_vertex(surface_mesh.next_polygon_vertex(next_polygon_edge));
                    vertices.push_back(
                        surface_mesh.vertex(next_vertex_id));
                }
                else
                {
                    next_polygon_edge = surface_mesh.prev_on_border({cur_polygon_id, cur_edge_id});
                    index_t next_vertex_id = surface_mesh.polygon_vertex(next_polygon_edge);
                    vertices.push_back(
                        surface_mesh.vertex(next_vertex_id));
                }
                cur_polygon_id = next_polygon_edge.polygon_id;
                cur_edge_id = next_polygon_edge.local_edge_id;
                line_last_vertex_is_a_corner = is_corner(geomodel,
                                                         vertices.back());
            }
        }

        /*!
     * @brief If a point correspond to a corner of the geomodel, returns the
     * index of the corner, else NO_ID is returned
     * @todo use an aabb tree (optimization)
     */
        index_t find_corner(const RINGMesh::GeoModel3D &geomodel, const vec3 &point)
        {
            for (index_t i = 0; i < geomodel.nb_corners(); ++i)
            {
                if (geomodel.corner(i).vertex(0) == point)
                {
                    return i;
                }
            }
            return NO_ID;
        }
    };

}

namespace SCAR
{

    TSurfLineParser::TSurfLineParser(
        GeoModelBuilderTS &gm_builder,
        RINGMesh::GeoModel3D &geomodel)
        : RINGMesh::GocadBaseParser(gm_builder, geomodel)
    {
    }

    TSurfLoadingStorage::TSurfLoadingStorage()
    {
        cur_surf_first_vertex_index_ = GOCAD_OFFSET;
    }

    void GeoModelBuilderTS::load_file()
    {
        // Read file and import information contained herein
        read_file();

        remove_duplicate_lines();
        retrieve_geomodel_topology();
    }

    void GeoModelBuilderTS::read_line()
    {
        std::string keyword = file_line().field(0);
        std::unique_ptr<TSurfLineParser> tsurf_parser(
            TSurfLineParserFactory::create(keyword, *this, geomodel_));
        if (tsurf_parser)
        {
            tsurf_parser->execute(file_line(), ts_load_storage_);
        }
        else
        {
            std::unique_ptr<RINGMesh::GocadLineParser> gocad_parser =
                RINGMesh::GocadLineFactory::create(keyword, *this, geomodel_);
            if (gocad_parser)
            {
                gocad_parser->execute(file_line(), ts_load_storage_);
            }
        }
    }

    void GeoModelBuilderTS::remove_duplicate_lines()
    {
        std::set<RINGMesh::gmme_id> duplicate_lines_to_remove;
        find_duplicate_lines(duplicate_lines_to_remove);
        remove.remove_mesh_entities(duplicate_lines_to_remove);
    }

    void GeoModelBuilderTS::find_duplicate_lines(
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

    bool GeoModelBuilderTS::is_line_duplicated(const index_t line_id)
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

    bool GeoModelBuilderTS::are_lines_identical(
        const RINGMesh::Line3D &line1,
        const RINGMesh::Line3D &line2)
    {
        // To be identical 2 lines must have same number of points, same boundaries
        // (corners) and same other points.
        return same_number_vertices(line1, line2) && same_boundaries(line1, line2) && same_internal_vertices(line1, line2);
    }

    bool GeoModelBuilderTS::same_number_vertices(
        const RINGMesh::Line3D &line1,
        const RINGMesh::Line3D &line2)
    {
        return line1.nb_vertices() == line2.nb_vertices();
    }

    bool GeoModelBuilderTS::same_boundaries(
        const RINGMesh::Line3D &line1,
        const RINGMesh::Line3D &line2)
    {
        if (line1.boundary_gmme(0).index() != line2.boundary_gmme(0).index() && line1.boundary_gmme(0).index() != line2.boundary_gmme(1).index())
        {
            return false;
        }
        if (line1.boundary_gmme(1).index() != line2.boundary_gmme(0).index() && line1.boundary_gmme(1).index() != line2.boundary_gmme(1).index())
        {
            return false;
        }
        return true;
    }

    bool GeoModelBuilderTS::same_internal_vertices(
        const RINGMesh::Line3D &line1,
        const RINGMesh::Line3D &line2)
    {
        bool same_order = true;
        if (line1.boundary_gmme(1).index() == line2.boundary_gmme(0).index() && line1.boundary_gmme(0).index() == line2.boundary_gmme(1).index())
        {
            same_order = false;
        }
        index_t nb_line_vertices = line1.nb_vertices();
        for (index_t v = 1; v < nb_line_vertices - 1; ++v)
        {
            // If same order, compare vertex of same id. Else, we must compare,
            // first of line1 with last (of index nb_line_vertices - 1) of line2,
            // and so on...
            if ((same_order && line1.vertex(v) != line2.vertex(v)) || (!same_order && line1.vertex(v) != line2.vertex((nb_line_vertices - 1) - v)))
            {
                return false;
            }
        }
        return true;
    }

    // @todo Split this function in several functions
    void GeoModelBuilderTS::retrieve_geomodel_topology()
    {
        // Check and retrieve Corner-Line relationships
        // In the case of TSurf loading, this part is not useful.
        for (index_t c = 0; c < geomodel_.nb_corners(); ++c)
        {
            index_t cur_corner_geomodel_id = geomodel_.mesh.vertices.geomodel_vertex_id(
                RINGMesh::gmme_id(RINGMesh::Corner3D::type_name_static(), c));
            std::vector<RINGMesh::GMEVertex> line_gme_vertices =
                geomodel_.mesh.vertices.gme_type_vertices(
                    RINGMesh::Line3D::type_name_static(), cur_corner_geomodel_id);
            for (index_t gmev = 0; gmev < line_gme_vertices.size(); ++gmev)
            {
                index_t cur_line_id = line_gme_vertices[gmev].gmme.index();
                const RINGMesh::Line3D &cur_line = geomodel_.line(cur_line_id);
                bool corner_in_line_boundaries = false;
                for (index_t boundary = 0; boundary < cur_line.nb_boundaries();
                     ++boundary)
                {
                    if (cur_line.boundary_gmme(boundary).index() == c)
                    {
                        corner_in_line_boundaries = true;
                        break;
                    }
                }
                if (!corner_in_line_boundaries)
                {
                    topology.add_line_corner_boundary_relation(cur_line.index(), c);
                }
            }
        }

        // Check and retrieve Line-Surface relationships
        for (index_t l = 0; l < geomodel_.nb_lines(); ++l)
        {
            const RINGMesh::Line3D &cur_line = geomodel_.line(l);

            scar_assert(cur_line.nb_vertices() > 1);
            std::vector<index_t> incident_surfaces;
            if (cur_line.nb_vertices() == 2)
            {
                index_t cur_boundary0_geomodel_id =
                    geomodel_.mesh.vertices.geomodel_vertex_id(
                        cur_line.boundary_gmme(0));
                index_t cur_boundary1_geomodel_id =
                    geomodel_.mesh.vertices.geomodel_vertex_id(
                        cur_line.boundary_gmme(1));

                std::vector<RINGMesh::GMEVertex> cur_boundary0_surface_gme_vertices =
                    geomodel_.mesh.vertices.gme_type_vertices(
                        RINGMesh::Surface3D::type_name_static(),
                        cur_boundary0_geomodel_id);
                std::vector<RINGMesh::GMEVertex> cur_boundary1_surface_gme_vertices =
                    geomodel_.mesh.vertices.gme_type_vertices(
                        RINGMesh::Surface3D::type_name_static(),
                        cur_boundary1_geomodel_id);

                // Find incident surfaces to both line boundaries
                for (index_t gmev_b0 = 0;
                     gmev_b0 < cur_boundary0_surface_gme_vertices.size(); ++gmev_b0)
                {
                    index_t cur_surface_id_b0 =
                        cur_boundary0_surface_gme_vertices[gmev_b0].gmme.index();
                    for (index_t gmev_b1 = 0;
                         gmev_b1 < cur_boundary1_surface_gme_vertices.size();
                         ++gmev_b1)
                    {
                        index_t cur_surface_id_b1 =
                            cur_boundary1_surface_gme_vertices[gmev_b1].gmme.index();
                        if (cur_surface_id_b1 == cur_surface_id_b0)
                        {
                            incident_surfaces.push_back(cur_surface_id_b0);
                            break;
                        }
                    }
                }
            }
            else
            {
                // Check which surfaces are incident to any internal point
                index_t line_vertex_geomodel_id =
                    geomodel_.mesh.vertices.geomodel_vertex_id(cur_line.gmme(), 1);
                std::vector<RINGMesh::GMEVertex> surface_gme_vertices =
                    geomodel_.mesh.vertices.gme_type_vertices(
                        RINGMesh::Surface3D::type_name_static(),
                        line_vertex_geomodel_id);
                for (index_t gmev = 0; gmev < surface_gme_vertices.size(); ++gmev)
                {
                    index_t cur_surface_id = surface_gme_vertices[gmev].gmme.index();
                    incident_surfaces.push_back(cur_surface_id);
                }
            }

            for (index_t surf = 0; surf < incident_surfaces.size(); ++surf)
            {
                const RINGMesh::Surface3D &cur_surface = geomodel_.surface(
                    incident_surfaces[surf]);
                index_t nb_times_line_is_boundary = static_cast<index_t>(std::count(
                    incident_surfaces.begin(), incident_surfaces.end(),
                    incident_surfaces[surf]));
                scar_assert(
                    nb_times_line_is_boundary == 1 || nb_times_line_is_boundary == 2);
                bool line_in_surface_boundaries = false;
                for (index_t boundary = 0; boundary < cur_surface.nb_boundaries();
                     ++boundary)
                {
                    if (cur_surface.boundary_gmme(boundary).index() == l)
                    {
                        line_in_surface_boundaries = true;
                        break;
                    }
                }
                if (!line_in_surface_boundaries)
                {
                    topology.add_surface_line_boundary_relation(cur_surface.index(),
                                                                cur_line.index());
                }
                if (nb_times_line_is_boundary == 2)
                {
                    // This line must be an inside border of the surface
                    // @todo WARNING Not sure it will works (porting to new ringmesh)
                    if (!cur_line.is_inside_border(cur_surface))
                    {
                        topology.add_surface_line_boundary_relation(cur_surface.index(),
                                                                    cur_line.index());
                    }
                }
            }
        }
    }

    void GeoModelBuilderTS::save_geomodel(
        const RINGMesh::GeoModel2D &geomodel,
        const std::string &filename)
    {
        std::ofstream out(filename.c_str());
        out.precision(16);

        // Print Model3d headers
        out << "GOCAD TSurf 1" << std::endl
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

        for (index_t surface_id = 0; surface_id < geomodel.nb_surfaces();
             ++surface_id)
        {
            out << "TFACE" << std::endl;
            const RINGMesh::Surface2D &cur_surface = geomodel.surface(surface_id);
            index_t begin_vertex_count = vertex_count;
            if (!cur_surface.is_meshed() || !cur_surface.is_simplicial())
            {
                throw SCARException("I/O", "Only model with triangular surfaces can "
                                           "be saved in .ts file format ");
            }
            for (index_t v = 0; v < cur_surface.nb_vertices(); ++v)
            {
                out << "VRTX " << ++vertex_count << " " << cur_surface.vertex(v).x
                    << " " << cur_surface.vertex(v).y << " " << 0 << std::endl;
            }
            for (index_t t = 0; t < cur_surface.nb_mesh_elements(); ++t)
            {
                out << "TRGL "
                    << cur_surface.mesh_element_vertex_index({t, 0}) + begin_vertex_count + 1 << " "
                    << cur_surface.mesh_element_vertex_index({t, 1}) + begin_vertex_count + 1 << " "
                    << cur_surface.mesh_element_vertex_index({t, 2}) + begin_vertex_count + 1 << std::endl;
            }
        }
        out << "END" << std::endl;
    }

    void GeoModelBuilderTS::save_geomodel(
        const RINGMesh::GeoModel3D &geomodel,
        const std::string &filename)
    {
        std::ofstream out(filename.c_str());
        out.precision(16);

        // Print Model3d headers
        out << "GOCAD TSurf 1" << std::endl
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

        for (index_t surface_id = 0; surface_id < geomodel.nb_surfaces();
             ++surface_id)
        {
            out << "TFACE" << std::endl;
            const RINGMesh::Surface3D &cur_surface = geomodel.surface(surface_id);
            index_t begin_vertex_count = vertex_count;
            if (!cur_surface.is_simplicial())
            {
                throw SCARException("I/O", "Only model with triangular surfaces can "
                                           "be saved in .ts file format ");
            }
            for (index_t v = 0; v < cur_surface.nb_vertices(); ++v)
            {
                out << "VRTX " << ++vertex_count << " " << cur_surface.vertex(v).x
                    << " " << cur_surface.vertex(v).y << " "
                    << cur_surface.vertex(v).z << std::endl;
            }
            for (index_t t = 0; t < cur_surface.nb_mesh_elements(); ++t)
            {
                out << "TRGL "
                    << cur_surface.mesh_element_vertex_index({t, 0}) + begin_vertex_count + 1 << " "
                    << cur_surface.mesh_element_vertex_index({t, 1}) + begin_vertex_count + 1 << " "
                    << cur_surface.mesh_element_vertex_index({t, 2}) + begin_vertex_count + 1 << std::endl;
            }
        }
        out << "END" << std::endl;
    }

    void tsurf_import_factory_initialize()
    {
        TSurfLineParserFactory::register_creator<LoadTSurfNewGObj>("GOCAD");
        TSurfLineParserFactory::register_creator<LoadTSurfName>("name:");
        TSurfLineParserFactory::register_creator<LoadTSurfSurface>("TFACE");
        TSurfLineParserFactory::register_creator<LoadTSurfEndInterface>("END");
        TSurfLineParserFactory::register_creator<LoadTSurfCorner>("BSTONE");
        TSurfLineParserFactory::register_creator<LoadTSurfTriangle>("TRGL");
        TSurfLineParserFactory::register_creator<LoadTSurfLine>("BORDER");
        TSurfLineParserFactory::register_creator<LoadTSurfAtom>("ATOM");
        TSurfLineParserFactory::register_creator<LoadTSurfAtom>("PATOM");
    }

}
