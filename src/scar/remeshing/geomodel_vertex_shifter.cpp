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

#include <scar/remeshing/geomodel_vertex_shifter.h>

#include <geogram/basic/attributes.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/io/io.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>

#include <ringpcl/correlation_map.h>

#include <scar/repair/topology_recovery.h>
#include <scar/tools/geometry.h>
#include <scar/tools/utils.h>

namespace
{
    using namespace SCAR;

    template <index_t DIMENSION>
    vecn<DIMENSION> generate_random_vertex(double max_coord)
    {
        vecn<DIMENSION> random_point;
        for (auto c : RINGMesh::range(DIMENSION))
        {
            random_point[c] = max_coord * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
        }
        return random_point;
    }

    bool check_polygon(const RINGMesh::Surface2D &surface, index_t polygon_id)
    {
        auto v0_to_v1 = GEO::normalize(
            surface.mesh_element_vertex({polygon_id, 1}) - surface.mesh_element_vertex({polygon_id, 0}));
        auto v0_to_v2 = GEO::normalize(
            surface.mesh_element_vertex({polygon_id, 2}) - surface.mesh_element_vertex({polygon_id, 0}));
        vec3 v01{v0_to_v1.x, v0_to_v1.y, 0};
        vec3 v02{v0_to_v2.x, v0_to_v2.y, 0};
        return cross(v01, v02).z > 0;
    }

    bool check_polygon(const RINGMesh::Surface3D &surface, index_t polygon_id)
    {
        scar_unused(surface);
        scar_unused(polygon_id);
        Logger::out("Check polygon", "Not implemented in 3D");
        return true;
    }
}

namespace SCAR
{

    template <index_t DIMENSION>
    GeoModelVertexShifter<DIMENSION>::GeoModelVertexShifter(
        RINGMesh::GeoModel<DIMENSION> &output_geomodel,
        const GeoModelTopologyRecoverer<DIMENSION> &topology_recoverer,
        ConvexShapeType convex_type,
        const double &angle,
        bool verbose)
        : geomodel_(output_geomodel),
          topology_recoverer_(topology_recoverer),
          shape_type_(convex_type),
          angle_(angle),
          verbose_(verbose)
    {
    }

    template <index_t DIMENSION>
    void GeoModelVertexShifter<DIMENSION>::shift_model_vertices()
    {
        edit_corners();
        repel_corners_from_lines();
        edit_lines();
    }

    template <index_t DIMENSION>
    class CornerMaxClique
    {
    public:
        CornerMaxClique(
            const std::vector<index_t> &corner_ids,
            const std::vector<double> &corner_excls,
            const vecn<DIMENSION> barycenter)
            : ids(corner_ids), excl(corner_excls), bary(barycenter)
        {
            scar_assert(corner_ids.size() == corner_excls.size());
            cardinality = corner_ids.size();
        }

        void set_processed()
        {
            scar_assert(status);
            status = false;
        }

    public:
        std::vector<index_t> ids;
        std::vector<double> excl;
        vecn<DIMENSION> bary;
        size_t cardinality;
        bool status = true;
    };

    template <index_t DIMENSION>
    std::vector<CornerMaxClique<DIMENSION>> build_corner_max_cliques(
        const std::vector<index_t> &connected_components,
        const CorrelationMap<bool> &connectivity_and_invalidity,
        const RINGMesh::GeoModel<DIMENSION> &geomodel,
        ConvexShapeType shape_type)
    {
        std::vector<CornerMaxClique<DIMENSION>> results;
        CorrelationMap<bool> submatrix(
            static_cast<index_t>(connected_components.size()), false);
        for (auto i : RINGMesh::range(connected_components.size() - 1))
        {
            for (auto j : RINGMesh::range(i + 1, connected_components.size()))
            {
                submatrix(i, j) = connectivity_and_invalidity(
                    connected_components[i], connected_components[j]);
            }
        }

        auto max_cliques = MaximalCliquesAPI::maximal_cliques_listing(submatrix);
        for (const auto &max_clique : max_cliques)
        {
            std::vector<index_t> corner_ids;
            corner_ids.reserve(max_clique.size());
            std::vector<double> exclusions;
            exclusions.reserve(max_clique.size());
            vecn<DIMENSION> bary;
            for (auto i_id : max_clique)
            {
                bary += geomodel.corner(connected_components[i_id]).vertex(0);
                corner_ids.push_back(connected_components[i_id]);
                auto corner_exclusion = get_exclusion_shape(shape_type,
                                                            geomodel.corner(connected_components[i_id]).mesh(), 0);
                exclusions.push_back(corner_exclusion->width(bary));
            }
            bary /= max_clique.size();
            results.emplace_back(corner_ids, exclusions, bary);
        }

        return results;
    }

    template <index_t DIMENSION>
    index_t get_next_max_clique(
        const std::vector<CornerMaxClique<DIMENSION>> &max_cliques)
    {
        index_t id = NO_ID;
        for (auto i : RINGMesh::range(max_cliques.size()))
        {
            if (!max_cliques[i].status)
            {
                continue;
            }
            if (id == NO_ID)
            { // Init
                id = i;
                continue;
            }
            if (max_cliques[i].cardinality < max_cliques[id].cardinality)
            {
                id = i;
            }
        }
        return id;
    }

    template <index_t DIMENSION>
    void process_max_clique(
        RINGMesh::GeoModel<DIMENSION> &geomodel,
        std::vector<CornerMaxClique<DIMENSION>> &max_cliques,
        index_t id)
    {
        std::vector<bool> corners_on_voi(geomodel.nb_corners(), false);
        for (const auto &line : geomodel.lines())
        {
            if (line.is_on_voi())
            {
                corners_on_voi[line.boundary(0).index()] = true;
                corners_on_voi[line.boundary(1).index()] = true;
            }
        }

        auto &cur_clique = max_cliques[id];
        // Compute freedom corners
        std::vector<bool> freedom(cur_clique.cardinality, true);
        for (auto i : RINGMesh::range(cur_clique.cardinality))
        {
            auto c_id = cur_clique.ids[i];
            if (corners_on_voi[c_id])
            {
                freedom[i] = false;
                continue;
            }

            for (auto j : RINGMesh::range(max_cliques.size()))
            {
                if (j == id)
                {
                    continue;
                }
                const auto &other_clique = max_cliques[j];
                if (other_clique.status)
                {
                    for (auto ji : RINGMesh::range(other_clique.cardinality))
                    {
                        if (other_clique.ids[ji] == c_id)
                        {
                            freedom[i] = false;
                            break;
                        }
                    }
                }
                if (!freedom[i])
                {
                    break;
                }
            }
        }

        // Move freedom corners and potentially associated cliques
        RINGMesh::GeoModelBuilder<DIMENSION> builder(geomodel);
        std::vector<double> max_penetrations(cur_clique.cardinality, 0);
        for (auto i : RINGMesh::range(cur_clique.cardinality))
        {
            if (!freedom[i])
            {
                continue;
            }
            for (auto j : RINGMesh::range(cur_clique.cardinality))
            {
                if (i == j)
                {
                    continue;
                }
                double cur_dist = (geomodel.corner(cur_clique.ids[i]).vertex(0) - geomodel.corner(cur_clique.ids[j]).vertex(0)).length();
                double cur_obj = cur_clique.excl[i] + cur_clique.excl[j];
                double diff = cur_obj - cur_dist;
                if (freedom[j])
                {
                    diff *= 0.5;
                }
                max_penetrations[i] = std::max(max_penetrations[i], diff);
            }
        }
        // For pairs of corners on VOI
        std::vector<index_t> pairs_on_voi; // not a good struct here
        for (auto i : RINGMesh::range(cur_clique.cardinality))
        {
            if (freedom[i])
            {
                continue;
            }
            for (auto j : RINGMesh::range(cur_clique.cardinality))
            {
                if (i == j)
                {
                    continue;
                }
                if (corners_on_voi[cur_clique.ids[i]] && corners_on_voi[cur_clique.ids[j]])
                {
                    double cur_dist =
                        (geomodel.corner(cur_clique.ids[i]).vertex(0) - geomodel.corner(cur_clique.ids[j]).vertex(0)).length();
                    double cur_obj = cur_clique.excl[i] + cur_clique.excl[j];
                    double diff = cur_obj - cur_dist;
                    max_penetrations[i] = std::max(max_penetrations[i], diff);
                    pairs_on_voi.push_back(i);
                }
            }
        }

        for (auto i : RINGMesh::range(cur_clique.cardinality))
        {
            vecn<DIMENSION> displacement;
            displacement = normalize(
                geomodel.corner(cur_clique.ids[i]).vertex(0) - cur_clique.bary);
            if (RINGMesh::contains(pairs_on_voi, i))
            {
                scar_assert(pairs_on_voi.size() == 2);
                displacement = normalize(
                    geomodel.corner(cur_clique.ids[pairs_on_voi[0]]).vertex(0) - geomodel.corner(cur_clique.ids[pairs_on_voi[1]]).vertex(0));
                if (dot(displacement,
                        geomodel.corner(cur_clique.ids[i]).vertex(0) - cur_clique.bary) < 0)
                {
                    displacement *= -1;
                }
            }
            displacement *= max_penetrations[i];

            std::set<index_t> corners_moved;
            std::queue<index_t> corners_to_move;
            corners_to_move.push(cur_clique.ids[i]);
            while (!corners_to_move.empty())
            {
                auto cur_c_id = corners_to_move.front();
                corners_to_move.pop();

                // Apply displacement
                auto new_location = geomodel.corner(cur_c_id).vertex(0) + displacement;
                builder.geometry.set_mesh_entity_vertex(
                    geomodel.corner(cur_c_id).gmme(), 0, new_location, true);
                corners_moved.insert(cur_c_id);

                // Propagate on alredy processed cliques (ie status == false)
                for (auto j : RINGMesh::range(max_cliques.size()))
                {
                    if (j == id)
                    {
                        continue;
                    }
                    if (max_cliques[j].status)
                    {
                        continue;
                    }
                    if (RINGMesh::contains(max_cliques[j].ids, cur_c_id))
                    {
                        for (auto k : RINGMesh::range(max_cliques[j].cardinality))
                        {
                            if (!set_contains(max_cliques[j].ids[k],
                                              corners_moved))
                            {
                                corners_to_move.push(max_cliques[j].ids[k]);
                            }
                        }
                    }
                }
            }
        }

        // Set process
        cur_clique.set_processed();
    }

    template <index_t DIMENSION>
    void GeoModelVertexShifter<DIMENSION>::edit_corners()
    {
        if (enum_contains(topology_recoverer_.geological_rules,
                          GeologicalRule::REMESH_DFN))
        {
            CorrelationMap<bool> connectivity(geomodel_.nb_corners(), false);
            for (const auto &line : geomodel_.lines())
            {
                connectivity(line.boundary(0).index(),
                             line.boundary(1).index()) = true;
            }

            CorrelationMap<bool> connectivity_and_invalidity(
                geomodel_.nb_corners(), false);
            for (auto c1_id : RINGMesh::range(geomodel_.nb_corners() - 1))
            {
                for (auto c2_id : RINGMesh::range(c1_id + 1,
                                                  geomodel_.nb_corners()))
                {
                    if (topology_recoverer_.invalidity_graph(c1_id, c2_id) != NO_ID && connectivity(c1_id, c2_id))
                    {
                        connectivity_and_invalidity(c1_id, c2_id) = true;
                    }
                }
            }

            std::vector<std::vector<index_t>> connected_components;
            CorrelationMapAPI::get_connected_components(connectivity_and_invalidity,
                                                        connected_components, true);
            for (auto &cc : connected_components)
            {
                if (cc.size() == 1)
                {
                    continue;
                }
                auto max_cliques = build_corner_max_cliques(cc,
                                                            connectivity_and_invalidity, geomodel_, shape_type_);
                auto clique_id = get_next_max_clique(max_cliques);
                while (clique_id != NO_ID)
                {
                    // Process the max clique
                    process_max_clique(geomodel_, max_cliques, clique_id);

                    // For next step
                    clique_id = get_next_max_clique(max_cliques);
                }
            }
        }

        std::set<index_t> modified_vertices;
        vecn<DIMENSION> ref_right;
        ref_right[0] = 1.;
        RINGMesh::GeoModelBuilder<DIMENSION> builder(geomodel_);

        // For each corner
        for (auto c_id : RINGMesh::range(geomodel_.nb_corners()))
        {
            const auto &corner = geomodel_.corner(c_id);
            auto A = geomodel_.corner(c_id).vertex(0);

            // Get incident lines and tags them if they are on voi
            std::set<index_t> incident_lines;
            std::set<index_t> incident_lines_on_voi;
            for (auto il_itr : RINGMesh::range(corner.nb_incident_entities()))
            {
                incident_lines.insert(corner.incident_entity(il_itr).index());
                if (corner.incident_entity(il_itr).is_on_voi())
                {
                    incident_lines_on_voi.insert(
                        corner.incident_entity(il_itr).index());
                }
            }

            // Ordering lines from a line on the voi (if any)
            std::vector<index_t> ordered_lines;
            CorrelationMap<bool> line_linked_by_surfaces(geomodel_.nb_lines(),
                                                         false);
            for (const auto &surf : geomodel_.surfaces())
            {
                std::set<index_t> boundary_lines_incident_to_corner;
                for (auto b : RINGMesh::range(surf.nb_boundaries()))
                {
                    if (set_contains(surf.boundary(b).index(),
                                     incident_lines))
                    {
                        boundary_lines_incident_to_corner.insert(
                            surf.boundary(b).index());
                    }
                }
                if (boundary_lines_incident_to_corner.size() != 2)
                {
                    continue;
                }

                line_linked_by_surfaces(*boundary_lines_incident_to_corner.begin(),
                                        *std::next(boundary_lines_incident_to_corner.begin())) = true;
            }
            if (!incident_lines_on_voi.empty())
            {
                ordered_lines.push_back(*incident_lines_on_voi.begin());
                while (true)
                {
                    auto linked = CorrelationMapAPI::directly_linked(
                        line_linked_by_surfaces, ordered_lines.back(), false,
                        false);
                    if (linked.empty())
                    {
                        break;
                    }
                    if (linked.size() == 1)
                    {
                        if (RINGMesh::contains(ordered_lines, linked[0]))
                        {
                            break;
                        }
                        ordered_lines.push_back(linked[0]);
                    }
                    else
                    {
                        for (auto l : linked)
                        {
                            if (!RINGMesh::contains(ordered_lines, l))
                            {
                                ordered_lines.push_back(l);
                                break;
                            }
                        }
                    }
                }
            }
            else
            {
                for (auto i : incident_lines)
                {
                    ordered_lines.push_back(i);
                }
            }

            // Compute angle between pairs of lines
            for (auto i : RINGMesh::range(ordered_lines.size()))
            {
                auto j = (i + 1) % ordered_lines.size();
                vecn<DIMENSION> B;
                index_t v_in_line0;
                if (geomodel_.line(ordered_lines[i]).boundary(0).gmme() == corner.gmme())
                {
                    v_in_line0 = 1;
                    B = geomodel_.line(ordered_lines[i]).vertex(v_in_line0);
                }
                else
                {
                    v_in_line0 = geomodel_.line(ordered_lines[i]).nb_vertices() - 2;
                    B = geomodel_.line(ordered_lines[i]).vertex(v_in_line0);
                }
                vecn<DIMENSION> C;
                index_t v_in_line1;
                if (geomodel_.line(ordered_lines[j]).boundary(0).gmme() == corner.gmme())
                {
                    v_in_line1 = 1;
                    C = geomodel_.line(ordered_lines[j]).vertex(v_in_line1);
                }
                else
                {
                    v_in_line1 = geomodel_.line(ordered_lines[j]).nb_vertices() - 2;
                    C = geomodel_.line(ordered_lines[j]).vertex(v_in_line1);
                }
                double cur_angle = angle(B - A, C - A);
                if (cur_angle >= angle_)
                {
                    continue;
                }
                double lAB = (B - A).length();
                double lAC = (C - A).length();
                double lBC = (B - C).length();
                double lBD = lAB * lBC / (lAB + lAC);
                auto D = B + lBD * normalize(C - B);

                double new_angle0 = angle_ / 2;
                double new_angle1 = angle_ / 2;
                if (enum_contains(topology_recoverer_.geological_rules,
                                  GeologicalRule::KEEP_BOUNDARIES_UNMOVED))
                {
                    if (geomodel_.line(ordered_lines[i]).is_on_voi())
                    {
                        new_angle1 = angle_ - cur_angle / 2;
                        new_angle0 = cur_angle / 2;
                    }
                    if (geomodel_.line(ordered_lines[j]).is_on_voi())
                    {
                        new_angle0 = angle_ - cur_angle / 2;
                        new_angle1 = cur_angle / 2;
                    }
                }
                auto Bj = A + normalize(D - A) * lAB * cos(cur_angle / 2);
                auto Bi = A + normalize(D - A) * lAB * cos(new_angle0);
                auto B_new = Bi + normalize(B - Bj) * lAB * sin(new_angle0);
                auto Cj = A + normalize(D - A) * lAC * cos(cur_angle / 2);
                auto Ci = A + normalize(D - A) * lAC * cos(new_angle1);
                auto C_new = Ci + normalize(C - Cj) * lAC * sin(new_angle1);

                builder.geometry.set_mesh_entity_vertex(
                    geomodel_.line(ordered_lines[i]).gmme(), v_in_line0, B_new,
                    true);
                builder.geometry.set_mesh_entity_vertex(
                    geomodel_.line(ordered_lines[j]).gmme(), v_in_line1, C_new,
                    true);

                // @todo check on triangles
                index_t modified_v0 = geomodel_.mesh.vertices.geomodel_vertex_id(
                    geomodel_.line(ordered_lines[i]).gmme(), v_in_line0);
                index_t modified_v1 = geomodel_.mesh.vertices.geomodel_vertex_id(
                    geomodel_.line(ordered_lines[j]).gmme(), v_in_line1);
                modified_vertices.insert(modified_v0);
                modified_vertices.insert(modified_v1);

                auto gme_vertices_v0 = geomodel_.mesh.vertices.gme_type_vertices(
                    RINGMesh::surface_type_name_static(), modified_v0);
                for (auto gmev : gme_vertices_v0)
                {
                    const auto &cur_surf = geomodel_.surface(gmev.gmme.index());
                    auto polygons_to_check = cur_surf.mesh().polygons_around_vertex(
                        gmev.v_index, false, NO_ID);
                    for (auto polygon : polygons_to_check)
                    {
                        bool ok = check_polygon(cur_surf, polygon);
                    }
                }
                auto gme_vertices_v1 = geomodel_.mesh.vertices.gme_type_vertices(
                    RINGMesh::surface_type_name_static(), modified_v1);
                for (auto gmev : gme_vertices_v1)
                {
                    const auto &cur_surf = geomodel_.surface(gmev.gmme.index());
                    auto polygons_to_check = cur_surf.mesh().polygons_around_vertex(
                        gmev.v_index, false, NO_ID);
                    for (auto polygon : polygons_to_check)
                    {
                        bool ok = check_polygon(cur_surf, polygon);
                    }
                }
            }

            // Do the same in the inverse order
            std::vector<index_t> inverse_ordered_lines(ordered_lines.size());
            for (auto i : RINGMesh::range(ordered_lines.size()))
            {
                inverse_ordered_lines[i] =
                    ordered_lines[ordered_lines.size() - 1 - i];
            }
            for (auto i : RINGMesh::range(inverse_ordered_lines.size()))
            {
                auto j = (i + 1) % inverse_ordered_lines.size();
                vecn<DIMENSION> B;
                index_t v_in_line0;
                if (geomodel_.line(inverse_ordered_lines[i]).boundary(0).gmme() == corner.gmme())
                {
                    v_in_line0 = 1;
                    B = geomodel_.line(inverse_ordered_lines[i]).vertex(v_in_line0);
                }
                else
                {
                    v_in_line0 =
                        geomodel_.line(inverse_ordered_lines[i]).nb_vertices() - 2;
                    B = geomodel_.line(inverse_ordered_lines[i]).vertex(v_in_line0);
                }
                vecn<DIMENSION> C;
                index_t v_in_line1;
                if (geomodel_.line(inverse_ordered_lines[j]).boundary(0).gmme() == corner.gmme())
                {
                    v_in_line1 = 1;
                    C = geomodel_.line(inverse_ordered_lines[j]).vertex(v_in_line1);
                }
                else
                {
                    v_in_line1 =
                        geomodel_.line(inverse_ordered_lines[j]).nb_vertices() - 2;
                    C = geomodel_.line(inverse_ordered_lines[j]).vertex(v_in_line1);
                }
                double cur_angle = angle(B - A, C - A);
                if (cur_angle >= angle_)
                {
                    continue;
                }
                double lAB = (B - A).length();
                double lAC = (C - A).length();
                double lBC = (B - C).length();
                double lBD = lAB * lBC / (lAB + lAC);
                auto D = B + lBD * normalize(C - B);

                double new_angle0 = angle_ / 2;
                double new_angle1 = angle_ / 2;
                if (enum_contains(topology_recoverer_.geological_rules,
                                  GeologicalRule::KEEP_BOUNDARIES_UNMOVED))
                {
                    if (geomodel_.line(inverse_ordered_lines[i]).is_on_voi())
                    {
                        new_angle1 = angle_ - cur_angle / 2;
                        new_angle0 = cur_angle / 2;
                    }
                    if (geomodel_.line(inverse_ordered_lines[j]).is_on_voi())
                    {
                        new_angle0 = angle_ - cur_angle / 2;
                        new_angle1 = cur_angle / 2;
                    }
                }
                auto Bj = A + normalize(D - A) * lAB * cos(cur_angle / 2);
                auto Bi = A + normalize(D - A) * lAB * cos(new_angle0);
                auto B_new = Bi + normalize(B - Bj) * lAB * sin(new_angle0);
                auto Cj = A + normalize(D - A) * lAC * cos(cur_angle / 2);
                auto Ci = A + normalize(D - A) * lAC * cos(new_angle1);
                auto C_new = Ci + normalize(C - Cj) * lAC * sin(new_angle1);

                builder.geometry.set_mesh_entity_vertex(
                    geomodel_.line(inverse_ordered_lines[i]).gmme(), v_in_line0,
                    B_new, true);
                builder.geometry.set_mesh_entity_vertex(
                    geomodel_.line(inverse_ordered_lines[j]).gmme(), v_in_line1,
                    C_new, true);

                // @todo check on triangles
                index_t modified_v0 = geomodel_.mesh.vertices.geomodel_vertex_id(
                    geomodel_.line(inverse_ordered_lines[i]).gmme(), v_in_line0);
                index_t modified_v1 = geomodel_.mesh.vertices.geomodel_vertex_id(
                    geomodel_.line(inverse_ordered_lines[j]).gmme(), v_in_line1);
                modified_vertices.insert(modified_v0);
                modified_vertices.insert(modified_v1);
                auto gme_vertices_v0 = geomodel_.mesh.vertices.gme_type_vertices(
                    RINGMesh::surface_type_name_static(), modified_v0);
                for (auto gmev : gme_vertices_v0)
                {
                    const auto &cur_surf = geomodel_.surface(gmev.gmme.index());
                    auto polygons_to_check = cur_surf.mesh().polygons_around_vertex(
                        gmev.v_index, false, NO_ID);
                    for (auto polygon : polygons_to_check)
                    {
                        bool ok = check_polygon(cur_surf, polygon);
                    }
                }
                auto gme_vertices_v1 = geomodel_.mesh.vertices.gme_type_vertices(
                    RINGMesh::surface_type_name_static(), modified_v1);
                for (auto gmev : gme_vertices_v1)
                {
                    const auto &cur_surf = geomodel_.surface(gmev.gmme.index());
                    auto polygons_to_check = cur_surf.mesh().polygons_around_vertex(
                        gmev.v_index, false, NO_ID);
                    for (auto polygon : polygons_to_check)
                    {
                        bool ok = check_polygon(cur_surf, polygon);
                    }
                }
            }
        }

        index_t nb_triangles_modified = 0;
        for (auto t : RINGMesh::range(geomodel_.mesh.polygons.nb()))
        {
            for (auto i : RINGMesh::range(3))
            {
                if (modified_vertices.find(
                        geomodel_.mesh.polygons.vertex({t, i})) != modified_vertices.end())
                {
                    ++nb_triangles_modified;
                    break;
                }
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelVertexShifter<DIMENSION>::repel_corners_from_lines()
    {
        RINGMesh::GeoModelBuilder<DIMENSION> builder(geomodel_);
        for (const auto &invalid_edge : topology_recoverer_.edge_information)
        {
            if (topology_recoverer_.get_edge_type(invalid_edge) != EdgeType::CornerLine)
            {
                continue;
            }
            const auto &corner = geomodel_.corner(invalid_edge.node1.id.index());
            const auto &line = geomodel_.line(invalid_edge.node2.id.index());
            index_t closest_edge_id;
            vecn<DIMENSION> closest_point;
            double distance;
            std::tie(closest_edge_id, closest_point, distance) =
                line.edge_aabb().closest_edge(corner.vertex(0));
            auto direction = normalize(closest_point - corner.vertex(0));
            double target_distance = 0;
            auto corner_exclusion = get_exclusion_shape(shape_type_, corner.mesh(),
                                                        0);
            target_distance += corner_exclusion->width(direction);
            auto rel_coord =
                PointOnLineCoordinates<DIMENSION>::relative_coordinates(
                    closest_point, line.mesh(), closest_edge_id);
            auto line_exclusion = get_exclusion_shape(shape_type_, line.mesh(),
                                                      rel_coord);
            target_distance += line_exclusion->width(-direction);

            double displacement = target_distance - distance;

            // Compute corner new location
            vecn<DIMENSION> new_location = corner.vertex(0) - (target_distance - distance) * 0.5 * direction;

            // Apply it
            builder.geometry.set_mesh_entity_vertex(corner.gmme(), 0, new_location,
                                                    true);

            // Propagate displacement along lines
            if (enum_contains(topology_recoverer_.geological_rules,
                              GeologicalRule::REMESH_DFN))
            {
                if (corner.nb_incident_entities() != 1)
                {
                    continue;
                }
                const auto &line = corner.incident_entity(0);
                if (line.nb_vertices() < 4)
                {
                    continue;
                }
                double displacement_init = displacement;
                bool increment = true;
                if (line.boundary_gmme(1) == corner.gmme())
                {
                    increment = false;
                }
                index_t cur_v_id = increment ? 1 : line.nb_vertices() - 2;
                while (displacement > displacement_init / 10.)
                {
                    vecn<DIMENSION> new_location = line.vertex(cur_v_id) - displacement * 0.5 * direction;
                    displacement /= 2.;
                    builder.geometry.set_mesh_entity_vertex(line.gmme(), cur_v_id,
                                                            new_location, true);
                    cur_v_id += increment ? 1 : -1;
                    if (cur_v_id <= 0 || cur_v_id >= line.nb_vertices())
                    {
                        break;
                    }
                }
            }
        }
    }

    template <index_t DIMENSION>
    void GeoModelVertexShifter<DIMENSION>::edit_lines()
    {
        RINGMesh::GeoModelBuilder<DIMENSION> builder(geomodel_);

        std::set<index_t> modified_vertices;

        for (auto l_id : RINGMesh::range(geomodel_.nb_lines()))
        {
            if (enum_contains(topology_recoverer_.geological_rules,
                              GeologicalRule::KEEP_BOUNDARIES_UNMOVED) &&
                geomodel_.line(l_id).is_on_voi())
            {
                // This line should not be moved
                continue;
            }

            index_t line_node_id = l_id + geomodel_.nb_corners();
            auto directly_linked = CorrelationMapAPI::directly_linked(
                topology_recoverer_.invalidity_graph, line_node_id, NO_ID, false);
            for (auto node_id : directly_linked)
            {
                if (topology_recoverer_.node_information[node_id].type != RINGMesh::line_type_name_static())
                {
                    continue;
                }
                EdgeInformation edge_info =
                    topology_recoverer_.edge_information[topology_recoverer_.invalidity_graph(
                        node_id, line_node_id)];
                const auto &line = geomodel_.line(l_id);
                GEO::Attribute<double> exclusion_attribute(
                    line.mesh().vertex_attribute_manager(), tolerance_att_name);
                index_t other_line_id;
                if (edge_info.node1.id.index() != l_id)
                {
                    other_line_id = edge_info.node1.id.index();
                }
                else
                {
                    other_line_id = edge_info.node2.id.index();
                }
                const auto &other_line = geomodel_.line(other_line_id);
                GEO::Attribute<double> exclusion_attribute_other(
                    other_line.mesh().vertex_attribute_manager(),
                    tolerance_att_name);
                for (const auto &part : edge_info.node1.parts)
                {
                    for (auto v_id : part.get_range())
                    {
                        if (v_id <= 1)
                        {
                            continue;
                        }
                        if (v_id >= line.nb_vertices() - 2)
                        {
                            continue;
                        }

                        // Find closest vertex in the other line and compute the distance
                        auto old_location = geomodel_.line(l_id).vertex(v_id);
                        auto exclusion_distance = exclusion_attribute[v_id];
                        index_t closest_edge_id;
                        vecn<DIMENSION> closest_point;
                        double old_distance;
                        std::tie(closest_edge_id, closest_point, std::ignore) =
                            other_line.edge_aabb().closest_edge(old_location);

                        double prop1_edge = 1. - ((closest_point - other_line.mesh_element_vertex({closest_edge_id, 0})).length() / other_line.mesh_element_size(closest_edge_id));
                        double prop2_edge = 1. - prop1_edge;
                        auto exclusion_distance_other1 =
                            exclusion_attribute_other[other_line.mesh_element_vertex_index(
                                {closest_edge_id, 0})];
                        auto exclusion_distance_other2 =
                            exclusion_attribute_other[other_line.mesh_element_vertex_index(
                                {closest_edge_id, 1})];
                        auto exclusion_distance_other = exclusion_distance_other1;
                        if (prop2_edge > prop1_edge)
                        {
                            exclusion_distance_other = exclusion_distance_other2;
                        }
                        index_t closest_vertex = closest_edge_id;
                        if (prop2_edge > prop1_edge)
                        {
                            closest_vertex++;
                        }
                        closest_point = other_line.vertex(closest_vertex);
                        old_distance = (closest_point - old_location).length();
                        // Condition is respected
                        if (old_distance >= (exclusion_distance + exclusion_distance_other))
                        {
                            continue;
                        }

                        // Find magnitude of displacement
                        double disp;
                        if (enum_contains(topology_recoverer_.geological_rules,
                                          GeologicalRule::KEEP_BOUNDARIES_UNMOVED) &&
                            other_line.is_on_voi())
                        {
                            // The other cannot move so we have to
                            // displaced this line for the whole missing distance
                            disp = (exclusion_distance + exclusion_distance_other) - old_distance;
                        }
                        else
                        {
                            if (other_line_id < l_id)
                            {
                                // We already moved the other line, so here we have to
                                // displaced this line for the whole missing distance
                                disp = (exclusion_distance + exclusion_distance_other) - old_distance;
                                disp = 0;
                            }
                            else
                            {
                                double prop1_exc =
                                    exclusion_distance / (exclusion_distance + exclusion_distance_other);
                                disp = exclusion_distance - (old_distance * prop1_exc);
                                double prop2_exc = 1 - prop1_exc;
                                auto disp2 = exclusion_distance_other - (old_distance * prop2_exc);

                                vecn<DIMENSION> new_location2 = closest_point - normalize(old_location - closest_point) * disp2;
                                builder.geometry.set_mesh_entity_vertex(
                                    other_line.gmme(), closest_vertex, new_location2,
                                    true);
                            }
                        }

                        vecn<DIMENSION> new_location = old_location + normalize(old_location - closest_point) * disp;
                        builder.geometry.set_mesh_entity_vertex(line.gmme(), v_id,
                                                                new_location, true);

                        index_t modified_vertex =
                            geomodel_.mesh.vertices.geomodel_vertex_id(line.gmme(),
                                                                       v_id);
                        modified_vertices.insert(modified_vertex);
                        auto gme_vertices =
                            geomodel_.mesh.vertices.gme_type_vertices(
                                RINGMesh::surface_type_name_static(),
                                modified_vertex);
                        for (auto gmev : gme_vertices)
                        {
                            const auto &cur_surf = geomodel_.surface(
                                gmev.gmme.index());
                            auto polygons_to_check =
                                cur_surf.mesh().polygons_around_vertex(gmev.v_index,
                                                                       false, NO_ID);
                        }
                    }
                }
            }
        }

        index_t nb_triangles_modified = 0;
        for (auto t : RINGMesh::range(geomodel_.mesh.polygons.nb()))
        {
            for (auto i : RINGMesh::range(3))
            {
                if (modified_vertices.find(
                        geomodel_.mesh.polygons.vertex({t, i})) != modified_vertices.end())
                {
                    ++nb_triangles_modified;
                    break;
                }
            }
        }
    }

    // Explicit class for DIM = 2 and DIM = 3
    template class scar_api GeoModelVertexShifter<2>;

    template class scar_api GeoModelVertexShifter<3>;
}
