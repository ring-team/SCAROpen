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

#include <scar/io/geomodel_builder_msh.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/mesh/mesh_index.h>

namespace SCAR
{

    ////////////////////////////////////////////////
    // This code is adapted from the RINGMesh one //
    ////////////////////////////////////////////////

    const char EOL = '\n';
    const char SPACE = ' ';

    // Vertices order for a vertex in GMSH
    static index_t vertices_in_vertex[1] = {0};
    // Vertices order for an edge in GMSH
    static index_t vertices_in_edge[2] = {0, 1};
    // Vertices order for a triangle in GMSH
    static index_t vertices_in_triangle[3] = {0, 1, 2};
    // Vertices order for a quad in GMSH
    static index_t vertices_in_quad[4] = {0, 1, 2, 3};
    // Vertices order for a tetrahedron in GMSH
    static index_t vertices_in_tetrahedron[4] = {0, 1, 2, 3};
    // Vertices order for an hexahedron in GMSH
    static index_t vertices_in_hexahedron[8] = {4, 0, 5, 1, 7, 3, 6, 2};
    // Vertices order for a prism in GMSH
    static index_t vertices_in_prism[6] = {0, 1, 2, 3, 4, 5};
    // Vertices order for a pyramid in GMSH
    static index_t vertices_in_pyramid[5] = {0, 1, 2, 3, 4};

    // GMSH count begin at 1
    index_t gmsh_offset = 1;

    // This is a tricky table that associate an unique id with the id of
    // elements
    // in GMSH. The unique id is computed as the sum of the MeshEntityType index
    // and the number of vertices of the elements
    // Vertex :      0 + 1 = 1
    // Edge :        1 + 2 = 3
    // Triangle :    2 + 3 = 5
    // Quad :        2 + 4 = 6
    // Tetrahedron : 3 + 4 = 7
    // Pyramid :     3 + 5 = 8
    // Prism         3 + 6 = 9
    // Hexahedron :  3 + 8 = 11
    index_t element_type[12] =
        {NO_ID, 15, NO_ID, 1, NO_ID, 2, 3, 4, 7, 6, NO_ID, 5};

    // This is a tricky table that associate an unique id with another table
    // containing the ordering of vertices inside the elementS
    index_t *vertices_in_elements[12] = {nullptr, vertices_in_vertex, nullptr,
                                         vertices_in_edge, nullptr,
                                         vertices_in_triangle, vertices_in_quad,
                                         vertices_in_tetrahedron,
                                         vertices_in_pyramid, vertices_in_prism,
                                         nullptr, vertices_in_hexahedron};

    // The physical id is a mandatory value for GMSH which is not used
    //    index_t physical_id = 0;

    // in GMSH, a tag is a physical id (equal to the geological parent index)
    // and a geometry id
    index_t nb_of_tags = 2;

    /*!
     * @brief Count all the elements in GMSH
     * an Element can be :
     * - A Corner
     * - An Edge
     * - A Triangle
     * (- A Quad
     * - A Tetrahedron
     * - A Pyramid
     * - A Prism
     * - An Hexaedron)
     */
    index_t count_elements(const RINGMesh::GeoModel2D &geomodel)
    {
        index_t nb_elements = 0;
        const auto &gmme_types =
            geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();
        for (const auto &cur_mesh_entity_type : gmme_types)
        {
            for (auto index_of_gmme_of_the_current_type : RINGMesh::range(
                     geomodel.nb_mesh_entities(cur_mesh_entity_type)))
            {
                RINGMesh::gmme_id cur_gmme_id = RINGMesh::gmme_id(
                    cur_mesh_entity_type, index_of_gmme_of_the_current_type);
                const auto &cur_gmme = geomodel.mesh_entity(cur_gmme_id);
                nb_elements += cur_gmme.nb_mesh_elements();
            }
        }
        return nb_elements;
    }

    /*!
     * @brief Find the gmsh type on an element using
     * the number of vertices and the mesh entity index
     * in which the element belong
     */
    index_t find_gmsh_element_type(
        index_t nb_vertices,
        index_t mesh_entity_type_index)
    {
        return element_type[nb_vertices + mesh_entity_type_index];
    }

    /*!
     * @brief Find the gmsh local vertex index of an element using
     * the number of vertices and the mesh entity index
     * in which the element belong
     */
    index_t find_gmsh_element_local_vertex_id(
        index_t nb_vertices,
        index_t mesh_entity_type_index,
        index_t local_vertex_index)
    {
        return vertices_in_elements[nb_vertices + mesh_entity_type_index][local_vertex_index];
    }

    /*!
     * @brief Export for the GMSH format 2.2 which is described here:
     * http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format
     * NB : Mesh entities are also exported
     */
    void GeoModelBuilderMSH::save_geomodel(
        const RINGMesh::GeoModel2D &geomodel,
        const std::string &filename)
    {
        geomodel.mesh.vertices.clear();
        std::ofstream out(filename.c_str());
        out.precision(16);

        out << "$MeshFormat" << EOL;
        out << "2.2 0 8" << EOL;
        out << "$EndMeshFormat" << EOL;

        out << "$Nodes" << EOL;
        out << geomodel.mesh.vertices.nb() << EOL;
        for (auto v : RINGMesh::range(geomodel.mesh.vertices.nb()))
        {
            out << v + gmsh_offset << SPACE << geomodel.mesh.vertices.vertex(v)
                << SPACE << 0 << EOL;
        }
        out << "$EndNodes" << EOL;

        out << "$Elements" << EOL;
        out << count_elements(geomodel) << EOL;
        const auto &gmme_types =
            geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();

        index_t element_index = 1;
        for (auto gmme_type_index : RINGMesh::range(
                 geomodel.entity_type_manager().mesh_entity_manager.nb_mesh_entity_types()))
        {
            auto cur_mesh_entity_type = gmme_types[gmme_type_index];
            for (auto index_of_gmme_of_the_current_type : RINGMesh::range(
                     geomodel.nb_mesh_entities(cur_mesh_entity_type)))
            {
                RINGMesh::gmme_id cur_gmme_id = RINGMesh::gmme_id(
                    cur_mesh_entity_type, index_of_gmme_of_the_current_type);
                const auto &cur_gmme = geomodel.mesh_entity(cur_gmme_id);

                for (auto elem_in_cur_gmme : RINGMesh::range(
                         cur_gmme.nb_mesh_elements()))
                {
                    index_t tag_to_set =
                        cur_gmme.has_parent() ? cur_gmme.parent(0).index() + gmsh_offset : 0;
                    if (cur_gmme_id.type() == RINGMesh::line_type_name_static())
                    {
                        if (!geomodel.line(cur_gmme_id.index()).is_on_voi())
                        {
                            tag_to_set = 0;
                        }
                    }

                    index_t nb_vertices_in_cur_element =
                        cur_gmme.nb_mesh_element_vertices(elem_in_cur_gmme);
                    index_t gmsh_element_type = find_gmsh_element_type(
                        nb_vertices_in_cur_element, gmme_type_index);
                    out << element_index++ << SPACE << gmsh_element_type << SPACE
                        << nb_of_tags << SPACE << tag_to_set << SPACE << tag_to_set
                        << SPACE /*<< nb_vertices_in_cur_element << SPACE*/;
                    for (auto v_index_in_cur_element : RINGMesh::range(
                             nb_vertices_in_cur_element))
                    {
                        out
                            << geomodel.mesh.vertices.geomodel_vertex_id(
                                   cur_gmme_id,
                                   cur_gmme.mesh_element_vertex_index(
                                       RINGMesh::ElementLocalVertex(elem_in_cur_gmme,
                                                                    find_gmsh_element_local_vertex_id(
                                                                        nb_vertices_in_cur_element,
                                                                        gmme_type_index,
                                                                        v_index_in_cur_element)))) +
                                   gmsh_offset
                            << SPACE;
                    }
                    out << EOL;
                }
            }
        }
        out << "$EndElements" << EOL;
        out << std::flush;
    }

}
