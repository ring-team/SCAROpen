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

#include <scar/tools/triangle_quality.h>

#include <geogram/basic/attributes.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/mesh_aabb.h>

#include <scar/tools/utils.h>

namespace
{
    using namespace SCAR;
    /*!
     * @param[in] mesh_qual_mode mesh quality number.
     * @return the property name associated to the mesh quality number
     * \p mesh_qual_mode.
     */
    std::string mesh_qual_mode_to_prop_name(TriangleQualityMode mesh_qual_mode)
    {
        std::string quality_name;
        switch (mesh_qual_mode)
        {
        case TriangleQualityMode::ALL:
            quality_name = "ALL";
            break;
        case TriangleQualityMode::MIN_ANGLE:
            quality_name = "MIN_ANGLE";
            break;
        case TriangleQualityMode::MAX_ANGLE:
            quality_name = "MAX_ANGLE";
            break;
        case TriangleQualityMode::ANGLE_DIFFERENCE:
            quality_name = "ANGLE_DIFFERENCE";
            break;
        case TriangleQualityMode::RADIUS_RATIO:
            quality_name = "RADIUS_RATIO";
            break;
        case TriangleQualityMode::EDGE_RATIO:
            quality_name = "EDGE_RATIO";
            break;
        case TriangleQualityMode::EDGE_TO_CIRCUMRADIUS:
            quality_name = "EDGE_TO_CIRCUMRADIUS";
            break;
        case TriangleQualityMode::EDGE_TO_INRADIUS:
            quality_name = "EDGE_TO_INRADIUS";
            break;
        case TriangleQualityMode::MIN_HEIGHT:
            quality_name = "MIN_HEIGHT";
            break;

        default:
            scar_assert_not_reached;
        }
        scar_assert(!quality_name.empty());
        return quality_name;
    }

    template <index_t DIMENSION>
    void triangle_angles(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2,
        double &alpha,
        double &beta,
        double &gamma)
    {
        // Side lengths
        double e0 = (v2 - v1).length();
        double e1 = (v2 - v0).length();
        double e2 = (v1 - v0).length();

        // From Cosine law
        alpha = std::acos((e1 * e1 + e2 * e2 - e0 * e0) / (2 * e1 * e2));
        beta = std::acos((e0 * e0 + e2 * e2 - e1 * e1) / (2 * e0 * e2));
        gamma = std::acos((e0 * e0 + e1 * e1 - e2 * e2) / (2 * e0 * e1));

        // Converting to degree
        alpha = alpha * 180 / M_PI;
        beta = beta * 180 / M_PI;
        gamma = gamma * 180 / M_PI;
        return;
    }

    template <index_t DIMENSION>
    double triangle_quality_min_angle(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2)
    {
        double alpha, beta, gamma = -1;
        triangle_angles(v0, v1, v2, alpha, beta, gamma);
        return std::min(std::min(alpha, beta), gamma);
    }

    template <index_t DIMENSION>
    double triangle_quality_max_angle(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2)
    {
        double alpha, beta, gamma = -1;
        triangle_angles(v0, v1, v2, alpha, beta, gamma);
        return std::max(std::max(alpha, beta), gamma);
    }

    template <index_t DIMENSION>
    double triangle_quality_angle_difference(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2)
    {
        double alpha, beta, gamma = -1;
        triangle_angles(v0, v1, v2, alpha, beta, gamma);
        double angle_min = std::min(std::min(alpha, beta), gamma);
        double angle_max = std::max(std::max(alpha, beta), gamma);
        return angle_max - angle_min;
    }

    template <index_t DIMENSION>
    double triangle_quality_radius_ratio(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2)
    {
        double e0 = (v2 - v1).length();
        double e1 = (v2 - v0).length();
        double e2 = (v1 - v0).length();
        double half_perim = (e0 + e1 + e2) / 2.;
        double area = std::sqrt(
            half_perim * (half_perim - e0) * (half_perim - e1) * (half_perim - e2));
        double incircle_radius = area / half_perim;
        double circum_radius = e0 * e1 * e2 / (4 * area);
        scar_assert(circum_radius > global_epsilon);
        return 2 * incircle_radius / circum_radius;
    }

    template <index_t DIMENSION>
    double triangle_quality_edge_ratio(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2)
    {
        // Side lengths
        double e0 = (v2 - v1).length();
        double e1 = (v2 - v0).length();
        double e2 = (v1 - v0).length();
        double edge_min = std::min(std::min(e0, e1), e2);
        double edge_max = std::max(std::max(e0, e1), e2);
        scar_assert(edge_min > global_epsilon) return edge_max / edge_min;
    }

    template <index_t DIMENSION>
    double triangle_quality_edge_to_circumradius(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2)
    {
        double e0 = (v2 - v1).length();
        double e1 = (v2 - v0).length();
        double e2 = (v1 - v0).length();
        double half_perim = (e0 + e1 + e2) / 2.;
        double area = std::sqrt(
            half_perim * (half_perim - e0) * (half_perim - e1) * (half_perim - e2));
        double circum_radius = e0 * e1 * e2 / (4 * area);
        return circum_radius / half_perim;
    }

    template <index_t DIMENSION>
    double triangle_quality_edge_to_inradius(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2)
    {
        double e0 = (v2 - v1).length();
        double e1 = (v2 - v0).length();
        double e2 = (v1 - v0).length();
        double half_perim = (e0 + e1 + e2) / 2.;
        double area = std::sqrt(
            half_perim * (half_perim - e0) * (half_perim - e1) * (half_perim - e2));
        double incircle_radius = area / half_perim;
        return incircle_radius / half_perim;
    }

    template <index_t DIMENSION>
    double triangle_quality_min_height(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2)
    {
        // Side lengths
        double e0 = (v2 - v1).length();
        double e1 = (v2 - v0).length();
        double e2 = (v1 - v0).length();
        double half_perim = (e0 + e1 + e2) / 2.;
        double area = std::sqrt(
            half_perim * (half_perim - e0) * (half_perim - e1) * (half_perim - e2));
        auto max_edge_length = std::max(std::max(e0, e1), e2);
        scar_assert(max_edge_length > global_epsilon) return 2 * area / max_edge_length;
    }

    template <index_t DIMENSION>
    double get_triangle_quality(
        const vecn<DIMENSION> &v0,
        const vecn<DIMENSION> &v1,
        const vecn<DIMENSION> &v2,
        TriangleQualityMode mesh_qual_mode)
    {
        double quality = -1;
        switch (mesh_qual_mode)
        {
        case TriangleQualityMode::MIN_ANGLE:
            quality = triangle_quality_min_angle(v0, v1, v2);
            break;
        case TriangleQualityMode::MAX_ANGLE:
            quality = triangle_quality_max_angle(v0, v1, v2);
            break;
        case TriangleQualityMode::ANGLE_DIFFERENCE:
            quality = triangle_quality_angle_difference(v0, v1, v2);
            break;
        case TriangleQualityMode::RADIUS_RATIO:
            quality = triangle_quality_radius_ratio(v0, v1, v2);
            break;
        case TriangleQualityMode::EDGE_RATIO:
            quality = triangle_quality_edge_ratio(v0, v1, v2);
            break;
        case TriangleQualityMode::EDGE_TO_CIRCUMRADIUS:
            quality = triangle_quality_edge_to_circumradius(v0, v1, v2);
            break;
        case TriangleQualityMode::EDGE_TO_INRADIUS:
            quality = triangle_quality_edge_to_inradius(v0, v1, v2);
            break;
        case TriangleQualityMode::MIN_HEIGHT:
            quality = triangle_quality_min_height(v0, v1, v2);
            break;
        default:
            scar_assert_not_reached;
        }
        return quality;
    }
}

namespace SCAR
{
    template <index_t DIMENSION>
    void compute_prop_triangle_mesh_quality(
        TriangleQualityMode mesh_qual_mode,
        const RINGMesh::GeoModel<DIMENSION> &geomodel)
    {
        scar_assert(geomodel.nb_surfaces() != 0);
        for (const auto &surface : geomodel.surfaces())
        {
            if (!surface.is_simplicial())
            {
                throw SCARException("Quality", "Surfaces should be triangulated.");
            }
            GEO::AttributesManager &surface_attr_mgr =
                surface.polygon_attribute_manager();
            std::string mode = mesh_qual_mode_to_prop_name(mesh_qual_mode);
            if (mode == "ALL")
            {
                for (auto m : RINGMesh::range(1, 9))
                {
                    GEO::Attribute<double> attr(surface_attr_mgr,
                                                mesh_qual_mode_to_prop_name(
                                                    static_cast<TriangleQualityMode>(m)));
                    for (auto polygon_itr : RINGMesh::range(
                             surface.nb_mesh_elements()))
                    {
                        attr[polygon_itr] = get_triangle_quality(
                            surface.mesh_element_vertex({polygon_itr, 0}),
                            surface.mesh_element_vertex({polygon_itr, 1}),
                            surface.mesh_element_vertex({polygon_itr, 2}),
                            static_cast<TriangleQualityMode>(m));
                    }
                }
            }
            else
            {
                index_t nb_bad = 0;
                double min_value = max_double();
                GEO::Attribute<double> attr(surface_attr_mgr,
                                            mesh_qual_mode_to_prop_name(mesh_qual_mode));
                for (auto polygon_itr : RINGMesh::range(surface.nb_mesh_elements()))
                {
                    auto value = get_triangle_quality(surface.mesh_element_vertex({polygon_itr, 0}), surface.mesh_element_vertex({polygon_itr, 1}), surface.mesh_element_vertex({polygon_itr, 2}), mesh_qual_mode);
                    attr[polygon_itr] = value;
                    min_value = std::min(min_value, value);
                    if (value < 0.2)
                    {
                        ++nb_bad;
                    }
                }
                Logger::out("MinValue", surface.gmme());
                Logger::out("MinValue", mode, " - min value = ", min_value);
                Logger::out("MinValue", mode, " - % bad = ",
                            ((double)nb_bad / (double)surface.nb_mesh_elements()));
            }
        }
        return;
    }

    template <index_t DIMENSION>
    void output_quality_values(
        const std::string &filename,
        TriangleQualityMode mesh_qual_mode,
        const RINGMesh::GeoModel<DIMENSION> &geomodel)
    {

        scar_assert(geomodel.nb_surfaces() != 0);

        std::string mode = mesh_qual_mode_to_prop_name(mesh_qual_mode);
        if (mode == "ALL")
        {
            for (auto m : RINGMesh::range(1, 8))
            {
                std::ofstream out;
                std::string cur_mode = mesh_qual_mode_to_prop_name(
                    static_cast<TriangleQualityMode>(m));
                std::string cur_filename = filename + cur_mode;
                out.open(cur_filename);

                for (const auto &surface : geomodel.surfaces())
                {
                    if (!surface.is_simplicial())
                    {
                        throw SCARException("Quality",
                                            "Surfaces should be triangulated.");
                    }
                    GEO::AttributesManager &surface_attr_mgr =
                        surface.polygon_attribute_manager();
                    GEO::Attribute<double> attr(surface_attr_mgr, cur_mode);
                    for (auto polygon_itr : RINGMesh::range(
                             surface.nb_mesh_elements()))
                    {
                        out << attr[polygon_itr] << "\n";
                    }
                }
                out << std::endl;
                out.close();
            }
        }
        else
        {
            std::ofstream out;
            out.open(filename);

            for (const auto &surface : geomodel.surfaces())
            {
                if (!surface.is_simplicial())
                {
                    throw SCARException("Quality",
                                        "Surfaces should be triangulated.");
                }
                GEO::AttributesManager &surface_attr_mgr =
                    surface.polygon_attribute_manager();
                GEO::Attribute<double> attr(surface_attr_mgr,
                                            mesh_qual_mode_to_prop_name(mesh_qual_mode));
                for (auto polygon_itr : RINGMesh::range(surface.nb_mesh_elements()))
                {
                    out << attr[polygon_itr] << "\n";
                }
            }
            out << std::endl;
        }
        return;
    }

    template void scar_api compute_prop_triangle_mesh_quality(
        TriangleQualityMode,
        const RINGMesh::GeoModel2D &);

    template void scar_api output_quality_values(
        const std::string &,
        TriangleQualityMode,
        const RINGMesh::GeoModel2D &);

    template void scar_api compute_prop_triangle_mesh_quality(
        TriangleQualityMode,
        const RINGMesh::GeoModel3D &);

    template void scar_api output_quality_values(
        const std::string &filename,
        TriangleQualityMode mesh_qual_mode,
        const RINGMesh::GeoModel3D &geomodel);
}
