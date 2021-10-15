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

#include <ringmesh/basic/geometry.h>
#include <ringmesh/mesh/mesh_base.h>
#include <ringmesh/mesh/surface_mesh.h>

#include <scar/tools/distance.h>
#include <scar/tools/utils.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModel);
    FORWARD_DECLARATION_DIMENSION_CLASS(LineMesh);
    FORWARD_DECLARATION_DIMENSION_CLASS(Box);
    FORWARD_DECLARATION_DIMENSION_CLASS(Frame);
    FORWARD_DECLARATION_DIMENSION_CLASS(ReferenceFrame);
}

namespace SCAR
{
    enum struct ConvexShapeType
    {
        NSphere,
        NEllipsoid,
        NCuboid,
        Undefined
    };

    template <index_t DIMENSION>
    class scar_api ConvexShape
    {
    public:
        ConvexShape(ConvexShapeType type)
            : type_(type)
        {
        }

        virtual ~ConvexShape() = default;

        ConvexShapeType type() const
        {
            return type_;
        }

        virtual const vecn<DIMENSION> &center() const = 0;

        virtual bool is_point_inside(const vecn<DIMENSION> &point) const = 0;

        virtual vecn<DIMENSION> farthest_point(
            const vecn<DIMENSION> &direction) const = 0;

        virtual double width(const vecn<DIMENSION> &direction) const = 0;

        virtual RINGMesh::Box<DIMENSION> aabb() const = 0;

        virtual std::tuple<bool, std::vector<vecn<DIMENSION>>> segment_intersections(
            const RINGMesh::Geometry::Segment<DIMENSION> &segment) const = 0;

        virtual void show() const = 0;

    protected:
        const ConvexShapeType type_;
    };

    ALIAS_2D_AND_3D(ConvexShape);

    template <index_t DIMENSION>
    using ShapeFactoryFromMeshVertex = RINGMesh::Factory<ConvexShapeType, ConvexShape<DIMENSION>, const RINGMesh::MeshBase<DIMENSION> &, index_t>;
    ALIAS_2D_AND_3D(ShapeFactoryFromMeshVertex);

    template <index_t DIMENSION>
    using ShapeFactoryFromInterpolation = RINGMesh::Factory<ConvexShapeType,
                                                            ConvexShape<DIMENSION>,
                                                            const RINGMesh::MeshBase<DIMENSION> &,
                                                            const std::vector<index_t> &,
                                                            const std::vector<double> &>;
    ALIAS_2D_AND_3D(ShapeFactoryFromInterpolation);

    template <index_t DIMENSION>
    std::unique_ptr<ConvexShape<DIMENSION>> get_exclusion_shape(
        const ConvexShapeType shape_type,
        const RINGMesh::MeshBase<DIMENSION> &mesh,
        const index_t vertex_id)
    {
        auto exclusion_shape = ShapeFactoryFromMeshVertex<DIMENSION>::create(
            shape_type, mesh, vertex_id);
        if (exclusion_shape == nullptr)
        {
            throw SCARException("Shape", "Exclusion shape can not be created");
        }
        return exclusion_shape;
    }

    // Interpolation on Line edge
    template <index_t DIMENSION>
    std::unique_ptr<ConvexShape<DIMENSION>> get_exclusion_shape(
        const ConvexShapeType shape_type,
        const RINGMesh::MeshBase<DIMENSION> &mesh,
        const LineRelativeCoordinates line_rcoord)
    {
        auto boundary_v_ids = double_integer_bounds(line_rcoord);
        std::vector<index_t> edge_vertex_ids{boundary_v_ids.first,
                                             boundary_v_ids.second};
        std::vector<double> linear_coefs{1 - line_rcoord.norm, line_rcoord.norm};
        auto exclusion_shape = ShapeFactoryFromInterpolation<DIMENSION>::create(
            shape_type, mesh, edge_vertex_ids, linear_coefs);

        if (exclusion_shape == nullptr)
        {
            throw SCARException("Shape", "Exclusion shape can not be created");
        }
        return exclusion_shape;
    }

    // Interpolation on Surface triangle
    std::unique_ptr<ConvexShape3D> get_exclusion_shape(
        const ConvexShapeType shape_type,
        const RINGMesh::SurfaceMesh3D &mesh,
        const index_t element_id,
        const std::array<double, 3> bcoords);

    template <index_t DIMENSION>
    class scar_api NSphere : public ConvexShape<DIMENSION>
    {
    public:
        NSphere(const vecn<DIMENSION> &center, double radius);

        NSphere(const RINGMesh::MeshBase<DIMENSION> &mesh, index_t vertex_id);

        NSphere(
            const RINGMesh::MeshBase<DIMENSION> &mesh,
            const std::vector<index_t> &element_v_id,
            const std::vector<double> &linear_coef);

        const vecn<DIMENSION> &center() const final
        {
            return center_;
        }

        double radius() const
        {
            return radius_;
        }

        vecn<DIMENSION> farthest_point(const vecn<DIMENSION> &direction) const final
        {
            return center_ + normalize(direction);
        }

        double width(const vecn<DIMENSION> &direction) const final
        {
            scar_unused(direction);
            return radius_;
        }

        RINGMesh::Box<DIMENSION> aabb() const final;

        bool is_point_inside(const vecn<DIMENSION> &point) const final
        {
            return (center_ - point).length2() <= std::pow(radius_, 2);
        }

        std::tuple<bool, std::vector<vecn<DIMENSION>>> segment_intersections(
            const RINGMesh::Geometry::Segment<DIMENSION> &segment) const final;

        void show() const final
        {
            Logger::out("Convex", "Radius = ", radius_);
        }

    private:
        vecn<DIMENSION> center_;
        double radius_;
    };

    using Circle = NSphere<2>;
    using Sphere = NSphere<3>;

    template <index_t DIMENSION>
    vecn<DIMENSION> get_local_coordinates(
        const RINGMesh::ReferenceFrame<DIMENSION> &ref_frame,
        const vecn<DIMENSION> &space_point);

    template <index_t DIMENSION>
    vecn<DIMENSION> get_space_coordinates(
        const RINGMesh::ReferenceFrame<DIMENSION> &ref_frame,
        const vecn<DIMENSION> &local_point);

    template <index_t DIMENSION>
    RINGMesh::ReferenceFrame<DIMENSION> normalize_frame(
        const RINGMesh::ReferenceFrame<DIMENSION> &ref_frame);

    template <index_t DIMENSION>
    vecn<DIMENSION> vec_multiplication(
        const vecn<DIMENSION> &vec1,
        const vecn<DIMENSION> &vec2)
    {
        vecn<DIMENSION> result;
        for (auto c : RINGMesh::range(DIMENSION))
        {
            result[c] = vec1[c] * vec2[c];
        }
        return result;
    }

    template <index_t DIMENSION>
    class scar_api NEllipsoid : public ConvexShape<DIMENSION>
    {
    public:
        NEllipsoid(
            const vecn<DIMENSION> &center,
            const RINGMesh::Frame<DIMENSION> &ppl_axis);

        NEllipsoid(
            const RINGMesh::MeshBase<DIMENSION> &mesh,
            index_t element_id);

        NEllipsoid(
            const RINGMesh::MeshBase<DIMENSION> &mesh,
            const std::vector<index_t> &element_v_id,
            const std::vector<double> &linear_coef);

        const vecn<DIMENSION> &center() const final
        {
            return center_;
        }

        vecn<DIMENSION> principal_axis_lengthes() const
        {
            vecn<DIMENSION> lengthes;
            for (auto c : RINGMesh::range(DIMENSION))
            {
                lengthes[c] = (ppl_axis_[c]).length();
            }
            return lengthes;
        }

        vecn<DIMENSION> farthest_point(const vecn<DIMENSION> &direction) const final
        {
            RINGMesh::ReferenceFrame<DIMENSION> ellipse_ref_frame(center_,
                                                                  ppl_axis_);
            // Rotation
            vecn<DIMENSION> transformed_direction = get_local_coordinates(
                ellipse_ref_frame, center_ + direction);
            transformed_direction = vec_multiplication(principal_axis_lengthes(),
                                                       transformed_direction);
            transformed_direction = normalize(transformed_direction);
            transformed_direction = vec_multiplication(principal_axis_lengthes(),
                                                       transformed_direction);
            // Rotation in the other sense + translation
            return get_space_coordinates(normalize_frame(ellipse_ref_frame),
                                         transformed_direction);
        }

        double width(const vecn<DIMENSION> &direction) const final
        {
            auto extremal_point = farthest_point(direction);
            return (extremal_point - center_).length2();
        }

        RINGMesh::Box<DIMENSION> aabb() const final;

        bool is_point_inside(const vecn<DIMENSION> &point) const final
        {
            vecn<DIMENSION> local_coord = get_local_coordinates(
                RINGMesh::ReferenceFrame<DIMENSION>(center_, ppl_axis_), point);
            return local_coord.length2() <= 1.;
        }

        std::tuple<bool, std::vector<vecn<DIMENSION>>> segment_intersections(
            const RINGMesh::Geometry::Segment<DIMENSION> &segment) const final;

        void show() const final
        {
            Logger::out("Convex", "Not implemented yet");
        }

    private:
        vecn<DIMENSION> center_;
        RINGMesh::Frame<DIMENSION> ppl_axis_;
    };

    using Ellipse = NEllipsoid<2>;
    using Ellipsoid = NEllipsoid<3>;

    void scar_api register_forms();

    template <index_t DIMENSION>
    bool are_convex_intersecting(
        const ConvexShape<DIMENSION> &convex1,
        const ConvexShape<DIMENSION> &convex2);

    template <index_t DIMENSION>
    vecn<DIMENSION> penetration_vector(
        const ConvexShape<DIMENSION> &convex1,
        const ConvexShape<DIMENSION> &convex2);
}
