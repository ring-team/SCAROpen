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

#include <scar/tools/convex_shape.h>

#include <geogram/basic/attributes.h>

#include <ringmesh/basic/box.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/mesh/mesh_index.h>

namespace SCAR
{

    template <index_t DIMENSION>
    double segment_point_normalized_coordinates(
        const RINGMesh::Geometry::Segment<DIMENSION> &oriented_segment,
        const vecn<DIMENSION> &point)
    {
        if (oriented_segment.length() < global_epsilon && (point - oriented_segment.p0).length() < global_epsilon)
        {
            return 0;
        }
        return std::sqrt(
            (point - oriented_segment.p0).length2() / (oriented_segment.p1 - oriented_segment.p0).length2());
    }

    template <index_t DIMENSION>
    bool intersection_along_center_segment(
        const std::unique_ptr<ConvexShape<DIMENSION>> &shape1,
        const std::unique_ptr<ConvexShape<DIMENSION>> &shape2)
    {
        const RINGMesh::Geometry::Segment<DIMENSION> center_segment(
            shape1->center(), shape2->center());

        double shape1_width = shape1->width(-center_segment.direction());
        double shape2_width = shape2->width(center_segment.direction());

        return shape1_width + shape2_width >= center_segment.length();
    }

    template <index_t DIMENSION>
    bool aabb_intersection_barycenter_in_both_shapes(
        const std::unique_ptr<ConvexShape<DIMENSION>> &shape1,
        const std::unique_ptr<ConvexShape<DIMENSION>> &shape2,
        const vecn<DIMENSION> shape_aabb_intersection)
    {
        return shape1->is_point_inside(shape_aabb_intersection) && shape2->is_point_inside(shape_aabb_intersection);
    }

    template <index_t DIMENSION>
    bool are_spheres_intersecting(
        const NSphere<DIMENSION> &sphere1,
        const NSphere<DIMENSION> &sphere2)
    {
        return (sphere1.center() - sphere2.center()).length() < (sphere1.radius() + sphere2.radius());
    }

    template <index_t DIMENSION>
    bool gjk_convexes_intersection(
        const ConvexShape<DIMENSION> &convex1,
        const ConvexShape<DIMENSION> &convex2)
    {
        scar_unused(convex1);
        scar_unused(convex2);
        throw SCARException("GJK", "TODO: implement GJK algorithm..");
        return false;
    }

    template <index_t DIMENSION>
    bool are_convex_intersecting(
        const ConvexShape<DIMENSION> &convex1,
        const ConvexShape<DIMENSION> &convex2)
    {
        if (convex1.type() == ConvexShapeType::NSphere && convex2.type() == ConvexShapeType::NSphere)
        {
            // Use fastest intersection algorithm for circles
            vecn<DIMENSION> arbitrary_direction;
            arbitrary_direction[0] = 1.;
            return are_spheres_intersecting(
                NSphere<DIMENSION>(convex1.center(),
                                   convex1.width(arbitrary_direction)),
                NSphere<DIMENSION>(convex2.center(),
                                   convex2.width(arbitrary_direction)));
        }
        return gjk_convexes_intersection(convex1, convex2);
    }

    template <index_t DIMENSION>
    vecn<DIMENSION> epa_convexes_penetration(
        const ConvexShape<DIMENSION> &convex1,
        const ConvexShape<DIMENSION> &convex2)
    {
        scar_unused(convex1);
        scar_unused(convex2);
        throw SCARException("GJK", "TODO: implement GJK-EPA algorithm..");
        return vecn<DIMENSION>();
    }

    template <index_t DIMENSION>
    vecn<DIMENSION> penetration_vector(
        const ConvexShape<DIMENSION> &convex1,
        const ConvexShape<DIMENSION> &convex2)
    {
        if (convex1.type() == ConvexShapeType::NSphere && convex2.type() == ConvexShapeType::NSphere)
        {
            // Use fastest computation for circles
            vecn<DIMENSION> direction = normalize(
                convex1.center() - convex2.center());
            double penetration_dist =
                (convex1.center() - convex2.center()).length() < (convex1.width(direction) + convex2.width(direction));

            return penetration_dist * direction;
        }
        return epa_convexes_penetration(convex1, convex2);
    }

    template <index_t DIMENSION>
    NSphere<DIMENSION>::NSphere(const vecn<DIMENSION> &center, double radius)
        : ConvexShape<DIMENSION>(ConvexShapeType::NSphere),
          center_(center),
          radius_(radius)
    {
    }

    template <index_t DIMENSION>
    NSphere<DIMENSION>::NSphere(
        const RINGMesh::MeshBase<DIMENSION> &mesh,
        index_t element_id)
        : ConvexShape<DIMENSION>(ConvexShapeType::NSphere)
    {
        center_ = mesh.vertex(element_id);
        GEO::Attribute<double> radius_attribute(mesh.vertex_attribute_manager(),
                                                tolerance_att_name);
        radius_ = radius_attribute[element_id];
    }

    template <index_t DIMENSION>
    NSphere<DIMENSION>::NSphere(
        const RINGMesh::MeshBase<DIMENSION> &mesh,
        const std::vector<index_t> &element_v_id,
        const std::vector<double> &linear_coef)
        : ConvexShape<DIMENSION>(ConvexShapeType::NSphere)
    {
        scar_assert(linear_coef.size() == element_v_id.size());
        radius_ = 0;
        GEO::Attribute<double> radius_attribute(mesh.vertex_attribute_manager(),
                                                tolerance_att_name);
        for (auto i : RINGMesh::range(element_v_id.size()))
        {
            center_ += linear_coef[i] * mesh.vertex(element_v_id[i]);
            radius_ += linear_coef[i] * radius_attribute[element_v_id[i]];
        }
    }

    template <index_t DIMENSION>
    std::tuple<bool, std::vector<vecn<DIMENSION>>> NSphere<DIMENSION>::segment_intersections(
        const RINGMesh::Geometry::Segment<DIMENSION> &segment) const
    {
        return RINGMesh::Intersection::segment_sphere(segment,
                                                      RINGMesh::Geometry::Sphere<DIMENSION>{this->center_, this->radius_});
    }

    template <index_t DIMENSION>
    RINGMesh::Box<DIMENSION> NSphere<DIMENSION>::aabb() const
    {
        RINGMesh::Box<DIMENSION> box;

        box.add_point(
            center_ + RINGMesh::initialize_vecn_coordinates<DIMENSION>(radius_));
        box.add_point(
            center_ - RINGMesh::initialize_vecn_coordinates<DIMENSION>(radius_));
        return box;
    }

    template <index_t DIMENSION>
    NEllipsoid<DIMENSION>::NEllipsoid(
        const vecn<DIMENSION> &center,
        const RINGMesh::Frame<DIMENSION> &ppl_axis)
        : ConvexShape<DIMENSION>(ConvexShapeType::NEllipsoid),
          center_(center),
          ppl_axis_(ppl_axis)
    {
    }

    template <index_t DIMENSION>
    NEllipsoid<DIMENSION>::NEllipsoid(
        const RINGMesh::MeshBase<DIMENSION> &mesh,
        index_t element_id)
        : ConvexShape<DIMENSION>(ConvexShapeType::NSphere)
    {
        center_ = mesh.vertex(element_id);
        GEO::Attribute<double> ellipse_major_ppl_attribute(
            mesh.vertex_attribute_manager(), ellipse_major_ppl_att_name);
        GEO::Attribute<double> ellipse_minor_ppl_attribute(
            mesh.vertex_attribute_manager(), ellipse_minor_ppl_att_name);
        // @todo frame is by default oriented on axis (read in the mesh orientation)
        ppl_axis_[0] *= ellipse_major_ppl_attribute[element_id];
        ppl_axis_[1] *= ellipse_minor_ppl_attribute[element_id];
        ppl_axis_[2] *= ellipse_minor_ppl_attribute[element_id];
        scar_assert_not_reached;
    }

    template <index_t DIMENSION>
    NEllipsoid<DIMENSION>::NEllipsoid(
        const RINGMesh::MeshBase<DIMENSION> &mesh,
        const std::vector<index_t> &element_v_id,
        const std::vector<double> &linear_coef)
        : ConvexShape<DIMENSION>(ConvexShapeType::NEllipsoid)
    {
        scar_unused(mesh);
        scar_unused(element_v_id);
        scar_unused(linear_coef);
        scar_assert(linear_coef.size() == element_v_id.size());
        Logger::err("NEllipsoid",
                    "Not able yet to build a NEllipsoid by linear interpolation.");
        scar_assert_not_reached;
    }

    template <index_t DIMENSION>
    std::tuple<bool, std::vector<vecn<DIMENSION>>> NEllipsoid<DIMENSION>::segment_intersections(
        const RINGMesh::Geometry::Segment<DIMENSION> &segment) const
    {
        scar_assert_not_reached;
        scar_unused(segment);
        return std::make_tuple(false, std::vector<vecn<DIMENSION>>());
    }

    template <index_t DIMENSION>
    RINGMesh::Box<DIMENSION> NEllipsoid<DIMENSION>::aabb() const
    {
        RINGMesh::Box<DIMENSION> box;

        for (auto coord : RINGMesh::range(DIMENSION))
        {
            scar_unused(coord);
            vecn<DIMENSION> direction;
            direction[DIMENSION] = 1.;
            box.add_point(this->farthest_point(direction));
            box.add_point(this->farthest_point(-direction));
        }

        return box;
    }

    std::unique_ptr<ConvexShape3D> get_exclusion_shape(
        const ConvexShapeType shape_type,
        const RINGMesh::SurfaceMesh3D &mesh,
        const index_t element_id,
        const std::array<double, 3> bcoords)
    {
        index_t v0_id, v1_id, v2_id;
        v0_id = mesh.polygon_vertex(RINGMesh::ElementLocalVertex(element_id, 0));
        v1_id = mesh.polygon_vertex(RINGMesh::ElementLocalVertex(element_id, 1));
        v2_id = mesh.polygon_vertex(RINGMesh::ElementLocalVertex(element_id, 2));
        std::vector<index_t> triangle_vertex_ids{v0_id, v1_id, v2_id};
        std::vector<double> linear_coefs{bcoords[0], bcoords[1], bcoords[2]};
        auto exclusion_shape = ShapeFactoryFromInterpolation3D::create(shape_type,
                                                                       mesh, triangle_vertex_ids, linear_coefs);

        if (exclusion_shape == nullptr)
        {
            throw SCARException("Shape", "Exclusion shape can not be created");
        }
        return exclusion_shape;
    }

    template <index_t DIMENSION>
    vecn<DIMENSION> get_local_coordinates(
        const RINGMesh::ReferenceFrame<DIMENSION> &ref_frame,
        const vecn<DIMENSION> &space_point)
    {
        vecn<DIMENSION> local_coords;
        for (auto c : RINGMesh::range(DIMENSION))
        {
            local_coords[c] = dot((space_point - ref_frame.origin()),
                                  ref_frame[c]) /
                              (ref_frame[c]).length();
        }
        return local_coords;
    }

    template <index_t DIMENSION>
    vecn<DIMENSION> get_space_coordinates(
        const RINGMesh::ReferenceFrame<DIMENSION> &ref_frame,
        const vecn<DIMENSION> &local_point)
    {
        vecn<DIMENSION> space_coords = ref_frame.origin();
        for (auto c : RINGMesh::range(DIMENSION))
        {
            space_coords += local_point[c] * ref_frame[c];
        }
        return space_coords;
    }

    template <index_t DIMENSION>
    RINGMesh::ReferenceFrame<DIMENSION> normalize_frame(
        const RINGMesh::ReferenceFrame<DIMENSION> &ref_frame)
    {
        RINGMesh::ReferenceFrame<DIMENSION> normalized_frame(ref_frame);
        for (auto c : RINGMesh::range(DIMENSION))
        {
            normalized_frame[c] /= (normalized_frame[c]).length();
        }
        return normalized_frame;
    }

    void scar_api register_forms()
    {
        ShapeFactoryFromMeshVertex2D::register_creator<Circle>(
            ConvexShapeType::NSphere);
        ShapeFactoryFromMeshVertex2D::register_creator<Ellipse>(
            ConvexShapeType::NEllipsoid);

        ShapeFactoryFromMeshVertex3D::register_creator<Sphere>(
            ConvexShapeType::NSphere);
        ShapeFactoryFromMeshVertex3D::register_creator<Ellipsoid>(
            ConvexShapeType::NEllipsoid);

        ShapeFactoryFromInterpolation2D::register_creator<Circle>(
            ConvexShapeType::NSphere);
        ShapeFactoryFromInterpolation2D::register_creator<Ellipse>(
            ConvexShapeType::NEllipsoid);

        ShapeFactoryFromInterpolation3D::register_creator<Sphere>(
            ConvexShapeType::NSphere);
        ShapeFactoryFromInterpolation3D::register_creator<Ellipsoid>(
            ConvexShapeType::NEllipsoid);
    }

    template class scar_api NSphere<2>;
    template class scar_api NEllipsoid<2>;

    template bool scar_api are_convex_intersecting(
        const ConvexShape2D &,
        const ConvexShape2D &);

    template vec2 scar_api penetration_vector(
        const ConvexShape2D &convex1,
        const ConvexShape2D &convex2);

    template class scar_api NSphere<3>;
    template class scar_api NEllipsoid<3>;

    template bool scar_api are_convex_intersecting(
        const ConvexShape3D &,
        const ConvexShape3D &);

    template vec3 scar_api penetration_vector(
        const ConvexShape3D &convex1,
        const ConvexShape3D &convex2);
}
