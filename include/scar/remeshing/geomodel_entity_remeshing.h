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

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModel);
    ALIAS_2D_AND_3D(GeoModel);
}

namespace SCAR
{
    FORWARD_DECLARATION_DIMENSION_CLASS(GeoModelTopologyMaker);
    ALIAS_2D_AND_3D(GeoModelTopologyMaker);
}

namespace SCAR
{

#ifdef SCAR_WITH_MMG
    void scar_api mesh_surfaces_using_mmg2d(
        RINGMesh::GeoModel2D &geomodel,
        const std::vector<double> &surface_resolutions);
#endif

    template <index_t DIMENSION>
    void remove_corners_with_valence_2(RINGMesh::GeoModel<DIMENSION> &geomodel);

    template <index_t DIMENSION>
    class scar_api GeoModelRemesherBase
    {
        scar_disable_copy(GeoModelRemesherBase);

    public:
        virtual ~GeoModelRemesherBase() = default;

        virtual void remesh_geomodel_entities();

    protected:
        GeoModelRemesherBase(
            RINGMesh::GeoModel<DIMENSION> &geomodel,
            const RINGMesh::GeoModel<DIMENSION> &init_geomodel,
            const GeoModelTopologyMaker<DIMENSION> &topology_maker,
            bool verbose);

    protected:
        RINGMesh::GeoModel<DIMENSION> &geomodel_;
        const RINGMesh::GeoModel<DIMENSION> &init_geomodel_;
        const GeoModelTopologyMaker<DIMENSION> &topology_maker_;
        bool verbose_;
    };

    template <index_t DIMENSION>
    class scar_api GeoModelRemesher final : public GeoModelRemesherBase<DIMENSION>
    {
    public:
        GeoModelRemesher(
            RINGMesh::GeoModel<DIMENSION> &geomodel,
            const RINGMesh::GeoModel<DIMENSION> &init_geomodel,
            const GeoModelTopologyMaker<DIMENSION> &topology_maker,
            bool verbose);
    };

    template <>
    class scar_api GeoModelRemesher<2> final : public GeoModelRemesherBase<2>
    {
    public:
        GeoModelRemesher(
            RINGMesh::GeoModel2D &geomodel,
            const RINGMesh::GeoModel2D &init_geomodel,
            const GeoModelTopologyMaker2D &topology_maker,
            bool verbose);
        ~GeoModelRemesher() override = default;
        void remesh_geomodel_entities() override;
    };

    template <>
    class scar_api GeoModelRemesher<3> final : public GeoModelRemesherBase<3>
    {
    public:
        GeoModelRemesher(
            RINGMesh::GeoModel3D &geomodel,
            const RINGMesh::GeoModel3D &init_geomodel,
            const GeoModelTopologyMaker3D &topology_maker,
            bool verbose);
        ~GeoModelRemesher() override = default;
        void remesh_geomodel_entities() override;
    };

}
