# Copyright (c) 2021 ASGA and Universite de Lorraine. All Rights Reserved.
# 
# This program (SCAR) was developed in the frame of the RING project managed by 
# ASGA and Universite de Lorraine. 
# 
# It is distributed under a dual licensing scheme:
# 
# 1. RING Consortium License
# Members of the RING-COCAD Consortium may only use this file in
# accordance with the terms of described in the GOCAD Advancement Agreement, 
# without the prior written authorization of the ASGA.
# Licencee agrees to attach or embed this Notice on all copies 
# of the program, including partial copies or modified versions thereof.
# Please use: contact at ring dash team dot org, for more information. 
# 
# 2. GNU General Public License Usage
# Alternatively, this file may be used under the terms of the 
# GNU General Public license version 3. The licenses are as published by 
# the Free Software Foundation and appearing in the file 
# https://www.gnu.org/licenses/gpl-3.0.html
# 


# Cmake for scar library that create the target: SCAR::scar

# target source list:
#--------------------------------------------------
# The following function adds every .cpp and .h file to the target named: ${target_name}
#
# You must not to edit variables ${target_name} and ${lib_source_dir} because they are 
# automatically set by a macro. 
#
# You need to complete the path "${lib_source_dir}" (this directory) with the 
# path to the header and source files to be included in this library. 
# NB: subdirectories may be used to organize source files if needed, but all files 
# should be listed here.
#--------------------------------------------------

target_sources(${target_name}
    PRIVATE
	#TO EDIT: fill in the list of sources as follow:
        "${lib_source_dir}/basic/common.cpp"
        "${lib_source_dir}/basic/command_line.cpp"
        "${lib_source_dir}/io/geomodel_builder_msh.cpp"
        "${lib_source_dir}/io/geomodel_builder_pl.cpp"
        "${lib_source_dir}/io/geomodel_builder_ts.cpp"
		"${lib_source_dir}/remeshing/geomodel_entity_remeshing.cpp"
		"${lib_source_dir}/remeshing/geomodel_vertex_shifter.cpp"
        "${lib_source_dir}/repair/geomodel_correspondence_map.cpp"
        "${lib_source_dir}/repair/repairer.cpp"
        "${lib_source_dir}/repair/topology_recovery.cpp"
        "${lib_source_dir}/tools/complexity_measures.cpp"
        "${lib_source_dir}/tools/convex_shape.cpp"
        "${lib_source_dir}/tools/distance.cpp"
        "${lib_source_dir}/tools/geometry.cpp"
        "${lib_source_dir}/tools/projection.cpp"
        "${lib_source_dir}/tools/triangle_quality.cpp"
        "${lib_source_dir}/tools/utils.cpp"

    PRIVATE # Could be PUBLIC from CMake 3.3
        "${lib_include_dir}/basic/common.h"
        "${lib_include_dir}/basic/command_line.h"
        "${lib_include_dir}/basic/types.h"
        "${lib_include_dir}/io/geomodel_builder_msh.h"
        "${lib_include_dir}/io/geomodel_builder_pl.h"
        "${lib_include_dir}/io/geomodel_builder_ts.h"
        "${lib_include_dir}/io/geomodel_io.h"
		"${lib_include_dir}/remeshing/geomodel_entity_remeshing.h"
		"${lib_include_dir}/remeshing/geomodel_vertex_shifter.h"
        "${lib_include_dir}/repair/additional_rules.h"
        "${lib_include_dir}/repair/geomodel_correspondence_map.h"
        "${lib_include_dir}/repair/repairer.h"
        "${lib_include_dir}/repair/topology_recovery.h"
        "${lib_include_dir}/tools/complexity_measures.h"
        "${lib_include_dir}/tools/convex_shape.h"
        "${lib_include_dir}/tools/distance.h"
        "${lib_include_dir}/tools/entity_analysis.h"
        "${lib_include_dir}/tools/geometry.h"
        "${lib_include_dir}/tools/projection.h"
        "${lib_include_dir}/tools/triangle_quality.h"
        "${lib_include_dir}/tools/utils.h"
	#TO EDIT: fill in the list of includes as follow:
	#	"${lib_include_dir}/my_file.h"

)

# List of linked targets: 
# This part is commented because the current library does not depend on any 
# other target.
#--------------------------------------------------
# The following function declares the dependency of the target (${target_name}) 
# to other target library.
#
# You must not to edit variables ${target_name} because it is automatically set 
# by a macro. 
#
# You need to complete the function adding dependent target name after the Keyword 
# PRIVATE and before the ")".
#
# The target name of a library created in a ringlab module is 
# "MODULE_NAMESPACE::library_name" (ex: RINGMesh::basic, RINGpcl::ringpcl) 
#--------------------------------------------------

set(RINGMesh_link_libraries RINGMesh::basic 
RINGMesh::geomodel_core 
RINGMesh::geomodel_builder 
RINGMesh::geomodel_tools 
RINGMesh::io 
RINGMesh::geogram_extension)
if(SCAR_WITH_APPLICATIONS)
	set(RINGMesh_link_libraries ${RINGMesh_link_libraries} RINGMesh::visualize)
endif()


target_link_libraries(${target_name} 
   PRIVATE
        ${RINGMesh_link_libraries}
		${MMG_LIBRARY}
		ringpcl
		
)


