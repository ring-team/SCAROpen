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


#----------------
# list all binary source files and dependencies
# use the macro: ringlab_add_binary(path_to/file.cpp target_deps)

# ex: ringlab_add_binary(${binary_source_dir}/file_name.cpp ModuleNamespace1::library_name1 ModuleNamespace2::library_name2 ...
# ${binary_source_dir} is set by cmake... just use it.
#----------------

set(link_libraries_exe 
    RINGMesh::basic 
    RINGMesh::geomodel_core 
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
	scar)

if(SCAR_WITH_MMG)
	set(link_libraries_exe ${link_libraries_exe} ${MMG_LIBRARY})
endif()
	
ringlab_add_binary(${binary_source_dir}/scar_2D_model_flip_triangles.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	
ringlab_add_binary(${binary_source_dir}/scar_2D_substract.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	
ringlab_add_binary(${binary_source_dir}/scar_add_entity_boundaries.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	
ringlab_add_binary(${binary_source_dir}/scar_add_info.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	
ringlab_add_binary(${binary_source_dir}/scar_clean.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	
ringlab_add_binary(${binary_source_dir}/scar_complexity_measures.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	
ringlab_add_binary(${binary_source_dir}/scar_convert.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	
ringlab_add_binary(${binary_source_dir}/scar_dilate.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	
ringlab_add_binary(${binary_source_dir}/scar_repair.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	
ringlab_add_binary(${binary_source_dir}/scar_triangle_quality.cpp 
    RINGMesh::basic
    RINGMesh::geomodel_core
    RINGMesh::geomodel_builder 
    RINGMesh::geomodel_tools 
    RINGMesh::io 
    RINGMesh::geogram_extension
    scar)
	