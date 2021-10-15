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
		"${lib_source_dir}/correlation_map.cpp"
		

    PRIVATE # Could be PUBLIC from CMake 3.3
        "${lib_include_dir}/utils.h"
        "${lib_include_dir}/correlation_map.h"
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

target_link_libraries(${target_name} 
   PRIVATE
        RINGMesh::basic 

)


