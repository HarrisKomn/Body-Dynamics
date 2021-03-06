###############################################################################
cmake_minimum_required (VERSION 3.0)
project (VVR_OGL_LABORATORY)

###############################################################################

find_package(OpenGL REQUIRED)

# for rdm (emacs)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# c++11, -g option is used to export debug symbols for gdb
list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -g ${CMAKE_CXX_FLAGS}")
# for grouping into folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory ! (and give it a clever name, like bin_Visual2012_64bits/)" )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
    message( "Your Source Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
    message( "Your Build Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()

# Compile external dependencies 
add_subdirectory (external)

# On Visual 2005 and above, this module can set the debug working directory
cmake_policy(SET CMP0026 OLD)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/rpavlik-cmake-modules-fe2273")
include(CreateLaunchers)
include(MSVCMultipleProcessCompile) # /MP

###############################################################################

include_directories(
    external/glfw-3.1.2/include/GLFW/
    external/glm-0.9.7.1/
    external/glew-1.13.0/include/
    external/Simple-OpenGL-Image-Library/include/
    external/tinyxml2/include/
    external/tinyobjloader/include/
    .
)

set(ALL_LIBS
    ${OPENGL_LIBRARY}
    glfw
    GLEW_1130
    SOIL
    TINYXML2
)

add_definitions(
    -DTW_STATIC
    -DTW_NO_LIB_PRAGMA
    -DTW_NO_DIRECT3D
    -DGLEW_STATIC
    -D_CRT_SECURE_NO_WARNINGS
)

###############################################################################
# lab07
add_executable(lab07 
    lab07/lab.cpp
    lab07/RigidBody.cpp
    lab07/RigidBody.h
    lab07/Cube.cpp
    lab07/Cube.h
    lab07/Sphere.cpp
    lab07/Sphere.h
    lab07/Box.cpp
    lab07/Box.h
    lab07/MassSpringDamper.cpp
    lab07/MassSpringDamper.h
    lab07/Collision.cpp
    lab07/Collision.h
         
    common/util.cpp
    common/util.h
    common/shader.cpp
    common/shader.h
    common/camera.cpp
    common/camera.h
    common/ModelLoader.cpp
    common/ModelLoader.h
    common/texture.cpp
    common/texture.h
    
    lab07/StandardShading.fragmentshader
    lab07/StandardShading.vertexshader
)
target_link_libraries(lab07
    ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(lab07 
    PROPERTIES 
    XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lab07/"
    PROJECT_LABEL "Lab 07 - Physics"
    FOLDER "Exercise"
)
create_target_launcher(lab07 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lab07/")
create_default_target_launcher(lab07 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lab07/") 

###############################################################################

SOURCE_GROUP(common REGULAR_EXPRESSION ".*/common/.*" )
SOURCE_GROUP(shaders REGULAR_EXPRESSION ".*/.*shader$" )
SOURCE_GROUP(obj REGULAR_EXPRESSION ".*/.*obj$" )

###############################################################################
# copy 
if (NOT ${CMAKE_GENERATOR} MATCHES "Xcode" )

add_custom_command(
   TARGET lab07 POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/lab07${CMAKE_EXECUTABLE_SUFFIX}" "${CMAKE_CURRENT_SOURCE_DIR}/lab07/"
)

elseif (${CMAKE_GENERATOR} MATCHES "Xcode" )

endif (NOT ${CMAKE_GENERATOR} MATCHES "Xcode" )

