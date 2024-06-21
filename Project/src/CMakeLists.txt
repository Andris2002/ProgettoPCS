list(APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/ImportExport.cpp")
list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/ImportExport.hpp")

list(APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/AlgFractures.cpp")
list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/AlgFractures.hpp")

list(APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.cpp")
list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.hpp")

list(APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/DFNstructs.cpp")
list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/DFNstructs.hpp")

list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)
