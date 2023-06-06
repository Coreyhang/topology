cmake_minimum_required(VERSION 3.5)

project (topology)

find_package ( OpenMP )
if ( OpenMP_FOUND )
    message("dkashdasojdnsalkdjn")
	list ( APPEND PRJ_COMPILE_DEF ENABLE_OPENMP )
	list ( APPEND PRJ_LIBRARIES ${OpenMP_CXX_LIBRARIES} )
	list ( APPEND PRJ_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS} )
else()
    message("sjasss")
endif ()

set(CMAKE_CXX_FLAGS "-fopenmp ${CMAKE_CXX_FLAGS} -g -Wall")

TARGET_LINK_LIBRARIES(${OpenMP_CXX_LIBRARIES})

set(SOURCES
    src/topology.c
    src/utils.c
    src/topo_gen.c
    src/cord.c
    src/floyd_warshall.c
)

add_executable(topology ${SOURCES})

target_include_directories(topology
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)
