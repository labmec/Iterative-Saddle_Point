set(MESHESDIR ${CMAKE_SOURCE_DIR}/targets/meshes)
set(CONFIG_H "${CMAKE_SOURCE_DIR}/targets/meshpath_config.h")
configure_file("${CMAKE_SOURCE_DIR}/targets/meshpath_config.h.in" ${CONFIG_H})

add_executable(iterative-no-condense-darcy iterative-no-condense-darcy.cpp)
target_link_libraries(iterative-no-condense-darcy PRIVATE libproject)

add_executable(iterative-condensed-darcy iterative-condensed-darcy.cpp)
target_link_libraries(iterative-condensed-darcy PRIVATE libproject)

add_executable(iterative-condensed-stokes iterative-condensed-stokes.cpp)
target_link_libraries(iterative-condensed-stokes PRIVATE libproject)
