set(MESHESDIR ${CMAKE_SOURCE_DIR}/targets/meshes)
set(CONFIG_H "${CMAKE_SOURCE_DIR}/targets/meshpath_config.h")
configure_file("${CMAKE_SOURCE_DIR}/targets/meshpath_config.h.in" ${CONFIG_H})

add_executable(iterative-no-condense-darcy iterative-no-condense-darcy.cpp)
target_link_libraries(iterative-no-condense-darcy PRIVATE libproject)
