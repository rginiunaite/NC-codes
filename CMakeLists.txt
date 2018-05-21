cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Aboria/cmake"
        ${CMAKE_MODULE_PATH})

# Boost
find_package(Boost 1.50.0 COMPONENTS python REQUIRED)
list(APPEND LIBRARIES ${Boost_LIBRARIES})
list(APPEND INCLUDES ${Boost_INCLUDE_DIRS})

# VTK
find_package(VTK REQUIRED)
if (VTK_FOUND)
    add_definitions(-DHAVE_VTK)
endif(VTK_FOUND)
list(APPEND LIBRARIES ${VTK_LIBRARIES})
list(APPEND INCLUDES ${VTK_INCLUDE_DIRS})

# Eigen
find_package(Eigen3 REQUIRED)
list(APPEND INCLUDES ${EIGEN3_INCLUDE_DIR})

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


# Aboria
set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
list(APPEND INCLUDES Aboria/src)
list(APPEND INCLUDES Aboria/third-party)

include_directories(src ${INCLUDES})



#add_executable(one_type_cells one_type_cells.cpp)
#add_executable(two_types_open_extra two_types_open_extra.cpp)
#add_executable(two_types_open_extra_main two_types_open_extra_main.cpp)
#add_executable(forced_gradient forced_gradient.cpp)
#add_executable(two_types_chains two_types_chains.cpp)
#add_executable(two_types_chains_non_uniform two_types_chains_non_uniform.cpp)
#add_executable(forced_gradient_non_uniform forced_gradient_non_uniform.cpp)
#add_executable(one_type_cells_non_uniform one_type_cells_non_uniform.cpp)
#add_executable(open_space_improved_version open_space_improved_version.cpp)
#add_executable(check2 two_types_open_extra_main.cpp)
#add_executable(fixed_gradient_3D fixed_gradient_3D.cpp)
#add_executable(touching touching.cpp)
#add_executable(combination_touching_open_space combination_touching_open_space.cpp)
#add_executable(combination_touching_open_space_main combination_touching_open_space_main.cpp)
#add_executable(touching_3angles touching_3angles.cpp)
#add_executable(persistent_movement persistent_movement.cpp)
#add_executable(persistent_movement_many_sim persistent_movement_many_sim.cpp)
#add_executable(touching_3angles_many_sim touching_3angles_many_sim.cpp)
add_executable(final final.cpp)




#target_link_libraries(one_type_cells ${LIBRARIES})
#target_link_libraries(two_types_open_extra ${LIBRARIES})
#target_link_libraries(two_types_open_extra_main ${LIBRARIES})
#target_link_libraries(forced_gradient ${LIBRARIES})
#target_link_libraries(two_types_chains ${LIBRARIES})
#target_link_libraries(two_types_chains_non_uniform ${LIBRARIES})
#target_link_libraries(forced_gradient_non_uniform ${LIBRARIES})
#target_link_libraries(one_type_cells_non_uniform ${LIBRARIES})
#target_link_libraries(open_space_improved_version ${LIBRARIES})
#target_link_libraries(check2 ${LIBRARIES})
#target_link_libraries(fixed_gradient_3D ${LIBRARIES})
#target_link_libraries(touching ${LIBRARIES})
#target_link_libraries(combination_touching_open_space ${LIBRARIES})
#target_link_libraries(combination_touching_open_space_main ${LIBRARIES})
#target_link_libraries(touching_3angles ${LIBRARIES})
#target_link_libraries(persistent_movement ${LIBRARIES})
#target_link_libraries(persistent_movement_many_sim ${LIBRARIES})
#target_link_libraries(touching_3angles_many_sim ${LIBRARIES})
target_link_libraries(final ${LIBRARIES})


