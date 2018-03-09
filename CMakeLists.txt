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

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11")


# Aboria
set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
list(APPEND INCLUDES Aboria/src)
list(APPEND INCLUDES Aboria/third-party)

include_directories(src ${INCLUDES})



#add_executable(one_type_cells one_type_cells.cpp)
#add_executable(two_types_open_extra two_types_open_extra.cpp)
#add_executable(forced_gradient forced_gradient.cpp)
add_executable(two_types_chains two_types_chains.cpp)


#target_link_libraries(one_type_cells ${LIBRARIES})
#target_link_libraries(two_types_open_extra ${LIBRARIES})
#target_link_libraries(forced_gradient ${LIBRARIES})
target_link_libraries(two_types_chains ${LIBRARIES})



