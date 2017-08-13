project(cpd-examples CXX)
cmake_minimum_required(VERSION 2.8.12)

find_package(Cpd REQUIRED)
find_package(PCL REQUIRED)

# dpkg -l *libqt5*
find_package(Qt5X11Extras QUIET)
if (Qt5X11Extras_FOUND)
  include_directories(${Qt5X11Extras_INCLUDE_DIRS})
  add_definitions(${Qt5X11Extras_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5X11Extras_EXECUTABLE_COMPILE_FLAGS}")
  list(APPEND QT_LIBRARIES ${Qt5X11Extras_LIBRARIES})
endif()

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# add_executable(voxel_grid voxel_grid.cpp)
# target_link_libraries(voxel_grid ${PCL_LIBRARIES})

add_executable(cpd-rigid rigid.cpp)
target_link_libraries(cpd-rigid PRIVATE Cpd::Library-C++ ${PCL_LIBRARIES})

# add_executable(cpd-version version.cpp)
# target_link_libraries(cpd-version PRIVATE Cpd::Library-C++)

# add_executable(cpd-transform transform.cpp)
# target_link_libraries(cpd-transform PRIVATE Cpd::Library-C++)
