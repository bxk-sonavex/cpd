project(cpd-examples CXX)
cmake_minimum_required(VERSION 2.8.12)

find_package(Cpd REQUIRED)
find_package(PCL REQUIRED COMPONENTS comman io filters search visualization)

# dpkg -l *libqt5*
find_package(Qt5X11Extras QUIET)
if (Qt5X11Extras_FOUND)
  include_directories(${Qt5X11Extras_INCLUDE_DIRS})
  add_definitions(${Qt5X11Extras_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5X11Extras_EXECUTABLE_COMPILE_FLAGS}")
  list(APPEND QT_LIBRARIES ${Qt5X11Extras_LIBRARIES})
endif()

# check for OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
else(OPENMP_FOUND)
  message (STATUS "OpenMP not found")
endif()

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(cpd-rigid rigid.cpp)
target_link_libraries(cpd-rigid PRIVATE Cpd::Library-C++ ${PCL_LIBRARIES} pcl_io_ply)