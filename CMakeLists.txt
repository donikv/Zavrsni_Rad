cmake_minimum_required (VERSION 3.2)
project (landau_vishkin)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -march=native")

# Define output directories for created binaries and libraries.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file(GLOB SOURCES "landau_vishkin/src/*.cpp")

add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

# Create libraries.
add_library(landau_vishkin SHARED ${SOURCES})

target_include_directories(landau_vishkin PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/landau_vishkin/include>
    $<INSTALL_INTERFACE:include>)

add_library(landau_vishkin_static STATIC ${SOURCES})

target_include_directories(landau_vishkin_static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/landau_vishkin/include>
    $<INSTALL_INTERFACE:include>)

# Build binaries.
add_executable(main main.cpp)
target_link_libraries(main landau_vishkin_static pthread)

# Create target 'install' for installing libraries.
install(TARGETS landau_vishkin landau_vishkin_static DESTINATION lib)
install(FILES landau_vishkin/include/Algorithms.hpp DESTINATION include)