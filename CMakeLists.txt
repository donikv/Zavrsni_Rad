cmake_minimum_required (VERSION 3.2)
project (landau_wishkin)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -march=native")

# Define output directories for created binaries and libraries.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file(GLOB SOURCES "landau_wishkin/src/*.cpp")

add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

# Create libraries.
add_library(landau_wishkin SHARED ${SOURCES})

target_include_directories(landau_wishkin PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/landau_wishkin/include>
    $<INSTALL_INTERFACE:include>)

add_library(landau_wishkin_static STATIC ${SOURCES})

target_include_directories(landau_wishkin_static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/landau_wishkin/include>
    $<INSTALL_INTERFACE:include>)

# Build binaries.
add_executable(main main.cpp)
target_link_libraries(main landau_wishkin_static pthread)

# Create target 'install' for installing libraries.
install(TARGETS landau_wishkin landau_wishkin_static DESTINATION lib)
install(FILES landau_wishkin/include/Algorithms.hpp DESTINATION include)