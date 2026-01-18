cmake_minimum_required(VERSION 3.16)
project(PILOTS VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PILOTS_ENABLE_OPENMP "Enable OpenMP parallelism" ON)

add_library(pilots_core
    src/LammpsDumpReader.cpp
  src/Runner.cpp
)

target_include_directories(pilots_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Version macro for results.json
target_compile_definitions(pilots_core PUBLIC PILOTS_VERSION_STR="${PROJECT_VERSION}")

# OpenMP (Level-2 optimization)
if(PILOTS_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(pilots_core PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(pilots_core PUBLIC PILOTS_HAS_OPENMP=1)
    else()
        message(STATUS "OpenMP not found; building without OpenMP")
        target_compile_definitions(pilots_core PUBLIC PILOTS_HAS_OPENMP=0)
    endif()
else()
    target_compile_definitions(pilots_core PUBLIC PILOTS_HAS_OPENMP=0)
endif()

# One-translation-unit-per-measure workflow:
# add a new measure by dropping a single .cpp file into src/measures/
# which both implements and registers the measure factory.
file(GLOB PILOTS_MEASURE_SOURCES
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/measures/*.cpp"
)

add_executable(pilots
  src/main.cpp
  ${PILOTS_MEASURE_SOURCES}
)

target_link_libraries(pilots PRIVATE pilots_core)

target_compile_definitions(pilots PRIVATE PILOTS_VERSION_STR="${PROJECT_VERSION}")

enable_testing()
add_test(
  NAME msd_config_smoke
  COMMAND pilots --config ${CMAKE_CURRENT_SOURCE_DIR}/tests/configs/msd.ini
)
