cmake_minimum_required(VERSION 3.15)

# Create main project.
project(awex VERSION 0.1.2 LANGUAGES Fortran C)

# Specify position-independent code since we will build a shared library.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_library(cawex
   src/airsea_fluxes.F90
   src/airsea_variables.F90
   src/humidity.F90
   src/longwave_radiation.F90
   src/fairall.F90
   src/t10m.F90
   src/kara.F90
   src/kondo.F90
   src/solar_zenith_angle.F90
   src/albedo_water.F90
   src/shortwave_radiation.F90
)
target_compile_definitions(cawex PRIVATE GUSTINESS)

if (MSVC)
  target_compile_options(cawex PRIVATE /libs:static)
endif()

# Find target Python interpreter.
find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_custom_command(OUTPUT _awex.c
  COMMENT "Creating ${CMAKE_CURRENT_BINARY_DIR}/_awex.c from ${CMAKE_CURRENT_SOURCE_DIR}/_awex.pyx"
  COMMAND Python3::Interpreter -m cython "${CMAKE_CURRENT_SOURCE_DIR}/_awex.pyx" -o "${CMAKE_CURRENT_BINARY_DIR}/_awex.c"
  DEPENDS _awex.pyx
  VERBATIM
)

Python3_add_library(_awex MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/_awex.c")
target_link_libraries(_awex PRIVATE cawex Python3::NumPy)

install(TARGETS _awex DESTINATION .)
