project(TROEFIT Fortran)
cmake_minimum_required(VERSION 2.8)

add_executable(troefit.x ${PROJECT_SOURCE_DIR}/troefit.f)

install(TARGETS troefit.x DESTINATION bin)
