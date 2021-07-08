project(TROEFIT Fortran)
cmake_minimum_required(VERSION 2.8)

add_executable(troefit.x ${PROJECT_SOURCE_DIR}/src/troefit.f)

<<<<<<< HEAD
install(TARGETS troefit.x)
=======
install(TARGETS troefit.x DESTINATION ${PROJECT_SOURCE_DIR}/../../bin)
>>>>>>> external build things
