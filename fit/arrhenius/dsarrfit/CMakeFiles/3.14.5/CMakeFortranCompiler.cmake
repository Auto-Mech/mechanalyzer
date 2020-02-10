set(CMAKE_Fortran_COMPILER "/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "7.3.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")


set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/usr/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/usr/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/lib/gcc/x86_64-conda_cos6-linux-gnu/7.3.0/finclude;/blues/gpfs/home/software/spack-0.10.1/opt/spack/linux-centos7-x86_64/intel-17.0.4/intel-mkl-2017.3.196-v7uuj6zmthzln35n2hb7i5u5ybncv5ev/compilers_and_libraries_2017.4.196/linux/mkl/include;/blues/gpfs/home/software/spack-0.10.1/opt/spack/linux-centos7-x86_64/gcc-4.8.5/intel-17.0.4-74uvhjiulyqgvsmywifbbuo46v5n42xc/tbb/include;/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/lib/gcc/x86_64-conda_cos6-linux-gnu/7.3.0/include;/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/lib/gcc/x86_64-conda_cos6-linux-gnu/7.3.0/include-fixed;/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/x86_64-conda_cos6-linux-gnu/include;/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/x86_64-conda_cos6-linux-gnu/sysroot/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gomp;gcc_s;gcc;quadmath;m;gcc_s;gcc;pthread;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/lib/gcc/x86_64-conda_cos6-linux-gnu/7.3.0;/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/lib/gcc;/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/x86_64-conda_cos6-linux-gnu/lib;/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/x86_64-conda_cos6-linux-gnu/sysroot/lib;/blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
