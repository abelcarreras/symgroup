# CMake project file for symgroup

cmake_minimum_required (VERSION 2.6)
project (symgroup)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
            "Choose the type of build, options are: None Debug Release."
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix "." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
        "Choose the installation directory; by default it installs in the WFNSYM directory."
        FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -ff2c -ffixed-line-length-0")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -fbacktrace -ff2c -ffixed-line-length-0")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    # ifort
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g -traceback -debug extended")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
    # g77
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# MKL library
if(WITH_MKL)
    message(STATUS "Use MKL library")
    link_directories("/opt/intel/mkl/lib/")
    include_directories("/opt/intel/compilers_and_libraries/mac/include/intel64/")
    #  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lmkl_def -lmkl_core -lmkl_lapack -lmkl_intel_thread  -lmkl_sequential -lmkl_intel_lp64 -lguide  -lpthread")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lmkl_core -lmkl_intel_thread -lmkl_sequential -lmkl_intel_lp64 -lpthread")

else(WITH_MKL)
    link_directories("/usr/lib")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack -lblas")
endif(WITH_MKL)


# build executables
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/EXE)
set( LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/LIB)
set (EXECUTABLES "symgroup")

add_executable ("symgroup"
        "src/symgroup.F"
        "src/radius.F"
        "src/connectivity.F"
        "src/jacobi.F"
        "src/linear_algebra.F"
        "src/mass.F"
        "src/measure.F"
        "src/operations.F"
        "src/read_input.F")

ADD_LIBRARY ("symgrouplib" "STATIC"
        "src/symgrouplib.F"
        "src/radius.F"
        "src/connectivity.F"
        "src/jacobi.F"
        "src/linear_algebra.F"
        "src/mass.F"
        "src/measure.F"
        "src/operations.F")

# install executables and scripts
install (TARGETS ${EXECUTABLES}
        RUNTIME DESTINATION "bin")
