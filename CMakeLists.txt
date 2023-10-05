cmake_minimum_required(VERSION 3.15)
project(minimal_plasma)

enable_language(C)

set(PLASMA_SOURCE
        control/control.c
        control/descriptor.c
        control/context.c
        control/auxiliary.c
        control/allocate.c
        control/async.c
        control/plasmaos.c
        control/plasma_threadsetting.c
        quark/quark.c
        quark/icl_hash.c
        quark/icl_list.c
        quark/quarkos.c
        core_blas/global.c
        core_blas/primes.c
        core_blas/gkkleader.c
        core_blas/descriptor.c
        core_blas-qwrapper/qwrapper_alloc.c

        control/pdtile.c
        core_blas/core_dpotrf.c
        core_blas/core_dsyrk.c
        core_blas/core_dtrsm.c
        core_blas/core_dgemm.c
        core_blas/core_dlacpy.c
        core_blas/core_dgetrip.c
        core_blas/core_dswpab.c
        core_blas/core_dshift.c
        core_blas/core_dplgsy.c
        core_blas-qwrapper/qwrapper_dlacpy.c
        core_blas-qwrapper/qwrapper_dtrsm.c
        core_blas-qwrapper/qwrapper_dsyrk.c
        core_blas-qwrapper/qwrapper_dpotrf.c
        core_blas-qwrapper/qwrapper_dgemm.c
        core_blas-qwrapper/qwrapper_dgetrip.c
        core_blas-qwrapper/qwrapper_dshift.c
        core_blas-qwrapper/qwrapper_dplgsy.c

        compute/dpotrf.c
        compute/pdpotrf.c
        compute/pdpack.c
        compute/pdgetmi2.c
        compute/pdshift.c
        compute/dgecfi.c
        compute/dgecfi2.c
        compute/dlacpy.c
        compute/dplgsy.c
        compute/pdlacpy.c
        compute/pdplgsy.c

        #  POTRI
        #  core_blas/core_dtrtri.c
        #  core_blas-qwrapper/qwrapper_dtrtri.c
        #  core_blas-qwrapper/qwrapper_dlauum.c
        #  core_blas-qwrapper/qwrapper_dtrmm.c
        #  compute/dpotri.c
        #  compute/pdtrtri.c
        #  compute/pdlauum.c
)

add_library(plasma STATIC ${PLASMA_SOURCE})
target_include_directories(plasma PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/quark>
        $<INSTALL_INTERFACE:include>
)

target_include_directories(plasma PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/control>
        $<INSTALL_INTERFACE:include>
)
target_compile_definitions(plasma PRIVATE $<$<COMPILE_LANGUAGE:C>:ADD_>)
find_package(MKL)
if (NOT ${MKL_FOUND})
    message(STATUS "MKL NOT FOUND")
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif ()
find_package(OpenMP REQUIRED)
if (${MKL_FOUND})
    target_link_libraries(plasma PRIVATE MKL::MKL)
    target_compile_definitions(plasma PUBLIC $<$<COMPILE_LANGUAGE:C>:PLASMA_WITH_MKL>)
else ()
    target_link_libraries(plasma PRIVATE BLAS::BLAS LAPACK::LAPACK lapacke)
endif ()
target_link_libraries(plasma PRIVATE OpenMP::OpenMP_C)
target_compile_options(plasma PRIVATE $<$<COMPILE_LANGUAGE:C>:-O3>)

#install(TARGETS plasma
#        EXPORT plasmaTargets
#        DESTINATION plasma)
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/control/ DESTINATION include)
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/quark/ DESTINATION include)
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/compute/ DESTINATION include)
#
#install(EXPORT plasmaTargets
#        FILE plasmaTargets.cmake
#        NAMESPACE PLASMA::
#        DESTINATION lib/cmake/PLASMA
#)

add_executable(example_dpotrf examples/example_dpotrf.c)
target_link_libraries(example_dpotrf PRIVATE plasma)
target_compile_options(example_dpotrf PRIVATE $<$<COMPILE_LANGUAGE:C>:-O3>)
