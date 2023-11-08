cmake_minimum_required(VERSION 3.5)

if(NOT CMAKE_VERSION VERSION_LESS 3.14)
  cmake_policy(SET CMP0083 NEW)
endif()

set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")


set(PACKAGE_VERSION 1.0)
project(orlqp VERSION ${PACKAGE_VERSION})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(osqp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(CLN REQUIRED)
find_package(GiNaC REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include)

list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

add_library(orlqp SHARED)
target_sources(orlqp PRIVATE
    src/osqp.cpp
    src/mpc.cpp
    src/qp_array.cpp
    src/symbolic.cpp
    )

target_include_directories(orlqp PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${GINAC_INCLUDE_DIRS}
    ${CLN_INCLUDE_DIRS}
)

target_link_libraries(orlqp PUBLIC 
    osqp::osqp
    Eigen3::Eigen
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES})

set_target_properties(orlqp PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;/usr/local/lib"
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    orlqpConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(TARGETS orlqp EXPORT orlqpTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(EXPORT orlqpTargets
    FILE orlqpTargets.cmake
    NAMESPACE orlqp::
    DESTINATION lib/cmake/orlqp
)

configure_file(cmake/orlqpConfig.cmake.in orlqpConfig.cmake @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/orlqpConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/orlqpConfigVersion.cmake"
        DESTINATION lib/cmake/orlqp
)

install(DIRECTORY include/ DESTINATION include)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    enable_testing()

    # Add the executable
    add_executable(simple_mpc_test tests/simple_mpc.cpp)
    target_include_directories(simple_mpc_test PUBLIC include)
    target_link_libraries(simple_mpc_test PUBLIC stdc++ m orlqp osqp::osqp)

    # Add the executable
    add_executable(simple_mpc_array_test tests/simple_mpc_array.cpp)
    target_include_directories(simple_mpc_array_test PUBLIC include)
    target_link_libraries(simple_mpc_array_test PUBLIC stdc++ m orlqp osqp::osqp)

    # Add the executable
    add_executable(tianze_lip tests/tianze_lip.cpp)
    target_include_directories(tianze_lip PUBLIC include)
    target_link_libraries(tianze_lip PUBLIC stdc++ m orlqp osqp::osqp ${GINAC_LIBRARIES} ${CLN_LIBRARIES})
endif()