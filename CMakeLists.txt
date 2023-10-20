cmake_minimum_required(VERSION 3.5)

set(PACKAGE_VERSION 1.0)
project(orlqp VERSION ${PACKAGE_VERSION})

find_package(osqp REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
include_directories(include)

add_library(orlqp SHARED)
target_sources(orlqp PRIVATE
    src/osqp.cpp
    src/mpc.cpp
    src/qp_array.cpp
    )
target_include_directories(orlqp PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)


# Generate the necessary configuration files for find_package
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    orlqpConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Define the locations that will be used by the <project-name>Config.cmake file
# This file helps users find the include directory after installation.
install(TARGETS orlqp EXPORT orlqpTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

# Export the targets to a script
install(EXPORT orlqpTargets
    FILE orlqpTargets.cmake
    NAMESPACE orlqp::
    DESTINATION lib/cmake/orlqp
)

# Create a Config file
configure_file(cmake/orlqpConfig.cmake.in orlqpConfig.cmake @ONLY)

# Install the configuration files
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/orlqpConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/orlqpConfigVersion.cmake"
        DESTINATION lib/cmake/orlqp
)

# Install public headers
install(DIRECTORY include/ DESTINATION include)

# Enable testing only if this project is being built directly
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    enable_testing()

    # Add the executable
    add_executable(simple_mpc_test tests/simple_mpc.cpp)
    target_include_directories(simple_mpc_test PUBLIC include)
    target_link_libraries(simple_mpc_test PUBLIC stdc++ stdc++fs m orlqp osqp::osqp)

    # Add the executable
    add_executable(simple_mpc_array_test tests/simple_mpc_array.cpp)
    target_include_directories(simple_mpc_array_test PUBLIC include)
    target_link_libraries(simple_mpc_array_test PUBLIC stdc++ stdc++fs m orlqp osqp::osqp)
endif()