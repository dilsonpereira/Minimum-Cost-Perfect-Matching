cmake_minimum_required(VERSION 3.7)

project(mincostmatching)

option(MCPM_BUILD_TESTS "Build mincost perfect matching test executable" ON)
option(MCPM_BUILD_EXAMPLE "Build mincost perfect matching example executable" ON)
enable_language(C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
set(MCPM_INCLUDE_LIBS src/ examples/ test/ include/)

# build gtest
add_subdirectory(external)


if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror -fdiagnostics-color=always")
endif()
if ( MSVC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4 /WX")
endif()

add_library(mincostmatching)


target_include_directories(mincostmatching PRIVATE ${MCPM_INCLUDE_LIBS})

add_subdirectory(src)

if ( ${MCPM_BUILD_TESTS} )
    add_executable(mincostmatching_test)
    add_subdirectory(test)
    target_include_directories(mincostmatching_test PRIVATE src ${MCPM_INCLUDE_LIBS})
    target_link_libraries(mincostmatching_test mincostmatching gtest_main)
endif()


if ( ${MCPM_BUILD_EXAMPLE} )
    add_executable(mincostmatching_example)
    add_subdirectory(example)
    target_include_directories(mincostmatching_example PRIVATE src ${MCPM_INCLUDE_LIBS})
    target_link_libraries(mincostmatching_example mincostmatching)
endif()
