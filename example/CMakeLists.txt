
file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_sources(mincostmatching_example PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Example.cpp
)
