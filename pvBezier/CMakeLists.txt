set(LIBRARY_NAME "pvBezierLib")
add_library(${LIBRARY_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/src/bernstein_basis.cpp
        include/pvBezier/partial_t.hpp)

target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(${LIBRARY_NAME} PUBLIC my_eigen)