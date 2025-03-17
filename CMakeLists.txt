cmake_minimum_required(VERSION 3.12)
project(SoftRaster)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
if(BUILD_TYPE_UPPER STREQUAL "RELEASE")
    add_compile_options(
        $<$<CXX_COMPILER_ID:MSVC>:/O2>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3>
        -DNDEBUG
    )
endif()

add_definitions(-DUNICODE -D_UNICODE)

add_subdirectory(funcs)
add_subdirectory(math)
add_subdirectory(renderPipeline)

set(MAIN_SRCS
    main.cpp
    renderPipeline/RS_Window.cpp
    renderPipeline/GPU.cpp
    math/VecMath.cpp
    funcs/frameBufferObject.cpp
)

add_executable(softRaster ${MAIN_SRCS})

target_link_libraries(softRaster
    PRIVATE
        funcs
        math
        renderPipeline
)