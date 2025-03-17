file(GLOB PIPELINE 
    *.cpp 
    *.hpp 
    GPU.cpp
    RW_Window.cpp
)
add_library(renderPipeline ${PIPELINE})

