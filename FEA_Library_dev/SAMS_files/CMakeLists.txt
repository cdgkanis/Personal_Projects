cmake_minimum_required(VERSION 3.29)
project(SAMS)

set(CMAKE_CXX_STANDARD 20)

add_executable(SAMS
        main.cpp
        Physics.cpp
        Commands.cpp
        Import_Handler.cpp
        Import_Handler.h
)

target_include_directories(SAMS PRIVATE
        ${PROJECT_SOURCE_DIR}
)

add_custom_command(TARGET SAMS POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_SOURCE_DIR}/Materials.json
        $<TARGET_FILE_DIR:SAMS>/Materials.json
)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0 -g1)
endif()

set_source_files_properties(Commands.cpp PROPERTIES COMPILE_OPTIONS "-g0")

add_custom_command(TARGET SAMS POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_SOURCE_DIR}/plot_transient.py
        $<TARGET_FILE_DIR:SAMS>/plot_transient.py
)