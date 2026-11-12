# The interactive explainer (see README.md in this directory).
#
# With Emscripten this builds the WebAssembly module which the page loads.
# With any other compiler it builds explainer_check, which runs the same
# interface natively so that it can be compared with euroc_benchmark.

set(explainer_include_dirs
    ${CMAKE_CURRENT_SOURCE_DIR}/../examples/euroc_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/../examples/attitude_benchmark)

if(EMSCRIPTEN)
    add_executable(ukf_explainer wasm/explainer_api.cpp)
    target_link_libraries(ukf_explainer PRIVATE ukf::double)
    target_include_directories(ukf_explainer PRIVATE ${explainer_include_dirs})
    target_compile_options(ukf_explainer PRIVATE -fno-rtti)
    # One self-contained .wasm file with no JavaScript glue: the page loads
    # it directly. Memory can grow, because recordings differ in length.
    target_link_options(ukf_explainer PRIVATE
        -sSTANDALONE_WASM=1
        --no-entry
        -sALLOW_MEMORY_GROWTH=1
        -sINITIAL_MEMORY=33554432
        -sSTACK_SIZE=1048576
        -sMALLOC=emmalloc
        -fno-rtti)
    set_target_properties(ukf_explainer PROPERTIES SUFFIX ".wasm")
    ukf_project_warnings(ukf_explainer)
else()
    add_executable(explainer_check wasm/explainer_check.cpp wasm/explainer_api.cpp)
    target_link_libraries(explainer_check PRIVATE ukf::double)
    target_include_directories(explainer_check PRIVATE ${explainer_include_dirs})
    ukf_project_warnings(explainer_check)
endif()
