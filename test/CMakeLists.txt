
add_executable(spline_test spline_test.cpp)
target_link_libraries(spline_test  liegroup)

add_executable(main main.cpp)
target_link_libraries(main  liegroup)

add_executable(expm expm.cpp)
target_link_libraries(expm  liegroup)

add_executable(moments moments.cpp)
target_link_libraries(moments  liegroup)

add_executable(se3_gen se3_gen.cpp)
target_link_libraries(se3_gen  liegroup)

add_executable(se3_moments se3_moments.cpp)
target_link_libraries(se3_moments  liegroup)

add_executable(extrinsics extrinsics.cpp)
target_link_libraries(extrinsics  liegroup)

