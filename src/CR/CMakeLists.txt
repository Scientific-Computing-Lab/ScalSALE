include_directories($ENV{SCR_PATH}/scr/include)
link_directories($ENV{SCR_PATH}/scr/lib64)
add_library(CR SHARED cr.f90)

target_link_libraries(CR scr)
target_link_libraries(CR scrf)