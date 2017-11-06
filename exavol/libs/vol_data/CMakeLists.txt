cmake_minimum_required(VERSION 2.6)
project(vol_data)

# The CGV framework is needed
find_package(cgv COMPONENTS  utils type reflect data media)

# Source files
file(GLOB_RECURSE SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.cxx")
# Header files
file(GLOB_RECURSE HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.h")

include_directories("${cgv_INCLUDE_DIRS}")

cgv_add_library(vol_data
	SOURCES ${SOURCES} 
	PUBLIC_HEADERS ${HEADERS}
	SHARED_DEFINITIONS VOL_DATA_EXPORTS)

	
target_link_libraries(vol_data ${cgv_LIBRARIES})

cgv_write_find_file(vol_data)

