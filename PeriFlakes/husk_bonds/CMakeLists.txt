
cmake_minimum_required(VERSION 2.8.9)

set(huskname bonds)

set(KERNEL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
set(KERNEL_FILES
   ${CMAKE_CURRENT_SOURCE_DIR}/line_intersection.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/bond_pressure.c 
)

# Required: Cornflakes
list(APPEND CMAKE_MODULE_PATH "$ENV{CORNFLAKES_DIR}/cmake")
include(cornflakes)
include_directories(${CORNFLAKES_INCLUDES})
# Required: GSL
find_package(GSL REQUIRED)
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS "")
find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${KERNEL_INCLUDES})

swig_add_module(${huskname}_lib python ${huskname}_swig.i ${KERNEL_FILES})
set_property(SOURCE ${KERNEL_FILES} APPEND_STRING PROPERTY COMPILE_FLAGS " -fPIC")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(link_to_cornflakes_for_swig lib_cornflakes)
else()
  set(link_to_cornflakes_for_swig swig_cornflakes)
endif()

swig_link_libraries(${huskname}_lib m ${GSL_LIBRARIES} ${PYTHON_LIBRARIES}  ${link_to_cornflakes_for_swig} ${CORNFLAKES_LIBRARIES})
