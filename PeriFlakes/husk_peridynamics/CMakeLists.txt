
cmake_minimum_required(VERSION 2.8.9)

set(huskname peridynamics)

set(KERNEL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
set(KERNEL_FILES
   ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_quarticA.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_quadr.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_quarticA.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_rational_2_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_quarticA.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_rational_2_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_rational_0_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_rational_1_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_quadr.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_quadr.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_rational_2_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_rational_2_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_rational_0_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_rational_0_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_inv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_rational_0_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_rational_0_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_rational_1_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_quadr.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_rational_0_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_rational_2_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_inv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_quarticA.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_cubic.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_linear.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_cubic.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_quadr.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_const.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_rational_0_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_rational_2_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_const.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_rational_0_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_rational_0_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_rational_0_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_cubic.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_cubic.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_rational_0_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_rational_1_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_quarticA.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_inv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_rational_1_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_rational_1_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_linear.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_rational_2_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_rational_1_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_const.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_rational_0_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_rational_2_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_rational_1_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_rational_2_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_rational_0_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_linear.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_rational_0_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_linear.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_inv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_rational_0_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_rational_2_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_const.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_rational_0_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_rational_2_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_rational_2_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_rational_2_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_inv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_const.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_rational_0_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_const.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_rational_2_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_quarticA.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_linear.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_quarticA.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_cubic.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_rational_1_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_quarticA.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_quadr.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_cubic.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_rational_1_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_inv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_rational_1_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_rational_0_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_rational_1_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_linear.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_inv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_rational_0_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_rational_0_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_rational_1_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_rational_2_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_rational_2_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_quadr.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_inv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_rational_1_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_linear.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_rational_1_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_rational_1_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_linear.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_rational_0_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_rational_0_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_rational_0_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_cubic.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_rational_2_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_rational_2_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_rational_2_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_rational_1_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_rational_1_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_rational_2_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_rational_1_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_quadr.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/smooth_cubic.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_rational_1_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_rational_1_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_rational_1_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_cubic.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_rational_0_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_rational_0_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_rational_1_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_rational_1_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_rational_2_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_const.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_linear.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_quadr.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_rational_2_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_rational_2_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_rational_2_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_rational_1_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_rational_0_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Silling_rational_1_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Oterkus2_rational_1_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_rational_2_1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood_rational_2_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011minus_const.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Silling_quarticA.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_rational_0_2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_const.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fstab_Littlewood2011_rational_2_0.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/BondBased_inv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/Fbased_rational_1_0.c 
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
