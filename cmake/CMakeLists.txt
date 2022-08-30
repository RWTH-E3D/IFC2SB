# CMake file for Project IFC2SB

cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(IFC2SB VERSION 1.0)

# FLAG for static OCC, IFCOS, QT5
OPTION(STATIC_OCC_IFCOS "Build IFC2SB and SBV with static OCC, IFCOpenShell and QT5" OFF)

# Create cache entries if absent for environment variables
MACRO(UNIFY_ENVVARS_AND_CACHE VAR)
    IF ((NOT DEFINED ${VAR}) AND ("$ENV{${VAR}}" STREQUAL ""))
        SET(${VAR} "$ENV{${VAR}}" CACHE STRING "${VAR}" FORCE)
    ENDIF ()
ENDMACRO()

# Qt5
set(CMAKE_AUTOMOC ON) # Instruct CMAKE to run MOC automatically when needed
find_package(Qt5 COMPONENTS Core Widgets Gui REQUIRED)

# messages ###########################################
message("\n### Qt5\n")
message(STATUS "Qt5Core include files found in ...")
foreach (line IN LISTS Qt5Core_INCLUDE_DIRS)
    message("\t${line}")
endforeach ()
message("\n")
message(STATUS "Qt5Widgets include files found in ...")
foreach (line IN LISTS Qt5Widgets_INCLUDE_DIRS)
    message("\t${line}")
endforeach ()
message("\n")
message(STATUS "Qt5Gui include files found in ...")
foreach (line IN LISTS Qt5Gui_INCLUDE_DIRS)
    message("\t${line}")
endforeach ()
###################################################### 	

# OpenCascade
SET(STATIC_OCC_INCLUDE_PATH "/home/fluid/Downloads/Libraries/opencascade-7.4.0_static/build/include/opencascade")
SET(STATIC_OCC_LIBRARY_PATH "/home/fluid/Downloads/Libraries/opencascade-7.4.0_static/build/lin64/gcc/lib")

message("\n### OpenCascade\n")
if (STATIC_OCC_IFCOS)
    if (NOT STATIC_OCC_INCLUDE_PATH STREQUAL "")
        #DO NOTHING
    else ()
        message(STATUS "STATIC_OCC_INCLUDE_PATH at ${STATIC_OCC_INCLUDE_PATH}")
        message(FATAL_ERROR "Undefined STATIC_OCC_INCLUDE_PATH env var.")
    endif ()
    if (NOT STATIC_OCC_LIBRARY_PATH STREQUAL "")
        link_directories(${STATIC_OCC_LIBRARY_PATH})
    else ()
        message(STATUS "STATIC_OCC_LIBRARY_PATH at ${STATIC_OCC_LIBRARY_PATH}")
        message(FATAL_ERROR "Undefined STATIC_OCC_LIBRARY_PATH env var.")
    endif ()

    SET(OPENCASCADE_LIBS "libTKTObjDRAW.a;libTKQADraw.a;libTKXDEDRAW.a;libTKDCAF.a;libTKXSDRAW.a;libTKViewerTest.a;libTKTopTest.a;libTKDraw.a;libTKRWMesh.a;libTKBinXCAF.a;libTKXmlXCAF.a;libTKVRML.a;libTKSTL.a;libTKXDESTEP.a;libTKXDEIGES.a;libTKXCAF.a;libTKIGES.a;libTKSTEP.a;libTKSTEP209.a;libTKSTEPAttr.a;libTKSTEPBase.a;libTKXSBase.a;libTKVCAF.a;libTKXmlTObj.a;libTKBinTObj.a;libTKTObj.a;libTKStd.a;libTKStdL.a;libTKXml.a;libTKBin.a;libTKXmlL.a;libTKBinL.a;libTKCAF.a;libTKLCAF.a;libTKCDF.a;libTKMeshVS.a;libTKOpenGl.a;libTKV3d.a;libTKService.a;libTKXMesh.a;libTKMesh.a;libTKFeat.a;libTKOffset.a;libTKFillet.a;libTKHLR.a;libTKBool.a;libTKBO.a;libTKShHealing.a;libTKPrim.a;libTKTopAlgo.a;libTKGeomAlgo.a;libTKBRep.a;libTKGeomBase.a;libTKG3d.a;libTKG2d.a;libTKMath.a;libTKernel.a;")
else ()
    find_package(OpenCASCADE)

    # messages ###########################################
    message("\n### OpenCASCADE\n")
    message(STATUS "OpenCASCADE include files found in ...")
    foreach (line IN LISTS OpenCASCADE_INCLUDE_DIR)
        message("\t${line}")
    endforeach ()
    # reduce linked libraries:
    SET(OPENCASCADE_LIBS "TKGeomBase;TKernel;TKG3d;TKBRep;TKTopAlgo;TKG2d;TKMath;TKPrim;TKFillet;TKBO;TKOffset;TKShHealing;TKV3d;TKSTEP;TKXSBase;TKService;TKOpenGl;TKIGES;TKVRML;TKSTL;TKSTEPBase;TKGeomAlgo;TKMesh;TKBool;")
    ######################################################
endif ()

message(STATUS "Used OpenCASCADE libraries:")
foreach (line IN LISTS OPENCASCADE_LIBS)
    message("\t${line}")
endforeach ()

# IfcOpenShell
message("### IfcOpenShell\n")
SET(STATIC_IFCOS_INCLUDE_PATH "/home/fluid/Downloads/Libraries/IfcOpenShell-0.6.0_2/src")
SET(STATIC_IFCOS_LIBRARY_PATH "/home/fluid/Downloads/Libraries/IfcOpenShell-0.6.0_2/build")
SET(SHARED_IFCOS_INCLUDE_PATH "/usr/local/include/")

if (STATIC_OCC_IFCOS)
    if (NOT STATIC_IFCOS_INCLUDE_PATH STREQUAL "")
        SET(IFCOS_INCLUDE_DIR "${STATIC_IFCOS_INCLUDE_PATH}")
    else ()
        message(STATUS "STATIC_IFCOS_INCLUDE_PATH at ${STATIC_IFCOS_INCLUDE_PATH}")
        message(FATAL_ERROR "Undefined STATIC_IFCOS_INCLUDE_PATH env var.")
    endif ()
    if (NOT STATIC_IFCOS_LIBRARY_PATH STREQUAL "")
        link_directories(${STATIC_IFCOS_LIBRARY_PATH})
    else ()
        message(STATUS "STATIC_IFCOS_LIBRARY_PATH at ${STATIC_IFCOS_LIBRARY_PATH}")
        message(FATAL_ERROR "Undefined STATIC_IFCOS_LIBRARY_PATH env var.")
    endif ()
    SET(IFCOS_LIBRARIES "libIfcGeom.a;libIfcGeom_ifc2x3.a;libIfcGeom_ifc4.a;libIfcGeom_ifc4x1.a;libIfcGeom_ifc4x2.a;libIfcGeom_ifc4x3_rc1.a;libIfcParse.a;libSerializers.a;libSerializers_ifc2x3.a;libSerializers_ifc4.a;libSerializers_ifc4x1.a;libSerializers_ifc4x2.a;libSerializers_ifc4x3_rc1.a;")
else (STATIC_OCC_IFCOS)
    if (NOT SHARED_IFCOS_INCLUDE_PATH STREQUAL "")
        SET(IFCOS_INCLUDE_DIR "${SHARED_IFCOS_INCLUDE_PATH}")
    else ()
        message(STATUS "SHARED_IFCOS_INCLUDE_PATH at ${SHARED_IFCOS_INCLUDE_PATH}")
        message(FATAL_ERROR "Undefined SHARED_IFCOS_INCLUDE_PATH env var.")
    endif ()
    SET(IFCOS_LIBRARIES "IfcParse;IfcGeom")
endif (STATIC_OCC_IFCOS)

IF ("${IFCOS_INCLUDE_DIR}" STREQUAL "")
    SET(IFCOS_INCLUDE_DIR "/usr/local/include" CACHE FILEPATH "IFC OpenShell header files")
    message(STATUS "Looking for IFC OpenShell include files in ${IFCOS_INCLUDE_DIR}")
    message(STATUS "Use IFCOS_INCLUDE_DIR to specify another directory")
ELSE ()
    SET(IFCOS_INCLUDE_DIR ${IFCOS_INCLUDE_DIR} CACHE FILEPATH "IFC OpehShell header files")
    message(STATUS "Looking for IFC OpenShell include files in: ${IFCOS_INCLUDE_DIR}")
ENDIF ()

# set ifcparse and ifcgeom include directories, libraries
SET(IFCOSP_INCLUDE_DIR "${IFCOS_INCLUDE_DIR}/ifcparse")
SET(IFCOSG_INCLUDE_DIR "${IFCOS_INCLUDE_DIR}/ifcgeom")

# check ifcparse
FIND_FILE(Ifc4_h "Ifc4.h" ${IFCOSP_INCLUDE_DIR})
IF (Ifc4_h)
    message(STATUS "IfcParse Header files found")
ELSE ()
    message(FATAL_ERROR "Unable to find header files, aborting")
ENDIF ()

# check ifcgeom
FIND_FILE(IfcGeom_h "IfcGeom.h" ${IFCOSG_INCLUDE_DIR})
IF (IfcGeom_h)
    message(STATUS "IFCGeom Header files found")
ELSE ()
    message(FATAL_ERROR "Unable to find header files, aborting")
ENDIF ()

# Intel tbb
find_package(TBB REQUIRED tbb)
message("\n### TBB\n")
message(STATUS "TBB library found (${TBB_VERSION}) in ${TBB_IMPORTED_TARGETS}")

# OpenMP
find_package(OpenMP REQUIRED)

# Boost
set(Boost_USE_STATIC_LIBS ON)  # only find static libs
set(Boost_USE_DEBUG_LIBS ON) # ignore debug libs and
set(Boost_USE_RELEASE_LIBS ON)  # only find release libs
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.71.0 COMPONENTS program_options filesystem regex REQUIRED)
message("\n### Boost\n")
message(STATUS "Boost include files found in ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries found in ${Boost_LIBRARY_DIRS}")
message("\n###\n")

# Clipper
set(CLIPPER_INCLUDE_PATH "/usr/local/include/polyclipping" CACHE STRING "Polyclipping include path")

# include all directories
if (STATIC_OCC_IFCOS)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${IFCOSG_INCLUDE_DIR} ${IFCOSP_INCLUDE_DIR} ${STATIC_OCC_INCLUDE_PATH} ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${CLIPPER_INCLUDE_PATH})
    SET(ADD_LIBS "pthread;dl;X11;xcb;Xau;Xdmcp;fontconfig;freetype;expat;png16;z;GL;xkbcommon;icui18n;icuuc;icudata;tbbmalloc;tbb;polyclipping;gmp;mpfr") # wouldn't compile without tbbmalloc;tbb stating explicite
    SET(STATIC_LIBS "${IFCOS_LIBRARIES};${OPENCASCADE_LIBS};")
    SET(TARGET_LIBS "${STATIC_LIBS};${STATIC_LIBS};${STATIC_LIBS};${Qt5Core_LIBRARIES};${Qt5Widgets_LIBRARIES};${Qt5Gui_LIBRARIES};${TBB_IMPORTED_TARGETS};${ADD_LIBS};${Boost_LIBRARIES};")
    SET(TARGET_LIBS_SBV "${STATIC_LIBS};${STATIC_LIBS};${Qt5Core_LIBRARIES};${Qt5Widgets_LIBRARIES};${Qt5Gui_LIBRARIES};${ADD_LIBS};")
else (STATIC_OCC_IFCOS)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${IFCOSG_INCLUDE_DIR} ${IFCOSP_INCLUDE_DIR} ${OpenCASCADE_INCLUDE_DIR} ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${CLIPPER_INCLUDE_PATH})
    SET(TARGET_LIBS "${Qt5Core_LIBRARIES};${Qt5Widgets_LIBRARIES};${Qt5Gui_LIBRARIES};${OpenCASCADE_LIBRARIES};${TBB_IMPORTED_TARGETS};${IFCOS_LIBRARIES};${Boost_LIBRARIES};polyclipping;gmp;mpfr")
    SET(TARGET_LIBS_SBV "${Qt5Core_LIBRARIES};${Qt5Widgets_LIBRARIES};${Qt5Gui_LIBRARIES};${OpenCASCADE_LIBRARIES};${IFCOS_LIBRARIES}")
endif (STATIC_OCC_IFCOS)

# source files name extraction
file(GLOB_RECURSE IFC2SB_SOURCE_CPP ../src/*.cpp)
file(GLOB_RECURSE IFC2SB_SOURCE_CXX ../src/*.cxx)

# add executable IFC2SB
add_executable(
        IFC2SB
        ${MOC_HEADERS}
        ${IFC2SB_SOURCE_CPP}
        ${IFC2SB_SOURCE_CXX}
)

target_link_libraries(
        IFC2SB
        ${TARGET_LIBS}
)

# SBV
file(GLOB_RECURSE SBV_SOURCE_CPP ../sb_viewer/*.cpp)
file(GLOB_RECURSE SBV_VIEWER_CXX ../src/viewer/*.cxx)
file(GLOB_RECURSE SBV_VIEWER_CPP ../src/viewer/*.cpp)

# add executable SBV
add_executable(
        SBV
        ${MOC_HEADERS}
        ${SBV_SOURCE_CPP}
        ${SBV_VIEWER_CXX}
        ${SBV_VIEWER_CPP}
)

target_link_libraries(
        SBV
        ${TARGET_LIBS_SBV}
)

# By default build type is set to Release
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -fopenmp -O3")

message(STATUS "Static build: ${STATIC_OCC_IFCOS}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message("\n###\n")
