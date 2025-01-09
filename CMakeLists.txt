cmake_minimum_required (VERSION 3.12)

project(Gemini++
   VERSION 2014.10.2
   DESCRIPTION "Gemini++ code to model the statistical decay of compound nuclei"
   LANGUAGES CXX
   )

list(APPEND CMAKE_MODULE_PATH
      	       ${PROJECT_SOURCE_DIR}/cmake/CMakeListsUtilities
	       )

include(GNUPackageInstallDirs)

include(StandardRpathHandling)

set(SOURCES
   Angle.cpp      Fus.cpp           LightP.cpp   Random.cpp      TlArray.cpp
   AngleDist.cpp  Gdr.cpp           Mass.cpp     Scission.cpp    TlBarDist.cpp
   Chart.cpp      History.cpp       Nucleus.cpp  SigBarDist.cpp  Weight.cpp
   Evap.cpp       LevelDensity.cpp  Nuclide.cpp  SigCharged.cpp  Yrast.cpp
   )

set(HEADERS
   CAngleDist.h  CFus.h           CLightP.h   CRandom.h    CSigBarDist.h  CWeight.h
   CAngle.h      CGdr.h           CMass.h     CRun.h       CSigCharged.h  CYrast.h
   CChart.h      CHistory.h       CNucleus.h  CRunThick.h  CTlArray.h     SStoreEvap.h
   CEvap.h       CLevelDensity.h  CNuclide.h  CScission.h  CTlBarDist.h
   )

add_library(Gemini++ SHARED ${SOURCES})

target_compile_definitions(Gemini++
      	             	      PUBLIC GINPUT="${CMAKE_INSTALL_FULL_PKGDATADIR}/")

install(TARGETS Gemini++
 EXPORT ${PROJECT_NAME}Exports
 LIBRARY DESTINATION ${CMAKE_INSTALL_PKGLIBDIR}
 INCLUDES DESTINATION ${CMAKE_INSTALL_PKGINCDIR}
)
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PKGINCDIR})
install(DIRECTORY tbl tl DESTINATION ${CMAKE_INSTALL_PKGDATADIR})

set(PROJECT_CONFIG_TEMPLATE_FILE
   ${PROJECT_SOURCE_DIR}/cmake/CMakeListsUtilities/projectConfig.cmake.in
   )
set(NO_PROJECT_TARGETS_NAMESPACE ON)
include(InstallConfigPackageVersion)
