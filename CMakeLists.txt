cmake_minimum_required (VERSION 3.12)

project(gemini++)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/CMakeListsUtilities)

include(GNUPackageInstallDirs)

include(StandardRpathHandling)

set(SOURCES Angle.cpp      Fus.cpp           LightP.cpp   Random.cpp      TlArray.cpp
AngleDist.cpp  Gdr.cpp           Mass.cpp     Scission.cpp    TlBarDist.cpp
Chart.cpp      History.cpp       Nucleus.cpp  SigBarDist.cpp  Weight.cpp
Evap.cpp       LevelDensity.cpp  Nuclide.cpp  SigCharged.cpp  Yrast.cpp)

set(HEADERS CAngleDist.h  CFus.h           CLightP.h   CRandom.h    CSigBarDist.h  CWeight.h
CAngle.h      CGdr.h           CMass.h     CRun.h       CSigCharged.h  CYrast.h
CChart.h      CHistory.h       CNucleus.h  CRunThick.h  CTlArray.h     SStoreEvap.h
CEvap.h       CLevelDensity.h  CNuclide.h  CScission.h  CTlBarDist.h)

add_library(gemini++ SHARED ${SOURCES})
target_compile_definitions(gemini++ PUBLIC GINPUT="${CMAKE_INSTALL_FULL_PKGDATADIR}/")

install(TARGETS gemini++ DESTINATION ${CMAKE_INSTALL_PKLIBDIR})
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PKGINCDIR})
install(DIRECTORY tbl tl DESTINATION ${CMAKE_INSTALL_PKGDATADIR})
