cmake_minimum_required (VERSION 2.8)

project(Gemini)

set(SOURCES Angle.cpp      Fus.cpp           LightP.cpp   Random.cpp      TlArray.cpp
AngleDist.cpp  Gdr.cpp           Mass.cpp     Scission.cpp    TlBarDist.cpp
Chart.cpp      History.cpp       Nucleus.cpp  SigBarDist.cpp  Weight.cpp
Evap.cpp       LevelDensity.cpp  Nuclide.cpp  SigCharged.cpp  Yrast.cpp)

set(HEADERS CAngleDist.h  CFus.h           CLightP.h   CRandom.h    CSigBarDist.h  CWeight.h
CAngle.h      CGdr.h           CMass.h     CRun.h       CSigCharged.h  CYrast.h
CChart.h      CHistory.h       CNucleus.h  CRunThick.h  CTlArray.h     SStoreEvap.h
CEvap.h       CLevelDensity.h  CNuclide.h  CScission.h  CTlBarDist.h)

add_library(Gemini SHARED ${SOURCES})

set(LIB_INSTALL ${CMAKE_INSTALL_PREFIX}/lib)
set(HEAD_INSTALL ${CMAKE_INSTALL_PREFIX}/include/gemini)
set(DATA_INSTALL ${CMAKE_INSTALL_PREFIX}/share/gemini)

install(TARGETS Gemini DESTINATION ${LIB_INSTALL})
install(FILES ${HEADERS} DESTINATION ${HEAD_INSTALL})
install(DIRECTORY tbl tl DESTINATION ${DATA_INSTALL})
