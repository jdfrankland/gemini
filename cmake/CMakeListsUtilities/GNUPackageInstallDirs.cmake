# for CMake versions < 3.21
include(CheckProjectIsTopLevel)

if(PROJECT_IS_TOP_LEVEL)
#-- if this a stand-alone installation, define our installataion layout
   include(GNUInstallDirs)
   set(CMAKE_INSTALL_PKGINCDIR ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGLIBDIR ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGCONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
else()
#-- when built as a subproject, we define our installation relative to the main project
   set(CMAKE_INSTALL_PKGINCDIR ${CMAKE_INSTALL_PKGINCDIR}/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGLIBDIR ${CMAKE_INSTALL_PKGLIBDIR}/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGCONFIGDIR ${CMAKE_INSTALL_PKGCONFIGDIR}/${PROJECT_NAME})
endif()

