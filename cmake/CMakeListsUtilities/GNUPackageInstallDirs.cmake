cmake_minimum_required(VERSION 3.7) # for GNUInstallDirs_get_absolute_install_dir macro

# for CMake versions < 3.21
include(CheckProjectIsTopLevel)

if(PROJECT_IS_TOP_LEVEL)
#-- if this a stand-alone installation, define our installataion layout
   include(GNUInstallDirs)
   set(CMAKE_INSTALL_PKGINCDIR ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGLIBDIR ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGCONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGDATADIR ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME})
else()
#-- when built as a subproject, we define our installation relative to the main project
   set(CMAKE_INSTALL_PKGINCDIR ${CMAKE_INSTALL_PKGINCDIR}/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGLIBDIR ${CMAKE_INSTALL_PKGLIBDIR}/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGCONFIGDIR ${CMAKE_INSTALL_PKGCONFIGDIR}/${PROJECT_NAME})
   set(CMAKE_INSTALL_PKGDATADIR ${CMAKE_INSTALL_PKGDATADIR}/${PROJECT_NAME})
endif()

#-- make absolute path variables using macro from GNUInstallDirs >=3.7
foreach(dir
   PKGINCDIR
   PKGLIBDIR
   PKGCONFIGDIR
   PKGDATADIR
   )
   if(CMAKE_VERSION VERSION_LESS 3.20)
      GNUInstallDirs_get_absolute_install_dir(CMAKE_INSTALL_FULL_${dir} CMAKE_INSTALL_${dir})
   else()
      GNUInstallDirs_get_absolute_install_dir(CMAKE_INSTALL_FULL_${dir} CMAKE_INSTALL_${dir} ${dir})
   endif()
endforeach()
 
