cmake_minimum_required(VERSION 3.7) # for GNUInstallDirs_get_absolute_install_dir macro

# for CMake versions < 3.21
include(CheckProjectIsTopLevel)

# by default, we lower-case the project name for installation directories
if(NOT DO_NOT_LOWER_CASE_PROJECT_INSTALL_DIRS)
   string(TOLOWER ${PROJECT_NAME} proj_name)
endif()

if(PROJECT_IS_TOP_LEVEL)
#-- if this a stand-alone installation, define our installataion layout
   include(GNUInstallDirs)
   set(CMAKE_INSTALL_PKGINCDIR ${CMAKE_INSTALL_INCLUDEDIR}/${proj_name})
   set(CMAKE_INSTALL_PKGLIBDIR ${CMAKE_INSTALL_LIBDIR}/${proj_name})
   set(CMAKE_INSTALL_PKGCONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${proj_name})
   set(CMAKE_INSTALL_PKGDATADIR ${CMAKE_INSTALL_DATADIR}/${proj_name})
else()
#-- when built as a subproject, we define our installation relative to the main project
   set(CMAKE_INSTALL_PKGINCDIR ${CMAKE_INSTALL_PKGINCDIR}/${proj_name})
   set(CMAKE_INSTALL_PKGLIBDIR ${CMAKE_INSTALL_PKGLIBDIR}/${proj_name})
   set(CMAKE_INSTALL_PKGCONFIGDIR ${CMAKE_INSTALL_PKGCONFIGDIR}/${proj_name})
   set(CMAKE_INSTALL_PKGDATADIR ${CMAKE_INSTALL_PKGDATADIR}/${proj_name})
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
 
