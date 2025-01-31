# CMakeListsUtilities

Various useful bits of code used for writing CMakeLists.txt files:

  - CheckProjectIsTopLevel: sets PROJECT_IS_TOP_LEVEL even if CMake version < 3.21
  - GNUPackageInstallDirs: extension of GNUInstallDirs which adds project-specific directories for libraries, header files, data files, and CMake package config files
  - StandardRpathHandling: make sure all libraries/executables have correct RPATH settings
  - InstallConfigPackage: installs Config files so other packages can find us (no version information)
  - InstallConfigPackageVersion: installs Config files so other packages can find us (with version information)
  - projectConfig.cmake.in: template to be used as base for package Config files

## StandardRpathHandling

Usage:
~~~{.cmake}
include(StandardRpathHandling)
~~~

  + Ensures the RPATH/RUNPATH for all executables & libraries is set to find our package's libraries at runtime, as well as any 3rd party libraries on which they depend, without needing to modify `LD_LIBRARY_PATH`
  + `CMAKE_INSTALL_PKGLIBDIR` must be set beforehand: it should contain the path (relative to `CMAKE_INSTALL_PREFIX`) which contains our package's libraries (e.g. `lib` or `lib/packageName`)

## InstallConfigPackage, InstallConfigPackageVersion

Usage:
~~~{.cmake}
include(InstallConfigPackage)
~~~
or 
~~~{.cmake}
include(InstallConfigPackageVersion)
~~~

  + Typically used at the end of the top level CMakeLists.txt in order to generate and install the `projectNameConfig.cmake` file which will allow others to find and use our package with a `find_package(projectName)` command
  + use `InstallConfigPackageVersion` if your project has a version number, otherwise use `InstallConfigPackage`
  + The project's exported targets will be written in a projectName-targets.cmake file, make sure that:
    - when exporting targets, use `${PROJECT_NAME}Exports` as the name of the target set;
    - if required, set `PROJECT_TARGETS_NAMESPACE` to the name that will be prefixed to targets (as `${PROJECT_TARGETS_NAMESPACE}::[target]`). By default, `${PROJECT_NAME}` will be used, unless the variable NO_PROJECT_TARGETS_NAMESPACE is set;;
    - if required, set `PROJECT_CONFIG_TEMPLATE_FILE` to contain the path to the template Config file (see projectConfig.cmake.in below). By default, we look for a file `${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in`;
  + set `CMAKE_INSTALL_PKGCONFIGDIR` beforehand to the relative path where all Config files will be installed (typicaly, `lib/cmake/projectName`)
  + any calls to `find_dependency` generated by `FindPackageAddDependency` (see https://gitlab.in2p3.fr/jdfcode/cmake/findpackageadddependency.git) will be automatically added to the Config script, as long as the template file contains the necessary code (see projectConfig.cmake.in below).
  + for any subprojects added with `add_subproject` (see https://gitlab.in2p3.fr/jdfcode/cmake/addsubproject) the list of subproject RPATH additions is retrieved and will be automatically added to the Config script, as long as the template file contains the necessary code (see projectConfig.cmake.in below).

## projectConfig.cmake.in

Template file to be used with `configure_package_config_file()` or with `InstallConfigPackage[Version]` (see above).

Some manual intervention may be required depending on the project, but it should work (minimally) as is.
  + any calls to `find_dependency` generated by `FindPackageAddDependency` (see https://gitlab.in2p3.fr/jdfcode/cmake/findpackageadddependency.git) will be automatically added to the Config script
  + all subprojects added with `add_subproject` (see https://gitlab.in2p3.fr/jdfcode/cmake/addsubproject) will have the paths to their libraries added to the install RPATH
  + any library or executable linking to our targets gets a full RPATH to allow it to find all the required shared libs

## GNUPackageInstallDirs
Just an extension of GNUInstallDirs which adds project-specific directories for installing
libraries, header files, data files and CMake package config files. The result depends on
whether the current project is being built as a top-level project or a subproject of
another:

~~~{.cmake}
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
~~~

PROJECT_NAME must be set before calling. By default, the project name is converted to all lower case to be used in the directory names, unless variable DO_NOT_LOWER_CASE_PROJECT_INSTALL_DIRS is set.

If being used as a subproject, the master project must have defined the CMAKE_INSTALL_* variables
(e.g. by itself calling `include(GNUPackageInstallDirs)`)
