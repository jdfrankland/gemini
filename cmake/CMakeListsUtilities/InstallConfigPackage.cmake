#---install cmake stuff for find_package for use by other projects
# -- (no version information)

#-- namespace by default = project name
if(NOT PROJECT_TARGETS_NAMESPACE)
    set(PROJECT_TARGETS_NAMESPACE ${PROJECT_NAME})
endif()

install(EXPORT ${PROJECT_NAME}Exports
    DESTINATION ${CMAKE_INSTALL_PKGCONFIGDIR}
    NAMESPACE ${PROJECT_TARGETS_NAMESPACE}::
    FILE ${PROJECT_NAME}-targets.cmake
)

include(CMakePackageConfigHelpers)

#---retrieve any package dependencies set by calls to find_package_add_dependency
get_property(PACKAGE_DEPENDENCIES GLOBAL PROPERTY ${PROJECT_NAME}_PACKAGE_DEPENDENCIES)

#---retrieve any package dependencies set by calls to add_subproject
get_property(SUBPROJECT_DEPENDENCIES GLOBAL PROPERTY ${PROJECT_NAME}_SUBPROJECT_DEPENDENCIES)

#-- config template file by default: ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in
if(NOT PROJECT_CONFIG_TEMPLATE_FILE)
    set(PROJECT_CONFIG_TEMPLATE_FILE ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in)
endif()

configure_package_config_file(${PROJECT_CONFIG_TEMPLATE_FILE}
                              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
                              INSTALL_DESTINATION ${CMAKE_INSTALL_PKGCONFIGDIR}
                              )

install(FILES
       ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
     DESTINATION
       ${CMAKE_INSTALL_PKGCONFIGDIR}
     )
