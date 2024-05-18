CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#===============================================================================
# Load macros and functions for Elements-based projects
#===============================================================================

find_package(ElementsProject)

#===============================================================================
# Declare project name and version
# Example with dependency:
#                         elements_project(MyProject 1.0 USE Element 3.9)
#===============================================================================

elements_project(EXT_PF1_dm4lsst 0.2 USE Elements 5.10.0 ST_DataModelBindings 8.0.3 ST_DataModelTools 8.0.3)
