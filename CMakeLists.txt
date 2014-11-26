########################################################
# cmake file for building MarlinTrk
# - based on $MARLIN/example/mymarlin/CMakeLists.txt
#   by @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( MarlinTrk )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 11)
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

cmake_policy(SET CMP0008 NEW)

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( LCIO REQUIRED )
INCLUDE_DIRECTORIES( ${LCIO_INCLUDE_DIRS} )
LINK_LIBRARIES( ${LCIO_LIBRARIES} )
ADD_DEFINITIONS( ${LCIO_DEFINITIONS} )

FIND_PACKAGE( GEAR REQUIRED )
INCLUDE_DIRECTORIES( ${GEAR_INCLUDE_DIRS} )
LINK_LIBRARIES( ${GEAR_LIBRARIES} )
ADD_DEFINITIONS( ${GEAR_DEFINITIONS} )

FIND_PACKAGE( streamlog REQUIRED )
FIND_PACKAGE( streamlog REQUIRED )
INCLUDE_DIRECTORIES( ${streamlog_INCLUDE_DIRS} )
LINK_LIBRARIES( ${streamlog_LIBRARIES} )
ADD_DEFINITIONS( ${streamlog_DEFINITIONS} )

FIND_PACKAGE( KalTest REQUIRED ) 
INCLUDE_DIRECTORIES( ${KalTest_INCLUDE_DIRS} )
LINK_LIBRARIES( ${KalTest_LIBRARIES} )
ADD_DEFINITIONS( ${KalTest_DEFINITIONS} )

FIND_PACKAGE( KalDet REQUIRED ) 
INCLUDE_DIRECTORIES( ${KalDet_INCLUDE_DIRS} )
LINK_LIBRARIES( ${KalDet_LIBRARIES} )
ADD_DEFINITIONS( ${KalDet_DEFINITIONS} )

FIND_PACKAGE( ROOT REQUIRED ) 
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

FIND_PACKAGE( DDKalTest REQUIRED ) 
INCLUDE_DIRECTORIES( ${DDKalTest_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DDKalTest_LIBRARIES} )
ADD_DEFINITIONS( ${DDKalTest_DEFINITIONS} )

FIND_PACKAGE( DD4hep REQUIRED ) 
INCLUDE_DIRECTORIES( ${DD4hep_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DD4hep_LIBRARIES} )
ADD_DEFINITIONS( ${DD4hep_DEFINITIONS} )


FIND_PACKAGE( CLHEP REQUIRED ) 
INCLUDE_DIRECTORIES( ${CLHEP_INCLUDE_DIRS} )
LINK_LIBRARIES( ${CLHEP_LIBRARIES} )
ADD_DEFINITIONS( ${CLHEP_DEFINITIONS} )



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( BEFORE ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )

# needed for adding header files to xcode project
IF(CMAKE_GENERATOR MATCHES "Xcode")
  FILE( GLOB_RECURSE library_headers "*.h" )
  ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} ${library_headers})
ELSE()
  ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
ENDIF()                

INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( MarlinTrkConfig.cmake MarlinTrkConfigVersion.cmake MarlinTrkLibDeps.cmake )

