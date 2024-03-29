########################################################
# cmake file for building MarlinTrk
# - based on $MARLIN/example/mymarlin/CMakeLists.txt
#   by @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( MarlinTrk )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 2 )
SET( ${PROJECT_NAME}_VERSION_MINOR 9)
SET( ${PROJECT_NAME}_VERSION_PATCH 2 )

cmake_policy(SET CMP0008 NEW)

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


OPTION( MARLINTRK_USE_GEAR "Set to OFF to not provide Gear backward compatibility " ON)
SET( MARLINTRK_USE_GEAR_CONFIG_VARIABLE 0 )
#--------------------------------------------


FIND_PACKAGE( LCIO REQUIRED )
INCLUDE_DIRECTORIES( SYSTEM ${LCIO_INCLUDE_DIRS} )
LINK_LIBRARIES( ${LCIO_LIBRARIES} )
ADD_DEFINITIONS( ${LCIO_DEFINITIONS} )

IF( MARLINTRK_USE_GEAR)
  FIND_PACKAGE( GEAR REQUIRED )
  INCLUDE_DIRECTORIES( SYSTEM  ${GEAR_INCLUDE_DIRS} )
  LINK_LIBRARIES( ${GEAR_LIBRARIES} )
  ADD_DEFINITIONS( ${GEAR_DEFINITIONS} )
  ADD_DEFINITIONS( -D MARLINTRK_BACKWARD_GEAR_WRAPPERS )
  SET( MARLINTRK_USE_GEAR_CONFIG_VARIABLE 1 )
ENDIF(MARLINTRK_USE_GEAR)


FIND_PACKAGE( streamlog REQUIRED )
INCLUDE_DIRECTORIES( SYSTEM  ${streamlog_INCLUDE_DIRS} )
LINK_LIBRARIES( ${streamlog_LIBRARIES} )
ADD_DEFINITIONS( ${streamlog_DEFINITIONS} )

FIND_PACKAGE( KalTest REQUIRED ) 
INCLUDE_DIRECTORIES( SYSTEM  ${KalTest_INCLUDE_DIRS} )
LINK_LIBRARIES( ${KalTest_LIBRARIES} )
ADD_DEFINITIONS( ${KalTest_DEFINITIONS} )

#---- this set kSdim=5 in KalTest ---
# ( should be exported by KalTest really )
ADD_DEFINITIONS( -D __NOT0__ )

FIND_PACKAGE( KalDet REQUIRED ) 
INCLUDE_DIRECTORIES( SYSTEM  ${KalDet_INCLUDE_DIRS} )
LINK_LIBRARIES( ${KalDet_LIBRARIES} )
ADD_DEFINITIONS( ${KalDet_DEFINITIONS} )

FIND_PACKAGE( ROOT REQUIRED ) 
INCLUDE_DIRECTORIES( SYSTEM  ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

FIND_PACKAGE( DDKalTest REQUIRED ) 
INCLUDE_DIRECTORIES( SYSTEM  ${DDKalTest_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DDKalTest_LIBRARIES} )
ADD_DEFINITIONS( ${DDKalTest_DEFINITIONS} )

# DD4hep depends on geant4. geant4 calls find_package clhep internally and also contains a findclhep.cmake (before 10.5)
# that findclhep cmake is not playing nice with other things so we call find clhep before find dd4hep
FIND_PACKAGE(CLHEP REQUIRED)
INCLUDE_DIRECTORIES( SYSTEM  ${CLHEP_INCLUDE_DIRS} )
LINK_LIBRARIES( ${CLHEP_LIBRARIES} )
ADD_DEFINITIONS( ${CLHEP_DEFINITIONS} )

FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDRec ) 
INCLUDE_DIRECTORIES( SYSTEM  ${DD4hep_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DD4hep_LIBRARIES}  ${DD4hep_COMPONENT_LIBRARIES} )
ADD_DEFINITIONS( ${DD4hep_DEFINITIONS} )

FIND_PACKAGE( aidaTT REQUIRED )
ADD_DEFINITIONS( "-D AIDATT_USE_DD4HEP -DUSE_GBL -DUSE_LCIO" ) 
INCLUDE_DIRECTORIES( SYSTEM  ${aidaTT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${aidaTT_LIBRARIES} )
ADD_DEFINITIONS( ${aidaTT_DEFINITIONS} )

FIND_PACKAGE( GSL REQUIRED)
INCLUDE_DIRECTORIES( SYSTEM  ${GSL_INCLUDE_DIRS} )
LINK_LIBRARIES( ${GSL_LIBRARIES} )
ADD_DEFINITIONS( ${GSL_DEFINITIONS} )

FIND_PACKAGE( GBL REQUIRED)
INCLUDE_DIRECTORIES( SYSTEM  ${GBL_INCLUDE_DIRS} )
LINK_LIBRARIES( ${GBL_LIBRARIES} )
ADD_DEFINITIONS( ${GBL_DEFINITIONS} )


### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
#ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#ADD_DEFINITIONS( "-Wno-long-long" )

ADD_DEFINITIONS( " -UNDEBUG " )

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
GENERATE_PACKAGE_CONFIGURATION_FILES( MarlinTrkConfig.cmake MarlinTrkConfigVersion.cmake )

