#include "MarlinTrk/Factory.h"

#include "MarlinTrk/MarlinKalTest.h"
#include "MarlinTrk/MarlinDDKalTest.h"

#include "streamlog/streamlog.h"

#include <sstream>

namespace MarlinTrk{
  
  

  IMarlinTrkSystem*  Factory::createMarlinTrkSystem(const std::string& systemType,  
						    const gear::GearMgr* gearMgr,
						    const std::string& options ){
    
    
    // check if we have already instantiated a tracking system of the requested type:
    
    TrkSystemMap::iterator tsI = instance()->_map.find( systemType ) ;

    if( tsI != instance()->_map.end() ) {

      streamlog_out( DEBUG4 ) << " Factory::createMarlinTrkSystem(): return already created IMarlinTrkSystem "
			      << " of type: " << systemType << std::endl ;

      return tsI->second ;
    }
    

    //--------------------------    
    
    IMarlinTrkSystem* trkSystem = 0 ;
    
    streamlog_out(  MESSAGE ) << " Factory::createMarlinTrkSystem:  creating IMarlinTrkSystem of type \"" 
			      << systemType << "\""  << std::endl ;
    

    if( systemType == std::string( "KalTest" ) ) {
      
      trkSystem = new MarlinKalTest( *gearMgr ) ;
      
    } else if( systemType == std::string( "DDKalTest" ) ) {
      
      trkSystem = new MarlinDDKalTest ;
    }
      
    if( ! trkSystem ) {
      
      std::stringstream log ;
      log << " Factory::createMarlinTrkSystem - cannot create IMarlinTrkSystem for type : " << systemType ;
      throw Exception( log.str() ) ;
    }
    
    instance()->_map.insert( std::make_pair( systemType, trkSystem ) ) ;

    instance()->_currentTrkSystem  = trkSystem ;

    return instance()->_currentTrkSystem ;
  }

  //-------------------------------------------------------------------------------------------------------------------
  
  IMarlinTrkSystem* Factory::getMarlinTrkSystem(const std::string& systemType) {  
    
    TrkSystemMap::iterator tsI = instance()->_map.find( systemType ) ;

    if( tsI == instance()->_map.end() ) {

      std::stringstream log ;
      log << " Factory::getMarlinTrkSystem called without a preceeding call to createMarlinTrkSystem() for type : " 
	  << systemType << std::endl ;

      throw Exception( log.str() ) ;
    }

    streamlog_out( DEBUG4 ) << " Factory::getMarlinTrkSystem(): return IMarlinTrkSystem "
			      << " of type: " << systemType << std::endl ;

    instance()->_currentTrkSystem  = tsI->second ;

    return  instance()->_currentTrkSystem ;
  }
    
  //-------------------------------------------------------------------------------------------------------------------
  IMarlinTrkSystem* Factory::getCurrentMarlinTrkSystem() {  
    
    IMarlinTrkSystem* current = instance()->_currentTrkSystem  ;

    if( current == 0  ){
    
      std::stringstream log ;
      log << " Factory::getCurrentMarlinTrkSystem called without a preceeding call to createMarlinTrkSystem() ot getMarlinTrkSystem() " ;

      throw Exception( log.str() ) ;
    } 

    streamlog_out( DEBUG4 ) << " Factory::getCurrentMarlinTrkSystem() called - return allready initialized IMarlinTrkSystem " << std::endl ;
    
    return current ; 
  }
  //-------------------------------------------------------------------------------------------------------------------

  Factory* Factory::instance() {
    static Factory _me ;    
    return &_me ;
  }
 //-------------------------------------------------------------------------------------------------------------------








}
