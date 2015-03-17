#include "MarlinTrk/Factory.h"

#include "MarlinTrk/MarlinKalTest.h"
#include "MarlinTrk/MarlinDDKalTest.h"

#include "streamlog/streamlog.h"

#include <sstream>

namespace MarlinTrk{
  
  

  IMarlinTrkSystem*  Factory::createMarlinTrkSystem(const std::string& systemType,  
						    const gear::GearMgr* gearMgr,
						    const std::string& options ){
    
    
    std::string& type = instance()->_myTrkSystemName  ;
    
    IMarlinTrkSystem* current = instance()->_currentTrkSystem  ;
    
    // check if we have already instantiated a different tracking system:
    
    if( current != 0 &&  type != systemType ){
      
      std::stringstream log ;
      log << " Factory::createMarlinTrkSystem - cannot create IMarlinTrkSystem for type : " << systemType  << "   already previously created with type : " << type << std::endl ;
      throw Exception( log.str() ) ;
    }
    
    //--------------------------    
    
    
    if( systemType == std::string( "KalTest" ) ) {
      
      static MarlinKalTest* kaltest_system = 0 ;
      
      if( ! kaltest_system ) {
	
	streamlog_out(  MESSAGE ) << " Factory::createMarlinTrkSystem:  creating IMarlinTrkSystem of type \"KalTest\"" << std::endl ;
	
	kaltest_system = new MarlinKalTest( *gearMgr ) ;
      }
      
      instance()->_currentTrkSystem = kaltest_system ;
      instance()->_myTrkSystemName = systemType ;

      return kaltest_system ;
      
    } else if( systemType == std::string( "DDKalTest" ) ) {
      
      static MarlinDDKalTest* ddkaltest_system = 0 ;

      if( ! ddkaltest_system ) {

	streamlog_out(  MESSAGE ) << " Factory::createMarlinTrkSystem:  creating IMarlinTrkSystem of type \"DDKalTest\"" << std::endl ;

	ddkaltest_system = new MarlinDDKalTest ;
      }
      
      instance()->_currentTrkSystem = ddkaltest_system ;
      instance()->_myTrkSystemName = systemType ;

      return ddkaltest_system ;
      
    } else {
      
      std::stringstream log ;
      log << " Factory::createMarlinTrkSystem - cannot create IMarlinTrkSystem for type : " << systemType ;
      throw Exception( log.str() ) ;
      
    }

    return 0 ;
  }

  //-------------------------------------------------------------------------------------------------------------------

  
  IMarlinTrkSystem* Factory::getCurrentMarlinTrkSystem() {  
    
    const std::string&  type = instance()->_myTrkSystemName  ;

    IMarlinTrkSystem* current = instance()->_currentTrkSystem  ;

    if( current == 0  ){
      
      std::stringstream log ;
      log << " Factory::getCurrentMarlinTrkSystem called without a preceeding call to createMarlinTrkSystem(const std::string& systemType, const gear::GearMgr* gearMgr, const std::string& options ) " ;

      throw Exception( log.str() ) ;
    } 
    
    streamlog_out( DEBUG6 ) << " Factory::getCurrentMarlinTrkSystem() called - return allready initialized IMarlinTrkSystem of type : " << type << std::endl ;

    return current ; 
  }

  //-------------------------------------------------------------------------------------------------------------------

  Factory* Factory::instance() {
    static Factory _me ;    
    return &_me ;
  }
 //-------------------------------------------------------------------------------------------------------------------








}
