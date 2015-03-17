#include "MarlinTrk/Factory.h"

#include "MarlinTrk/MarlinKalTest.h"
#include "MarlinTrk/MarlinDDKalTest.h"

#include "streamlog/streamlog.h"

#include <sstream>

namespace MarlinTrk{
  
  
  IMarlinTrkSystem*  Factory::createMarlinTrkSystem(const std::string& systemType,  
						    const gear::GearMgr* gearMgr,
						    const std::string& options ){
    
    
    if( systemType == std::string( "KalTest" ) ) {
      
      streamlog_out( DEBUG5 ) << " create IMarlinTrkSystem of type \"KalTest\"" << std::endl ;
      
      return new MarlinKalTest( *gearMgr ) ;
      
      
    } else if( systemType == std::string( "DDKalTest" ) ) {
      
      static MarlinDDKalTest* ddkaltest_system = 0 ;

      if( ! ddkaltest_system ) {

	streamlog_out(  MESSAGE ) << " Factory::createMarlinTrkSystem:  creating IMarlinTrkSystem of type \"DDKalTest\"" << std::endl ;

	ddkaltest_system = new MarlinDDKalTest ;
      }
      
      return ddkaltest_system ;
      //    return new MarlinDDKalTest ;
      

    }else{
      
      std::stringstream log ;
      log << " Factory::createMarlinTrkSystem - cannot create IMarlinTrkSystem for type : " << systemType ;
      throw Exception( log.str() ) ;
      
    }

    return 0 ;
  }
  
  
  
}
