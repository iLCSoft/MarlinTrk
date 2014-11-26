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
      
      streamlog_out(  DEBUG5 ) << " create IMarlinTrkSystem of type \"DDKalTest\"" << std::endl ;
      
      return new MarlinDDKalTest ;
      

    }else{
      
      std::stringstream log ;
      log << " Factory::createMarlinTrkSystem - cannot create IMarlinTrkSystem for type : " << systemType ;
      throw Exception( log.str() ) ;
      
    }

    return 0 ;
  }
  
  
  
}
