#include "MarlinTrk/Factory.h"

#include "MarlinTrk/MarlinKalTest.h"



namespace MarlinTrk{
  
  
  IMarlinTrkSystem*  Factory::createMarlinTrkSystem(const std::string& systemType,  
                                                    const gear::GearMgr* mgr , 
                                                    const std::string& options ){
    
    
    if( systemType == std::string( "KalTest" ) )
      return new MarlinKalTest( *mgr ) ;
    
    
    //FIXME:    add other implemetations here ....
    //
    //   ....
    
    return 0 ;
    
  }
  
  
  
}
