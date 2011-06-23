#include "Factory.h"

#include "MarlinKalTest.h"



namespace MarlinTrk{


  IMarlinTrkSystem*  Factory::createMarlinTrkSystem( const std::string& systemType,  
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
