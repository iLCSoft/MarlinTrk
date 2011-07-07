#ifndef MarlinTrk_Factory_h
#define MarlinTrk_Factory_h

#include "IMarlinTrkSystem.h"
#include "gear/GEAR.h"
#include "gear/GearMgr.h"

namespace MarlinTrk{
  
  /** Simple Factory for creating the MarlinTrkSystem of a certain type:
   *  KalTest, LEPTracking, GenFit,...<br>
   *
   *  For now only KalTest.
   *
   * @author F.Gaede, DESY
   * @version $ID:$
   */
  class Factory {
    
  public:
    
    virtual ~Factory() {};
    
    /** Create the MarlinTrkSystem instance of the specified type:<br>
     *  KalTest, LEPTracking, GenFit,...<br>
     *  Returns null if type not implemented...
     * 
     *  For now only KalTest.
     */
    static IMarlinTrkSystem* createMarlinTrkSystem( const std::string& systemType,  
						    const gear::GearMgr* mgr , 
						    const std::string& options ) ;
    
    
  } ;
  
} // end of MarlinTrk namespace 

#endif

