#ifndef MarlinTrk_Factory_h
#define MarlinTrk_Factory_h

#include "IMarlinTrkSystem.h"
#include "gear/GEAR.h"
#include "gear/GearMgr.h"

#include <string>

namespace MarlinTrk{
  
  /** Simple Factory for creating the MarlinTrkSystem of a certain type:
   *  KalTest, LEPTracking, GenFit,...<br>
   *
   *  For now only KalTest.
   *
   * @author F.Gaede, DESY
   * @version $Id$
   */
  class Factory {
    
  public:
    
    virtual ~Factory() {}
    
    /** Create the MarlinTrkSystem instance of the specified type:<br>
     *  KalTest, DDKalTest, aidaTT,...<br>
     *  Returns 0 if type not implemented...
     * 
     *  For now only KalTest and DDKalTest.
     */
    static IMarlinTrkSystem* createMarlinTrkSystem(const std::string& systemType,  
                                                   const gear::GearMgr* gearMgr,
                                                   const std::string& options ) ;
    
    
    /** Return the current MarlinTrkSystem - only valid after a preceeding call
     *  to  createMarlinTrkSystem(), otherwise an exception is thrown.
     *  This is useful for several modules (e.g. Marlin processors) using the same 
     *  IMarlinTrkSystem. It is the users responsibility to make sure one module
     *  has created it.
     */
    static IMarlinTrkSystem* getCurrentMarlinTrkSystem() ;


    static Factory* instance() ;

  protected:

    Factory() : _currentTrkSystem(0) , _myTrkSystemName("UNKNOWN") {}

    IMarlinTrkSystem* _currentTrkSystem ;

    std::string _myTrkSystemName ;
    
  } ;
  
} // end of MarlinTrk namespace 

#endif

