#ifndef MarlinTrk_Factory_h
#define MarlinTrk_Factory_h

#include "IMarlinTrkSystem.h"
#include "gear/GEAR.h"
#include "gear/GearMgr.h"

#include <string>
#include <map>

namespace MarlinTrk{
  
  /** Factory methods for creating the MarlinTrkSystem of a certain type:
   *  KalTest, DDKalTest, aidaTT,...<br>
   *  Currently implemented: KalTest, DDKalTest.
   *  The returned instance for a given type is cached, thus: <br>
   * 
   *  DO NOT DELETE THE POINTER at the end of your software module (Marlin processor) ! 
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
     * 
     *  DO NOT DELETE THE POINTER at the end of your software module (Marlin processor) ! 
     * 
     *  For now only KalTest or DDKalTest.
     */
    static IMarlinTrkSystem* createMarlinTrkSystem(const std::string& systemType,  
                                                   const gear::GearMgr* gearMgr,
                                                   const std::string& options ) ;
    
    
    /** Return the MarlinTrkSystem of the given type - only valid after a preceeding call
     *  to  createMarlinTrkSystem() for that type, otherwise an exception is thrown.
     *  This is useful for several modules (e.g. Marlin processors) using the same 
     *  IMarlinTrkSystem. It is the users responsibility to make sure one module
     *  has created it.
     * 
     *  DO NOT DELETE THE POINTER at the end of your software module (Marlin processor) ! 
     * 
     */
    static IMarlinTrkSystem* getMarlinTrkSystem(const std::string& systemType) ;


    /** Return the current MarlinTrkSystem, i.e. the one returned in the last
     *  call to createMarlinTrkSystem() or getMarlinTrkSystem() preceeding this
     *  call. To be used in cases where the concrete type of the tracking system
     *  does not matter - or when it is savely initialized in the enclosing
     *  software module (e.e. the Marlin processor).
     * 
     *  DO NOT DELETE THE POINTER at the end of your software module (Marlin processor) ! 
     * 
     */
    static IMarlinTrkSystem* getCurrentMarlinTrkSystem() ;


    static Factory* instance() ;

  protected:

    typedef std::map< std::string, IMarlinTrkSystem*> TrkSystemMap ;

    Factory() : _currentTrkSystem(0){}

    IMarlinTrkSystem* _currentTrkSystem ;

    TrkSystemMap _map ;

  } ;
  
} // end of MarlinTrk namespace 

#endif

