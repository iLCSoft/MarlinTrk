#ifndef MarlinKalTest_h
#define MarlinKalTest_h

#include "MarlinTrk/IMarlinTrkSystem.h"

#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif

#include "gear/GearMgr.h"

//LCIO:
#include "lcio.h"
#include "UTIL/BitField64.h" 
#include "UTIL/LCTOOLS.h"
#include <LCRTRelations.h>

#include "streamlog/streamlog.h"

#include "TObjArray.h"
#include "TVector3.h"

#include <cmath>
#include <vector>


class TKalDetCradle ;
class TVKalDetector ;
class ILDVMeasLayer ;
class THelicalTrack ;

class ILDCylinderMeasLayer;

namespace EVENT{
  class TrackerHit ;
}

namespace MarlinTrk{
  
  /** Interface to KaltTest Kalman fitter - instantiates and holds the detector geometry.
   */
  class MarlinKalTest : public MarlinTrk::IMarlinTrkSystem {
    
  public:
    
    friend class MarlinKalTestTrack;
    
    // define some configuration constants
    static const bool FitBackward   = kIterBackward ;
    static const bool FitForward    = kIterForward ;
    static const bool OrderOutgoing  = true ;
    static const bool OrderIncoming  = false ;
    
    
    /** Default c'tor, initializes the geometry from GEAR. */
    MarlinKalTest( const gear::GearMgr& gearMgr) ;
    
    /** d'tor */
    ~MarlinKalTest() ;
    
    /** initialise track fitter system */
    void init() ; 
    
    /** instantiate its implementation of the IMarlinTrack */
    MarlinTrk::IMarlinTrack* createTrack()  ;
    
    
  protected:
    
    /** take multiple scattering into account during the fit */
    void includeMultipleScattering( bool on )  ;
    
    /** take energy loss into account during the fit */
    void includeEnergyLoss( bool on )  ;
    
    /** Store active measurement module IDs for a given TVKalDetector needed for navigation  */
    void storeActiveMeasurementModuleIDs(TVKalDetector* detector);  
    
    /** Store active measurement module IDs needed for navigation  */
    void getSensitiveMeasurementModules( int detElementID, std::vector< const ILDVMeasLayer *>& measmodules) const; 
    
    /** Store active measurement module IDs needed for navigation  */
    void getSensitiveMeasurementModulesForLayer( int layerID, std::vector<const ILDVMeasLayer *>& measmodules) const;
    
    //  void init(bool MSOn, bool EnergyLossOn) ;
    bool is_initialised ;
    
    //** find the measurment layer for a given hit 
    const ILDVMeasLayer* findMeasLayer( EVENT::TrackerHit * trkhit) const ; 
    //** find the measurment layer for a given det element ID and point in space 
    const ILDVMeasLayer* findMeasLayer( int detElementID, const TVector3& point) const ;
    
    // get the last layer crossed by the helix when extrapolating from the present position to the pca to point
    const ILDVMeasLayer* getLastMeasLayer(THelicalTrack const& helix, TVector3 const& point) const ;
    
    const ILDCylinderMeasLayer* getIPLayer() const { return _ipLayer; }
    
    const ILDCylinderMeasLayer* _ipLayer ;
    
    const gear::GearMgr* _gearMgr ;
    
    TKalDetCradle* _det ;            // the detector cradle
    
    std::multimap< int,const ILDVMeasLayer *> _active_measurement_modules;
    
    std::multimap< int,const ILDVMeasLayer *> _active_measurement_modules_by_layer;
    
    
#ifdef MARLINTRK_DIAGNOSTICS_ON

  private:    
    MarlinTrk::DiagnosticsController _diagnostics;

  public:    

    /** Return the pointer to the Diagnositics Object. Forseen for internal diagnostics, only available when complied with MARLINTRK_DIAGNOSTICS_ON defined. 
     */
    virtual void * getDiagnositicsPointer() { return &_diagnostics ; }
            
#endif

    
  } ;
  
} // end of namespace MarlinTrk

#endif
