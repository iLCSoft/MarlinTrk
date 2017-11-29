#ifndef MarlinDDKalTest_h
#define MarlinDDKalTest_h

#include "MarlinTrk/IMarlinTrkSystem.h"

#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif

#include "DDKalTest/DDKalDetector.h"

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
class DDVMeasLayer ;
class THelicalTrack ;

class DDCylinderMeasLayer;

namespace EVENT{
  class TrackerHit ;
}

namespace MarlinTrk{
  
  /** Interface to KaltTest Kalman fitter - instantiates and holds the detector geometry.
   */
  class MarlinDDKalTest : public MarlinTrk::IMarlinTrkSystem {
    
  public:
    
    friend class MarlinDDKalTestTrack;
    
    // define some configuration constants
    static const bool FitBackward   = kIterBackward ;
    static const bool FitForward    = kIterForward ;
    static const bool OrderOutgoing  = true ;
    static const bool OrderIncoming  = false ;
    
    /// Default c'tor
    MarlinDDKalTest() ;
    MarlinDDKalTest(const MarlinDDKalTest&) = delete;
    MarlinDDKalTest const& operator=(const MarlinDDKalTest&) = delete;

    /** d'tor */
    ~MarlinDDKalTest() ;
    
    /** Sets the specified option ( one of the constants defined in IMarlinTrkSystem::CFG )
     *  to the given value. Override here to re-configure E-loss and QMS
     *  after the initialization.
     */
    virtual void setOption(unsigned CFGOption, bool val) ;

    /** initialise track fitter system */
    void init() ; 
    
   /// the name of the implementation 
    virtual std::string name() { return "DDKalTest" ; }

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
    void getSensitiveMeasurementModules( int detElementID, std::vector< const DDVMeasLayer *>& measmodules) const; 
    
    /** Store active measurement module IDs needed for navigation  */
    void getSensitiveMeasurementModulesForLayer( int layerID, std::vector<const DDVMeasLayer *>& measmodules) const;
    
    //  void init(bool MSOn, bool EnergyLossOn) ;
    bool is_initialised=false;
    
    //** find the measurment layer for a given hit 
    const DDVMeasLayer* findMeasLayer( EVENT::TrackerHit * trkhit) const ; 
    //** find the measurment layer for a given det element ID and point in space 
    const DDVMeasLayer* findMeasLayer( int detElementID, const TVector3& point) const ;
    
    // get the last layer crossed by the helix when extrapolating from the present position to the pca to point
    const DDVMeasLayer* getLastMeasLayer(THelicalTrack const& helix, TVector3 const& point) const ;
    
    const DDCylinderMeasLayer* getIPLayer() const { return _ipLayer; }
    

    // members:

    const DDCylinderMeasLayer* _ipLayer=nullptr;
    
    TKalDetCradle* _det=nullptr;         // the detector cradle
    
    std::multimap< int,const DDVMeasLayer *> _active_measurement_modules{};
    
    std::multimap< int,const DDVMeasLayer *> _active_measurement_modules_by_layer{};

    std::vector< DDKalDetector* > _detectors{};
    
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
