#ifndef MarlinAidaTT_h
#define MarlinAidaTT_h

#include "MarlinTrk/IMarlinTrkSystem.h"

#include "DDRec/SurfaceManager.h"

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

#include "aidaTT/AidaTT.hh"
#include "aidaTT/IGeometry.hh"

namespace EVENT{
  class TrackerHit ;
}

namespace aidaTT{
  class IGeometry ;
  class IBField ;
  class IFittingAlgorithm ;
  class IPropagation ;
}

namespace MarlinTrk{
  
  typedef std::multimap< long, const aidaTT::ISurface* > SurfMap ;


  /** Interface to KaltTest Kalman fitter - instantiates and holds the detector geometry.
   */
  class MarlinAidaTT : public MarlinTrk::IMarlinTrkSystem {
    
  public:
    
    friend class MarlinAidaTTTrack;
    
    // define some configuration constants
    static const bool FitBackward   = kIterBackward ;
    static const bool FitForward    = kIterForward ;
    static const bool OrderOutgoing  = true ;
    static const bool OrderIncoming  = false ;
    
    /// Default c'tor
    MarlinAidaTT() ;
    
    /** d'tor */
    ~MarlinAidaTT() ;
    
    /** initialise track fitter system */
    void init() ; 
    
   /// the name of the implementation 
    virtual std::string name() { return "AidaTT" ; }

    /** instantiate its implementation of the IMarlinTrack */
    MarlinTrk::IMarlinTrack* createTrack()  ;
    
    
  protected:
    
    // /** take multiple scattering into account during the fit */
    // void includeMultipleScattering( bool on )  ;
    
    // /** take energy loss into account during the fit */
    // void includeEnergyLoss( bool on )  ;
    
    // /** Store active measurement module IDs for a given TVKalDetector needed for navigation  */
    // void storeActiveMeasurementModuleIDs(TVKalDetector* detector);  
    
    // /** Store active measurement module IDs needed for navigation  */
    // void getSensitiveMeasurementModules( int detElementID, std::vector< const DDVMeasLayer *>& measmodules) const; 
    
    // /** Store active measurement module IDs needed for navigation  */
    // void getSensitiveMeasurementModulesForLayer( int layerID, std::vector<const DDVMeasLayer *>& measmodules) const;
    
    // //  void init(bool MSOn, bool EnergyLossOn) ;
    // bool is_initialised ;
    
    // //** find the measurment layer for a given hit 
    // const DDVMeasLayer* findMeasLayer( EVENT::TrackerHit * trkhit) const ; 
    // //** find the measurment layer for a given det element ID and point in space 
    // const DDVMeasLayer* findMeasLayer( int detElementID, const TVector3& point) const ;
    
    // // get the last layer crossed by the helix when extrapolating from the present position to the pca to point
    // const DDVMeasLayer* getLastMeasLayer(THelicalTrack const& helix, TVector3 const& point) const ;
    
    // std::multimap< int,const DDVMeasLayer *> _active_measurement_modules;
    // std::multimap< int,const DDVMeasLayer *> _active_measurement_modules_by_layer;
        
    
    bool _useQMS ;
    bool _usedEdx ;
    bool _is_initialised ;
    
    /// multi-map of surfaces
    SurfMap _surfMap ;
    
    const aidaTT::IGeometry*         _geom ;
    aidaTT::IBField*           _bfield ;
    aidaTT::IFittingAlgorithm* _fitter ;
    aidaTT::IPropagation*      _propagation  ;

 } ;
  
} // end of namespace MarlinTrk

#endif
