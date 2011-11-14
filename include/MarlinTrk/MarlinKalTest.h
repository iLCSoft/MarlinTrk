#ifndef INCLUDE_MarlinKalTest
#define INCLUDE_MarlinKalTest 1

#include "MarlinTrk/IMarlinTrkSystem.h"


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
  void getSensitiveMeasurementModules( int detElementID, std::vector<ILDVMeasLayer*>& measmodules);
  
  /** Store active measurement module IDs needed for navigation  */
  void getSensitiveMeasurementModulesForLayer( int layerID, std::vector<ILDVMeasLayer*>& measmodules);
  
  //  void init(bool MSOn, bool EnergyLossOn) ;
  bool is_initialised ;
  
  //** find the measurment layer for a given hit 
  const ILDVMeasLayer* findMeasLayer( EVENT::TrackerHit * trkhit) ;
  //** find the measurment layer for a given det element ID and point in space 
  const ILDVMeasLayer* findMeasLayer( int detElementID, const TVector3& point) ;
  
  // get the last layer crossed by the helix when extrapolating from the present position to the pca to point
  const ILDVMeasLayer* getLastMeasLayer(THelicalTrack const& helix, TVector3 const& point) ;
  
  const ILDCylinderMeasLayer* getIPLayer() { return _ipLayer; }
  
  const ILDCylinderMeasLayer* _ipLayer ;
  
  const gear::GearMgr* _gearMgr ;
  
  TKalDetCradle* _det ;            // the detector cradle
  
  std::multimap< Int_t, ILDVMeasLayer*> _active_measurement_modules;
  
  std::multimap< Int_t, ILDVMeasLayer*> _active_measurement_modules_by_layer;
  
  
  
} ;

#endif
