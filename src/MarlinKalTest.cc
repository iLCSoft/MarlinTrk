#include "MarlinKalTest.h"

#include "MarlinKalTestTrack.h"

#include "kaltest/TKalDetCradle.h"
#include "kaltest/TVKalDetector.h"

#include "kaldet/ILDVMeasLayer.h"

//#include "kaldet/ILDIPKalDetector.h"
#include "kaldet/ILDVXDKalDetector.h"
//#include "kaldet/ILDSITKalDetector.h"
#include "kaldet/ILDTPCKalDetector.h"

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"

#include <ILDDetectorIDs.h>

#include <math.h>
#include <cmath>

#include "streamlog/streamlog.h"


MarlinKalTest::MarlinKalTest( const gear::GearMgr& gearMgr, bool MSOn, bool EnergyLossOn) :  
  _gearMgr( &gearMgr )  
{
  
  streamlog_out( DEBUG4 ) << "  MarlinKalTest - initializing the detector ..." << std::endl ;
  
  _det = new TKalDetCradle ; // from kaltest. TKalDetCradle inherits from TObjArray ... 
  _det->SetOwner( true ) ; // takes care of deleting subdetector in the end ...
  
  // this could be made a public init() method taking options ....
  streamlog_out( DEBUG4 ) << "  MarlinKalTest - call init " << std::endl ;
  
//  this->includeMultipleScattering(MSOn) ;  
//  this->includeEnergyLoss(EnergyLossOn) ;  
//
//  init() ;
  
  streamlog_out( DEBUG4 ) << "  MarlinKalTest - established " << std::endl ;

}

MarlinKalTest::~MarlinKalTest(){
  
  delete _det ;
}


void MarlinKalTest::init() {
  

//  Double_t ip_radius  = 0.5     ; //mm
//  Double_t ip_halfl   = 2500.0  ; //mm

//  const Double_t bz = _gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;

  //  ILDIPKalDetector* ipdet = new ILDIPKalDetector( ip_radius, ip_halfl, bz )  ; 
  // now store the measurement layer id's for the active layers 
  //  this->storeActiveMeasurementLayerIDs(ipdet);

  ILDVXDKalDetector* vxddet = new ILDVXDKalDetector( *_gearMgr )  ;
  // now store the measurement layer id's for the active layers 
  this->storeActiveMeasurementLayerIDs(vxddet);
  
//  ILDSITKalDetector* sitdet = new ILDSITKalDetector( *_gearMgr )  ;
//  // now store the measurement layer id's for the active layers 
//  this->storeActiveMeasurementLayerIDs(sitdet);

  streamlog_out( DEBUG4 ) << "  MarlinKalTest - create an ILDTPCKalDetector " << std::endl ;

  ILDTPCKalDetector* tpcdet = new ILDTPCKalDetector( *_gearMgr )  ;
  // now store the measurement layer id's for the active layers 
  this->storeActiveMeasurementLayerIDs(tpcdet);


//  _det->Install( *ipdet ) ;  
  _det->Install( *vxddet ) ;  
//  _det->Install( *sitdet ) ;    
  _det->Install( *tpcdet ) ;  

  
  _det->Close() ;          // close the cradle
  _det->Sort() ;           // sort meas. layers from inside to outside
    
}

MarlinTrk::IMarlinTrack* MarlinKalTest::createTrack()  {

  return new MarlinKalTestTrack(this) ;

}

void MarlinKalTest::includeMultipleScattering( bool msOn ) {
  
  if( msOn == true ) {
    _det->SwitchOnMS();
  }
  else{
    _det->SwitchOffMS();
  } 

} 

void MarlinKalTest::includeEnergyLoss( bool energyLossOn ) {
  
  if( energyLossOn == true ) {
    _det->SwitchOnDEDX();
  }
  else{
    _det->SwitchOffDEDX();
  } 

} 



ILDVMeasLayer* MarlinKalTest::getSensitiveMeasurementLayer( Int_t layerID ){
  std::map< Int_t, ILDVMeasLayer*>::iterator it;
  it = this->_active_measurement_layer.find(layerID);
  if( it != this->_active_measurement_layer.end() ){
    return it->second;
  }
  else{
    return NULL; // SJA:FIXME: this should be an exception
  }
    
}


void MarlinKalTest::storeActiveMeasurementLayerIDs(TVKalDetector* detector){
  
  Int_t nLayers = detector->GetEntriesFast() ;
  
  for( int i=0  ; i < nLayers ; ++i ){
    
    ILDVMeasLayer* ml = dynamic_cast<ILDVMeasLayer*>( detector->At( i ) ); 
    
    if( ml->IsActive() ) {
      streamlog_out(DEBUG) << "MarlinKalTest::storeActiveMeasurementLayerIDs added active layer with "
			   << " LayerID = " << ml->getLayerID()
			   << std::endl ;
      this->_active_measurement_layer[ ml->getLayerID() ] = ml;
    }

  }
  
}
