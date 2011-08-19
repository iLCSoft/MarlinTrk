#include "MarlinTrk/MarlinKalTest.h"
#include "MarlinTrk/MarlinKalTestTrack.h"

#include "kaltest/TKalDetCradle.h"
#include "kaltest/TVKalDetector.h"
#include "kaltest/THelicalTrack.h"

#include "kaldet/ILDVMeasLayer.h"

#include "kaldet/ILDSupportKalDetector.h"
#include "kaldet/ILDVXDKalDetector.h"
//#include "kaldet/ILDSITKalDetector.h"
#include "kaldet/ILDTPCKalDetector.h"

//SJA:FIXME: only needed for storing the modules in the layers map
#include <UTIL/BitField64.h>
#include "UTIL/ILDConf.h"

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"

#include <math.h>
#include <cmath>

#include <utility>

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
  //  now store the measurement layer id's for the active layers 
  //  this->storeActiveMeasurementLayerIDs(ipdet);

  ILDSupportKalDetector* supportdet = new ILDSupportKalDetector( *_gearMgr )  ;  

  ILDVXDKalDetector* vxddet = new ILDVXDKalDetector( *_gearMgr )  ;
  // now store the measurement layer id's for the active layers 
  this->storeActiveMeasurementModuleIDs(vxddet);
  
//  ILDSITKalDetector* sitdet = new ILDSITKalDetector( *_gearMgr )  ;
//  // now store the measurement layer id's for the active layers 
//  this->storeActiveMeasurementLayerIDs(sitdet);

  streamlog_out( DEBUG4 ) << "  MarlinKalTest - create an ILDTPCKalDetector " << std::endl ;

  ILDTPCKalDetector* tpcdet = new ILDTPCKalDetector( *_gearMgr )  ;
  // now store the measurement layer id's for the active layers 
  this->storeActiveMeasurementModuleIDs(tpcdet);


  _det->Install( *supportdet ) ;  
  _det->Install( *vxddet ) ;  
//  _det->Install( *sitdet ) ;    
  _det->Install( *tpcdet ) ;  

  
  _det->Close() ;          // close the cradle
  _det->Sort() ;           // sort meas. layers from inside to outside

  streamlog_out( DEBUG4 ) << "  MarlinKalTest - number of layers = " << _det->GetEntriesFast() << std::endl ;

//  //  int ilayer =  _det->GetEntriesFast()-1 ;
//  int ilayer = -1 ;
//  const ILDVMeasLayer   &ml  = *dynamic_cast<ILDVMeasLayer *>(_det->At(ilayer)); 
//  streamlog_out( DEBUG4 ) << "  MarlinKalTest - name of last layer = " << ml.GetMLName() << std::endl ;

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

void MarlinKalTest::getSensitiveMeasurementModulesForLayer( int layerID, std::vector<ILDVMeasLayer*>& measmodules){

  if( ! measmodules.empty() ) {
    
    std::stringstream errorMsg;
    errorMsg << "MarlinKalTest::getSensitiveMeasurementModulesForLayer vector passed as second argument is not empty " << std::endl ; 
    throw MarlinTrk::Exception(errorMsg.str());
    
  }
  
  std::pair<std::multimap<Int_t, ILDVMeasLayer*>::iterator, std::multimap<Int_t, ILDVMeasLayer*>::iterator> ii;
  ii = this->_active_measurement_modules_by_layer.equal_range(layerID); // set the first and last entry in ii;

  std::multimap<Int_t, ILDVMeasLayer*>::iterator it; //Iterator to be used along with ii


  for(it = ii.first; it != ii.second; ++it)
    {
      //      std::cout<<"Key = "<<it->first<<"    Value = "<<it->second << std::endl ;
      measmodules.push_back( it->second ) ; 
    }
      
}

void MarlinKalTest::getSensitiveMeasurementModules( int moduleID , std::vector<ILDVMeasLayer*>& measmodules ){

  if( ! measmodules.empty() ) {
    
    std::stringstream errorMsg;
    errorMsg << "MarlinKalTest::getSensitiveMeasurementLayer vector passed as second argument is not empty " << std::endl ; 
    throw MarlinTrk::Exception(errorMsg.str());

  }

  std::pair<std::multimap<Int_t, ILDVMeasLayer*>::iterator, std::multimap<Int_t, ILDVMeasLayer*>::iterator> ii;
  ii = this->_active_measurement_modules.equal_range(moduleID); // set the first and last entry in ii;

  std::multimap<Int_t, ILDVMeasLayer*>::iterator it; //Iterator to be used along with ii


  for(it = ii.first; it != ii.second; ++it)
    {
      //      std::cout<<"Key = "<<it->first<<"    Value = "<<it->second << std::endl ;
      measmodules.push_back( it->second ) ; 
    }
      
}


void MarlinKalTest::storeActiveMeasurementModuleIDs(TVKalDetector* detector){
  
  Int_t nLayers = detector->GetEntriesFast() ;
  
  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  
  for( int i=0; i < nLayers; ++i ){
    
    ILDVMeasLayer* ml = dynamic_cast<ILDVMeasLayer*>( detector->At( i ) ); 

    if( ! ml ) {
      std::stringstream errorMsg;
      errorMsg << "MarlinKalTest::storeActiveMeasurementLayerIDs  dynamic_cast to ILDVMeasLayer* failed " << std::endl ; 
      throw MarlinTrk::Exception(errorMsg.str());
    }

    if( ml->IsActive() ) {

      this->_active_measurement_modules.insert(std::pair<int,ILDVMeasLayer*>(ml->getLayerID(),ml));
      
      
      encoder.setValue( ml->getLayerID() );

      // set the module field value to 0 leaving only sub_det, side and layer set.
      encoder[ILDCellID0::module] = 0 ;
      encoder[ILDCellID0::sensor] = 0 ;

      int subdet_layer_id = encoder.lowWord() ;

      this->_active_measurement_modules_by_layer.insert(std::pair<int,ILDVMeasLayer*>(subdet_layer_id,ml));
      
      streamlog_out(DEBUG) << "MarlinKalTest::storeActiveMeasurementLayerIDs added active layer with "
			   << " ModuleID = " << ml->getLayerID()
			   << " LayerID = " << subdet_layer_id
			   << std::endl ;

      

    }

  }
  
}

const ILDVMeasLayer*  MarlinKalTest::getLastMeasLayer(THelicalTrack const& hel, TVector3 const& point) {

  THelicalTrack helix = hel;

  TMatrixD covK(5,5) ;
  double deflection_to_point = 0 ;
  helix.MoveTo(  point, deflection_to_point , 0 , 0) ;

  bool isfwd = ((helix.GetKappa() > 0 && deflection_to_point < 0) || (helix.GetKappa() <= 0 && deflection_to_point > 0)) ? true : false;
  
  int mode = isfwd ? -1 : +1 ;

//  streamlog_out( DEBUG4 ) << "  MarlinKalTest - getLastMeasLayer deflection to point = " << deflection_to_point << " kappa = " << helix.GetKappa()  << "  mode = " << mode << std::endl ;
//  streamlog_out( DEBUG4 ) << " Point to move to:" << std::endl;
//  point.Print();

  int nsufaces =  _det->GetEntriesFast();

  const ILDVMeasLayer* ml_retval = NULL;
  double min_deflection = DBL_MAX;

  for(int i=0; i<nsufaces; ++i){
   
    const ILDVMeasLayer   &ml  = *dynamic_cast<ILDVMeasLayer *>(_det->At(i)); 

    double defection_angle = 0 ;
    TVector3 crossing_point ;   

    const TVSurface *sfp = dynamic_cast<const TVSurface *>(&ml);  // surface at destination       


    int does_cross = sfp->CalcXingPointWith(helix, crossing_point, defection_angle, mode) ;

    if( does_cross ) {

      const double deflection = fabs( deflection_to_point - defection_angle ) ;

      if( deflection < min_deflection ) {

//	streamlog_out( DEBUG4 ) << "  MarlinKalTest - crossing found for suface = " << ml.GetMLName() 
//				<< std::endl
//				<< "  min_deflection = " << min_deflection
//				<< "  deflection = " << deflection
//				<< "  deflection angle = " << defection_angle 
//				<< std::endl 
//				<< " x = " << crossing_point.X() 
//				<< " y = " << crossing_point.Y() 
//				<< " z = " << crossing_point.Z() 
//				<< " r = " << crossing_point.Perp() 
//				<< std::endl ;
	
	min_deflection = deflection ;
	ml_retval = &ml ;
      }

     
    }

  }
  
  return ml_retval;
}
