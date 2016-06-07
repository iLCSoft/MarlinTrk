#include "MarlinTrk/MarlinDDKalTest.h"
#include "MarlinTrk/MarlinDDKalTestTrack.h"

#include "kaltest/TKalDetCradle.h"
#include "kaltest/TVKalDetector.h"
#include "kaltest/THelicalTrack.h"

#include "DDKalTest/DDVMeasLayer.h"
#include "DDKalTest/DDKalDetector.h"
#include "DDKalTest/DDCylinderMeasLayer.h"

// #include "DDKalTest/DDSupportKalDetector.h"
// #include "DDKalTest/DDVXDKalDetector.h"

//SJA:FIXME: only needed for storing the modules in the layers map
#include <UTIL/BitField64.h>
#include "UTIL/ILDConf.h"

#include "DD4hep/LCDD.h"
#include "DDRec/SurfaceManager.h"

#include <algorithm>
#include <string>
#include <math.h>
#include <cmath>
#include <fstream>

#include <utility>

#include "streamlog/streamlog.h"

//#include "DDKalTest/DDMeasurementSurfaceStoreFiller.h"

namespace MarlinTrk{
  
  
  MarlinDDKalTest::MarlinDDKalTest()  :
    _ipLayer(NULL) {
    
    streamlog_out( DEBUG4 ) << "  MarlinDDKalTest - initializing the detector ..." << std::endl ;
    
    _det = new TKalDetCradle ;  // from kaltest. TKalDetCradle inherits from TObjArray ... 
    _det->SetOwner( true ) ;    // takes care of deleting subdetector in the end ...
    
    is_initialised = false; 
    
    this->registerOptions() ;
    
    streamlog_out( DEBUG4 ) << "  MarlinDDKalTest - established " << std::endl ;
  }

  MarlinDDKalTest::~MarlinDDKalTest(){
    
#ifdef MARLINTRK_DIAGNOSTICS_ON
    _diagnostics.end();
#endif
    
    delete _det ;
  }
  
  
  void MarlinDDKalTest::init() {
    
     
    this->includeMultipleScattering( getOption(IMarlinTrkSystem::CFG::useQMS) ) ;  

    this->includeEnergyLoss( getOption(IMarlinTrkSystem::CFG::usedEdx) ) ; 

    streamlog_out( DEBUG5 ) << " -------------------------------------------------------------------------------- " << std::endl ;
    streamlog_out( DEBUG5 ) << "  MarlinDDKalTest::init() called with the following options :                     " << std::endl ;
    streamlog_out( DEBUG5 ) <<    this->getOptions() ;
    streamlog_out( DEBUG5 ) << " -------------------------------------------------------------------------------- " << std::endl ;
    
    if( is_initialised ) {
      
      streamlog_out( DEBUG5 ) << "  MarlinDDKalTest::init()  - already initialized - only options are set .. " << std::endl ;
      
      return ;
    }
    

    streamlog_out( DEBUG5 ) << " ##################### MarlinDDKalTest::init()  - initializing  " << std::endl ;

    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();

    double minS = 1.e99; 
    DDCylinderMeasLayer* ipLayer = 0 ;


    // for the tracking we get all tracking detectors and all passive detectors (beam pipe,...)

    std::vector< DD4hep::Geometry::DetElement>        detectors   = lcdd.detectors( "tracker" ) ;
    const std::vector< DD4hep::Geometry::DetElement>& passiveDets = lcdd.detectors( "passive" ) ;
    const std::vector< DD4hep::Geometry::DetElement>& calos       = lcdd.detectors( "calorimeter" ) ;

    detectors.reserve( detectors.size() + passiveDets.size() + calos.size() ) ;

    std::copy( passiveDets.begin() , passiveDets.end() , std::back_inserter( detectors )  ) ;

    for ( std::vector< DD4hep::Geometry::DetElement>::const_iterator it=calos.begin() ; it != calos.end() ; ++it ){

    // for ( std::vector< DD4hep::Geometry::DetElement>::const_iterator it=passiveDets.begin() ; it != passiveDets.end() ; ++it ){

      std::string name = it->name() ;
      std::transform( name.begin() , name.end() , name.begin() , ::tolower ) ;
      if( name.find( "ecal" ) != std::string::npos ){

	detectors.push_back( *it ) ;
      }
    }  



    for ( std::vector< DD4hep::Geometry::DetElement>::iterator it=detectors.begin() ; it != detectors.end() ; ++it ){

      DD4hep::Geometry::DetElement det = *it ;

      streamlog_out( DEBUG5 ) << "  MarlinDDKalTest::init() - creating DDKalDetector for : " << det.name() << std::endl ;

      DDKalDetector* kalDet = new DDKalDetector( det ) ;

      this->storeActiveMeasurementModuleIDs( kalDet ) ;

      _det->Install( *kalDet ) ;


      Int_t nLayers = kalDet->GetEntriesFast() ;
    

      // --- keep the cylinder meas layer with smallest sorting policy (radius) as ipLayer
      // fixme: this should be implemented in a more explicit way ...
      for( int i=0; i < nLayers; ++i ) {
	const TVSurface* tvs = dynamic_cast<const TVSurface*>( kalDet->At( i ) ); 

	double s = tvs->GetSortingPolicy() ;
	if( s < minS &&  dynamic_cast< DDCylinderMeasLayer* > (  kalDet->At( i) )  ) {
	  minS = s  ;
	  ipLayer = dynamic_cast< DDCylinderMeasLayer* > (  kalDet->At( i) ) ;
	}
      }

      if( streamlog_level( DEBUG5 ) ) {   // dump surfaces to text file 
	
	std::map< double, DDVMeasLayer*> smap ;
	std::ofstream file ;
	std::stringstream s ; s << "DDKalTest_" <<  det.name() << "_surfaces.txt" ;
	file.open( s.str() , std::ofstream::out  ) ; 
	lcio::BitField64 bf(  UTIL::ILDCellID0::encoder_string ) ;
	
	for( unsigned i=0,N=kalDet->GetEntriesFast() ; i<N ;++i){
	  DDVMeasLayer* ml = dynamic_cast<DDVMeasLayer*> ( kalDet->At( i ) ) ;
	  TVSurface* s =  dynamic_cast<TVSurface*> ( kalDet->At( i ) ) ;
	  smap[ s->GetSortingPolicy() ] = ml ;
	}
	for( std::map<double,DDVMeasLayer*>::iterator it=smap.begin() ; it!=smap.end() ; ++it){
	  bf.setValue( it->second->getCellIDs()[0] ) ;
	  file << " "  <<  std::scientific << std::setw(10) << it->first  <<  "\t" << bf.valueString() << *it->second->surface()  << "\n"  ;
	}
	file.close() ;
      }

    }



    if( ipLayer) {

      _ipLayer = ipLayer ;

      streamlog_out( MESSAGE ) << " MarlinDDKalTest: install IP layer at radius : " << minS << std::endl ;
    }
    //-------------------------------------------------------------------------------


    _det->Close() ;          // close the cradle
    //done in Close()    _det->Sort() ;           // sort meas. layers from inside to outside
    
    streamlog_out( DEBUG4 ) << "  MarlinDDKalTest - number of layers = " << _det->GetEntriesFast() << std::endl ;
           

    if( streamlog_level( DEBUG ) ) {

      lcio::BitField64 bf(  UTIL::ILDCellID0::encoder_string ) ;
      
      for( unsigned i=0,N=_det->GetEntriesFast() ; i<N ;++i){

	DDVMeasLayer* ml = dynamic_cast<DDVMeasLayer*> ( _det->At( i ) ) ;
	
	bf.setValue( ml->getLayerID() ) ;

	TVSurface* s =  dynamic_cast<TVSurface*> ( _det->At( i ) ) ;

	streamlog_out( DEBUG ) << " *** meas. layer : " << bf.valueString() << "  sorting: " <<  s->GetSortingPolicy()  << std::endl ;
      }

    }

    is_initialised = true; 
    
  }
  
  MarlinTrk::IMarlinTrack* MarlinDDKalTest::createTrack()  {
    
    if ( ! is_initialised ) {
      
      std::stringstream errorMsg;
      errorMsg << "MarlinDDKalTest::createTrack: Fitter not initialised. MarlinDDKalTest::init() must be called before MarlinDDKalTest::createTrack()" << std::endl ; 
      throw MarlinTrk::Exception(errorMsg.str());
      
    }
    
    return new MarlinDDKalTestTrack(this) ;
    
  }
  
  void MarlinDDKalTest::includeMultipleScattering( bool msOn ) {
    
    if( msOn == true ) {
      _det->SwitchOnMS();
    }
    else{
      _det->SwitchOffMS();
    } 
    
  } 
  
  void MarlinDDKalTest::includeEnergyLoss( bool energyLossOn ) {
    
    if( energyLossOn == true ) {
      _det->SwitchOnDEDX();
    }
    else{
      _det->SwitchOffDEDX();
    } 
    
  } 
  
  void MarlinDDKalTest::getSensitiveMeasurementModulesForLayer( int layerID, std::vector< const DDVMeasLayer *>& measmodules) const {
    
    if( ! measmodules.empty() ) {
      
      std::stringstream errorMsg;
      errorMsg << "MarlinDDKalTest::getSensitiveMeasurementModulesForLayer vector passed as second argument is not empty " << std::endl ; 
      throw MarlinTrk::Exception(errorMsg.str());
      
    }
    
    streamlog_out( DEBUG0 ) << "MarlinDDKalTest::getSensitiveMeasurementModulesForLayer: layerID = " << layerID << std::endl;
    
    std::multimap<Int_t, const DDVMeasLayer *>::const_iterator it; //Iterator to be used along with ii
    
    
    
    //  for(it = _active_measurement_modules_by_layer.begin(); it != _active_measurement_modules_by_layer.end(); ++it) {
    //    streamlog_out( DEBUG0 ) << "Key = "<< ttdecodeILD(it->first) <<"    Value = "<<it->second << std::endl ;
    //  }
    
    
    std::pair<std::multimap<Int_t, const DDVMeasLayer *>::const_iterator, std::multimap<Int_t, const DDVMeasLayer *>::const_iterator> ii;  
    
    // set the module and sensor bit ranges to zero as these are not used in the map 
    lcio::BitField64 bf(  UTIL::ILDCellID0::encoder_string ) ;
    bf.setValue( layerID ) ;
    bf[lcio::ILDCellID0::module] = 0 ;
    bf[lcio::ILDCellID0::sensor] = 0 ;
    layerID = bf.lowWord();
    
    ii = this->_active_measurement_modules_by_layer.equal_range(layerID); // set the first and last entry in ii;
    
    for(it = ii.first; it != ii.second; ++it) {
      //    streamlog_out( DEBUG0 ) <<"Key = "<< it->first <<"    Value = "<<it->second << std::endl ;
      measmodules.push_back( it->second ) ; 
    }
    
  }
  
  void MarlinDDKalTest::getSensitiveMeasurementModules( int moduleID , std::vector< const DDVMeasLayer *>& measmodules ) const {
    
    if( ! measmodules.empty() ) {
      
      std::stringstream errorMsg;
      errorMsg << "MarlinDDKalTest::getSensitiveMeasurementLayer vector passed as second argument is not empty " << std::endl ; 
      throw MarlinTrk::Exception(errorMsg.str());
      
    }
    
    std::pair<std::multimap<int, const DDVMeasLayer *>::const_iterator, std::multimap<Int_t, const DDVMeasLayer *>::const_iterator> ii;
    ii = this->_active_measurement_modules.equal_range(moduleID); // set the first and last entry in ii;
    
    std::multimap<int,const DDVMeasLayer *>::const_iterator it; //Iterator to be used along with ii
    
    
    for(it = ii.first; it != ii.second; ++it) {
      //      std::cout<<"Key = "<<it->first<<"    Value = "<<it->second << std::endl ;
      measmodules.push_back( it->second ) ; 
    }
  }
  
  
  void MarlinDDKalTest::storeActiveMeasurementModuleIDs(TVKalDetector* detector) {
    
    Int_t nLayers = detector->GetEntriesFast() ;
    
    for( int i=0; i < nLayers; ++i ) {
      
      const DDVMeasLayer* ml = dynamic_cast<const DDVMeasLayer*>( detector->At( i ) ); 
      
      if( ! ml ) {
        std::stringstream errorMsg;
        errorMsg << "MarlinDDKalTest::storeActiveMeasurementLayerIDs dynamic_cast to DDVMeasLayer* failed " << std::endl ; 
        throw MarlinTrk::Exception(errorMsg.str());
      }
      
      if( ml->IsActive() ) {
        
        // then get all the sensitive element id's assosiated with this DDVMeasLayer and store them in the map 
        std::vector<int>::const_iterator it = ml->getCellIDs().begin();
        
        while ( it!=ml->getCellIDs().end() ) {
          
          int sensitive_element_id = *it;
          this->_active_measurement_modules.insert(std::pair<int,const DDVMeasLayer*>( sensitive_element_id, ml ));        
          ++it;
          
        }
        
        int subdet_layer_id = ml->getLayerID() ;
        
        this->_active_measurement_modules_by_layer.insert(std::pair<int ,const DDVMeasLayer*>(subdet_layer_id,ml));
        
        streamlog_out(DEBUG0) << "MarlinDDKalTest::storeActiveMeasurementLayerIDs added active layer with "
        << " LayerID = " << subdet_layer_id << " and DetElementIDs  " ;
        
        for (it = ml->getCellIDs().begin(); it!=ml->getCellIDs().end(); ++it) {
          
          streamlog_out(DEBUG0) << " : " << *it ;
          
        }
        
        streamlog_out(DEBUG0) << std::endl;
        
        
        
        
      }
      
    }
    
  }
  
  const DDVMeasLayer*  MarlinDDKalTest::getLastMeasLayer(THelicalTrack const& hel, TVector3 const& point) const {
    
    THelicalTrack helix = hel;
    
    double deflection_to_point = 0 ;
    helix.MoveTo(  point, deflection_to_point , 0 , 0) ;
    
    bool isfwd = ((helix.GetKappa() > 0 && deflection_to_point < 0) || (helix.GetKappa() <= 0 && deflection_to_point > 0)) ? true : false;
    
    int mode = isfwd ? -1 : +1 ;
    
    //  streamlog_out( DEBUG4 ) << "  MarlinDDKalTest - getLastMeasLayer deflection to point = " << deflection_to_point << " kappa = " << helix.GetKappa()  << "  mode = " << mode << std::endl ;
    //  streamlog_out( DEBUG4 ) << " Point to move to:" << std::endl;
    //  point.Print();
    
    int nsufaces =  _det->GetEntriesFast();
    
    const DDVMeasLayer* ml_retval = 0;
    double min_deflection = DBL_MAX;
    
    for(int i=0; i<nsufaces; ++i) {
      
      const DDVMeasLayer   &ml  = *dynamic_cast< const DDVMeasLayer *>(_det->At(i)); 
      
      double defection_angle = 0 ;
      TVector3 crossing_point ;   
      
      const TVSurface *sfp = dynamic_cast<const TVSurface *>(&ml);  // surface at destination       
      
      
      int does_cross = sfp->CalcXingPointWith(helix, crossing_point, defection_angle, mode) ;
      
      if( does_cross ) {
        
        const double deflection = fabs( deflection_to_point - defection_angle ) ;
        
        if( deflection < min_deflection ) {
          
          //      streamlog_out( DEBUG4 ) << "  MarlinDDKalTest - crossing found for suface = " << ml.GetMLName() 
          //                              << std::endl
          //                              << "  min_deflection = " << min_deflection
          //                              << "  deflection = " << deflection
          //                              << "  deflection angle = " << defection_angle 
          //                              << std::endl 
          //                              << " x = " << crossing_point.X() 
          //                              << " y = " << crossing_point.Y() 
          //                              << " z = " << crossing_point.Z() 
          //                              << " r = " << crossing_point.Perp() 
          //                              << std::endl ;
          
          min_deflection = deflection ;
          ml_retval = &ml ;
        }
        
      }
      
    }
    
    return ml_retval;
  }
  
  const DDVMeasLayer* MarlinDDKalTest::findMeasLayer( EVENT::TrackerHit * trkhit) const {
    
    const TVector3 hit_pos( trkhit->getPosition()[0], trkhit->getPosition()[1], trkhit->getPosition()[2]) ;
    
    return this->findMeasLayer( trkhit->getCellID0(), hit_pos ) ;
    
  }
  
  const DDVMeasLayer* MarlinDDKalTest::findMeasLayer( int detElementID, const TVector3& point) const {
    
    const DDVMeasLayer* ml = 0; // return value 
    
    std::vector<const DDVMeasLayer*> meas_modules ;
    
    // search for the list of measurement layers associated with this CellID
    this->getSensitiveMeasurementModules( detElementID, meas_modules ) ; 
    
    if( meas_modules.size() == 0 ) { // no measurement layers found 
      
      UTIL::BitField64 encoder( UTIL::ILDCellID0::encoder_string ) ; 
      encoder.setValue(detElementID) ;
      
      std::stringstream errorMsg;
      errorMsg << "MarlinDDKalTest::findMeasLayer module id unkown: moduleID = " << detElementID 
	       << " [" << encoder.valueString() << "]" << std::endl ; 
       throw MarlinTrk::Exception(errorMsg.str());
      
    } 
    else if (meas_modules.size() == 1) { // one to one mapping 
      
      ml = meas_modules[0] ;
      
    }
    else { // layer has been split 
      
      bool surf_found(false);
      
      // loop over the measurement layers associated with this CellID and find the correct one using the position of the hit
      for( unsigned int i=0; i < meas_modules.size(); ++i) {
        
        
        
        const TVSurface* surf = 0;
        
        if( ! (surf = dynamic_cast<const TVSurface*> (  meas_modules[i] )) ) {
          std::stringstream errorMsg;
          errorMsg << "MarlinDDKalTest::findMeasLayer dynamic_cast failed for surface type: moduleID = " << detElementID << std::endl ; 
          throw MarlinTrk::Exception(errorMsg.str());
        }
        
        bool hit_on_surface = surf->IsOnSurface(point);
        
        if( (!surf_found) && hit_on_surface ){
          
          ml = meas_modules[i] ;
          surf_found = true ;
          
        }
        else if( surf_found && hit_on_surface ) {  // only one surface should be found, if not throw 
          
          std::stringstream errorMsg;
          errorMsg << "MarlinDDKalTest::findMeasLayer point found to be on two surfaces: moduleID = " << detElementID << std::endl ; 
          throw MarlinTrk::Exception(errorMsg.str());
        }      
        
      }
      if( ! surf_found ){ // print out debug info
        streamlog_out(DEBUG1) << "MarlinDDKalTest::findMeasLayer point not found to be on any surface matching moduleID = "
        << detElementID
        << ": x = " << point.x()
        << " y = " << point.y()
        << " z = " << point.z()
        << std::endl ;
      }
      else{
        streamlog_out(DEBUG1) << "MarlinDDKalTest::findMeasLayer point found to be on surface matching moduleID = "
        << detElementID
        << ": x = " << point.x()
        << " y = " << point.y()
        << " z = " << point.z()
        << std::endl ;
      }
    }
    
    return ml ;
    
  }
  
} // end of namespace MarlinTrk
