#include "MarlinTrk/MarlinAidaTT.h"
#include "MarlinTrk/MarlinAidaTTTrack.h"

// // //SJA:FIXME: only needed for storing the modules in the layers map
// #include <UTIL/BitField64.h>
// #include "UTIL/ILDConf.h"

#include "aidaTT/IBField.hh"
#include "aidaTT/ConstantSolenoidBField.hh"
#include "aidaTT/analyticalPropagation.hh"
#include "aidaTT/simplifiedPropagation.hh"
#include "aidaTT/GBLInterface.hh"
#include "aidaTT/fitResults.hh"
#include "aidaTT/Vector5.hh"
#include "aidaTT/utilities.hh"
#include "aidaTT/LCIOPersistency.hh"
#include "aidaTT/Vector3D.hh"
#include "aidaTT/IGeometry.hh"


#include "DD4hep/LCDD.h"
#include "DD4hep/Fields.h"
#include "DD4hep/DD4hepUnits.h"

#include <algorithm>
#include <string>
#include <math.h>
#include <cmath>

#include <utility>

#include "streamlog/streamlog.h"

namespace MarlinTrk{
  
  
  MarlinAidaTT::MarlinAidaTT() : _useQMS(false), _usedEdx(false) ,_is_initialised(false){
    
    this->registerOptions() ;
    
    streamlog_out( DEBUG4 ) << "  MarlinAidaTT - constructed " << std::endl ;
  }
  
  MarlinAidaTT::~MarlinAidaTT(){
    delete  _geom ;
    delete  _bfield ;
    delete  _fitter ;
    delete  _propagation  ;
  }
  
  
  void MarlinAidaTT::init() {
    
    _useQMS =  getOption( IMarlinTrkSystem::CFG::useQMS )  ;  
    
    _usedEdx =  getOption( IMarlinTrkSystem::CFG::usedEdx ) ; 
    
    streamlog_out( DEBUG5 ) << " -------------------------------------------------------------------------------- " << std::endl ;
    streamlog_out( DEBUG5 ) << "  MarlinAidaTT::init() called with the following options :                        " << std::endl ;
    streamlog_out( DEBUG5 ) <<    this->getOptions() ;
    streamlog_out( DEBUG5 ) << " -------------------------------------------------------------------------------- " << std::endl ;
    
    if( _is_initialised ) {
      
      streamlog_out( DEBUG5 ) << "  MarlinAidaTT::init()  - already initialized - only options are set .. " << std::endl ;
      
      return ;
    }
    
   
    streamlog_out( DEBUG5 ) << " ##################### MarlinAidaTT::init()  - initializing  " << std::endl ;
    
    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
    
    
    double origin[3] = { 0., 0., 0. }, bfield[3] ;
    lcdd.field().magneticField( origin , bfield  ) ;

    _bfield = new aidaTT::ConstantSolenoidBField( bfield[2] / dd4hep::tesla ) ; 

    _propagation = new aidaTT::analyticalPropagation();
    //_propagation = new aidaTT::simplifiedPropagation();


    _fitter = new aidaTT::GBLInterface();

    
    // get the all surfaces in the detector
    DD4hep::Geometry::DetElement world = lcdd.world() ;
    
    _geom = & aidaTT::IGeometry::instance() ;
    
    const std::vector<const aidaTT::ISurface*>& surfaces = _geom->getSurfaces() ;
    
    for(std::vector<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf) {
      _surfMap.insert( std::make_pair( (*surf)->id(), *surf ) ) ;
    }
    
    streamlog_out( DEBUG5 ) << "  MarlinAidaTT - number of surfaces = " << _surfMap.size() << std::endl ;

    _is_initialised = true; 

  }
  
  MarlinTrk::IMarlinTrack* MarlinAidaTT::createTrack()  {
    
    if ( ! _is_initialised ) {
      
      std::stringstream errorMsg;
      
      errorMsg << "MarlinAidaTT::createTrack: Fitter not initialised. MarlinAidaTT::init() must be called before MarlinAidaTT::createTrack()" << std::endl ; 
      throw MarlinTrk::Exception(errorMsg.str());
      
    }
    return new MarlinAidaTTTrack(this) ;
  
  }
  
  
  
} // end of namespace MarlinTrk
