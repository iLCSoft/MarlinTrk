#include "MarlinTrk/MarlinAidaTTTrack.h"

#include "MarlinTrk/MarlinAidaTT.h"
#include "MarlinTrk/IMarlinTrkSystem.h"


#include <lcio.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>

#include <UTIL/BitField64.h>
#include <UTIL/Operators.h>
#include <UTIL/ILDConf.h>

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
#include "aidaTT/DD4hepGeometry.hh"

#include "DD4hep/DD4hepUnits.h"

#include <sstream>

#include "streamlog/streamlog.h"


namespace MarlinTrk {
  
  //---------------------------------------------------------------------------------------------------------------
  
  namespace{ 
    std::string cellIDString( int detElementID ) {
      lcio::BitField64 bf(  UTIL::ILDCellID0::encoder_string ) ;
      bf.setValue( detElementID ) ;
      return bf.valueString() ;
    }
  }
  //---------------------------------------------------------------------------------------------------------------
  
  
  MarlinAidaTTTrack::MarlinAidaTTTrack( MarlinAidaTT* mAidaTT) 
    : _aidaTT( mAidaTT ) , _initialised( false )  {
    
  }
  
  
  MarlinAidaTTTrack::~MarlinAidaTTTrack(){
    delete _fitTrajectory ;
  }
  
  int MarlinAidaTTTrack::addHit( EVENT::TrackerHit * trkhit) {
    _lcioHits.push_back( trkhit ) ;
    return success ;
  } 
  
  int MarlinAidaTTTrack::initialise( bool fitDirection ) {; 
    
    if ( _initialised ) {
      throw MarlinTrk::Exception("Track fit already initialised");   
    }
    unsigned nHits =  _lcioHits.size() ;
    if( nHits < 3) {
      
      streamlog_out( ERROR) << "<<<<<< MarlinAidaTTTrack::initialise: Shortage of Hits! nhits = "  
			    << _lcioHits.size() << " >>>>>>>" << std::endl;
      return error ;
    }
    
      //--------- get the start helix from three points
    bool backwards = false ;
    
    lcio::TrackerHit* h1 = ( backwards ?  _lcioHits[ nHits-1 ] : _lcioHits[    0    ] ) ;
    lcio::TrackerHit* h2 =  _lcioHits[ (nHits+1) / 2 ] ;
    lcio::TrackerHit* h3 = ( backwards ?  _lcioHits[    0    ] : _lcioHits[ nHits-1 ] ) ;

    const double* pos1 = h1->getPosition() ;
    const double* pos2 = h2->getPosition() ;
    const double* pos3 = h3->getPosition() ;

    aidaTT::Vector3D x1( pos1[0] * dd4hep::mm, pos1[1] * dd4hep::mm , pos1[2] * dd4hep::mm ) ;
    aidaTT::Vector3D x2( pos2[0] * dd4hep::mm, pos2[1] * dd4hep::mm , pos2[2] * dd4hep::mm ) ;
    aidaTT::Vector3D x3( pos3[0] * dd4hep::mm, pos3[1] * dd4hep::mm , pos3[2] * dd4hep::mm ) ;

    calculateStartHelix( x1, x2,  x3 , _initialTrackParams , backwards ) ;
             
    streamlog_out( DEBUG3 )  << "  start helix from three points : " << _initialTrackParams << std::endl ;

    return myInit() ;
  }
  

  int MarlinAidaTTTrack::initialise(  const EVENT::TrackState& ts, double /*bfield_z*/, bool fitDirection ) {

    _initialTrackParams = aidaTT::readLCIO( &ts ) ;

    return myInit() ;
  }

  int MarlinAidaTTTrack::myInit() {
    
    moveHelixTo( _initialTrackParams, aidaTT::Vector3D()  ) ; // move to origin
    
    // --- set some large errors to the covariance matrix ( might not be needed for GBL ) 
    _initialTrackParams.covarianceMatrix().Unit() ;
    _initialTrackParams.covarianceMatrix()( aidaTT::OMEGA, aidaTT::OMEGA ) = 1.e-2 ;
    _initialTrackParams.covarianceMatrix()( aidaTT::TANL , aidaTT::TANL  ) = 1.e2 ;
    _initialTrackParams.covarianceMatrix()( aidaTT::PHI0 , aidaTT::PHI0  ) = 1.e2 ;
    _initialTrackParams.covarianceMatrix()( aidaTT::D0   , aidaTT::D0    ) = 1.e5 ;
    _initialTrackParams.covarianceMatrix()( aidaTT::Z0   , aidaTT::Z0    ) = 1.e5 ;

    _fitTrajectory = new aidaTT::trajectory( _initialTrackParams, _aidaTT->_fitter, _aidaTT->_bfield, 
					     _aidaTT->_propagation, _aidaTT->_geom );

    // add the Interaction Point as the first element of the trajectory
    int ID = 1;
    _fitTrajectory->addElement( aidaTT::Vector3D(), &ID);

    _initialised = true ;

    return success ;
  } 
  
  int MarlinAidaTTTrack::addAndFit( EVENT::TrackerHit* trkhit, double& chi2increment, double maxChi2Increment) {
    
    streamlog_out( DEBUG ) << "MarlinAidaTTTrack::addAndFit: nothing to be done ... " << std::endl ;
    return success ;
  }
  
  int MarlinAidaTTTrack::testChi2Increment( EVENT::TrackerHit* trkhit, double& chi2increment ) {

    streamlog_out( DEBUG ) << "MarlinAidaTTTrack::testChi2Increment: nothing to be done ... " << std::endl ;
    return success ;
  }
  
  int MarlinAidaTTTrack::fit( double maxChi2Increment ) {
    
    streamlog_out(DEBUG2) << "MarlinAidaTTTrack::fit() called " << std::endl ;
    
    if ( ! _initialised ) {
      throw MarlinTrk::Exception("Track fit not initialised");   
    }
    
    unsigned nHits =  _lcioHits.size() ;
    if( nHits < 3) {
      streamlog_out( ERROR) << "<<<<<< MarlinAidaTTTrack::initialise: Shortage of Hits! nhits = "  
			    << _lcioHits.size() << " >>>>>>>" << std::endl;
      return error ;
    }
    
    std::vector<double> precision(2) ;
    
    // ==== loop over hits and prepare the corresponding trajectory elements
    for(unsigned i=0 ; i < nHits ; ++i){
      
      EVENT::TrackerHit* hit = _lcioHits[i] ;
      
      long id = hit->getCellID0() ;
      
      SurfMap::iterator it = _aidaTT->_surfMap.find( id ) ;
      
      if( it == _aidaTT->_surfMap.end() ){
	
	streamlog_out( ERROR ) << "  MarlinAidaTTTrack::fit - no surface found for hit with id : " << cellIDString( id ) 
			       << std::endl ;
	continue ; // ignore hit
      }
      const aidaTT::ISurface* surf = it->second ;
      
      // get the hit position in dd4hep/aidaTT units
      double hitpos[3] ;
      for(unsigned int i = 0; i < 3; ++i) hitpos[i] = hit->getPosition()[i] * dd4hep::mm;
      

      //---- compute the precision from the hit errors

      double du,dv ;

      TrackerHitPlane* planarhit = dynamic_cast<TrackerHitPlane*>( hit );
      if( planarhit != 0 ) {

	du = planarhit->getdU() * dd4hep::mm  ;
	dv = planarhit->getdV() * dd4hep::mm  ;

      } else { // we have a TPC hit which is not yet using the CylinderTrackerHit ...

	const FloatVec& cov = hit->getCovMatrix();
	
	du = sqrt( cov[0] + cov[2] ) * dd4hep::mm  ;
	dv = sqrt( cov[5]          ) * dd4hep::mm  ;
      }

      precision[0] = 1./ (du*du) ;
      precision[1] = 1./ (dv*dv) ;

      
      _fitTrajectory->addMeasurement( hitpos, precision, *surf, hit );
    }
    

    // do the fit:
    _fitTrajectory->prepareForFitting();
    
    int fit_ok = _fitTrajectory->fit();

    //    _fitResult = _fitTrajectory->getFitResults() ;

    return ( fit_ok ? success : error  )  ;
  }
  
  
  /** smooth all track states 
   */
  int MarlinAidaTTTrack::smooth(){

    streamlog_out( DEBUG2 )  << "MarlinAidaTTTrack::smooth() - nothing to do ...." << std::endl ;
    return success ;
  }
  
  
  /** smooth track states from the last filtered hit back to the measurement site associated with the given hit 
   */
  int MarlinAidaTTTrack::smooth( EVENT::TrackerHit* trkhit ) {
    
    streamlog_out( DEBUG2 )  << "MarlinAidaTTTrack::smooth() - nothing to do ...." << std::endl ;
    return success ;
  }
  
  
  int  MarlinAidaTTTrack::getTrackState( IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {
    
    const aidaTT::fitResults& result = *_fitTrajectory->getFitResults();

    ts = *aidaTT::createLCIO( result.estimatedParameters() );
        
    ts.setLocation(lcio::TrackState::AtIP);
    
    chi2 = result.chiSquare() ;

    ndf = result.ndf() ;

    return success ;
  }
  
  
  int MarlinAidaTTTrack::getTrackState( EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {
    
    const aidaTT::fitResults* result = _fitTrajectory->getFitResults( trkhit->getCellID0() );

    if( result == 0 ) 
      return error ;

    ts = *aidaTT::createLCIO( result->estimatedParameters() );
        
    chi2 = result->chiSquare() ;

    ndf = result->ndf() ;

    return success ;
  }
  
  
  int MarlinAidaTTTrack::getHitsInFit( std::vector<std::pair<EVENT::TrackerHit*, double> >& hits ) {
    
    streamlog_out( DEBUG2 )  << "MarlinAidaTTTrack::getHitsInFit() - returning all hits used ..." << std::endl ;
    
    hits.reserve( _lcioHits.size() ) ;

    for(unsigned i=0,N=_lcioHits.size(); i<N ; ++i){

      hits.push_back(  std::make_pair( _lcioHits[i] , 0. ) );
    }
    
    return success ;
  }
  
  int MarlinAidaTTTrack::getOutliers( std::vector<std::pair<EVENT::TrackerHit*, double> >& hits ) {
    streamlog_out( DEBUG2 )  << "MarlinAidaTTTrack::getOutliers() - nothing to do ..." << std::endl ;
    return success ;
  }
  
  
  int MarlinAidaTTTrack::getNDF( int& ndf ){

    const aidaTT::fitResults& result = *_fitTrajectory->getFitResults();
    ndf = result.ndf();
    return success;
  }
  
  
  
  int MarlinAidaTTTrack::getTrackerHitAtPositiveNDF( EVENT::TrackerHit*& trkhit ) {

    trkhit = _lcioHits[0] ; // ????
    return success;    
  }
  
  
  int MarlinAidaTTTrack::extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){  
    
    // const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    // return this->extrapolate( point, site, ts, chi2, ndf ) ;
    return success ;
  }
  
  int MarlinAidaTTTrack::extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {
    
    // TKalTrackSite* site = 0 ;
    // int error_code = getSiteFromLCIOHit(trkhit, site);
    
    // if( error_code != success ) return error_code;
    
    // return this->extrapolate( point, *site, ts, chi2, ndf ) ;
    return success ;
    
  }
  
  
  int MarlinAidaTTTrack::extrapolateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    // const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    // return this->extrapolateToLayer( layerID, site, ts, chi2, ndf, detElementID, mode ) ;
    return success ;
  }


  int MarlinAidaTTTrack::extrapolateToLayer( int layerID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
  
    // TKalTrackSite* site = 0;
    // int error_code = getSiteFromLCIOHit(trkhit, site);
    
    //   if( error_code != success ) return error_code ;
    
    //   return this->extrapolateToLayer( layerID, *site, ts, chi2, ndf, detElementID, mode ) ;
    return success ;
  }
  
  
  
  
  int MarlinAidaTTTrack::extrapolateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    // const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    // return this->extrapolateToDetElement( detElementID, site, ts, chi2, ndf, mode ) ;
    return success ;
  }
  
  
  int MarlinAidaTTTrack::extrapolateToDetElement( int detElementID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    // TKalTrackSite* site = 0;
    // int error_code = getSiteFromLCIOHit(trkhit, site);
    
    // if( error_code != success ) return error_code ;
    
    // return this->extrapolateToDetElement( detElementID, *site, ts, chi2, ndf, mode ) ;
    return success ;
  }
  
  
  
  int MarlinAidaTTTrack::propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){

    if( point[0] == 0.0 && point[1] == 0.0 && point[2] == 0.0 ) {

      return getTrackState(  ts, chi2, ndf ) ;


    }else{

      streamlog_out( WARNING )  << "MarlinAidaTTTrack::propagate not yet implemented for point otherthan IP " 
				<< std::endl ;





    }
    return success ;
  }
  
  int MarlinAidaTTTrack::propagate( const gear::Vector3D& point, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){
    
    return this->propagate( point, ts, chi2, ndf ) ;
  }
  
  
  
  int MarlinAidaTTTrack::propagateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    // const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    // return this->propagateToLayer( layerID, site, ts, chi2, ndf, detElementID, mode ) ;
    return success ;
    
  }
  
  
  int MarlinAidaTTTrack::propagateToLayer( int layerID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    // TKalTrackSite* site = 0;
    // int error_code = getSiteFromLCIOHit(trkhit, site);
    
    // if( error_code != success ) return error_code ;
    
    // return this->propagateToLayer( layerID, *site, ts, chi2, ndf, detElementID, mode ) ;
    return success ;
    
  }
  
  
  int MarlinAidaTTTrack::propagateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    // const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    // return this->propagateToDetElement( detElementID, site, ts, chi2, ndf, mode ) ;
    return success ;
    
  }
  
  
  int MarlinAidaTTTrack::propagateToDetElement( int detElementID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    // TKalTrackSite* site = 0;
    // int error_code = getSiteFromLCIOHit(trkhit, site);
    
    // if( error_code != success ) return error_code ;
    
    // return this->propagateToDetElement( detElementID, *site, ts, chi2, ndf, mode ) ;
    return success ;
    
  }
  
    
  int MarlinAidaTTTrack::intersectionWithDetElement( int detElementID, gear::Vector3D& point, int mode ) {  
    
    // const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    // const DDVMeasLayer* ml = 0;
    // return this->intersectionWithDetElement( detElementID, site, point, ml, mode ) ;
    return success ;
    
  }
  
  
  int MarlinAidaTTTrack::intersectionWithDetElement( int detElementID,  EVENT::TrackerHit* trkhit, gear::Vector3D& point, int mode ) {  
    
    // TKalTrackSite* site = 0;
    // int error_code = getSiteFromLCIOHit(trkhit, site);
    
    // if( error_code != success ) return error_code ;
    
    // const DDVMeasLayer* ml = 0;
    // return this->intersectionWithDetElement( detElementID, *site, point, ml, mode ) ;
    return success ;
    
  }
  
  
  int MarlinAidaTTTrack::intersectionWithLayer( int layerID, gear::Vector3D& point, int& detElementID, int mode ) {  
    
    // const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    // const DDVMeasLayer* ml = 0;
    // return this->intersectionWithLayer( layerID, site, point, detElementID, ml,  mode ) ;
    return success ;
    
  }
  
  
  int MarlinAidaTTTrack::intersectionWithLayer( int layerID,  EVENT::TrackerHit* trkhit, gear::Vector3D& point, int& detElementID, int mode ) {  
    
    // TKalTrackSite* site = 0;
    // int error_code = getSiteFromLCIOHit(trkhit, site);
    
    // if( error_code != success ) return error_code ;
    
    // const DDVMeasLayer* ml = 0;
    // return this->intersectionWithLayer( layerID, *site, point, detElementID, ml, mode ) ;
    return success ;
    
  }
  
  std::string MarlinAidaTTTrack::toString() {

    return std::string(" AidaTTTrack - to String ... " ) ;
  }

  
  
  
} // end of namespace MarlinTrk 
