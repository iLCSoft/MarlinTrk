#include "MarlinTrk/MarlinAidaTTTrack.h"

#include "MarlinTrk/MarlinAidaTT.h"
#include "MarlinTrk/IMarlinTrkSystem.h"



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

#include <lcio.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>

#include <UTIL/BitField64.h>
#include <UTIL/Operators.h>
#include <UTIL/ILDConf.h>

#include <sstream>

#include "streamlog/streamlog.h"

  
using namespace UTIL ;

namespace MarlinTrk {

  //---------------------------------------------------------------------------------------------------------------
  
  namespace{ 
    std::string cellIDString( int detElementID) {
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
    
    // ==== store hits in a map ===============
    std::map< int, EVENT::TrackerHit*> hitMap ;
    for(unsigned i=0 ; i < nHits ; ++i){
      hitMap[ _lcioHits[i]->getCellID0()  ] = _lcioHits[i] ;
    }
    
    //==== compute _all_ surface intersections ====================== 
    _intersections = &_fitTrajectory->getIntersectionsWithSurfaces( _aidaTT->_geom->getSurfaces() ) ;
    

    //========= loop over all intersections =========
    int pointLabel = 0 ; 
    for( std::vector<std::pair<double, const aidaTT::ISurface*> >::const_iterator it =  
	   _intersections->begin() ; it != _intersections->end() ; ++it ){
      
      const aidaTT::ISurface* surf = it->second ;

      _indexMap[ surf->id() ] = ++pointLabel ;  // label 0 is for the IP point 

      EVENT::TrackerHit* hit = hitMap[ surf->id() ] ;
      
      streamlog_out(DEBUG) << "MarlinAidaTTTrack::fit() - intersection " << pointLabel <<": << at s = " << it->first <<  " surface id : " << cellIDString( surf->id()  ) << std::endl ;

      if( hit != 0 ){ //-------- we have to add a measurement   
	
       	// get the hit position in dd4hep/aidaTT units
	double hitpos[3] ;
	for(unsigned int i = 0; i < 3; ++i) hitpos[i] = hit->getPosition()[i] * dd4hep::mm;
	
	//---- compute the precision from the hit errors
	double du,dv ;
	std::vector<double> precision ;
    
	
	TrackerHitPlane* planarhit = dynamic_cast<TrackerHitPlane*>( hit );
	if( planarhit != 0 ) {
	  
	  du = planarhit->getdU() * dd4hep::mm  ;
	  dv = planarhit->getdV() * dd4hep::mm  ;
	  
	} else { // we have a TPC hit which is not yet using the CylinderTrackerHit ...
	  
	  const FloatVec& cov = hit->getCovMatrix();

	  du = sqrt( cov[0] + cov[2] ) * dd4hep::mm  ;
	  dv = sqrt( cov[5]          ) * dd4hep::mm  ;
	}

	precision.push_back( 1./ (du*du) );

	if( ! surf->type().isMeasurement1D()  )
	  precision.push_back( 1./ (dv*dv) );
	
	_fitTrajectory->addMeasurement( hitpos, precision, *surf, hit , _aidaTT->_useQMS );

      } else  { // we just add a scatterer

	if (_aidaTT->_useQMS )
	  _fitTrajectory->addScatterer( *surf ) ;
      }
    }
      

    // do the fit:
    _fitTrajectory->prepareForFitting();
    
    int fit_ok = _fitTrajectory->fit();
    
    streamlog_out( DEBUG4 )  << "MarlinAidaTTTrack::fit() - fit worked  " << fit_ok  << std::endl ;
   
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
    
    std::map<int,int>::iterator it =  _indexMap.find( trkhit->getCellID0()  ) ;

    if( it == _indexMap.end() ) {
      
      streamlog_out( DEBUG2 )  << " MarlinAidaTTTrack::getTrackState(): " 
			       << " no surface intersection for given hit " 
			       << cellIDString( trkhit->getCellID0()) << std::endl  ;
      return error ; 
    }

    return getTrackState( aidaTT::Vector3D( trkhit->getPosition() ), it->second, ts, chi2, ndf ) ;
  }
  

  int MarlinAidaTTTrack::getTrackState( const aidaTT::Vector3D& refPoint, int label, 
					IMPL::TrackStateImpl& ts, double& chi2, int& ndf  ){

    const aidaTT::fitResults* result = _fitTrajectory->getFitResults( label );
    
    if( result == 0 ) {
      
      streamlog_out( DEBUG2 )  << " MarlinAidaTTTrack::getTrackState(): " 
			       << " no result at label " <<  label 
			       << " close to " << refPoint << std::endl  ;
      return error ;
    }
    
    // results are returned with origin as reference point
    aidaTT::trackParameters resTS = result->estimatedParameters() ;
    

    const double* pos = refPoint ;
    aidaTT::Vector3D hitPos( pos[0] * dd4hep::mm, pos[1] * dd4hep::mm , pos[2] * dd4hep::mm ) ; 
    //    hitPos = dd4hep::mm * hitPos ;
    
    streamlog_out( DEBUG2 )  << " MarlinAidaTTTrack::getTrackState(): " 
			     << " moving helix to " << hitPos  << std::endl ;
    
    bool calcCovMat = true ;
    moveHelixTo( resTS , hitPos , calcCovMat ) ;
    
    ts = *aidaTT::createLCIO( resTS );
        
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

    trkhit = _lcioHits[0] ; // is this correct ???
    return success;    
  }
  
  
  int MarlinAidaTTTrack::extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){  
    
    return propagate( point, ts, chi2, ndf ) ;
  }
  
  int MarlinAidaTTTrack::extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {
    
    return propagate( point, trkhit, ts, chi2, ndf ) ;
  }
  
  
  int MarlinAidaTTTrack::extrapolateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    return propagateToLayer( layerID, ts, chi2, ndf, detElementID, mode  ) ;
  }
  
  
  int MarlinAidaTTTrack::extrapolateToLayer( int layerID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    return propagateToLayer( layerID, trkhit, ts, chi2, ndf, detElementID, mode  ) ;
  }
  
  
  
  
  int MarlinAidaTTTrack::extrapolateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    return propagateToDetElement( detElementID, ts, chi2, ndf, mode  ) ;
  }
  
  
  int MarlinAidaTTTrack::extrapolateToDetElement( int detElementID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    return propagateToDetElement( detElementID, trkhit, ts, chi2, ndf, mode  ) ;
  }
  
  
  
  int MarlinAidaTTTrack::propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){

    if( point[0] == 0.0 && point[1] == 0.0 && point[2] == 0.0 ) {

      return getTrackState(  ts, chi2, ndf ) ;


    }else{

      aidaTT::Vector3D thePoint( point[0]*dd4hep::mm,point[1]*dd4hep::mm, point[2]*dd4hep::mm ) ;

      //========= loop over all intersections and find the one closest to given point
      double minDist2 = 1.e99 ;
      unsigned index = 0, count = 0 ;
      for( std::vector<std::pair<double, const aidaTT::ISurface*> >::const_iterator it =  
	     _intersections->begin() ; it != _intersections->end() ; ++it ){
	
	++count ;

	// get the (cached) intersection point

	double s ; aidaTT::Vector2D uv ; aidaTT::Vector3D position ;
	_fitTrajectory->_calculateIntersectionWithSurface( it->second, s , &uv, &position );
	
	aidaTT::Vector3D dv = position - thePoint ;
	double dist2 = dv.r2() ;
	
	if( dist2 < minDist2 ){
	  minDist2 = dist2 ;
	  index = count ;
	}
      }

      streamlog_out( DEBUG2 )  << "MarlinAidaTTTrack::propagate(): found closest intersection "
			       << " to point " << point << " at surface " 
			       << *(*_intersections)[index].second << std::endl ;



      return getTrackState( aidaTT::Vector3D( point[0], point[1], point[2] ), 
			    index, ts, chi2, ndf ) ;

    }
  }
  
  int MarlinAidaTTTrack::propagate( const gear::Vector3D& point, EVENT::TrackerHit*, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){
    
    return this->propagate( point, ts, chi2, ndf ) ;
  }
  
  
  
  int MarlinAidaTTTrack::propagateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int ) { 
    
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
    encoder.reset() ;  // reset to 0
    

    // compute a mask for the layerid
    int mask=0 ;
    mask |= encoder[lcio::ILDCellID0::subdet].mask() ;
    mask |= encoder[lcio::ILDCellID0::side  ].mask() ;
    mask |= encoder[lcio::ILDCellID0::layer ].mask() ;

    // loop over intersections to find the (first) intersection w/ the given layerid
    
    double s ; aidaTT::Vector2D uv ; aidaTT::Vector3D position ;

    //========= loop over all intersections and find on that matches layerID
    int theID = -1 ;
    unsigned index = 0, count = 0 ;
    for( std::vector<std::pair<double, const aidaTT::ISurface*> >::const_iterator it =  
	   _intersections->begin() ; it != _intersections->end() ; ++it ){
      ++count ;

      int id = it->second->id() ;

      if( layerID == ( id & mask ) ){

	index = count ;
	theID = id ;
	
	_fitTrajectory->_calculateIntersectionWithSurface( it->second, s , &uv, &position );
	
	break ;
      }
    }

    if( theID == -1 ){
      streamlog_out( ERROR )  << "MarlinAidaTTTrack::propagate(): no intersection "
			      << " found for layerID "  << cellIDString( layerID ) 
			      << std::endl ;
      return error ;
    }


    detElementID = theID ;

    // need to convert intersection position back to mm ...
    aidaTT::Vector3D point( position[0]/dd4hep::mm, position[1]/dd4hep::mm,position[2]/dd4hep::mm ) ;

    int res = getTrackState( point, index, ts, chi2, ndf ) ;

    return res ;
  }
  
  

  int MarlinAidaTTTrack::propagateToLayer( int layerID, EVENT::TrackerHit*, IMPL::TrackStateImpl& ts, 
					   double& chi2, int& ndf, int& detElementID, int  ) { 
    
    return propagateToLayer( layerID, ts, chi2, ndf, detElementID ) ;
  }
  
  
  int MarlinAidaTTTrack::propagateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int ) { 
    
    std::map<int,int>::iterator it =  _indexMap.find( detElementID ) ;
    
    if( it == _indexMap.end() ) {
      
      streamlog_out( DEBUG2 )  << " MarlinAidaTTTrack::propagateToDetElement(): " 
			       << " no surface intersection for given DetElementID " 
			       << cellIDString( detElementID ) << std::endl  ;
      return error ; 
    }
    
    
    const aidaTT::ISurface* surf = (*_intersections)[ it->second - 1  ].second ;

    // sanity check:
    if( surf->id() != detElementID ){
      
      std::stringstream s ; s << "MarlinAidaTTTrack::propagateToDetElement() - inconsistent ids: detElementID = " 
			      << cellIDString( detElementID ) << " and surf.id() " <<   surf->id() << std::endl ;

      throw MarlinTrk::Exception( s.str() );   
    }

    double s ; aidaTT::Vector2D uv ; aidaTT::Vector3D position ;

    _fitTrajectory->_calculateIntersectionWithSurface( surf, s , &uv, &position );
    
    // need to convert intersection position back to mm ...
    aidaTT::Vector3D point( position[0]/dd4hep::mm, position[1]/dd4hep::mm,position[2]/dd4hep::mm ) ;
    
    return getTrackState( point , it->second, ts, chi2, ndf ) ;
  }
  
  
  int MarlinAidaTTTrack::propagateToDetElement( int detElementID, EVENT::TrackerHit*, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int ) { 
    
    return propagateToDetElement(detElementID, ts, chi2, ndf ) ;
  }
  
    
  int MarlinAidaTTTrack::intersectionWithDetElement( int detElementID, gear::Vector3D& point, int) {  
    
    SurfMap::iterator it = _aidaTT->_surfMap.find( detElementID ) ;

    if( it == _aidaTT->_surfMap.end() ){
      streamlog_out( DEBUG2 )  << " MarlinAidaTTTrack::intersectionWithDetElement() - no surface found for " << cellIDString( detElementID ) << std::endl ;
      return error ;
    }

    const aidaTT::ISurface* surf = it->second ;
    
    double s ; aidaTT::Vector2D uv ; aidaTT::Vector3D position ;

    bool intersects = _fitTrajectory->_calculateIntersectionWithSurface( surf, s , &uv, &position );

    // need to convert intersection position back to mm ...
    if( intersects) 
      point =  aidaTT::Vector3D( position[0]/dd4hep::mm, position[1]/dd4hep::mm,position[2]/dd4hep::mm ) ;

    return (intersects ? success : error ) ;
  }
  
  
  int MarlinAidaTTTrack::intersectionWithDetElement( int detElementID,  EVENT::TrackerHit*, gear::Vector3D& point, int mode ) {  
    
    return intersectionWithDetElement( detElementID, point, mode ) ;
  }
  
  
  int MarlinAidaTTTrack::intersectionWithLayer( int layerID, gear::Vector3D& point, int& detElementID, int mode ) {  
    
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
    encoder.reset() ;  // reset to 0
    
    // compute a mask for the layerid
    int mask=0 ;
    mask |= encoder[lcio::ILDCellID0::subdet].mask() ;
    mask |= encoder[lcio::ILDCellID0::side  ].mask() ;
    mask |= encoder[lcio::ILDCellID0::layer ].mask() ;

    int theID = -1 ;
    for( std::vector<std::pair<double, const aidaTT::ISurface*> >::const_iterator it =  
	   _intersections->begin() ; it != _intersections->end() ; ++it ){
      
      int id = it->second->id() ;
      
      if( layerID == ( id & mask ) ){
	theID = id ;
	break ;
      }
    }
    
    if( theID == -1 ){
      streamlog_out( ERROR )  << "MarlinAidaTTTrack::intersectionWithLayer(): no intersection "
			      << " found for layerID "  << cellIDString( layerID ) 
			      << std::endl ;
      return error ;
    }

    detElementID = theID ;

    return intersectionWithDetElement( detElementID, point, mode  ) ;
  }
  
  
  int MarlinAidaTTTrack::intersectionWithLayer( int layerID,  EVENT::TrackerHit*, gear::Vector3D& point, int& detElementID, int ) {  

    return intersectionWithLayer( layerID, point, detElementID ) ;
  }
  
  std::string MarlinAidaTTTrack::toString() {

    return std::string(" AidaTTTrack - to String ... " ) ;
  }

  
  
  
} // end of namespace MarlinTrk 
