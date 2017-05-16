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
#include "aidaTT/materialUtils.hh"

#include "DD4hep/DD4hepUnits.h"

#include <lcio.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>

#include <UTIL/BitField64.h>
#include <UTIL/Operators.h>
#include <UTIL/ILDConf.h>
#include "UTIL/LCTrackerConf.h"

#include <sstream>

#include "streamlog/streamlog.h"

  
using namespace UTIL ;

namespace MarlinTrk {

  //---------------------------------------------------------------------------------------------------------------
  
  namespace{ 
    std::string cellIDString( int detElementID) {
      lcio::BitField64 bf(  UTIL::LCTrackerCellID::encoding_string() ) ;
      bf.setValue( detElementID ) ;
      return bf.valueString() ;
    }
  }
  //---------------------------------------------------------------------------------------------------------------
  
  
  MarlinAidaTTTrack::MarlinAidaTTTrack( MarlinAidaTT* mAidaTT) 
    : _aidaTT( mAidaTT ) , _initialised( false ) , _mass( aidaTT::pionMass )  {
    
  }
  
  
  MarlinAidaTTTrack::~MarlinAidaTTTrack(){
    delete _fitTrajectory ;
  }
  

  void MarlinAidaTTTrack::setMass(double mass) { _mass =  mass ;  } 
  
  double MarlinAidaTTTrack::getMass() { return _mass ; }


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

    
    aidaTT::trackParameters tp ;
    calculateStartHelix( x1, x2,  x3 , tp , backwards ) ;
             
    streamlog_out( DEBUG3 )  << "  start helix from three points : " <<  tp << std::endl ;


#if 1 // create prefit from this helix
    _initialTrackParams = createPreFit( tp ) ;
#else
    _initialTrackParams = tp ;
#endif

    return myInit() ;
  }
  

  int MarlinAidaTTTrack::initialise(  const EVENT::TrackState& ts, double, bool dir) {
    
    if ( _initialised ) {
      throw MarlinTrk::Exception("Track fit already initialised");   
    }

#if 0 // debug code - initialize with hits and a prefit ....
    
    return initialise( dir ) ;

#else
    _initialTrackParams = aidaTT::readLCIO( &ts ) ;

    return myInit() ;
#endif

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

    _fitTrajectory = new aidaTT::trajectory( _initialTrackParams, _aidaTT->_fitter, //_aidaTT->_bfield, 
					     _aidaTT->_propagation, _aidaTT->_geom );

    _fitTrajectory->setMass( _mass ) ;

    // add the Interaction Point as the first element of the trajectory
    // int ID = 1;
    // _fitTrajectory->addElement( aidaTT::Vector3D(), &ID);
    //done internally in trajectory ...

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

      EVENT::TrackerHit* hit = hitMap[ surf->id() ] ;
      
      streamlog_out(DEBUG7) << "MarlinAidaTTTrack::fit() - intersection - current pointLabel : " << pointLabel  
			   << ":  at s = " << it->first <<  " surface id : " 
			   << cellIDString( surf->id()  ) << std::endl 
			   << *surf << std::endl ;

      if( hit != 0 ){ //-------- we have to add a measurement   
	
	double hitpos[3] ;
	std::vector<double> precision ;
	getHitInfo( hit, hitpos, precision , surf) ;

	_fitTrajectory->addMeasurement( hitpos, precision, *surf, hit , _aidaTT->_useQMS );
	_indexMap[ surf->id() ] = ++pointLabel ;  // label 0 is for the IP point 

	streamlog_out(DEBUG) << "MarlinAidaTTTrack::fit()  addMeasurement called for pointLabel : " << pointLabel << std::endl ;

      } else  { // we just add a scatterer

	if (_aidaTT->_useQMS ){
	  
	  // ignore virtual surface with no material (e.g. inside the beam pipe )

	  if( ! ( surf->innerMaterial().density() < 1e-6  && 
		  surf->outerMaterial().density() < 1e-6 )  ) {

	    _fitTrajectory->addScatterer( *surf ) ;
	    _indexMap[ surf->id() ] = ++pointLabel ;  // label 0 is for the IP point 

	    streamlog_out(DEBUG) << "MarlinAidaTTTrack::fit()  addScatterer called for pointLabel : " << pointLabel << std::endl ;
	  }
	}
      }
    }
      

    // do the fit:
    _fitTrajectory->prepareForFitting();
    
    int fit_ok = _fitTrajectory->fit();
    
    streamlog_out( DEBUG4 )  << "MarlinAidaTTTrack::fit() - fit worked  " << fit_ok  
			     << std::endl ;
   

#if 0 // ---------------- try to run a refit --- does not work really ....
    // refit one more time w/ new start parameters
    const aidaTT::fitResults* result = _fitTrajectory->getFitResults();
    const aidaTT::trackParameters& tp =  result->estimatedParameters()  ;

    streamlog_out( MESSAGE ) << " --- first fit result : " << tp << std::endl ;
 
    _fitTrajectory->setInitialTrackParameters( tp );
    _fitTrajectory->prepareForFitting();
    fit_ok = _fitTrajectory->fit();
    streamlog_out( DEBUG4 )  << "MarlinAidaTTTrack::fit() - refit worked  " << fit_ok  << std::endl ;

    result = _fitTrajectory->getFitResults();
    const aidaTT::trackParameters& tp1 =  result->estimatedParameters()  ;

    streamlog_out( MESSAGE ) << " --- second fit result : " << tp1 << std::endl ;
#endif

    return ( fit_ok ? success : error  )  ;
  }
  
  
  aidaTT::trackParameters MarlinAidaTTTrack::createPreFit(aidaTT::trackParameters& tp ){
    
    // create a prefit from the hits w/o QMS and dEdx
    
    moveHelixTo( tp, aidaTT::Vector3D()  ) ; // move to origin

    aidaTT::trajectory traj(  tp , _aidaTT->_fitter, //_aidaTT->_bfield, 
			      _aidaTT->_propagation, _aidaTT->_geom ) ;
    
    traj.setMass( _mass ) ;

    // Add the Interaction Point as the first element of the trajectory
    // int ID = 1;
    // aidaTT::Vector3D IntPoint(0,0,0);
    // traj.addElement(IntPoint, &ID);
    // done internally in trajectory    

    // try to use up to 25 or so hits....
    unsigned nHits=_lcioHits.size() ;

    int step = ( nHits <= 25  ? 1 : int( 1.*nHits/25. )  ) ; 

    streamlog_out( DEBUG1 ) << " MarlinAidaTTTrack::createPreFit() :  will use every " 
			    <<  step << "-th hit for prefit ! " << std::endl ; 

    for(unsigned i=0 ; i < nHits ; i+= step ){
      
      EVENT::TrackerHit* hit = _lcioHits[i] ;
      
      long hitid = hit->getCellID0() ;

      SurfMap::iterator it = _aidaTT->_surfMap.find( hitid ) ;

      if( it == _aidaTT->_surfMap.end() ){

	streamlog_out( DEBUG1 ) << " MarlinAidaTTTrack::createPreFit() : no surface found for id : " 
				<< cellIDString( hitid ) << std::endl ;
	continue;
      }

      const aidaTT::ISurface* surf = it->second ;
      
      double hitpos[3] ;
      std::vector<double> precision ;
      getHitInfo( hit, hitpos, precision , surf) ;
      
      streamlog_out( DEBUG1 ) << " MarlinAidaTTTrack::createPreFit() : adding hit for "
			      <<  cellIDString( hitid ) << " at " 
			      << aidaTT::Vector3D( hit->getPosition() ) << std::endl ;
	
      

      traj.addMeasurement( hitpos, precision, *surf, hit );
    }

    traj.prepareForFitting();
	      
    int success = traj.fit();
	      

    if( success ) { 

      const aidaTT::fitResults* result = traj.getFitResults();
      
      const aidaTT::trackParameters& newTP = result->estimatedParameters() ;

      streamlog_out( DEBUG4 ) << " MarlinAidaTTTrack::createPreFit() : prefit tp: " 
			      << newTP << std::endl ;


      return newTP ;

    } else {

      streamlog_out( WARNING ) << " MarlinAidaTTTrack::createPreFit() : prefit failed for tp: " 
			       << tp << std::endl ;
      
      return tp ;
    }

  }


  void MarlinAidaTTTrack::getHitInfo( const EVENT::TrackerHit* hit, double* hitpos, 
				      std::vector<double>& precision, const aidaTT::ISurface* surf){
    
    // get the hit position in dd4hep/aidaTT units
    for(unsigned int i = 0; i < 3; ++i) hitpos[i] = hit->getPosition()[i] * dd4hep::mm;
    
    //---- compute the precision from the hit errors
    double du,dv ;
    
    const TrackerHitPlane* planarhit = dynamic_cast<const TrackerHitPlane*>( hit );
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
    else
      precision.push_back( 0. );
    
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
  
  
  int MarlinAidaTTTrack::extrapolate( const Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){  
    
    return propagate( point, ts, chi2, ndf ) ;
  }
  
  int MarlinAidaTTTrack::extrapolate( const Vector3D& point, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {
    
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
  
  
  
  int MarlinAidaTTTrack::propagate( const Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){

    if( point[0] == 0.0 && point[1] == 0.0 && point[2] == 0.0 ) {

      return getTrackState(  ts, chi2, ndf ) ;


    }else{

      aidaTT::Vector3D thePoint( point[0]*dd4hep::mm,point[1]*dd4hep::mm, point[2]*dd4hep::mm ) ;

      //========= loop over all intersections and find the one closest to given point
      // double minDist2 = 1.e99 ;
      // unsigned index = 0, count = 0 ;
      // for( std::vector<std::pair<double, const aidaTT::ISurface*> >::const_iterator it =  
      // 	     _intersections->begin() ; it != _intersections->end() ; ++it ){
	
      // 	// get the (cached) intersection point

      // 	double s ; aidaTT::Vector2D uv ; aidaTT::Vector3D position ;
      // 	_fitTrajectory->_calculateIntersectionWithSurface( it->second, s , &uv, &position );
	
      // 	aidaTT::Vector3D dv = position - thePoint ;
      // 	double dist2 = dv.r2() ;
	
      // 	if( dist2 < minDist2 ){
      // 	  minDist2 = dist2 ;
      // 	  index = count ;
      // 	}
      // 	++count ;
      // }
      // streamlog_out( DEBUG2 )  << "MarlinAidaTTTrack::propagate(): found closest intersection "
      // 			       << " to point " << point << " at surface " 
      // 			       << *(*_intersections)[index].second << std::endl ;
      
      //      const aidaTT::ISurface* s = _intersections->at(index).second ;
      //int label = _indexMap[ s->id() ] ;

      const std::vector<aidaTT::trajectoryElement*>& elemVec = 	_fitTrajectory->trajectoryElements() ;
      double minDist2 = 1.e99 ;
      unsigned label = 0 ;
      for( unsigned i=0,N=elemVec.size() ; i < N ; ++i ){
	const aidaTT::Vector3D& rp  = elemVec[i]->getTrackParameters()->referencePoint() ; 
       	aidaTT::Vector3D dv = rp - thePoint ;
       	double dist2 = dv.r2() ;
       	if( dist2 < minDist2 ){
       	  minDist2 = dist2 ;
       	  label = i ;
       	}
      }

      return getTrackState( aidaTT::Vector3D( point[0], point[1], point[2] ), 
			    label, ts, chi2, ndf ) ;

    }
  }
  
  int MarlinAidaTTTrack::propagate( const Vector3D& point, EVENT::TrackerHit*, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){
    
    return this->propagate( point, ts, chi2, ndf ) ;
  }
  
  
  
  int MarlinAidaTTTrack::propagateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int ) { 
    
    UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ; 
    encoder.reset() ;  // reset to 0
    

    // compute a mask for the layerid
    int mask=0 ;
    mask |= encoder[lcio::LCTrackerCellID::subdet()].mask() ;
    mask |= encoder[lcio::LCTrackerCellID::side()  ].mask() ;
    mask |= encoder[lcio::LCTrackerCellID::layer() ].mask() ;

    // loop over intersections to find the (first) intersection w/ the given layerid
    //    double s ; aidaTT::Vector2D uv ; aidaTT::Vector3D position ;
    // //========= loop over all intersections and find on that matches layerID
    // int theID = -1 ;
    // unsigned index = 0, count = 0 ;
    // for( std::vector<std::pair<double, const aidaTT::ISurface*> >::const_iterator it =  
    // 	   _intersections->begin() ; it != _intersections->end() ; ++it ){
    //   ++count ;

    //   int id = it->second->id() ;

    //   if( layerID == ( id & mask ) ){

    // 	index = count ;
    // 	theID = id ;
	
    // 	_fitTrajectory->_calculateIntersectionWithSurface( it->second, s , &uv, &position );
	
    // 	break ;
    //   }
    // }

    streamlog_out( DEBUG )  << "MarlinAidaTTTrack::propagate():  looking for cellID " << cellIDString( layerID )  << std::endl ;

    const std::vector<aidaTT::trajectoryElement*>& elemVec = _fitTrajectory->trajectoryElements() ;
    unsigned label = 0 ;
    int theID = -1 ;
    unsigned index = 0 ;
    // first element has not surface ...
    for( unsigned i=1,N=elemVec.size() ; i < N ; ++i ){
      
      int id =  elemVec[i]->surface().id() ;
      
      
      streamlog_out( DEBUG )  << "MarlinAidaTTTrack::propagate():  comparing to cellID " << cellIDString( id ) << std::endl ; 

      if( layerID == ( id & mask ) ){
	theID = id ;
	label = i ;
	break ;
      }
    }

    if( theID == -1 ){
      streamlog_out( ERROR )  << "MarlinAidaTTTrack::propagate(): no intersection "
			      << " found for layerID "  << cellIDString( layerID ) 
			      << std::endl ;
      return no_intersection ;
    }
    
    detElementID = theID ;

    const aidaTT::Vector3D& position  = elemVec[label]->getTrackParameters()->referencePoint() ; 

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
    
    
    // const aidaTT::ISurface* surf = (*_intersections)[ it->second - 1  ].second ;
    // // sanity check:
    // if( surf->id() != detElementID ){
    //   std::stringstream s ; s << "MarlinAidaTTTrack::propagateToDetElement() - inconsistent ids: detElementID = " 
    // 			      << cellIDString( detElementID ) << " and surf.id() " <<   surf->id() << std::endl ;
    //   throw MarlinTrk::Exception( s.str() );   
    // }
    // double s ; aidaTT::Vector2D uv ; aidaTT::Vector3D position ;
    // _fitTrajectory->_calculateIntersectionWithSurface( surf, s , &uv, &position );
    
    int label = it->second ;
    const std::vector<aidaTT::trajectoryElement*>& elemVec = 	_fitTrajectory->trajectoryElements() ;
    const aidaTT::Vector3D& position  = elemVec[label]->getTrackParameters()->referencePoint() ; 

    // need to convert intersection position back to mm ...
    aidaTT::Vector3D point( position[0]/dd4hep::mm, position[1]/dd4hep::mm,position[2]/dd4hep::mm ) ;
    
    return getTrackState( point , it->second, ts, chi2, ndf ) ;
  }
  
  
  int MarlinAidaTTTrack::propagateToDetElement( int detElementID, EVENT::TrackerHit*, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int ) { 
    
    return propagateToDetElement(detElementID, ts, chi2, ndf ) ;
  }
  
    
  int MarlinAidaTTTrack::intersectionWithDetElement( int detElementID, Vector3D& point, int) {  
    
 //    SurfMap::iterator it = _aidaTT->_surfMap.find( detElementID ) ;
 //    if( it == _aidaTT->_surfMap.end() ){
 //      streamlog_out( DEBUG2 )  << " MarlinAidaTTTrack::intersectionWithDetElement() - no surface found for " << cellIDString( detElementID ) << std::endl ;
 //      return error ;
 //    }
 //    const aidaTT::ISurface* surf = it->second ;
 //    double s ; aidaTT::Vector2D uv ; aidaTT::Vector3D position ;
 //    bool intersects = _fitTrajectory->_calculateIntersectionWithSurface( surf, s , &uv, &position );
 //    // need to convert intersection position back to mm ...
 // if( intersects) 
    std::map<int,int>::iterator it =  _indexMap.find( detElementID ) ;
    
    if( it == _indexMap.end() ) {
      
      streamlog_out( DEBUG2 )  << " MarlinAidaTTTrack::intersectionWithDetElement(): " 
			       << " no surface intersection for given DetElementID " 
			       << cellIDString( detElementID ) << std::endl  ;
      return error ; 
    }
    int label = it->second ;
    const std::vector<aidaTT::trajectoryElement*>& elemVec = 	_fitTrajectory->trajectoryElements() ;
    const aidaTT::Vector3D& position  = elemVec[label]->getTrackParameters()->referencePoint() ; 


      point =  aidaTT::Vector3D( position[0]/dd4hep::mm, position[1]/dd4hep::mm,position[2]/dd4hep::mm ) ;

      return success ;
      //    return (intersects ? success : error ) ;
  }
  
  
  int MarlinAidaTTTrack::intersectionWithDetElement( int detElementID,  EVENT::TrackerHit*, Vector3D& point, int mode ) {  
    
    return intersectionWithDetElement( detElementID, point, mode ) ;
  }
  
  
  int MarlinAidaTTTrack::intersectionWithLayer( int layerID, Vector3D& point, int& detElementID, int mode ) {  
    
    UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ; 
    encoder.reset() ;  // reset to 0
    
    // compute a mask for the layerid
    int mask=0 ;
    mask |= encoder[lcio::LCTrackerCellID::subdet()].mask() ;
    mask |= encoder[lcio::LCTrackerCellID::side()  ].mask() ;
    mask |= encoder[lcio::LCTrackerCellID::layer() ].mask() ;

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
  
  
  int MarlinAidaTTTrack::intersectionWithLayer( int layerID,  EVENT::TrackerHit*, Vector3D& point, int& detElementID, int ) {  

    return intersectionWithLayer( layerID, point, detElementID ) ;
  }
  
  std::string MarlinAidaTTTrack::toString() {

    return std::string(" AidaTTTrack - to String ... " ) ;
  }

  
  
  
} // end of namespace MarlinTrk 
