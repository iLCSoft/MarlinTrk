
#include "MarlinTrk/MarlinKalTestTrack.h"

#include "MarlinTrk/MarlinKalTest.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

#include <kaltest/TKalDetCradle.h>
#include <kaltest/TKalTrack.h>
#include <kaltest/TKalTrackState.h>
#include "kaltest/TKalTrackSite.h"
#include "TKalFilterCond.h"

#include <lcio.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>

#include <UTIL/BitField64.h>

#include <ILDCellIDEncoding.h>

#include "kaldet/ILDCylinderMeasLayer.h"
#include "kaldet/ILDCylinderHit.h"

#include "kaldet/ILDPlanarMeasLayer.h"
#include "kaldet/ILDPlanarHit.h"

#include "gear/GEAR.h"
#include "gear/BField.h"

#include "streamlog/streamlog.h"


/** Helper class for defining a filter condition based on the delta chi2 in the AddAndFilter step.
 */
class KalTrackFilter : public TKalFilterCond{

public:
  
  /** C'tor - takes as optional argument the maximum allowed delta chi2 for adding the hit (in IsAccepted() )
   */
  KalTrackFilter(double maxDeltaChi2 = DBL_MAX) : _maxDeltaChi2( maxDeltaChi2 ) {
  } 
  virtual ~KalTrackFilter() {} 
  
  virtual Bool_t IsAccepted(const TKalTrackSite &site) {
    
    double deltaChi2 = site.GetDeltaChi2();
    
    streamlog_out( DEBUG3 ) << " KalTrackFilter::IsAccepted called  !  deltaChi2 = "  <<  deltaChi2  << std::endl;

    return ( deltaChi2 < _maxDeltaChi2 )   ; 
  }

protected:

  double _maxDeltaChi2 ;

} ;
//---------------------------------------------------------------------------------------------------------------



MarlinKalTestTrack::MarlinKalTestTrack( MarlinKalTest* ktest) 
  : _ktest(ktest)
{

  _kaltrack = new TKalTrack() ;
  _kaltrack->SetOwner() ;

  _kalhits = new TObjArray() ;

  _initialised = false ;
  _smoothed = false ;

}


MarlinKalTestTrack::~MarlinKalTestTrack(){
  delete _kaltrack ;
  delete _kalhits ;
}



int MarlinKalTestTrack::addHit( EVENT::TrackerHit * trkhit) 
{

  return this->addHit( trkhit, _ktest->findMeasLayer( trkhit )) ;
 
} 

int MarlinKalTestTrack::addHit( EVENT::TrackerHit * trkhit, const ILDVMeasLayer* ml) 
{
  if( trkhit && ml ) {
    return this->addHit( trkhit, ml->ConvertLCIOTrkHit(trkhit), ml) ;
  }
  else {
    return 1 ;
  }

}

int MarlinKalTestTrack::addHit( EVENT::TrackerHit* trkhit, ILDVTrackHit* kalhit, const ILDVMeasLayer* ml) 
{

  if( kalhit && ml ) {
    _kalhits->Add(kalhit ) ;  // Add hit and set surface found 
    _lcio_hits_to_kaltest_hits[trkhit] = kalhit ; // add hit to map relating lcio and kaltest hits
    _kaltest_hits_to_lcio_hits[kalhit] = trkhit ; // add hit to map relating kaltest and lcio hits
  }
  else{
    return 1 ;
  }

  streamlog_out(DEBUG3) << "MarlinKalTestTrack::MarlinKalTestTrack hit added " 
			<< "number of hits for track = " << _kalhits->GetEntries() 
			<< std::endl ;

  return 0 ;

}


int MarlinKalTestTrack::initialise( bool direction ) {; 

  //SJA:FIXME: check here if the track is already initialised, and for now don't allow it to be re-initialised
  //           if the track is going to be re-initialised then we would need to do it directly on the first site
  if ( _initialised ) {
    
    throw MarlinTrk::Exception("Track fit already initialised");   
    
  }

  if (_kalhits->GetEntries() < 3) {
    
    streamlog_out( ERROR) << "<<<<<< MarlinKalTestTrack::initialise: Shortage of Hits! nhits = "  
			  << _kalhits->GetEntries() << " >>>>>>>" << std::endl;
    return 1;

  }

  const bool gkDir = direction ; 

  // establish the hit order
  Int_t i1, i2, i3; // (i1,i2,i3) = (1st,mid,last) hit to filter
  if (gkDir == kIterBackward) {
    i3 = 0 ; // fg: first index is 0 and not 1 
    i1 = _kalhits->GetEntries() - 1;
    i2 = i1 / 2;
  } else {
    i1 = 0 ; 
    i3 = _kalhits->GetEntries() - 1;
    i2 = i3 / 2;
  }

  TVTrackHit *startingHit = dynamic_cast<TVTrackHit *>(_kalhits->At(i1));
  
  // ---------------------------
  //  Create an initial start site for the track using the first hit
  // ---------------------------
  // set up a dummy hit needed to create initial site  

  TVTrackHit* pDummyHit = NULL;

  if ( (pDummyHit = dynamic_cast<ILDCylinderHit *>( startingHit )) ) {
    pDummyHit = (new ILDCylinderHit(*static_cast<ILDCylinderHit*>( startingHit )));
  }
  else if ( (pDummyHit = dynamic_cast<ILDPlanarHit *>( startingHit )) ) {
    pDummyHit = (new ILDPlanarHit(*static_cast<ILDPlanarHit*>( startingHit )));
  }
  else {
    streamlog_out( ERROR) << "<<<<<<<<< MarlinKalTestTrack::initialise: dynamic_cast failed for hit type >>>>>>>" << std::endl;
    return false;
  }

  TVTrackHit& dummyHit = *pDummyHit;

  //SJA:FIXME: this constants should go in a header file
  // give the dummy hit huge errors so that it does not contribute to the fit
  dummyHit(0,1) = 1.e6;   // give a huge error to d
  dummyHit(1,1) = 1.e6;   // give a huge error to z   
  
  // use dummy hit to create initial site
  TKalTrackSite& initialSite = *new TKalTrackSite(dummyHit);
  
  initialSite.SetHitOwner();// site owns hit
  initialSite.SetOwner();   // site owns states

 // ---------------------------
  //  Create initial helix
  // ---------------------------

  TVTrackHit &h1 = *dynamic_cast<TVTrackHit *>(_kalhits->At(i1)); // first hit
  TVTrackHit &h2 = *dynamic_cast<TVTrackHit *>(_kalhits->At(i2)); // middle hit
  TVTrackHit &h3 = *dynamic_cast<TVTrackHit *>(_kalhits->At(i3)); // last hit
  TVector3    x1 = h1.GetMeasLayer().HitToXv(h1);
  TVector3    x2 = h2.GetMeasLayer().HitToXv(h2);
  TVector3    x3 = h3.GetMeasLayer().HitToXv(h3);

  // create helix using 3 global space points 
  THelicalTrack helstart(x1, x2, x3, h1.GetBfield(), gkDir); // initial helix 

  // ---------------------------
  //  Set up initial track state ... could try to use lcio track parameters ...
  // ---------------------------

  static TKalMatrix initialState(kSdim,1) ;
  initialState(0,0) = 0.0 ;                       // dr
  initialState(1,0) = helstart.GetPhi0() ;        // phi0
  initialState(2,0) = helstart.GetKappa() ;       // kappa
  initialState(3,0) = 0.0 ;                       // dz
  initialState(4,0) = helstart.GetTanLambda() ;   // tan(lambda)
  if (kSdim == 6) initialState(5,0) = 0.;         // t0


  // ---------------------------
  //  Set up initial Covariance Matrix with very large errors 
  // ---------------------------
  
  static TKalMatrix Cov(kSdim,kSdim);
  for (Int_t i=0; i<kSdim; i++) {
    Cov(i,i) = 1.e6;   // initialise diagonal elements of dummy error matrix
  }


  // Add initial states to the site 
  initialSite.Add(new TKalTrackState(initialState,Cov,initialSite,TVKalSite::kPredicted));
  initialSite.Add(new TKalTrackState(initialState,Cov,initialSite,TVKalSite::kFiltered));

  // add the initial site to the track: that is, give the track initial parameters and covariance 
  // matrix at the starting measurement layer
  _kaltrack->Add(&initialSite);

  _initialised = true ;

  return 0 ;

}

int MarlinKalTestTrack::initialise( const IMPL::TrackStateImpl& ts, double bfield_z, bool initalise_at_end ) {

  //SJA:FIXME: check here if the track is already initialised, and for now don't allow it to be re-initialised
  //           if the track is going to be re-initialised then we would need to do it directly on the first site
  if ( _initialised ) {
    
    throw MarlinTrk::Exception("Track fit already initialised");   
    
  }
  
  // for GeV, Tesla, R in mm  
  double alpha = 2.99792458E-4 ;
  double kappa = ts.getOmega() * bfield_z * alpha ;

  THelicalTrack helix( ts.getD0(),
		       ts.getPhi(),
		       kappa,
		       ts.getZ0(),
		       ts.getTanLambda(),
		       ts.getReferencePoint()[0],
		       ts.getReferencePoint()[1],
		       ts.getReferencePoint()[2],
		       bfield_z );

  TMatrixD cov(5,5) ;	
  EVENT::FloatVec covLCIO( 15 )  ; 

  cov( 0 , 0 )  =   covLCIO[ 0] ; //   d0,   d0
      
  cov( 1 , 0 )  = - covLCIO[ 1] ; //   phi0, d0
  cov( 1 , 1 )  =   covLCIO[ 2] ; //   phi0, phi 
      
  cov( 2 , 0 ) = - covLCIO[ 3] / alpha ;           //   omega, d0
  cov( 2 , 1 ) =   covLCIO[ 4] / alpha ;           //   omega, phi
  cov( 2 , 2 ) =   covLCIO[ 5] / (alpha * alpha) ; //   omega, omega
      
  cov( 3 , 0 ) = - covLCIO[ 6] ;         //   z0  , d0
  cov( 3 , 1 ) =   covLCIO[ 7] ;         //   z0  , phi
  cov( 3 , 2 ) =   covLCIO[ 8] / alpha ; //   z0  , omega
  cov( 3 , 3 ) =   covLCIO[ 9] ;         //   z0  , z0
      
  cov( 4 , 0 ) = - covLCIO[10] ;         //   tanl, d0 
  cov( 4 , 1 ) =   covLCIO[11] ;         //   tanl, phi
  cov( 4 , 2 ) =   covLCIO[12] / alpha ; //   tanl, omega    
  cov( 4 , 3 ) =   covLCIO[13] ;         //   tanl, z0
  cov( 4 , 4 ) =   covLCIO[14] ;         //   tanl, tanl

  // move the helix to either the position of the last hit or the first depending on initalise_at_end

  // default case initalise_at_end
  int index = _kalhits->GetEntries() - 1 ;
  // or initialise at start 
  if( ! initalise_at_end ){
    index = 0 ;
  }

  TVTrackHit* kalhit = dynamic_cast<TVTrackHit *>(_kalhits->At(index)); 

  TVector3 initial_pivot = kalhit->GetMeasLayer().HitToXv(*kalhit);
  
  double dphi;
  helix.MoveTo( initial_pivot, dphi, NULL, &cov );

  // ---------------------------
  //  Create an initial start site for the track using the  hit
  // ---------------------------
  // set up a dummy hit needed to create initial site  

  TVTrackHit* pDummyHit = NULL;

  if ( (pDummyHit = dynamic_cast<ILDCylinderHit *>( kalhit )) ) {
    pDummyHit = (new ILDCylinderHit(*static_cast<ILDCylinderHit*>( kalhit )));
  }
  else if ( (pDummyHit = dynamic_cast<ILDPlanarHit *>( kalhit )) ) {
    pDummyHit = (new ILDPlanarHit(*static_cast<ILDPlanarHit*>( kalhit )));
  }
  else {
    streamlog_out( ERROR) << "<<<<<<<<< MarlinKalTestTrack::initialise: dynamic_cast failed for hit type >>>>>>>" << std::endl;
    return false;
  }

  TVTrackHit& dummyHit = *pDummyHit;

  //SJA:FIXME: this constants should go in a header file
  // give the dummy hit huge errors so that it does not contribute to the fit
  dummyHit(0,1) = 1.e6;   // give a huge error to d
  dummyHit(1,1) = 1.e6;   // give a huge error to z   
  
  // use dummy hit to create initial site
  TKalTrackSite& initialSite = *new TKalTrackSite(dummyHit);
  
  initialSite.SetHitOwner();// site owns hit
  initialSite.SetOwner();   // site owns states

  // ---------------------------
  //  Set up initial track state 
  // ---------------------------
  
  static TKalMatrix initialState(kSdim,1) ;
  initialState(0,0) = helix.GetDrho() ;        // d0
  initialState(1,0) = helix.GetPhi0() ;        // phi0
  initialState(2,0) = helix.GetKappa() ;       // kappa
  initialState(3,0) = helix.GetDz();           // dz
  initialState(4,0) = helix.GetTanLambda() ;   // tan(lambda)
  if (kSdim == 6) initialState(5,0) = 0.;      // t0


  // ---------------------------
  //  Set up initial Covariance Matrix
  // ---------------------------
  
  TKalMatrix covK(kSdim,kSdim) ;  for(int i=0;i<5;++i) for(int j=0;j<5;++j) covK[i][j] = cov[i][j] ;
  if (kSdim == 6) covK(5,5) = 1.e6; // t0
  

  // Add initial states to the site 
  initialSite.Add(new TKalTrackState(initialState,covK,initialSite,TVKalSite::kPredicted));
  initialSite.Add(new TKalTrackState(initialState,covK,initialSite,TVKalSite::kFiltered));

  // add the initial site to the track: that is, give the track initial parameters and covariance 
  // matrix at the starting measurement layer
  _kaltrack->Add(&initialSite);
  
  _initialised = true ;
  return 0 ;

} 

int MarlinKalTestTrack::addAndFit( ILDVTrackHit* kalhit, double& chi2increment, TKalTrackSite*& site, double maxChi2Increment) {

  streamlog_out(DEBUG3) << "MarlinKalTestTrack::addAndFit called " << std::endl ;

  if ( ! _initialised ) {
    
    throw MarlinTrk::Exception("Track fit not initialised");   
    
  }

  TKalTrackSite* temp_site = new TKalTrackSite(*kalhit); // create new site for this hit

  KalTrackFilter filter( maxChi2Increment );
  
  temp_site->SetFilterCond( &filter ) ;

  if (!_kaltrack->AddAndFilter(*temp_site)) {        

    chi2increment = temp_site->GetDeltaChi2() ;
    // get the measurement layer of the current hit
    const ILDVMeasLayer* ml =  dynamic_cast<const ILDVMeasLayer*>( &(kalhit->GetMeasLayer() ) ) ;
    streamlog_out( DEBUG4 )  << "Kaltrack::fit : site discarded! at index : " << ml->GetIndex() << " for type " << ml->GetMLName() << " layer ID " << ml->getLayerID() << std::endl ;
    
    delete temp_site;  // delete site if filter step failed      
    
    return 1 ;

  }

  site = temp_site;
  chi2increment = site->GetDeltaChi2() ;

  return 0 ;

}

int MarlinKalTestTrack::addAndFit( EVENT::TrackerHit* trkhit, double& chi2increment, double maxChi2Increment) {
    
  const ILDVMeasLayer* ml = _ktest->findMeasLayer( trkhit ) ;
  ILDVTrackHit* kalhit = ml->ConvertLCIOTrkHit(trkhit) ;

  TKalTrackSite* site;
  int error = this->addAndFit( kalhit, chi2increment, site, maxChi2Increment);

  if( error != 0 ){
    return 1 ;
  }
  else {
    this->addHit( trkhit, kalhit, ml ) ; 
    _hit_used_for_sites[trkhit] = site ;
  }

  return 0 ;

}

int MarlinKalTestTrack::fit( bool fitDirection ) {

  // SJA:FIXME: what do we do about calling fit after we have already added hits and filtered
  // I guess this would created new sites when addAndFit is called 
  // one option would be to remove the sites 
  // need to check where the sites are stored ...  probably in the KalTrackSystem
  // 

  streamlog_out(DEBUG4) << "MarlinKalTestTrack::fit() called " << std::endl ;
  
  const bool gkDir = fitDirection ; 

  if ( ! _initialised ) {

    throw MarlinTrk::Exception("Track fit not initialised");   

  }

  // ---------------------------
  //  Prepare hit iterrator for adding hits to kaltrack
  // ---------------------------

  TIter next(_kalhits, gkDir); // fit inwards to IP, if gkDir = kIterBackward

  // ---------------------------
  //  Start Kalman Filter
  // ---------------------------

  ILDVTrackHit *kalhit = 0;
  
  while ( (kalhit = dynamic_cast<ILDVTrackHit *>( next() ) ) ) {

    double chi2increment;
    TKalTrackSite* site;
    int error = this->addAndFit( kalhit, chi2increment, site);
    
    // here do dynamic cast repeatedly in DEBUG statement as this will be stripped out any way for production code
    // otherwise we have to do the cast outside of the DEBUG statement and it won't be stripped out 
    streamlog_out( DEBUG3 )  << "Kaltrack::fit :  add site to track at index : " 
			     << (dynamic_cast<const ILDVMeasLayer*>( &(kalhit->GetMeasLayer() ) ))->GetIndex() 
			     << " for type " 
			     << dynamic_cast<const ILDVMeasLayer*>( &(kalhit->GetMeasLayer() ) )->GetMLName() 
			     << " layer ID " 
			     << dynamic_cast<const ILDVMeasLayer*>( &(kalhit->GetMeasLayer() ) )->getLayerID() 
			     << " error = " << error 
			     << std::endl ;

    // find the lcio hit for this kaltest hit
    std::map<ILDVTrackHit*,EVENT::TrackerHit*>::iterator it;
    it = _kaltest_hits_to_lcio_hits.find(kalhit) ;

    if( it == _kaltest_hits_to_lcio_hits.end() ) { // something went wrong and this kalhit has no lcio trkhit associated
	  
      std::stringstream errorMsg;
      errorMsg << "MarlinKalTestTrack::fit hit pointer " << kalhit << " not stored in _kaltest_hits_to_lcio_hits map" << std::endl ; 
      throw MarlinTrk::Exception(errorMsg.str());
	
    }

    EVENT::TrackerHit* trkhit = it->second;  
    
    if( error == 0 ){ // add trkhit to map associating trkhits and sites
      _hit_used_for_sites[it->second] = site;
    } 
    else { // hit rejected by the filter, so store in the list of rejected hits
      _hit_not_used_for_sites.push_back(trkhit) ;
    }
        
  } // end of Kalman filter
  
  if( _ktest->getOption(  MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing ) ){
    streamlog_out( DEBUG3 )  << "Perform Smoothing for All Previous Measurement Sites " << std::endl ;
    _kaltrack->SmoothAll() ;
  }

  return 0;

}


int  MarlinKalTestTrack::getTrackState( IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {

  streamlog_out( DEBUG4 )  << "MarlinKalTestTrack::getTrackState( IMPL::TrackStateImpl& ts ) " << std::endl ;
  
  // use the last filtered track state 
  const TVKalSite& site = *(dynamic_cast<const TVKalSite*>(_kaltrack->Last())) ;

  this->ToLCIOTrackState( site, ts, chi2, ndf );
    
  return 0 ;

}


int MarlinKalTestTrack::getTrackState( EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {

  streamlog_out( DEBUG4 )  << "MarlinKalTestTrack::getTrackState( EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts ) using hit: " << trkhit << std::endl ;

  std::map<EVENT::TrackerHit*,TKalTrackSite*>::iterator it;

  int error = getSiteFromLCIOHit(trkhit, it);

  if( error != 0 ) return error;

  const TVKalSite& site = *(it->second) ; 
  
  this->ToLCIOTrackState( site, ts, chi2, ndf );

  return 0 ;
}

int MarlinKalTestTrack::extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){  

  const TVKalSite& site = *(dynamic_cast<const TVKalSite*>(_kaltrack->Last())) ;

  return this->extrapolate( point, site, ts, chi2, ndf ) ;

}

int MarlinKalTestTrack::extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {

  std::map<EVENT::TrackerHit*,TKalTrackSite*>::iterator it;
  int error = getSiteFromLCIOHit( trkhit, it ) ; 
  if( error !=0 ) return error;

  const TVKalSite& site = *(it->second);

  return this->extrapolate( point, site, ts, chi2, ndf ) ;

}

int MarlinKalTestTrack::extrapolate( const gear::Vector3D& point, const TVKalSite& site ,IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){  
  
  streamlog_out(DEBUG4) << "MarlinKalTestTrack::extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) called " << std::endl ;
  
  TKalTrackState& trkState = (TKalTrackState&) site.GetCurState(); // this segfaults if no hits are present
  
  THelicalTrack helix = trkState.GetHelix() ;
  double dPhi ;

  // convert the gear point supplied to TVector3
  const TVector3 tpoint( point.x(), point.y(), point.z() ) ;

  Int_t sdim = trkState.GetDimension();  // dimensions of the track state, it will be 5 or 6
  TKalMatrix sv(sdim,1);

  // now move to the point
  TKalMatrix  DF(sdim,sdim);  
  DF.UnitMatrix();                           
  helix.MoveTo(  tpoint , dPhi , &DF , 0) ;  // move helix to desired point, and get propagator matrix

  TMatrixD c0(trkState.GetCovMat());  

  TKalMatrix DFt  = TKalMatrix(TMatrixD::kTransposed, DF);
  c0 = DF * c0 * DFt ;                 // update the covariance matrix 

  this->ToLCIOTrackState( helix, c0, ts, chi2, ndf );

  return 0;
  
} 


int MarlinKalTestTrack::extrapolateToLayer( bool direction, int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) { 

  const TVKalSite& site = *(dynamic_cast<const TVKalSite*>(_kaltrack->Last())) ;
  
  return this->extrapolateToLayer( direction, layerID, site, ts, chi2, ndf, detElementID ) ;

}


int MarlinKalTestTrack::extrapolateToLayer( bool direction, int layerID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) { 

 std::map<EVENT::TrackerHit*,TKalTrackSite*>::iterator it;
 int error = getSiteFromLCIOHit( trkhit, it ) ; 
 if( error !=0 ) return error;
 
 const TVKalSite& site = *(it->second);
 
 return this->extrapolateToLayer( direction, layerID, site, ts, chi2, ndf, detElementID ) ;
 
}


int MarlinKalTestTrack::extrapolateToLayer( bool direction, int layerID, const TVKalSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) { 

  streamlog_out(DEBUG4) << "MarlinKalTestTrack::extrapolateToLayer( bool direction, int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) called " << std::endl ;

  gear::Vector3D crossing_point ;

  int error = this->intersectionWithLayer( direction, layerID, site, crossing_point, detElementID ) ;
  
  if( error != 0 ) return error ;

  return this->extrapolate( crossing_point, site, ts, chi2, ndf ) ;

} 



int MarlinKalTestTrack::propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){

  const TVKalSite& site = *(dynamic_cast<const TVKalSite*>(_kaltrack->Last())) ;

  return this->propagate( point, site, ts, chi2, ndf ) ;

}

int MarlinKalTestTrack::propagate( const gear::Vector3D& point, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){

  std::map<EVENT::TrackerHit*,TKalTrackSite*>::iterator it;
  int error = getSiteFromLCIOHit( trkhit, it ) ; 
  if( error !=0 ) return error;

  const TVKalSite& site = *(it->second);

  return this->propagate( point, site, ts, chi2, ndf ) ;

}

int MarlinKalTestTrack::propagate( const gear::Vector3D& point, const TVKalSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){

  streamlog_out(DEBUG4) << "MarlinKalTestTrack::propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) called " << std::endl ;

  // convert the gear point supplied to TVector3
  const TVector3 tpoint( point.x(), point.y(), point.z() ) ;

  TKalTrackState& trkState = (TKalTrackState&) site.GetCurState(); // this segfaults if no hits are present

  THelicalTrack helix = trkState.GetHelix() ;
  double dPhi ;

  // the last layer crossed by the track before point 
  const ILDVMeasLayer* ml = _ktest->getLastMeasLayer(helix, tpoint);
  
  Int_t sdim = trkState.GetDimension();  // dimensions of the track state, it will be 5 or 6
  TKalMatrix sv(sdim,1);

  TKalMatrix  F(sdim,sdim);              // propagator matrix to be returned by transport function
  F.UnitMatrix();                        // set the propagator matrix to the unit matrix

  TKalMatrix  Q(sdim,sdim);              // noise matrix to be returned by transport function 
  Q.Zero();        
  TVector3    x0;                        // intersection point to be returned by transport

  const TVMeasLayer& tvml   = dynamic_cast<const TVMeasLayer&>(*ml) ;         // cast from ILDVMeasurementLayer
  const TKalTrackSite& track_site = dynamic_cast<const TKalTrackSite&>(site) ;   // cast from TVKalSite
  
  _ktest->_det->Transport(track_site, tvml, x0, sv, F, Q ) ;      // transport to last layer cross before point 
  
  // given that we are sure to have intersected the layer tvml as this was provided via getLastMeasLayer, x0 will lie on the layer
  // this could be checked with the method isOnSurface 
  // so F will be the propagation matrix from the current location to the last surface and Q will be the noise matrix up to this point 
  TMatrixD c0(trkState.GetCovMat());  

  TKalMatrix Ft  = TKalMatrix(TMatrixD::kTransposed, F);
  c0 = F * c0 * Ft + Q; // update covaraince matrix and add the MS assosiated with moving to tvml
  
  helix.MoveTo(  x0 , dPhi , 0 , 0 ) ;  // move the helix to tvml

  // get whether the track is incomming or outgoing at the last surface
  const TVSurface *sfp = dynamic_cast<const TVSurface *>(ml);   // last surface
  TMatrixD dxdphi = helix.CalcDxDphi(0);                        // tangent vector at last surface                       
  TVector3 dxdphiv(dxdphi(0,0),dxdphi(1,0),dxdphi(2,0));        // convert matirix diagonal to vector
  Double_t cpa = helix.GetKappa();                              // get pt 
  
  Bool_t isout = -cpa*dxdphiv.Dot(sfp->GetOutwardNormal(x0)) < 0 ? kTRUE : kFALSE;  // out-going or in-coming at the destination surface

  // now move to the point
  TKalMatrix  DF(sdim,sdim);  
  DF.UnitMatrix();                           
  helix.MoveTo(  tpoint , dPhi , &DF , 0) ;  // move helix to desired point, and get propagator matrix
  
  TKalMatrix Qms(sdim, sdim);                                       
  tvml.CalcQms(isout, helix, dPhi, Qms);     // calculate MS for the final step through the present material 

  TKalMatrix DFt  = TKalMatrix(TMatrixD::kTransposed, DF);
  c0 = DF * c0 * DFt + Qms ;                 // update the covariance matrix 

  this->ToLCIOTrackState( helix, c0, ts, chi2, ndf );

  return 0;

}


int MarlinKalTestTrack::propagateToLayer( bool direction, int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) { 

  const TVKalSite& site = *(dynamic_cast<const TVKalSite*>(_kaltrack->Last())) ;
  
  return this->propagateToLayer( direction, layerID, site, ts, chi2, ndf, detElementID) ;

}


int MarlinKalTestTrack::propagateToLayer( bool direction, int layerID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) { 

 std::map<EVENT::TrackerHit*,TKalTrackSite*>::iterator it;
 int error = getSiteFromLCIOHit( trkhit, it ) ; 
 if( error !=0 ) return error;
 
 const TVKalSite& site = *(it->second);
 
 return this->propagateToLayer( direction, layerID, site, ts, chi2, ndf, detElementID) ;
 
}


int MarlinKalTestTrack::propagateToLayer( bool direction, int layerID, const TVKalSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) { 

  streamlog_out(DEBUG4) << "MarlinKalTestTrack::propagateToLayer( bool direction, int layerID, const TVKalSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) called " << std::endl;

  gear::Vector3D crossing_point ;

  int error = this->intersectionWithLayer( direction, layerID, site, crossing_point, detElementID) ;
  
  if( error != 0 ) return error ;

  return this->propagate( crossing_point, site, ts, chi2, ndf) ;

} 



int MarlinKalTestTrack::intersectionWithLayer( bool direction, int layerID, gear::Vector3D& point, int& detElementID) {  

  const TVKalSite& site = *(dynamic_cast<const TVKalSite*>(_kaltrack->Last())) ;

  return this->intersectionWithLayer( direction, layerID, site, point, detElementID) ;

}

int MarlinKalTestTrack::intersectionWithLayer( bool direction, int layerID,  EVENT::TrackerHit* trkhit, gear::Vector3D& point, int& detElementID) {  

  std::map<EVENT::TrackerHit*,TKalTrackSite*>::iterator it;
  int error = getSiteFromLCIOHit( trkhit, it ) ; 
  if( error !=0 ) return error;
  
  const TVKalSite& site = *(it->second);
  return this->intersectionWithLayer( direction, layerID, site, point, detElementID) ;

}

int MarlinKalTestTrack::intersectionWithLayer( bool direction, int layerID, const TVKalSite& site, gear::Vector3D& point, int& detElementID) {  

  streamlog_out(DEBUG4) << "MarlinKalTestTrack::intersectionWithLayer( bool direction, int layerID, const TVKalSite& site, gear::Vector3D& point, int& detElementID) called " << std::endl;

  std::vector<ILDVMeasLayer*> meas_modules ;
  _ktest->getSensitiveMeasurementModulesForLayer( layerID, meas_modules ) ;  

  TVSurface* surf = NULL;
  TVector3 xto;       // reference point at destination to be returned by CalcXinPointWith  
  Double_t dphi_min = DBL_MAX;  // use to store the min deflection angle found so that can avoid the crossing on the far side of the layer

  if( meas_modules.size() == 0 ) {
    
    std::stringstream errorMsg;
    errorMsg << "MarlinKalTestTrack::MarlinKalTestTrack layer id unkown: layerID = " << layerID << std::endl ; 
    throw MarlinTrk::Exception(errorMsg.str());
    
  } 
  else { 
    
    bool surf_found(false);

    for( unsigned int i=0; i < meas_modules.size(); ++i) {
     
      streamlog_out(DEBUG2) << "MarlinKalTestTrack::intersectionWithLayer try and find intersection for detElementID = " << (meas_modules[i])->getLayerID() << " layerID " << layerID << std::endl ;

      if( ! (surf = dynamic_cast<TVSurface*> (  meas_modules[i] )) )   {
      
	std::stringstream errorMsg;
	errorMsg << "MarlinKalTestTrack::MarlinKalTestTrack dynamic_cast failed for surface type: layerID = " << layerID << std::endl ; 
	throw MarlinTrk::Exception(errorMsg.str());

      }
      
      TKalTrackState& trkState = (TKalTrackState&) site.GetCurState(); // this segfaults if no hits are present
      
      THelicalTrack helix = trkState.GetHelix() ;
           
      TVector3 pivot = helix.GetPivot() ;

      Double_t dphi = 0;  // deflection angle to destination to be returned by CalcXingPointWith

      int crossing_exist = surf->CalcXingPointWith(helix, xto, dphi, direction) ;

      streamlog_out(DEBUG2) << "MarlinKalTestTrack::intersectionWithLayer crossing_exist = " << crossing_exist << " dphi " << dphi << std::endl ;
      
      // make sure we get the next crossing 
      if( crossing_exist > 0 && dphi < dphi_min ) { 
	
	surf_found = true ;
	dphi_min = dphi ;
	point[0] = xto.X();
	point[1] = xto.Y();
	point[2] = xto.Z();
	detElementID = meas_modules[i]->getLayerID() ;
      }

    }

    if( ! surf_found ){
      streamlog_out(DEBUG3) << "MarlinKalTestTrack::intersectionWithLayer no intersection found with layerID = "
			    << layerID
			    << std::endl ;
      return -1;
    }
   

    streamlog_out(DEBUG3) << "MarlinKalTestTrack::intersectionWithLayer intersection with layerID = "
			  << layerID
			  << ": at x = " << point.x()
			  << " y = "     << point.y()
			  << " z = "     << point.z()
			  << " detElementID = " << detElementID 
			  << std::endl ;

  }

  return 0;

} 



void MarlinKalTestTrack::ToLCIOTrackState( const THelicalTrack& helix, const TMatrixD& cov, IMPL::TrackStateImpl& ts, double& chi2, int& ndf){

  chi2 = _kaltrack->GetChi2();
  ndf  = _kaltrack->GetNDF();
 
  //============== convert parameters to LCIO convention ====
  
  // fill 5x5 covariance matrix from the 6x6 covariance matrix return by trkState.GetCovMat()  above
  TMatrixD covK(5,5) ;  for(int i=0;i<5;++i) for(int j=0;j<5;++j) covK[i][j] = cov[i][j] ;
  
  //  this is for incomming tracks ...
  double phi       =    toBaseRange( helix.GetPhi0() + M_PI/2. ) ;
  double omega     =    1. /helix.GetRho()  ;              
  double d0        =  - helix.GetDrho() ; 
  double z0        =    helix.GetDz()   ;
  double tanLambda =    helix.GetTanLambda()  ;

  ts.setD0( d0 ) ;  
  ts.setPhi( phi  ) ; // fi0  - M_PI/2.  ) ;  
  ts.setOmega( omega  ) ;
  ts.setZ0( z0  ) ;  
  ts.setTanLambda( tanLambda ) ;  
    
  Double_t cpa  = helix.GetKappa();
  double alpha = omega / cpa  ; // conversion factor for omega (1/R) to kappa (1/Pt) 

  EVENT::FloatVec covLCIO( 15 )  ; 
  covLCIO[ 0] =   covK( 0 , 0 )   ; //   d0,   d0

  covLCIO[ 1] = - covK( 1 , 0 )   ; //   phi0, d0
  covLCIO[ 2] =   covK( 1 , 1 )   ; //   phi0, phi

  covLCIO[ 3] = - covK( 2 , 0 ) * alpha   ; //   omega, d0
  covLCIO[ 4] =   covK( 2 , 1 ) * alpha   ; //   omega, phi
  covLCIO[ 5] =   covK( 2 , 2 ) * alpha * alpha  ; //   omega, omega

  covLCIO[ 6] = - covK( 3 , 0 )   ; //   z0  , d0
  covLCIO[ 7] =   covK( 3 , 1 )   ; //   z0  , phi
  covLCIO[ 8] =   covK( 3 , 2 ) * alpha   ; //   z0  , omega
  covLCIO[ 9] =   covK( 3 , 3 )   ; //   z0  , z0

  covLCIO[10] = - covK( 4 , 0 )   ; //   tanl, d0
  covLCIO[11] =   covK( 4 , 1 )   ; //   tanl, phi
  covLCIO[12] =   covK( 4 , 2 ) * alpha  ; //   tanl, omega
  covLCIO[13] =   covK( 4 , 3 )   ; //   tanl, z0
  covLCIO[14] =   covK( 4 , 4 )   ; //   tanl, tanl

 
  ts.setCovMatrix( covLCIO ) ;


  float pivot[3] ;

  pivot[0] =  helix.GetPivot().X() ;
  pivot[1] =  helix.GetPivot().Y() ;
  pivot[2] =  helix.GetPivot().Z() ;

  ts.setReferencePoint( pivot ) ;

  streamlog_out( DEBUG4 ) << " kaltest track parameters: "
			 << " chi2/ndf " << chi2 / ndf  
    			 << " chi2 " <<  chi2 << std::endl 
    
    			 << "\t D0 "          <<  d0         <<  "[+/-" << sqrt( covLCIO[0] ) << "] " 
			 << "\t Phi :"        <<  phi        <<  "[+/-" << sqrt( covLCIO[2] ) << "] " 
			 << "\t Omega "       <<  omega      <<  "[+/-" << sqrt( covLCIO[5] ) << "] " 
			 << "\t Z0 "          <<  z0         <<  "[+/-" << sqrt( covLCIO[9] ) << "] " 
			 << "\t tan(Lambda) " <<  tanLambda  <<  "[+/-" << sqrt( covLCIO[14]) << "] " 
    
			 << "\t pivot : [" << pivot[0] << ", " << pivot[1] << ", "  << pivot[2] 
			 << " - r: " << std::sqrt( pivot[0]*pivot[0]+pivot[1]*pivot[1] ) << "]" 
			 << std::endl ;


}


void MarlinKalTestTrack::ToLCIOTrackState( const TVKalSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {

  TKalTrackState& trkState = (TKalTrackState&) site.GetCurState(); // GetCutState will return the last added state to this site
  // Assuming everything has proceeded as expected 
  // this will be Predicted -> Filtered -> Smoothed 
  
  THelicalTrack helix = trkState.GetHelix() ;
  
  TMatrixD c0(trkState.GetCovMat());  
  
  this->ToLCIOTrackState( helix, c0, ts, chi2, ndf );

}


int MarlinKalTestTrack::getSiteFromLCIOHit( EVENT::TrackerHit* trkhit, std::map<EVENT::TrackerHit*,TKalTrackSite*>::iterator& it ){

  it = _hit_used_for_sites.find(trkhit) ;  

  if( it == _hit_used_for_sites.end() ) { // hit not associated with any site
    
    bool found = false;
    
    for( unsigned int i = 0; i < _hit_not_used_for_sites.size(); ++i) {
      if( trkhit == _hit_not_used_for_sites[i] ) found = true ;
    }

    if( found ) {
      streamlog_out( DEBUG4 )  << "MarlinKalTestTrack::getSiteFromLCIOHit: hit was rejected during filtering" << std::endl ;
      return 1 ;
    }
    else {
      return 2 ; // 1 hit was rejected by fit, 2 hit not associated with this track 
    }
  } 

  streamlog_out( DEBUG3 )  << "MarlinKalTestTrack::getSiteFromLCIOHit: site found" << std::endl ;
  return 0 ;
  
}
