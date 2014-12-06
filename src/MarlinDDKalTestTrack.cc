#include "MarlinTrk/MarlinDDKalTestTrack.h"

#include "MarlinTrk/MarlinDDKalTest.h"
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
#include <UTIL/Operators.h>
#include <UTIL/ILDConf.h>

#include "DDKalTest/DDCylinderMeasLayer.h" // needed for dedicated IP Layer
#include "DDKalTest/DDCylinderHit.h"
#include "DDKalTest/DDPlanarHit.h"
//#include "DDKalTest/DDPlanarStripHit.h"

#include "gear/GEAR.h"
#include "gear/BField.h"

#include "streamlog/streamlog.h"


/** Helper class for defining a filter condition based on the delta chi2 in the AddAndFilter step.
 */
class KalTrackFilter : public TKalFilterCond{
  
public:
  
  /** C'tor - takes as optional argument the maximum allowed delta chi2 for adding the hit (in IsAccepted() )
   */
  KalTrackFilter(double maxDeltaChi2 = DBL_MAX) : _maxDeltaChi2( maxDeltaChi2 ), _passed_last_filter_step(true) {
  } 
  virtual ~KalTrackFilter() {} 
  
  virtual Bool_t IsAccepted(const TKalTrackSite &site) {
    
    double deltaChi2 = site.GetDeltaChi2();
    
    streamlog_out( DEBUG1 ) << " KalTrackFilter::IsAccepted called  !  deltaChi2 = "  << std::scientific <<  deltaChi2  << " _maxDeltaChi2 = " << _maxDeltaChi2 << std::endl;

    _passed_last_filter_step = deltaChi2 < _maxDeltaChi2;
        
    return ( _passed_last_filter_step )   ; 
  }
  
  void resetFilterStatus() { _passed_last_filter_step = true; }
  bool passedLastFilterStep() const { return _passed_last_filter_step; }
  
protected:
  
  double _maxDeltaChi2 ;
  bool _passed_last_filter_step;
  
} ;

namespace MarlinTrk {
  
  //---------------------------------------------------------------------------------------------------------------
  
  std::string toString( int detElementID ) {
    lcio::BitField64 bf(  UTIL::ILDCellID0::encoder_string ) ;
    bf.setValue( detElementID ) ;
    return bf.valueString() ;
  }
  
  //---------------------------------------------------------------------------------------------------------------
  
  
  MarlinDDKalTestTrack::MarlinDDKalTestTrack( MarlinDDKalTest* ktest) 
  : _ktest(ktest) {
    
    _kaltrack = new TKalTrack() ;
    _kaltrack->SetOwner() ;
    
    _kalhits = new TObjArray() ;
    _kalhits->SetOwner() ;
    
    _initialised = false ;
    _fitDirection = false ;
    _smoothed = false ;
    
    _trackHitAtPositiveNDF = 0;
    _hitIndexAtPositiveNDF = 0;
    
#ifdef MARLINTRK_DIAGNOSTICS_ON
      _ktest->_diagnostics.new_track(this) ;
#endif
    
    
  }
  
  
  MarlinDDKalTestTrack::~MarlinDDKalTestTrack(){
    
#ifdef MARLINTRK_DIAGNOSTICS_ON    
    _ktest->_diagnostics.end_track() ;
#endif
    
    delete _kaltrack ;
    delete _kalhits ;
  }
  
  
  
  int MarlinDDKalTestTrack::addHit( EVENT::TrackerHit * trkhit) {

    return this->addHit( trkhit, _ktest->findMeasLayer( trkhit )) ;

  } 
  
  int MarlinDDKalTestTrack::addHit( EVENT::TrackerHit * trkhit, const DDVMeasLayer* ml) {

    streamlog_out(DEBUG1) << "MarlinDDKalTestTrack::addHit: trkhit = "  << trkhit->id() << " addr: " << trkhit << " ml = " << ml << std::endl ;
    
    if( trkhit && ml ) {
      return this->addHit( trkhit, ml->ConvertLCIOTrkHit(trkhit), ml) ;
    }
    else {
      streamlog_out( ERROR ) << " MarlinDDKalTestTrack::addHit - bad inputs " <<  trkhit << " ml : " << ml << std::endl ;
       return bad_intputs ;
    }
    
  }
  
  int MarlinDDKalTestTrack::addHit( EVENT::TrackerHit* trkhit, DDVTrackHit* kalhit, const DDVMeasLayer* ml) {
    
    if( kalhit && ml ) {
      _kalhits->Add(kalhit ) ;  // Add hit and set surface found 
      _lcio_hits_to_kaltest_hits[trkhit] = kalhit ; // add hit to map relating lcio and kaltest hits
                                                    //    _kaltest_hits_to_lcio_hits[kalhit] = trkhit ; // add hit to map relating kaltest and lcio hits
    }
    else {
      delete kalhit;
      return bad_intputs ;
    }
    
    streamlog_out(DEBUG1) << "MarlinDDKalTestTrack::addHit: hit added " 
    << "number of hits for track = " << _kalhits->GetEntries() 
    << std::endl ;
    
    return success ;
    
  }
  
  
  int MarlinDDKalTestTrack::initialise( bool fitDirection ) {; 
    
    
    //SJA:FIXME: check here if the track is already initialised, and for now don't allow it to be re-initialised
    //           if the track is going to be re-initialised then we would need to do it directly on the first site
    if ( _initialised ) {
      
      throw MarlinTrk::Exception("Track fit already initialised");   
      
    }
    
    if (_kalhits->GetEntries() < 3) {
      
      streamlog_out( ERROR) << "<<<<<< MarlinDDKalTestTrack::initialise: Shortage of Hits! nhits = "  
      << _kalhits->GetEntries() << " >>>>>>>" << std::endl;
      return error ;
      
    }
    
    _fitDirection =  fitDirection ; 
    
    // establish the hit order
    Int_t i1, i2, i3; // (i1,i2,i3) = (1st,mid,last) hit to filter
    if (_fitDirection == kIterBackward) {
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
    
    TVTrackHit* pDummyHit = 0;
    
    if ( (pDummyHit = dynamic_cast<DDCylinderHit *>( startingHit )) ) {
      pDummyHit = (new DDCylinderHit(*static_cast<DDCylinderHit*>( startingHit )));
    }
    else if ( (pDummyHit = dynamic_cast<DDPlanarHit *>( startingHit )) ) {
      pDummyHit = (new DDPlanarHit(*static_cast<DDPlanarHit*>( startingHit )));
    }
    // else if ( DDPlanarStripHit_DIM == 2 && (pDummyHit = dynamic_cast<DDPlanarStripHit *>( startingHit )) ) {
    //   pDummyHit = (new DDPlanarStripHit(*static_cast<DDPlanarStripHit*>( startingHit )));
    // }
    else {
      streamlog_out( ERROR) << "<<<<<<<<< MarlinDDKalTestTrack::initialise: dynamic_cast failed for hit type >>>>>>>" << std::endl;
      return error ;
    }
    
    TVTrackHit& dummyHit = *pDummyHit;
    
    //SJA:FIXME: this constants should go in a header file
    // give the dummy hit huge errors so that it does not contribute to the fit
    dummyHit(0,1) = 1.e16;   // give a huge error to d
    dummyHit(1,1) = 1.e16;   // give a huge error to z   
    
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
    
    if ( h1.GetDimension() == 1 || h2.GetDimension() == 1 || h3.GetDimension() == 1  ) {
      
      throw MarlinTrk::Exception("Track fit cannot be initialised from 1 Dimentional hits. Use method MarlinDDKalTestTrack::initialise(  const EVENT::TrackState& ts, double bfield_z, bool fitDirection )");   
      
    }
    
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::initialise Create initial helix from hits: \n "
    << "P1 x = " << x1.x() << " y = " << x1.y() << " z = " << x1.z() << " r = " << x1.Perp() << "\n "
    << "P2 x = " << x2.x() << " y = " << x2.y() << " z = " << x2.z() << " r = " << x2.Perp() << "\n "
    << "P3 x = " << x3.x() << " y = " << x3.y() << " z = " << x3.z() << " r = " << x3.Perp() << "\n "
    << "Bz = " << h1.GetBfield() << " direction = " << _fitDirection
    << std::endl;
    
    // create helix using 3 global space points 
    THelicalTrack helstart(x1, x2, x3, h1.GetBfield(), _fitDirection); // initial helix 
    
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
    
    TKalMatrix Cov(kSdim,kSdim);

    // make sure everything is initialised to zero
    for (int i=0; i<kSdim*kSdim; ++i) {
      Cov.GetMatrixArray()[i] = 0.0;
    }
    
//    for (Int_t i=0; i<kSdim; i++) {
//      // fg: if the error is too large the initial helix parameters might be changed extremely by the first three (or so) hits,
//      //     such that the fit will not work because the helix curls away and does not hit the next layer !!!
//      Cov(i,i) = 1.e2 ;   // initialise diagonal elements of dummy error matrix
//    }

    // prefer translation over rotation of the trackstate early in the fit 
    
    Cov(0,0) = 1.e6 ; // d0
    Cov(1,1) = 1.e2 ; // dphi0
    Cov(2,2) = 1.e2 ; // dkappa
    Cov(3,3) = 1.e6 ; // dz
    Cov(4,4) = 1.e2 ; // dtanL
    if (kSdim == 6) Cov(5,5) = 1.e2;  // t0
          
    // Add initial states to the site 
    initialSite.Add(new TKalTrackState(initialState,Cov,initialSite,TVKalSite::kPredicted));
    initialSite.Add(new TKalTrackState(initialState,Cov,initialSite,TVKalSite::kFiltered));
    
    // add the initial site to the track: that is, give the track initial parameters and covariance 
    // matrix at the starting measurement layer
    _kaltrack->Add(&initialSite);
    
    _initialised = true ;
    
    streamlog_out( DEBUG2 ) << " track parameters used for init : " << std::scientific << std::setprecision(6) 
			    << "\t D0 "          <<  0.0
			    << "\t Phi :"        <<  toBaseRange( helstart.GetPhi0() + M_PI/2. )
			    << "\t Omega "       <<  1. /helstart.GetRho() 
			    << "\t Z0 "          <<  0.0
			    << "\t tan(Lambda) " <<  helstart.GetTanLambda()
      
			    << "\t pivot : [" << helstart.GetPivot().X() << ", " << helstart.GetPivot().Y() << ", "  << helstart.GetPivot().Z()
			    << " - r: " << std::sqrt( helstart.GetPivot().X()*helstart.GetPivot().X()+helstart.GetPivot().Y()*helstart.GetPivot().Y() ) << "]"
			    << std::endl ;
    
#ifdef MARLINTRK_DIAGNOSTICS_ON
    
    // convert to LICO parameters first
    
    double d0        =    0.0 ;
    double phi       =    toBaseRange( helstart.GetPhi0() + M_PI/2. );
    double omega     =    1. /helstart.GetRho()  ;              
    double z0        =    0.0 ;
    double tanLambda =    helstart.GetTanLambda()  ;
    
//    Cov.Print();
    
    _ktest->_diagnostics.set_intial_track_parameters(d0,
                                                     phi,
                                                     omega,
                                                     z0,
                                                     tanLambda,
                                                     helstart.GetPivot().X(),
                                                     helstart.GetPivot().Y(),
                                                     helstart.GetPivot().Z(),
                                                     Cov);

    
    
#endif
    

    return success ;
    
  }
  
  int MarlinDDKalTestTrack::initialise(  const EVENT::TrackState& ts, double bfield_z, bool fitDirection ) {

    if (_kalhits->GetEntries() == 0) {
      
      streamlog_out( ERROR) << "<<<<<< MarlinDDKalTestTrack::Initialise: Number of Hits is Zero. Cannot Initialise >>>>>>>" << std::endl;
      return error ;
      
    }
    
    //SJA:FIXME: check here if the track is already initialised, and for now don't allow it to be re-initialised
    //           if the track is going to be re-initialised then we would need to do it directly on the first site
    if ( _initialised ) {
      
      throw MarlinTrk::Exception("Track fit already initialised");   
      
    }
        
    streamlog_out( DEBUG2 ) << "MarlinDDKalTestTrack::initialise using TrackState: track parameters used for init : "
    << "\t D0 "          <<  ts.getD0()         
    << "\t Phi :"        <<  ts.getPhi()        
    << "\t Omega "       <<  ts.getOmega()      
    << "\t Z0 "          <<  ts.getZ0()         
    << "\t tan(Lambda) " <<  ts.getTanLambda()  
    
    << "\t pivot : [" << ts.getReferencePoint()[0] << ", " << ts.getReferencePoint()[1] << ", "  << ts.getReferencePoint()[2] 
    << " - r: " << std::sqrt( ts.getReferencePoint()[0]*ts.getReferencePoint()[0]+ts.getReferencePoint()[1]*ts.getReferencePoint()[1] ) << "]" 
    << std::endl ;

    
    _fitDirection = fitDirection ;
    
    // for GeV, Tesla, R in mm  
    double alpha = bfield_z * 2.99792458E-4 ;
    double kappa;
    if ( bfield_z == 0.0 )
      kappa = DBL_MAX;
    else kappa = ts.getOmega() / alpha ;
    
    THelicalTrack helix( -ts.getD0(),
                        toBaseRange( ts.getPhi() - M_PI/2. ) ,
                        kappa,
                        ts.getZ0(),
                        ts.getTanLambda(),
                        ts.getReferencePoint()[0],
                        ts.getReferencePoint()[1],
                        ts.getReferencePoint()[2],
                        bfield_z );
    
    TMatrixD cov(5,5) ;   

    EVENT::FloatVec covLCIO = ts.getCovMatrix();
    
    cov( 0 , 0 ) =   covLCIO[ 0] ;                   //   d0, d0
    cov( 0 , 1 ) = - covLCIO[ 1] ;                   //   d0, phi
    cov( 0 , 2 ) = - covLCIO[ 3] / alpha ;           //   d0, kappa
    cov( 0 , 3 ) = - covLCIO[ 6] ;                   //   d0, z0
    cov( 0 , 4 ) = - covLCIO[10] ;                   //   d0, tanl
    
    cov( 1 , 0 ) = - covLCIO[ 1] ;                   //   phi, d0
    cov( 1 , 1 ) =   covLCIO[ 2] ;                   //   phi, phi
    cov( 1 , 2 ) =   covLCIO[ 4] / alpha ;           //   phi, kappa
    cov( 1 , 3 ) =   covLCIO[ 7] ;                   //   phi, z0
    cov( 1 , 4 ) =   covLCIO[11] ;                   //   tanl, phi
    
    cov( 2 , 0 ) = - covLCIO[ 3] / alpha ;           //   kappa, d0
    cov( 2 , 1 ) =   covLCIO[ 4] / alpha ;           //   kappa, phi
    cov( 2 , 2 ) =   covLCIO[ 5] / (alpha * alpha) ; //   kappa, kappa
    cov( 2 , 3 ) =   covLCIO[ 8] / alpha ;           //   kappa, z0
    cov( 2 , 4 ) =   covLCIO[12] / alpha ;           //   kappa, tanl
    
    cov( 3 , 0 ) = - covLCIO[ 6] ;                   //   z0, d0
    cov( 3 , 1 ) =   covLCIO[ 7] ;                   //   z0, phi
    cov( 3 , 2 ) =   covLCIO[ 8] / alpha ;           //   z0, kappa
    cov( 3 , 3 ) =   covLCIO[ 9] ;                   //   z0, z0
    cov( 3 , 4 ) =   covLCIO[13] ;                   //   z0, tanl
    
    cov( 4 , 0 ) = - covLCIO[10] ;                   //   tanl, d0 
    cov( 4 , 1 ) =   covLCIO[11] ;                   //   tanl, phi
    cov( 4 , 2 ) =   covLCIO[12] / alpha ;           //   tanl, kappa    
    cov( 4 , 3 ) =   covLCIO[13] ;                   //   tanl, z0
    cov( 4 , 4 ) =   covLCIO[14] ;                   //   tanl, tanl
    
//    cov.Print();
    
    // move the helix to either the position of the last hit or the first depending on initalise_at_end
    
    // default case initalise_at_end
    int index = _kalhits->GetEntries() - 1 ;
    // or initialise at start 
    if( _fitDirection == IMarlinTrack::forward ){
      index = 0 ;
    }
    
    TVTrackHit* kalhit = dynamic_cast<TVTrackHit *>(_kalhits->At(index)); 
    
    double dphi;
   
    TVector3 initial_pivot ;

    // Leave the pivot at the origin for a 1-dim hit
    if (kalhit->GetDimension() > 1) {
      
      initial_pivot = kalhit->GetMeasLayer().HitToXv(*kalhit) ;
    }
    else{
      initial_pivot =  TVector3(0.0,0.0,0.0);
    }

    
    // ---------------------------
    //  Create an initial start site for the track using the  hit
    // ---------------------------
    // set up a dummy hit needed to create initial site  
    
    TVTrackHit* pDummyHit = 0;
  
    if ( (pDummyHit = dynamic_cast<DDCylinderHit *>( kalhit )) ) {
      pDummyHit = (new DDCylinderHit(*static_cast<DDCylinderHit*>( kalhit )));
      
    }
    else if ( (pDummyHit = dynamic_cast<DDPlanarHit *>( kalhit )) ) {
      pDummyHit = (new DDPlanarHit(*static_cast<DDPlanarHit*>( kalhit )));
      //}
      //else if ( (pDummyHit = dynamic_cast<DDPlanarStripHit *>( kalhit )) ) {
      //pDummyHit = (new DDPlanarStripHit(*static_cast<DDPlanarStripHit*>( kalhit )));
      if( pDummyHit->GetDimension() == 1 ) {
      
	const TVMeasLayer *ml = &pDummyHit->GetMeasLayer();
	
	const TVSurface* surf = dynamic_cast<const TVSurface*>(ml);
	
	if (surf) {
	  double phi;
	  
	  surf->CalcXingPointWith(helix, initial_pivot, phi);        
	  
	} else {
	  streamlog_out( ERROR) << "<<<<<<<<< MarlinDDKalTestTrack::initialise: dynamic_cast failed for TVSurface  >>>>>>>" << std::endl;
	  return error ;        
	}
	
      }
    }
    else {
      streamlog_out( ERROR) << "<<<<<<<<< MarlinDDKalTestTrack::initialise: dynamic_cast failed for hit type >>>>>>>" << std::endl;
      return error ;
    }
    
    TVTrackHit& dummyHit = *pDummyHit;
        
    //SJA:FIXME: this constants should go in a header file
    // give the dummy hit huge errors so that it does not contribute to the fit
    dummyHit(0,1) = 1.e16;   // give a huge error to d
    
    if(dummyHit.GetDimension()>1) dummyHit(1,1) = 1.e16;   // give a huge error to z   
    
    // use dummy hit to create initial site
    TKalTrackSite& initialSite = *new TKalTrackSite(dummyHit);
    
    initialSite.SetHitOwner();// site owns hit
    initialSite.SetOwner();   // site owns states
    
    // ---------------------------
    //  Set up initial track state 
    // ---------------------------

    helix.MoveTo( initial_pivot, dphi, 0, &cov );  
    
    static TKalMatrix initialState(kSdim,1) ;
    initialState(0,0) = helix.GetDrho() ;        // d0
    initialState(1,0) = helix.GetPhi0() ;        // phi0
    initialState(2,0) = helix.GetKappa() ;       // kappa
    initialState(3,0) = helix.GetDz();           // dz
    initialState(4,0) = helix.GetTanLambda() ;   // tan(lambda)
    if (kSdim == 6) initialState(5,0) = 0.;      // t0
    
    // make sure that the pivot is in the right place 
    initialSite.SetPivot(initial_pivot);
    
    // ---------------------------
    //  Set up initial Covariance Matrix
    // ---------------------------
    
    TKalMatrix covK(kSdim,kSdim) ;  
    for(int i=0;i<5;++i) {
      for(int j=0;j<5;++j) {
        covK[i][j] = cov[i][j] ;  
      }
    }
    if (kSdim == 6) covK(5,5) = 1.e6; // t0
        
//    covK.Print();
    
    // Add initial states to the site 
    initialSite.Add(new TKalTrackState(initialState,covK,initialSite,TVKalSite::kPredicted));
    initialSite.Add(new TKalTrackState(initialState,covK,initialSite,TVKalSite::kFiltered));
    
    // add the initial site to the track: that is, give the track initial parameters and covariance 
    // matrix at the starting measurement layer
    _kaltrack->Add(&initialSite);
    
    _initialised = true ;

    
#ifdef MARLINTRK_DIAGNOSTICS_ON
    
    // convert to LICO parameters first
    
    double d0        =  - helix.GetDrho() ;
    double phi       =    toBaseRange( helix.GetPhi0() + M_PI/2. );
    double omega     =    1. /helix.GetRho()  ;              
    double z0        =    helix.GetDz()   ;
    double tanLambda =    helix.GetTanLambda()  ;
        
    _ktest->_diagnostics.set_intial_track_parameters(d0,
                                                     phi,
                                                     omega,
                                                     z0,
                                                     tanLambda,
                                                     helix.GetPivot().X(),
                                                     helix.GetPivot().Y(),
                                                     helix.GetPivot().Z(),
                                                     covK);

    

#endif

    return success ;
    
  } 
  
  int MarlinDDKalTestTrack::addAndFit( DDVTrackHit* kalhit, double& chi2increment, TKalTrackSite*& site, double maxChi2Increment) {
    
    streamlog_out(DEBUG1) << "MarlinDDKalTestTrack::addAndFit called : maxChi2Increment = "  << std::scientific << maxChi2Increment << std::endl ;
    
    if ( ! _initialised ) {
      
      throw MarlinTrk::Exception("Track fit not initialised");   
      
    }
    
    // here do dynamic cast repeatedly in DEBUG statement as this will be stripped out any way for production code
    // otherwise we have to do the cast outside of the DEBUG statement and it won't be stripped out 
    streamlog_out( DEBUG1 )  << "Kaltrack::addAndFit :  add site to track at index : " 
    << (dynamic_cast<const DDVMeasLayer*>( &(kalhit->GetMeasLayer() ) ))->GetIndex() 
    << " for type " 
    << dynamic_cast<const DDVMeasLayer*>( &(kalhit->GetMeasLayer() ) )->GetName() ;
    streamlog_out( DEBUG0 ) << " with CellIDs:";
    
    for (unsigned int i = 0; i < (dynamic_cast<const DDVMeasLayer*>( &(kalhit->GetMeasLayer() ) )->getNCellIDs());++i) {
      streamlog_out( DEBUG0 )  << " : " 
      << dynamic_cast<const DDVMeasLayer*>( &(kalhit->GetMeasLayer() ) )->getCellIDs()[i] ;
      
    }
    
    streamlog_out( DEBUG1 ) << std::endl ;
    
    TKalTrackSite* temp_site = new TKalTrackSite(*kalhit); // create new site for this hit
    
    KalTrackFilter filter( maxChi2Increment );
    filter.resetFilterStatus();
    
    temp_site->SetFilterCond( &filter ) ;
    
    
    // this is the only point at which a hit is actually filtered 
    // and it is here that we can get the GetDeltaChi2 vs the maxChi2Increment
    // it will always be possible to get the delta chi2 so long as we have a link to the sites ...
    // although calling smooth will natrually update delta chi2.
    
    
    if (!_kaltrack->AddAndFilter(*temp_site)) {        
      
      chi2increment = temp_site->GetDeltaChi2() ;
      // get the measurement layer of the current hit
      const DDVMeasLayer* ml =  dynamic_cast<const DDVMeasLayer*>( &(kalhit->GetMeasLayer() ) ) ;
      TVector3 pos = ml->HitToXv(*kalhit);
      streamlog_out( DEBUG1 )  << "Kaltrack::addAndFit : site discarded! at index : " << ml->GetIndex()
      << " for type " << ml->GetName() 
      << " chi2increment = " << chi2increment
      << " maxChi2Increment = " << maxChi2Increment
      << " x = " << pos.x()
      << " y = " << pos.y()
      << " z = " << pos.z()
      << " with CellIDs: " << std::endl;
      
      for (unsigned int i = 0; i < (dynamic_cast<const DDVMeasLayer*>( &(kalhit->GetMeasLayer() ) )->getNCellIDs());++i) {
        streamlog_out( DEBUG1 )  << " CellID = " 
        << dynamic_cast<const DDVMeasLayer*>( &(kalhit->GetMeasLayer() ) )->getCellIDs()[i] 
        << std::endl ;
      }
      

#ifdef MARLINTRK_DIAGNOSTICS_ON
      _ktest->_diagnostics.record_rejected_site(kalhit, temp_site); 
#endif
      
      delete temp_site;  // delete site if filter step failed      
      
      
      // compiling the code below with the cmake option CMAKE_BUILD_TYPE=Debug
      // and with LDFLAGS=-Wl,--no-undefined
      // causes an undefined reference error
      // the problem gets fixed using the if/else statement below
      
      // this fails
      //return filter.usedForLastFilterStep() ? site_fails_chi2_cut : site_discarded ;
      
      // this also fails..
      //bool rc = filter.usedForLastFilterStep() ;
      //return (rc ? site_fails_chi2_cut : site_discarded);
      
      // but this works ?!! 
      //return ( true ? site_fails_chi2_cut : site_discarded);
      
      // and this also works..
      streamlog_out(DEBUG2) << " addAndFit : Site Fails Chi2 cut ? " << filter.passedLastFilterStep() << std::endl;

      if( filter.passedLastFilterStep() == false ) {
        return site_fails_chi2_cut ;
      } else {
        return site_discarded ;
      }
      
    }
    
    site = temp_site;
    chi2increment = site->GetDeltaChi2() ;

#ifdef MARLINTRK_DIAGNOSTICS_ON
    _ktest->_diagnostics.record_site(kalhit, site);  
#endif
    
    return success ;
    
  }
  
  int MarlinDDKalTestTrack::addAndFit( EVENT::TrackerHit* trkhit, double& chi2increment, double maxChi2Increment) {
    
    if( ! trkhit ) { 
      streamlog_out( ERROR) << "MarlinDDKalTestTrack::addAndFit( EVENT::TrackerHit* trkhit, double& chi2increment, double maxChi2Increment): trkhit == 0"  << std::endl;
      return bad_intputs ; 
    }
    
    const DDVMeasLayer* ml = _ktest->findMeasLayer( trkhit ) ;
    
    if( ml == 0 ){  
      // fg: not sure if ml should ever be 0 - but it seems to happen, 
      //     if point is not on surface and more than one surface exists ...
      
      streamlog_out( ERROR ) << ">>>>>>>>>>>  no measurment layer found for trkhit cellid0 : " 
      << toString( trkhit->getCellID0() ) << " at " 
      << gear::Vector3D( trkhit->getPosition() ) << std::endl ;
      
      return  IMarlinTrack::bad_intputs ; 
    }
    
    DDVTrackHit* kalhit = ml->ConvertLCIOTrkHit(trkhit) ;
    
    if( kalhit == 0 ){  //fg: ml->ConvertLCIOTrkHit returns 0 if hit not on surface !!!
      return IMarlinTrack::bad_intputs ;
    }
    
    TKalTrackSite* site = 0 ;
    int error_code = this->addAndFit( kalhit, chi2increment, site, maxChi2Increment);
    
    if( error_code != success ){

      delete kalhit;

      // if the hit fails for any reason other than the Chi2 cut record the Chi2 contibution as DBL_MAX
      if( error_code != site_fails_chi2_cut ) {
        chi2increment = DBL_MAX;
      }

      _outlier_chi2_values.push_back(std::make_pair(trkhit, chi2increment));

      streamlog_out( DEBUG2 ) << ">>>>>>>>>>>  addAndFit Number of Outliers : "
      << _outlier_chi2_values.size() << std::endl;
      
      return error_code ;
    }
    else {
      this->addHit( trkhit, kalhit, ml ) ; 
      _hit_used_for_sites[trkhit] = site ;
      _hit_chi2_values.push_back(std::make_pair(trkhit, chi2increment));
    }
    
    // set the values for the point at which the fit becomes constained 
    if( _trackHitAtPositiveNDF == 0 && _kaltrack->GetNDF() >= 0){

      _trackHitAtPositiveNDF = trkhit;
      _hitIndexAtPositiveNDF = _kaltrack->IndexOf( site );
      
      streamlog_out( DEBUG2 ) << ">>>>>>>>>>>  Fit is now constrained at : "
      << toString( trkhit->getCellID0() ) 
      << " pos " << gear::Vector3D( trkhit->getPosition() )
      << " trkhit = " << _trackHitAtPositiveNDF
      << " index of kalhit = " << _hitIndexAtPositiveNDF
      << " NDF = " << _kaltrack->GetNDF() 
      <<  std::endl; 
    
    }

    return success ;
    
  }
  
  
  
  int MarlinDDKalTestTrack::testChi2Increment( EVENT::TrackerHit* trkhit, double& chi2increment ) {
    
    if( ! trkhit ) { 
      streamlog_out( ERROR) << "MarlinDDKalTestTrack::addAndFit( EVENT::TrackerHit* trkhit, double& chi2increment, double maxChi2Increment): trkhit == 0"  << std::endl;
      return IMarlinTrack::bad_intputs ; 
    }
    
    const DDVMeasLayer* ml = _ktest->findMeasLayer( trkhit ) ;
    
    if( ml == 0 ){  
      // fg: not sure if ml should ever be 0 - but it seems to happen, 
      //     if point is not on surface and more than one surface exists ...
      
      streamlog_out( ERROR ) << ">>>>>>>>>>>  no measurment layer found for trkhit cellid0 : " 
      << toString( trkhit->getCellID0() ) << " at " 
      << gear::Vector3D( trkhit->getPosition() ) << std::endl ;
      
      return  IMarlinTrack::bad_intputs ; 
      
    }
    
    DDVTrackHit* kalhit = ml->ConvertLCIOTrkHit(trkhit) ;
    
    if( kalhit == 0 ){  //fg: ml->ConvertLCIOTrkHit returns 0 if hit not on surface !!!
      return IMarlinTrack::bad_intputs ;
    }
    
    
    TKalTrackSite* site = 0 ;
    int error_code = this->addAndFit( kalhit, chi2increment, site, -DBL_MAX); // using -DBL_MAX here ensures the hit will never be added to the fit
    
    delete kalhit;  
    
    return error_code;
    
  }
  
  
  
  int MarlinDDKalTestTrack::fit( double maxChi2Increment ) {
    
    // SJA:FIXME: what do we do about calling fit after we have already added hits and filtered
    // I guess this would created new sites when addAndFit is called 
    // one option would be to remove the sites 
    // need to check where the sites are stored ...  probably in the KalTrackSystem
    // 
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::fit() called " << std::endl ;
    
    if ( ! _initialised ) {
      
      throw MarlinTrk::Exception("Track fit not initialised");   
      
    }
    
    // ---------------------------
    //  Prepare hit iterrator for adding hits to kaltrack
    // ---------------------------
    
    TIter next(_kalhits, _fitDirection); 
    
    // ---------------------------
    //  Start Kalman Filter
    // ---------------------------
    
    DDVTrackHit *kalhit = 0;
    
    while ( (kalhit = dynamic_cast<DDVTrackHit *>( next() ) ) ) {
      
      double chi2increment;
      TKalTrackSite* site;
      int error_code = this->addAndFit( kalhit, chi2increment, site, maxChi2Increment );
      
      
      EVENT::TrackerHit* trkhit = kalhit->getLCIOTrackerHit();
      
      if( error_code == 0 ){ // add trkhit to map associating trkhits and sites
        _hit_used_for_sites[trkhit] = site;
        _hit_chi2_values.push_back(std::make_pair(trkhit, chi2increment));

        // set the values for the point at which the fit becomes constained 
        if( _trackHitAtPositiveNDF == 0 && _kaltrack->GetNDF() >= 0){
          
          _trackHitAtPositiveNDF = trkhit;
          _hitIndexAtPositiveNDF = _kaltrack->IndexOf( site );
          
          streamlog_out( DEBUG2 ) << ">>>>>>>>>>>  Fit is now constrained at : "
          << toString( trkhit->getCellID0() ) 
          << " pos " << gear::Vector3D( trkhit->getPosition() )
          << " trkhit = " << _trackHitAtPositiveNDF
          << " index of kalhit = " << _hitIndexAtPositiveNDF
          << " NDF = " << _kaltrack->GetNDF() 
          <<  std::endl; 
          
        }
            
      } 
      else { // hit rejected by the filter, so store in the list of rejected hits

        // if the hit fails for any reason other than the Chi2 cut record the Chi2 contibution as DBL_MAX
        if( error_code != site_fails_chi2_cut ) {
          chi2increment = DBL_MAX;
        }
        
        _outlier_chi2_values.push_back(std::make_pair(trkhit, chi2increment));
        streamlog_out( DEBUG2 ) << ">>>>>>>>>>>  fit(): Number of Outliers : "
        << _outlier_chi2_values.size() << std::endl;

        _hit_not_used_for_sites.push_back(trkhit) ;
        
      }
      
    } // end of Kalman filter
    
    if( _ktest->getOption(  MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing ) ){
      streamlog_out( DEBUG2 )  << "Perform Smoothing for All Previous Measurement Sites " << std::endl ;
      int error = this->smooth() ;
      
      if( error != success ) return error ;
      
    }
    
    //return _hit_used_for_sites.empty() == false ? success : all_sites_fail_fit ;
    if( _hit_used_for_sites.empty() == false )
    {
        return success ;
    }
    else{
        return all_sites_fail_fit ;
    }
    
  }
  
  
  /** smooth all track states 
   */
  int MarlinDDKalTestTrack::smooth(){
    
    streamlog_out( DEBUG2 )  << "MarlinDDKalTestTrack::smooth() " << std::endl ;
    
    //fg: we should actually smooth all sites - it is then up to the user which smoothed tracks state to take 
    //    for any furthter extrapolation/propagation ...
 
    if( !_smoothed ) 
      _kaltrack->SmoothAll() ;
    
    //SJA:FIXME: in the current implementation it is only possible to smooth back to the 4th site.
    // This is due to the fact that the covariance matrix is not well defined at the first 3 measurement sites filtered.
    
    //    _kaltrack->SmoothBackTo( _hitIndexAtPositiveNDF + 1 ) ;
    
   _smoothed = true ;

    return success ;
    
  }
  
  
  /** smooth track states from the last filtered hit back to the measurement site associated with the given hit 
   */
  int MarlinDDKalTestTrack::smooth( EVENT::TrackerHit* trkhit ) {
    
    streamlog_out( DEBUG2 )  << "MarlinDDKalTestTrack::smooth( EVENT::TrackerHit* " << trkhit << "  ) " << std::endl ;

    
    if ( !trkhit ) {
      return bad_intputs ;
    }
    
    
    std::map<EVENT::TrackerHit*,TKalTrackSite*>::const_iterator it;
    
    
    TKalTrackSite* site = 0 ;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code ;
    
    int index = _kaltrack->IndexOf( site );
    
    _kaltrack->SmoothBackTo( index ) ;
    
    _smoothed = true ;
    
    return success ;
    
  }
  
  
  int  MarlinDDKalTestTrack::getTrackState( IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {
    
    streamlog_out( DEBUG2 )  << "MarlinDDKalTestTrack::getTrackState( IMPL::TrackStateImpl& ts ) " << std::endl ;
    
    // use the last filtered track state 
    const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    this->ToLCIOTrackState( site, ts, chi2, ndf );
    
    return success ;
    
    
  }
  
  
  int MarlinDDKalTestTrack::getTrackState( EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {
    
    streamlog_out( DEBUG2 )  << "MarlinDDKalTestTrack::getTrackState( EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts ) using hit: " << trkhit << " with cellID0 = " << trkhit->getCellID0() << std::endl ;
    
    TKalTrackSite* site = 0 ;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code ;
    
    streamlog_out( DEBUG1 )  << "MarlinDDKalTestTrack::getTrackState: site " << site << std::endl;
    
    this->ToLCIOTrackState( *site, ts, chi2, ndf );
    
    return success ;
  }
  
  
  int MarlinDDKalTestTrack::getHitsInFit( std::vector<std::pair<EVENT::TrackerHit*, double> >& hits ) {
    
    std::copy( _hit_chi2_values.begin() , _hit_chi2_values.end() , std::back_inserter(  hits  )  ) ;

    // this needs more thought. What about when the hits are added using addAndFit?

    // need to check the order so that we can return the list ordered in time
    // as they will be added to _hit_chi2_values in the order of fitting 
    // not in the order of time
//    
//    if( _fitDirection == IMarlinTrack::backward ){    
//      std::reverse_copy( _hit_chi2_values.begin() , _hit_chi2_values.end() , std::back_inserter(  hits  )  ) ;
//    } else {
//      std::copy( _hit_chi2_values.begin() , _hit_chi2_values.end() , std::back_inserter(  hits  )  ) ;
//    }
    
    return success ;
    
  }
  
  int MarlinDDKalTestTrack::getOutliers( std::vector<std::pair<EVENT::TrackerHit*, double> >& hits ) {

    std::copy( _outlier_chi2_values.begin() , _outlier_chi2_values.end() , std::back_inserter(  hits  )  ) ;

    
    // this needs more thought. What about when the hits are added using addAndFit?
//    // need to check the order so that we can return the list ordered in time
//    // as they will be added to _hit_chi2_values in the order of fitting 
//    // not in the order of time
//
//    if( _fitDirection == IMarlinTrack::backward ){    
//      std::reverse_copy( _outlier_chi2_values.begin() , _outlier_chi2_values.end() , std::back_inserter(  hits  )  ) ;
//    } else {
//      std::copy( _outlier_chi2_values.begin() , _outlier_chi2_values.end() , std::back_inserter(  hits  )  ) ;
//    }

        
    return success ;
    
  }
  
  
  int MarlinDDKalTestTrack::getNDF( int& ndf ){
    
    if( _initialised == false ) { 
      return error;
    } else {
      
      ndf = _kaltrack->GetNDF();
      return success;
      
    }
    
  }
  
  int MarlinDDKalTestTrack::getTrackerHitAtPositiveNDF( EVENT::TrackerHit*& trkhit ) {

      trkhit = _trackHitAtPositiveNDF;
      return success;    

  }
  
  
  int MarlinDDKalTestTrack::extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){  
    
    const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    return this->extrapolate( point, site, ts, chi2, ndf ) ;
    
  }
  
  int MarlinDDKalTestTrack::extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) {
    
    TKalTrackSite* site = 0 ;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code;
    
    return this->extrapolate( point, *site, ts, chi2, ndf ) ;
    
  }
  
  int MarlinDDKalTestTrack::extrapolate( const gear::Vector3D& point, const TKalTrackSite& site ,IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){  
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) called " << std::endl ;
    
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
    
    return success;
    
  } 
  
  
  int MarlinDDKalTestTrack::extrapolateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    return this->extrapolateToLayer( layerID, site, ts, chi2, ndf, detElementID, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::extrapolateToLayer( int layerID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    TKalTrackSite* site = 0;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code ;
    
    return this->extrapolateToLayer( layerID, *site, ts, chi2, ndf, detElementID, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::extrapolateToLayer( int layerID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::extrapolateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) called " << std::endl ;
    
    gear::Vector3D crossing_point ;
    const DDVMeasLayer* ml = 0;
    
    int error_code = this->intersectionWithLayer( layerID, site, crossing_point, detElementID, ml, mode ) ;
    
    if( error_code != 0 ) return error_code ;
    
    return this->extrapolate( crossing_point, site, ts, chi2, ndf ) ;
    
  } 
  
  
  int MarlinDDKalTestTrack::extrapolateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    return this->extrapolateToDetElement( detElementID, site, ts, chi2, ndf, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::extrapolateToDetElement( int detElementID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    TKalTrackSite* site = 0;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code ;
    
    return this->extrapolateToDetElement( detElementID, *site, ts, chi2, ndf, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::extrapolateToDetElement( int detElementID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::extrapolateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) called " << std::endl ;
    
    gear::Vector3D crossing_point ;
    
    const DDVMeasLayer* ml = 0;
    int error_code = this->intersectionWithDetElement( detElementID, site, crossing_point, ml, mode ) ;
    
    if( error_code != 0 ) return error_code ;
    
    return this->extrapolate( crossing_point, site, ts, chi2, ndf ) ;
    
  } 
  
  
  
  int MarlinDDKalTestTrack::propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){
    
    const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    // check if the point is inside the beampipe
    // SJA:FIXME: really we should also check if the PCA to the point is also less than R
    const DDVMeasLayer* ml = (  (_ktest->getIPLayer() ) && point.r() < _ktest->getIPLayer()->GetR()) ? _ktest->getIPLayer() : 0;
    
    return this->propagate( point, site, ts, chi2, ndf, ml ) ;
    
  }
  
  int MarlinDDKalTestTrack::propagate( const gear::Vector3D& point, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){
    
    TKalTrackSite* site = 0;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code ;
    
    // check if the point is inside the beampipe
    // SJA:FIXME: really we should also check if the PCA to the point is also less than R

    const DDVMeasLayer* ml = _ktest->getIPLayer();

    if ( ml )
      if (point.r() > _ktest->getIPLayer()->GetR()) ml = NULL;

//    const DDVMeasLayer* ml = (point.r() < _ktest->getIPLayer()->GetR()) ? _ktest->getIPLayer() : 0;
    
    return this->propagate( point, *site, ts, chi2, ndf, ml ) ;
    
  }
  
  int MarlinDDKalTestTrack::propagate( const gear::Vector3D& point, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, const DDVMeasLayer* ml ){
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::propagate( const gear::Vector3D& point, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) called " << std::endl ;
    
    // convert the gear point supplied to TVector3
    const TVector3 tpoint( point.x(), point.y(), point.z() ) ;
    
    TKalTrackState& trkState = (TKalTrackState&) site.GetCurState(); // this segfaults if no hits are present
    
    THelicalTrack helix = trkState.GetHelix() ;
    double dPhi = 0.0;
    
    
    Int_t sdim = trkState.GetDimension();  // dimensions of the track state, it will be 5 or 6
    TKalMatrix sv(sdim,1);
    
    TKalMatrix  F(sdim,sdim);              // propagator matrix to be returned by transport function
    F.UnitMatrix();                        // set the propagator matrix to the unit matrix
    
    TKalMatrix  Q(sdim,sdim);              // noise matrix to be returned by transport function 
    Q.Zero();        
    TVector3    x0;                        // intersection point to be returned by transport
    
    TMatrixD c0(trkState.GetCovMat());  
    
    // the last layer crossed by the track before point 
    if( ! ml ){
      ml = _ktest->getLastMeasLayer(helix, tpoint);
    }
    
    if ( ml ) {
      
      _ktest->_det->Transport(site, *ml, x0, sv, F, Q ) ;      // transport to last layer cross before point 
      
      // given that we are sure to have intersected the layer ml as this was provided via getLastMeasLayer, x0 will lie on the layer
      // this could be checked with the method isOnSurface 
      // so F will be the propagation matrix from the current location to the last surface and Q will be the noise matrix up to this point 
      
      
      TKalMatrix Ft  = TKalMatrix(TMatrixD::kTransposed, F);
      c0 = F * c0 * Ft + Q; // update covaraince matrix and add the MS assosiated with moving to tvml
      
      helix.MoveTo(  x0 , dPhi , 0 , 0 ) ;  // move the helix to tvml
      
      
    }
    else { // the current site is at the last surface before the point to propagate to 
      
      ml = dynamic_cast<const DDVMeasLayer*>(&(site.GetHit().GetMeasLayer())) ;
      
    }
    
    // get whether the track is incomming or outgoing at the last surface
    const TVSurface *sfp = dynamic_cast<const TVSurface *>(ml);   // last surface
    
    TMatrixD dxdphi = helix.CalcDxDphi(0);                        // tangent vector at last surface                       
    TVector3 dxdphiv(dxdphi(0,0),dxdphi(1,0),dxdphi(2,0));        // convert matirix diagonal to vector
//    Double_t cpa = helix.GetKappa();                              // get pt 

    Bool_t isout = -dPhi*dxdphiv.Dot(sfp->GetOutwardNormal(x0)) < 0 ? kTRUE : kFALSE;  // out-going or in-coming at the destination surface
    
    // now move to the point
    TKalMatrix  DF(sdim,sdim);  
    DF.UnitMatrix();                           
    helix.MoveTo(  tpoint , dPhi , &DF , 0) ;  // move helix to desired point, and get propagator matrix

    TKalMatrix Qms(sdim, sdim);
    ml->CalcQms(isout, helix, dPhi, Qms);     // calculate MS for the final step through the present material 
    
    TKalMatrix DFt  = TKalMatrix(TMatrixD::kTransposed, DF);
    c0 = DF * c0 * DFt + Qms ;                 // update the covariance matrix 
    
    
    this->ToLCIOTrackState( helix, c0, ts, chi2, ndf );
    
    return success;
    
  }
  
  
  int MarlinDDKalTestTrack::propagateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    return this->propagateToLayer( layerID, site, ts, chi2, ndf, detElementID, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::propagateToLayer( int layerID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    TKalTrackSite* site = 0;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code ;
    
    return this->propagateToLayer( layerID, *site, ts, chi2, ndf, detElementID, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::propagateToLayer( int layerID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode ) { 
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::propagateToLayer( int layerID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) called " << std::endl;
    
    gear::Vector3D crossing_point ;
    
    const DDVMeasLayer* ml = 0;
    
    int error_code = this->intersectionWithLayer( layerID, site, crossing_point, detElementID, ml, mode) ;
    
    if( error_code != success ) return error_code ;
    
    return this->propagate( crossing_point, site, ts, chi2, ndf , ml) ;
    
  } 
  
  
  int MarlinDDKalTestTrack::propagateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    return this->propagateToDetElement( detElementID, site, ts, chi2, ndf, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::propagateToDetElement( int detElementID, EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    TKalTrackSite* site = 0;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code ;
    
    return this->propagateToDetElement( detElementID, *site, ts, chi2, ndf, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::propagateToDetElement( int detElementID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) { 
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::propagateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode ) called " << std::endl ;
    
    gear::Vector3D crossing_point ;
    
    const DDVMeasLayer* ml = 0;
    int error_code = this->intersectionWithDetElement( detElementID, site, crossing_point, ml, mode ) ;
    
    if( error_code != 0 ) return error_code ;
    
    return this->propagate( crossing_point, site, ts, chi2, ndf, ml ) ;
    
  } 
  
  
  int MarlinDDKalTestTrack::intersectionWithDetElement( int detElementID, gear::Vector3D& point, int mode ) {  
    
    const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    
    const DDVMeasLayer* ml = 0;
    return this->intersectionWithDetElement( detElementID, site, point, ml, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::intersectionWithDetElement( int detElementID,  EVENT::TrackerHit* trkhit, gear::Vector3D& point, int mode ) {  
    
    TKalTrackSite* site = 0;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code ;
    
    const DDVMeasLayer* ml = 0;
    return this->intersectionWithDetElement( detElementID, *site, point, ml, mode ) ;
    
  }
  
  int MarlinDDKalTestTrack::intersectionWithDetElement( int detElementID, const TKalTrackSite& site, gear::Vector3D& point, const DDVMeasLayer*& ml, int mode ) {
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::intersectionWithDetElement( int detElementID, const TKalTrackSite& site, gear::Vector3D& point, const DDVMeasLayer*& ml, int mode) called " << std::endl;
    
    std::vector<const DDVMeasLayer*> meas_modules ;
    _ktest->getSensitiveMeasurementModules( detElementID, meas_modules ) ;  
    
    if( meas_modules.size() == 0 ) {
      
      std::stringstream errorMsg;
      errorMsg << "MarlinDDKalTestTrack::intersectionWithDetElement detector element id unkown: detElementID = " 
      << toString( detElementID )  << std::endl ; 
      
      throw MarlinTrk::Exception(errorMsg.str());
      
    } 
    
    int dummy_detElementID; // not returned here as this is a single element as far as the outside world is concerned. Could check they are equal if we wanted ...
    int error_code = this->findIntersection( meas_modules, site, point, dummy_detElementID, ml, mode ) ;
    
    if( error_code == success ){
      
      
      streamlog_out(DEBUG1) << "MarlinDDKalTestTrack::intersectionWithDetElement intersection with detElementID = "
      <<  toString( detElementID ) 
      << ": at x = " << point.x()
      << " y = "     << point.y()
      << " z = "     << point.z()
      << std::endl ;
    }
    
    else if( error_code == no_intersection ) {
      
      ml = 0;
      
      streamlog_out(DEBUG1) << "MarlinDDKalTestTrack::intersectionWithDetElement No intersection with detElementID = "
      << toString( detElementID )
      << std::endl ;
      
    }
    
    return error_code ;
    
  }
  
  int MarlinDDKalTestTrack::intersectionWithLayer( int layerID, gear::Vector3D& point, int& detElementID, int mode ) {  
    
    const TKalTrackSite& site = *(dynamic_cast<const TKalTrackSite*>(_kaltrack->Last())) ;
    const DDVMeasLayer* ml = 0;
    return this->intersectionWithLayer( layerID, site, point, detElementID, ml,  mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::intersectionWithLayer( int layerID,  EVENT::TrackerHit* trkhit, gear::Vector3D& point, int& detElementID, int mode ) {  
    
    TKalTrackSite* site = 0;
    int error_code = getSiteFromLCIOHit(trkhit, site);
    
    if( error_code != success ) return error_code ;
    
    const DDVMeasLayer* ml = 0;
    return this->intersectionWithLayer( layerID, *site, point, detElementID, ml, mode ) ;
    
  }
  
  
  int MarlinDDKalTestTrack::intersectionWithLayer( int layerID, const TKalTrackSite& site, gear::Vector3D& point, int& detElementID, const DDVMeasLayer*& ml, int mode ) {  
    
    streamlog_out(DEBUG2) << "MarlinDDKalTestTrack::intersectionWithLayer( int layerID, const TKalTrackSite& site, gear::Vector3D& point, int& detElementID, int mode) called layerID = " << layerID << std::endl;
    
    std::vector<DDVMeasLayer const*> meas_modules ;
    _ktest->getSensitiveMeasurementModulesForLayer( layerID, meas_modules ) ;  
    
    if( meas_modules.size() == 0 ) {
      
      streamlog_out(DEBUG5)<< "MarlinDDKalTestTrack::intersectionWithLayer layer id unknown: layerID = " << toString( layerID ) << std::endl ;
      return no_intersection;
      
    } 
    
    //  int index_of_intersected;
    int error_code = this->findIntersection( meas_modules, site, point, detElementID, ml, mode ) ;
    
    if( error_code == success ){
      
      
      streamlog_out(DEBUG1) << "MarlinDDKalTestTrack::intersectionWithLayer intersection with layerID = "
      << layerID
      << ": at x = " << point.x()
      << " y = "     << point.y()
      << " z = "     << point.z()
      << " r = "     << point.rho()
      << " detElementID = " << detElementID 
      << " " << toString( detElementID )
      << std::endl ;
      
    }
    else if( error_code == no_intersection ) {
      
      ml = 0;
      streamlog_out(DEBUG1) << "MarlinDDKalTestTrack::intersectionWithLayer No intersection with layerID = "
      << layerID 
      << " " << toString( layerID )
      << std::endl ;
      
    }
    
    return error_code ;
    
    
  } 
  
  
  int MarlinDDKalTestTrack::findIntersection( const DDVMeasLayer& meas_module, const TKalTrackSite& site, gear::Vector3D& point, double& dphi, int& detElementID, int mode ) {
    
    
    TKalTrackState& trkState = (TKalTrackState&) site.GetCurState(); 
    
    
    //--------- DEBUG --------------
    // TKalTrackState* tsSmoothed = (  &((TVKalSite&)site).GetState(TVKalSite::kSmoothed) != 0 ? 
    //                              &(TKalTrackState&) ((TVKalSite&)site).GetState( TVKalSite::kSmoothed )  : 0 ) ;
    // if( tsSmoothed == &trkState ) 
    //   streamlog_out(DEBUG2) << "************ MarlinDDKalTestTrack::intersectionWithLayer : using smoothed TrackState !!!!! " << std::endl ;
    
    // TKalTrackState* tsFiltered = (  &((TVKalSite&)site).GetState(TVKalSite::kFiltered) != 0 ? 
    //                              &(TKalTrackState&) ((TVKalSite&)site).GetState( TVKalSite::kFiltered )  : 0 ) ;
    // if( tsFiltered == &trkState ) 
    //   streamlog_out(DEBUG2) << "************ MarlinDDKalTestTrack::intersectionWithLayer : using filtered TrackState !!!!! " << std::endl ;
    // //------------------------------
    
    
    THelicalTrack helix = trkState.GetHelix() ;
    
    TVector3 xto;       // reference point at destination to be returned by CalcXinPointWith  
    
    int crossing_exist = meas_module.getIntersectionAndCellID(helix, xto, dphi, detElementID, mode);
    //  int crossing_exist = surf->CalcXingPointWith(helix, xto, dphi, mode) ;
    
    streamlog_out(DEBUG1) << "MarlinDDKalTestTrack::intersectionWithLayer crossing_exist = " << crossing_exist << " dphi " << dphi << " with detElementIDs: " <<  detElementID ;
    
    streamlog_out(DEBUG1) << std::endl ;
    
    
    if( crossing_exist == 0 ) { 
      return no_intersection ;
    }
    else {
      
      point[0] = xto.X();
      point[1] = xto.Y();
      point[2] = xto.Z();
      
    }
    
    return success ;
    
  }
  
  
  
  int MarlinDDKalTestTrack::findIntersection( std::vector<DDVMeasLayer const*>& meas_modules, const TKalTrackSite& site, gear::Vector3D& point, int& detElementID, const DDVMeasLayer*& ml, int mode ) {
    
    unsigned int n_modules = meas_modules.size() ;
    
    double dphi_min = DBL_MAX;  // use to store the min deflection angle found so that can avoid the crossing on the far side of the layer
    bool surf_found(false);
    
    for( unsigned int i = 0 ; i < n_modules ; ++i ){
      
      double dphi = 0;
      // need to send a temporary point as we may get the crossing point with the layer on the oposite side of the layer 
      gear::Vector3D point_temp ;
      
      int temp_detElementID;
      
      int error_code = findIntersection( *meas_modules[i], site, point_temp, dphi, temp_detElementID, mode ) ;
      
      if( error_code == success ) {
        
        // make sure we get the next crossing 
        if( fabs(dphi) < dphi_min ) {     
          
          dphi_min = fabs(dphi) ;
          surf_found = true ;
          ml = meas_modules[i];
          detElementID = temp_detElementID;
          point = point_temp ;
        }
        
      }
      else if( error_code != no_intersection ) { // in which case error_code is an error rather than simply a lack of intersection, so return  
        
        return error_code ;
        
      }
      
    }
    
    // check if the surface was found and return accordingly 
    if ( surf_found ) {
      return success ;
    }
    else {
      return no_intersection ;
    }
    
    
  }
  
  
  
  void MarlinDDKalTestTrack::ToLCIOTrackState( const THelicalTrack& helix, const TMatrixD& cov, IMPL::TrackStateImpl& ts, double& chi2, int& ndf) const {
    
    chi2 = _kaltrack->GetChi2();
    ndf  = _kaltrack->GetNDF();
    
    //============== convert parameters to LCIO convention ====
    
    // fill 5x5 covariance matrix from the 6x6 covariance matrix above
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
    
    streamlog_out( DEBUG2 ) << " kaltest track parameters: "
    << " chi2/ndf " << chi2 / ndf  
    << " chi2 " <<  chi2 
    << " ndf " <<  ndf
    << " prob " <<  TMath::Prob(chi2, ndf)
    << std::endl 
    
    << "\t D0 "          <<  d0         <<  "[+/-" << sqrt( covLCIO[0] ) << "] " 
    << "\t Phi :"        <<  phi        <<  "[+/-" << sqrt( covLCIO[2] ) << "] " 
    << "\t Omega "       <<  omega      <<  "[+/-" << sqrt( covLCIO[5] ) << "] " 
    << "\t Z0 "          <<  z0         <<  "[+/-" << sqrt( covLCIO[9] ) << "] " 
    << "\t tan(Lambda) " <<  tanLambda  <<  "[+/-" << sqrt( covLCIO[14]) << "] " 
    
    << "\t pivot : [" << pivot[0] << ", " << pivot[1] << ", "  << pivot[2] 
    << " - r: " << std::sqrt( pivot[0]*pivot[0]+pivot[1]*pivot[1] ) << "]" 
    << std::endl ;
    
    
  }
  
  
  void MarlinDDKalTestTrack::ToLCIOTrackState( const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) const {
    
    TKalTrackState& trkState = (TKalTrackState&) site.GetCurState(); // GetCutState will return the last added state to this site
                                                                     // Assuming everything has proceeded as expected 
                                                                     // this will be Predicted -> Filtered -> Smoothed 
    
    THelicalTrack helix = trkState.GetHelix() ;
    
    TMatrixD c0(trkState.GetCovMat());  
    
    this->ToLCIOTrackState( helix, c0, ts, chi2, ndf );
    
  }
  
  
  int MarlinDDKalTestTrack::getSiteFromLCIOHit( EVENT::TrackerHit* trkhit, TKalTrackSite*& site ) const {
    
    std::map<EVENT::TrackerHit*,TKalTrackSite*>::const_iterator it;
    
    it = _hit_used_for_sites.find(trkhit) ;  
    
    if( it == _hit_used_for_sites.end() ) { // hit not associated with any site
      
      bool found = false;
      
      for( unsigned int i = 0; i < _hit_not_used_for_sites.size(); ++i) {
        if( trkhit == _hit_not_used_for_sites[i] ) found = true ;
      }
      
      if( found ) {
        streamlog_out( DEBUG2 )  << "MarlinDDKalTestTrack::getSiteFromLCIOHit: hit was rejected during filtering" << std::endl ;
        return site_discarded ;
      }
      else {
        streamlog_out( DEBUG2 )  << "MarlinDDKalTestTrack::getSiteFromLCIOHit: hit " << trkhit << " not in list of supplied hits" << std::endl ;
        return bad_intputs ; 
      }
    } 
    
    site = it->second;
    
    
    streamlog_out( DEBUG1 )  << "MarlinDDKalTestTrack::getSiteFromLCIOHit: site " << site << " found for hit " << trkhit << std::endl ;
    return success ;
    
  }
  
} // end of namespace MarlinTrk 
