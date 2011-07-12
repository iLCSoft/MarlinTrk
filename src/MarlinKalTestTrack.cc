
#include "MarlinTrk/MarlinKalTestTrack.h"

#include "MarlinTrk/MarlinKalTest.h"

#include <kaltest/TKalTrack.h>
#include <kaltest/TKalTrackState.h>
#include "kaltest/TKalTrackSite.h"

#include <lcio.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>

#include <UTIL/BitField64.h>

#include <ILDCellIDEncoding.h>

#include "kaldet/ILDCylinderMeasLayer.h"
#include "kaldet/ILDCylinderHit.h"

#include "kaldet/ILDPlanarMeasLayer.h"
#include "kaldet/ILDPlanarHit.h"

#include "streamlog/streamlog.h"


MarlinKalTestTrack::MarlinKalTestTrack( MarlinKalTest* ktest) 
  : _ktest(ktest)
{

  _kaltrack = new TKalTrack() ;
  _kaltrack->SetOwner() ;

  _kalhits = new TObjArray() ;

}


MarlinKalTestTrack::~MarlinKalTestTrack(){
  delete _kaltrack ;
  delete _kalhits ;
}


void MarlinKalTestTrack::addHit( EVENT::TrackerHit * trkhit) 
{
  
  //  const ILDVMeasLayer* ml = _ktest->getSensitiveMeasurementLayer( trkhit->getCellID0() ) ;
  
  std::vector<ILDVMeasLayer*> measlayers ;
  _ktest->getSensitiveMeasurementLayer( trkhit->getCellID0(), measlayers ) ;

  if( measlayers.size() == 0 ) {
    
    std::stringstream errorMsg;
    errorMsg << "MarlinKalTestTrack::MarlinKalTestTrack hit layer id unkown: layerID = " << trkhit->getCellID0() << std::endl ; 
    throw MarlinTrk::Exception(errorMsg.str());
    
  } 
  else if (measlayers.size() == 1) {
    _kalhits->Add( measlayers[0]->ConvertLCIOTrkHit(trkhit) ) ; 
  }
  else { // layer has been split 
    
    bool surf_found(false);

    for( unsigned int i=0; i < measlayers.size(); ++i) {
     
      const TVector3 hit( trkhit->getPosition()[0], trkhit->getPosition()[1], trkhit->getPosition()[2]) ;
      
      TVSurface* surf = NULL;

      if( (surf = dynamic_cast<ILDCylinderMeasLayer*> (  measlayers[i] )) )   {
	// surf = dynamic_cast<ILDCylinderMeasLayer*> (  measlayers[i] )
      }
      else if ( (surf = dynamic_cast<ILDPlanarMeasLayer *>( measlayers[i] )) )  {
	// surf = dynamic_cast<ILDPlanarMeasLayer*> (  measlayers[i] )
      }
      else {
	std::stringstream errorMsg;
	errorMsg << "MarlinKalTestTrack::MarlinKalTestTrack dynamic_cast failed for surface type: layerID = " << trkhit->getCellID0() << std::endl ; 
	throw MarlinTrk::Exception(errorMsg.str());
      }
      
      bool hit_on_surface = surf->IsOnSurface(hit);

      if( (!surf_found) && hit_on_surface ){
	_kalhits->Add( measlayers[i]->ConvertLCIOTrkHit(trkhit) ) ;  // Add hit and set surface found 
	surf_found = true ;
      }
      else if( surf_found && hit_on_surface ) {  // only one surface should be found, if not throw 
	std::stringstream errorMsg;
	errorMsg << "MarlinKalTestTrack::MarlinKalTestTrack hit found to be on two surfaces: layerID = " << trkhit->getCellID0() << std::endl ; 
	throw MarlinTrk::Exception(errorMsg.str());
      }      
    }
    if( ! surf_found ){
      streamlog_out(DEBUG3) << "MarlinKalTestTrack::MarlinKalTestTrack hit not found to be on any surface matching layerID = "
			    << trkhit->getCellID0()
			    << ": x = " << trkhit->getPosition()[0]
			    << " y = " << trkhit->getPosition()[1]
			    << " z = " << trkhit->getPosition()[2]
			    << std::endl ;
    }

    streamlog_out(DEBUG3) << "MarlinKalTestTrack::MarlinKalTestTrack hit found to be on surface matching layerID = "
			  << trkhit->getCellID0()
			  << ": x = " << trkhit->getPosition()[0]
			  << " y = " << trkhit->getPosition()[1]
			  << " z = " << trkhit->getPosition()[2]
			  << std::endl ;
  }


  streamlog_out(DEBUG3) << "MarlinKalTestTrack::MarlinKalTestTrack hit added " 
			<< "number of hits for track = " << _kalhits->GetEntries() 
			<< std::endl ;

} 




bool MarlinKalTestTrack::fit( Bool_t fitDirection ) {

  //SJA:FIXME: do we need to sort the hits here ... ?

  const Bool_t gkDir = fitDirection ; 

  streamlog_out(DEBUG4) << "MarlinKalTestTrack::fit() called " << std::endl ;

  if (_kalhits->GetEntries() < 3) {
    
    streamlog_out( ERROR) << "<<<<<< KalTrack::fitTrack(): Shortage of Hits! nhits = "  
			  << _kalhits->GetEntries() << " >>>>>>>" << std::endl;
    return false;

  }

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
    //    pDummyHit = (new ILDCylinderHit(*static_cast<ILDTPCHit*>( startingHit )));
  }
  else if ( (pDummyHit = dynamic_cast<ILDPlanarHit *>( startingHit )) ) {
    //    pDummyHit = (new ILDPlanarHit(*static_cast<ILDTPCHit*>( startingHit )));
  }
  else {
    streamlog_out( ERROR) << "<<<<<< KalTrack::fitTrack(): dynamic_cast failed for hit type >>>>>>>" << std::endl;
    return false;
  }

  TVTrackHit &dummyHit = *pDummyHit;

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
    Cov(i,i) = 1.e6;   // dummy error matrix
  }


  // Add initial states to the site 
  initialSite.Add(new TKalTrackState(initialState,Cov,initialSite,TVKalSite::kPredicted));
  initialSite.Add(new TKalTrackState(initialState,Cov,initialSite,TVKalSite::kFiltered));

  // add the initial site to the track: that is, give the track initial parameters and covariance 
  // matrix at the starting measurement layer
  _kaltrack->Add(&initialSite);

  // ---------------------------
  //  Prepare hit iterrator for adding hits to kaltrack
  // ---------------------------

  TIter next(_kalhits, gkDir); // come in to IP, if gkDir = kIterBackward

  // ---------------------------
  //  Start Kalman Filter
  // ---------------------------

  TVTrackHit *hitp = 0;
  
  while ( (hitp = dynamic_cast<TVTrackHit *>( next() ) ) ) {
   
    TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create new site for this hit

    // get the measurement layer of the current hit
    const ILDVMeasLayer* ml =  dynamic_cast<const ILDVMeasLayer*>( &(hitp->GetMeasLayer() ) ) ;
    
    streamlog_out( DEBUG3 )  << "Kaltrack::fit :  add site to track at index : " << ml->GetIndex() << " for type " << ml->GetMLName() << " layer ID " << ml->getLayerID() << std::endl ;
    
    // try to add the site and filter 
    if (!_kaltrack->AddAndFilter(site)) {        
      streamlog_out( DEBUG4 )  << "Kaltrack::fit : site discarded! at index : " << ml->GetIndex() << " for type " << ml->GetMLName() << " layer ID " << ml->getLayerID() << std::endl ;
      delete &site;                        // delete it if failed      
    }
  } // end of Kalman filter
  
  return true;

}



IMPL::TrackStateImpl* MarlinKalTestTrack::propagateToIP(){

  streamlog_out(DEBUG4) << "MarlinKalTestTrack::PropagateToIP() called " << std::endl ;

  // create track state to be returned
  IMPL::TrackStateImpl* trk = new IMPL::TrackStateImpl();

  // get the current site. This will be the IP for now as that is the last hit we added ... assuming it did not fail SJA:FIXME: how should this be enforced?  
  TVKalSite& cursite = _kaltrack->GetCurSite();

  TKalTrackState& trkState = (TKalTrackState&) cursite.GetCurState(); // this segfaults if no hits are present

  Int_t    ndf  = _kaltrack->GetNDF();
  Double_t chi2 = _kaltrack->GetChi2();

  //============== convert parameters to LCIO convention ====
  
  //  ---- get parameters at origin 
  
  THelicalTrack helix = trkState.GetHelix() ;
  double dPhi ;

  // need to get the 5x5 sub matrix of the covariance matrix
  const TMatrixD& c0 =  trkState.GetCovMat() ;

  // fill 5x5 covariance matrix from the 6x6 covariance matrix return by trkState.GetCovMat()  above
  TMatrixD covK(5,5) ;  for(int i=0;i<5;++i) for(int j=0;j<5;++j) covK[i][j] = c0[i][j] ;

  helix.MoveTo(  TVector3( 0., 0., 0. ) , dPhi , 0 , &covK ) ;

  //  this is for incomming tracks ...
  double phi       =    toBaseRange( helix.GetPhi0() + M_PI/2. ) ;
  double omega     =    1. /helix.GetRho()  ;              
  double d0        =  - helix.GetDrho() ; 
  double z0        =    helix.GetDz()   ;
  double tanLambda =    helix.GetTanLambda()  ;

  trk->setD0( d0 ) ;  
  trk->setPhi( phi  ) ; // fi0  - M_PI/2.  ) ;  
  trk->setOmega( omega  ) ;
  trk->setZ0( z0  ) ;  
  trk->setTanLambda( tanLambda ) ;  
    
  Double_t cpa  = trkState(2, 0);
  double alpha = omega / cpa  ; // conversion factor for omega (1/R) to kappa (1/Pt) 

  EVENT::FloatVec cov( 15 )  ; 
  cov[ 0] =   covK( 0 , 0 )   ; //   d0,   d0

  cov[ 1] = - covK( 1 , 0 )   ; //   phi0, d0
  cov[ 2] =   covK( 1 , 1 )   ; //   phi0, phi

  cov[ 3] = - covK( 2 , 0 ) * alpha   ; //   omega, d0
  cov[ 4] =   covK( 2 , 1 ) * alpha   ; //   omega, phi
  cov[ 5] =   covK( 2 , 2 ) * alpha * alpha  ; //   omega, omega

  cov[ 6] = - covK( 3 , 0 )   ; //   z0  , d0
  cov[ 7] =   covK( 3 , 1 )   ; //   z0  , phi
  cov[ 8] =   covK( 3 , 2 ) * alpha   ; //   z0  , omega
  cov[ 9] =   covK( 3 , 3 )   ; //   z0  , z0

  cov[10] = - covK( 4 , 0 )   ; //   tanl, d0
  cov[11] =   covK( 4 , 1 )   ; //   tanl, phi
  cov[12] =   covK( 4 , 2 ) * alpha  ; //   tanl, omega
  cov[13] =   covK( 4 , 3 )   ; //   tanl, z0
  cov[14] =   covK( 4 , 4 )   ; //   tanl, tanl

 
  trk->setCovMatrix( cov ) ;


  float pivot[3] ;

//  SJA:the code below looks suspisious as I am not sure that this is the point of closest approach to the IP rathet it seams to be the position of the last measurement, in this case the pseudo measurement inside the beampipe ILDIPMeasL
//  pivot[0] =  ((TKalTrackSite&) cursite).GetPivot()(0) ;
//  pivot[1] =  ((TKalTrackSite&) cursite).GetPivot()(1) ;
//  pivot[2] =  ((TKalTrackSite&) cursite).GetPivot()(2) ;

// the reference point is in fact 0,0,0 as this was used for helix.MoveTo
  pivot[0] =  0.0 ;
  pivot[1] =  0.0 ;
  pivot[2] =  0.0 ;

  trk->setReferencePoint( pivot ) ;

  streamlog_out( DEBUG4 ) << " kaltest track parameters: "
			 << " chi2/ndf " << chi2 / ndf  
    			 << " chi2 " <<  chi2 << std::endl 
    
    			 << "\t D0 "          <<  d0         <<  "[+/-" << sqrt( cov[0] ) << "] " 
			 << "\t Phi :"        <<  phi        <<  "[+/-" << sqrt( cov[2] ) << "] " 
			 << "\t Omega "       <<  omega      <<  "[+/-" << sqrt( cov[5] ) << "] " 
			 << "\t Z0 "          <<  z0         <<  "[+/-" << sqrt( cov[9] ) << "] " 
			 << "\t tan(Lambda) " <<  tanLambda  <<  "[+/-" << sqrt( cov[14]) << "] " 
    
			 << "\t pivot : [" << pivot[0] << ", " << pivot[1] << ", "  << pivot[2] 
			 << " - r: " << std::sqrt( pivot[0]*pivot[0]+pivot[1]*pivot[1] ) << "]" 
			 << std::endl ;
  

  return trk;

}

