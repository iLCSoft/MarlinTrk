
#include "MarlinTrk/MarlinKalTestTrack.h"

#include "MarlinTrk/MarlinKalTest.h"

#include <kaltest/TKalDetCradle.h>
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

#include "gear/GEAR.h"
#include "gear/BField.h"

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
 
  //SJA:FIXME: what happens if the fit is already performed and then we add a hit which is in-between the hits already fitted. 
 
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

      if( (surf = dynamic_cast<TVSurface*> (  measlayers[i] )) )   {
      	// surf = dynamic_cast<TVSurface*> (  measlayers[i] )
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




int MarlinKalTestTrack::fit( bool fitDirection ) {

  //SJA:FIXME: do we need to sort the hits here ... ?

  const bool gkDir = fitDirection ; 

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
    pDummyHit = (new ILDCylinderHit(*static_cast<ILDCylinderHit*>( startingHit )));
  }
  else if ( (pDummyHit = dynamic_cast<ILDPlanarHit *>( startingHit )) ) {
    pDummyHit = (new ILDPlanarHit(*static_cast<ILDPlanarHit*>( startingHit )));
  }
  else {
    streamlog_out( ERROR) << "<<<<<< KalTrack::fitTrack(): dynamic_cast failed for hit type >>>>>>>" << std::endl;
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

  TIter next(_kalhits, gkDir); // fit inwards to IP, if gkDir = kIterBackward

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
  
  return 0;

}



int MarlinKalTestTrack::propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts ){

  streamlog_out(DEBUG4) << "MarlinKalTestTrack::Propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts ) called " << std::endl ;

  // get the current site. SJA:FIXME: should we check here if this is valid?
  // it would be better to get the site closest to the point in s ...
  // here we assume that we want to take the last filtered site
  const TVKalSite& cursite = _kaltrack->GetCurSite();

  TKalTrackState& trkState = (TKalTrackState&) cursite.GetCurState(); // this segfaults if no hits are present

  Int_t    ndf  = _kaltrack->GetNDF();
  Double_t chi2 = _kaltrack->GetChi2();
  
  
  THelicalTrack helix = trkState.GetHelix() ;
  double dPhi ;


  // convert the gear point supplied to TVector3
  const TVector3 tpoint( point.x(), point.y(), point.z() ) ;

  // the last layer crossed by the track before point 
  const ILDVMeasLayer* ml = _ktest->getLastMeasLayer(helix, tpoint);
  
  streamlog_out( DEBUG4 ) << "  MarlinKalTestTrack - last surface before point = " << ml->GetMLName() << std::endl;

  Int_t sdim = trkState.GetDimension();  // dimensions of the track state, it will be 5 or 6
  TKalMatrix sv(sdim,1);

  TKalMatrix  F(sdim,sdim);              // propagator matrix to be returned by transport function
  F.UnitMatrix();                        // set the propagator matrix to the unit matrix

  TKalMatrix  Q(sdim,sdim);                     // noise matrix to be returned by transport function 
  Q.Zero();        
  TVector3    x0;                        // intersection point to be returned by transport

  const TVMeasLayer& tvml   = dynamic_cast<const TVMeasLayer&>(*ml) ;         // cast from ILDVMeasurementLayer
  const TKalTrackSite& site = dynamic_cast<const TKalTrackSite&>(cursite) ;   // cast from TVKalSite
  
  int return_code = _ktest->_det->Transport(site, tvml, x0, sv, F, Q ) ;      // transport to last layer cross before point 
  
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
  c0 = DF * c0 + Qms * DFt ;                 // update the covariance matrix 


  //============== convert parameters to LCIO convention ====

  // fill 5x5 covariance matrix from the 6x6 covariance matrix return by trkState.GetCovMat()  above
  TMatrixD covK(5,5) ;  for(int i=0;i<5;++i) for(int j=0;j<5;++j) covK[i][j] = c0[i][j] ;

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
    
  //  Double_t cpa  = trkState(2, 0);
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

 
  ts.setCovMatrix( cov ) ;


  float pivot[3] ;

  pivot[0] =  point.x() ;
  pivot[1] =  point.y() ;
  pivot[2] =  point.z() ;

  ts.setReferencePoint( pivot ) ;

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
  

  return 0;

}

