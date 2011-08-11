#ifndef MarlinKalTestTrack_h
#define MarlinKalTestTrack_h

#include "IMarlinTrack.h"
#include "IMarlinTrkSystem.h"

#include <TObjArray.h>

#include <cmath>

#include "TMatrixD.h"


class TKalTrack ;
class THelicalTrack ;

class MarlinKalTest;


namespace EVENT{
  class TrackerHit ;
}


class MarlinKalTestTrack : public MarlinTrk::IMarlinTrack {

 public:

  MarlinKalTestTrack(MarlinKalTest* ktest) ;

  ~MarlinKalTestTrack() ;

 protected:
  
 private:

  MarlinKalTestTrack(const MarlinKalTestTrack&) ;                 // Prevent copy-construction
  MarlinKalTestTrack& operator=(const MarlinKalTestTrack&) ;      // Prevent assignment

  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange( double phi){
    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
    return phi ;
  }

  // make member functions private to force use through interface


  // add hit to track
   void addHit(EVENT::TrackerHit* hit) ;
  
  // perform the fit of all current hits, return code via int
   int fit( bool fitDirection ) ;
  
  
  // propagate the fit to the point of closest approach to the nominal IP=(0.0,0.0,0.0). 
   //   int propagateToIP( IMPL::TrackStateImpl& ts) ;
  
  // propagate the fit to the point of closest approach to the given point. 
   int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts) ;
  
    // propagate to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
   int propagateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, int& layerNumber) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  // propagate to numbered sensitive layer, returning TrackState via provided reference 
   int propagateToLayer( bool direction, int layerNumber, IMPL::TrackStateImpl& ts) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  
  // extrapolate the fit to the point of closest approach to the given point. 
   int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts) ;
  
  // extrapolate to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
   int extrapolateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, int& layerNumber) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  // extrapolate to numbered sensitive layer, returning TrackState via provided reference 
   int extrapolateToLayer( bool direction, int layerNumber, IMPL::TrackStateImpl& ts) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  
  // extrapolate to next sensitive layer, returning intersection point in global coordinates, layer number of sensitive layer returned via layerNumber reference 
   int intersectionWithNextLayer( bool direction, int& layerNumber, gear::Vector3D& point) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  // extrapolate to numbered sensitive layer, returning intersection point in global coordinates 
   int intersectionWithLayer( bool direction, int layerNumber, gear::Vector3D& point)  ;

   // fill LCIO Track State with parameters from TKalTrack
   void ToLCIOTrackState( IMPL::TrackStateImpl& ts ) ;
  
   // fill LCIO Track State with parameters from helix and cov matrix 
   void ToLCIOTrackState( const THelicalTrack& helix, const TMatrixD& cov, IMPL::TrackStateImpl& ts ) ;


  // memeber variables 

  TKalTrack* _kaltrack;


  EVENT::TrackerHitVec _lcioHits ; 

  TObjArray* _kalhits;

  MarlinKalTest* _ktest;

} ;

#endif
