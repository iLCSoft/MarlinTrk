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

  // make member functions private to force use through interface
  
  //** add hit to track
  void addHit(EVENT::TrackerHit* hit) ;

  //** initialise the fit using the supplied hits only, using the given order to determine the direction of the track
  // SJA::FIXME: replace bool with type specifying the order. For now direction determines if the the hits were added odered in time of in reverse time.
  int initialise( bool direction ) ; 

  //** initialise the fit with a track state 
  int initialise( const IMPL::TrackStateImpl& ts) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  //** perform the fit of all current hits, return code via int
  int fit( bool fitDirection ) ;  
  
  //** update the current fit using the supplied hit, return code via int. Provides the Chi2 increment to the fit from adding the hit via reference. 
  int addAndFit( EVENT::TrackerHit* hit, double& chi2increment, double maxChi2Increment=DBL_MAX ) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); };


  //** get track state, return code via int
  int getTrackState( IMPL::TrackStateImpl& ts ) ;
  
  //** get track state at measurement associated with the given hit, return code via int
  int getTrackState( EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts ) ;
  
  
  //** propagate the fit to the point of closest approach to the given point. 
  int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts) ;
  
  //** propagate track state at measurement associated with the given hit, the fit to the point of closest approach to the given point. 
  int propagate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts)  {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;

  //** propagate to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
  int propagateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, int& layerNumber) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;

  //** propagate track state at measurement associated with the given hit, to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
  int propagateToNextLayer( bool direction, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, int& layerNumber) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  //** propagate to numbered sensitive layer, returning TrackState via provided reference 
  int propagateToLayer( bool direction, int layerNumber, IMPL::TrackStateImpl& ts) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;

  //** propagate track state at measurement associated with the given hit, to numbered sensitive layer, returning TrackState via provided reference 
  int propagateToLayer( bool direction, int layerNumber, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ; 
 
 
  //** extrapolate the fit to the point of closest approach to the given point. 
  int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts) ;
  
  //** extrapolate track state at measurement associated with the given hit, to the point of closest approach to the given point. 
  int extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
    
  //** extrapolate to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
  int extrapolateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, int& layerNumber) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;

  //** extrapolate track state at measurement associated with the given hit, to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
  int extrapolateToNextLayer( bool direction, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, int& layerNumber) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  //** extrapolate to numbered sensitive layer, returning TrackState via provided reference 
  int extrapolateToLayer( bool direction, int layerNumber, IMPL::TrackStateImpl& ts) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  //** extrapolate track state at measurement associated with the given hit, to numbered sensitive layer, returning TrackState via provided reference 
  int extrapolateToLayer( bool direction, int layerNumber, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts)  {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  

  //** extrapolate to next sensitive layer, returning intersection point in global coordinates, layer number of sensitive layer returned via layerNumber reference 
  int intersectionWithNextLayer( bool direction, int& layerNumber, gear::Vector3D& point) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  //** extrapolate track state at measurement associated with the given hit, to next sensitive layer, returning intersection point in global coordinates, layer number of sensitive layer returned via layerNumber reference 
  int intersectionWithNextLayer( bool direction, EVENT::TrackerHit* hit, int& layerNumber, gear::Vector3D& point) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;

  //** extrapolate to numbered sensitive layer, returning intersection point in global coordinates 
  int intersectionWithLayer( bool direction, int layerNumber, gear::Vector3D& point)  ;

  //** extrapolate track state at measurement associated with the given hit, to numbered sensitive layer, returning intersection point in global coordinates 
  int intersectionWithLayer( bool direction, int layerNumber, EVENT::TrackerHit* hit, gear::Vector3D& point) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;

  //** end of memeber functions from IMarlinTrack interface

  
  //** fill LCIO Track State with parameters from TKalTrack
  void ToLCIOTrackState( IMPL::TrackStateImpl& ts ) ;
  
  //** fill LCIO Track State with parameters from helix and cov matrix 
  void ToLCIOTrackState( const THelicalTrack& helix, const TMatrixD& cov, IMPL::TrackStateImpl& ts ) ;


  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange( double phi){
    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
    return phi ;
  }
  
  
  // memeber variables 
  
  TKalTrack* _kaltrack;
  
  EVENT::TrackerHitVec _lcioHits ; 

  TObjArray* _kalhits;

  MarlinKalTest* _ktest;

  //** used to store whether initial track state has been supplied or created 
  bool _initialised ;

  //** used to store whether smoothing has been performed
  bool _smoothed ;
  
} ;

#endif
