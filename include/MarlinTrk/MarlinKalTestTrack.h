#ifndef MarlinKalTestTrack_h
#define MarlinKalTestTrack_h

#include "IMarlinTrack.h"
#include "IMarlinTrkSystem.h"

#include <TObjArray.h>

#include <cmath>

class TKalTrack ;
class TKalTrackState ;

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
  
  // perform the fit of all current hits, return true is fit succeeds and false if it fails
  bool fit( bool fitDirection ) ;
  
  // Propagate the fit to the point of closest approach to the given point. The responsiblitiy for deletion lies with the caller.
  IMPL::TrackStateImpl* propagate(gear::Vector3D* point) {  throw MarlinTrk::Exception("Function Not Implemented Yet "); } ; 

  // Propagate the fit to the point of closest approach to the nominal IP=(0.0,0.,0.0). The responsiblitiy for deletion lies with the caller.
  IMPL::TrackStateImpl* propagateToIP() ;
  
  // extrapolate to next sensitive layer, returning intersection point in global coordinates, layer number of sensitive layer returned via layerNumber reference 
  gear::Vector3D intersectionWithNextSensitiveLayer( bool direction, int& layerNumber) { throw MarlinTrk::Exception("Function Not Implemented Yet "); } ;
  
  // extrapolate to numbered sensitive layer, returning intersection point in global coordinates 
  gear::Vector3D intersectionWithSensitiveLayer( int layerNumber) { throw MarlinTrk::Exception("Function Not Implemented Yet"); } ; 
  

  // memeber variables 

  TKalTrack* _kaltrack;


  EVENT::TrackerHitVec _lcioHits ; 

  TObjArray* _kalhits;

  MarlinKalTest* _ktest;

} ;

#endif
