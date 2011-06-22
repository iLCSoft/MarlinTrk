#ifndef MarlinKalTestTrack_h
#define MarlinKalTestTrack_h

#include "IMarlinTrack.h"

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

  //  MarlinKalTestTrack( EVENT::Track* trk, MarlinKalTest* ktest) ;
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
  bool fit( Bool_t fitDirection );
  
 // returns a pointer to a New LCIO Track with fit parameters determinded at the IP. The responsiblitiy for deletion lies with the caller.
  IMPL::TrackImpl* getIPFit();

  // returns a pointer to an LCIO Track whose referece point is the closest to the specified point.
  virtual IMPL::TrackImpl* getNearestFitToPoint(float* point) {  MarlinTrk::IMarlinTrackException exp; throw; } ; // Not Implemented Yet !
  
  // returns a pointer to an LCIO Track whose referece point is the closest to a cylinder of radius r which is centered at the origin parallel to the z axis.
  virtual IMPL::TrackImpl* getNearestFitToCylinder(float r) {  MarlinTrk::IMarlinTrackException exp; throw; } ; // Not Implemented Yet !
  
  // returns a pointer to an LCIO Track whose referece point is the closest to a plane normal to the z axis.
  virtual IMPL::TrackImpl* getNearestFitToZPlane(float z) {  MarlinTrk::IMarlinTrackException exp; throw; } ; // Not Implemented Yet !

   // add hit to track
  void addHit(EVENT::TrackerHit* hit)  ;  

  // memeber variables 
  EVENT::Track* _initialLCTrack ;

  TKalTrack* _kaltrack;


  EVENT::TrackerHitVec _lcioHits ; 

  TObjArray* _kalhits;

  MarlinKalTest* _ktest;

} ;

#endif
