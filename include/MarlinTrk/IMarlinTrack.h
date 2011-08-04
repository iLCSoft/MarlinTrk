#ifndef IMarlinTrack_h
#define IMarlinTrack_h

#include "lcio.h"

#include "EVENT/TrackerHit.h"
#include "IMPL/TrackStateImpl.h"

#include "gearimpl/Vector3D.h"


#include <exception>

namespace MarlinTrk{
    
  
  class IMarlinTrack {
    
  public:
    
    virtual ~IMarlinTrack() {};

    // add hit to track
    virtual void addHit(EVENT::TrackerHit* hit) = 0 ;

    // perform the fit of all current hits, return code via int
    virtual int fit( bool fitDirection ) = 0 ;


//    // propagate the fit to the point of closest approach to the nominal IP=(0.0,0.0,0.0). 
//    virtual int propagateToIP( IMPL::TrackStateImpl& ts) = 0 ;

    // propagate the fit to the point of closest approach to the given point. 
    virtual int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts) = 0 ;
    
    // propagate to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
    virtual int propagateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, int& layerNumber) = 0 ;

    // propagate to numbered sensitive layer, returning TrackState via provided reference 
    virtual int propagateToLayer( bool direction, int layerNumber, IMPL::TrackStateImpl& ts) = 0 ;


    // extrapolate the fit to the point of closest approach to the given point. 
    virtual int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts) = 0 ;
  
    // extrapolate to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
    virtual int extrapolateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, int& layerNumber) = 0 ;

    // extrapolate to numbered sensitive layer, returning TrackState via provided reference 
    virtual int extrapolateToLayer( bool direction, int layerNumber, IMPL::TrackStateImpl& ts) = 0 ;


    // extrapolate to next sensitive layer, returning intersection point in global coordinates, layer number of sensitive layer returned via layerNumber reference 
    virtual int intersectionWithNextLayer( bool direction, int& layerNumber, gear::Vector3D& point) = 0 ;
    
    // extrapolate to numbered sensitive layer, returning intersection point in global coordinates 
    virtual int intersectionWithLayer( int layerNumber, bool direction, gear::Vector3D& point) = 0 ;


       
    
  protected:
    
  private:
    
    IMarlinTrack& operator=( const IMarlinTrack&) ; // disallow assignment operator 
    
  } ;
  
} // end of MarlinTrk namespace 

#endif

