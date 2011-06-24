#ifndef IMarlinTrack_h
#define IMarlinTrack_h

#include "lcio.h"

#include "EVENT/TrackerHit.h"
#include "IMPL/TrackStateImpl.h"

#include "gearimpl/Vector3D.h"


#include <exception>

namespace MarlinTrk{
  
  class IMarlinTrackException: public std::exception
    {
      virtual const char* what() const throw()
      {
	return "IMarlinTrackException occurred";
      }
    } ;
  
  
  class IMarlinTrack {
    
  public:
    
    virtual ~IMarlinTrack() {};

    // add hit to track
    virtual void addHit(EVENT::TrackerHit* hit) = 0 ;

    // perform the fit of all current hits, return true is fit succeeds and false if it fails
    virtual bool fit( bool fitDirection ) = 0 ;

    // Propagate the fit to the point of closest approach to the given point. The responsiblitiy for deletion lies with the caller.
    virtual IMPL::TrackStateImpl* propagate(gear::Vector3D* point) = 0 ;
    
    // Propagate the fit to the point of closest approach to the nominal IP=(0.0,0.0,0.0). The responsiblitiy for deletion lies with the caller.
    virtual IMPL::TrackStateImpl* propagateToIP() = 0 ;
  
    // extrapolate to next sensitive layer, returning intersection point in global coordinates, layer number of sensitive layer returned via layerNumber reference 
    virtual gear::Vector3D intersectionWithNextSensitiveLayer( bool direction, int& layerNumber) = 0 ;

    // extrapolate to numbered sensitive layer, returning intersection point in global coordinates 
    virtual gear::Vector3D intersectionWithSensitiveLayer( int layerNumber) = 0 ;

       
    
  protected:
    
  private:
    
    IMarlinTrack& operator=( const IMarlinTrack&) ; // disallow assignment operator 
    
  } ;
  
} // end of MarlinTrk namespace 

#endif

