#ifndef IMarlinTrack_h
#define IMarlinTrack_h

#include "lcio.h"
#include "IMPL/TrackImpl.h"
#include "EVENT/TrackerHit.h"

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
    
    virtual bool fit( bool fitDirection ) = 0 ;
    
    // returns a pointer to a New LCIO Track with fit parameters determinded at the IP. The responsiblitiy for deletion lies with the caller.
    virtual IMPL::TrackImpl* getIPFit() = 0 ;
    
    // returns a pointer to an LCIO Track whose referece point is the closest to the specified point.
    virtual IMPL::TrackImpl* getNearestFitToPoint(float* point) = 0 ;
    
    // returns a pointer to an LCIO Track whose referece point is the closest to a cylinder of radius r which is centered at the origin parallel to the z axis.
    virtual IMPL::TrackImpl* getNearestFitToCylinder(float r) = 0 ;
    
    // returns a pointer to an LCIO Track whose referece point is the closest to a plane normal to the z axis.
    virtual IMPL::TrackImpl* getNearestFitToZPlane(float z) = 0 ;
    
    // add hit to track
    virtual void addHit(EVENT::TrackerHit* hit) = 0 ;
    
  protected:
    
  private:
    
    IMarlinTrack& operator=( const IMarlinTrack&) ; // disallow assignment operator 
    
  } ;
  
} // end of MarlinTrk namespace 

#endif

