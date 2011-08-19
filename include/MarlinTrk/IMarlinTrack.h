#ifndef IMarlinTrack_h
#define IMarlinTrack_h

#include <cfloat>

#include "lcio.h"

#include "EVENT/TrackerHit.h"
#include "IMPL/TrackStateImpl.h"

#include "gearimpl/Vector3D.h"


#include <exception>

namespace MarlinTrk{
    
  
  class IMarlinTrack {
    
  public:
    
    virtual ~IMarlinTrack() {};

    //** add hit to track
    virtual int addHit(EVENT::TrackerHit* hit) = 0 ;

    //** initialise the fit using the supplied hits only, using the given order to determine the direction of the track
    // SJA::FIXME: replace bool with type specifying the order. For now direction determines if the the hits were added odered in time of in reverse time.
    virtual int initialise( bool direction ) = 0 ; 

    //** initialise the fit with a track state, and z component of the B field in Tesla. The default for initalise_at_end is set to true as it is expected that the most common case will be to fit backwards in time. If it is desired to fit in the opposite direction then initalise_at_end should be set to false and the intialisation will be done at the first hit. 
    virtual int initialise( const IMPL::TrackStateImpl& ts, double bfield_z, bool initalise_at_end = true ) = 0 ;

    //** perform the fit of all current hits, return code via int
    virtual int fit( bool fitDirection ) = 0 ;

    //** update the current fit using the supplied hit, return code via int. Provides the Chi2 increment to the fit from adding the hit via reference. 
    virtual int addAndFit( EVENT::TrackerHit* hit, double& chi2increment, double maxChi2Increment=DBL_MAX ) = 0 ;


    //** get track state, return code via int
    virtual int getTrackState( IMPL::TrackStateImpl& ts ) = 0 ;

    //** get track state at measurement associated with the given hit, return code via int
    virtual int getTrackState( EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts ) = 0 ;

    
    //** propagate the fit to the point of closest approach to the given point. 
    virtual int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts) = 0 ;

    //** propagate track state at measurement associated with the given hit, the fit to the point of closest approach to the given point. 
    virtual int propagate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts) = 0 ;
    
    //** propagate the fit to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
    virtual int propagateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, int& layerNumber) = 0 ;

    //** propagate track state at measurement associated with the given hit, to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
    virtual int propagateToNextLayer( bool direction, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, int& layerNumber) = 0 ;

    //** propagate fit to numbered sensitive layer, returning TrackState via provided reference 
    virtual int propagateToLayer( bool direction, int layerNumber, IMPL::TrackStateImpl& ts) = 0 ;

    //** propagate track state at measurement associated with the given hit, to numbered sensitive layer, returning TrackState via provided reference 
    virtual int propagateToLayer( bool direction, int layerNumber, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts) = 0 ;


    //** extrapolate the fit to the point of closest approach to the given point. 
    virtual int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts) = 0 ;

    //** extrapolate track state at measurement associated with the given hit, to the point of closest approach to the given point. 
    virtual int extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts) = 0 ;

    //** extrapolate the fit to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
    virtual int extrapolateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, int& layerNumber) = 0 ;

    //** extrapolate track state at measurement associated with the given hit, to next sensitive layer, returning TrackState via provided reference, layer number of sensitive layer returned via layerNumber reference 
    virtual int extrapolateToNextLayer( bool direction, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, int& layerNumber) = 0 ;

    //** extrapolate the fit to numbered sensitive layer, returning TrackState via provided reference 
    virtual int extrapolateToLayer( bool direction, int layerNumber, IMPL::TrackStateImpl& ts) = 0 ;

    //** extrapolate track state at measurement associated with the given hit, to numbered sensitive layer, returning TrackState via provided reference 
    virtual int extrapolateToLayer( bool direction, int layerNumber, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts) = 0 ;


    //** extrapolate the fit to next sensitive layer, returning intersection point in global coordinates, layer number of sensitive layer returned via layerNumber reference 
    virtual int intersectionWithNextLayer( bool direction, int& layerNumber, gear::Vector3D& point) = 0 ;

    //** extrapolate track state at measurement associated with the given hit, to next sensitive layer, returning intersection point in global coordinates, layer number of sensitive layer returned via layerNumber reference 
    virtual int intersectionWithNextLayer( bool direction, EVENT::TrackerHit* hit, int& layerNumber, gear::Vector3D& point) = 0 ;
    
    //** extrapolate the fit to numbered sensitive layer, returning intersection point in global coordinates 
    virtual int intersectionWithLayer( bool direction, int layerNumber, gear::Vector3D& point) = 0 ;

    //** extrapolate track state at measurement associated with the given hit, to numbered sensitive layer, returning intersection point in global coordinates 
    virtual int intersectionWithLayer( bool direction, int layerNumber, EVENT::TrackerHit* hit, gear::Vector3D& point) = 0 ;

    
  protected:
    
  private:
    
    IMarlinTrack& operator=( const IMarlinTrack&) ; // disallow assignment operator 
    
  } ;
  
} // end of MarlinTrk namespace 

#endif

