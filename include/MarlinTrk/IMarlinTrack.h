#ifndef IMarlinTrack_h
#define IMarlinTrack_h

#include <cfloat>

#include "lcio.h"

#include "EVENT/TrackerHit.h"
#include "IMPL/TrackStateImpl.h"

#include "gearimpl/Vector3D.h"


#include <exception>

namespace MarlinTrk{
  
  
  /** Interface for generic tracks in MarlinTrk. The interface should provide the functionality to
   *  perform track finding and fitting. It is asssumed that the underlying implemetation will by 
   *  a Kalman Filter or a similar algorithm.
   *
   * @version $Id$
   * @author S.Aplin, F. Gaede DESY
   */
  
  class IMarlinTrack {
    
  public:
    
    /** boolean constant for defining backward direction - to be used for intitialize, fit, propagate etc. */
    static const bool backward = false ;
    
    /** boolean constant for defining backward direction - to be used for intitialize, fit, propagate etc. */
    static const bool forward  = ! backward  ;
    
    /**default d'tor*/
    virtual ~IMarlinTrack() {};
    
    /** add hit to track - hits have to be given ordered with in time (typically outgoing)
     */
    virtual int addHit(EVENT::TrackerHit* hit) = 0 ;

    /** initialise the fit using the supplied hits only - use IMarlinTrack::backward or IMarlinTrack::forward to specify the direction
     *  with respect to the order used in addHit(EVENT::TrackerHit* hit).
     */
    virtual int initialise( bool direction ) = 0 ; 

    /** same as initialise( bool direction )*/
    inline int  initialize( bool direction ) { return initialise( direction ) ; } ; 
    
    /** initialise the fit with a track state, and z component of the B field in Tesla.
     *  The default for initalise_at_end is set to true as it is expected that the most common case will be to fit backwards in time.
     *  If it is desired to fit in the opposite direction then initalise_at_end should be set to false and the intialisation will be done at the first hit. 
     */
    virtual int initialise( const IMPL::TrackStateImpl& ts, double bfield_z, bool initalise_at_end = true ) = 0 ;

    /** same as initialise(const IMPL::TrackStateImpl& ts, double bfield_z, bool initalise_at_end = true ) */
    inline int  initialize( const IMPL::TrackStateImpl& ts, double bfield_z, bool initalise_at_end = true ) { 

      return initialise( ts, bfield_z, initalise_at_end ) ;
    } ; 

    /** perform the fit of all current hits, return code via int - use IMarlinTrack::backward or IMarlinTrack::forward to specify the direction
     *  with respect to the order used in addHit(EVENT::TrackerHit* hit)
     */
    virtual int fit( bool fitDirection ) = 0 ;
    
    /** update the current fit using the supplied hit, return code via int. Provides the Chi2 increment to the fit from adding the hit via reference. 
     */
    virtual int addAndFit( EVENT::TrackerHit* hit, double& chi2increment, double maxChi2Increment=DBL_MAX ) = 0 ;
    
    
    // Track State Accessesors
    
    /** get track state, returning TrackState, chi2 and ndf via reference 
     */
    virtual int getTrackState( IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

    /** get track state at measurement associated with the given hit, returning TrackState, chi2 and ndf via reference 
     */
    virtual int getTrackState( EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;


    // PROPAGATORS 
    
    /** propagate the fit to the point of closest approach to the given point, returning TrackState, chi2 and ndf via reference    
     */
    virtual int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

    /** propagate track state at measurement associated with the given hit, the fit to the point of closest approach to the given point,
     *  returning TrackState, chi2 and ndf via reference   
     */
    virtual int propagate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
    
    /** propagate the fit to next sensitive layer, returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference 
     */
    virtual int propagateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& layerID ) = 0 ;

    /** propagate track state at measurement associated with the given hit, to next sensitive layer, 
     *	returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference
     */
    virtual int propagateToNextLayer( bool direction, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& layerID ) = 0 ;
    
    /** propagate fit to numbered sensitive layer, returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference 
     */
    virtual int propagateToLayer( bool direction, int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) = 0 ;

    /** propagate track state at measurement associated with the given hit, to numbered sensitive layer, 
     *  returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference 
     */
    virtual int propagateToLayer( bool direction, int layerID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) = 0 ;


    // EXTRAPOLATORS

    /** extrapolate the fit to the point of closest approach to the given point, returning TrackState, chi2 and ndf via reference   
     */
    virtual int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

    /** extrapolate track state at measurement associated with the given hit, to the point of closest approach to the given point, 
     *	returning TrackState, chi2 and ndf via reference   
     */
    virtual int extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

    /** extrapolate the fit to next sensitive layer, returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference 
     */
    virtual int extrapolateToNextLayer( bool direction, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) = 0 ;

    /** extrapolate track state at measurement associated with the given hit, to next sensitive layer, 
     *	returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference 
     */
    virtual int extrapolateToNextLayer( bool direction, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) = 0 ;

    /** extrapolate the fit to numbered sensitive layer, returning TrackState via provided reference 
     */
    virtual int extrapolateToLayer( bool direction, int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) = 0 ;

    /** extrapolate track state at measurement associated with the given hit, to numbered sensitive layer, 
     *  returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference 
     */
    virtual int extrapolateToLayer( bool direction, int layerID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) = 0 ;


    // INTERSECTORS

    /** extrapolate the fit to next sensitive layer, returning intersection point in global 
     *  coordinates and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int intersectionWithNextLayer( bool direction, gear::Vector3D& point, int& detElementID ) = 0 ;

    /** extrapolate track state at measurement associated with the given hit, to next sensitive layer, 
     *  returning intersection point in global coordinates and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int intersectionWithNextLayer( bool direction, EVENT::TrackerHit* hit, gear::Vector3D& point, int& detElementID ) = 0 ;
    
    /** extrapolate the fit to numbered sensitive layer, returning intersection point in global coordinates and integer ID of the 
     *  intersected sensitive detector element via reference 
     */
    virtual int intersectionWithLayer( bool direction, int layerID, gear::Vector3D& point, int& detElementID ) = 0 ;
    
    /** extrapolate track state at measurement associated with the given hit, to numbered sensitive layer,
     *  returning intersection point in global coordinates and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int intersectionWithLayer( bool direction, int layerID, EVENT::TrackerHit* hit, gear::Vector3D& point, int& detElementID ) = 0 ;
    
    
  protected:
    
  private:
    
    IMarlinTrack& operator=( const IMarlinTrack&) ; // disallow assignment operator 
    
  } ;
  
} // end of MarlinTrk namespace 

#endif

