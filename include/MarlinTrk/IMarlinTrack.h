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
    
    /** boolean constant for defining backward direction - to be used for intitialise */
    static const bool backward = false ;
    
    /** boolean constant for defining backward direction - to be used for intitialise */
    static const bool forward  = ! backward  ;
    
    
    
    static  const int modeBackward = - 1 ;
    static  const int modeClosest  =   0 ;
    static  const int modeForward  = + 1 ;
    
    
    static const int success  = 0 ;  // no error
    static const int error = 1 ;
    static const int bad_intputs = 3 ;
    static const int no_intersection = 4 ; // no intersection found
    static const int site_discarded = 5 ;  // measurement discarded by the fitter
    static const int site_fails_chi2_cut = 6 ;  // measurement discarded by the fitter due to chi2 cut
    
    
    
    
    /**default d'tor*/
    virtual ~IMarlinTrack() {};
    
    /** add hit to track - the hits have to be added ordered in time ( i.e. typically outgoing )
     *  this order will define the direction of the energy loss used in the fit
     */
    virtual int addHit(EVENT::TrackerHit* hit) = 0 ;
    
    /** initialise the fit using the hits added up to this point -
     *  the fit direction has to be specified using IMarlinTrack::backward or IMarlinTrack::forward. 
     *  this is the order  wrt the order used in addHit() that will be used in the fit() 
     */
    virtual int initialise( bool fitDirection ) = 0 ; 
    
    /** initialise the fit with a track state, and z component of the B field in Tesla.
     *  the fit direction has to be specified using IMarlinTrack::backward or IMarlinTrack::forward. 
     *  this is the order that will be used in the fit().
     *  it is the users responsibility that the track state is consistent with the order
     *  of the hits used in addHit() ( i.e. the direction of energy loss )
     */
    virtual int initialise(  const EVENT::TrackState& ts, double bfield_z, bool fitDirection ) = 0 ;
    
    
    /** perform the fit of all current hits, returns error code ( IMarlinTrack::success if no error ) .
     *  the fit will be performed  in the order specified at initialise() wrt the order used in addHit(), i.e.
     *  IMarlinTrack::backward implies fitting from the outside to the inside for tracks comming from the IP.
     */
    virtual int fit() = 0 ;
    
    
    /** update the current fit using the supplied hit, return code via int. Provides the Chi2 increment to the fit from adding the hit via reference. 
     *  the given hit will not be added if chi2increment > maxChi2Increment. 
     */
    virtual int addAndFit( EVENT::TrackerHit* hit, double& chi2increment, double maxChi2Increment=DBL_MAX ) = 0 ;

    
    /** obtain the chi2 increment which would result in adding the hit to the fit. This method will not alter the current fit, and the hit will not be stored in the list of hits or outliers
     */
    virtual int testChi2Increment( EVENT::TrackerHit* hit, double& chi2increment ) = 0 ;

    
    /** smooth all track states 
     */
    virtual int smooth() = 0 ;
    
    
    /** smooth track states from the last filtered hit back to the measurement site associated with the given hit 
     */
    virtual int smooth( EVENT::TrackerHit* hit ) = 0 ;
    
    
    
    // Track State Accessesors
    
    /** get track state, returning TrackState, chi2 and ndf via reference 
     */
    virtual int getTrackState( IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
    
    
    /** get track state at measurement associated with the given hit, returning TrackState, chi2 and ndf via reference 
     */
    virtual int getTrackState( EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
    
    /** get the list of hits included in the fit, together with the chi2 contributions of the hits. 
     *  Pointers to the hits together with their chi2 contribution will be filled into a vector of 
     *  pairs consitining of the pointer as the first part of the pair and the chi2 contribution as
     *  the second.
     */
    virtual int getHitsInFit( std::vector<std::pair<EVENT::TrackerHit*, double> >& hits ) = 0 ;

    /** get the list of hits which have been rejected by from the fit due to the a chi2 increment greater than threshold,
     *  Pointers to the hits together with their chi2 contribution will be filled into a vector of 
     *  pairs consitining of the pointer as the first part of the pair and the chi2 contribution as
     *  the second.
     */
    virtual int getOutliers( std::vector<std::pair<EVENT::TrackerHit*, double> >& hits ) = 0 ;
    
    // PROPAGATORS 
    
    /** propagate the fit to the point of closest approach to the given point, returning TrackState, chi2 and ndf via reference    
     */
    virtual int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
    
    
    /** propagate the fit at the measurement site associated with the given hit, to the point of closest approach to the given point,
     *  returning TrackState, chi2 and ndf via reference   
     */
    virtual int propagate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
    
    
    /** propagate fit to numbered sensitive layer, returning TrackState, chi2, ndf and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int propagateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** propagate the fit at the measurement site associated with the given hit, to numbered sensitive layer, 
     *  returning TrackState, chi2, ndf and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int propagateToLayer( int layerID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** propagate the fit to sensitive detector element, returning TrackState, chi2 and ndf via reference
     */
    virtual int propagateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
    
    /** propagate the fit at the measurement site associated with the given hit, to sensitive detector element, 
     *  returning TrackState, chi2 and ndf via reference 
     */
    virtual int propagateToDetElement( int detEementID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
    
    
    
    // EXTRAPOLATORS
    
    /** extrapolate the fit to the point of closest approach to the given point, returning TrackState, chi2 and ndf via reference   
     */
    virtual int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to the point of closest approach to the given point, 
     *  returning TrackState, chi2 and ndf via reference   
     */
    virtual int extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
    
    /** extrapolate the fit to numbered sensitive layer, returning TrackState, chi2, ndf and integer ID of the intersected sensitive detector element via reference
     */
    virtual int extrapolateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to numbered sensitive layer, 
     *  returning TrackState, chi2, ndf and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int extrapolateToLayer( int layerID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit to sensitive detector element, returning TrackState, chi2 and ndf via reference
     */
    virtual int extrapolateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to sensitive detector element, 
     *  returning TrackState, chi2 and ndf via reference 
     */
    virtual int extrapolateToDetElement( int detEementID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
    
    
    // INTERSECTORS
    
    /** extrapolate the fit to numbered sensitive layer, returning intersection point in global coordinates and integer ID of the 
     *  intersected sensitive detector element via reference 
     */
    virtual int intersectionWithLayer( int layerID, gear::Vector3D& point, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to numbered sensitive layer,
     *  returning intersection point in global coordinates and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int intersectionWithLayer( int layerID, EVENT::TrackerHit* hit, gear::Vector3D& point, int& detElementID, int mode=modeClosest ) = 0  ;
    
    
    /** extrapolate the fit to numbered sensitive detector element, returning intersection point in global coordinates via reference 
     */
    virtual int intersectionWithDetElement( int detElementID, gear::Vector3D& point, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to sensitive detector element,
     *  returning intersection point in global coordinates via reference 
     */
    virtual int intersectionWithDetElement( int detEementID, EVENT::TrackerHit* hit, gear::Vector3D& point, int mode=modeClosest ) = 0  ;
    
    
  protected:
    
  private:
    
    IMarlinTrack& operator=( const IMarlinTrack&) ; // disallow assignment operator 
    
  } ;
  
  
  
  
  /** Helper function to convert error return code to string */
  inline std::string errorCode( int error ){
    
    switch ( error ){ 
      case IMarlinTrack::success           : return "IMarlinTrack::success";         break;
      case IMarlinTrack::error             : return "IMarlinTrack::error";           break;
      case IMarlinTrack::bad_intputs       : return "IMarlinTrack::bad_intputs";     break;
      case IMarlinTrack::no_intersection   : return "IMarlinTrack::no_intersection"; break;
      case IMarlinTrack::site_discarded    : return "IMarlinTrack::site_discarded";  break;
      default: return "UNKNOWN" ; 
    }
  }
  
  
  
  
} // end of MarlinTrk namespace 

#endif

