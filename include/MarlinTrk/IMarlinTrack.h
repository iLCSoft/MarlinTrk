#ifndef IMarlinTrack_h
#define IMarlinTrack_h

#include <cfloat>

#include "lcio.h"

#include "EVENT/TrackerHit.h"
#include "IMPL/TrackStateImpl.h"

#include "DDRec/Vector3D.h"

#include <exception>
#include <string>


#ifdef MARLINTRK_BACKWARD_GEAR_WRAPPERS
#include "gearimpl/Vector3D.h"
#endif



namespace MarlinTrk{
  
  /// the Vector3D used for the tracking interface
  typedef dd4hep::rec::Vector3D Vector3D ;
  
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
    static const bool backward ;
    
    /** boolean constant for defining backward direction - to be used for intitialise */
    static const bool forward ;
    
    
    
    static const int modeBackward ;
    static const int modeClosest  ;
    static const int modeForward  ;
    
    
    static const int success ;  // no error
    static const int error ;
    static const int bad_intputs ;
    static const int no_intersection ; // no intersection found
    static const int site_discarded ;  // measurement discarded by the fitter
    static const int site_fails_chi2_cut ;  // measurement discarded by the fitter due to chi2 cut
    static const int all_sites_fail_fit ;   // no single measurement added to the fit
    
    
    /**default d'tor*/
    virtual ~IMarlinTrack() {};
    
    /** set the mass of the charged particle (GeV) that is used for energy loss and multiple scattering -
     * default value if this method is not called is the pion mass. 
     */
    virtual void setMass(double mass) = 0 ;

    /** return the of the charged particle (GeV) that is used for energy loss and multiple scattering.
     */
    virtual double getMass() = 0 ;

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
    virtual int fit( double maxChi2Increment=DBL_MAX ) = 0 ;
    
    
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
    
    /** get the current number of degrees of freedom for the fit.
     */
    virtual int getNDF( int& ndf ) = 0 ;
    
    /** get TrackeHit at which fit became constrained, i.e. ndf >= 0
     */
    virtual int getTrackerHitAtPositiveNDF( EVENT::TrackerHit*& trkhit ) = 0 ;
    
    // PROPAGATORS 
    
    /** propagate the fit to the point of closest approach to the given point, returning TrackState, chi2 and ndf via reference    
     */
    virtual int propagate( const Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

    
    /** propagate the fit at the measurement site associated with the given hit, to the point of closest approach to the given point,
     *  returning TrackState, chi2 and ndf via reference   
     */
    virtual int propagate( const Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;


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
    virtual int extrapolate( const Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to the point of closest approach to the given point, 
     *  returning TrackState, chi2 and ndf via reference   
     */
    virtual int extrapolate( const Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
    
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
    virtual int intersectionWithLayer( int layerID, Vector3D& point, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to numbered sensitive layer,
     *  returning intersection point in global coordinates and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int intersectionWithLayer( int layerID, EVENT::TrackerHit* hit, Vector3D& point, int& detElementID, int mode=modeClosest ) = 0  ;

    
    /** extrapolate the fit to numbered sensitive detector element, returning intersection point in global coordinates via reference 
     */
    virtual int intersectionWithDetElement( int detElementID, Vector3D& point, int mode=modeClosest ) = 0  ;


    /** extrapolate the fit at the measurement site associated with the given hit, to sensitive detector element,
     *  returning intersection point in global coordinates via reference 
     */
    virtual int intersectionWithDetElement( int detEementID, EVENT::TrackerHit* hit, Vector3D& point, int mode=modeClosest ) = 0  ;
    
    /** Dump this track to a string for debugging - implementation dependant.
     */
    virtual std::string toString() ;


    //-------------------------------------------------------------------------------------------------------------------------------

#ifdef MARLINTRK_BACKWARD_GEAR_WRAPPERS

    int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){
      Vector3D v( point.x(), point.y(), point.z() ) ;
      return propagate( v, ts, chi2, ndf ) ;
    }
    int propagate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){
      Vector3D v( point.x(), point.y(), point.z() ) ;
      return propagate( v, hit, ts, chi2, ndf ) ;
    }
    int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){
      Vector3D v( point.x(), point.y(), point.z() ) ;
      return extrapolate( v, ts, chi2, ndf ) ;
    }
    int extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ){
      Vector3D v( point.x(), point.y(), point.z() ) ;
      return extrapolate( v, hit, ts, chi2, ndf ) ;
    }
    int intersectionWithLayer( int layerID, gear::Vector3D& point, int& detElementID, int mode=modeClosest ){
      Vector3D v( point.x(), point.y(), point.z() ) ;
      int ret =  intersectionWithLayer( layerID, v, detElementID, mode ) ;
      point = v ;
      return ret ;
    }
    int intersectionWithLayer( int layerID, EVENT::TrackerHit* hit, gear::Vector3D& point, int& detElementID, int mode=modeClosest ){
      Vector3D v( point.x(), point.y(), point.z() ) ;
      int ret = intersectionWithLayer( layerID, hit, v, detElementID, mode ) ;
      point = v ;
      return ret ;
    }
    int intersectionWithDetElement( int detElementID, gear::Vector3D& point, int mode=modeClosest ){
      Vector3D v( point.x(), point.y(), point.z() ) ;
      int ret = intersectionWithDetElement( detElementID, v, mode ) ;
      point = v ;
      return ret ;
    }   
    int intersectionWithDetElement( int detElementID, EVENT::TrackerHit* hit, gear::Vector3D& point, int mode=modeClosest ){
      Vector3D v( point.x(), point.y(), point.z() ) ;
      int ret = intersectionWithDetElement( detElementID, hit, v, mode ) ;
      point = v ;
      return ret ;
    }

#endif

    //-------------------------------------------------------------------------------------------------------------------------------

    
  protected:
    
  private:
    
    IMarlinTrack& operator=( const IMarlinTrack&) ; // disallow assignment operator 
                
  } ;
  
  /** Helper function to convert error return code to string */
  std::string errorCode( int error );
  
  
  
  
} // end of MarlinTrk namespace 

#endif

