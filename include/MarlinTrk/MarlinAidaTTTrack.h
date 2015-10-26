#ifndef MarlinAidaTTTrack_h
#define MarlinAidaTTTrack_h

#include "IMarlinTrack.h"
#include "IMarlinTrkSystem.h"

#include "aidaTT/AidaTT.hh"

#include <cmath>

namespace MarlinTrk {
  class MarlinAidaTT;
}

namespace EVENT{
  class TrackerHit ;
}



/** Implementation of the IMarlinTrack interface, using KalTest and KalDet to provide 
 *  the needed functionality for a Kalman Filter.
 *
 * @version $Id: $
 * @autor F. Gaede DESY, Oct 2016
 */

namespace MarlinTrk{

class MarlinAidaTTTrack : public MarlinTrk::IMarlinTrack {
  
public:

    
  MarlinAidaTTTrack(MarlinAidaTT* aidaTT) ;
  
  ~MarlinAidaTTTrack() ;
  
protected:
  
private:
  
  MarlinAidaTTTrack(const MarlinAidaTTTrack&) ;                 // Prevent copy-construction
  MarlinAidaTTTrack& operator=(const MarlinAidaTTTrack&) ;      // Prevent assignment
  
  // make member functions private to force use through interface
  
  /** add hit to track - the hits have to be added ordered in time ( i.e. typically outgoing )
   *  this order will define the direction of the energy loss used in the fit
   */
  int addHit(EVENT::TrackerHit* hit) ;
  
  /** initialise the fit using the hits added up to this point -
   *  the fit direction has to be specified using IMarlinTrack::backward or IMarlinTrack::forward. 
   *  this is the order  wrt the order used in addHit() that will be used in the fit() 
   */
  int initialise( bool fitDirection ); 
  
  /** initialise the fit with a track state
   *  the fit direction has to be specified using IMarlinTrack::backward or IMarlinTrack::forward. 
   *  this is the order that will be used in the fit().
   *  it is the users responsibility that the track state is consistent with the order
   *  of the hits used in addHit() ( i.e. the direction of energy loss )
   *  Note: the bfield_z is not taken from the argument but from the first hit 
   *  should consider changing the interface ...
   */
  int initialise( const EVENT::TrackState& ts, double /*bfield_z*/, bool fitDirection ) ;
  
  
  /** perform the fit of all current hits, returns error code ( IMarlinTrack::success if no error ) .
   *  the fit will be performed  in the order specified at initialise() wrt the order used in addHit(), i.e.
   *  IMarlinTrack::backward implies fitting from the outside to the inside for tracks comming from the IP.
   */
  int fit( double maxChi2Increment=DBL_MAX ) ;
  
  
  /** smooth all track states 
   */
  int smooth() ;
  
  
  /** smooth track states from the last filtered hit back to the measurement site associated with the given hit 
   */
  int smooth( EVENT::TrackerHit* hit )  ;
  
  
  /** update the current fit using the supplied hit, return code via int. Provides the Chi2 increment to the fit from adding the hit via reference. 
   *  the given hit will not be added if chi2increment > maxChi2Increment. 
   */
  int addAndFit( EVENT::TrackerHit* hit, double& chi2increment, double maxChi2Increment=DBL_MAX ) ;
  
  
  /** obtain the chi2 increment which would result in adding the hit to the fit. This method will not alter the current fit, and the hit will not be stored in the list of hits or outliers
   */
  int testChi2Increment( EVENT::TrackerHit* hit, double& chi2increment ) ;

  
  // Track State Accessesors
  
  /** get track state, returning TrackState, chi2 and ndf via reference 
   */
  int getTrackState( IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
  
  
  /** get track state at measurement associated with the given hit, returning TrackState, chi2 and ndf via reference 
   */
  int getTrackState( EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
  
  
  /** get the list of hits included in the fit, together with the chi2 contributions of the hits. 
   *  Pointers to the hits together with their chi2 contribution will be filled into a vector of 
   *  pairs consitining of the pointer as the first part of the pair and the chi2 contribution as
   *  the second.
   */
  int getHitsInFit( std::vector<std::pair<EVENT::TrackerHit*, double> >& hits ) ;
  
  /** get the list of hits which have been rejected by from the fit due to the a chi2 increment greater than threshold,
   *  Pointers to the hits together with their chi2 contribution will be filled into a vector of 
   *  pairs consitining of the pointer as the first part of the pair and the chi2 contribution as
   *  the second.
   */
  int getOutliers( std::vector<std::pair<EVENT::TrackerHit*, double> >& hits ) ;


  /** get the current number of degrees of freedom for the fit.
   */
  int getNDF( int& ndf ) ;
  
  /** get TrackeHit at which fit became constrained, i.e. ndf >= 0
   */
  int getTrackerHitAtPositiveNDF( EVENT::TrackerHit*& trkhit ) ;
  
  // PROPAGATORS 
  
  /** propagate the fit to the point of closest approach to the given point, returning TrackState, chi2 and ndf via reference    
   */
  int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
  
  
  /** propagate the fit at the measurement site associated with the given hit, to the point of closest approach to the given point,
   *  returning TrackState, chi2 and ndf via reference   
   */
  int propagate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
  
  
  /** propagate the fit to the numbered sensitive layer, returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference 
   */
  int propagateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest )  ;
  
  /** propagate the fit at the measurement site associated with the given hit, to numbered sensitive layer, 
   *  returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference 
   */
  int propagateToLayer( int layerID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest )  ;
  
  /** propagate the fit to sensitive detector element, returning TrackState, chi2 and ndf via reference
   */
  int propagateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
  
  /** propagate the fit at the measurement site associated with the given hit, to sensitive detector element, 
   *  returning TrackState, chi2 and ndf via reference 
   */
  int propagateToDetElement( int detEementID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
  
  
  
  // EXTRAPOLATORS
  
  /** extrapolate the fit to the point of closest approach to the given point, returning TrackState, chi2 and ndf via reference   
   */
  int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
  
  /** extrapolate the fit at the measurement site associated with the given hit, to the point of closest approach to the given point, 
   *    returning TrackState, chi2 and ndf via reference   
   */
  int extrapolate( const gear::Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
  
  
  /** extrapolate the fit to numbered sensitive layer, returning TrackState via provided reference 
   */
  int extrapolateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest )  ;
  
  /** extrapolate the fit at the measurement site associated with the given hit, to numbered sensitive layer, 
   *  returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference 
   */
  int extrapolateToLayer( int layerID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest )  ;
  
  /** extrapolate the fit to sensitive detector element, returning TrackState, chi2 and ndf via reference
   */
  int extrapolateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
  
  /** extrapolate the fit at the measurement site associated with the given hit, to sensitive detector element, 
   *  returning TrackState, chi2 and ndf via reference 
   */
  int extrapolateToDetElement( int detEementID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
  
  
  /** Dump this track to a string for debugging.
   */
  std::string toString() ;
 
  // INTERSECTORS
  
  
  /** extrapolate the fit to numbered sensitive layer, returning intersection point in global coordinates and integer ID of the 
   *  intersected sensitive detector element via reference 
   */
  int intersectionWithLayer( int layerID, gear::Vector3D& point, int& detElementID, int mode=modeClosest )  ;
  
  /** extrapolate the fit at the measurement site associated with the given hit, to numbered sensitive layer,
   *  returning intersection point in global coordinates and integer ID of the intersected sensitive detector element via reference 
   */
  int intersectionWithLayer( int layerID, EVENT::TrackerHit* hit, gear::Vector3D& point, int& detElementID, int mode=modeClosest )  ;
  
  /** extrapolate the fit to numbered sensitive detector element, returning intersection point in global coordinates via reference 
   */
  int intersectionWithDetElement( int detElementID, gear::Vector3D& point, int mode=modeClosest )  ;
  
  /** extrapolate the fit at the measurement site associated with the given hit, to sensitive detector element,
   *  returning intersection point in global coordinates via reference 
   */ 
  int intersectionWithDetElement( int detElementID, EVENT::TrackerHit* hit, gear::Vector3D& point, int mode=modeClosest )  ;
  
  
  //** end of memeber functions from IMarlinTrack interface
  
  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange( double phi) const {
    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
    return phi ;
  }
  
protected:
  
  aidaTT::trackParameters createPreFit(aidaTT::trackParameters& tp ) ;

  void getHitInfo( const EVENT::TrackerHit* hit, double* hitpos, std::vector<double>& precision, const aidaTT::ISurface* surf ) ;

  int getTrackState( const aidaTT::Vector3D& refPoint, int label, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;

  /// common initialization
  int myInit() ;

  // memeber variables 
  
  MarlinAidaTT* _aidaTT;
  
  /** used to store whether initial track state has been supplied or created 
   */
  bool _initialised ;
  
  /** used to store the fit direction supplied to intialise 
   */
  bool _fitDirection ;
  
  /** used to store whether smoothing has been performed
   */
  bool _smoothed ;
  
  std::vector< EVENT::TrackerHit*> _lcioHits ;
 
  aidaTT::trackParameters _initialTrackParams ;

  aidaTT::trajectory* _fitTrajectory ;
  
  const std::vector<std::pair<double, const aidaTT::ISurface*> >* _intersections ;

  std::map< int, int > _indexMap ;
} ;

} // end of namespace MarlinTrk
  
#endif
