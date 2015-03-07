#ifndef INCLUDE_MarlinTrkUtils
#define INCLUDE_MarlinTrkUtils 1

#include <vector>

#include <cfloat>

#include <LCIOSTLTypes.h>

namespace IMPL {
  class TrackImpl ;
  class TrackStateImpl ;
}

namespace EVENT {
  class TrackerHit;
  class TrackState;
}

namespace UTIL {
  class BitField64;
}


namespace MarlinTrk{
  class IMarlinTrack ;
}


namespace MarlinTrk{

  /** Takes a list of hits and uses the IMarlinTrack inferface to fit them using a supplied prefit containing
   *  a covariance matrix for the initialisation. The TrackImpl will have the 4 trackstates added to
   *  it @IP, @First_Hit, @Last_Hit and @CaloFace */
  int createFinalisedLCIOTrack(
      IMarlinTrack* marlinTrk,
      std::vector<EVENT::TrackerHit*>& hit_list,
      IMPL::TrackImpl* track,
      bool fit_direction,
      EVENT::TrackState* pre_fit,
      float bfield_z,
      double maxChi2Increment=DBL_MAX);
  
  /** Takes a list of hits and uses the IMarlinTrack inferface to fit them using a supplied covariance matrix
   *  for the initialisation. The TrackImpl will have the 4 trackstates added to
   *  it @IP, @First_Hit, @Last_Hit and @CaloFace */
  int createFinalisedLCIOTrack(
      IMarlinTrack* marlinTrk,
      std::vector<EVENT::TrackerHit*>& hit_list,
      IMPL::TrackImpl* track,
      bool fit_direction,
      const EVENT::FloatVec& initial_cov_for_prefit,
      float bfield_z,
      double maxChi2Increment=DBL_MAX);
  
  /** Provides the values of a track state from the first, middle and last hits in the hit_list. */
  int createPrefit( std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackStateImpl* pre_fit, float bfield_z, bool fit_direction );

  /** Takes a list of hits and uses the IMarlinTrack inferface to fit them using a supplied prefit containing a covariance matrix for the initialisation. */  
  int createFit( std::vector<EVENT::TrackerHit*>& hit_list, IMarlinTrack* marlinTrk, EVENT::TrackState* pre_fit, float bfield_z, bool fit_direction, double maxChi2Increment=DBL_MAX );

  /** Takes a fitted MarlinTrack, TrackImpl to record the fit and the hits which have been added to the fit.
   *  The TrackImpl will have the 4 trackstates added to it @IP, @First_Hit, @Last_Hit and @CaloFace.
   *  Note: the hit list is needed as the IMarlinTrack only contains the hits used in the fit, not the spacepoints
   *  (if any have been included) so as the strip hits cannot point to the space points we need to have the list so
   *  that they can be recorded in the LCIO TrackImpl */
  int finaliseLCIOTrack(
      IMarlinTrack* marlinTrk,
      IMPL::TrackImpl* track,
      std::vector<EVENT::TrackerHit*>& hit_list,
      bool fit_direction,
      IMPL::TrackStateImpl* atLastHit=0,
      IMPL::TrackStateImpl* atCaloFace=0);
  
  /** Set the subdetector hit numbers for the TrackImpl */
  void addHitNumbersToTrack(IMPL::TrackImpl* track, std::vector<EVENT::TrackerHit*>& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder);

  /** Set the subdetector hit numbers for the TrackImpl */
  void addHitNumbersToTrack(IMPL::TrackImpl* track, std::vector<std::pair<EVENT::TrackerHit* , double> >& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder);
  
}

#endif
