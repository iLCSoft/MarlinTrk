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
}

namespace UTIL {
  class BitField64;
}


namespace MarlinTrk{
  class IMarlinTrack ;
}


namespace MarlinTrk{

  /** Takes a list of hits and uses the IMarlinTrack inferface to fit them using a supplied prefit containing a covariance matrix for the initialisation. The TrackImpl will have the 4 trackstates added to it @IP, @First_Hit, @Last_Hit and @CaloFace */
  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackImpl* track, bool fit_backwards, IMPL::TrackStateImpl* pre_fit, float bfield_z, double maxChi2Increment=DBL_MAX);
  
  /** Takes a list of hits and uses the IMarlinTrack inferface to fit them using a supplied covariance matrix for the initialisation. The TrackImpl will have the 4 trackstates added to it @IP, @First_Hit, @Last_Hit and @CaloFace */
  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackImpl* track, bool fit_backwards, const EVENT::FloatVec& initial_cov_for_prefit, float bfield_z, double maxChi2Increment=DBL_MAX);
  
  /** Provides the values of a track state from the first, middle and last hits in the hit_list. */
  int createPrefit( std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackStateImpl* pre_fit, float bfield_z, bool fit_backwards );

  /** Takes a list of hits and uses the IMarlinTrack inferface to fit them using a supplied prefit containing a covariance matrix for the initialisation. */  
  int createFit( std::vector<EVENT::TrackerHit*>& hit_list, IMarlinTrack* marlinTrk, IMPL::TrackStateImpl* pre_fit, float bfield_z, bool fit_backwards, double maxChi2Increment=DBL_MAX );

  /** Set the subdetector hit numbers for the TrackImpl */
  void addHitNumbersToTrack(IMPL::TrackImpl* track, std::vector<EVENT::TrackerHit*>& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder);

  /** Set the subdetector hit numbers for the TrackImpl */
  void addHitNumbersToTrack(IMPL::TrackImpl* track, std::vector<std::pair<EVENT::TrackerHit* , double> >& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder);
  
}

#endif
