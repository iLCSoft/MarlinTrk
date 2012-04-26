#ifndef INCLUDE_MarlinTrkUtils
#define INCLUDE_MarlinTrkUtils 1

#include <vector>

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

  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackImpl* track, bool fit_backwards, IMPL::TrackStateImpl* pre_fit, float bfield_z);
  
  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackImpl* track, bool fit_backwards, const EVENT::FloatVec& initial_cov_for_prefit, float bfield_z);
  
  int createPrefit( std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackStateImpl* pre_fit, float bfield_z, bool fit_backwards );
  
  int createFit( std::vector<EVENT::TrackerHit*>& hit_list, IMarlinTrack* marlinTrk, IMPL::TrackStateImpl* pre_fit, float bfield_z, bool fit_backwards );
  
  void addHitsToTrack(IMPL::TrackImpl* track, std::vector<EVENT::TrackerHit*>& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder);

  void addHitsToTrack(IMPL::TrackImpl* track, std::vector<std::pair<EVENT::TrackerHit* , double> >& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder);
  
}

#endif
