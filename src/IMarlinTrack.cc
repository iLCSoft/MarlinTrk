
#include "MarlinTrk/IMarlinTrack.h"

namespace MarlinTrk{

  const int IMarlinTrack::modeBackward = - 1 ;
  const int IMarlinTrack::modeClosest  =   0 ;
  const int IMarlinTrack::modeForward  = + 1 ;


 const int IMarlinTrack::success  = 0 ;  // no error
 const int IMarlinTrack::error = 1 ;
 const int IMarlinTrack::bad_intputs = 3 ;
 const int IMarlinTrack::no_intersection = 4 ; // no intersection found
 const int IMarlinTrack::site_discarded = 5 ;  // measurement discarded by the fitter
 const int IMarlinTrack::site_fails_chi2_cut = 6 ;  // measurement discarded by the fitter due to chi2 cut
 const int IMarlinTrack::all_sites_fail_fit = 7 ;   // no single measurement added to the fit

}
