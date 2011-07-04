#ifndef INCLUDE_ILDPlanarHitExt
#define INCLUDE_ILDPlanarHitExt 1

namespace ILDPlanarHitExt{

  // helper class to assign additional parameters to TrackerHits
  struct PlanarHitExtStruct{
  PlanarHitExtStruct() : u(0.0), du(0.0), v(0.0), dv(0.0), ladder_r(0), ladder_phi(0), ladder_offset(0), ladder_number(-1), ladder_length(0.0), ladder_width(0.0)  {}
    double u ;
    double du ;
    double v ;
    double dv ;
    double ladder_r ;
    double ladder_phi ;
    double ladder_offset ;
    int    ladder_number ;
    double ladder_length ;
    double ladder_width  ;

  } ;
  struct PlanarHitExt : lcio::LCOwnedExtension<PlanarHitExt, PlanarHitExtStruct> {} ;
  
}

#endif
