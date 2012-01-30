#ifndef ICOORDINATESYSTEM_h
#define ICOORDINATESYSTEM_h

#include "TVector3.h"
#include "TRotation.h"

namespace GearExtensions{
  
  /** An abstract base class for coordinate systems
   */
  class ICoordinateSystem{
    
    
  public:
    
    virtual ~ICoordinateSystem() { /* no-op */ }
    
    /** @return the local coordinates of the point */
    virtual TVector3 getLocalPoint( TVector3 globalPoint ) = 0;
    
    /** @return the global coordinates of the point */
    virtual TVector3 getGlobalPoint( TVector3 localPoint ) = 0;
    
    /** @return the global coordinates of the origin of the coordinate system */
    virtual TVector3 getOrigin() = 0;
    
    
    
  };
  
} // end namespace

#endif

