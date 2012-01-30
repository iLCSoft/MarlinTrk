#include "CartesianCoordinateSystem.h"

using namespace GearExtensions;

TVector3 CartesianCoordinateSystem::getLocalPoint( TVector3 globalPoint ){
   
   // First we get us a new origin via translation
   TVector3 x = globalPoint  - _T;
   
   // Then we do the rotation
   TRotation R_inv = _R.Inverse();
   
   return R_inv*x;   
   
}

TVector3 CartesianCoordinateSystem::getGlobalPoint( TVector3 localPoint ){
   
   
   // The point in global coordinates is the place of the origin in gloabl coordinates (=T) plus
   // the local coordinates rotated into the global ones   
   return _T + _R*localPoint;
   
}
