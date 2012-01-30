#ifndef CARTESIANCOORDINATESYSTEM_h
#define CARTESIANCOORDINATESYSTEM_h

#include "ICoordinateSystem.h"

namespace GearExtensions{
   
/** Cartesian coordinate system class */
class CartesianCoordinateSystem : public ICoordinateSystem{
   
   
public:
   
   CartesianCoordinateSystem( TVector3 T , TRotation R ): _T(T), _R(R){}
   
   /** @return the local coordinates of the point */
   virtual TVector3 getLocalPoint( TVector3 globalPoint );
   
   /** @return the global coordinates of the point */
   virtual TVector3 getGlobalPoint( TVector3 localPoint );
   
   /** @return the global coordinates of the origin of the coordinate system */
   virtual TVector3 getOrigin(){ return _T; }
   
   /** @return a rotation Matrix. local = R*global (after the translation)
      */
   TRotation getR(){ return _R;}
   
private:
   
   /** The translation vector (= the origin of the Coordinate System ) */
   TVector3 _T;
   
   /** Rotation Matrix
    * Definition: global = R* local
    */
   TRotation _R;
   
};
   
} // end namespace

#endif
   
