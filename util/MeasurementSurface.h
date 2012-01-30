#ifndef MEASUREMENTSURFACE_H
#define MEASUREMENTSURFACE_H

//#include "ICoordinateSystem.h"

namespace MarlinTrk {
  
  /** Measurement Surface Class */
  namespace GearExtensions{
    
    class ICoordinateSystem;
    
    class MeasurementSurface{
      
    public:
      
      MeasurementSurface( int ID, ICoordinateSystem* csys ): _ID( ID ), _coordinateSystem( csys ){}
      
      ~MeasurementSurface() ;
      
      int getID(){ return _ID; }
      ICoordinateSystem* getCoordinateSystem(){ return _coordinateSystem; }
      
    private:
      
      int _ID;
      ICoordinateSystem* _coordinateSystem;
      
      
      MeasurementSurface(const MeasurementSurface& m){};   // copy constructor is private --> no cpoying allowed
      MeasurementSurface& operator= (MeasurementSurface const& m); // assignment not allowed either
      
      
    };
    
    
  } //end of GearExtensions namespace 
  
} //end of MarlinTrk namespace 

#endif
