#ifndef MEASUREMENTSURFACE_H
#define MEASUREMENTSURFACE_H

#include <vector>
#include "IBoundary.h"

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
      
      /** adds a boundary */
      void addBoundary( IBoundary* boundary ){ _boundaries.push_back( boundary ); }
      
      /** Checks if a point in local coordinates is within the boundaries */
      bool isLocalInBoundary( CLHEP::Hep3Vector local );
      
      /** Checks if a point in global coordinates is within the boundaries */
      bool isGlobalInBoundary( CLHEP::Hep3Vector global );
      
    private:
      
      int _ID;
      ICoordinateSystem* _coordinateSystem;
      
      std::vector<IBoundary*> _boundaries;
      
      MeasurementSurface(const MeasurementSurface& m){};   // copy constructor is private --> no cpoying allowed
      MeasurementSurface& operator= (MeasurementSurface const& m); // assignment not allowed either
      
      
    };
    
    
  } //end of GearExtensions namespace 
  
} //end of MarlinTrk namespace 

#endif
