
#include "MeasurementSurface.h"
#include "ICoordinateSystem.h"

namespace MarlinTrk {

  namespace GearExtensions {
    MeasurementSurface::~MeasurementSurface(){ 
      
      delete _coordinateSystem; 
      for( unsigned i=0; i< _boundaries.size(); i++) delete _boundaries[i];
      _boundaries.clear();
      
    }
    
    
    bool MeasurementSurface::isLocalInBoundary( CLHEP::Hep3Vector local ){
      
      
      for( unsigned i=0; i<_boundaries.size(); i++ ){
        
        if ( _boundaries[i]->isInBoundary( local ) ) return true;  // if it's within one of the surfaces, return true
        
      }
      
      return false; // was within no surface
      
    }
    
    
    
    bool MeasurementSurface::isGlobalInBoundary( CLHEP::Hep3Vector global ){
      
      
      CLHEP::Hep3Vector local = _coordinateSystem->getLocalPoint( global );
      
      return isLocalInBoundary( local );
      
    }
    
    
    
  }
  
}
