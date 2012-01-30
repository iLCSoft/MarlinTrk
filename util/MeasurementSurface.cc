
#include "MeasurementSurface.h"
#include "ICoordinateSystem.h"

namespace MarlinTrk {

  namespace GearExtensions {
    MeasurementSurface::~MeasurementSurface(){ delete _coordinateSystem; }
  }
  
}
