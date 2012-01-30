#ifndef MeasurementSurfaceStore_h
#define MeasurementSurfaceStore_h

/** MeasurementSurfaceStore: Class to hold and manage a collection of MeasurementSufaces  
 *
 * @author S.Aplin DESY
 */

#include <map>
#include <exception>

namespace GearExtensions{
  class MeasurementSurface;
}

namespace gear{
  class GearMgr;
  class ZPlanarParameters;
  class FTDParameters;
}


namespace GearExtensions{
  
  
  class MeasurementSurfaceStoreException: public std::exception {
    virtual const char* what() const throw() {
      return "MeasurementSurfaceStoreException occurred";
    }
  } ;
  
  
  class MeasurementSurfaceStore {
    
  public:
    
    /** Accessor Method */
    static MeasurementSurfaceStore& Instance() {
      
      static MeasurementSurfaceStore singleton;
      
      return singleton;
      
    }
    
    // Other non-static member functions
    
  public:
    
    /** Destructor */
    ~MeasurementSurfaceStore();   
    
    void initialise(gear::GearMgr* gear_mgr) ;
    
    /** Get Measurement Surface via ID */
    MeasurementSurface* GetMeasurementSurface( int ID ) const ;  
    
    
  private:
    
    MeasurementSurfaceStore() { _measurement_surface_map.clear() ;}                               // Private constructor
    
    
    
    MeasurementSurfaceStore(const MeasurementSurfaceStore&) ;                 // Prevent copy-construction
    MeasurementSurfaceStore& operator=(const MeasurementSurfaceStore&) ;      // Prevent assignment
    
    void addMeasurementSurface(MeasurementSurface* ms); 
    
    void createStore(gear::GearMgr* gear_mgr);
    
    /** adds MeasurementSufaces to the store
     * @param param: the ZPlanarParameters pointer of the detector, of which the measurement surfaces shall be added
     * 
     * @param det_id: the detector id (as in ILDConf)
     */
    void storeZPlanar( const gear::ZPlanarParameters* param , int det_id );
    
    void storeFTD( const gear::FTDParameters* param );
    
    // private member variables
    std::map<int,MeasurementSurface* > _measurement_surface_map;
    
    typedef std::map<int, MeasurementSurface*>::const_iterator ms_map_it ; 
    
    static bool _isInitialised;
    
  };
  
}

#endif
