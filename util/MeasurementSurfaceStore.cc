
#include "MeasurementSurfaceStore.h"

#include <stdexcept>
#include <vector>
#include <algorithm>

#include "MeasurementSurface.h"

#include <gear/ZPlanarParameters.h>
#include <gear/ZPlanarLayerLayout.h>
#include <gear/FTDLayerLayout.h>
#include <gear/FTDParameters.h>
#include <gear/GEAR.h>
#include <gear/GearMgr.h>

#include "UTIL/ILDConf.h"

#include "CartesianCoordinateSystem.h"

#include "streamlog/streamlog.h"

namespace MarlinTrk {
  
  namespace GearExtensions{
    
    
    void printRotation( const CLHEP::HepRotation& R ){
      
      streamlog_out(DEBUG4).precision(2);
      streamlog_out(DEBUG4).setf( std::ios::fixed , std::ios::floatfield );
      streamlog_out(DEBUG4) << "R:\n"
      << "| " << R.xx() << "  " << R.xy() << "  " << R.xz() << " |\n"
      << "| " << R.yx() << "  " << R.yy() << "  " << R.yz() << " |\n"
      << "| " << R.zx() << "  " << R.zy() << "  " << R.zz() << " |\n";
    }
    
    bool MeasurementSurfaceStore::_isInitialised = false ;
    
    MeasurementSurfaceStore::~MeasurementSurfaceStore(){
      
      ms_map_it it = _measurement_surface_map.begin();
      std::vector<MeasurementSurface*> deleted_objects;
      
      for( /**/; it!=_measurement_surface_map.end(); ++it) 
        
        if( std::find( deleted_objects.begin(), deleted_objects.end(), (*it).second ) != deleted_objects.end() ) {
          delete (*it).second ;
          deleted_objects.push_back((*it).second) ;
        }
    }
    
    MeasurementSurface* MeasurementSurfaceStore::GetMeasurementSurface(int ID) const {
      
      ms_map_it it = _measurement_surface_map.find(ID) ;        
      
      if ( it == _measurement_surface_map.end() ) { 
        MeasurementSurfaceStoreException exp;
        throw exp ; 
      } 
      else { 
        return (*it).second ; 
      }
      
    }
    
    void MeasurementSurfaceStore::initialise(gear::GearMgr* gear_mgr){
      
      if ( ! _isInitialised) {
        this->createStore(gear_mgr); 
        _isInitialised = true ;
      }
      
      
    }
    
    void MeasurementSurfaceStore::addMeasurementSurface(MeasurementSurface* ms) {
      
      int ID = ms->getID();
      
      ms_map_it it = _measurement_surface_map.find(ID) ; 
      
      MeasurementSurfaceStoreException exp;
      
      if ( it != _measurement_surface_map.end() ) { 
        throw exp; 
      } 
      else { 
        _measurement_surface_map[ID] = ms  ; 
      }
    }
    
    
    void MeasurementSurfaceStore::createStore(gear::GearMgr* gear_mgr){
      
      
      const gear::ZPlanarParameters* paramVXD = &(gear_mgr->getVXDParameters());
      storeZPlanar( paramVXD , UTIL::ILDDetID::VXD );
      
      const gear::ZPlanarParameters* paramSIT = &(gear_mgr->getSITParameters());
      storeZPlanar( paramSIT , UTIL::ILDDetID::SIT );
      
      const gear::ZPlanarParameters* paramSET = &(gear_mgr->getSETParameters());
      storeZPlanar( paramSET , UTIL::ILDDetID::SET );
      
      const gear::FTDParameters* paramFTD = &(gear_mgr->getFTDParameters());
      storeFTD( paramFTD );
      
    }
    
    
    
    
    
    void MeasurementSurfaceStore::storeZPlanar( const gear::ZPlanarParameters* param , int det_id ){
      
      
      const gear::ZPlanarLayerLayout& layerLayout = param->getZPlanarLayerLayout();
      
      unsigned nLayers = layerLayout.getNLayers();
      
      for( unsigned layerNumber = 0; layerNumber < nLayers; layerNumber++ ){
        
        unsigned nLadders = layerLayout.getNLadders( layerNumber );
        double ladder_r            = layerLayout.getSensitiveDistance(layerNumber); // the distance of the ladders from (0,0,0)
        double sensitive_offset    = layerLayout.getSensitiveOffset(layerNumber); // the offset, see ZPlanarLayerLayout.h for more details
        double deltaPhi            = ( 2 * M_PI ) / nLadders ; // the phi difference between two ladders
        double phi0                = layerLayout.getPhi0( layerNumber );
        double stripAngle = 0.; // TODO: implement
        
        
        
        streamlog_out(DEBUG4) << "NLADDERS = " << nLadders << "\n";
        
        for( unsigned ladderNumber = 0; ladderNumber < nLadders; ladderNumber++ ){
          
          // determine the CellID0 for the  ladder
          UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
          cellID[ lcio::ILDCellID0::subdet ] = det_id ;
          cellID[ lcio::ILDCellID0::side   ] = 0 ;
          cellID[ lcio::ILDCellID0::layer  ] = layerNumber ;
          cellID[ lcio::ILDCellID0::module ] = ladderNumber ;
          cellID[ lcio::ILDCellID0::sensor ] = 0 ;
          int cellID0 = cellID.lowWord();
          
          streamlog_out(DEBUG4) << "layer = " << layerNumber << "\tladder = " << ladderNumber << "\n";
          
          // Let's start with the translation T: the new center of coordinates:
          // The center of the first ladder (when we ignore an offset and phi0 for now) is (R,0,0)
          // If we include the offset, the center gets shifted to (R,offset,0)
          CLHEP::Hep3Vector T( ladder_r, sensitive_offset, 0 );
          // Now we have to take into account phi0 and that the number of the ladder.
          // Together the center is rotated by phi0 + ladderNumber*deltaPhi around the z axis
          CLHEP::HepRotation rot;
          rot.rotateZ( deltaPhi * ladderNumber + phi0 );
          
          T = rot * T;
          
          // Next, we want to determinte the rotation matrix R
          // We start with u,v,w alligned with x,y,z.
          // As u is perpendicular to the strip orientation it looks like this.
          //               y
          //               |     
          //            ---|---
          //            |  |  |
          //            |  |  |
          //      <--------|--------> x
          //            |  |  |
          //            |  |  |
          //            ---|---
          //               |
          // To get this sensor in place we have to do a few rotations:
          // First we'll rotate around the z axis. With a strip angle of 0,
          // we would just rotate by 90°, but with a strip angle by
          // 90°-stripAngle in clockwise direction. 
          CLHEP::HepRotation R;
          printRotation( R );        
          R.rotateZ( stripAngle - M_PI/2. );
          printRotation( R );
          
          // Next we rotate 90° clockwise around y, so the strip now points in z direction (if strip angle == 0)
          R.rotateY( -M_PI/2. );
          
          // Finally we have to get the ladder in place w.r. to its number and the resulting phi angle
          R.rotateZ( deltaPhi * ladderNumber + phi0 );
          
          CartesianCoordinateSystem* cartesian = new CartesianCoordinateSystem( T, R );
          
          
          MeasurementSurface* ms = new MeasurementSurface( cellID0, cartesian );
          addMeasurementSurface( ms );
          
        }
        
      }
      
    }
    
    
    void MeasurementSurfaceStore::storeFTD( const gear::FTDParameters* param ){
      
      const gear::FTDLayerLayout& ftdLayers = param->getFTDLayerLayout() ;
      unsigned nLayers = ftdLayers.getNLayers();
      
      UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
      cellID[ lcio::ILDCellID0::subdet ] = UTIL::ILDDetID::FTD ;
      
      for( unsigned layer = 0; layer < nLayers; layer++ ){
        
        
        cellID[ lcio::ILDCellID0::layer  ] = layer ;
        
        unsigned nPetals = ftdLayers.getNPetals( layer );
        double deltaPhi = ( 2 * M_PI ) / nPetals;
        double phi0 = ftdLayers.getPhi0( layer );
        
        for( unsigned petal=0; petal< nPetals; petal++ ){
          
          
          cellID[ lcio::ILDCellID0::module ] = petal ;
          
          unsigned nSensors = 4; // FIXME: this should come from gear, also each sensor should get it's own CoordinateSystem!!!
          // TODO: and why are we starting at 1 anyway?
          
          for ( unsigned sensor = 1; sensor <=nSensors; sensor++ ){
            
            double stripAngle = 0.; // TODO: implement
            
            cellID[ lcio::ILDCellID0::side   ] = -1 ;                    
            cellID[ lcio::ILDCellID0::sensor ] = sensor ;
            int cellID0 = cellID.lowWord();
            
            // We start with the Translation Vector T (=origin of new CoordinateSystem )
            // the first petal (if phi0 = 0) sits with its symmetry axis on the x axis.
            // It starts at x= rInner and reaches rInner+ width on the x axis
            // So far the center of Mass will be on the x axis.
            // If we take into account that there are two sensors on each side, they all have their 
            // own centers of mass. We therefore just divide the big trapezoid into two smaller ones.
            double xmin = ftdLayers.getSensitiveLengthMin(layer);
            double xmax = ftdLayers.getSensitiveLengthMax(layer);
            
            // taking only half the petal:
            if( sensor%2 == 0 )xmax = (xmax+xmin)/2.;      // sensor 2 or 4, means the ones closer to the IP
            else               xmin = (xmax+xmin)/2.;      // sensor 1 or 3, the outer ones
            
            double x = (xmin+xmax) / 2.;
            double y = 0.;
            double z = -ftdLayers.getSensitiveZposition( layer, petal, sensor );
            CLHEP::Hep3Vector T( x , y, z);
            
            // Now we only have to rotate the petal around the z-axis into its place
            CLHEP::HepRotation rot;
            rot.rotateZ( petal * deltaPhi + phi0 );
            T = rot * T;
            
            
            //On to the rotation matrix
            CLHEP::HepRotation R;
            R.rotateZ( petal * deltaPhi + phi0 + stripAngle - M_PI/2. );
            
            
            CartesianCoordinateSystem* cartesian = new CartesianCoordinateSystem( T, R );
            MeasurementSurface* ms = new MeasurementSurface( cellID0, cartesian );
            addMeasurementSurface( ms );
            
            // Once more for the other side
            cellID[ lcio::ILDCellID0::side   ] = 1 ;   
            cellID0 = cellID.lowWord();
            
            T.setZ( -T.z() ); // switch 
                              // R is pretty much the same as the strip orientation will be the same,
                              // but as (by chosen definition) w should point towards the IP,
                              // we have to flip u around. So we acutally have to rotate 180° around v
                              // So first we get the vector v
            CLHEP::Hep3Vector v = R*CLHEP::Hep3Vector(0,1,0);
            // Then we rotate around it
            R.rotate( M_PI , v );
            
            CartesianCoordinateSystem* cartesian2 = new CartesianCoordinateSystem( T, R );
            MeasurementSurface* ms2 = new MeasurementSurface( cellID0, cartesian2 );
            addMeasurementSurface( ms2 );
            
            
            
            //TODO: do we need alpha in here?
            
            
          }
          
        }      
        
      }
      
    }
    
  } // end of GearExtensions namespace 
  
} // end of MarlinTrk namespace 











