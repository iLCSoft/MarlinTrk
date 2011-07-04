#ifndef TrackerHitLayerNumberProcessor_h
#define TrackerHitLayerNumberProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>


using namespace lcio ;
using namespace marlin ;


/** ======= TrackerHitLayerNumberProcessor ========== <br>
 * Calculates Layer Numbers for TrackerHit collections <br> 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collections of TrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of digitized TrackerHits in the vertex detector and SIT <br>
 * @param VTXHitCollection The name of output collection of digitized VTX TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param SITHitCollection The name of output collection of digitized SIT TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param TPCHitCollection The name of output collection of digitized TPC TrackerHits <br>
 * (default name TPCTrackerHits) <br>
 * @param FTDHitCollection The name of output collection of digitized FTD TrackerHits <br>
 * (default name TPCTrackerHits) <br>
 * @param SETHitCollection The name of output collection of digitized SET TrackerHits <br>
 * (default name ETDTrackerHits) <br>
 * @param ETDHitCollection The name of output collection of digitized ETD TrackerHits <br>
 * (default name ETDTrackerHits) <br>
 * <br>
 * 
 */
class TrackerHitLayerNumberProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackerHitLayerNumberProcessor ; }
  
  
  TrackerHitLayerNumberProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

  
 protected:

  std::string _colNameVTX ;
  std::string _colNameSIT ;
  std::string _colNameTPC ;
  std::string _colNameFTD ;
  std::string _colNameSET ;
  std::string _colNameETD ;

  int _nRun ;
  int _nEvt ;


} ;

#endif



