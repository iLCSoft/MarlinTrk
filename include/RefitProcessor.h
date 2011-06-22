#ifndef RefitProcessor_h
#define RefitProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <UTIL/LCRelationNavigator.h>

using namespace lcio ;
using namespace marlin ;

class MarlinKalTest ;


/**  Track Refitter processor for marlin. Refits an input track collection, producing a new collection of tracks
 * 

 *  <h4>Input - Prerequisites</h4>
 *  Needs a collection of LCIO Tracks. 
 *
 *  <h4>Output</h4> 
 *  Refitted LCIO Track Collection
 * 
 * @param InputTrackCollectionName Name of the Track collection to be refitted
 * @param OutputTrackCollectionName Name of the refitted Track collection to be refitted
 * 
 * @author S. J. Aplin, DESY
 */

class RefitProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new RefitProcessor ; }
  
  
  RefitProcessor() ;
  
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

  /* helper function to get collection using try catch block */
  LCCollection* GetCollection( LCEvent * evt, std::string colName ) ;
  
  /* helper function to get relations using try catch block */
  LCRelationNavigator* GetRelations(LCEvent * evt, std::string RelName ) ;

  /** Input track collection name for refitting.
   */
  std::string _input_track_col_name ;

 /** Input track relations name for refitting.
  */
  std::string _input_track_rel_name ;

  /** refitted track collection name.
   */
  std::string _output_track_col_name ;
  
 /** Output track relations name for refitting.
  */
  std::string _output_track_rel_name ;

  /** pointer to the MarlinKalTest instance which sets up the geometery
   */
  MarlinKalTest* _kaltest ;

  bool _MSOn ;
  bool _ElossOn ;

  int _n_run ;
  int _n_evt ;



} ;

#endif



