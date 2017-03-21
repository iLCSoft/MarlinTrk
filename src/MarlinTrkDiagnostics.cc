
#include "MarlinTrk/MarlinTrkDiagnostics.h"

#include "EVENT/LCObject.h"
#include "UTIL/BitSet32.h"
#include "UTIL/LCTrackerConf.h"

#ifdef MARLINTRK_DIAGNOSTICS_ON

namespace MarlinTrk{


  void getMCParticlesForTrackerHit(EVENT::TrackerHit* trkhit, std::vector<EVENT::MCParticle*>& mcps){

    if ( !trkhit ) {
      return;
    }
    
    // make sure there is nothing in the vector we wish to return
    mcps.clear();
    
    // first check if this is a composite space point
    if(UTIL::BitSet32( trkhit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]){
      
      const EVENT::LCObjectVec rawObjects = trkhit->getRawHits();

      for (unsigned iraw = 0; iraw < rawObjects.size(); ++iraw) {

        EVENT::TrackerHit* rawHit = dynamic_cast< EVENT::TrackerHit* >( rawObjects[iraw] );

        if( rawHit && rawHit->ext<MarlinTrk::MCTruth4HitExt>()){
          
          EVENT::MCParticle* mcp = rawHit->ext<MarlinTrk::MCTruth4HitExt>()->simhit->getMCParticle();
          bool found = false;
          // check that it is not already in the vector
          for (unsigned imcp=0; imcp<mcps.size(); ++imcp) {
            if (mcp == mcps[imcp]) {
              found = true;
              break;
            }
          }
          
          if( found == false ) mcps.push_back(mcp);
          
        }
        
        
      } // end of loop over rawObjects

      // end if COMPOSITE_SPACEPOINT
    } else {
     
      
      if( trkhit->ext<MarlinTrk::MCTruth4HitExt>()){
        mcps.push_back(trkhit->ext<MarlinTrk::MCTruth4HitExt>()->simhit->getMCParticle());
      }

      
    }
  }
}

#endif
