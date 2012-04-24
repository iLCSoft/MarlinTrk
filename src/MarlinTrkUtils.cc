
#include "MarlinTrk/MarlinTrkUtils.h"

#include <vector>

#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/HelixTrack.h"

#include "lcio.h"
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <EVENT/TrackerHit.h>

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>

#include "streamlog/streamlog.h"

#define MIN_NDF 6

namespace MarlinTrk {
  
  
  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackImpl* track, bool fit_backwards, const EVENT::FloatVec& initial_cov_for_prefit, float bfield_z){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return -1 ;
    
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ")  ) ;
    }
    
    int return_error = 0;
    
    
    ///////////////////////////////////////////////////////
    // produce prefit parameters
    ///////////////////////////////////////////////////////
    
    IMPL::TrackStateImpl* pre_fit = new IMPL::TrackStateImpl();
    
    pre_fit->setCovMatrix(initial_cov_for_prefit);
    
    return_error = createPrefit(hit_list, pre_fit, bfield_z, fit_backwards);
    
    
    ///////////////////////////////////////////////////////
    // use prefit parameters to produce Finalised track
    ///////////////////////////////////////////////////////
    
    if( return_error == 0 ) {
      
      return_error = createFinalisedLCIOTrack( marlinTrk, hit_list, track, fit_backwards, pre_fit, bfield_z);
      
    } else {
      streamlog_out(DEBUG3) << "MarlinTrk::createFinalisedLCIOTrack : Prefit failed error = " << return_error << std::endl;
    }
    
    delete pre_fit;
    
    return return_error;
    
  }
  
  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackImpl* track, bool fit_backwards, IMPL::TrackStateImpl* pre_fit, float bfield_z){
    
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return -1 ;
    
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ")  ) ;
    }
    
    if( pre_fit == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackStateImpl == NULL ")  ) ;
    }
    
    
    int fit_status = createFit(hit_list, marlinTrk, pre_fit, bfield_z, fit_backwards);
    
    if( fit_status != IMarlinTrack::success ){ 
      
      streamlog_out(DEBUG3) << "MarlinTrk::createFinalisedLCIOTrack fit failed: fit_status = " << fit_status << std::endl; 
      
      return fit_status;
      
    } 
    
    int error = finaliseLCIOTrack(marlinTrk, track);
    
    
    return error;
    
  }
  
  
  
  
  int createFit( std::vector<EVENT::TrackerHit*>& hit_list, IMarlinTrack* marlinTrk, IMPL::TrackStateImpl* pre_fit, float bfield_z, bool fit_backwards ){
    
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return -1 ;
    
    if( marlinTrk == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::createFit: IMarlinTrack == NULL ")  ) ;
    }
    
    if( pre_fit == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::createFit: TrackStateImpl == NULL ")  ) ;
    }
    
    int return_error = 0;
    
    
    ///////////////////////////////////////////////////////
    // check that the prefit has the reference point at the correct location 
    ///////////////////////////////////////////////////////
    
    if (( fit_backwards == IMarlinTrack::backward && pre_fit->getLocation() != lcio::TrackState::AtLastHit ) 
        ||  
        ( fit_backwards == IMarlinTrack::forward && pre_fit->getLocation() != lcio::TrackState::AtFirstHit )) {            
      std::stringstream ss ;
      
      ss << "MarlinTrk::createFinalisedLCIOTrack track state must be set at either first or last hit. Location = ";
      ss << pre_fit->getLocation();
      
      throw EVENT::Exception( ss.str() );
      
    } 
    
    ///////////////////////////////////////////////////////
    // add hits to IMarlinTrk  
    ///////////////////////////////////////////////////////
    
    EVENT::TrackerHitVec::iterator it = hit_list.begin();
    
    //  start by trying to add the hits to the track we want to finally use. 
    streamlog_out(DEBUG2) << "MarlinTrk::createFit Start Fit: AddHits: number of hits to fit " << hit_list.size() << std::endl;
    
    EVENT::TrackerHitVec added_hits;
    unsigned int ndof_added = 0;
    
    for( it = hit_list.begin() ; it != hit_list.end() ; ++it ) {
      
      EVENT::TrackerHit* trkHit = *it;
      bool isSuccessful = false; 
      
      if( UTIL::BitSet32( trkHit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
        
        //Split it up and add both hits to the MarlinTrk
        const EVENT::LCObjectVec rawObjects = trkHit->getRawHits();                    
        
        for( unsigned k=0; k< rawObjects.size(); k++ ){
          
          EVENT::TrackerHit* rawHit = dynamic_cast< EVENT::TrackerHit* >( rawObjects[k] );
          
          if( marlinTrk->addHit( rawHit ) == IMarlinTrack::success ){
            
            isSuccessful = true; //if at least one hit from the spacepoint gets added
            ++ndof_added;
//            streamlog_out(DEBUG4) << "MarlinTrk::createFit ndof_added = " << ndof_added << std::endl;
          }
        }
      }
      else { // normal non composite hit
        
        if (marlinTrk->addHit( trkHit ) == IMarlinTrack::success ) {
          isSuccessful = true;
          ndof_added += 2;
//          streamlog_out(DEBUG4) << "MarlinTrk::createFit ndof_added = " << ndof_added << std::endl;          
        }
      }
      
      if (isSuccessful) {
        added_hits.push_back(trkHit);
      }        
      else{
        streamlog_out(DEBUG4) << "Hit " << it - hit_list.begin() << " Dropped " << std::endl;          
      }
      
    }
      
    if( ndof_added < MIN_NDF ) {
      streamlog_out(DEBUG4) << "MarlinTrk::createFit : Cannot fit less with less than " << MIN_NDF << " degrees of freedom. Number of hits =  " << added_hits.size() << " ndof = " << ndof_added << std::endl;
      return -1;
    }
      
    
    
    ///////////////////////////////////////////////////////
    // set the initial track parameters  
    ///////////////////////////////////////////////////////
    
    return_error = marlinTrk->initialise( *pre_fit, bfield_z, IMarlinTrack::backward ) ;
    if (return_error != IMarlinTrack::success) {
      
      streamlog_out(DEBUG5) << "MarlinTrk::createFit Initialisation of track fit failed with error : " << return_error << std::endl;
      
      return return_error;
      
    }
    
    ///////////////////////////////////////////////////////
    // try fit and return error
    ///////////////////////////////////////////////////////
    
    return marlinTrk->fit() ;
    
  }
  
  
  
  int createPrefit( std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackStateImpl* pre_fit, float bfield_z, bool fit_backwards){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return -1 ;
    
    if( pre_fit == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackStateImpl == NULL ")  ) ;
    }
    
    ///////////////////////////////////////////////////////
    // loop over all the hits and create a list consisting only 2D hits 
    ///////////////////////////////////////////////////////
    
    EVENT::TrackerHitVec twoD_hits;
    
    for (int ihit=0; ihit < hit_list.size(); ++ihit) {
      
      // check if this a space point or 2D hit 
      if(UTIL::BitSet32( hit_list[ihit]->getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] == false ){
        // then add to the list 
        twoD_hits.push_back(hit_list[ihit]);
        
      }
    }
    
    ///////////////////////////////////////////////////////
    // check that there are enough 2-D hits to create a helix 
    ///////////////////////////////////////////////////////
    
    if (twoD_hits.size() < 3) { // no chance to initialise print warning and return
      streamlog_out(WARNING) << "MarlinTrk::createFinalisedLCIOTrack Cannot create helix from less than 3 2-D hits" << std::endl;
      return -1;
    }
    
    ///////////////////////////////////////////////////////
    // make a helix from 3 hits to get a trackstate
    ///////////////////////////////////////////////////////
    
    // SJA:FIXME: this may not be the optimal 3 hits to take in certain cases where the 3 hits are not well spread over the track length 
    const double* x1 = twoD_hits[0]->getPosition();
    const double* x2 = twoD_hits[ twoD_hits.size()/2 ]->getPosition();
    const double* x3 = twoD_hits.back()->getPosition();
    
    HelixTrack helixTrack( x1, x2, x3, bfield_z, fit_backwards );
    
    helixTrack.moveRefPoint(hit_list.back()->getPosition()[0], hit_list.back()->getPosition()[1], hit_list.back()->getPosition()[2]);
    
    const float referencePoint[3] = { helixTrack.getRefPointX() , helixTrack.getRefPointY() , helixTrack.getRefPointZ() };
    
    pre_fit->setD0(helixTrack.getD0()) ;
    pre_fit->setPhi(helixTrack.getPhi0()) ;
    pre_fit->setOmega(helixTrack.getOmega()) ;
    pre_fit->setZ0(helixTrack.getZ0()) ;
    pre_fit->setTanLambda(helixTrack.getTanLambda()) ;
    
    pre_fit->setReferencePoint(referencePoint) ;
    
    if ( fit_backwards == IMarlinTrack::backward ) {
      pre_fit->setLocation(lcio::TrackState::AtLastHit);
    } else {
      pre_fit->setLocation(lcio::TrackState::AtFirstHit);
    }
    
    
    return 0;
    
  }
  
  int finaliseLCIOTrack( IMarlinTrack* marlintrk, IMPL::TrackImpl* track ){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( marlintrk == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: IMarlinTrack == NULL ")  ) ;
    }
    
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ")  ) ;
    }
    
    
    ///////////////////////////////////////////////////////
    // error to return if any
    ///////////////////////////////////////////////////////
    int return_error = 0;
    
    
    int ndf = 0;
    double chi2 = -DBL_MAX;
    
    ///////////////////////////////////////////////////////
    // first create trackstate at IP
    ///////////////////////////////////////////////////////
    const gear::Vector3D point(0.,0.,0.); // nominal IP
    
    IMPL::TrackStateImpl* trkStateIP = new IMPL::
    TrackStateImpl() ;
    
    
    ///////////////////////////////////////////////////////
    // make sure that the track state can be propagated to the IP and that the NDF is not less than 0
    ///////////////////////////////////////////////////////
    
    return_error = marlintrk->propagate(point, *trkStateIP, chi2, ndf ) ;
    
    if ( return_error != 0 || ndf < 0 ) { 
      streamlog_out(DEBUG4) << "MarlinTrk::finaliseLCIOTrack: return_code for propagation = " << return_error << " NDF = " << ndf << std::endl;
      delete trkStateIP;
      return return_error;
    }
    
    trkStateIP->setLocation(  lcio::TrackState::AtIP ) ;
    track->trackStates().push_back(trkStateIP);
    track->setChi2(chi2);
    track->setNdf(ndf);
    
    
    ///////////////////////////////////////////////////////
    // set the track states at the first and last hits 
    ///////////////////////////////////////////////////////
    
    std::vector<std::pair<EVENT::TrackerHit*, double> > hits_in_fit;
    std::vector<std::pair<EVENT::TrackerHit*, double> > outliers;
    
    marlintrk->getHitsInFit(hits_in_fit);
    marlintrk->getOutliers(outliers);
    
    for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      track->addHit(hits_in_fit[ihit].first);
    }

    for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      track->addHit(outliers[ihit].first);
    }
    
    ///////////////////////////////////////////////////////
    // first hit
    ///////////////////////////////////////////////////////
    
    IMPL::TrackStateImpl* trkStateAtFirstHit = new IMPL::TrackStateImpl() ;
    EVENT::TrackerHit* firstHit = hits_in_fit.front().first;
    
    
    return_error = marlintrk->getTrackState(firstHit, *trkStateAtFirstHit, chi2, ndf ) ;
    
    if ( return_error == 0 ) {
      trkStateAtFirstHit->setLocation(  lcio::TrackState::AtFirstHit ) ;
      track->trackStates().push_back(trkStateAtFirstHit);
    } else {
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at First Hit " << firstHit << std::endl ;
      delete trkStateAtFirstHit;
    }
    
    ///////////////////////////////////////////////////////
    // last hit
    ///////////////////////////////////////////////////////  
    
    IMPL::TrackStateImpl* trkStateAtLastHit = new IMPL::TrackStateImpl() ;
    EVENT::TrackerHit* lastHit = hits_in_fit.back().first;
    
    return_error = marlintrk->getTrackState(lastHit, *trkStateAtLastHit, chi2, ndf ) ;
    
    if ( return_error == 0 ) {
      trkStateAtLastHit->setLocation(  lcio::TrackState::AtLastHit ) ;
      track->trackStates().push_back(trkStateAtLastHit);
    } else {
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at Last Hit " << lastHit << std::endl ;
      delete trkStateAtLastHit;
    }
    
    
    ///////////////////////////////////////////////////////
    // set the track state at Calo Face 
    ///////////////////////////////////////////////////////
    
    IMPL::TrackStateImpl* trkStateCalo = new IMPL::TrackStateImpl;
    bool tanL_is_positive = trkStateIP->getTanLambda()>0 ;
    
    return_error = createTrackStateAtCaloFace(marlintrk, trkStateCalo, lastHit, tanL_is_positive);
    
    if ( return_error == 0 ) {
      trkStateCalo->setLocation(  lcio::TrackState::AtCalorimeter ) ;
      track->trackStates().push_back(trkStateCalo);
    } else {
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at Calo Face "  << std::endl ;
      delete trkStateCalo;
    }
    
    
    ///////////////////////////////////////////////////////
    // done
    ///////////////////////////////////////////////////////
    
    
    return return_error;
    
  }
  
  
  
  
  
  int createTrackStateAtCaloFace( IMarlinTrack* marlintrk, IMPL::TrackStateImpl* trkStateCalo, EVENT::TrackerHit* trkhit, bool tanL_is_positive){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( marlintrk == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::createTrackStateAtCaloFace: IMarlinTrack == NULL ")  ) ;
    }
    
    if( trkStateCalo == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::createTrackStateAtCaloFace: TrackImpl == NULL ")  ) ;
    }
    
    if( trkhit == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::createTrackStateAtCaloFace: TrackHit == NULL ")  ) ;
    }
    
    
    int return_error = 0;
    
    double chi2 = -DBL_MAX;
    int ndf = 0;
    
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
    encoder.reset() ;  // reset to 0
    
    encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::ECAL ;
    encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::barrel;
    encoder[lcio::ILDCellID0::layer]  = 0 ;
    
    int detElementID = 0;
    
    return_error = marlintrk->propagateToLayer(encoder.lowWord(), trkhit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
    
    if (return_error == IMarlinTrack::no_intersection ) { // try forward or backward
      if (tanL_is_positive) {
        encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::fwd;
      }
      else{
        encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::bwd;
      }
      return_error = marlintrk->propagateToLayer(encoder.lowWord(), trkhit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
    }
    
    if (return_error !=IMarlinTrack::success ) {
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> FinalRefit :  could not get TrackState at Calo Face: return_error = " << return_error << std::endl ;
    }
    
    
    
    
    return return_error;
    
    
  }
  
  void addHitsToTrack(IMPL::TrackImpl* track, std::vector<EVENT::TrackerHit*>& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::createTrackStateAtCaloFace: TrackImpl == NULL ")  ) ;
    }

    std::map<int, int> hitNumbers; 
    
    hitNumbers[lcio::ILDDetID::VXD] = 0;
    hitNumbers[lcio::ILDDetID::SIT] = 0;
    hitNumbers[lcio::ILDDetID::FTD] = 0;
    hitNumbers[lcio::ILDDetID::TPC] = 0;
    hitNumbers[lcio::ILDDetID::SET] = 0;
    hitNumbers[lcio::ILDDetID::ETD] = 0;
    
    for(unsigned int j=0; j<hit_list.size(); ++j) {
      
      track->addHit(hit_list.at(j)) ;
      
      cellID_encoder.setValue(hit_list.at(j)->getCellID0()) ;
      int detID = cellID_encoder[UTIL::ILDCellID0::subdet];
      ++hitNumbers[detID];
      //    streamlog_out( DEBUG1 ) << "Hit from Detector " << detID << std::endl;     
    }
    
    int offset = 2 ;
    if ( hits_in_fit == false ) { // all hit atributed by patrec
      offset = 1 ;
    }
    
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - offset ] = hitNumbers[lcio::ILDDetID::VXD];
    track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - offset ] = hitNumbers[lcio::ILDDetID::FTD];
    track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - offset ] = hitNumbers[lcio::ILDDetID::SIT];
    track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - offset ] = hitNumbers[lcio::ILDDetID::TPC];
    track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - offset ] = hitNumbers[lcio::ILDDetID::SET];
    track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - offset ] = hitNumbers[lcio::ILDDetID::ETD];

    
    
  }
  
}
