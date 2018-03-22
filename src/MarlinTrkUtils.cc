#include "MarlinTrk/MarlinTrkUtils.h"

#include <vector>
#include <algorithm>
#include <memory>
#include <math.h>

#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/Factory.h"

#include "MarlinTrk/MarlinDDKalTest.h"

#include "lcio.h"
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <EVENT/TrackerHit.h>

#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
#include <UTIL/Operators.h>

#include "streamlog/streamlog.h"

#include "TMatrixD.h"

#define MIN_NDF 6


namespace MarlinTrk {
  
  using namespace lcio ;
  using namespace UTIL ;

//  // Check if a square matrix is Positive Definite 
//  bool Matrix_Is_Positive_Definite(const EVENT::FloatVec& matrix){
//    
//    std::cout << "\n MarlinTrk::Matrix_Is_Positive_Definite(EVENT::FloatVec& matrix): " << std::endl;
//    
//    int icol,irow;
//    
//    int nrows = 5;
//    
//    TMatrixD cov(nrows,nrows) ; 
//    
//    bool matrix_is_positive_definite = true;
//    
//    int icov = 0;
//    for(int irow=0; irow<nrows; ++irow ){
//      for(int jcol=0; jcol<irow+1; ++jcol){
//        cov(irow,jcol) = matrix[icov];
//        cov(jcol,irow) = matrix[icov];
////        std::cout << " " << matrix[icov] ;
//        ++icov ;
//      }
////      std::cout << std::endl;
//    }
//    
////    cov.Print();
//    
//    double *pU = cov.GetMatrixArray();
//    
//    for (icol = 0; icol < nrows; icol++) {
//      const int rowOff = icol * nrows;
//      
//      //Compute fU(j,j) and test for non-positive-definiteness.
//      double ujj = pU[rowOff+icol];
//      double diagonal = ujj;
////      std::cout << "ERROR: diagonal = " << diagonal << std::endl;
//
//      for (irow = 0; irow < icol; irow++) {
//        const int pos_ij = irow*nrows+icol;
//        std::cout << " " << pU[pos_ij] ;
//        ujj -= pU[pos_ij]*pU[pos_ij];
//      }
//      std::cout  << " " << diagonal << std::endl;
//
//      
//      if (ujj <= 0) {
//        matrix_is_positive_definite = false;
//      }
//    }
//    
//      std::cout  << std::endl;
//    
//    if ( matrix_is_positive_definite == false ) {
//      std::cout << "******************************************************" << std::endl;
//      std::cout << "** ERROR:  matrix shown not to be positive definite **" << std::endl;
//      std::cout << "******************************************************" << std::endl;
//    }
//    
//    return matrix_is_positive_definite;
//    
//  }
  
  
  
  int createTrackStateAtCaloFace( IMarlinTrack* marlinTrk, IMPL::TrackStateImpl* track, EVENT::TrackerHit* trkhit, bool tanL_is_positive );
  
  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackImpl* track, bool fit_direction, const EVENT::FloatVec& initial_cov_for_prefit, float bfield_z, double maxChi2Increment){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return IMarlinTrack::bad_intputs ;
    
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ")  ) ;
    }
    
    int return_error = 0;
    
    
    ///////////////////////////////////////////////////////
    // produce prefit parameters
    ///////////////////////////////////////////////////////
    
    IMPL::TrackStateImpl pre_fit ;
    
    return_error = createPrefit(hit_list, &pre_fit, bfield_z, fit_direction);
    
    pre_fit.setCovMatrix(initial_cov_for_prefit);

    streamlog_out( DEBUG3 ) << " **** createFinalisedLCIOTrack - created pre-fit: " <<  toString( &pre_fit )  << std::endl ;


    ///////////////////////////////////////////////////////
    // use prefit parameters to produce Finalised track
    ///////////////////////////////////////////////////////
    
    if( return_error == 0 ) {
      
      return_error = createFinalisedLCIOTrack( marlinTrk, hit_list, track, fit_direction, &pre_fit, bfield_z, maxChi2Increment);
      
    } else {
      streamlog_out(DEBUG3) << "MarlinTrk::createFinalisedLCIOTrack : Prefit failed error = " << return_error << std::endl;
    }
    
    
    return return_error;
    
  }
  
  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackImpl* track, bool fit_direction, EVENT::TrackState* pre_fit, float bfield_z, double maxChi2Increment){
    
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return IMarlinTrack::bad_intputs ;
    
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ")  ) ;
    }
    
    // if( pre_fit == 0 ){
    //   throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackStateImpl == NULL ")  ) ;
    // }
    
    
    int fit_status = createFit(hit_list, marlinTrk, pre_fit, bfield_z, fit_direction, maxChi2Increment);
    
    if( fit_status != IMarlinTrack::success ){ 
      
      streamlog_out(DEBUG3) << "MarlinTrk::createFinalisedLCIOTrack fit failed: fit_status = " << fit_status << std::endl; 
      
      return fit_status;
      
    } 
    
    int error = finaliseLCIOTrack(marlinTrk, track, hit_list, fit_direction );
    
    
    return error;
    
  }
  
  
  
  
  int createFit( std::vector<EVENT::TrackerHit*>& hit_list, IMarlinTrack* marlinTrk, EVENT::TrackState* pre_fit, float bfield_z, bool fit_direction, double maxChi2Increment){
    
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return IMarlinTrack::bad_intputs;
    
    if( marlinTrk == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::createFit: IMarlinTrack == NULL ")  ) ;
    }
    
    // if( pre_fit == 0 ){
    //   throw EVENT::Exception( std::string("MarlinTrk::createFit: TrackStateImpl == NULL ")  ) ;
    // }
    
    int return_error = 0;
    
    
//    ///////////////////////////////////////////////////////
//    // check that the prefit has the reference point at the correct location 
//    ///////////////////////////////////////////////////////
//    
//    if (( fit_direction == IMarlinTrack::backward && pre_fit->getLocation() != lcio::TrackState::AtLastHit ) 
//        ||  
//        ( fit_direction == IMarlinTrack::forward && pre_fit->getLocation() != lcio::TrackState::AtFirstHit )) {            
//      std::stringstream ss ;
//      
//      ss << "MarlinTrk::createFinalisedLCIOTrack track state must be set at either first or last hit. Location = ";
//      ss << pre_fit->getLocation();
//      
//      throw EVENT::Exception( ss.str() );
//      
//    } 
    
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
            streamlog_out(DEBUG4) << "MarlinTrk::createFit ndof_added = " << ndof_added << std::endl;
          }
        }
      }
      else { // normal non composite hit
        
        if (marlinTrk->addHit( trkHit ) == IMarlinTrack::success ) {
          isSuccessful = true;
          ndof_added += 2;
	  streamlog_out(DEBUG4) << "MarlinTrk::createFit ndof_added = " << ndof_added << std::endl;          
        }
      }
      
      if (isSuccessful) {
        added_hits.push_back(trkHit);
      }        
      else{
        streamlog_out(DEBUG2) << "Hit " << it - hit_list.begin() << " Dropped " << std::endl;
      }
      
    }
      
    if( ndof_added < MIN_NDF ) {
      streamlog_out(DEBUG2) << "MarlinTrk::createFit : Cannot fit less with less than " << MIN_NDF << " degrees of freedom. Number of hits =  " << added_hits.size() << " ndof = " << ndof_added << std::endl;
      return IMarlinTrack::bad_intputs;
    }
      
    
    
    ///////////////////////////////////////////////////////
    // set the initial track parameters  
    ///////////////////////////////////////////////////////
    

    if( pre_fit == 0 ) {

      streamlog_out(DEBUG5) << "MarlinTrk::createFit : null pointer for pre_fit given - will fall back "
			    << " to default initialisation ..." << std::endl ;

      return_error = marlinTrk->initialise( fit_direction ) ; 

    } else {

      return_error = marlinTrk->initialise( *pre_fit, bfield_z, fit_direction ) ;//IMarlinTrack::backward ) ;
    }


    if (return_error != IMarlinTrack::success) {
      
      streamlog_out(DEBUG5) << "MarlinTrk::createFit Initialisation of track fit failed with error : " << return_error << std::endl;
      
      return return_error;
      
    }
    
    
#if 0 // DEBUG code:

    IMPL::TrackStateImpl ts ; 
    double chi2(0) ;
    int ndf(0) ;
    int ii = marlinTrk->propagate( Vector3D(),  ts, chi2, ndf ) ;
    
    streamlog_out(DEBUG5) << " MarlinTrk::createFit:   pre-fit, propagated to the IP : " << ts << std::endl ;

    
    return IMarlinTrack::success ;
#endif


    ///////////////////////////////////////////////////////
    // try fit and return error
    ///////////////////////////////////////////////////////
    
    return marlinTrk->fit(maxChi2Increment) ;
    
  }
  
  
  
  int createPrefit( std::vector<EVENT::TrackerHit*>& hit_list, IMPL::TrackStateImpl* pre_fit, float bfield_z, bool fit_direction){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return IMarlinTrack::bad_intputs ;
    
    if( pre_fit == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackStateImpl == NULL ")  ) ;
    }
    
    ///////////////////////////////////////////////////////
    // loop over all the hits and create a list consisting only 2D hits 
    ///////////////////////////////////////////////////////
    
    EVENT::TrackerHitVec twoD_hits;
    
    for (unsigned ihit=0; ihit < hit_list.size(); ++ihit) {
      
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
      return IMarlinTrack::bad_intputs;
    }
    
    ///////////////////////////////////////////////////////
    // make a helix from 3 hits to get a trackstate
    ///////////////////////////////////////////////////////
    
    // SJA:FIXME: this may not be the optimal 3 hits to take in certain cases where the 3 hits are not well spread over the track length 
    const double* x1 = twoD_hits[0]->getPosition();
    const double* x2 = twoD_hits[ twoD_hits.size()/2 ]->getPosition();
    const double* x3 = twoD_hits.back()->getPosition();
    
    HelixTrack helixTrack( x1, x2, x3, bfield_z, HelixTrack::forwards );

    helixTrack.moveRefPoint(0.0, 0.0, 0.0);

    // if ( fit_direction == IMarlinTrack::backward ) {
    //   pre_fit->setLocation(lcio::TrackState::AtLastHit);
    //   helixTrack.moveRefPoint(hit_list.back()->getPosition()[0], hit_list.back()->getPosition()[1], hit_list.back()->getPosition()[2]);
    // } else {
    //   pre_fit->setLocation(lcio::TrackState::AtFirstHit);
    //   helixTrack.moveRefPoint(hit_list.front()->getPosition()[0], hit_list.front()->getPosition()[1], hit_list.front()->getPosition()[2]);
    // }

    const float referencePoint[3] = { float(helixTrack.getRefPointX()) ,  float(helixTrack.getRefPointY()) ,  float(helixTrack.getRefPointZ() )};
    
    pre_fit->setD0(helixTrack.getD0()) ;
    pre_fit->setPhi(helixTrack.getPhi0()) ;
    pre_fit->setOmega(helixTrack.getOmega()) ;
    pre_fit->setZ0(helixTrack.getZ0()) ;
    pre_fit->setTanLambda(helixTrack.getTanLambda()) ;
    
    pre_fit->setReferencePoint(referencePoint) ;
    
    return IMarlinTrack::success;
    
  }
  
  int finaliseLCIOTrack( IMarlinTrack* marlintrk, IMPL::TrackImpl* track, std::vector<EVENT::TrackerHit*>& hit_list, bool fit_direction, IMPL::TrackStateImpl* atLastHit, IMPL::TrackStateImpl* atCaloFace){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( marlintrk == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: IMarlinTrack == NULL ")  ) ;
    }
    
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ")  ) ;
    }
    
    if( atCaloFace && atLastHit == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: atLastHit == NULL ")  ) ;
    }

    if( atLastHit && atCaloFace == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: atCaloFace == NULL ")  ) ;
    }

    
    
    ///////////////////////////////////////////////////////
    // error to return if any
    ///////////////////////////////////////////////////////
    int return_error = 0;
    
    int ndf = 0;
    double chi2 = -DBL_MAX;
    
    /////////////////////////////////////////////////////////////
    // First check NDF to see if it make any sense to continue.
    // The track will be dropped if the NDF is less than 0 
    /////////////////////////////////////////////////////////////
    
    return_error = marlintrk->getNDF(ndf);

    if ( return_error != IMarlinTrack::success) {
      streamlog_out(DEBUG3) << "MarlinTrk::finaliseLCIOTrack: getNDF returns " << return_error << std::endl;
      return return_error;
    } else if( ndf < 0 ) {
      streamlog_out(DEBUG8) << "MarlinTrk::finaliseLCIOTrack: number of degrees of freedom less than 0 track dropped : NDF = " << ndf << std::endl;
      return IMarlinTrack::error;
    } else {
      streamlog_out(DEBUG4) << "MarlinTrk::finaliseLCIOTrack: NDF = " << ndf << std::endl;
    }
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // get the list of hits used in the fit
    // add these to the track, add spacepoints as long as at least on strip hit is used.  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::vector<std::pair<EVENT::TrackerHit*, double> > hits_in_fit;
    std::vector<std::pair<EVENT::TrackerHit*, double> > outliers;
    std::vector<EVENT::TrackerHit*> used_hits;
        
    hits_in_fit.reserve(300);
    outliers.reserve(300);
    
    marlintrk->getHitsInFit(hits_in_fit);
    marlintrk->getOutliers(outliers);
    
    ///////////////////////////////////////////////
    // now loop over the hits provided for fitting 
    // we do this so that the hits are added in the
    // order in which they have been fitted
    ///////////////////////////////////////////////
    
    for ( unsigned ihit = 0; ihit < hit_list.size(); ++ihit) {
      
      EVENT::TrackerHit* trkHit = hit_list[ihit];
      
      if( UTIL::BitSet32( trkHit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
        
        // get strip hits 
        const EVENT::LCObjectVec rawObjects = trkHit->getRawHits();                    
        
        for( unsigned k=0; k< rawObjects.size(); k++ ){
          
          EVENT::TrackerHit* rawHit = dynamic_cast< EVENT::TrackerHit* >( rawObjects[k] );
          
          bool is_outlier = false;
          
          // here we loop over outliers as this will be faster than looping over the used hits
          for ( unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
            
            if ( rawHit == outliers[ohit].first ) { 
              is_outlier = true;                                    
              break; // break out of loop over outliers
            }
          }
          
          if (is_outlier == false) {
            used_hits.push_back(hit_list[ihit]);
            track->addHit(used_hits.back());
            break; // break out of loop over rawObjects
          }          
        }
      } else {
        
        bool is_outlier = false;
        
        // here we loop over outliers as this will be faster than looping over the used hits
        for ( unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
          
          if ( trkHit == outliers[ohit].first ) { 
            is_outlier = true;                                    
            break; // break out of loop over outliers
          }
        }
        
        if (is_outlier == false) {
          used_hits.push_back(hit_list[ihit]);
          track->addHit(used_hits.back());
        }          
        
        
      }
      
    }

    
//    ///////////////////////////////////////////////////////////////////////////
//    // We now need to find out at which point the fit is constrained 
//    // and therefore be able to provide well formed (pos. def.) cov. matrices
//    ///////////////////////////////////////////////////////////////////////////
//    
    
        
    ///////////////////////////////////////////////////////
    // first hit
    ///////////////////////////////////////////////////////
    
    IMPL::TrackStateImpl* trkStateAtFirstHit = new IMPL::TrackStateImpl() ;
    EVENT::TrackerHit* firstHit = ( fit_direction == IMarlinTrack::backward ? hits_in_fit.back().first : hits_in_fit.front().first ) ;

     ///////////////////////////////////////////////////////
    // last hit
    ///////////////////////////////////////////////////////
    
    EVENT::TrackerHit* lastHit =  ( fit_direction == IMarlinTrack::backward ? hits_in_fit.front().first :  hits_in_fit.back().first ) ;
          
    EVENT::TrackerHit* last_constrained_hit = 0 ;     
    marlintrk->getTrackerHitAtPositiveNDF(last_constrained_hit);


    streamlog_out(DEBUG3) << "MarlinTrk::finaliseLCIOTrack: firstHit : " << toString( firstHit ) 
			  << " lastHit:                                " << toString( lastHit )
			  << " last constrained hit:                   " << toString( last_constrained_hit )
			  << " fit direction is forward : " << fit_direction << std::endl ;

    //fgx    return_error = marlintrk->smooth(lastHit);
    return_error = marlintrk->smooth( last_constrained_hit );
    
    streamlog_out(DEBUG4) << "MarlinTrk::finaliseLCIOTrack: return_code for smoothing to last constrained hit " 
			  << last_constrained_hit << " = " << return_error << " NDF = " << ndf << std::endl;

    if ( return_error != IMarlinTrack::success ) { 
      delete trkStateAtFirstHit;
      //      delete trkStateAtLastHit;
      return return_error ;
    }

    
    ///////////////////////////////////////////////////////
    // first create trackstate at IP
    ///////////////////////////////////////////////////////
    const Vector3D point(0.,0.,0.); // nominal IP
    
    IMPL::TrackStateImpl* trkStateIP = new IMPL::
    TrackStateImpl() ;
      
    

    streamlog_out(DEBUG4) << "MarlinTrk::finaliseLCIOTrack: finalised kaltest track  : "
			  << marlintrk->toString() << std::endl ;



    ///////////////////////////////////////////////////////
    // make sure that the track state can be propagated to the IP 
    ///////////////////////////////////////////////////////
    
    MarlinTrk::IMarlinTrkSystem* trksystem =  MarlinTrk::Factory::getCurrentMarlinTrkSystem() ;

    bool usingAidaTT = ( trksystem->name() == "AidaTT" ) ;

    // if we fitted backwards, the firstHit is the last one used in the fit and we simply propagate to the IP:

    if(  fit_direction == IMarlinTrack::backward ||  usingAidaTT ) {   

      return_error = marlintrk->propagate(point, firstHit, *trkStateIP, chi2, ndf ) ;

    } else { 

      // if we fitted forward, we start from the last_constrained hit
      // and then add the last inner hits with a Kalman step ... 
      
      // create a temporary IMarlinTrack 
      
      auto mTrk = std::shared_ptr<MarlinTrk::IMarlinTrack>( trksystem->createTrack() ) ;
      
      IMPL::TrackStateImpl ts;

      double chi2Tmp = 0 ;
      int ndfTmp = 0 ;
      return_error = marlintrk->getTrackState( last_constrained_hit, ts ,  chi2 , ndf ) ;
      
      streamlog_out( DEBUG3  )  << "  MarlinTrk::finaliseLCIOTrack:--  TrackState at last constrained hit : " << std::endl
				<< toString( &ts )    << std::endl ;
      
      //need to add a dummy hit to the track
      mTrk->addHit( last_constrained_hit ) ;
      
      double _bfield = 42.0 ;   
      // fixme: the implementation for DDKalTest does no longer need this value but the IMarlinTrk interface is not yet changed
      mTrk->initialise( ts ,  _bfield ,  fit_direction ) ;

      std::vector<std::pair<EVENT::TrackerHit*, double> >::reverse_iterator hI = hits_in_fit.rbegin() ;

      while( (*hI).first != last_constrained_hit ){

	streamlog_out( DEBUG0  )  << "  MarlinTrk::finaliseLCIOTrack:--  hit in reverse_iterator : "  << std::endl
				  << toString( (*hI).first ) << std::endl ; 
	++hI ;
      }

      ++hI ;

      while( hI != hits_in_fit.rend()  ){

	EVENT::TrackerHit* h =  (*hI).first ;


	double deltaChi ;
	double maxChi2Increment = 1e10 ; // ???

	int addHit = mTrk->addAndFit(  h , deltaChi, maxChi2Increment ) ; 


	streamlog_out( DEBUG3 ) << " MarlinTrk::finaliseLCIOTrack: hit " << toString( h )   
				<< "  added : " << MarlinTrk::errorCode( addHit )
				<< "  deltaChi2: " << deltaChi 
				<< std::endl ;

	if( addHit !=  MarlinTrk::IMarlinTrack::success ){

	  streamlog_out( ERROR ) << " ****  MarlinTrk::finaliseLCIOTrack:  could not add inner hit to track !!! " << std::endl ; 
	}

	++hI ;

      }//------------------------------------

      streamlog_out(DEBUG4) << "MarlinTrk::finaliseLCIOTrack: temporary kaltest track for track state at the IP: "
			    <<  mTrk->toString() << std::endl ;

      // now propagate the temporary track to the IP
      return_error = mTrk->propagate( point, firstHit, *trkStateIP, chi2Tmp, ndfTmp ) ;

      
      streamlog_out( DEBUG4 ) << " ***  MarlinTrk::finaliseLCIOTrack: - propagated temporary track fromfirst hit to IP : " <<  toString( trkStateIP ) << std::endl ;
    }


    if ( return_error != IMarlinTrack::success ) { 
      streamlog_out(DEBUG4) << "MarlinTrk::finaliseLCIOTrack: return_code for propagation = " << return_error << " NDF = " << ndf << std::endl;
      delete trkStateIP;
      delete trkStateAtFirstHit;
      //      delete trkStateAtLastHit;

      return return_error ;
    }
    
    trkStateIP->setLocation(  lcio::TrackState::AtIP ) ;
    track->trackStates().push_back(trkStateIP);
    track->setChi2(chi2);
    track->setNdf(ndf);
    
              
    ///////////////////////////////////////////////////////
    // set the track states at the first and last hits 
    ///////////////////////////////////////////////////////    
    
    ///////////////////////////////////////////////////////
    // @ first hit
    ///////////////////////////////////////////////////////
    
    streamlog_out( DEBUG5 ) << "  >>>>>>>>>>>MarlinTrk::finaliseLCIOTrack:  create TrackState AtFirstHit" << std::endl ;

    
    return_error = marlintrk->getTrackState(firstHit, *trkStateAtFirstHit, chi2, ndf ) ;
    
    if ( return_error == IMarlinTrack::success ) {
      trkStateAtFirstHit->setLocation(  lcio::TrackState::AtFirstHit ) ;
      track->trackStates().push_back(trkStateAtFirstHit);
    } else {
      streamlog_out( WARNING ) << "  >>>>>>>>>>>MarlinTrk::finaliseLCIOTrack:  MarlinTrk::finaliseLCIOTrack:  could not get TrackState at First Hit " << firstHit << std::endl ;
      delete trkStateAtFirstHit;
    }
    
    double r_first = firstHit->getPosition()[0]*firstHit->getPosition()[0] + firstHit->getPosition()[1]*firstHit->getPosition()[1];
    
    track->setRadiusOfInnermostHit(sqrt(r_first));
    
    if ( atLastHit == 0 && atCaloFace == 0 ) {
    
      ///////////////////////////////////////////////////////
      // @ last hit
      ///////////////////////////////////////////////////////  
      
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack: create TrackState AtLastHit : using trkhit " << last_constrained_hit << std::endl ;
      
      Vector3D last_hit_pos(lastHit->getPosition());
      
      IMPL::TrackStateImpl* trkStateAtLastHit = new IMPL::TrackStateImpl() ;

      return_error = marlintrk->propagate(last_hit_pos, last_constrained_hit, *trkStateAtLastHit, chi2, ndf);
            
//      return_error = marlintrk->getTrackState(lastHit, *trkStateAtLastHit, chi2, ndf ) ;
      
      if ( return_error == IMarlinTrack::success ) {
        trkStateAtLastHit->setLocation(  lcio::TrackState::AtLastHit ) ;
        track->trackStates().push_back(trkStateAtLastHit);
      } else {
        streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at Last Hit " << last_constrained_hit << std::endl ;
        delete trkStateAtLastHit;
      }
      
//      const EVENT::FloatVec& ma = trkStateAtLastHit->getCovMatrix();
//      
//      Matrix_Is_Positive_Definite( ma );



      ///////////////////////////////////////////////////////
      // set the track state at Calo Face 
      ///////////////////////////////////////////////////////
      
      IMPL::TrackStateImpl* trkStateCalo = new IMPL::TrackStateImpl;
      bool tanL_is_positive = trkStateIP->getTanLambda()>0 ;
      
      return_error = createTrackStateAtCaloFace(marlintrk, trkStateCalo, last_constrained_hit, tanL_is_positive);
      //      return_error = createTrackStateAtCaloFace(marlintrk, trkStateCalo, lastHit, tanL_is_positive);
      
      if ( return_error == IMarlinTrack::success ) {
        trkStateCalo->setLocation(  lcio::TrackState::AtCalorimeter ) ;
        track->trackStates().push_back(trkStateCalo);
      } else {
        streamlog_out( DEBUG9 ) << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at Calo Face "  << std::endl ;
        delete trkStateCalo;

	//FIXME: ignore track state at Calo face for debugging new tracking ...
#if 0
	return_error = IMarlinTrack::success ;
        streamlog_out( DEBUG9 ) << "     MarlinTrk::finaliseLCIOTrack:            ignore missing TrackState at Calo Face  for debugging " << std::endl ;
#endif    

      }

      
    } else {
      track->trackStates().push_back(atLastHit);
      track->trackStates().push_back(atCaloFace);
    }
    
    
    
    ///////////////////////////////////////////////////////
    // done
    ///////////////////////////////////////////////////////
    
    
    return return_error;
    
  }
  
  
  
  
  
  int createTrackStateAtCaloFace( IMarlinTrack* marlintrk, IMPL::TrackStateImpl* trkStateCalo, EVENT::TrackerHit* trkhit, bool tanL_is_positive ){
    
    streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> createTrackStateAtCaloFace : using trkhit " 
			    << UTIL::toString( trkhit ) << " tanL_is_positive = " << tanL_is_positive << std::endl ;
    
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
    
    UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ; 
    encoder.reset() ;  // reset to 0
    
    // ================== need to get the correct ID(s) for the calorimeter face  ============================

    unsigned ecal_barrel_face_ID = lcio::ILDDetID::ECAL ;
    unsigned ecal_endcap_face_ID = lcio::ILDDetID::ECAL_ENDCAP ;

    //=========================================================================================================

    encoder[lcio::LCTrackerCellID::subdet()] = ecal_barrel_face_ID ;
    encoder[lcio::LCTrackerCellID::side()]   = lcio::ILDDetID::barrel;
    encoder[lcio::LCTrackerCellID::layer()]  = 0 ;
    
    int detElementID = 0;
    
    return_error = marlintrk->propagateToLayer(encoder.lowWord(), trkhit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
    
    if (return_error == IMarlinTrack::no_intersection ) { // try forward or backward

      encoder[lcio::LCTrackerCellID::subdet()] = ecal_endcap_face_ID ;

      if (tanL_is_positive) {
        encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::fwd;
      }
      else{
        encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::bwd;
      }
      return_error = marlintrk->propagateToLayer(encoder.lowWord(), trkhit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
    }
    
    //fg: for curling tracks the propagated track has the wrong z0 whereas it should be 0. really 
    if( std::abs( trkStateCalo->getZ0() ) > std::abs( 2.*M_PI/trkStateCalo->getOmega() * trkStateCalo->getTanLambda() ) ){

      streamlog_out( DEBUG2 ) << "  >>>>>>>>>>> createTrackStateAtCaloFace : setting z0 to 0. for track state at calorimeter : " 
			       << toString(trkStateCalo ) << std::endl ;

      trkStateCalo->setZ0( 0. ) ;
    } 

    if (return_error !=IMarlinTrack::success ) {
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> createTrackStateAtCaloFace :  could not get TrackState at Calo Face: return_error = " << return_error << std::endl ;
    }
    
    
    return return_error;
    
    
  }
  
  void addHitNumbersToTrack(IMPL::TrackImpl* track, std::vector<EVENT::TrackerHit*>& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::addHitsToTrack: TrackImpl == NULL ")  ) ;
    }

    std::map<int, int> hitNumbers; 
    
    for(unsigned int j=0; j<hit_list.size(); ++j) {
      
      cellID_encoder.setValue(hit_list.at(j)->getCellID0()) ;
      int detID = cellID_encoder[UTIL::LCTrackerCellID::subdet()];
      ++hitNumbers[detID];
    }
    
    int offset = 2 ;
    if ( hits_in_fit == false ) { // all hit atributed by patrec
      offset = 1 ;
    }
    
    // this assumes that there is no tracker with an index larger than the ecal ...
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ECAL);


    for(  std::map<int, int>::iterator it = hitNumbers.begin() ; 
	  it != hitNumbers.end() ; ++it ){

      int detIndex = it->first ;
      track->subdetectorHitNumbers().at( 2 * detIndex - offset ) = it->second ;

    }

  }
  
  void addHitNumbersToTrack(IMPL::TrackImpl* track, std::vector<std::pair<EVENT::TrackerHit* , double> >& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::addHitsToTrack: TrackImpl == NULL ")  ) ;
    }
    
    std::map<int, int> hitNumbers; 
    
    for(unsigned int j=0; j<hit_list.size(); ++j) {
      
      cellID_encoder.setValue(hit_list.at(j).first->getCellID0()) ;
      int detID = cellID_encoder[UTIL::LCTrackerCellID::subdet()];
      ++hitNumbers[detID];
      //    streamlog_out( DEBUG1 ) << "Hit from Detector " << detID << std::endl;     
    }
    
    int offset = 2 ;
    if ( hits_in_fit == false ) { // all hit atributed by patrec
      offset = 1 ;
    }
    

    // this assumes that there is no tracker with an index larger than the ecal ...
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ECAL);
    
    
    for(  std::map<int, int>::iterator it = hitNumbers.begin() ; 
	  it != hitNumbers.end() ; ++it ){
      
      int detIndex = it->first ;
      track->subdetectorHitNumbers().at( 2 * detIndex - offset ) = it->second ;
    }
    
  }
  
}
