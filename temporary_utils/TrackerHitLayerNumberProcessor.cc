/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "TrackerHitLayerNumberProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/MCParticle.h>

#include <gsl/gsl_randist.h>
#include "marlin/ProcessorEventSeeder.h"

#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>

#include <ILDDetectorIDs.h>

#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>
#include <sstream>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

TrackerHitLayerNumberProcessor aTrackerHitLayerNumberProcessor ;


TrackerHitLayerNumberProcessor::TrackerHitLayerNumberProcessor() : Processor("TrackerHitLayerNumberProcessor") {
  
  // modify processor description
  _description = "TrackerHitLayerNumberProcessor should create layer Numbers for TrackerHits" ;
  

  // register steering parameters: name, description, class-variable, default value

  
  // Input collections
  registerInputCollection( LCIO::TRACKERHIT,
                            "VTXHitCollection" , 
                            "Name of the vxd TrackerHit collection"  ,
                            _colNameVTX ,
                            std::string("VTXTrackerHits") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                            "SITHitCollection" , 
                            "Name of the sit TrackerHit collection"  ,
                            _colNameSIT ,
                            std::string("SITTrackerHits") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                            "TPCHitCollection" , 
                            "Name of the tpc TrackerHit collection"  ,
                            _colNameTPC ,
                            std::string("TPCTrackerHits") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                            "FTDHitCollection" , 
                            "Name of the ftd TrackerHit collection"  ,
                            _colNameFTD ,
                            std::string("FTDTrackerHits") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                            "SETHitCollection" , 
                            "Name of the set TrackerHit collection"  ,
                            _colNameSET ,
                            std::string("SETTrackerHits") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                            "ETDHitCollection" , 
                            "Name of the etd TrackerHit collection"  ,
                            _colNameETD ,
                            std::string("ETDTrackerHits") ) ;
  
  
}


void TrackerHitLayerNumberProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

}

void TrackerHitLayerNumberProcessor::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void TrackerHitLayerNumberProcessor::processEvent( LCEvent * evt ) { 
	
  for (int iColl=0;iColl<6;++iColl) {

    LCCollection* THcol = 0 ;
    try{
      if (iColl==0)
        THcol = evt->getCollection( _colNameVTX ) ;
      else if (iColl==1)
        THcol = evt->getCollection( _colNameSIT ) ;
      else if (iColl==2)
        THcol = evt->getCollection( _colNameTPC ) ;
      else if (iColl==3)
        THcol = evt->getCollection( _colNameFTD ) ;
      else if (iColl==4)
        THcol = evt->getCollection( _colNameSET ) ;
      else if (iColl==5)
        THcol = evt->getCollection( _colNameETD ) ;

    }
    catch(DataNotAvailableException &e){

      if (iColl==0)
        streamlog_out(DEBUG) << "Collection " << _colNameVTX.c_str() << " is unavailable in event " << _nEvt << std::endl;
      else if (iColl==1)
        streamlog_out(DEBUG) << "Collection " << _colNameSIT.c_str() << " is unavailable in event " << _nEvt << std::endl;
      else if (iColl==2)
        streamlog_out(DEBUG) << "Collection " << _colNameTPC.c_str() << " is unavailable in event " << _nEvt << std::endl;
      else if (iColl==3)
        streamlog_out(DEBUG) << "Collection " << _colNameFTD.c_str() << " is unavailable in event " << _nEvt << std::endl;
      else if (iColl==4)
        streamlog_out(DEBUG) << "Collection " << _colNameSET.c_str() << " is unavailable in event " << _nEvt << std::endl;
      else if (iColl==5)
        streamlog_out(DEBUG) << "Collection " << _colNameETD.c_str() << " is unavailable in event " << _nEvt << std::endl;
    }

    if( THcol != 0 ){    
    
      int nSimHits = THcol->getNumberOfElements()  ;

      if (iColl==0) {
        //VXD
        streamlog_out( DEBUG ) << " processing collection " << _colNameVTX 
                               << " with " <<  nSimHits  << " hits ... " << std::endl ;
      }            
      else if (iColl==1) {
        //SIT
        streamlog_out( DEBUG ) << " processing collection " << _colNameSIT 
                               << " with " <<  nSimHits  << " hits ... " << std::endl ;
      }
      else if (iColl==2) {
        //TPC
        streamlog_out( DEBUG ) << " processing collection " << _colNameTPC 
                               << " with " <<  nSimHits  << " hits ... " << std::endl ;
      }
      else if (iColl==3) {
        //TPC
        streamlog_out( DEBUG ) << " processing collection " << _colNameFTD 
                               << " with " <<  nSimHits  << " hits ... " << std::endl ;
      }
      else if (iColl==4) {
        //TPC
        streamlog_out( DEBUG ) << " processing collection " << _colNameSET 
                               << " with " <<  nSimHits  << " hits ... " << std::endl ;
      }
      else if (iColl==5) {
        //TPC
        streamlog_out( DEBUG ) << " processing collection " << _colNameETD 
                               << " with " <<  nSimHits  << " hits ... " << std::endl ;
      }
    
      for(int i=0; i< nSimHits; ++i){
        
        TrackerHit* trkHit = dynamic_cast<TrackerHit*>( THcol->getElementAt( i ) ) ;
                
        
        const double *pos ;
        pos =  trkHit->getPosition() ;  
        
        gear::Vector3D hitvec(pos[0],pos[1],pos[2]);
        
        streamlog_out(DEBUG) <<"Position of hit = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
        
        CLHEP::Hep3Vector point(pos[0], pos[1], pos[2]);
                
        int layerID = -1;

        if (iColl==0) { // VXD
          streamlog_out( DEBUG ) << " processing collection " << _colNameVTX << std::endl;
          int layerNumber = trkHit->getType() % 100;
          layerID = ILDDetectorIDs::DetID::VXD * ILDDetectorIDs::DetID::Factor + layerNumber - 1 ; // SJA:NOTE: take care of the -1
        }
        else if (iColl==1) { // SIT
          streamlog_out( DEBUG ) << " processing collection " << _colNameSIT << std::endl;
          int layerNumber = trkHit->getType() % 400;
          layerID = ILDDetectorIDs::DetID::SIT * ILDDetectorIDs::DetID::Factor + layerNumber - 1 ; // SJA:NOTE: take care of the -1
        }
        else if (iColl==2) { // TPC
          streamlog_out( DEBUG ) << " processing collection " << _colNameTPC << std::endl;
          const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
          const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

          int padIndex = padLayout.getNearestPad(point.perp(),point.phi());
          int iRowHit = padLayout.getRowNumber(padIndex);

          layerID = ILDDetectorIDs::DetID::TPC * ILDDetectorIDs::DetID::Factor + iRowHit ; 
        }
        else if (iColl==3) { // FTD
          streamlog_out( DEBUG ) << " processing collection " << _colNameFTD << std::endl;
          int layerNumber = trkHit->getType() % 200;
          layerID = ILDDetectorIDs::DetID::FTD * ILDDetectorIDs::DetID::Factor + layerNumber - 1 ; // SJA:NOTE: take care of the -1
        }
        else if (iColl==4) { // SET
          streamlog_out( DEBUG ) << " processing collection " << _colNameSET << std::endl;
          int layerNumber = trkHit->getType() % 400;
          layerID = ILDDetectorIDs::DetID::SET * ILDDetectorIDs::DetID::Factor + layerNumber - 1 ; // SJA:NOTE: take care of the -1
        }
        else if (iColl==5) { // ETD
          streamlog_out( DEBUG ) << " processing collection " << _colNameETD << std::endl;
          int layerNumber = trkHit->getType() % 200;
          layerID = ILDDetectorIDs::DetID::ETD * ILDDetectorIDs::DetID::Factor + layerNumber - 1 ; // SJA:NOTE: take care of the -1
        }


  
        trkHit->ext<ILDDetectorIDs::HitInfo>() = new ILDDetectorIDs::HitInfoStruct ;
        trkHit->ext<ILDDetectorIDs::HitInfo>()->layerID = ( layerID ) ; 

        streamlog_out( DEBUG ) << " hit has layer id " << layerID << std::endl ;

      }          
    }
  }
  _nEvt ++ ;
}



void TrackerHitLayerNumberProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrackerHitLayerNumberProcessor::end(){ 
  
  streamlog_out(MESSAGE) << " end()  " << name() 
                         << " processed " << _nEvt << " events in " << _nRun << " runs "
                         << std::endl ;

}




