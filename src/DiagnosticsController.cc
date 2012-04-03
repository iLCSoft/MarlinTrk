
#include "MarlinTrk/MarlinTrkDiagnostics.h"

#ifdef MARLINTRK_DIAGNOSTICS_ON

#include "MarlinTrk/DiagnosticsController.h"

#include "MarlinTrk/MarlinKalTestTrack.h"

#include "streamlog/streamlog.h"

#include "kaldet/ILDVTrackHit.h"
#include "kaltest/TKalTrackSite.h"
#include "kaltest/TKalTrack.h"

#include "TFile.h"
#include "TTree.h"

#include "MarlinTrk/MarlinTrkNtuple.h"
#include "MarlinTrk/HelixTrack.h"

#include "EVENT/MCParticle.h"
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>




namespace MarlinTrk{
  
  
  DiagnosticsController::DiagnosticsController()   {
    
    _initialised = false;
    _recording_on = false; // recording is set to off by default so that processors do not have to perform any action if they are not interested in diagnostics, i.e. no need to call init
    _current_track = 0;
    _currentMCP = 0;
    _mcpInfoStored=false;
    _skip_track = false;    
    
  } // constructor
  
  
  
  void DiagnosticsController::init(std::string root_file_name, std::string root_tree_name, bool recording_on){
    
    streamlog_out(DEBUG4) << " DiagnosticsController::init called " << "\n"
    << "\t root file name = " << root_file_name << "\n"
    << "\t root tree name = " << root_tree_name << "\n"
    << "\t recording on = " << recording_on << "\n"
    << std::endl;

    _recording_on = recording_on;
    
    if ( _recording_on == false ) {      
      _initialised = true;
      return;
    }
    
    _root_file_name = root_file_name+".root";
    _root_tree_name = root_tree_name;
    
    _root_file = new TFile(_root_file_name.c_str(),"RECREATE");
    
    _tree = new TTree( _root_tree_name.c_str(), _root_tree_name.c_str());
          
    _tree->Branch("nsites", &_track_record.nsites ,"nsites/I" );
    
    _tree->Branch("nsites_vxd", &_track_record.nsites_vxd ,"nsites_vxd/I" );
    _tree->Branch("nsites_sit", &_track_record.nsites_sit ,"nsites_sit§/I" );
    _tree->Branch("nsites_ftd", &_track_record.nsites_ftd ,"nsites_ftd/I" );
    _tree->Branch("nsites_tpc", &_track_record.nsites_tpc ,"nsites_tpc/I" );
    _tree->Branch("nsites_set", &_track_record.nsites_set ,"nsites_set/I" );
    
    _tree->Branch("px_mcp", &_track_record.px_mcp ,"px_mcp/F" );
    _tree->Branch("py_mcp", &_track_record.py_mcp ,"py_mcp/F" );
    _tree->Branch("pz_mcp", &_track_record.pz_mcp ,"pz_mcp/F" );
    _tree->Branch("p_mcp", &_track_record.p_mcp ,"p_mcp/F" );
    _tree->Branch("theta_mcp", &_track_record.theta_mcp ,"theta_mcp/F" );
    _tree->Branch("phi_mcp", &_track_record.phi_mcp ,"phi_mcp/F" );
    _tree->Branch("pdg_mcp", &_track_record.pdg_mcp ,"pdg_mcp/I" );
    
    _tree->Branch("d0_mcp", &_track_record.d0_mcp ,"d0_mcp/F" );
    _tree->Branch("phi0_mcp", &_track_record.phi0_mcp ,"phi0_mcp/F" );
    _tree->Branch("omega_mcp", &_track_record.omega_mcp ,"omega_mcp/F" );
    _tree->Branch("z0_mcp", &_track_record.z0_mcp ,"z0_mcp/F" );
    _tree->Branch("tanL_mcp", &_track_record.tanL_mcp ,"tanL_mcp/F" );
    
    _tree->Branch("d0_seed", &_track_record.d0_seed ,"d0_seed/F" );
    _tree->Branch("phi0_seed", &_track_record.phi0_seed ,"phi0_seed/F" );
    _tree->Branch("omega_seed", &_track_record.omega_seed ,"omega_seed/F" );
    _tree->Branch("z0_seed", &_track_record.z0_seed ,"z0_seed/F" );
    _tree->Branch("tanL_seed", &_track_record.tanL_seed ,"tanL_seed/F" );
    
    _tree->Branch("d0_ip", &_track_record.d0_ip ,"d0_ip/F" );
    _tree->Branch("phi0_ip", &_track_record.phi0_ip ,"phi0_ip/F" );
    _tree->Branch("omega_ip", &_track_record.omega_ip ,"omega_ip/F" );
    _tree->Branch("z0_ip", &_track_record.z0_ip ,"z0_ip/F" );
    _tree->Branch("tanL_ip", &_track_record.tanL_ip ,"tanL_ip/F" );
    
    _tree->Branch("cov_ip_d0d0", &_track_record.cov_ip_d0d0 ,"cov_ip_d0d0/F" );
    _tree->Branch("cov_ip_phi0d0", &_track_record.cov_ip_phi0d0 ,"cov_ip_phi0d0/F" );
    _tree->Branch("cov_ip_phi0phi0", &_track_record.cov_ip_phi0phi0 ,"cov_ip_phi0phi0/F" );
    _tree->Branch("cov_ip_omegad0", &_track_record.cov_ip_omegad0 ,"cov_ip_omegad0/F" );
    _tree->Branch("cov_ip_omegaphi0", &_track_record.cov_ip_omegaphi0 ,"cov_ip_omegaphi0/F" );
    _tree->Branch("cov_ip_omegaomega", &_track_record.cov_ip_omegaomega ,"cov_ip_omegaomega/F" );
    _tree->Branch("cov_ip_z0phi0", &_track_record.cov_ip_z0phi0 ,"cov_ip_z0phi0/F" );
    _tree->Branch("cov_ip_z0omega", &_track_record.cov_ip_z0omega ,"cov_ip_z0omega/F" );
    _tree->Branch("cov_ip_z0z0", &_track_record.cov_ip_z0z0 ,"cov_ip_z0z0/F" );
    _tree->Branch("cov_ip_z0tanL", &_track_record.cov_ip_z0tanL ,"cov_ip_z0tanL/F" );
    _tree->Branch("cov_ip_tanLd0", &_track_record.cov_ip_tanLd0 ,"cov_ip_tanLd0/F" );
    _tree->Branch("cov_ip_tanLphi0", &_track_record.cov_ip_tanLphi0 ,"cov_ip_tanLphi0/F" );
    _tree->Branch("cov_ip_tanLomega", &_track_record.cov_ip_tanLomega ,"cov_ip_tanLomega/F" );
    _tree->Branch("cov_ip_tanLtanL", &_track_record.cov_ip_tanLtanL ,"cov_ip_tanLtanL/F" );

    
    _tree->Branch("ndf", &_track_record.ndf ,"ndf/I" );
    _tree->Branch("chi2", &_track_record.chi2 ,"chi2/F" );
    _tree->Branch("prob", &_track_record.prob ,"prob/F" );
    
    _tree->Branch("CellID0", _track_record.CellID0 ,"CellID0[nsites]/I" );
    _tree->Branch("rejected", _track_record.rejected ,"rejected[nsites]/I" );
    
    
    _tree->Branch("site_x", _track_record.site_x ,"site_x[nsites]/F" );
    _tree->Branch("site_y", _track_record.site_y ,"site_y[nsites]/F" );
    _tree->Branch("site_z", _track_record.site_z ,"site_z[nsites]/F" );
    
    _tree->Branch("ref_point_x", _track_record.ref_point_x ,"ref_point_x[nsites]/F" );
    _tree->Branch("ref_point_y", _track_record.ref_point_y ,"ref_point_y[nsites]/F" );
    _tree->Branch("ref_point_z", _track_record.ref_point_z ,"ref_point_z[nsites]/F" );
    
    _tree->Branch("d0_mc", _track_record.d0_mc ,"d0_mc[nsites]/F" );
    _tree->Branch("phi0_mc", _track_record.phi0_mc ,"phi0_mc[nsites]/F" );
    _tree->Branch("omega_mc", _track_record.omega_mc ,"omega_mc[nsites]/F" );
    _tree->Branch("z0_mc", _track_record.z0_mc ,"z0_mc[nsites]/F" );
    _tree->Branch("tanL_mc", _track_record.tanL_mc ,"tanL_mc[nsites]/F" );
    
    _tree->Branch("d0_predicted", _track_record.d0_predicted ,"d0_predicted[nsites]/F" );
    _tree->Branch("phi0_predicted", _track_record.phi0_predicted ,"phi0_predicted[nsites]/F" );
    _tree->Branch("omega_predicted", _track_record.omega_predicted ,"omega_predicted[nsites]/F" );
    _tree->Branch("z0_predicted", _track_record.z0_predicted ,"z0_predicted[nsites]/F" );
    _tree->Branch("tanL_predicted", _track_record.tanL_predicted ,"tanL_predicted[nsites]/F" );
    
    _tree->Branch("d0_filtered", _track_record.d0_filtered ,"d0_filtered[nsites]/F" );
    _tree->Branch("phi0_filtered", _track_record.phi0_filtered ,"phi0_filtered[nsites]/F" );
    _tree->Branch("omega_filtered", _track_record.omega_filtered ,"omega_filtered[nsites]/F" );
    _tree->Branch("z0_filtered", _track_record.z0_filtered ,"z0_filtered[nsites]/F" );
    _tree->Branch("tanL_filtered", _track_record.tanL_filtered ,"tanL_filtered[nsites]/F" );
    
    _tree->Branch("d0_smoothed", _track_record.d0_smoothed ,"d0_smoothed[nsites]/F" );
    _tree->Branch("phi0_smoothed", _track_record.phi0_smoothed ,"phi0_smoothed[nsites]/F" );
    _tree->Branch("omega_smoothed", _track_record.omega_smoothed ,"omega_smoothed[nsites]/F" );
    _tree->Branch("z0_smoothed", _track_record.z0_smoothed ,"z0_smoothed[nsites]/F" );
    _tree->Branch("tanL_smoothed", _track_record.tanL_smoothed ,"tanL_smoothed[nsites]/F" );
    
    
    _tree->Branch("chi2_inc_filtered", _track_record.chi2_inc_filtered ,"chi2_inc_filtered[nsites]/F" );
    _tree->Branch("chi2_inc_smoothed", _track_record.chi2_inc_smoothed ,"chi2_inc_smoothed[nsites]/F" );
    
    _tree->Branch("dim", _track_record.dim ,"dim[nsites]/I" );
    
    _tree->Branch("cov_smoothed_d0d0", _track_record.cov_smoothed_d0d0 ,"cov_smoothed_d0d0[nsites]/F" );
    _tree->Branch("cov_smoothed_phi0d0", _track_record.cov_smoothed_phi0d0 ,"cov_smoothed_phi0d0[nsites]/F" );
    _tree->Branch("cov_smoothed_phi0phi0", _track_record.cov_smoothed_phi0phi0 ,"cov_smoothed_phi0phi0[nsites]/F" );
    _tree->Branch("cov_smoothed_omegad0", _track_record.cov_smoothed_omegad0 ,"cov_smoothed_omegad0[nsites]/F" );
    _tree->Branch("cov_smoothed_omegaphi0", _track_record.cov_smoothed_omegaphi0 ,"cov_smoothed_omegaphi0[nsites]/F" );
    _tree->Branch("cov_smoothed_omegaomega", _track_record.cov_smoothed_omegaomega ,"cov_smoothed_omegaomega[nsites]/F" );
    _tree->Branch("cov_smoothed_z0phi0", _track_record.cov_smoothed_z0phi0 ,"cov_smoothed_z0phi0[nsites]/F" );
    _tree->Branch("cov_smoothed_z0omega", _track_record.cov_smoothed_z0omega ,"cov_smoothed_z0omega[nsites]/F" );
    _tree->Branch("cov_smoothed_z0z0", _track_record.cov_smoothed_z0z0 ,"cov_smoothed_z0z0[nsites]/F" );
    _tree->Branch("cov_smoothed_z0tanL", _track_record.cov_smoothed_z0tanL ,"cov_smoothed_z0tanL[nsites]/F" );
    _tree->Branch("cov_smoothed_tanLd0", _track_record.cov_smoothed_tanLd0 ,"cov_smoothed_tanLd0[nsites]/F" );
    _tree->Branch("cov_smoothed_tanLphi0", _track_record.cov_smoothed_tanLphi0 ,"cov_smoothed_tanLphi0[nsites]/F" );
    _tree->Branch("cov_smoothed_tanLomega", _track_record.cov_smoothed_tanLomega ,"cov_smoothed_tanLomega[nsites]/F" );
    _tree->Branch("cov_smoothed_tanLtanL", _track_record.cov_smoothed_tanLtanL ,"cov_smoothed_tanLtanL[nsites]/F" );
    
    _tree->Branch("cov_predicted_d0d0", _track_record.cov_predicted_d0d0 ,"cov_predicted_d0d0[nsites]/F" );
    _tree->Branch("cov_predicted_phi0d0", _track_record.cov_predicted_phi0d0 ,"cov_predicted_phi0d0[nsites]/F" );
    _tree->Branch("cov_predicted_phi0phi0", _track_record.cov_predicted_phi0phi0 ,"cov_predicted_phi0phi0[nsites]/F" );
    _tree->Branch("cov_predicted_omegad0", _track_record.cov_predicted_omegad0 ,"cov_predicted_omegad0[nsites]/F" );
    _tree->Branch("cov_predicted_omegaphi0", _track_record.cov_predicted_omegaphi0 ,"cov_predicted_omegaphi0[nsites]/F" );
    _tree->Branch("cov_predicted_omegaomega", _track_record.cov_predicted_omegaomega ,"cov_predicted_omegaomega[nsites]/F" );
    _tree->Branch("cov_predicted_z0phi0", _track_record.cov_predicted_z0phi0 ,"cov_predicted_z0phi0[nsites]/F" );
    _tree->Branch("cov_predicted_z0omega", _track_record.cov_predicted_z0omega ,"cov_predicted_z0omega[nsites]/F" );
    _tree->Branch("cov_predicted_z0z0", _track_record.cov_predicted_z0z0 ,"cov_predicted_z0z0[nsites]/F" );
    _tree->Branch("cov_predicted_z0tanL", _track_record.cov_predicted_z0tanL ,"cov_predicted_z0tanL[nsites]/F" );
    _tree->Branch("cov_predicted_tanLd0", _track_record.cov_predicted_tanLd0 ,"cov_predicted_tanLd0[nsites]/F" );
    _tree->Branch("cov_predicted_tanLphi0", _track_record.cov_predicted_tanLphi0 ,"cov_predicted_tanLphi0[nsites]/F" );
    _tree->Branch("cov_predicted_tanLomega", _track_record.cov_predicted_tanLomega ,"cov_predicted_tanLomega[nsites]/F" );
    _tree->Branch("cov_predicted_tanLtanL", _track_record.cov_predicted_tanLtanL ,"cov_predicted_tanLtanL[nsites]/F" );
    
    _tree->Branch("cov_filtered_d0d0", _track_record.cov_filtered_d0d0 ,"cov_filtered_d0d0[nsites]/F" );
    _tree->Branch("cov_filtered_phi0d0", _track_record.cov_filtered_phi0d0 ,"cov_filtered_phi0d0[nsites]/F" );
    _tree->Branch("cov_filtered_phi0phi0", _track_record.cov_filtered_phi0phi0 ,"cov_filtered_phi0phi0[nsites]/F" );
    _tree->Branch("cov_filtered_omegad0", _track_record.cov_filtered_omegad0 ,"cov_filtered_omegad0[nsites]/F" );
    _tree->Branch("cov_filtered_omegaphi0", _track_record.cov_filtered_omegaphi0 ,"cov_filtered_omegaphi0[nsites]/F" );
    _tree->Branch("cov_filtered_omegaomega", _track_record.cov_filtered_omegaomega ,"cov_filtered_omegaomega[nsites]/F" );
    _tree->Branch("cov_filtered_z0phi0", _track_record.cov_filtered_z0phi0 ,"cov_filtered_z0phi0[nsites]/F" );
    _tree->Branch("cov_filtered_z0omega", _track_record.cov_filtered_z0omega ,"cov_filtered_z0omega[nsites]/F" );
    _tree->Branch("cov_filtered_z0z0", _track_record.cov_filtered_z0z0 ,"cov_filtered_z0z0[nsites]/F" );
    _tree->Branch("cov_filtered_z0tanL", _track_record.cov_filtered_z0tanL ,"cov_filtered_z0tanL[nsites]/F" );
    _tree->Branch("cov_filtered_tanLd0", _track_record.cov_filtered_tanLd0 ,"cov_filtered_tanLd0[nsites]/F" );
    _tree->Branch("cov_filtered_tanLphi0", _track_record.cov_filtered_tanLphi0 ,"cov_filtered_tanLphi0[nsites]/F" );
    _tree->Branch("cov_filtered_tanLomega", _track_record.cov_filtered_tanLomega ,"cov_filtered_tanLomega[nsites]/F" );
    _tree->Branch("cov_filtered_tanLtanL", _track_record.cov_filtered_tanLtanL ,"cov_filtered_tanLtanL[nsites]/F" );
    
    
    _initialised = true;
    
  }
  
    
  void DiagnosticsController::clear_track_record(){
    
    streamlog_out(DEBUG) << " DiagnosticsController::clear_track_record called " << std::endl;
    
    _track_record.nsites= 0 ;
    
    _track_record.nsites_vxd= 0 ;
    _track_record.nsites_sit= 0 ;
    _track_record.nsites_ftd = 0 ;
    _track_record.nsites_tpc = 0 ;
    _track_record.nsites_set = 0 ;
    
    _track_record.px_mcp = 0 ;
    _track_record.py_mcp = 0 ;
    _track_record.pz_mcp = 0 ;
    _track_record.p_mcp = 0 ;
    _track_record.theta_mcp = 0 ;
    _track_record.phi_mcp = 0 ;
    _track_record.pdg_mcp = 0 ;
    
    _track_record.d0_mcp = 0 ;
    _track_record.phi0_mcp = 0 ;
    _track_record.omega_mcp = 0 ;
    _track_record.z0_mcp = 0 ;
    _track_record.tanL_mcp = 0 ;
    
    _track_record.ndf = 0 ;
    _track_record.chi2 = 0 ;
    _track_record.prob = 0 ;
    
    _track_record.d0_seed = 0 ;
    _track_record.phi0_seed = 0 ;
    _track_record.omega_seed = 0 ;
    _track_record.z0_seed = 0 ;
    _track_record.tanL_seed = 0 ; 
    
    _track_record.d0_ip = 0 ;
    _track_record.phi0_ip = 0 ;
    _track_record.omega_ip = 0 ;
    _track_record.z0_ip = 0 ;
    _track_record.tanL_ip = 0 ; 
    
    _track_record.cov_ip_d0d0 = 0 ;      
    _track_record.cov_ip_phi0d0 = 0 ;      
    _track_record.cov_ip_phi0phi0 = 0 ;      
    _track_record.cov_ip_omegad0 = 0 ;      
    _track_record.cov_ip_omegaphi0 = 0 ;      
    _track_record.cov_ip_omegaomega = 0 ;      
    _track_record.cov_ip_z0d0 = 0 ;      
    _track_record.cov_ip_z0phi0 = 0 ;      
    _track_record.cov_ip_z0omega = 0 ;      
    _track_record.cov_ip_z0z0 = 0 ;      
    _track_record.cov_ip_z0tanL = 0 ;      
    _track_record.cov_ip_tanLd0 = 0 ;      
    _track_record.cov_ip_tanLphi0 = 0 ;      
    _track_record.cov_ip_tanLomega = 0 ;      
    _track_record.cov_ip_tanLtanL = 0 ;  
    
    for ( int i = 0 ; i<MAX_SITES; ++i) {
      
      _track_record.CellID0[i] = 0 ;
      _track_record.rejected[i] = 0 ;
      
      _track_record.site_x[i] = 0 ;
      _track_record.site_y[i] = 0 ;
      _track_record.site_z[i] = 0 ;
      
      _track_record.ref_point_x[i] = 0 ;
      _track_record.ref_point_y[i] = 0 ;
      _track_record.ref_point_z[i] = 0 ;
      
      _track_record.d0_mc[i] = 0 ;
      _track_record.phi0_mc[i] = 0 ;
      _track_record.omega_mc[i] = 0 ;
      _track_record.z0_mc[i] = 0 ;
      _track_record.tanL_mc[i] = 0 ;
      
      _track_record.d0_predicted[i] = 0 ;
      _track_record.phi0_predicted[i] = 0 ;
      _track_record.omega_predicted[i] = 0 ;
      _track_record.z0_predicted[i] = 0 ;
      _track_record.tanL_predicted[i] = 0 ; 
      
      _track_record.d0_filtered[i] = 0 ;
      _track_record.phi0_filtered[i] = 0 ;
      _track_record.omega_filtered[i] = 0 ;
      _track_record.z0_filtered[i] = 0 ;
      _track_record.tanL_filtered[i] = 0 ; 
      
      _track_record.d0_smoothed[i] = 0 ;
      _track_record.phi0_smoothed[i] = 0 ;
      _track_record.omega_smoothed[i] = 0 ;
      _track_record.z0_smoothed[i] = 0 ;
      _track_record.tanL_smoothed[i] = 0 ; 
      
      
      _track_record.chi2_inc_filtered[i] = 0 ;
      _track_record.chi2_inc_smoothed[i] = 0 ;
      _track_record.dim[i] = 0 ;
      
      _track_record.cov_predicted_d0d0[i] = 0 ;      
      _track_record.cov_predicted_phi0d0[i] = 0 ;      
      _track_record.cov_predicted_phi0phi0[i] = 0 ;      
      _track_record.cov_predicted_omegad0[i] = 0 ;      
      _track_record.cov_predicted_omegaphi0[i] = 0 ;      
      _track_record.cov_predicted_omegaomega[i] = 0 ;      
      _track_record.cov_predicted_z0d0[i] = 0 ;      
      _track_record.cov_predicted_z0phi0[i] = 0 ;      
      _track_record.cov_predicted_z0omega[i] = 0 ;      
      _track_record.cov_predicted_z0z0[i] = 0 ;      
      _track_record.cov_predicted_z0tanL[i] = 0 ;      
      _track_record.cov_predicted_tanLd0[i] = 0 ;      
      _track_record.cov_predicted_tanLphi0[i] = 0 ;      
      _track_record.cov_predicted_tanLomega[i] = 0 ;      
      _track_record.cov_predicted_tanLtanL[i] = 0 ;      
      
      _track_record.cov_filtered_d0d0[i] = 0 ;      
      _track_record.cov_filtered_phi0d0[i] = 0 ;      
      _track_record.cov_filtered_phi0phi0[i] = 0 ;      
      _track_record.cov_filtered_omegad0[i] = 0 ;      
      _track_record.cov_filtered_omegaphi0[i] = 0 ;      
      _track_record.cov_filtered_omegaomega[i] = 0 ;      
      _track_record.cov_filtered_z0d0[i] = 0 ;      
      _track_record.cov_filtered_z0phi0[i] = 0 ;      
      _track_record.cov_filtered_z0omega[i] = 0 ;      
      _track_record.cov_filtered_z0z0[i] = 0 ;      
      _track_record.cov_filtered_z0tanL[i] = 0 ;      
      _track_record.cov_filtered_tanLd0[i] = 0 ;      
      _track_record.cov_filtered_tanLphi0[i] = 0 ;      
      _track_record.cov_filtered_tanLomega[i] = 0 ;      
      _track_record.cov_filtered_tanLtanL[i] = 0 ;      
      
      _track_record.cov_smoothed_d0d0[i] = 0 ;      
      _track_record.cov_smoothed_phi0d0[i] = 0 ;      
      _track_record.cov_smoothed_phi0phi0[i] = 0 ;      
      _track_record.cov_smoothed_omegad0[i] = 0 ;      
      _track_record.cov_smoothed_omegaphi0[i] = 0 ;      
      _track_record.cov_smoothed_omegaomega[i] = 0 ;      
      _track_record.cov_smoothed_z0d0[i] = 0 ;      
      _track_record.cov_smoothed_z0phi0[i] = 0 ;      
      _track_record.cov_smoothed_z0omega[i] = 0 ;      
      _track_record.cov_smoothed_z0z0[i] = 0 ;      
      _track_record.cov_smoothed_z0tanL[i] = 0 ;      
      _track_record.cov_smoothed_tanLd0[i] = 0 ;      
      _track_record.cov_smoothed_tanLphi0[i] = 0 ;      
      _track_record.cov_smoothed_tanLomega[i] = 0 ;      
      _track_record.cov_smoothed_tanLtanL[i] = 0 ;      
      
      
      
    }
    
  }
  
  void DiagnosticsController::new_track(MarlinKalTestTrack* trk){
    
    if ( _recording_on == false ) {
      return;
    }
    
    if ( _initialised == false ){
      streamlog_out(ERROR) << "DiagnosticsController::new_track: Diagnostics not initialised call init(std::string root_file_name, std::string root_tree_name, bool recording_off) first : exit(1) called from file " 
      << __FILE__
      << " line "
      << __LINE__
      << std::endl;
      
      exit(1);
    }
    
    streamlog_out(DEBUG3) << " DiagnosticsController::new_track called " << std::endl;
    
    this->clear_track_record();
    
    _current_track = trk;
        
    _skip_track = false;    
    _currentMCP = NULL;
    _mcpInfoStored=false;
    
  }
  
  void DiagnosticsController::end_track(){
    
    if ( _recording_on == false ) {
      return;
    }
    
    if ( _initialised == false ){
      streamlog_out(ERROR) << "DiagnosticsController::end_track: Diagnostics not initialised call init(std::string root_file_name, std::string root_tree_name, bool recording_off) first : exit(1) called from file " 
      << __FILE__
      << " line "
      << __LINE__
      << std::endl;
      
      exit(1);
    }

    streamlog_out(DEBUG3) << " DiagnosticsController::end_track called " << std::endl;
    
    if ( _skip_track ) { // just clear the track buffers and return.
      this->clear_track_record();
      return;
    } else {
    
      _track_record.chi2 = _current_track->_kaltrack->GetChi2();
      _track_record.ndf = _current_track->_kaltrack->GetNDF();
      _track_record.prob = TMath::Prob(_track_record.chi2, _track_record.ndf);
      
      TIter it(_current_track->_kaltrack,kIterForward);
      
      Int_t nsites =  _current_track->_kaltrack->GetEntries();
      
      if(_current_track->_smoothed){
        
        for (Int_t isite=1; isite<nsites; isite++) {
          
          UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
          encoder.setValue( _track_record.CellID0[isite] );
          
          
          TVKalSite* site = static_cast<TVKalSite *>( _current_track->_kaltrack->At(isite));
          
          if ( _track_record.rejected[isite] == 0 && encoder[lcio::ILDCellID0::subdet] != 0 ) {
            
            
            _track_record.chi2_inc_smoothed[isite] = site->GetDeltaChi2();
            
            TKalTrackState& trkState_smoothed = (TKalTrackState&) site->GetState(TVKalSite::kSmoothed); // this segfaults if no hits are present
            
            THelicalTrack helix_smoothed = trkState_smoothed.GetHelix() ;
            TMatrixD c0_smoothed(trkState_smoothed.GetCovMat());  
            
            Int_t sdim = trkState_smoothed.GetDimension();  // dimensions of the track state, it will be 5 or 6
            
            // move track state to the sim hit for comparison 
            const TVector3 tpoint( _track_record.ref_point_x[isite], _track_record.ref_point_y[isite], _track_record.ref_point_z[isite] ) ;
            
            // now move to the point
            TKalMatrix  DF(sdim,sdim);  
            DF.UnitMatrix();                           
            
            double dPhi=0;
            
            helix_smoothed.MoveTo(  tpoint , dPhi , &DF , &c0_smoothed) ;  // move helix to desired point, and get propagator matrix
            
            IMPL::TrackStateImpl ts;
            int ndf;
            double chi2;
            
            _current_track->ToLCIOTrackState( helix_smoothed, c0_smoothed, ts, chi2, ndf );
            
            
            _track_record.d0_smoothed[isite] = ts.getD0() ;
            _track_record.phi0_smoothed[isite] = ts.getPhi() ;
            _track_record.omega_smoothed[isite] = ts.getOmega() ;
            _track_record.z0_smoothed[isite] = ts.getZ0() ;
            _track_record.tanL_smoothed[isite] = ts.getTanLambda() ;
            
            
            _track_record.cov_smoothed_d0d0[isite] = ts.getCovMatrix()[0];      
            _track_record.cov_smoothed_phi0d0[isite] = ts.getCovMatrix()[1];      
            _track_record.cov_smoothed_phi0phi0[isite] = ts.getCovMatrix()[2];      
            _track_record.cov_smoothed_omegad0[isite] = ts.getCovMatrix()[3];      
            _track_record.cov_smoothed_omegaphi0[isite] = ts.getCovMatrix()[4];      
            _track_record.cov_smoothed_omegaomega[isite] = ts.getCovMatrix()[5];      
            _track_record.cov_smoothed_z0d0[isite] = ts.getCovMatrix()[6];      
            _track_record.cov_smoothed_z0phi0[isite] = ts.getCovMatrix()[7];      
            _track_record.cov_smoothed_z0omega[isite] = ts.getCovMatrix()[8];      
            _track_record.cov_smoothed_z0z0[isite] = ts.getCovMatrix()[9];      
            _track_record.cov_smoothed_z0tanL[isite] = ts.getCovMatrix()[10];      
            _track_record.cov_smoothed_tanLd0[isite] = ts.getCovMatrix()[11];      
            _track_record.cov_smoothed_tanLphi0[isite] = ts.getCovMatrix()[12];      
            _track_record.cov_smoothed_tanLomega[isite] = ts.getCovMatrix()[13];      
            _track_record.cov_smoothed_tanLtanL[isite] = ts.getCovMatrix()[14];  
            
          }
          
        }
      }
      
      IMPL::TrackStateImpl ts_at_ip;
      int ndf;
      double chi2;
      
      gear::Vector3D point(0.0,0.0,0.0);
      _current_track->propagate( point, ts_at_ip, chi2, ndf );
      
      _track_record.d0_ip = ts_at_ip.getD0() ;
      _track_record.phi0_ip = ts_at_ip.getPhi() ;
      _track_record.omega_ip = ts_at_ip.getOmega() ;
      _track_record.z0_ip = ts_at_ip.getZ0() ;
      _track_record.tanL_ip = ts_at_ip.getTanLambda() ;
      
      
      _track_record.cov_ip_d0d0 = ts_at_ip.getCovMatrix()[0];      
      _track_record.cov_ip_phi0d0 = ts_at_ip.getCovMatrix()[1];      
      _track_record.cov_ip_phi0phi0 = ts_at_ip.getCovMatrix()[2];      
      _track_record.cov_ip_omegad0 = ts_at_ip.getCovMatrix()[3];      
      _track_record.cov_ip_omegaphi0 = ts_at_ip.getCovMatrix()[4];      
      _track_record.cov_ip_omegaomega = ts_at_ip.getCovMatrix()[5];      
      _track_record.cov_ip_z0d0 = ts_at_ip.getCovMatrix()[6];      
      _track_record.cov_ip_z0phi0 = ts_at_ip.getCovMatrix()[7];      
      _track_record.cov_ip_z0omega = ts_at_ip.getCovMatrix()[8];      
      _track_record.cov_ip_z0z0 = ts_at_ip.getCovMatrix()[9];      
      _track_record.cov_ip_z0tanL = ts_at_ip.getCovMatrix()[10];      
      _track_record.cov_ip_tanLd0 = ts_at_ip.getCovMatrix()[11];      
      _track_record.cov_ip_tanLphi0 = ts_at_ip.getCovMatrix()[12];      
      _track_record.cov_ip_tanLomega = ts_at_ip.getCovMatrix()[13];      
      _track_record.cov_ip_tanLtanL = ts_at_ip.getCovMatrix()[14];  
      
      
      _tree->Fill();
      
    }
    
  }
  
  void DiagnosticsController::end(){
    
    if ( _recording_on == false ) {
      return;
    }
    
    if ( _initialised == false ){
      streamlog_out(ERROR) << "DiagnosticsController::end: Diagnostics not initialised call init(std::string root_file_name, std::string root_tree_name, bool recording_off) first : exit(1) called from file " 
      << __FILE__
      << " line "
      << __LINE__
      << std::endl;
      
      exit(1);
    }
          
    streamlog_out(DEBUG4) << " DiagnosticsController::end() called " << std::endl;
    
    //    _tree->Print();
    _root_file->Write();
    _root_file->Close();
    
  }
  
  
  void DiagnosticsController::set_intial_track_parameters(double d0, double phi0, double omega, double z0, double tanL){
    
    if ( _recording_on == false ) {
      return;
    }
    
    if ( _initialised == false ){
      streamlog_out(ERROR) << "DiagnosticsController::set_intial_track_parameters: Diagnostics not initialised call init(std::string root_file_name, std::string root_tree_name, bool recording_off) first : exit(1) called from file " 
      << __FILE__
      << " line "
      << __LINE__
      << std::endl;
      
      exit(1);
    }
    
    streamlog_out(DEBUG3) << " DiagnosticsController::set_intial_track_parameters called " << std::endl;
    
    _track_record.d0_seed= d0;
    _track_record.phi0_seed= phi0;
    _track_record.omega_seed= omega;
    _track_record.z0_seed= z0;
    _track_record.tanL_seed= tanL;
    
    streamlog_out(DEBUG3) << " $#$#$# Initial Track Parameters: " 
    << "d0 = " << d0 << " "  
    << "phi0 = " << phi0 << " "  
    << "omega = " << omega << " "  
    << "z0 = " << z0 << " "  
    << "tanL = " << tanL << " "  
    << std::endl;
    
  }
  
  
  void DiagnosticsController::record_rejected_site(ILDVTrackHit* hit, TKalTrackSite* site){

    if ( _recording_on == false ) {
      return;
    }
    
    if ( _initialised == false ){
      streamlog_out(ERROR) << "DiagnosticsController::record_rejected_site: Diagnostics not initialised call init(std::string root_file_name, std::string root_tree_name, bool recording_off) first : exit(1) called from file " 
      << __FILE__
      << " line "
      << __LINE__
      << std::endl;
      
      exit(1);
    }

      
    if(_skip_track) return;
    
    _track_record.rejected[_track_record.nsites] = 1;
    _track_record.CellID0[_track_record.nsites] = hit->getLCIOTrackerHit()->getCellID0();
    
    ++_track_record.nsites;
    streamlog_out(DEBUG2) << "record_rejected_site _track_record.nsites = " << _track_record.nsites << std::endl;

  }
  
  
  void DiagnosticsController::record_site(ILDVTrackHit* hit, TKalTrackSite* site){
        
    if ( _recording_on == false ) {
      return;
    }
    
    if ( _initialised == false ){
      streamlog_out(ERROR) << "DiagnosticsController::record_site: Diagnostics not initialised call init(std::string root_file_name, std::string root_tree_name, bool recording_off) first : exit(1) called from file " 
      << __FILE__
      << " line "
      << __LINE__
      << std::endl;
      
      exit(1);
    }
    
    streamlog_out(DEBUG2) << "DiagnosticsController::record_site called " << std::endl;
    
    if(_skip_track) return;
    
    EVENT::TrackerHit* trkhit = hit->getLCIOTrackerHit();
    EVENT::SimTrackerHit* simhit = trkhit->ext<MarlinTrk::MCTruth4HitExt>()->simhit;
    
    if ( ! simhit ) {
      
      streamlog_out(ERROR) << "SimTrackerHit not attached to TrackerHit using MCTruth4HitExt: exit(1) called from file " 
      << __FILE__
      << " line "
      << __LINE__
      << std::endl;
      
      exit(1);
      
    }
    
    EVENT::MCParticle* mcp = simhit->getMCParticle() ; 
    
    
    if( _track_record.nsites>-1 ){
      
      if ( _mcpInfoStored == false ) {
        _currentMCP = mcp;
        _track_record.px_mcp = mcp->getMomentum()[0];
        _track_record.py_mcp = mcp->getMomentum()[1];
        _track_record.pz_mcp = mcp->getMomentum()[2];
        
        double pt = sqrt(_track_record.px_mcp*_track_record.px_mcp + 
                         _track_record.py_mcp*_track_record.py_mcp ) ;
        
        _track_record.p_mcp = sqrt( pt*pt + _track_record.pz_mcp*_track_record.pz_mcp );
        
        
        _track_record.theta_mcp = atan2( pt, _track_record.pz_mcp );
        _track_record.phi_mcp   = atan2( _track_record.py_mcp, _track_record.px_mcp );
        
        _track_record.pdg_mcp = mcp->getPDG();
        
        //    HelixTrack helixMC( sim_pos, sim_mom, mcp->getCharge(), ml.GetBz() ) ;
        HelixTrack helixMC( mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), site->GetBfield() ) ;
        
        _track_record.d0_mcp= helixMC.getD0();
        _track_record.phi0_mcp = helixMC.getPhi0();
        _track_record.omega_mcp = helixMC.getOmega();
        _track_record.z0_mcp = helixMC.getZ0();
        _track_record.tanL_mcp = helixMC.getTanLambda();
        
        
      }
      else if( _currentMCP != mcp ) {
        _skip_track = true ; // do not store tracks formed from more than one MCParticle
        streamlog_out(WARNING) << "DiagnosticsController::record_site: Track skipped. Not storing tracks formed from more than one MCParticle " << std::endl;
        return ;
      }
      
      double sim_pos[3];
      sim_pos[0] = simhit->getPosition()[0];
      sim_pos[1] = simhit->getPosition()[1];
      sim_pos[2] = simhit->getPosition()[2];
      
      double sim_mom[3];
      sim_mom[0] = simhit->getMomentum()[0];
      sim_mom[1] = simhit->getMomentum()[1];
      sim_mom[2] = simhit->getMomentum()[2];
      
      if ( fabs(sim_mom[0]) < 1.e-09 && fabs(sim_mom[1]) < 1.e-09 && fabs(sim_mom[2]) < 1.e-09 ) {
        // then the momentum is sub eV and therefore the momentum was certainly not set
        streamlog_out(ERROR) << "Momentum not set in SimTrackerHit exit(1) called from file " 
        << __FILE__
        << " line "
        << __LINE__
        << std::endl;
        
        exit(1);
        
      }
      
                  
      _track_record.CellID0[_track_record.nsites] = trkhit->getCellID0() ;
      
      UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
      encoder.setValue( trkhit->getCellID0() );
      
      if (encoder[lcio::ILDCellID0::subdet] == lcio::ILDDetID::VXD) {
        ++_track_record.nsites_vxd;
      }
      else if (encoder[lcio::ILDCellID0::subdet] == lcio::ILDDetID::SIT) {
        ++_track_record.nsites_sit;
      }
      else if (encoder[lcio::ILDCellID0::subdet] == lcio::ILDDetID::FTD) {
        ++_track_record.nsites_ftd;
      }
      else if (encoder[lcio::ILDCellID0::subdet] == lcio::ILDDetID::TPC) {
        ++_track_record.nsites_tpc;
      }
      else if (encoder[lcio::ILDCellID0::subdet] == lcio::ILDDetID::SET) {
        ++_track_record.nsites_set;
      }
      
      _track_record.chi2_inc_filtered[_track_record.nsites] = site->GetDeltaChi2();
      _track_record.dim[_track_record.nsites] = site->GetDimension();
      
      // create helix from MCTruth at sim hit
      //    HelixTrack helixMC( sim_pos, sim_mom, mcp->getCharge(), ml.GetBz() ) ;
      HelixTrack helixMC( sim_pos, sim_mom, mcp->getCharge(), site->GetBfield() ) ;

      // here perhaps we should move the helix to the hit to calculate the deltas or though this could still be done in the analysis code, as both sim and rec hit positions are stored ?
      
      streamlog_out(DEBUG2) << " $#$#$# SimHit Track Parameters: " 
      << "d0 = " << helixMC.getD0() << " "  
      << "phi0 = " << helixMC.getPhi0() << " "  
      << "omega = " << helixMC.getOmega() << " "  
      << "z0 = " << helixMC.getZ0() << " "  
      << "tanL = " << helixMC.getTanLambda() << " "  
      << std::endl;
      
      _track_record.d0_mc[_track_record.nsites] = helixMC.getD0();
      _track_record.phi0_mc[_track_record.nsites] = helixMC.getPhi0();
      _track_record.omega_mc[_track_record.nsites] = helixMC.getOmega();
      _track_record.z0_mc[_track_record.nsites] = helixMC.getZ0();
      _track_record.tanL_mc[_track_record.nsites] = helixMC.getTanLambda();
      

      
      double rec_x = trkhit->getPosition()[0];
      double rec_y = trkhit->getPosition()[1];
      double rec_z = trkhit->getPosition()[2];
      
      _track_record.site_x[_track_record.nsites] = rec_x;
      _track_record.site_y[_track_record.nsites] = rec_y;
      _track_record.site_z[_track_record.nsites] = rec_z;
      
      
      // move track state to the sim hit for comparison 
      const TVector3 tpoint( sim_pos[0], sim_pos[1], sim_pos[2] ) ;
      
      IMPL::TrackStateImpl ts;
      int ndf;
      double chi2;
      double dPhi ;
      
      
      TKalTrackState& trkState_predicted = (TKalTrackState&) site->GetState(TVKalSite::kPredicted); // this segfaults if no hits are present
      
      THelicalTrack helix_predicted = trkState_predicted.GetHelix() ;
      TMatrixD c0_predicted(trkState_predicted.GetCovMat());  
      
      Int_t sdim = trkState_predicted.GetDimension();  // dimensions of the track state, it will be 5 or 6
      
      // now move to the point
      TKalMatrix  DF(sdim,sdim);  
      DF.UnitMatrix();                           
      helix_predicted.MoveTo(  tpoint , dPhi , &DF , &c0_predicted) ;  // move helix to desired point, and get propagator matrix
      
      _current_track->ToLCIOTrackState( helix_predicted, c0_predicted, ts, chi2, ndf );
      
      _track_record.d0_predicted[_track_record.nsites] = ts.getD0() ;
      _track_record.phi0_predicted[_track_record.nsites] = ts.getPhi() ;
      _track_record.omega_predicted[_track_record.nsites] = ts.getOmega() ;
      _track_record.z0_predicted[_track_record.nsites] = ts.getZ0() ;
      _track_record.tanL_predicted[_track_record.nsites] = ts.getTanLambda() ;
      
      
      _track_record.cov_predicted_d0d0[_track_record.nsites] = ts.getCovMatrix()[0];      
      _track_record.cov_predicted_phi0d0[_track_record.nsites] = ts.getCovMatrix()[1];      
      _track_record.cov_predicted_phi0phi0[_track_record.nsites] = ts.getCovMatrix()[2];      
      _track_record.cov_predicted_omegad0[_track_record.nsites] = ts.getCovMatrix()[3];      
      _track_record.cov_predicted_omegaphi0[_track_record.nsites] = ts.getCovMatrix()[4];      
      _track_record.cov_predicted_omegaomega[_track_record.nsites] = ts.getCovMatrix()[5];      
      _track_record.cov_predicted_z0d0[_track_record.nsites] = ts.getCovMatrix()[6];      
      _track_record.cov_predicted_z0phi0[_track_record.nsites] = ts.getCovMatrix()[7];      
      _track_record.cov_predicted_z0omega[_track_record.nsites] = ts.getCovMatrix()[8];      
      _track_record.cov_predicted_z0z0[_track_record.nsites] = ts.getCovMatrix()[9];      
      _track_record.cov_predicted_z0tanL[_track_record.nsites] = ts.getCovMatrix()[10];      
      _track_record.cov_predicted_tanLd0[_track_record.nsites] = ts.getCovMatrix()[11];      
      _track_record.cov_predicted_tanLphi0[_track_record.nsites] = ts.getCovMatrix()[12];      
      _track_record.cov_predicted_tanLomega[_track_record.nsites] = ts.getCovMatrix()[13];      
      _track_record.cov_predicted_tanLtanL[_track_record.nsites] = ts.getCovMatrix()[14];  
      
      
      TKalTrackState& trkState_filtered = (TKalTrackState&) site->GetState(TVKalSite::kFiltered); // this segfaults if no hits are present
      
      THelicalTrack helix_filtered = trkState_filtered.GetHelix() ;
      TMatrixD c0_filtered(trkState_filtered.GetCovMat());  
      
      DF.UnitMatrix();                           
      helix_filtered.MoveTo(  tpoint , dPhi , &DF , &c0_filtered) ;  // move helix to desired point, and get propagator matrix
      
      IMPL::TrackStateImpl ts_f;
      
      _current_track->ToLCIOTrackState( helix_filtered, c0_filtered, ts_f, chi2, ndf );
      
      _track_record.d0_filtered[_track_record.nsites] = ts_f.getD0() ;
      _track_record.phi0_filtered[_track_record.nsites] = ts_f.getPhi() ;
      _track_record.omega_filtered[_track_record.nsites] = ts_f.getOmega() ;
      _track_record.z0_filtered[_track_record.nsites] = ts_f.getZ0() ;
      _track_record.tanL_filtered[_track_record.nsites] = ts_f.getTanLambda() ;
      
      
      _track_record.cov_filtered_d0d0[_track_record.nsites] = ts_f.getCovMatrix()[0];      
      _track_record.cov_filtered_phi0d0[_track_record.nsites] = ts_f.getCovMatrix()[1];      
      _track_record.cov_filtered_phi0phi0[_track_record.nsites] = ts_f.getCovMatrix()[2];      
      _track_record.cov_filtered_omegad0[_track_record.nsites] = ts_f.getCovMatrix()[3];      
      _track_record.cov_filtered_omegaphi0[_track_record.nsites] = ts_f.getCovMatrix()[4];      
      _track_record.cov_filtered_omegaomega[_track_record.nsites] = ts_f.getCovMatrix()[5];      
      _track_record.cov_filtered_z0d0[_track_record.nsites] = ts_f.getCovMatrix()[6];      
      _track_record.cov_filtered_z0phi0[_track_record.nsites] = ts_f.getCovMatrix()[7];      
      _track_record.cov_filtered_z0omega[_track_record.nsites] = ts_f.getCovMatrix()[8];      
      _track_record.cov_filtered_z0z0[_track_record.nsites] = ts_f.getCovMatrix()[9];      
      _track_record.cov_filtered_z0tanL[_track_record.nsites] = ts_f.getCovMatrix()[10];      
      _track_record.cov_filtered_tanLd0[_track_record.nsites] = ts_f.getCovMatrix()[11];      
      _track_record.cov_filtered_tanLphi0[_track_record.nsites] = ts_f.getCovMatrix()[12];      
      _track_record.cov_filtered_tanLomega[_track_record.nsites] = ts_f.getCovMatrix()[13];      
      _track_record.cov_filtered_tanLtanL[_track_record.nsites] = ts_f.getCovMatrix()[14];   
      
      
      _track_record.ref_point_x[_track_record.nsites] = ts_f.getReferencePoint()[0];
      _track_record.ref_point_y[_track_record.nsites] = ts_f.getReferencePoint()[1];
      _track_record.ref_point_z[_track_record.nsites] = ts_f.getReferencePoint()[2];
      
    }
    
    ++_track_record.nsites;
    streamlog_out(DEBUG2) << "_track_record.nsites = " << _track_record.nsites << std::endl;
  }
  
}

#endif
