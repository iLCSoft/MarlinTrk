//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 14 16:01:15 2011 by ROOT version 5.28/00f
// from TTree mytree/a simple Tree with simple variables
// found on file: mytree.root
//////////////////////////////////////////////////////////

#ifndef TrkAnalysisTree_h
#define TrkAnalysisTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define MAX_MCPARTICLES 500
#define MAX_NTRACKS 500
#define MAX_MC_TRACK_RELS 1500 
#define MAX_TRACK_MC_RELS 1500 
#define MAX_ALLOWED_TRACK_MC_LINKS 3

class TrkAnalysisTree {
  public :
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leaf types
  
  Int_t nmcp;
  Int_t runNumber;
  Int_t eventNumber;
  
  Int_t mcp_pdg[MAX_MCPARTICLES];
  Int_t mcp_generator_status[MAX_MCPARTICLES];
  Int_t mcp_simulator_status[MAX_MCPARTICLES];
  Int_t mcp_ndaughters[MAX_MCPARTICLES];
  Int_t mcp_nparents[MAX_MCPARTICLES];
  Float_t mcp_mass[MAX_MCPARTICLES];
  Float_t mcp_charge[MAX_MCPARTICLES];
  Float_t mcp_energy[MAX_MCPARTICLES];
  Float_t mcp_px[MAX_MCPARTICLES];
  Float_t mcp_py[MAX_MCPARTICLES];
  Float_t mcp_pz[MAX_MCPARTICLES];
  Float_t mcp_vx[MAX_MCPARTICLES];
  Float_t mcp_vy[MAX_MCPARTICLES];
  Float_t mcp_vz[MAX_MCPARTICLES];
  Float_t mcp_d0[MAX_MCPARTICLES];
  Float_t mcp_phi0[MAX_MCPARTICLES];
  Float_t mcp_omega[MAX_MCPARTICLES];
  Float_t mcp_z0[MAX_MCPARTICLES];
  Float_t mcp_tanL[MAX_MCPARTICLES];
  Int_t mcp_nhit_vxd[MAX_MCPARTICLES];
  Int_t mcp_nhit_sit[MAX_MCPARTICLES];
  Int_t mcp_nhit_ftd[MAX_MCPARTICLES];
  Int_t mcp_nhit_tpc[MAX_MCPARTICLES];
  Int_t mcp_nhit_set[MAX_MCPARTICLES];
  Int_t mcp_nhit_etd[MAX_MCPARTICLES];
  Int_t mcp_ntrk_linked[MAX_MCPARTICLES];
  Int_t mcp_track_index_first[MAX_MCPARTICLES];
  Int_t mcp_track_index_second[MAX_MCPARTICLES];
  Int_t mcp_track_index_third[MAX_MCPARTICLES];
  Float_t mcp_track_weight_first[MAX_MCPARTICLES];
  Float_t mcp_track_weight_second[MAX_MCPARTICLES];
  Float_t mcp_track_weight_third[MAX_MCPARTICLES];
  Float_t mcp_track_weight_recp_first[MAX_MCPARTICLES];
  Float_t mcp_track_weight_recp_second[MAX_MCPARTICLES];
  Float_t mcp_track_weight_recp_third[MAX_MCPARTICLES];

  
//  Int_t mcp_link_bank_start[MAX_MCPARTICLES];
  
//  Int_t mc_track_link_bank_nrels;
//  Int_t mc_track_link_bank_track_index[MAX_MC_TRACK_RELS];
//  Float_t mc_track_link_bank_rel_weight[MAX_MC_TRACK_RELS]; 
  
  Int_t   ntracks;
  Float_t tracks_d0[MAX_NTRACKS];
  Float_t tracks_phi0[MAX_NTRACKS];
  Float_t tracks_omega[MAX_NTRACKS];
  Float_t tracks_z0[MAX_NTRACKS];
  Float_t tracks_tanL[MAX_NTRACKS];
  Float_t tracks_ref_point_x[MAX_NTRACKS];
  Float_t tracks_ref_point_y[MAX_NTRACKS];
  Float_t tracks_ref_point_z[MAX_NTRACKS];
  Int_t   tracks_type[MAX_NTRACKS];
  Float_t tracks_chi2[MAX_NTRACKS];
  Float_t tracks_prob[MAX_NTRACKS];
  Int_t   tracks_ndf[MAX_NTRACKS];
  Float_t tracks_radius_innermost_hit[MAX_NTRACKS];
  Int_t tracks_nhit_vxd[MAX_NTRACKS];
  Int_t tracks_nhit_sit[MAX_NTRACKS];
  Int_t tracks_nhit_ftd[MAX_NTRACKS];
  Int_t tracks_nhit_tpc[MAX_NTRACKS];
  Int_t tracks_nhit_set[MAX_NTRACKS];
  Int_t tracks_nhit_etd[MAX_NTRACKS];
  Float_t tracks_cov_d0d0[MAX_NTRACKS];      
  Float_t tracks_cov_phi0d0[MAX_NTRACKS];      
  Float_t tracks_cov_phi0phi0[MAX_NTRACKS];      
  Float_t tracks_cov_omegad0[MAX_NTRACKS];      
  Float_t tracks_cov_omegaphi0[MAX_NTRACKS];      
  Float_t tracks_cov_omegaomega[MAX_NTRACKS];      
  Float_t tracks_cov_z0d0[MAX_NTRACKS];      
  Float_t tracks_cov_z0phi0[MAX_NTRACKS];      
  Float_t tracks_cov_z0omega[MAX_NTRACKS];      
  Float_t tracks_cov_z0z0[MAX_NTRACKS];      
  Float_t tracks_cov_tanLd0[MAX_NTRACKS];      
  Float_t tracks_cov_tanLphi0[MAX_NTRACKS];      
  Float_t tracks_cov_tanLomega[MAX_NTRACKS];      
  Float_t tracks_cov_tanLz0[MAX_NTRACKS];      
  Float_t tracks_cov_tanLtanL[MAX_NTRACKS];
  Int_t tracks_nmcp_linked[MAX_NTRACKS];

  Int_t tracks_mcp_index_first[MAX_NTRACKS];
  Int_t tracks_mcp_index_second[MAX_NTRACKS];
  Int_t tracks_mcp_index_third[MAX_NTRACKS];
  Float_t tracks_mcp_weight_first[MAX_NTRACKS];
  Float_t tracks_mcp_weight_second[MAX_NTRACKS];
  Float_t tracks_mcp_weight_third[MAX_NTRACKS];
  Float_t tracks_mcp_weight_recp_first[MAX_NTRACKS];
  Float_t tracks_mcp_weight_recp_second[MAX_NTRACKS];
  Float_t tracks_mcp_weight_recp_third[MAX_NTRACKS];

//  Int_t tracks_link_bank_start[MAX_NTRACKS];
//  
//  Int_t track_mc_link_bank_nrels;
//  Int_t track_mc_link_bank_mc_index[MAX_TRACK_MC_RELS];
//  Float_t track_mc_link_bank_rel_weight[MAX_TRACK_MC_RELS]; 
  
  // List of branches
  
  TBranch        *b_eventNumber;    //!  
  TBranch        *b_runNumber;    //!  
  
  TBranch        *b_nmcp;    //!
  TBranch        *b_mcp_pdg;    //!  
  TBranch        *b_mcp_generator_status; //!
  TBranch        *b_mcp_simulator_status; //!
  TBranch        *b_mcp_ndaughters; //!
  TBranch        *b_mcp_nparents; //!  
  TBranch        *b_mcp_mass;    //! 
  TBranch        *b_mcp_charge;    //!
  TBranch        *b_mcp_energy;    //!
  TBranch        *b_mcp_px;    //!
  TBranch        *b_mcp_py;    //!
  TBranch        *b_mcp_pz;    //!
  TBranch        *b_mcp_vx;    //!
  TBranch        *b_mcp_vy;    //!
  TBranch        *b_mcp_vz;    //!
  TBranch        *b_mcp_d0;    //!
  TBranch        *b_mcp_phi0;    //!
  TBranch        *b_mcp_omega;    //!
  TBranch        *b_mcp_z0;    //!
  TBranch        *b_mcp_tanL;    //!
  TBranch        *b_mcp_nhit_vxd;    //!
  TBranch        *b_mcp_nhit_sit;    //!
  TBranch        *b_mcp_nhit_ftd;    //!
  TBranch        *b_mcp_nhit_tpc;    //!
  TBranch        *b_mcp_nhit_set;    //!
  TBranch        *b_mcp_nhit_etd;    //!
  TBranch        *b_mcp_ntrk_linked;    //!


  TBranch        *b_mcp_track_index_first;    //!
  TBranch        *b_mcp_track_index_second;    //!
  TBranch        *b_mcp_track_index_third;    //!
  TBranch        *b_mcp_track_weight_first;    //!
  TBranch        *b_mcp_track_weight_second;    //!
  TBranch        *b_mcp_track_weight_third;    //!
  TBranch        *b_mcp_track_weight_recp_first;    //!
  TBranch        *b_mcp_track_weight_recp_second;    //!
  TBranch        *b_mcp_track_weight_recp_third;    //!

  
//  TBranch        *b_mcp_link_bank_start;    //!
//
//  TBranch        *b_mc_track_link_bank_nrels; //!
//  TBranch        *b_mc_track_link_bank_track_index; //!
//  TBranch        *b_mc_track_link_bank_rel_weight; //!

  TBranch        *b_ntracks;     //!
  TBranch        *b_tracks_d0;     //!
  TBranch        *b_tracks_phi0;   //!
  TBranch        *b_tracks_omega;  //!
  TBranch        *b_tracks_z0;  //! 
  TBranch        *b_tracks_tanL;  //!
  TBranch        *b_tracks_ref_point_x;  //!
  TBranch        *b_tracks_ref_point_y;  //!
  TBranch        *b_tracks_ref_point_z;  //!
  TBranch        *b_tracks_type;  //!
  TBranch        *b_tracks_chi2;  //!
  TBranch        *b_tracks_prob;  //!
  TBranch        *b_tracks_ndf;  //!
  TBranch        *b_tracks_radius_innermost_hit;  //!
  TBranch        *b_tracks_nhit_vxd;  //!
  TBranch        *b_tracks_nhit_sit;  //!
  TBranch        *b_tracks_nhit_ftd;  //!
  TBranch        *b_tracks_nhit_tpc;  //!
  TBranch        *b_tracks_nhit_set;  //!
  TBranch        *b_tracks_nhit_etd;  //!
  TBranch        *b_tracks_cov_d0d0;  //! 
  TBranch        *b_tracks_cov_phi0d0;  //!      
  TBranch        *b_tracks_cov_phi0phi0;  //!   
  TBranch        *b_tracks_cov_omegad0;  //!      
  TBranch        *b_tracks_cov_omegaphi0;  //!      
  TBranch        *b_tracks_cov_omegaomega;  //!      
  TBranch        *b_tracks_cov_z0d0;  //!      
  TBranch        *b_tracks_cov_z0phi0;  //!      
  TBranch        *b_tracks_cov_z0omega;  //!      
  TBranch        *b_tracks_cov_z0z0;  //!      
  TBranch        *b_tracks_cov_tanLd0;  //!      
  TBranch        *b_tracks_cov_tanLphi0;  //!      
  TBranch        *b_tracks_cov_tanLomega;  //!      
  TBranch        *b_tracks_cov_tanLz0;  //!      
  TBranch        *b_tracks_cov_tanLtanL;  //! 
  TBranch        *b_tracks_nmcp_linked;  //! 
  
  TBranch        *b_tracks_mcp_index_first;    //!
  TBranch        *b_tracks_mcp_index_second;    //!
  TBranch        *b_tracks_mcp_index_third;    //!
  TBranch        *b_tracks_mcp_weight_first;    //!
  TBranch        *b_tracks_mcp_weight_second;    //!
  TBranch        *b_tracks_mcp_weight_third;    //!
  TBranch        *b_tracks_mcp_weight_recp_first;    //!
  TBranch        *b_tracks_mcp_weight_recp_second;    //!
  TBranch        *b_tracks_mcp_weight_recp_third;    //!

//  TBranch        *b_tracks_link_bank_start;  //! 
//  
//  TBranch        *b_track_mc_link_bank_nrels; //!
//  TBranch        *b_track_mc_link_bank_mc_index; //!
//  TBranch        *b_track_mc_link_bank_rel_weight; //!
  
  TrkAnalysisTree(TTree *tree=0);
  
  ~TrkAnalysisTree(); // no longer virtual 
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     CreateBranches(TTree *tree);
  // virtual void     Loop(); // removed as in NTuple.h
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  virtual void     Clear();

};

#endif

#ifdef TrkAnalysisTree_cxx
TrkAnalysisTree::TrkAnalysisTree(TTree *tree)
{
  if(tree) Init(tree);
}

TrkAnalysisTree::~TrkAnalysisTree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t TrkAnalysisTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t TrkAnalysisTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void TrkAnalysisTree::CreateBranches(TTree *tree)
{

  tree->Branch("nmcp", &nmcp ,"nmcp/I" );

  tree->Branch("mcp_pdg", mcp_pdg ,"mcp_pdg[nmcp]/I" );
  tree->Branch("mcp_generator_status", mcp_generator_status ,"mcp_generator_status[nmcp]/I" );
  tree->Branch("mcp_simulator_status", mcp_simulator_status ,"mcp_simulator_status[nmcp]/I" );
  tree->Branch("mcp_ndaughters", mcp_ndaughters ,"mcp_ndaughters[nmcp]/I" );
  tree->Branch("mcp_nparents", mcp_nparents ,"mcp_nparents[nmcp]/I" );
  tree->Branch("mcp_mass", mcp_mass ,"mcp_mass[nmcp]/F" );
  tree->Branch("mcp_charge", mcp_charge ,"mcp_charge[nmcp]/F" );
  tree->Branch("mcp_energy", mcp_energy ,"mcp_energy[nmcp]/F" );
  tree->Branch("mcp_px", mcp_px ,"mcp_px[nmcp]/F" );
  tree->Branch("mcp_py", mcp_py ,"mcp_py[nmcp]/F" );
  tree->Branch("mcp_pz", mcp_pz ,"mcp_pz[nmcp]/F" );
  tree->Branch("mcp_vx", mcp_vx ,"mcp_vx[nmcp]/F" );
  tree->Branch("mcp_vy", mcp_vy ,"mcp_vy[nmcp]/F" );
  tree->Branch("mcp_vz", mcp_vz ,"mcp_vz[nmcp]/F" );
  tree->Branch("mcp_d0", mcp_d0 ,"mcp_d0[nmcp]/F" );
  tree->Branch("mcp_phi0", mcp_phi0 ,"mcp_phi0[nmcp]/F" );
  tree->Branch("mcp_omega", mcp_omega ,"mcp_omega[nmcp]/F" );
  tree->Branch("mcp_z0", mcp_z0 ,"mcp_z0[nmcp]/F" );
  tree->Branch("mcp_tanL", mcp_tanL ,"mcp_tanL[nmcp]/F" );
  tree->Branch("mcp_nhit_vxd", mcp_nhit_vxd ,"mcp_nhit_vxd[nmcp]/I" );
  tree->Branch("mcp_nhit_sit", mcp_nhit_sit ,"mcp_nhit_sit[nmcp]/I" );
  tree->Branch("mcp_nhit_ftd", mcp_nhit_ftd ,"mcp_nhit_ftd[nmcp]/I" );
  tree->Branch("mcp_nhit_tpc", mcp_nhit_tpc ,"mcp_nhit_tpc[nmcp]/I" );
  tree->Branch("mcp_nhit_set", mcp_nhit_set ,"mcp_nhit_set[nmcp]/I" );
  tree->Branch("mcp_nhit_etd", mcp_nhit_etd ,"mcp_nhit_etd[nmcp]/I" );  
  tree->Branch("mcp_ntrk_linked", mcp_ntrk_linked ,"mcp_ntrk_linked[nmcp]/I" );
  
  tree->Branch("mcp_track_index_first", mcp_track_index_first ,"mcp_track_index_first[nmcp]/I" );
  tree->Branch("mcp_track_index_second", mcp_track_index_second ,"mcp_track_index_second[nmcp]/I" );
  tree->Branch("mcp_track_index_third", mcp_track_index_third ,"mcp_track_index_third[nmcp]/I" );
  tree->Branch("mcp_track_weight_first", mcp_track_weight_first ,"mcp_track_weight_first[nmcp]/F" );
  tree->Branch("mcp_track_weight_second", mcp_track_weight_second ,"mcp_track_weight_second[nmcp]/F" );
  tree->Branch("mcp_track_weight_third", mcp_track_weight_third ,"mcp_track_weight_third[nmcp]/F" );
  tree->Branch("mcp_track_weight_recp_first", mcp_track_weight_recp_first ,"mcp_track_weight_recp_first[nmcp]/F" );
  tree->Branch("mcp_track_weight_recp_second", mcp_track_weight_recp_second ,"mcp_track_weight_recp_second[nmcp]/F" );
  tree->Branch("mcp_track_weight_recp_third", mcp_track_weight_recp_third ,"mcp_track_weight_recp_third[nmcp]/F" );

  
//  tree->Branch("mcp_link_bank_start", mcp_link_bank_start ,"mcp_link_bank_start[nmcp]/I" );
  
//  tree->Branch("mc_track_link_bank_nrels", &mc_track_link_bank_nrels ,"mc_track_link_bank_nrels/I" );
//  tree->Branch("mc_track_link_bank_track_index", mc_track_link_bank_track_index ,"mc_track_link_bank_track_index[mc_track_link_bank_nrels]/I" );
//  tree->Branch("mc_track_link_bank_rel_weight", mc_track_link_bank_rel_weight ,"mc_track_link_bank_rel_weight[mc_track_link_bank_nrels]/F" );
  
  tree->Branch("ntracks", &ntracks ,"ntracks/I" );
  tree->Branch("tracks_d0", tracks_d0 ,"tracks_d0[ntracks]/F" );
  tree->Branch("tracks_phi0", tracks_phi0 ,"tracks_phi0[ntracks]/F" );
  tree->Branch("tracks_omega", tracks_omega ,"tracks_omega[ntracks]/F" );
  tree->Branch("tracks_z0", tracks_z0 ,"tracks_z0[ntracks]/F" );
  tree->Branch("tracks_tanL", tracks_tanL ,"tracks_tanL[ntracks]/F" );
  tree->Branch("tracks_ref_point_x", tracks_ref_point_x ,"tracks_ref_point_x[ntracks]/F" );
  tree->Branch("tracks_ref_point_y", tracks_ref_point_y ,"tracks_ref_point_y[ntracks]/F" );
  tree->Branch("tracks_ref_point_z", tracks_ref_point_z ,"tracks_ref_point_z[ntracks]/F" );
  tree->Branch("tracks_type", tracks_type ,"tracks_type[ntracks]/I" );  
  tree->Branch("tracks_chi2", tracks_chi2 ,"tracks_chi2[ntracks]/F" );
  tree->Branch("tracks_prob", tracks_prob ,"tracks_prob[ntracks]/F" );
  tree->Branch("tracks_ndf", tracks_ndf ,"tracks_ndf[ntracks]/I" );
  tree->Branch("tracks_radius_innermost_hit", tracks_radius_innermost_hit ,"tracks_radius_innermost_hit[ntracks]/F" );
  tree->Branch("tracks_nhit_vxd", tracks_nhit_vxd ,"tracks_nhit_vxd[ntracks]/I" );
  tree->Branch("tracks_nhit_sit", tracks_nhit_sit ,"tracks_nhit_sit[ntracks]/I" );
  tree->Branch("tracks_nhit_ftd", tracks_nhit_ftd ,"tracks_nhit_ftd[ntracks]/I" );
  tree->Branch("tracks_nhit_tpc", tracks_nhit_tpc ,"tracks_nhit_tpc[ntracks]/I" );
  tree->Branch("tracks_nhit_set", tracks_nhit_set ,"tracks_nhit_set[ntracks]/I" );
  tree->Branch("tracks_nhit_etd", tracks_nhit_etd ,"tracks_nhit_etd[ntracks]/I" );
  tree->Branch("tracks_cov_d0d0", tracks_cov_d0d0 ,"tracks_cov_d0d0[ntracks]/F" );
  tree->Branch("tracks_cov_phi0d0", tracks_cov_phi0d0 ,"tracks_cov_phi0d0[ntracks]/F" );
  tree->Branch("tracks_cov_phi0phi0", tracks_cov_phi0phi0 ,"tracks_cov_phi0phi0[ntracks]/F" );
  tree->Branch("tracks_cov_omegad0", tracks_cov_omegad0 ,"tracks_cov_omegad0[ntracks]/F" );
  tree->Branch("tracks_cov_omegaphi0", tracks_cov_omegaphi0 ,"tracks_cov_omegaphi0[ntracks]/F" );
  tree->Branch("tracks_cov_omegaomega", tracks_cov_omegaomega ,"tracks_cov_omegaomega[ntracks]/F" );
  tree->Branch("tracks_cov_z0d0", tracks_cov_z0d0 ,"tracks_cov_z0d0[ntracks]/F" );
  tree->Branch("tracks_cov_z0phi0", tracks_cov_z0phi0 ,"tracks_cov_z0phi0[ntracks]/F" );
  tree->Branch("tracks_cov_z0omega", tracks_cov_z0omega ,"tracks_cov_z0omega[ntracks]/F" );
  tree->Branch("tracks_cov_z0z0", tracks_cov_z0z0 ,"tracks_cov_z0z0[ntracks]/F" );
  tree->Branch("tracks_cov_tanLd0", tracks_cov_tanLd0 ,"tracks_cov_tanLd0[ntracks]/F" );
  tree->Branch("tracks_cov_tanLphi0", tracks_cov_tanLphi0 ,"tracks_cov_tanLphi0[ntracks]/F" );
  tree->Branch("tracks_cov_tanLomega", tracks_cov_tanLomega ,"tracks_cov_tanLomega[ntracks]/F" );
  tree->Branch("tracks_cov_tanLz0", tracks_cov_tanLz0 ,"tracks_cov_tanLz0[ntracks]/F" );
  tree->Branch("tracks_cov_tanLtanL", tracks_cov_tanLtanL ,"tracks_cov_tanLtanL[ntracks]/F" );
  tree->Branch("tracks_nmcp_linked", tracks_nmcp_linked ,"tracks_nmcp_linked[ntracks]/I" );
  
  tree->Branch("tracks_mcp_index_first", tracks_mcp_index_first ,"tracks_mcp_index_first[ntracks]/I" );
  tree->Branch("tracks_mcp_index_second", tracks_mcp_index_second ,"tracks_mcp_index_second[ntracks]/I" );
  tree->Branch("tracks_mcp_index_third", tracks_mcp_index_third ,"tracks_mcp_index_third[ntracks]/I" );  
  tree->Branch("tracks_mcp_weight_first", tracks_mcp_weight_first ,"tracks_mcp_weight_first[ntracks]/F" );
  tree->Branch("tracks_mcp_weight_second", tracks_mcp_weight_second ,"tracks_mcp_weight_second[ntracks]/F" );
  tree->Branch("tracks_mcp_weight_third", tracks_mcp_weight_third ,"tracks_mcp_weight_third[ntracks]/F" );
  tree->Branch("tracks_mcp_weight_recp_first", tracks_mcp_weight_recp_first ,"tracks_mcp_weight_recp_first[ntracks]/F" );
  tree->Branch("tracks_mcp_weight_recp_second", tracks_mcp_weight_recp_second ,"tracks_mcp_weight_recp_second[ntracks]/F" );
  tree->Branch("tracks_mcp_weight_recp_third", tracks_mcp_weight_recp_third ,"tracks_mcp_weight_recp_third[ntracks]/F" );

//  tree->Branch("tracks_link_bank_start", tracks_link_bank_start ,"tracks_link_bank_start[ntracks]/I" );
//
//  tree->Branch("track_mc_link_bank_nrels", &track_mc_link_bank_nrels ,"track_mc_link_bank_nrels/I" );
//  tree->Branch("track_mc_link_bank_mc_index", track_mc_link_bank_mc_index ,"track_mc_link_bank_mc_index[track_mc_link_bank_nrels]/I" );
//  tree->Branch("track_mc_link_bank_rel_weight", track_mc_link_bank_rel_weight ,"track_mc_link_bank_rel_weight[track_mc_link_bank_nrels]/F" );



}

void TrkAnalysisTree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  std::cout << " tree address " << tree << std::endl;
  
  // Set branch addresses and branch pointers
  if (!tree) { 
    std::cout << " tree address is NULL exit 1 called from " << __FILE__ << " line " <<  __LINE__ << std::endl;
    exit(1); 
    
  }

  tree->Print();
  
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
    
  fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber );
  fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber );
  
  fChain->SetBranchAddress("nmcp", &nmcp, &b_nmcp );

  fChain->SetBranchAddress("mcp_pdg", mcp_pdg, &b_mcp_pdg );
  fChain->SetBranchAddress("mcp_generator_status", mcp_generator_status, &b_mcp_generator_status );
  fChain->SetBranchAddress("mcp_simulator_status", mcp_simulator_status, &b_mcp_simulator_status );
  fChain->SetBranchAddress("mcp_ndaughters", mcp_ndaughters, &b_mcp_ndaughters );
  fChain->SetBranchAddress("mcp_nparents", mcp_nparents, &b_mcp_nparents );
  fChain->SetBranchAddress("mcp_mass", mcp_mass, &b_mcp_mass );
  fChain->SetBranchAddress("mcp_charge", mcp_charge, &b_mcp_charge );
  fChain->SetBranchAddress("mcp_energy", mcp_energy, &b_mcp_energy );
  fChain->SetBranchAddress("mcp_px", mcp_px, &b_mcp_px );
  fChain->SetBranchAddress("mcp_py", mcp_py, &b_mcp_py );
  fChain->SetBranchAddress("mcp_pz", mcp_pz, &b_mcp_pz );
  fChain->SetBranchAddress("mcp_vx", mcp_vx, &b_mcp_vx );
  fChain->SetBranchAddress("mcp_vy", mcp_vy, &b_mcp_vy );
  fChain->SetBranchAddress("mcp_vz", mcp_vz, &b_mcp_vz );
  fChain->SetBranchAddress("mcp_d0", mcp_d0, &b_mcp_d0 );
  fChain->SetBranchAddress("mcp_phi0", mcp_phi0, &b_mcp_phi0 );
  fChain->SetBranchAddress("mcp_omega", mcp_omega, &b_mcp_omega );
  fChain->SetBranchAddress("mcp_z0", mcp_z0, &b_mcp_z0 );
  fChain->SetBranchAddress("mcp_tanL", mcp_tanL, &b_mcp_tanL );
  fChain->SetBranchAddress("mcp_nhit_vxd", mcp_nhit_vxd, &b_mcp_nhit_vxd );
  fChain->SetBranchAddress("mcp_nhit_sit", mcp_nhit_sit, &b_mcp_nhit_sit );
  fChain->SetBranchAddress("mcp_nhit_ftd", mcp_nhit_ftd, &b_mcp_nhit_ftd );
  fChain->SetBranchAddress("mcp_nhit_tpc", mcp_nhit_tpc, &b_mcp_nhit_tpc );
  fChain->SetBranchAddress("mcp_nhit_set", mcp_nhit_set, &b_mcp_nhit_set );
  fChain->SetBranchAddress("mcp_nhit_etd", mcp_nhit_etd, &b_mcp_nhit_etd );
  fChain->SetBranchAddress("mcp_ntrk_linked", mcp_ntrk_linked, &b_mcp_ntrk_linked );
  fChain->SetBranchAddress("mcp_track_index_first", mcp_track_index_first, &b_mcp_track_index_first );
  fChain->SetBranchAddress("mcp_track_index_second", mcp_track_index_second, &b_mcp_track_index_second );
  fChain->SetBranchAddress("mcp_track_index_third", mcp_track_index_third, &b_mcp_track_index_third );
  fChain->SetBranchAddress("mcp_track_weight_first", mcp_track_weight_first, &b_mcp_track_weight_first );
  fChain->SetBranchAddress("mcp_track_weight_second", mcp_track_weight_second, &b_mcp_track_weight_second );
  fChain->SetBranchAddress("mcp_track_weight_third", mcp_track_weight_third, &b_mcp_track_weight_third );
  fChain->SetBranchAddress("mcp_track_weight_recp_first", mcp_track_weight_recp_first, &b_mcp_track_weight_recp_first );
  fChain->SetBranchAddress("mcp_track_weight_recp_second", mcp_track_weight_recp_second, &b_mcp_track_weight_recp_second );
  fChain->SetBranchAddress("mcp_track_weight_recp_third", mcp_track_weight_recp_third, &b_mcp_track_weight_recp_third );

//  fChain->SetBranchAddress("mcp_link_bank_start", mcp_link_bank_start, &b_mcp_link_bank_start );
//  
//  fChain->SetBranchAddress("mc_track_link_bank_nrels", &mc_track_link_bank_nrels, &b_mc_track_link_bank_nrels );
//  fChain->SetBranchAddress("mc_track_link_bank_track_index", mc_track_link_bank_track_index, &b_mc_track_link_bank_track_index );
//  fChain->SetBranchAddress("mc_track_link_bank_rel_weight", mc_track_link_bank_rel_weight, &b_mc_track_link_bank_rel_weight );
  
  
  fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks );
  fChain->SetBranchAddress("tracks_d0", tracks_d0, &b_tracks_d0);
  fChain->SetBranchAddress("tracks_phi0", tracks_phi0, &b_tracks_phi0 );
  fChain->SetBranchAddress("tracks_omega", tracks_omega, &b_tracks_omega);
  fChain->SetBranchAddress("tracks_z0", tracks_z0, &b_tracks_z0);
  fChain->SetBranchAddress("tracks_tanL", tracks_tanL, &b_tracks_tanL);
  fChain->SetBranchAddress("tracks_ndf", tracks_ndf, &b_tracks_ndf );
  fChain->SetBranchAddress("tracks_chi2", tracks_chi2, &b_tracks_chi2);
  fChain->SetBranchAddress("tracks_prob", tracks_prob, &b_tracks_prob);
  fChain->SetBranchAddress("tracks_type", tracks_type, &b_tracks_type);
  fChain->SetBranchAddress("tracks_ref_point_x", tracks_ref_point_x, &b_tracks_ref_point_x);
  fChain->SetBranchAddress("tracks_ref_point_y", tracks_ref_point_y, &b_tracks_ref_point_y);
  fChain->SetBranchAddress("tracks_ref_point_z", tracks_ref_point_z, &b_tracks_ref_point_z);
  fChain->SetBranchAddress("tracks_type", tracks_type ,&b_tracks_type );
  fChain->SetBranchAddress("tracks_chi2", tracks_chi2, &b_tracks_chi2 );
  fChain->SetBranchAddress("tracks_ndf", tracks_ndf, &b_tracks_ndf );
  fChain->SetBranchAddress("tracks_radius_innermost_hit", tracks_radius_innermost_hit, &b_tracks_radius_innermost_hit );
  fChain->SetBranchAddress("tracks_nhit_vxd", tracks_nhit_vxd, &b_tracks_nhit_vxd );
  fChain->SetBranchAddress("tracks_nhit_sit", tracks_nhit_sit, &b_tracks_nhit_sit );
  fChain->SetBranchAddress("tracks_nhit_ftd", tracks_nhit_ftd, &b_tracks_nhit_ftd );
  fChain->SetBranchAddress("tracks_nhit_tpc", tracks_nhit_tpc, &b_tracks_nhit_tpc );
  fChain->SetBranchAddress("tracks_nhit_set", tracks_nhit_set, &b_tracks_nhit_set );
  fChain->SetBranchAddress("tracks_nhit_etd", tracks_nhit_etd, &b_tracks_nhit_etd );
  fChain->SetBranchAddress("tracks_cov_d0d0", tracks_cov_d0d0, &b_tracks_cov_d0d0 );      
  fChain->SetBranchAddress("tracks_cov_phi0d0", tracks_cov_phi0d0, &b_tracks_cov_phi0d0 );      
  fChain->SetBranchAddress("tracks_cov_phi0phi0", tracks_cov_phi0phi0, &b_tracks_cov_phi0phi0 );      
  fChain->SetBranchAddress("tracks_cov_omegad0", tracks_cov_omegad0, &b_tracks_cov_omegad0 );
  fChain->SetBranchAddress("tracks_cov_omegaphi0", tracks_cov_omegaphi0, &b_tracks_cov_omegaphi0 );      
  fChain->SetBranchAddress("tracks_cov_omegaomega", tracks_cov_omegaomega, &b_tracks_cov_omegaomega );      
  fChain->SetBranchAddress("tracks_cov_z0d0", tracks_cov_z0d0, &b_tracks_cov_z0d0 );      
  fChain->SetBranchAddress("tracks_cov_z0phi0", tracks_cov_z0phi0, &b_tracks_cov_z0phi0 );      
  fChain->SetBranchAddress("tracks_cov_z0omega", tracks_cov_z0omega, &b_tracks_cov_z0omega );      
  fChain->SetBranchAddress("tracks_cov_z0z0", tracks_cov_z0z0, &b_tracks_cov_z0z0) ;      
  fChain->SetBranchAddress("tracks_cov_tanLd0", tracks_cov_tanLd0, &b_tracks_cov_tanLd0 );      
  fChain->SetBranchAddress("tracks_cov_tanLphi0", tracks_cov_tanLphi0, &b_tracks_cov_tanLphi0 );      
  fChain->SetBranchAddress("tracks_cov_tanLomega", tracks_cov_tanLomega, &b_tracks_cov_tanLomega );      
  fChain->SetBranchAddress("tracks_cov_tanLz0", tracks_cov_tanLz0, &b_tracks_cov_tanLz0 );      
  fChain->SetBranchAddress("tracks_cov_tanLtanL", tracks_cov_tanLtanL, &b_tracks_cov_tanLtanL );
  fChain->SetBranchAddress("tracks_nmcp_linked", tracks_nmcp_linked, &b_tracks_nmcp_linked );

  fChain->SetBranchAddress("tracks_mcp_index_first", tracks_mcp_index_first, &b_tracks_mcp_index_first );
  fChain->SetBranchAddress("tracks_mcp_index_second", tracks_mcp_index_second, &b_tracks_mcp_index_second );
  fChain->SetBranchAddress("tracks_mcp_index_third", tracks_mcp_index_third, &b_tracks_mcp_index_third );  
  fChain->SetBranchAddress("tracks_mcp_weight_first", tracks_mcp_weight_first, &b_tracks_mcp_weight_first );
  fChain->SetBranchAddress("tracks_mcp_weight_second", tracks_mcp_weight_second, &b_tracks_mcp_weight_second );
  fChain->SetBranchAddress("tracks_mcp_weight_third", tracks_mcp_weight_third, &b_tracks_mcp_weight_third );
  fChain->SetBranchAddress("tracks_mcp_weight_recp_first", tracks_mcp_weight_recp_first, &b_tracks_mcp_weight_recp_first );
  fChain->SetBranchAddress("tracks_mcp_weight_recp_second", tracks_mcp_weight_recp_second, &b_tracks_mcp_weight_recp_second );
  fChain->SetBranchAddress("tracks_mcp_weight_recp_third", tracks_mcp_weight_recp_third, &b_tracks_mcp_weight_recp_third );

//  fChain->SetBranchAddress("tracks_link_bank_start", tracks_link_bank_start, &b_tracks_link_bank_start );
//  
//  fChain->SetBranchAddress("track_mc_link_bank_nrels", &track_mc_link_bank_nrels, &b_track_mc_link_bank_nrels );
//  fChain->SetBranchAddress("track_mc_link_bank_mc_index", track_mc_link_bank_mc_index, &b_track_mc_link_bank_mc_index );
//  fChain->SetBranchAddress("track_mc_link_bank_rel_weight", track_mc_link_bank_rel_weight, &b_track_mc_link_bank_rel_weight );
  
  
  
  
  
  
  
  Notify();
}

Bool_t TrkAnalysisTree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void TrkAnalysisTree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t TrkAnalysisTree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void TrkAnalysisTree::Clear(){
  
  eventNumber = 0;
  runNumber = 0;
  
  nmcp = 0;
  
  for (int i=0; i<MAX_MCPARTICLES; ++i) {
    
    mcp_pdg[i]=0;    
    mcp_generator_status[i]=0;
    mcp_simulator_status[i]=0;
    mcp_ndaughters[i]=0;
    mcp_nparents[i]=0;      
    mcp_mass[i]=0;
    mcp_charge[i]=0;
    mcp_energy[i]=0;
    mcp_px[i]=0;
    mcp_py[i]=0;
    mcp_pz[i]=0;
    mcp_vx[i]=0;
    mcp_vy[i]=0;
    mcp_vz[i]=0;
    mcp_d0[i]=0;
    mcp_phi0[i]=0;
    mcp_omega[i]=0;
    mcp_z0[i]=0;
    mcp_tanL[i]=0;
    mcp_nhit_vxd[i]=0;
    mcp_nhit_sit[i]=0;
    mcp_nhit_ftd[i]=0;
    mcp_nhit_tpc[i]=0;
    mcp_nhit_set[i]=0;
    mcp_nhit_etd[i]=0;
    mcp_ntrk_linked[i]=0;
    mcp_track_index_first[i]=0;
    mcp_track_index_second[i]=0;
    mcp_track_index_third[i]=0;
    mcp_track_weight_first[i]=0;
    mcp_track_weight_second[i]=0;
    mcp_track_weight_third[i]=0;
    mcp_track_weight_recp_first[i]=0;
    mcp_track_weight_recp_second[i]=0;
    mcp_track_weight_recp_third[i]=0;
    
    //    mcp_link_bank_start[i]=0;
    
  }
  
//  mc_track_link_bank_nrels = 0;
//  
//  for (int i=0; i<MAX_MC_TRACK_RELS; ++i) {
//    
//
//    mc_track_link_bank_track_index[i] = 0;
//    mc_track_link_bank_rel_weight[i] = 0.0;
//    
//  }
  
  ntracks = 0;
  
  for (int i=0; i<MAX_NTRACKS; ++i) {
    
    tracks_d0[i]=0;
    tracks_phi0[i]=0;
    tracks_omega[i]=0;
    tracks_z0[i]=0;
    tracks_tanL[i]=0;
    tracks_ref_point_x[i]=0;
    tracks_ref_point_y[i]=0;
    tracks_ref_point_z[i]=0;
    tracks_type[i]=0;
    tracks_chi2[i]=0;
    tracks_prob[i]=0;
    tracks_ndf[i]=0;
    tracks_radius_innermost_hit[i]=0;
    tracks_nhit_vxd[i]=0;
    tracks_nhit_sit[i]=0;
    tracks_nhit_ftd[i]=0;
    tracks_nhit_tpc[i]=0;
    tracks_nhit_set[i]=0;
    tracks_nhit_etd[i]=0;
    tracks_cov_d0d0[i]=0;      
    tracks_cov_phi0d0[i]=0;      
    tracks_cov_phi0phi0[i]=0;      
    tracks_cov_omegad0[i]=0;      
    tracks_cov_omegaphi0[i]=0;      
    tracks_cov_omegaomega[i]=0;      
    tracks_cov_z0d0[i]=0;      
    tracks_cov_z0phi0[i]=0;      
    tracks_cov_z0omega[i]=0;      
    tracks_cov_z0z0[i]=0;      
    tracks_cov_tanLd0[i]=0;      
    tracks_cov_tanLphi0[i]=0;      
    tracks_cov_tanLomega[i]=0;      
    tracks_cov_tanLz0[i]=0;      
    tracks_cov_tanLtanL[i]=0;
    tracks_nmcp_linked[i]=0;
    tracks_mcp_index_first[i]=0;
    tracks_mcp_index_second[i]=0;
    tracks_mcp_index_third[i]=0;
    tracks_mcp_weight_first[i]=0;
    tracks_mcp_weight_second[i]=0;
    tracks_mcp_weight_third[i]=0;
    tracks_mcp_weight_recp_first[i]=0;
    tracks_mcp_weight_recp_second[i]=0;
    tracks_mcp_weight_recp_third[i]=0;
//    tracks_link_bank_start[i]=0;
    
  }  
  

//  track_mc_link_bank_nrels=0;
//  
//  for (int i=0; i<MAX_TRACK_MC_RELS; ++i) {
//    
//    track_mc_link_bank_mc_index[i]=0;
//    track_mc_link_bank_rel_weight[i]=0; 
//    
//  }
  
  
}

#endif // #ifdef TrkAnalysisTree_cxx
