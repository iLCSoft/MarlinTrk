//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 14 12:04:41 2011 by ROOT version 5.30/04
// from TTree truth_testing/truth_testing
// found on file: MarlinTrkNtuple.root
//////////////////////////////////////////////////////////

#ifndef MarlinTrkNtuple_ROOT_h
#define MarlinTrkNtuple_ROOT_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

#define MAX_SITES 300

class MarlinTrkNtuple {
public:
  TTree* fChain;  //! pointer to the analyzed TTree or TChain
  Int_t fCurrent; //! current Tree number in a TChain

  // Declaration of leaf types

  Int_t error_code;
  Int_t nsites;
  Int_t nsites_vxd;
  Int_t nsites_sit;
  Int_t nsites_ftd;
  Int_t nsites_tpc;
  Int_t nsites_set;
  Float_t x_mcp;
  Float_t y_mcp;
  Float_t z_mcp;
  Float_t px_mcp;
  Float_t py_mcp;
  Float_t pz_mcp;
  Float_t p_mcp;
  Float_t theta_mcp;
  Float_t phi_mcp;
  Int_t pdg_mcp;
  Float_t d0_mcp;
  Float_t phi0_mcp;
  Float_t omega_mcp;
  Float_t z0_mcp;
  Float_t tanL_mcp;
  Float_t d0_seed;
  Float_t phi0_seed;
  Float_t omega_seed;
  Float_t z0_seed;
  Float_t tanL_seed;
  Float_t seed_ref_point_x;
  Float_t seed_ref_point_y;
  Float_t seed_ref_point_z;
  Float_t cov_seed_d0d0;
  Float_t cov_seed_phi0d0;
  Float_t cov_seed_phi0phi0;
  Float_t cov_seed_kappad0;
  Float_t cov_seed_kappaphi0;
  Float_t cov_seed_kappakappa;
  Float_t cov_seed_z0d0;
  Float_t cov_seed_z0phi0;
  Float_t cov_seed_z0kappa;
  Float_t cov_seed_z0z0;
  Float_t cov_seed_tanLd0;
  Float_t cov_seed_tanLphi0;
  Float_t cov_seed_tanLkappa;
  Float_t cov_seed_tanLz0;
  Float_t cov_seed_tanLtanL;
  Int_t ndf;
  Float_t chi2;
  Float_t prob;
  Float_t d0_ip;
  Float_t phi0_ip;
  Float_t omega_ip;
  Float_t z0_ip;
  Float_t tanL_ip;
  Float_t cov_ip_d0d0;
  Float_t cov_ip_phi0d0;
  Float_t cov_ip_phi0phi0;
  Float_t cov_ip_omegad0;
  Float_t cov_ip_omegaphi0;
  Float_t cov_ip_omegaomega;
  Float_t cov_ip_z0d0;
  Float_t cov_ip_z0phi0;
  Float_t cov_ip_z0omega;
  Float_t cov_ip_z0z0;
  Float_t cov_ip_tanLd0;
  Float_t cov_ip_tanLphi0;
  Float_t cov_ip_tanLomega;
  Float_t cov_ip_tanLz0;
  Float_t cov_ip_tanLtanL;
  Int_t CellID0[MAX_SITES];                    //[nsites]
  Int_t rejected[MAX_SITES];                   //[nsites]
  Float_t site_x[MAX_SITES];                   //[nsites]
  Float_t site_y[MAX_SITES];                   //[nsites]
  Float_t site_z[MAX_SITES];                   //[nsites]
  Float_t ref_point_x[MAX_SITES];              //[nsites]
  Float_t ref_point_y[MAX_SITES];              //[nsites]
  Float_t ref_point_z[MAX_SITES];              //[nsites]
  Float_t d0_mc[MAX_SITES];                    //[nsites]
  Float_t phi0_mc[MAX_SITES];                  //[nsites]
  Float_t omega_mc[MAX_SITES];                 //[nsites]
  Float_t z0_mc[MAX_SITES];                    //[nsites]
  Float_t tanL_mc[MAX_SITES];                  //[nsites]
  Float_t d0_predicted[MAX_SITES];             //[nsites]
  Float_t phi0_predicted[MAX_SITES];           //[nsites]
  Float_t omega_predicted[MAX_SITES];          //[nsites]
  Float_t z0_predicted[MAX_SITES];             //[nsites]
  Float_t tanL_predicted[MAX_SITES];           //[nsites]
  Float_t d0_filtered[MAX_SITES];              //[nsites]
  Float_t phi0_filtered[MAX_SITES];            //[nsites]
  Float_t omega_filtered[MAX_SITES];           //[nsites]
  Float_t z0_filtered[MAX_SITES];              //[nsites]
  Float_t tanL_filtered[MAX_SITES];            //[nsites]
  Float_t d0_smoothed[MAX_SITES];              //[nsites]
  Float_t phi0_smoothed[MAX_SITES];            //[nsites]
  Float_t omega_smoothed[MAX_SITES];           //[nsites]
  Float_t z0_smoothed[MAX_SITES];              //[nsites]
  Float_t tanL_smoothed[MAX_SITES];            //[nsites]
  Float_t chi2_inc_filtered[MAX_SITES];        //[nsites]
  Float_t chi2_inc_smoothed[MAX_SITES];        //[nsites]
  Int_t dim[MAX_SITES];                        //[nsites]
  Float_t cov_smoothed_d0d0[MAX_SITES];        //[nsites]
  Float_t cov_smoothed_phi0d0[MAX_SITES];      //[nsites]
  Float_t cov_smoothed_phi0phi0[MAX_SITES];    //[nsites]
  Float_t cov_smoothed_omegad0[MAX_SITES];     //[nsites]
  Float_t cov_smoothed_omegaphi0[MAX_SITES];   //[nsites]
  Float_t cov_smoothed_omegaomega[MAX_SITES];  //[nsites]
  Float_t cov_smoothed_z0d0[MAX_SITES];        //[nsites]
  Float_t cov_smoothed_z0phi0[MAX_SITES];      //[nsites]
  Float_t cov_smoothed_z0omega[MAX_SITES];     //[nsites]
  Float_t cov_smoothed_z0z0[MAX_SITES];        //[nsites]
  Float_t cov_smoothed_tanLd0[MAX_SITES];      //[nsites]
  Float_t cov_smoothed_tanLphi0[MAX_SITES];    //[nsites]
  Float_t cov_smoothed_tanLomega[MAX_SITES];   //[nsites]
  Float_t cov_smoothed_tanLz0[MAX_SITES];      //[nsites]
  Float_t cov_smoothed_tanLtanL[MAX_SITES];    //[nsites]
  Float_t cov_predicted_d0d0[MAX_SITES];       //[nsites]
  Float_t cov_predicted_phi0d0[MAX_SITES];     //[nsites]
  Float_t cov_predicted_phi0phi0[MAX_SITES];   //[nsites]
  Float_t cov_predicted_omegad0[MAX_SITES];    //[nsites]
  Float_t cov_predicted_omegaphi0[MAX_SITES];  //[nsites]
  Float_t cov_predicted_omegaomega[MAX_SITES]; //[nsites]
  Float_t cov_predicted_z0d0[MAX_SITES];       //[nsites]
  Float_t cov_predicted_z0phi0[MAX_SITES];     //[nsites]
  Float_t cov_predicted_z0omega[MAX_SITES];    //[nsites]
  Float_t cov_predicted_z0z0[MAX_SITES];       //[nsites]
  Float_t cov_predicted_tanLd0[MAX_SITES];     //[nsites]
  Float_t cov_predicted_tanLphi0[MAX_SITES];   //[nsites]
  Float_t cov_predicted_tanLomega[MAX_SITES];  //[nsites]
  Float_t cov_predicted_tanLz0[MAX_SITES];     //[nsites]
  Float_t cov_predicted_tanLtanL[MAX_SITES];   //[nsites]
  Float_t cov_filtered_d0d0[MAX_SITES];        //[nsites]
  Float_t cov_filtered_phi0d0[MAX_SITES];      //[nsites]
  Float_t cov_filtered_phi0phi0[MAX_SITES];    //[nsites]
  Float_t cov_filtered_omegad0[MAX_SITES];     //[nsites]
  Float_t cov_filtered_omegaphi0[MAX_SITES];   //[nsites]
  Float_t cov_filtered_omegaomega[MAX_SITES];  //[nsites]
  Float_t cov_filtered_z0d0[MAX_SITES];        //[nsites]
  Float_t cov_filtered_z0phi0[MAX_SITES];      //[nsites]
  Float_t cov_filtered_z0omega[MAX_SITES];     //[nsites]
  Float_t cov_filtered_z0z0[MAX_SITES];        //[nsites]
  Float_t cov_filtered_tanLd0[MAX_SITES];      //[nsites]
  Float_t cov_filtered_tanLphi0[MAX_SITES];    //[nsites]
  Float_t cov_filtered_tanLomega[MAX_SITES];   //[nsites]
  Float_t cov_filtered_tanLz0[MAX_SITES];      //[nsites]
  Float_t cov_filtered_tanLtanL[MAX_SITES];    //[nsites]

  // List of branches
  TBranch* b_error_code;               //!
  TBranch* b_nsites;                   //!
  TBranch* b_nsites_vxd;               //!
  TBranch* b_nsites_sit;               //!
  TBranch* b_nsites_ftd;               //!
  TBranch* b_nsites_tpc;               //!
  TBranch* b_nsites_set;               //!
  TBranch* b_x_mcp;                    //!
  TBranch* b_y_mcp;                    //!
  TBranch* b_z_mcp;                    //!
  TBranch* b_px_mcp;                   //!
  TBranch* b_py_mcp;                   //!
  TBranch* b_pz_mcp;                   //!
  TBranch* b_p_mcp;                    //!
  TBranch* b_theta_mcp;                //!
  TBranch* b_phi_mcp;                  //!
  TBranch* b_pdg_mcp;                  //!
  TBranch* b_d0_mcp;                   //!
  TBranch* b_phi0_mcp;                 //!
  TBranch* b_omega_mcp;                //!
  TBranch* b_z0_mcp;                   //!
  TBranch* b_tanL_mcp;                 //!
  TBranch* b_d0_seed;                  //!
  TBranch* b_phi0_seed;                //!
  TBranch* b_omega_seed;               //!
  TBranch* b_z0_seed;                  //!
  TBranch* b_tanL_seed;                //!
  TBranch* b_seed_ref_point_x;         //!
  TBranch* b_seed_ref_point_y;         //!
  TBranch* b_seed_ref_point_z;         //!
  TBranch* b_cov_seed_d0d0;            //!
  TBranch* b_cov_seed_phi0d0;          //!
  TBranch* b_cov_seed_phi0phi0;        //!
  TBranch* b_cov_seed_kappad0;         //!
  TBranch* b_cov_seed_kappaphi0;       //!
  TBranch* b_cov_seed_kappakappa;      //!
  TBranch* b_cov_seed_z0d0;            //!
  TBranch* b_cov_seed_z0phi0;          //!
  TBranch* b_cov_seed_z0kappa;         //!
  TBranch* b_cov_seed_z0z0;            //!
  TBranch* b_cov_seed_tanLd0;          //!
  TBranch* b_cov_seed_tanLphi0;        //!
  TBranch* b_cov_seed_tanLkappa;       //!
  TBranch* b_cov_seed_tanLz0;          //!
  TBranch* b_cov_seed_tanLtanL;        //!
  TBranch* b_d0_ip;                    //!
  TBranch* b_phi0_ip;                  //!
  TBranch* b_omega_ip;                 //!
  TBranch* b_z0_ip;                    //!
  TBranch* b_tanL_ip;                  //!
  TBranch* b_cov_ip_d0d0;              //!
  TBranch* b_cov_ip_phi0d0;            //!
  TBranch* b_cov_ip_phi0phi0;          //!
  TBranch* b_cov_ip_omegad0;           //!
  TBranch* b_cov_ip_omegaphi0;         //!
  TBranch* b_cov_ip_omegaomega;        //!
  TBranch* b_cov_ip_z0d0;              //!
  TBranch* b_cov_ip_z0phi0;            //!
  TBranch* b_cov_ip_z0omega;           //!
  TBranch* b_cov_ip_z0z0;              //!
  TBranch* b_cov_ip_tanLd0;            //!
  TBranch* b_cov_ip_tanLphi0;          //!
  TBranch* b_cov_ip_tanLomega;         //!
  TBranch* b_cov_ip_tanLz0;            //!
  TBranch* b_cov_ip_tanLtanL;          //!
  TBranch* b_ndf;                      //!
  TBranch* b_chi2;                     //!
  TBranch* b_prob;                     //!
  TBranch* b_CellID0;                  //!
  TBranch* b_rejected;                 //!
  TBranch* b_site_x;                   //!
  TBranch* b_site_y;                   //!
  TBranch* b_site_z;                   //!
  TBranch* b_ref_point_x;              //!
  TBranch* b_ref_point_y;              //!
  TBranch* b_ref_point_z;              //!
  TBranch* b_d0_mc;                    //!
  TBranch* b_phi0_mc;                  //!
  TBranch* b_omega_mc;                 //!
  TBranch* b_z0_mc;                    //!
  TBranch* b_tanL_mc;                  //!
  TBranch* b_d0_predicted;             //!
  TBranch* b_phi0_predicted;           //!
  TBranch* b_omega_predicted;          //!
  TBranch* b_z0_predicted;             //!
  TBranch* b_tanL_predicted;           //!
  TBranch* b_d0_filtered;              //!
  TBranch* b_phi0_filtered;            //!
  TBranch* b_omega_filtered;           //!
  TBranch* b_z0_filtered;              //!
  TBranch* b_tanL_filtered;            //!
  TBranch* b_d0_smoothed;              //!
  TBranch* b_phi0_smoothed;            //!
  TBranch* b_omega_smoothed;           //!
  TBranch* b_z0_smoothed;              //!
  TBranch* b_tanL_smoothed;            //!
  TBranch* b_chi2_inc_filtered;        //!
  TBranch* b_chi2_inc_smoothed;        //!
  TBranch* b_dim;                      //!
  TBranch* b_cov_smoothed_d0d0;        //!
  TBranch* b_cov_smoothed_phi0d0;      //!
  TBranch* b_cov_smoothed_phi0phi0;    //!
  TBranch* b_cov_smoothed_omegad0;     //!
  TBranch* b_cov_smoothed_omegaphi0;   //!
  TBranch* b_cov_smoothed_omegaomega;  //!
  TBranch* b_cov_smoothed_z0d0;        //!
  TBranch* b_cov_smoothed_z0phi0;      //!
  TBranch* b_cov_smoothed_z0omega;     //!
  TBranch* b_cov_smoothed_z0z0;        //!
  TBranch* b_cov_smoothed_tanLd0;      //!
  TBranch* b_cov_smoothed_tanLphi0;    //!
  TBranch* b_cov_smoothed_tanLomega;   //!
  TBranch* b_cov_smoothed_tanLz0;      //!
  TBranch* b_cov_smoothed_tanLtanL;    //!
  TBranch* b_cov_predicted_d0d0;       //!
  TBranch* b_cov_predicted_phi0d0;     //!
  TBranch* b_cov_predicted_phi0phi0;   //!
  TBranch* b_cov_predicted_omegad0;    //!
  TBranch* b_cov_predicted_omegaphi0;  //!
  TBranch* b_cov_predicted_omegaomega; //!
  TBranch* b_cov_predicted_z0d0;       //!
  TBranch* b_cov_predicted_z0phi0;     //!
  TBranch* b_cov_predicted_z0omega;    //!
  TBranch* b_cov_predicted_z0z0;       //!
  TBranch* b_cov_predicted_tanLd0;     //!
  TBranch* b_cov_predicted_tanLphi0;   //!
  TBranch* b_cov_predicted_tanLomega;  //!
  TBranch* b_cov_predicted_tanLz0;     //!
  TBranch* b_cov_predicted_tanLtanL;   //!
  TBranch* b_cov_filtered_d0d0;        //!
  TBranch* b_cov_filtered_phi0d0;      //!
  TBranch* b_cov_filtered_phi0phi0;    //!
  TBranch* b_cov_filtered_omegad0;     //!
  TBranch* b_cov_filtered_omegaphi0;   //!
  TBranch* b_cov_filtered_omegaomega;  //!
  TBranch* b_cov_filtered_z0d0;        //!
  TBranch* b_cov_filtered_z0phi0;      //!
  TBranch* b_cov_filtered_z0omega;     //!
  TBranch* b_cov_filtered_z0z0;        //!
  TBranch* b_cov_filtered_tanLd0;      //!
  TBranch* b_cov_filtered_tanLphi0;    //!
  TBranch* b_cov_filtered_tanLomega;   //!
  TBranch* b_cov_filtered_tanLz0;      //!
  TBranch* b_cov_filtered_tanLtanL;    //!

  MarlinTrkNtuple(TTree* tree = 0);
  virtual ~MarlinTrkNtuple();
  virtual Int_t Cut(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree* tree);
  virtual void CreateBranches(TTree* tree);
  //   virtual void     Loop(); // removed as in NTuple.h
  virtual Bool_t Notify();
  virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef MarlinTrkNtuple_cxx
MarlinTrkNtuple::MarlinTrkNtuple(TTree* tree) { Init(tree); }

MarlinTrkNtuple::~MarlinTrkNtuple() {
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
}

Int_t MarlinTrkNtuple::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain)
    return 0;
  return fChain->GetEntry(entry);
}
Long64_t MarlinTrkNtuple::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void MarlinTrkNtuple::CreateBranches(TTree* tree) {

  tree->Branch("error_code", &error_code, "error_code/I");

  tree->Branch("nsites", &nsites, "nsites/I");

  tree->Branch("nsites_vxd", &nsites_vxd, "nsites_vxd/I");
  tree->Branch("nsites_sit", &nsites_sit, "nsites_sit§/I");
  tree->Branch("nsites_ftd", &nsites_ftd, "nsites_ftd/I");
  tree->Branch("nsites_tpc", &nsites_tpc, "nsites_tpc/I");
  tree->Branch("nsites_set", &nsites_set, "nsites_set/I");

  tree->Branch("x_mcp", &x_mcp, "x_mcp/F");
  tree->Branch("y_mcp", &y_mcp, "y_mcp/F");
  tree->Branch("z_mcp", &z_mcp, "z_mcp/F");

  tree->Branch("px_mcp", &px_mcp, "px_mcp/F");
  tree->Branch("py_mcp", &py_mcp, "py_mcp/F");
  tree->Branch("pz_mcp", &pz_mcp, "pz_mcp/F");
  tree->Branch("p_mcp", &p_mcp, "p_mcp/F");
  tree->Branch("theta_mcp", &theta_mcp, "theta_mcp/F");
  tree->Branch("phi_mcp", &phi_mcp, "phi_mcp/F");
  tree->Branch("pdg_mcp", &pdg_mcp, "pdg_mcp/I");

  tree->Branch("d0_mcp", &d0_mcp, "d0_mcp/F");
  tree->Branch("phi0_mcp", &phi0_mcp, "phi0_mcp/F");
  tree->Branch("omega_mcp", &omega_mcp, "omega_mcp/F");
  tree->Branch("z0_mcp", &z0_mcp, "z0_mcp/F");
  tree->Branch("tanL_mcp", &tanL_mcp, "tanL_mcp/F");

  tree->Branch("d0_seed", &d0_seed, "d0_seed/F");
  tree->Branch("phi0_seed", &phi0_seed, "phi0_seed/F");
  tree->Branch("omega_seed", &omega_seed, "omega_seed/F");
  tree->Branch("z0_seed", &z0_seed, "z0_seed/F");
  tree->Branch("tanL_seed", &tanL_seed, "tanL_seed/F");

  tree->Branch("seed_ref_point_x", &seed_ref_point_x, "seed_ref_point_x/F");
  tree->Branch("seed_ref_point_y", &seed_ref_point_y, "seed_ref_point_y/F");
  tree->Branch("seed_ref_point_z", &seed_ref_point_z, "seed_ref_point_z/F");

  tree->Branch("cov_seed_d0d0", &cov_seed_d0d0, "cov_seed_d0d0/F");
  tree->Branch("cov_seed_phi0d0", &cov_seed_phi0d0, "cov_seed_phi0d0/F");
  tree->Branch("cov_seed_phi0phi0", &cov_seed_phi0phi0, "cov_seed_phi0phi0/F");
  tree->Branch("cov_seed_kappad0", &cov_seed_kappad0, "cov_seed_kappad0/F");
  tree->Branch("cov_seed_kappaphi0", &cov_seed_kappaphi0, "cov_seed_kappaphi0/F");
  tree->Branch("cov_seed_kappakappa", &cov_seed_kappakappa, "cov_seed_kappakappa/F");
  tree->Branch("cov_seed_z0phi0", &cov_seed_z0phi0, "cov_seed_z0phi0/F");
  tree->Branch("cov_seed_z0kappa", &cov_seed_z0kappa, "cov_seed_z0kappa/F");
  tree->Branch("cov_seed_z0z0", &cov_seed_z0z0, "cov_seed_z0z0/F");
  tree->Branch("cov_seed_tanLd0", &cov_seed_tanLd0, "cov_seed_tanLd0/F");
  tree->Branch("cov_seed_tanLphi0", &cov_seed_tanLphi0, "cov_seed_tanLphi0/F");
  tree->Branch("cov_seed_tanLkappa", &cov_seed_tanLkappa, "cov_seed_tanLkappa/F");
  tree->Branch("cov_seed_tanLz0", &cov_seed_tanLz0, "cov_seed_tanLz0/F");
  tree->Branch("cov_seed_tanLtanL", &cov_seed_tanLtanL, "cov_seed_tanLtanL/F");

  tree->Branch("d0_ip", &d0_ip, "d0_ip/F");
  tree->Branch("phi0_ip", &phi0_ip, "phi0_ip/F");
  tree->Branch("omega_ip", &omega_ip, "omega_ip/F");
  tree->Branch("z0_ip", &z0_ip, "z0_ip/F");
  tree->Branch("tanL_ip", &tanL_ip, "tanL_ip/F");

  tree->Branch("cov_ip_d0d0", &cov_ip_d0d0, "cov_ip_d0d0/F");
  tree->Branch("cov_ip_phi0d0", &cov_ip_phi0d0, "cov_ip_phi0d0/F");
  tree->Branch("cov_ip_phi0phi0", &cov_ip_phi0phi0, "cov_ip_phi0phi0/F");
  tree->Branch("cov_ip_omegad0", &cov_ip_omegad0, "cov_ip_omegad0/F");
  tree->Branch("cov_ip_omegaphi0", &cov_ip_omegaphi0, "cov_ip_omegaphi0/F");
  tree->Branch("cov_ip_omegaomega", &cov_ip_omegaomega, "cov_ip_omegaomega/F");
  tree->Branch("cov_ip_z0phi0", &cov_ip_z0phi0, "cov_ip_z0phi0/F");
  tree->Branch("cov_ip_z0omega", &cov_ip_z0omega, "cov_ip_z0omega/F");
  tree->Branch("cov_ip_z0z0", &cov_ip_z0z0, "cov_ip_z0z0/F");
  tree->Branch("cov_ip_tanLd0", &cov_ip_tanLd0, "cov_ip_tanLd0/F");
  tree->Branch("cov_ip_tanLphi0", &cov_ip_tanLphi0, "cov_ip_tanLphi0/F");
  tree->Branch("cov_ip_tanLomega", &cov_ip_tanLomega, "cov_ip_tanLomega/F");
  tree->Branch("cov_ip_tanLz0", &cov_ip_tanLz0, "cov_ip_tanLz0/F");
  tree->Branch("cov_ip_tanLtanL", &cov_ip_tanLtanL, "cov_ip_tanLtanL/F");

  tree->Branch("ndf", &ndf, "ndf/I");
  tree->Branch("chi2", &chi2, "chi2/F");
  tree->Branch("prob", &prob, "prob/F");

  tree->Branch("CellID0", CellID0, "CellID0[nsites]/I");
  tree->Branch("rejected", rejected, "rejected[nsites]/I");

  tree->Branch("site_x", site_x, "site_x[nsites]/F");
  tree->Branch("site_y", site_y, "site_y[nsites]/F");
  tree->Branch("site_z", site_z, "site_z[nsites]/F");

  tree->Branch("ref_point_x", ref_point_x, "ref_point_x[nsites]/F");
  tree->Branch("ref_point_y", ref_point_y, "ref_point_y[nsites]/F");
  tree->Branch("ref_point_z", ref_point_z, "ref_point_z[nsites]/F");

  tree->Branch("d0_mc", d0_mc, "d0_mc[nsites]/F");
  tree->Branch("phi0_mc", phi0_mc, "phi0_mc[nsites]/F");
  tree->Branch("omega_mc", omega_mc, "omega_mc[nsites]/F");
  tree->Branch("z0_mc", z0_mc, "z0_mc[nsites]/F");
  tree->Branch("tanL_mc", tanL_mc, "tanL_mc[nsites]/F");

  tree->Branch("d0_predicted", d0_predicted, "d0_predicted[nsites]/F");
  tree->Branch("phi0_predicted", phi0_predicted, "phi0_predicted[nsites]/F");
  tree->Branch("omega_predicted", omega_predicted, "omega_predicted[nsites]/F");
  tree->Branch("z0_predicted", z0_predicted, "z0_predicted[nsites]/F");
  tree->Branch("tanL_predicted", tanL_predicted, "tanL_predicted[nsites]/F");

  tree->Branch("d0_filtered", d0_filtered, "d0_filtered[nsites]/F");
  tree->Branch("phi0_filtered", phi0_filtered, "phi0_filtered[nsites]/F");
  tree->Branch("omega_filtered", omega_filtered, "omega_filtered[nsites]/F");
  tree->Branch("z0_filtered", z0_filtered, "z0_filtered[nsites]/F");
  tree->Branch("tanL_filtered", tanL_filtered, "tanL_filtered[nsites]/F");

  tree->Branch("d0_smoothed", d0_smoothed, "d0_smoothed[nsites]/F");
  tree->Branch("phi0_smoothed", phi0_smoothed, "phi0_smoothed[nsites]/F");
  tree->Branch("omega_smoothed", omega_smoothed, "omega_smoothed[nsites]/F");
  tree->Branch("z0_smoothed", z0_smoothed, "z0_smoothed[nsites]/F");
  tree->Branch("tanL_smoothed", tanL_smoothed, "tanL_smoothed[nsites]/F");

  tree->Branch("chi2_inc_filtered", chi2_inc_filtered, "chi2_inc_filtered[nsites]/F");
  tree->Branch("chi2_inc_smoothed", chi2_inc_smoothed, "chi2_inc_smoothed[nsites]/F");

  tree->Branch("dim", dim, "dim[nsites]/I");

  tree->Branch("cov_smoothed_d0d0", cov_smoothed_d0d0, "cov_smoothed_d0d0[nsites]/F");
  tree->Branch("cov_smoothed_phi0d0", cov_smoothed_phi0d0, "cov_smoothed_phi0d0[nsites]/F");
  tree->Branch("cov_smoothed_phi0phi0", cov_smoothed_phi0phi0, "cov_smoothed_phi0phi0[nsites]/F");
  tree->Branch("cov_smoothed_omegad0", cov_smoothed_omegad0, "cov_smoothed_omegad0[nsites]/F");
  tree->Branch("cov_smoothed_omegaphi0", cov_smoothed_omegaphi0, "cov_smoothed_omegaphi0[nsites]/F");
  tree->Branch("cov_smoothed_omegaomega", cov_smoothed_omegaomega, "cov_smoothed_omegaomega[nsites]/F");
  tree->Branch("cov_smoothed_z0phi0", cov_smoothed_z0phi0, "cov_smoothed_z0phi0[nsites]/F");
  tree->Branch("cov_smoothed_z0omega", cov_smoothed_z0omega, "cov_smoothed_z0omega[nsites]/F");
  tree->Branch("cov_smoothed_z0z0", cov_smoothed_z0z0, "cov_smoothed_z0z0[nsites]/F");
  tree->Branch("cov_smoothed_tanLd0", cov_smoothed_tanLd0, "cov_smoothed_tanLd0[nsites]/F");
  tree->Branch("cov_smoothed_tanLphi0", cov_smoothed_tanLphi0, "cov_smoothed_tanLphi0[nsites]/F");
  tree->Branch("cov_smoothed_tanLomega", cov_smoothed_tanLomega, "cov_smoothed_tanLomega[nsites]/F");
  tree->Branch("cov_smoothed_tanLz0", cov_smoothed_tanLz0, "cov_smoothed_tanLz0[nsites]/F");
  tree->Branch("cov_smoothed_tanLtanL", cov_smoothed_tanLtanL, "cov_smoothed_tanLtanL[nsites]/F");

  tree->Branch("cov_predicted_d0d0", cov_predicted_d0d0, "cov_predicted_d0d0[nsites]/F");
  tree->Branch("cov_predicted_phi0d0", cov_predicted_phi0d0, "cov_predicted_phi0d0[nsites]/F");
  tree->Branch("cov_predicted_phi0phi0", cov_predicted_phi0phi0, "cov_predicted_phi0phi0[nsites]/F");
  tree->Branch("cov_predicted_omegad0", cov_predicted_omegad0, "cov_predicted_omegad0[nsites]/F");
  tree->Branch("cov_predicted_omegaphi0", cov_predicted_omegaphi0, "cov_predicted_omegaphi0[nsites]/F");
  tree->Branch("cov_predicted_omegaomega", cov_predicted_omegaomega, "cov_predicted_omegaomega[nsites]/F");
  tree->Branch("cov_predicted_z0phi0", cov_predicted_z0phi0, "cov_predicted_z0phi0[nsites]/F");
  tree->Branch("cov_predicted_z0omega", cov_predicted_z0omega, "cov_predicted_z0omega[nsites]/F");
  tree->Branch("cov_predicted_z0z0", cov_predicted_z0z0, "cov_predicted_z0z0[nsites]/F");
  tree->Branch("cov_predicted_tanLd0", cov_predicted_tanLd0, "cov_predicted_tanLd0[nsites]/F");
  tree->Branch("cov_predicted_tanLphi0", cov_predicted_tanLphi0, "cov_predicted_tanLphi0[nsites]/F");
  tree->Branch("cov_predicted_tanLomega", cov_predicted_tanLomega, "cov_predicted_tanLomega[nsites]/F");
  tree->Branch("cov_predicted_tanLz0", cov_predicted_tanLz0, "cov_predicted_tanLz0[nsites]/F");
  tree->Branch("cov_predicted_tanLtanL", cov_predicted_tanLtanL, "cov_predicted_tanLtanL[nsites]/F");

  tree->Branch("cov_filtered_d0d0", cov_filtered_d0d0, "cov_filtered_d0d0[nsites]/F");
  tree->Branch("cov_filtered_phi0d0", cov_filtered_phi0d0, "cov_filtered_phi0d0[nsites]/F");
  tree->Branch("cov_filtered_phi0phi0", cov_filtered_phi0phi0, "cov_filtered_phi0phi0[nsites]/F");
  tree->Branch("cov_filtered_omegad0", cov_filtered_omegad0, "cov_filtered_omegad0[nsites]/F");
  tree->Branch("cov_filtered_omegaphi0", cov_filtered_omegaphi0, "cov_filtered_omegaphi0[nsites]/F");
  tree->Branch("cov_filtered_omegaomega", cov_filtered_omegaomega, "cov_filtered_omegaomega[nsites]/F");
  tree->Branch("cov_filtered_z0phi0", cov_filtered_z0phi0, "cov_filtered_z0phi0[nsites]/F");
  tree->Branch("cov_filtered_z0omega", cov_filtered_z0omega, "cov_filtered_z0omega[nsites]/F");
  tree->Branch("cov_filtered_z0z0", cov_filtered_z0z0, "cov_filtered_z0z0[nsites]/F");
  tree->Branch("cov_filtered_tanLd0", cov_filtered_tanLd0, "cov_filtered_tanLd0[nsites]/F");
  tree->Branch("cov_filtered_tanLphi0", cov_filtered_tanLphi0, "cov_filtered_tanLphi0[nsites]/F");
  tree->Branch("cov_filtered_tanLomega", cov_filtered_tanLomega, "cov_filtered_tanLomega[nsites]/F");
  tree->Branch("cov_filtered_tanLz0", cov_filtered_tanLz0, "cov_filtered_tanLz0[nsites]/F");
  tree->Branch("cov_filtered_tanLtanL", cov_filtered_tanLtanL, "cov_filtered_tanLtanL[nsites]/F");
}

void MarlinTrkNtuple::Init(TTree* tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  std::cout << " tree address " << tree << std::endl;

  // Set branch addresses and branch pointers
  if (!tree)
    return;

  tree->Print();

  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("error_code", &error_code, &b_error_code);
  fChain->SetBranchAddress("nsites", &nsites, &b_nsites);
  fChain->SetBranchAddress("nsites_vxd", &nsites_vxd, &b_nsites_vxd);
  fChain->SetBranchAddress("nsites_sit", &nsites_sit, &b_nsites_sit);
  fChain->SetBranchAddress("nsites_ftd", &nsites_ftd, &b_nsites_ftd);
  fChain->SetBranchAddress("nsites_tpc", &nsites_tpc, &b_nsites_tpc);
  fChain->SetBranchAddress("nsites_set", &nsites_set, &b_nsites_set);
  fChain->SetBranchAddress("x_mcp", &x_mcp, &b_x_mcp);
  fChain->SetBranchAddress("y_mcp", &y_mcp, &b_y_mcp);
  fChain->SetBranchAddress("z_mcp", &z_mcp, &b_z_mcp);
  fChain->SetBranchAddress("px_mcp", &px_mcp, &b_px_mcp);
  fChain->SetBranchAddress("py_mcp", &py_mcp, &b_py_mcp);
  fChain->SetBranchAddress("pz_mcp", &pz_mcp, &b_pz_mcp);
  fChain->SetBranchAddress("p_mcp", &p_mcp, &b_p_mcp);
  fChain->SetBranchAddress("theta_mcp", &theta_mcp, &b_theta_mcp);
  fChain->SetBranchAddress("phi_mcp", &phi_mcp, &b_phi_mcp);
  fChain->SetBranchAddress("pdg_mcp", &pdg_mcp, &b_pdg_mcp);
  fChain->SetBranchAddress("d0_mcp", &d0_mcp, &b_d0_mcp);
  fChain->SetBranchAddress("phi0_mcp", &phi0_mcp, &b_phi0_mcp);
  fChain->SetBranchAddress("omega_mcp", &omega_mcp, &b_omega_mcp);
  fChain->SetBranchAddress("z0_mcp", &z0_mcp, &b_z0_mcp);
  fChain->SetBranchAddress("tanL_mcp", &tanL_mcp, &b_tanL_mcp);
  fChain->SetBranchAddress("d0_seed", &d0_seed, &b_d0_seed);
  fChain->SetBranchAddress("phi0_seed", &phi0_seed, &b_phi0_seed);
  fChain->SetBranchAddress("omega_seed", &omega_seed, &b_omega_seed);
  fChain->SetBranchAddress("z0_seed", &z0_seed, &b_z0_seed);
  fChain->SetBranchAddress("tanL_seed", &tanL_seed, &b_tanL_seed);
  fChain->SetBranchAddress("seed_ref_point_x", &seed_ref_point_x, &b_seed_ref_point_x);
  fChain->SetBranchAddress("seed_ref_point_y", &seed_ref_point_y, &b_seed_ref_point_y);
  fChain->SetBranchAddress("seed_ref_point_z", &seed_ref_point_z, &b_seed_ref_point_z);
  fChain->SetBranchAddress("cov_seed_d0d0", &cov_seed_d0d0, &b_cov_seed_d0d0);
  fChain->SetBranchAddress("cov_seed_phi0d0", &cov_seed_phi0d0, &b_cov_seed_phi0d0);
  fChain->SetBranchAddress("cov_seed_phi0phi0", &cov_seed_phi0phi0, &b_cov_seed_phi0phi0);
  fChain->SetBranchAddress("cov_seed_kappad0", &cov_seed_kappad0, &b_cov_seed_kappad0);
  fChain->SetBranchAddress("cov_seed_kappaphi0", &cov_seed_kappaphi0, &b_cov_seed_kappaphi0);
  fChain->SetBranchAddress("cov_seed_kappakappa", &cov_seed_kappakappa, &b_cov_seed_kappakappa);
  fChain->SetBranchAddress("cov_seed_z0d0", &cov_seed_z0d0, &b_cov_seed_z0d0);
  fChain->SetBranchAddress("cov_seed_z0phi0", &cov_seed_z0phi0, &b_cov_seed_z0phi0);
  fChain->SetBranchAddress("cov_seed_z0kappa", &cov_seed_z0kappa, &b_cov_seed_z0kappa);
  fChain->SetBranchAddress("cov_seed_z0z0", &cov_seed_z0z0, &b_cov_seed_z0z0);
  fChain->SetBranchAddress("cov_seed_tanLd0", &cov_seed_tanLd0, &b_cov_seed_tanLd0);
  fChain->SetBranchAddress("cov_seed_tanLphi0", &cov_seed_tanLphi0, &b_cov_seed_tanLphi0);
  fChain->SetBranchAddress("cov_seed_tanLomega", &cov_seed_tanLkappa, &b_cov_seed_tanLkappa);
  fChain->SetBranchAddress("cov_seed_tanLz0", &cov_seed_tanLz0, &b_cov_seed_tanLz0);
  fChain->SetBranchAddress("cov_seed_tanLtanL", &cov_seed_tanLtanL, &b_cov_seed_tanLtanL);
  fChain->SetBranchAddress("d0_ip", &d0_ip, &b_d0_ip);
  fChain->SetBranchAddress("phi0_ip", &phi0_ip, &b_phi0_ip);
  fChain->SetBranchAddress("omega_ip", &omega_ip, &b_omega_ip);
  fChain->SetBranchAddress("z0_ip", &z0_ip, &b_z0_ip);
  fChain->SetBranchAddress("tanL_ip", &tanL_ip, &b_tanL_ip);
  fChain->SetBranchAddress("cov_ip_d0d0", &cov_ip_d0d0, &b_cov_ip_d0d0);
  fChain->SetBranchAddress("cov_ip_phi0d0", &cov_ip_phi0d0, &b_cov_ip_phi0d0);
  fChain->SetBranchAddress("cov_ip_phi0phi0", &cov_ip_phi0phi0, &b_cov_ip_phi0phi0);
  fChain->SetBranchAddress("cov_ip_omegad0", &cov_ip_omegad0, &b_cov_ip_omegad0);
  fChain->SetBranchAddress("cov_ip_omegaphi0", &cov_ip_omegaphi0, &b_cov_ip_omegaphi0);
  fChain->SetBranchAddress("cov_ip_omegaomega", &cov_ip_omegaomega, &b_cov_ip_omegaomega);
  fChain->SetBranchAddress("cov_ip_z0d0", &cov_ip_z0d0, &b_cov_ip_z0d0);
  fChain->SetBranchAddress("cov_ip_z0phi0", &cov_ip_z0phi0, &b_cov_ip_z0phi0);
  fChain->SetBranchAddress("cov_ip_z0omega", &cov_ip_z0omega, &b_cov_ip_z0omega);
  fChain->SetBranchAddress("cov_ip_z0z0", &cov_ip_z0z0, &b_cov_ip_z0z0);
  fChain->SetBranchAddress("cov_ip_tanLd0", &cov_ip_tanLd0, &b_cov_ip_tanLd0);
  fChain->SetBranchAddress("cov_ip_tanLphi0", &cov_ip_tanLphi0, &b_cov_ip_tanLphi0);
  fChain->SetBranchAddress("cov_ip_tanLomega", &cov_ip_tanLomega, &b_cov_ip_tanLomega);
  fChain->SetBranchAddress("cov_ip_tanLz0", &cov_ip_tanLz0, &b_cov_ip_tanLz0);
  fChain->SetBranchAddress("cov_ip_tanLtanL", &cov_ip_tanLtanL, &b_cov_ip_tanLtanL);
  fChain->SetBranchAddress("ndf", &ndf, &b_ndf);
  fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
  fChain->SetBranchAddress("prob", &prob, &b_prob);
  fChain->SetBranchAddress("CellID0", CellID0, &b_CellID0);
  fChain->SetBranchAddress("rejected", rejected, &b_rejected);
  fChain->SetBranchAddress("site_x", site_x, &b_site_x);
  fChain->SetBranchAddress("site_y", site_y, &b_site_y);
  fChain->SetBranchAddress("site_z", site_z, &b_site_z);
  fChain->SetBranchAddress("ref_point_x", ref_point_x, &b_ref_point_x);
  fChain->SetBranchAddress("ref_point_y", ref_point_y, &b_ref_point_y);
  fChain->SetBranchAddress("ref_point_z", ref_point_z, &b_ref_point_z);
  fChain->SetBranchAddress("d0_mc", d0_mc, &b_d0_mc);
  fChain->SetBranchAddress("phi0_mc", phi0_mc, &b_phi0_mc);
  fChain->SetBranchAddress("omega_mc", omega_mc, &b_omega_mc);
  fChain->SetBranchAddress("z0_mc", z0_mc, &b_z0_mc);
  fChain->SetBranchAddress("tanL_mc", tanL_mc, &b_tanL_mc);
  fChain->SetBranchAddress("d0_predicted", d0_predicted, &b_d0_predicted);
  fChain->SetBranchAddress("phi0_predicted", phi0_predicted, &b_phi0_predicted);
  fChain->SetBranchAddress("omega_predicted", omega_predicted, &b_omega_predicted);
  fChain->SetBranchAddress("z0_predicted", z0_predicted, &b_z0_predicted);
  fChain->SetBranchAddress("tanL_predicted", tanL_predicted, &b_tanL_predicted);
  fChain->SetBranchAddress("d0_filtered", d0_filtered, &b_d0_filtered);
  fChain->SetBranchAddress("phi0_filtered", phi0_filtered, &b_phi0_filtered);
  fChain->SetBranchAddress("omega_filtered", omega_filtered, &b_omega_filtered);
  fChain->SetBranchAddress("z0_filtered", z0_filtered, &b_z0_filtered);
  fChain->SetBranchAddress("tanL_filtered", tanL_filtered, &b_tanL_filtered);
  fChain->SetBranchAddress("d0_smoothed", d0_smoothed, &b_d0_smoothed);
  fChain->SetBranchAddress("phi0_smoothed", phi0_smoothed, &b_phi0_smoothed);
  fChain->SetBranchAddress("omega_smoothed", omega_smoothed, &b_omega_smoothed);
  fChain->SetBranchAddress("z0_smoothed", z0_smoothed, &b_z0_smoothed);
  fChain->SetBranchAddress("tanL_smoothed", tanL_smoothed, &b_tanL_smoothed);
  fChain->SetBranchAddress("chi2_inc_filtered", chi2_inc_filtered, &b_chi2_inc_filtered);
  fChain->SetBranchAddress("chi2_inc_smoothed", chi2_inc_smoothed, &b_chi2_inc_smoothed);
  fChain->SetBranchAddress("dim", dim, &b_dim);
  fChain->SetBranchAddress("cov_smoothed_d0d0", cov_smoothed_d0d0, &b_cov_smoothed_d0d0);
  fChain->SetBranchAddress("cov_smoothed_phi0d0", cov_smoothed_phi0d0, &b_cov_smoothed_phi0d0);
  fChain->SetBranchAddress("cov_smoothed_phi0phi0", cov_smoothed_phi0phi0, &b_cov_smoothed_phi0phi0);
  fChain->SetBranchAddress("cov_smoothed_omegad0", cov_smoothed_omegad0, &b_cov_smoothed_omegad0);
  fChain->SetBranchAddress("cov_smoothed_omegaphi0", cov_smoothed_omegaphi0, &b_cov_smoothed_omegaphi0);
  fChain->SetBranchAddress("cov_smoothed_omegaomega", cov_smoothed_omegaomega, &b_cov_smoothed_omegaomega);
  fChain->SetBranchAddress("cov_smoothed_z0d0", cov_smoothed_z0d0, &b_cov_smoothed_z0d0);
  fChain->SetBranchAddress("cov_smoothed_z0phi0", cov_smoothed_z0phi0, &b_cov_smoothed_z0phi0);
  fChain->SetBranchAddress("cov_smoothed_z0omega", cov_smoothed_z0omega, &b_cov_smoothed_z0omega);
  fChain->SetBranchAddress("cov_smoothed_z0z0", cov_smoothed_z0z0, &b_cov_smoothed_z0z0);
  fChain->SetBranchAddress("cov_smoothed_tanLd0", cov_smoothed_tanLd0, &b_cov_smoothed_tanLd0);
  fChain->SetBranchAddress("cov_smoothed_tanLphi0", cov_smoothed_tanLphi0, &b_cov_smoothed_tanLphi0);
  fChain->SetBranchAddress("cov_smoothed_tanLomega", cov_smoothed_tanLomega, &b_cov_smoothed_tanLomega);
  fChain->SetBranchAddress("cov_smoothed_tanLz0", cov_smoothed_tanLz0, &b_cov_smoothed_tanLz0);
  fChain->SetBranchAddress("cov_smoothed_tanLtanL", cov_smoothed_tanLtanL, &b_cov_smoothed_tanLtanL);
  fChain->SetBranchAddress("cov_predicted_d0d0", cov_predicted_d0d0, &b_cov_predicted_d0d0);
  fChain->SetBranchAddress("cov_predicted_phi0d0", cov_predicted_phi0d0, &b_cov_predicted_phi0d0);
  fChain->SetBranchAddress("cov_predicted_phi0phi0", cov_predicted_phi0phi0, &b_cov_predicted_phi0phi0);
  fChain->SetBranchAddress("cov_predicted_omegad0", cov_predicted_omegad0, &b_cov_predicted_omegad0);
  fChain->SetBranchAddress("cov_predicted_omegaphi0", cov_predicted_omegaphi0, &b_cov_predicted_omegaphi0);
  fChain->SetBranchAddress("cov_predicted_omegaomega", cov_predicted_omegaomega, &b_cov_predicted_omegaomega);
  fChain->SetBranchAddress("cov_predicted_z0d0", cov_predicted_z0d0, &b_cov_predicted_z0d0);
  fChain->SetBranchAddress("cov_predicted_z0phi0", cov_predicted_z0phi0, &b_cov_predicted_z0phi0);
  fChain->SetBranchAddress("cov_predicted_z0omega", cov_predicted_z0omega, &b_cov_predicted_z0omega);
  fChain->SetBranchAddress("cov_predicted_z0z0", cov_predicted_z0z0, &b_cov_predicted_z0z0);
  fChain->SetBranchAddress("cov_predicted_tanLd0", cov_predicted_tanLd0, &b_cov_predicted_tanLd0);
  fChain->SetBranchAddress("cov_predicted_tanLphi0", cov_predicted_tanLphi0, &b_cov_predicted_tanLphi0);
  fChain->SetBranchAddress("cov_predicted_tanLomega", cov_predicted_tanLomega, &b_cov_predicted_tanLomega);
  fChain->SetBranchAddress("cov_predicted_tanLz0", cov_predicted_tanLz0, &b_cov_predicted_tanLz0);
  fChain->SetBranchAddress("cov_predicted_tanLtanL", cov_predicted_tanLtanL, &b_cov_predicted_tanLtanL);
  fChain->SetBranchAddress("cov_filtered_d0d0", cov_filtered_d0d0, &b_cov_filtered_d0d0);
  fChain->SetBranchAddress("cov_filtered_phi0d0", cov_filtered_phi0d0, &b_cov_filtered_phi0d0);
  fChain->SetBranchAddress("cov_filtered_phi0phi0", cov_filtered_phi0phi0, &b_cov_filtered_phi0phi0);
  fChain->SetBranchAddress("cov_filtered_omegad0", cov_filtered_omegad0, &b_cov_filtered_omegad0);
  fChain->SetBranchAddress("cov_filtered_omegaphi0", cov_filtered_omegaphi0, &b_cov_filtered_omegaphi0);
  fChain->SetBranchAddress("cov_filtered_omegaomega", cov_filtered_omegaomega, &b_cov_filtered_omegaomega);
  fChain->SetBranchAddress("cov_filtered_z0d0", cov_filtered_z0d0, &b_cov_filtered_z0d0);
  fChain->SetBranchAddress("cov_filtered_z0phi0", cov_filtered_z0phi0, &b_cov_filtered_z0phi0);
  fChain->SetBranchAddress("cov_filtered_z0omega", cov_filtered_z0omega, &b_cov_filtered_z0omega);
  fChain->SetBranchAddress("cov_filtered_z0z0", cov_filtered_z0z0, &b_cov_filtered_z0z0);
  fChain->SetBranchAddress("cov_filtered_tanLd0", cov_filtered_tanLd0, &b_cov_filtered_tanLd0);
  fChain->SetBranchAddress("cov_filtered_tanLphi0", cov_filtered_tanLphi0, &b_cov_filtered_tanLphi0);
  fChain->SetBranchAddress("cov_filtered_tanLomega", cov_filtered_tanLomega, &b_cov_filtered_tanLomega);
  fChain->SetBranchAddress("cov_filtered_tanLz0", cov_filtered_tanLz0, &b_cov_filtered_tanLz0);
  fChain->SetBranchAddress("cov_filtered_tanLtanL", cov_filtered_tanLtanL, &b_cov_filtered_tanLtanL);

  Notify();
}

Bool_t MarlinTrkNtuple::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void MarlinTrkNtuple::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain)
    return;
  fChain->Show(entry);
}
Int_t MarlinTrkNtuple::Cut(Long64_t entry) {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef MarlinTrkNtuple_cxx
