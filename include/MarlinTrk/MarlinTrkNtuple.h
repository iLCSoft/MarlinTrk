#ifndef INCLUDE_MarlinTrkNtuple
#define INCLUDE_MarlinTrkNtuple 1

#define MAX_SITES 300
struct {
  
  int error_code;
  
  int nsites;
  
  int nsites_vxd;
  int nsites_sit;
  int nsites_ftd;
  int nsites_tpc;
  int nsites_set;
  int nsites_etd;

  float x_mcp;
  float y_mcp;
  float z_mcp;

  float px_mcp;
  float py_mcp;
  float pz_mcp;

  float p_mcp;

  float theta_mcp;
  float phi_mcp;
  int   pdg_mcp;
  
  float d0_mcp;
  float phi0_mcp;
  float omega_mcp;
  float z0_mcp;
  float tanL_mcp;

  int   ndf;
  float chi2;
  float prob;
  
  float d0_seed;
  float phi0_seed;
  float omega_seed;
  float z0_seed;
  float tanL_seed; 

  float seed_ref_point_x;
  float seed_ref_point_y;
  float seed_ref_point_z;
  
  float cov_seed_d0d0 ;      
  float cov_seed_phi0d0 ;      
  float cov_seed_phi0phi0 ;      
  float cov_seed_kappad0 ;      
  float cov_seed_kappaphi0 ;      
  float cov_seed_kappakappa ;      
  float cov_seed_z0d0 ;      
  float cov_seed_z0phi0 ;      
  float cov_seed_z0kappa ;      
  float cov_seed_z0z0 ;      
  float cov_seed_tanLz0 ;      
  float cov_seed_tanLd0 ;      
  float cov_seed_tanLphi0 ;      
  float cov_seed_tanLkappa ;      
  float cov_seed_tanLtanL ;      
  
  float d0_ip;
  float phi0_ip;
  float omega_ip;
  float z0_ip;
  float tanL_ip; 

  float cov_ip_d0d0 ;      
  float cov_ip_phi0d0 ;      
  float cov_ip_phi0phi0 ;      
  float cov_ip_omegad0 ;      
  float cov_ip_omegaphi0 ;      
  float cov_ip_omegaomega ;      
  float cov_ip_z0d0 ;      
  float cov_ip_z0phi0 ;      
  float cov_ip_z0omega ;      
  float cov_ip_z0z0 ;      
  float cov_ip_tanLz0 ;      
  float cov_ip_tanLd0 ;      
  float cov_ip_tanLphi0 ;      
  float cov_ip_tanLomega ;      
  float cov_ip_tanLtanL ;      
  
  int   CellID0[MAX_SITES];
  int   rejected[MAX_SITES];

  float site_x[MAX_SITES];
  float site_y[MAX_SITES];
  float site_z[MAX_SITES];
  
  float ref_point_x[MAX_SITES];
  float ref_point_y[MAX_SITES];
  float ref_point_z[MAX_SITES];

  float d0_mc[MAX_SITES];
  float phi0_mc[MAX_SITES];
  float omega_mc[MAX_SITES];
  float z0_mc[MAX_SITES];
  float tanL_mc[MAX_SITES];
  
  float d0_predicted[MAX_SITES];
  float phi0_predicted[MAX_SITES];
  float omega_predicted[MAX_SITES];
  float z0_predicted[MAX_SITES];
  float tanL_predicted[MAX_SITES]; 

  float d0_filtered[MAX_SITES];
  float phi0_filtered[MAX_SITES];
  float omega_filtered[MAX_SITES];
  float z0_filtered[MAX_SITES];
  float tanL_filtered[MAX_SITES]; 

  float d0_smoothed[MAX_SITES];
  float phi0_smoothed[MAX_SITES];
  float omega_smoothed[MAX_SITES];
  float z0_smoothed[MAX_SITES];
  float tanL_smoothed[MAX_SITES]; 
  

  float chi2_inc_filtered[MAX_SITES];
  float chi2_inc_smoothed[MAX_SITES];
  int   dim[MAX_SITES];

  float cov_predicted_d0d0[MAX_SITES];      
  float cov_predicted_phi0d0[MAX_SITES];      
  float cov_predicted_phi0phi0[MAX_SITES];      
  float cov_predicted_omegad0[MAX_SITES];      
  float cov_predicted_omegaphi0[MAX_SITES];      
  float cov_predicted_omegaomega[MAX_SITES];      
  float cov_predicted_z0d0[MAX_SITES];      
  float cov_predicted_z0phi0[MAX_SITES];      
  float cov_predicted_z0omega[MAX_SITES];      
  float cov_predicted_z0z0[MAX_SITES];      
  float cov_predicted_tanLz0[MAX_SITES];      
  float cov_predicted_tanLd0[MAX_SITES];      
  float cov_predicted_tanLphi0[MAX_SITES];      
  float cov_predicted_tanLomega[MAX_SITES];      
  float cov_predicted_tanLtanL[MAX_SITES];      

  float cov_filtered_d0d0[MAX_SITES];      
  float cov_filtered_phi0d0[MAX_SITES];      
  float cov_filtered_phi0phi0[MAX_SITES];      
  float cov_filtered_omegad0[MAX_SITES];      
  float cov_filtered_omegaphi0[MAX_SITES];      
  float cov_filtered_omegaomega[MAX_SITES];      
  float cov_filtered_z0d0[MAX_SITES];      
  float cov_filtered_z0phi0[MAX_SITES];      
  float cov_filtered_z0omega[MAX_SITES];      
  float cov_filtered_z0z0[MAX_SITES];      
  float cov_filtered_tanLz0[MAX_SITES];      
  float cov_filtered_tanLd0[MAX_SITES];      
  float cov_filtered_tanLphi0[MAX_SITES];      
  float cov_filtered_tanLomega[MAX_SITES];      
  float cov_filtered_tanLtanL[MAX_SITES];      

  float cov_smoothed_d0d0[MAX_SITES];      
  float cov_smoothed_phi0d0[MAX_SITES];      
  float cov_smoothed_phi0phi0[MAX_SITES];      
  float cov_smoothed_omegad0[MAX_SITES];      
  float cov_smoothed_omegaphi0[MAX_SITES];      
  float cov_smoothed_omegaomega[MAX_SITES];      
  float cov_smoothed_z0d0[MAX_SITES];      
  float cov_smoothed_z0phi0[MAX_SITES];      
  float cov_smoothed_z0omega[MAX_SITES];      
  float cov_smoothed_z0z0[MAX_SITES];      
  float cov_smoothed_tanLz0[MAX_SITES];      
  float cov_smoothed_tanLd0[MAX_SITES];      
  float cov_smoothed_tanLphi0[MAX_SITES];      
  float cov_smoothed_tanLomega[MAX_SITES];      
  float cov_smoothed_tanLtanL[MAX_SITES];      


} _track_record;



#endif
