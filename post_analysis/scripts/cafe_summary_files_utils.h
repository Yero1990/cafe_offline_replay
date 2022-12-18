/* 
   ------------------
   Author: C. Yero
   
   Brief: A compilation of CaFe utility functions
   for reading summary files and performing numerical 
   operations on them
   
   ------------------
*/ 

#include "../../UTILS_CAFE/UTILS/parse_utils.h"
#include "../../UTILS_CAFE/UTILS/read_csv.h"
#include "../../UTILS_CAFE/UTILS/vector_operations.h"


double get(string var="", string target="", string kin=""){

  /* brief: returns pre-determined variable for particular (target, kin)
     
     ----------
     arguments:
     ----------
     
     get the integrated (either total or average of the following quantities_)
     var: total_charge, real_yield, real_yield_err, hms_trk_eff, hms_trk_eff_err, shms_trk_eff, shms_trk_eff_err


   */
 
  
  // set generic .csv file name
  string file_csv = Form("summary_files_pass1/EmissCut_100MeV/cafe_prod_%s_%s_report_summary.csv", target.c_str(), kin.c_str());


 
  // return total charge (sums over charge for each run)
  if( var.compare("total_charge")==0 ){
    vector<double> v_charge         = read_csv(file_csv.c_str(), "charge"); // mC
    double charge = vsum(v_charge);
    
    return charge;
  }

  // return total yield (sums over real yield for each run)
  if( var.compare("real_yield")==0 ){
    
    vector<double> v_real_Yield     = read_csv(file_csv, "real_Yield");
    double real_Yield = vsum(v_real_Yield);
    
    return real_Yield;
  }

  // return total yield err ( root of the sum of errors ^{2} ) - basic error propagation for sum of variables  
  if( var.compare("real_yield_err")==0 ){
    
    vector<double> v_real_Yield_err = read_csv(file_csv, "real_Yield_err");
    double real_Yield_err = sqrt( vsum ( vpow(v_real_Yield_err, 2) ) ); // sqrt [  err_1^{2} + err_2^{2} + . .  . err_n^{2} ]
    
    return real_Yield_err;
    
  }

  // return weighted avg HMS tracking (and error) efficiency (presumably tracking efficiency is ~ same for every run of the same kinematic)
  if( var.compare("hms_trk_eff")==0 || var.compare("hms_trk_eff_err")==0 ){

    vector<double> v_hTrkEff        = read_csv(file_csv, "hTrkEff");
    vector<double> v_hTrkEff_err    = read_csv(file_csv, "hTrkEff_err");

    //error in weighted average is passed by reference
    double hms_trk_eff_err = 0;          

    // calculated weighted average of a vector of elements
    double hms_trk_eff              = vavgw(v_hTrkEff, v_hTrkEff_err, hms_trk_eff_err);  

    
    if( var.compare("hms_trk_eff")==0 ){
      return hms_trk_eff;
    }

    if( var.compare("hms_trk_eff_err")==0 ){
      return hms_trk_eff_err;
    }    

  }

  // return weighted avg SHMS tracking (and error) efficiency (presumably tracking efficiency is ~ same for every run of the same kinematic)
  if( var.compare("shms_trk_eff")==0 || var.compare("shms_trk_eff_err")==0 ){

    vector<double> v_pTrkEff        = read_csv(file_csv, "pTrkEff");
    vector<double> v_pTrkEff_err    = read_csv(file_csv, "pTrkEff_err");

    //error in weighted average is passed by reference
    double shms_trk_eff_err = 0;          

    // calculated weighted average of a vector of elements
    double shms_trk_eff              = vavgw(v_pTrkEff, v_pTrkEff_err, shms_trk_eff_err);  

    
    if( var.compare("shms_trk_eff")==0 ){
      return shms_trk_eff;
    }

    if( var.compare("shms_trk_eff_err")==0 ){
      return shms_trk_eff_err;
    }    

  }


  // return weighted avg total live time (presumably live time is ~ same for every run of the same kinematic)
  if( var.compare("total_live_time")==0 || var.compare("total_live_time_err")==0 ){

   
    vector<double> v_tLT            = read_csv(file_csv, "tLT");
    vector<double> v_tLT_err        = read_csv(file_csv, "tLT_err_Bi");

    //error in weighted average is passed by reference
    double total_live_time_err = 0;          

    // calculated weighted average of a vector of elements
    double total_live_time              = vavgw(v_tLT, v_tLT_err, total_live_time_err);  

    
    if( var.compare("total_live_time")==0 ){
      return total_live_time;
    }

    if( var.compare("total_live_time_err")==0 ){
      return total_live_time_err;
    }    

  }

  
  return 0;
  
}


double get_param( string var="", string target="", string kin="" ){
  
  // brief: function to read parameters from the .csv files


  // set generic .csv file name
  string file_csv = Form("summary_files_pass1/EmissCut_100MeV/cafe_prod_%s_%s_report_summary.csv", target.c_str(), kin.c_str());

  // find parameters (commented variables)

  // NOTE:  need to parse string and conver to double
  //double transparency = FindString("transparency:", file_csv.c_str(), false, -1, true);
  //double tgt_area_density = FindString("target_areal_density", file_csv.c_str(), false, -1, true);

  //double Z = FindString("Z:", file_csv.c_str(), false, -1, true);
  //double N = FindString("N:", file_csv.c_str(), false, -1, true);
  //double A = FindString("A:", file_csv.c_str(), false, -1, true);


  
}

