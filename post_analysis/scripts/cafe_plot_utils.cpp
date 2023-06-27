#include "cafe_plot_utils.h"

void cafe_plot_utils(){
  
  cout << "" << endl;
  cout << "" << endl;
  cout << "--------------------------------" << endl;
  cout << "" << endl;
  cout << "Welcome to CaFe Plot Utils" << endl;
  cout << "" << endl;
  cout << "Please call a plotting utility \nfunction from this script " << endl;
  cout << "--------------------------------" << endl;

  // data/simc yield comparisons (yields scaled by charge, efficiencies, etc)
  
  //user may define scale factor within function)
  //cout <<  "hist_id: \n 0: Pm, 1: W, 2: Q2, 3: xbj, 4: nu, 5: Em, 6: th_rq, 7: th_pq \n
  //                      8: th_e, 9: kf, 10: th_p, 11: Pf, 12: |q|, 13: th_q, 14: Pmx_lab, \n
  //                      15: Pmy_lab, 16: Pmz_lab " << endl;

  TString target = "C12";
  TString kin = "MF";
  int run= 17098;

  
  yield_comparison(target.Data(), kin.Data(), run, 0);
  
  yield_comparison(target.Data(), kin.Data(), run, 1);
  
  yield_comparison(target.Data(), kin.Data(), run, 2);
  yield_comparison(target.Data(), kin.Data(), run, 3);
  yield_comparison(target.Data(), kin.Data(), run, 4);
  yield_comparison(target.Data(), kin.Data(), run, 5);
  yield_comparison(target.Data(), kin.Data(), run, 6);
  yield_comparison(target.Data(), kin.Data(), run, 7);
 
  
  yield_comparison(target.Data(), kin.Data(), run, 8);
  yield_comparison(target.Data(), kin.Data(), run, 9);
  yield_comparison(target.Data(), kin.Data(), run, 10);
  yield_comparison(target.Data(), kin.Data(), run, 11);
  yield_comparison(target.Data(), kin.Data(), run, 12);
  yield_comparison(target.Data(), kin.Data(), run, 13);
  
  yield_comparison(target.Data(), kin.Data(), run, 14);
  yield_comparison(target.Data(), kin.Data(), run, 15);
  yield_comparison(target.Data(), kin.Data(), run, 16);
  
    
}

  //=====================
  // EXAMPLES FOR USERS
  //=====================


  //--------------------------------------------
  // READ PARAMETERS FROM SUMMARY (.csv) FILES
  //--------------------------------------------

  
  // get summary file columns combined (either total or average or per run)
  /*
  double total_charge    = get_header("total_charge", "B11", "SRC"   );
  double total_yield     = get_header("real_yield", "Be9", "SRC"     );
   double total_yield_err = get_header("real_yield_err", "B10", "MF" );

  double hms_trk_eff     = get_header("hms_trk_eff", "B11", "SRC" );
  double hms_trk_err_err = get_header("hms_trk_eff_err", "LD2", "SRC" );

  double shms_trk_eff     = get_header("shms_trk_eff", "C12", "SRC" );
  double shms_trk_err_err = get_header("shms_trk_eff_err", "B11", "SRC" );
  
  double_t live_time      = get_header("total_live_time", "Fe54",  "MF");

  double transparency     = get_param("transparency","C12", "MF" );
  double tgt_area_density = get_param("tgt_area_density","Fe54", "SRC" );
  
  double N    = get_param("N","C12", "MF" );  # get number of neutrons
  double Z    = get_param("Z", "Fe54", "MF" );  # get number of protons
  double A    = get_param("A","Be9", "MF" );  # get mass number (N+Z)

  
  */


  

  // brief: this function returns an overlay of any two histogram objects from any two files (or histograms) input by the user
  //        the user also provides labels, titles and legend text (this is basically a quick way to make quality plot comparisons, without the
  //        hassle of dealing with ROOT ), there is also an optional flag ( bool norm ) to draw normalized histograms to an areal of 1.
  
  // version 0: returns overlay of histogram objects from a single input file

  /*
  //Example: overlyaing same histograms with different cuts (sequential cuts study)
  vector<TString> hist_name_Q2_mf ={"quality_plots/ACCP+PID+CTIME_CUTS/H_Q2_ACCP_PID_CTIME_CUTS",
				    "quality_plots/ACCP+PID+CTIME+Q2_CUTS/H_Q2_ACCP_PID_CTIME_Q2_CUTS",
				    "quality_plots/ACCP+PID+CTIME+Q2+Em_CUTS/H_Q2_ACCP_PID_CTIME_Q2_Em_CUTS",
				    "quality_plots/ACCP+PID+CTIME+Q2+Em+Pm_CUTS/H_Q2_ACCP_PID_CTIME_Q2_Em_Pm_CUTS",
  };

  vector<TString> hist_name_Q2_src ={"quality_plots/ACCP+PID+CTIME_CUTS/H_Q2_ACCP_PID_CTIME_CUTS",
				     "quality_plots/ACCP+PID+CTIME+Q2_CUTS/H_Q2_ACCP_PID_CTIME_Q2_CUTS",
				     "quality_plots/ACCP+PID+CTIME+Q2+Xbj_CUTS/H_Q2_ACCP_PID_CTIME_Q2_Xbj_CUTS",
				     "quality_plots/ACCP+PID+CTIME+Q2+Xbj+thrq_CUTS/H_Q2_ACCP_PID_CTIME_Q2_Xbj_thrq_CUTS",
				     "quality_plots/ACCP+PID+CTIME+Q2+Xbj+thrq+Pm_CUTS/H_Q2_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm_CUTS",
				     
  };


  vector<TString> hist_name_xbj_mf ={"quality_plots/ACCP+PID+CTIME_CUTS/H_xbj_ACCP_PID_CTIME_CUTS",
				    "quality_plots/ACCP+PID+CTIME+Q2_CUTS/H_xbj_ACCP_PID_CTIME_Q2_CUTS",
				    "quality_plots/ACCP+PID+CTIME+Q2+Em_CUTS/H_xbj_ACCP_PID_CTIME_Q2_Em_CUTS",
				    "quality_plots/ACCP+PID+CTIME+Q2+Em+Pm_CUTS/H_xbj_ACCP_PID_CTIME_Q2_Em_Pm_CUTS",
  };

  vector<TString> hist_name_xbj_src ={"quality_plots/ACCP+PID+CTIME_CUTS/H_xbj_ACCP_PID_CTIME_CUTS",
				     "quality_plots/ACCP+PID+CTIME+Q2_CUTS/H_xbj_ACCP_PID_CTIME_Q2_CUTS",
				     "quality_plots/ACCP+PID+CTIME+Q2+Xbj_CUTS/H_xbj_ACCP_PID_CTIME_Q2_Xbj_CUTS",
				     "quality_plots/ACCP+PID+CTIME+Q2+Xbj+thrq_CUTS/H_xbj_ACCP_PID_CTIME_Q2_Xbj_thrq_CUTS",
				     "quality_plots/ACCP+PID+CTIME+Q2+Xbj+thrq+Pm_CUTS/H_xbj_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm_CUTS",
				     
  };
   
  
  vector<int> clr_mf           ={2, 4, 6, 8};
  vector<int> clr_src          ={2, 4, 6, 8, 9};

  vector<TString> hist_leg_MF  ={"accp+pid+ctime", "Q^{2}", "Q^{2}+Em", "Q^{2}+Em+Pm"};
  vector<TString> hist_leg_SRC  ={"accp+pid+ctime", "Q^{2}", "Q^{2}+X_{bj}", "Q^{2}+X_{bj}+#theta_{rq}",  "Q^{2}+X_{bj}+#theta_{rq}+P_{m}"};
  */

  
  // --- plot Q2 MF for all targets ----
  //compare_histos("analyzed_files/combined/cafe_prod_LD2_MF_combined.root", hist_name_Q2_mf, clr_mf, hist_leg_MF, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: LD2 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Be9_MF_combined.root", hist_name_Q2_mf, clr_mf, hist_leg_MF, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: Be9 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_B10_MF_combined.root", hist_name_Q2_mf, clr_mf, hist_leg_MF, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: B10 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_B11_MF_combined.root", hist_name_Q2_mf, clr_mf, hist_leg_MF, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: B11 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_C12_MF_combined.root", hist_name_Q2_mf, clr_mf, hist_leg_MF, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: C12 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Ca40_MF_combined.root", hist_name_Q2_mf, clr_mf, hist_leg_MF, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: Ca40 MF", false);
  //compare_histos("analyzed_files/individual/pass1/cafe_prod_Ca48_MF_17096_-1_histos.root", hist_name_Q2_mf, clr_mf, hist_leg_MF, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: Ca48 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Fe54_MF_combined.root", hist_name_Q2_mf, clr_mf, hist_leg_MF, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: Fe54 MF", false);

  // ---- plot Q2 SRC for all targets -----
  //compare_histos("analyzed_files/combined/cafe_prod_LD2_SRC_combined.root", hist_name_Q2_src, clr_src, hist_leg_SRC, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: LD2 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Be9_SRC_combined.root", hist_name_Q2_src, clr_src, hist_leg_SRC, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: Be9 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_B10_SRC_combined.root", hist_name_Q2_src, clr_src, hist_leg_SRC, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: B10 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_B11_SRC_combined.root", hist_name_Q2_src, clr_src, hist_leg_SRC, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: B11 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_C12_SRC_combined.root", hist_name_Q2_src, clr_src, hist_leg_SRC, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: C12 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Ca40_SRC_combined.root", hist_name_Q2_src, clr_src, hist_leg_SRC, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: Ca40 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Ca48_SRC_combined.root", hist_name_Q2_src, clr_src, hist_leg_SRC, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: Ca48 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Fe54_SRC_combined.root", hist_name_Q2_src, clr_src, hist_leg_SRC, "Q2 [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}: Fe54 SRC", false);


   // --- plot Xbj MF for all targets ----
  //compare_histos("analyzed_files/combined/cafe_prod_LD2_MF_combined.root", hist_name_xbj_mf, clr_mf, hist_leg_MF, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: LD2 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Be9_MF_combined.root", hist_name_xbj_mf, clr_mf, hist_leg_MF, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: Be9 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_B10_MF_combined.root", hist_name_xbj_mf, clr_mf, hist_leg_MF, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: B10 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_B11_MF_combined.root", hist_name_xbj_mf, clr_mf, hist_leg_MF, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: B11 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_C12_MF_combined.root", hist_name_xbj_mf, clr_mf, hist_leg_MF, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: C12 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Ca40_MF_combined.root", hist_name_xbj_mf, clr_mf, hist_leg_MF, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: Ca40 MF", false);
  //compare_histos("analyzed_files/individual/pass1/cafe_prod_Ca48_MF_17096_-1_histos.root", hist_name_xbj_mf, clr_mf, hist_leg_MF, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: Ca48 MF", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Fe54_MF_combined.root", hist_name_xbj_mf, clr_mf, hist_leg_MF, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: Fe54 MF", false);


  // ---- plot Xbj SRC for all targets -----
  //compare_histos("analyzed_files/combined/cafe_prod_LD2_SRC_combined.root", hist_name_xbj_src, clr_src, hist_leg_SRC, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: LD2 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Be9_SRC_combined.root", hist_name_xbj_src, clr_src, hist_leg_SRC, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: Be9 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_B10_SRC_combined.root", hist_name_xbj_src, clr_src, hist_leg_SRC, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: B10 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_B11_SRC_combined.root", hist_name_xbj_src, clr_src, hist_leg_SRC, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: B11 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_C12_SRC_combined.root", hist_name_xbj_src, clr_src, hist_leg_SRC, "X_{bj}", "Counts", "x-Bjorken, X_{bj}: C12 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Ca40_SRC_combined.root", hist_name_xbj_src, clr_src, hist_leg_SRC, "X_{bj}", "Counts","x-Bjorken, X_{bj}: Ca40 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Ca48_SRC_combined.root", hist_name_xbj_src, clr_src, hist_leg_SRC, "X_{bj}", "Counts","x-Bjorken, X_{bj}: Ca48 SRC", false);
  //compare_histos("analyzed_files/combined/cafe_prod_Fe54_SRC_combined.root", hist_name_xbj_src, clr_src, hist_leg_SRC, "X_{bj}", "Counts","x-Bjorken, X_{bj}: Fe54 SRC", false);
  
  /*
  // version 1: input file path and histogram objects by user
  //compare_histos("path/to/file1.root", "kin_plots/H_Pm", "path/to/file2.root", "kin_plots/H_Pm", "X-label [units]", "Y-label [units]", "title",
  //		 "hist1_legend_title", "hist2_legend_title", true);

  // version 2: input histogram objects by user (assumed that two histograms have the same binning)
  // example: use dummy histograms to test
  TH1F *h1_test = new TH1F("h1_test", "hist 1", 100,-5,5);
  TH1F *h2_test = new TH1F("h2_test", "hist 2", 100,-5,5);
  h1_test->FillRandom("gaus"); h2_test->FillRandom("landau");
  
  compare_histos(h1_test, h2_test, "X-label [units]", "Y-label [units]", "title",
  		 "hist1_legend_title", "hist2_legend_title", true);

  */
  

  

  // brief: this function returns an overlay of the targets input by the user, with corresponding colors for any given histogram
  //        The histograms are read from analyzed .root files in a fixed location specified within the function.
  
  // NOTE: All .root files and summary (.csv) files used by these functions are assumed to be in a fixed location (set within the function)

  // histo names: randSub_plots/H_Em_nuc_rand_sub, randSub_plots/H_Q2_rand_sub, randSub_plots/H_xbj_rand_sub, randSub_plots/H_Pm_rand_sub,
  
  //vector<string> tgt = {"LD2", "Be9", "B10", "B11", "C12"};
  //vector<int> clr     = {2, 4, 6, 8, 9};    //  root cern color code ---> 2 : red,  4: blue, 6: magenta, 8: green, 9: purple

  //vector<string> tgt = {"Be9", "B10", "B11", "C12"};
  //vector<int> clr     = {1, 2, 4, 8};    //  root cern color code ---> 1: black,  2 : red,  4: blue, 6: magenta, 8: green, 9: purple

  //overlay_nuclei(tgt, clr, "SRC", "randSub_plots/H_Q2_rand_sub",     "Q^{2} [GeV^{2}]", "Normalized Counts",  "4-Momentum Transfer (light nuclei)");
  // overlay_nuclei(tgt, clr, "SRC", "randSub_plots/H_xbj_rand_sub",    "x_{Bj}",          "Normalized Counts",  "x-Bjorken (light nuclei)");
  //overlay_nuclei(tgt, clr, "SRC", "randSub_plots/H_Em_nuc_rand_sub", "E_{m} [GeV]",     "Normalized Counts",  "Missing Energy (light nuclei)");
  //overlay_nuclei(tgt, clr, "SRC", "randSub_plots/H_Pm_rand_sub",     "P_{m} [GeV/c]",   "Normalized Counts",  "Missing Momentum (light nuclei)");
 
  

  
  
  //-------------------------
  // HISTOGRAM SINGLE RATIOS
  //-------------------------
  // brief: this function returns any two histograms (A,B) and ratio, A/B. 
  //         The histograms (A,B) are read from analyzed combined .root files in a fixed location specified within the function.
  //         and A/B is calculated within the ratio function

  // NOTE: All .root files and summary (.csv) files used by these functions are assumed to be in a fixed location (set within the function)
    
  // Example of Use: 

  /*
  //declare histogram vector to retrieve histograms A: Ca48 MF,  B: Ca40 MF, and R: A/B, binned in Pmiss
  vector<TH1F*> hvec_be9;
  vector<TH1F*> hvec_b10;
  vector<TH1F*> hvec_b11;
  

  // call function to get ratio
  // arguments:("targetA", "kinematicsA", "targetB", "kinematicsB", "histogram object", bool show_histos? )
  hvec_be9 = get_single_ratios("Be9",       "MF",        "C12",       "MF",       "randSub_plots/H_Pm_rand_sub",     true);  //hvec[0]: histA, hvec[1]: histB, hvec[3]: R
  hvec_b10 = get_single_ratios("B10",       "MF",        "C12",       "MF",       "randSub_plots/H_Pm_rand_sub",     true); 
  hvec_b11 = get_single_ratios("B11",       "MF",        "C12",       "MF",       "randSub_plots/H_Pm_rand_sub",     true); 

  // access the histograms from the vector
  TH1F *be9 = &*hvec_be9[0];  // get histA
  TH1F *c12 = &*hvec_be9[1];  // get histB
  TH1F *b10 = &*hvec_b10[0];  
  TH1F *b11 = &*hvec_b11[0];

  be9->SetLineColor(kBlack);
  be9->SetLineWidth(3);
  b10->SetLineColor(kRed);
  b10->SetLineWidth(3);
  b11->SetLineColor(kBlue);
  b11->SetLineWidth(3);
  c12->SetLineColor(kGreen);
  c12->SetLineWidth(3);
  
  //access ratios
  TH1F *be9r = &*hvec_be9[2];
  TH1F *b10r = &*hvec_b10[2];
  TH1F *b11r = &*hvec_b11[2];
  be9r->SetLineColor(kBlack);
  be9r->SetMarkerColor(kBlack);
  b10r->SetLineColor(kRed);
  b10r->SetMarkerColor(kRed);
  b11r->SetLineColor(kBlue);
  b11r->SetMarkerColor(kBlue);
  
  
  // check that the retrieved histograms are as expected
  TCanvas *c = new TCanvas("c1","",800,1200);
  c->Divide(1,2);
  c->cd(1);
  be9->Draw("histE0");
  b10->Draw("sameshistE0");
  b11->Draw("sameshistE0");
  c12->Draw("sameshistE0");


  c->cd(2);
  be9r->Draw("E0");
  b10r->Draw("E0sames");
  b11r->Draw("E0sames");
  */
  
 
  
//}



/*
  For data:
  We know that nucleus A must have transparency T_A and H has transparency T_H=1
  
  R_data =  Y_A / Y_H = N_A * T_A / (N_H * T_H) = 1
  --> T_A = N_H * T_H / N_A = N_H / N_A

  T = Y_DATA / Y_PWIA

 For
  
 */
