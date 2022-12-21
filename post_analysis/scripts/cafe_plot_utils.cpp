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


  //=====================
  // EXAMPLES FOR USERS
  //=====================


  //--------------------------------------------
  // READ PARAMETERS FROM SUMMARY (.csv) FILES
  //--------------------------------------------

  /*
  // get summary file columns combined (either total or average)
  
  double total_charge    = get_header("total_charge", "B11", "SRC"   );
  double total_yield     = get_header("real_yield", "Be9", "SRC"     );
  double total_yield_err = get_header("real_yield_err", "B10", "MF" );

  double hms_trk_err     = get_header("hms_trk_eff", "B11", "SRC" );
  double hms_trk_err_err = get_header("hms_trk_eff_err", "LD2", "SRC" );

  double shms_trk_err     = get_header("shms_trk_eff", "C12", "SRC" );
  double shms_trk_err_err = get_header("shms_trk_eff_err", "B11", "SRC" );
  
  double_t live_time      = get_header("total_live_time", "Fe54,  "MF");

  double transparency     = get_param("transparency","C12", "MF" );
  double tgt_area_density = get_param("transparency","Fe54", "SRC" );
  
  double N    = get_param("N","C12", "MF" );  # get number of neutrons
  double Z    = get_param("Z, "Fe54", "MF" );  # get number of protons
  double A    = get_param("A","Be9", "MF" );  # get mass number (N+Z)

*/
  
  
  
  //-------------------------------
  // COMPARE (OVERLAY) HISTOGRAMS
  //-------------------------------
  // brief: this function returns an overlay of any two histogram objects from any two files (or histograms) input by the user
  //        the user also provides labels, titles and legend text (this is basically a quick way to make quality plot comparisons, without the
  //        hassle of dealing with ROOT ), there is also an optional flag ( bool norm ) to draw normalized histograms to an areal of 1.
  
  // version 0: returns overlay of histogram objects from a single input file

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
  compare_histos("analyzed_files/combined/cafe_prod_Fe54_SRC_combined.root", hist_name_xbj_src, clr_src, hist_leg_SRC, "X_{bj}", "Counts","x-Bjorken, X_{bj}: Fe54 SRC", false);
  
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
  

  
  /*
  //--------------------------------
  // OVERLAY NUCLEI (1D HISTOGRAMS)
  //--------------------------------
  // brief: this function returns an overlay of the targets input by the user, with corresponding colors for any given histogram
  //        The histograms are read from analyzed .root files in a fixed location specified within the function.

  // NOTE: All .root files and summary (.csv) files used by these functions are assumed to be in a fixed location (set within the function)

  vector<string> tgt = {"LD2", "Be9", "B10", "B11", "C12"};
  vector<int> clr     = {2, 4, 6, 8, 9};    //  root cern color code ---> 2 : red,  4: blue, 6: magenta, 8: green, 9: purple
  overlay_nuclei(tgt, clr, "MF", "kin_plots/H_Pm", "Missing Momentum [GeV/c]", "Normalized Counts", "Missing Momentum (light nuclei)");
 
  */

  
  /*
  //-------------------------
  // HISTOGRAM SINGLE RATIOS
  //-------------------------
  // brief: this function returns any two histograms (A,B) and ratio, A/B. 
  //         The histograms (A,B) are read from analyzed combined .root files in a fixed location specified within the function.
  //         and A/B is calculated within the ratio function

  // NOTE: All .root files and summary (.csv) files used by these functions are assumed to be in a fixed location (set within the function)
    
  // Example of Use: 

  // declare histogram vector to retrieve histograms A: Ca48 MF,  B: Ca40 MF, and R: A/B, binned in Pmiss
  vector<TH1F*> hvec;

  // call function to get ratio
  // arguments:("targetA", "kinematicsA", "targetB", "kinematicsB", "histogram object", bool show_histos? )
  hvec = get_single_ratios("Ca48",       "MF",        "Ca40",       "MF",       "kin_plots/H_Pm",     true); 

  // access the histograms from the vector
  TH1F *A = &*hvec[0];  // get histA
  TH1F *B = &*hvec[1];  // get histB
  TH1F *R = &*hvec[2];  // get ratio of two histograms, A/B
  */
  
  /*
  // check that the retrieved histograms are as expected
  TCanvas *c = new TCanvas("c1","",800,1200);
  c->Divide(1,2);
  c->cd(1);
  A->Draw("histE0");
  B->Draw("sameshistE0");
  cout << "A Integral: " << A->Integral() << endl;
  cout << "B Integral: " << B->Integral() << endl;
  c->cd(2);
  R->Draw("E0");
  cout << "R Integral: " << R->Integral() << endl;
  */
 
  
}
