/* 
   ------------------
   Author: C. Yero
   
   Brief: A compilation of CaFe utility functions
   for plotting histograms
   
   ------------------
*/ 

#include "cafe_summary_files_utils.h"

//____________________________________________________________________________________________________
void compare_histos(
		    TString file1_path="path/to/file1.root", TString file2_path="path/to/file2.root", 
		    TString hist1="kin_plots/H_Pm", TString hist2="kin_plots/H_Pm",
		    TString xlabel="X-label [units]", TString ylabel="Y-label [units]", TString title="title",
		    TString hist1_leg="hist1_legend_title", TString hist2_leg="hist2_legend_title", bool norm=true) {

  /* brief: generic function to compare (overlay) two 1D histograms from different files (assuming they have same binning)
     
     The arguments are as follows:

     file1_path, file2_path: ROOTfile paths ( /path/to/file.root )
     hist1, hist2          : complete path to histogram objects in file (for example if they are in a sub-direcotry, then it must be specified ("path/to/hist_object")
     xlabel, ylabel, title : self-explanatory (axis labels and plot title)
     hist1_leg, hist2_leg  : histograms legend names that can help identify what the histogram being plotted is
     norm                  : boolean flag that if set to true, draws the histograms normalized to an area of 1

   */

  gStyle->SetOptStat(0);

  int font_type = 132;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetTitleFont(font_type, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(font_type);
  gStyle->SetLegendTextSize(0.03);
  

  //Open  ROOT files;
  TFile *file1 = NULL;
  TFile *file2 = NULL;

  file1 = new TFile(file1_path.Data());
  file2 = new TFile(file2_path.Data());

  // declare 1D histos
  TH1F *H_hist1 = 0;
  TH1F *H_hist2 = 0;

  // get histogram objects
  file1->cd();
  file1->GetObject(hist1.Data(), H_hist1);
  file2->cd();
  file2->GetObject(hist2.Data(), H_hist2);

  // set histos aethetics
  H_hist1->SetLineColor(kRed);
  H_hist1->SetFillColorAlpha(kRed, 0.40);
  H_hist1->SetFillStyle(3004);

  H_hist2->SetLineColor(kBlue);
  H_hist2->SetFillColorAlpha(kBlue, 0.40);
  H_hist2->SetFillStyle(3005);

  // set y-range
  H_hist1->GetYaxis()->SetRangeUser(0, H_hist1->GetMaximum()+0.6*H_hist1->GetMaximum());

  // set histogram titles/labels/font
  H_hist1->SetTitle(title);
  
  H_hist1->GetXaxis()->SetLabelSize(0.04);
  H_hist1->GetYaxis()->SetLabelSize(0.04);
  
  H_hist1->GetYaxis()->SetTitle(ylabel);
  H_hist1->GetXaxis()->SetTitle(xlabel);

  H_hist1->GetYaxis()->CenterTitle();
  H_hist1->GetXaxis()->CenterTitle();

  H_hist1->SetLabelFont(font_type, "XY");
  H_hist1->SetTitleFont(font_type, "XY");
  H_hist1->SetTitleSize(0.05, "XY");
  H_hist1->SetTitleOffset(1., "XY");


  TCanvas *c = new TCanvas("c", "c", 900, 700);

  H_hist1->Draw("histE0");
  H_hist2->Draw("sameshistE0");


  if(norm) {
    H_hist1->DrawNormalized("histE0");
    H_hist2->DrawNormalized("sameshistE0");
  }
  
  // create legend ( displays hist legend label and integral counts)
  TLegend *leg = new TLegend(0.14,0.89,0.25,0.78);
  double h1_I, h2_I;
  double h1_Ierr, h2_Ierr;
  double nbins = H_hist1->GetNbinsX();  //Get total number of bins (excluding overflow)
  h1_I = H_hist1->IntegralAndError(1, nbins, h1_Ierr);
  h2_I = H_hist2->IntegralAndError(1, nbins, h2_Ierr);
  
  leg->AddEntry(H_hist1,Form("%s | Integral: %.3f", hist1_leg.Data(), h1_I),"f");
  leg->AddEntry(H_hist2,Form("%s | Integral: %.3f", hist2_leg.Data(), h2_I));
  // draw legend
  leg->Draw();
  
  
}

//____________________________________________________________________________________________________
void compare_histos(
		    TH1F *H_hist1=0, TH1F *H_hist2 = 0, 
		    TString xlabel="X-label [units]", TString ylabel="Y-label [units]", TString title="title",
		    TString hist1_leg="hist1_legend_title", TString hist2_leg="hist2_legend_title", bool norm=true) {

  /* brief: generic function to compare (overlay) two 1D histograms (assuming they have same binning)
     (note this has the same name as function that takes file paths, this is referred to as function overloading, and allows
     the user to utilize different variations of the same function)

     The arguments are as follows:

     H_hist1, H_hist2      : pre-defined user histogram objects (must have same binning)
     xlabel, ylabel, title : self-explanatory (axis labels and plot title)
     hist1_leg, hist2_leg  : histograms legend names that can help identify what the histogram being plotted is
     norm                  : boolean flag that if set to true, draws the histograms normalized to an area of 1

   */

  gStyle->SetOptStat(0);

  int font_type = 132;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetTitleFont(font_type, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(font_type);
  gStyle->SetLegendTextSize(0.03);

  // set histos aethetics
  H_hist1->SetLineColor(kRed);
  H_hist1->SetFillColorAlpha(kRed, 0.40);
  H_hist1->SetFillStyle(3004);

  H_hist2->SetLineColor(kBlue);
  H_hist2->SetFillColorAlpha(kBlue, 0.40);
  H_hist2->SetFillStyle(3005);

  // set y-range
  H_hist1->GetYaxis()->SetRangeUser(0, H_hist1->GetMaximum()+0.6*H_hist1->GetMaximum());

  // set histogram titles/labels/font
  H_hist1->SetTitle(title);
  
  H_hist1->GetXaxis()->SetLabelSize(0.04);
  H_hist1->GetYaxis()->SetLabelSize(0.04);
  
  H_hist1->GetYaxis()->SetTitle(ylabel);
  H_hist1->GetXaxis()->SetTitle(xlabel);

  H_hist1->GetYaxis()->CenterTitle();
  H_hist1->GetXaxis()->CenterTitle();

  H_hist1->SetLabelFont(font_type, "XY");
  H_hist1->SetTitleFont(font_type, "XY");
  H_hist1->SetTitleSize(0.05, "XY");
  H_hist1->SetTitleOffset(1., "XY");


  TCanvas *c = new TCanvas("c", "c", 900, 700);

  H_hist1->Draw("histE0");
  H_hist2->Draw("sameshistE0");


  if(norm) {
    H_hist1->DrawNormalized("histE0");
    H_hist2->DrawNormalized("sameshistE0");
  }
  
  // create legend ( displays hist legend label and integral counts)
  TLegend *leg = new TLegend(0.14,0.89,0.25,0.78);
  double h1_I, h2_I;
  double h1_Ierr, h2_Ierr;
  double nbins = H_hist1->GetNbinsX();  //Get total number of bins (excluding overflow)
  h1_I = H_hist1->IntegralAndError(1, nbins, h1_Ierr);
  h2_I = H_hist2->IntegralAndError(1, nbins, h2_Ierr);
  
  leg->AddEntry(H_hist1,Form("%s | Integral: %.3f", hist1_leg.Data(), h1_I),"f");
  leg->AddEntry(H_hist2,Form("%s | Integral: %.3f", hist2_leg.Data(), h2_I));
  // draw legend
  leg->Draw();
  
  
}

//____________________________________________________________________________________________________
void overlay_nuclei(const int n, TString tgt[], int clr[], TString kin="", TString hist_name="",
		    TString xlabel="X-label [units]", TString ylabel="Y-label [units]", TString title="title"){

  /* 
     brief: cafe-specific plotting utility to overlay n nuclei histograms, assuming a generic filename that varies only with (target, kin)
     the user may change this generic file-name accordingly

     -----------
     arguments:
     -----------
     n : number of nucei to overlay  
     tgt []: string array to hold target names (must be consistent with target as found in file name)
     clr[] : integer array to hold plotting colors of each target ( colors represented by integers, see https://root.cern.ch/doc/master/classTColor.html )
     kin   : single string to select kinematic setting for specified targets ("MF" or "SRC")
     hist_name : histogram object name (if histogram is in ROOT file sub-directory, then it must also be specified)
     xlabel, ylabel, title : strings to set histogram axis and title labels

     -----------------------
     example of code usage:
     -----------------------

     TString tgt[5] = {"LD2", "Be9", "B10", "B11", "C12"};   // create array of targets to plot
     int clr[5]     = {2, 4, 6, 8, 9};  // 
     overlay_nuclei(5, tgt, clr, "MF", "kin_plots/H_Pm", "Missing Momentum [GeV/c]", "Normalized Counts", "Missing Momentum (light nuclei)");

  */
  
  // dont show stats box
  gStyle->SetOptStat(0);

  // set global plotting style
  int font_type = 132;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetTitleFont(font_type, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(font_type);
  gStyle->SetLegendTextSize(0.03);
  
  
  TString fname[n];
  TCanvas *c = new TCanvas(Form("C_%s", kin.Data()), "", 900, 700);
  TLegend *leg = new TLegend(0.14,0.89,0.25,0.78);

  // loop over each file name
  for (int i=0; i<n; i++){

    // generic file name with specific target, kinematic
    fname[i] = Form("analyzed_files_combined_pass1/cafe_prod_%s_%s_combined.root", tgt[i].Data(), kin.Data());

    
    cout << "fname = " << fname[i] << endl;
    cout << clr[i] << endl;
    // read TFile
    TFile *file = NULL;
    file = new TFile(fname[i].Data());
  
    TH1F *H_hist =0;
    
    file->cd();
    file->GetObject(hist_name.Data(), H_hist);

    // set histos aethetics
    H_hist->SetLineColor(clr[i]);
    H_hist->SetLineWidth(2);
    H_hist->SetMarkerStyle(8);
    H_hist->SetMarkerSize(0.8);
    H_hist->SetMarkerColor(clr[i]);
	

    //H_hist->SetFillColorAlpha(clr[i], 0.40);
    //H_hist->SetFillStyle(3002);
    
    // set y-range
    H_hist->GetYaxis()->SetRangeUser(0.1, H_hist->GetMaximum()+0.6*H_hist->GetMaximum());
    
    // set histogram titles/labels/font
    H_hist->SetTitle(title);
    
    H_hist->GetXaxis()->SetLabelSize(0.04);
    H_hist->GetYaxis()->SetLabelSize(0.04);
    
    H_hist->GetYaxis()->SetTitle(ylabel);
    H_hist->GetXaxis()->SetTitle(xlabel);
    
    H_hist->GetYaxis()->CenterTitle();
    H_hist->GetXaxis()->CenterTitle();
    
    H_hist->SetLabelFont(font_type, "XY");
    H_hist->SetTitleFont(font_type, "XY");
    H_hist->SetTitleSize(0.05, "XY");
    H_hist->SetTitleOffset(1., "XY");

    // changed to canvas and draw
    c->cd();
    if(n==0) { H_hist->DrawNormalized("histE0") ;}
    else{
      H_hist->DrawNormalized("sameshistE0");
    }

    // add legend entry
    leg->AddEntry(H_hist,Form("%s %s", tgt[i].Data(), kin.Data()),"f");

  }

  cout <<"show canvas " << endl;
  leg->Draw();
  c->Show();
  
}




void single_ratios(vector<string> tgt={}, string hist_name="" ){
  
  // brief: calculate single ratios of SRC / MF per target nuclei for a given historam
  // using the analyzed CaFe *_combined.root files

  // define filename string arrays
  string fname_MF[tgt.size()];
  string fname_SRC[tgt.size()];

  string fname_csv_MF[tgt.size()];
  string fname_csv_SRC[tgt.size()];
  
  // define TFile
  TFile *file_MF = NULL;
  TFile *file_SRC = NULL;

  // define histogram objects
  TH1F *H_hist_MF =0;
  TH1F *H_hist_SRC =0;
  
  // define hsitogram scale variables 
  double scale_factor_MF,  Q_MF,  hms_trk_eff_MF,  shms_trk_eff_MF,  total_LT_MF;
  double scale_factor_SRC, Q_SRC, hms_trk_eff_SRC, shms_trk_eff_SRC, total_LT_SRC;

    
  // loop over each file name
  for (int i=0; i<tgt.size(); i++){


    // set summary file names (for getting charge, eff, etc to be used in scaling the histos)
    fname_csv_MF[i]  = Form("summary_files_pass1/EmissCut_100MeV/cafe_prod_%s_MF_report_summary.csv", tgt[i].c_str());
    fname_csv_SRC[i] = Form("summary_files_pass1/EmissCut_100MeV/cafe_prod_%s_SRC_report_summary.csv", tgt[i].c_str());

    // get info from summary files (for scaling histograms)
    Q_MF  = get("total_charge", tgt[i].c_str(),  "MF");    
    Q_SRC = get("total_charge", tgt[i].c_str(),  "SRC");

    hms_trk_eff_MF    = get("hms_trk_eff",  tgt[i].c_str(),  "MF");
    hms_trk_eff_SRC   = get("hms_trk_eff",  tgt[i].c_str(),  "SRC");    

    shms_trk_eff_MF   = get("shms_trk_eff", tgt[i].c_str(),  "MF");
    shms_trk_eff_SRC  = get("shms_trk_eff", tgt[i].c_str(),  "SRC");
    
    total_LT_MF       = get("total_live_time", tgt[i].c_str(),  "MF");
    total_LT_SRC      = get("total_live_time", tgt[i].c_str(),  "SRC");

    scale_factor_MF    = 1. / ( Q_MF *  hms_trk_eff_MF  * shms_trk_eff_MF  * total_LT_MF ) ;
    scale_factor_SRC   = 1. / ( Q_SRC * hms_trk_eff_SRC * shms_trk_eff_SRC * total_LT_SRC ) ;

    
    
    // set .root file names 
    fname_MF[i] = Form("analyzed_files_combined_pass1/cafe_prod_%s_MF_combined.root", tgt[i].c_str());
    fname_SRC[i] = Form("analyzed_files_combined_pass1/cafe_prod_%s_SRC_combined.root", tgt[i].c_str());

    cout << fname_MF[i] << endl;
    cout << fname_SRC[i] << endl;
    
    // read TFile
    file_MF = new TFile(fname_MF[i].c_str());
    file_SRC = new TFile(fname_SRC[i].c_str());

    // get histogram objects
    file_MF->cd();
    file_MF->GetObject(hist_name.c_str(), H_hist_MF);

    file_SRC->cd();
    file_SRC->GetObject(hist_name.c_str(), H_hist_SRC);

	
    
  }

}



void cafe_plot_utils(){

  cout << "" << endl;
  cout << "" << endl;
  cout << "--------------------------------" << endl;
  cout << "" << endl;
  cout << "Welcome to CaFe Plot Utils" << endl;
  cout << "" << endl;
  cout << "Please call a plotting utility \nfunction from this script " << endl;
  cout << "--------------------------------" << endl;

  
  //  compare_histos("kin_plots/H_Pm", "kin_plots/H_Pm", "X-label [units]", "Y-label [units]", "title",
  //		 "hist1_legend_title", "hist2_legend_title", true);

  /* root cern color code;
     2 : red,  4: blue, 6: magenta, 8: green, 9: purple
   */

  
  // select which nuclei to overlay and at which kinematics
  TString tgt[5] = {"LD2", "Be9", "B10", "B11", "C12"};
  int clr[5]     = {2, 4, 6, 8, 9}; 
  overlay_nuclei(5, tgt, clr, "MF", "kin_plots/H_Pm", "Missing Momentum [GeV/c]", "Normalized Counts", "Missing Momentum (light nuclei)");
  overlay_nuclei(5, tgt, clr, "SRC", "kin_plots/H_Pm", "Missing Momentum [GeV/c]", "Normalized Counts", "Missing Momentum (light nuclei)");

    
}
