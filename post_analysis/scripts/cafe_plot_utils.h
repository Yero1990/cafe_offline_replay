/* 
   ------------------
   Author: C. Yero
   
   Brief: A compilation of CaFe utility functions
   for plotting histograms
   
   ------------------
*/ 

#include "cafe_summary_files_utils.h"


//____________________________________________________________________________________________________
void compare_histos( TString file1_path="path/to/file1.root", TString file2_path="path/to/file2.root", 
		    TString hist1="kin_plots/H_Pm", TString hist2="kin_plots/H_Pm",
		    TString xlabel="X-label [units]", TString ylabel="Y-label [units]", TString title="title",
		    TString hist1_leg="hist1_legend_title", TString hist2_leg="hist2_legend_title", bool norm=true)
{

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




TH1F* single_ratios(string tgtA="",  string kinA="", string tgtB="", string kinB="", string hist_name="", bool show_histos=false ){
  
  // brief: calculate single ratios of binned histograms (hist_name) for  target, kinematics): R = tgtA_kinA / tgtB_kinB  
  // using the analyzed CaFe *_combined.root files .  Example:  Ca48 MF / Ca40 MF,  Fe54 SRC / Ca48 SRC, etc.
  
 
  // define histogram scale variables 
  double scale_factor_A,  Q_A,  hms_trk_eff_A,  shms_trk_eff_A,  total_LT_A;
  double scale_factor_B, Q_B, hms_trk_eff_B, shms_trk_eff_B, total_LT_B;

  
  // get info from summary files (for scaling histograms)
  Q_A  = get("total_charge", tgtA.c_str(),  kinA.c_str());    
  Q_B = get("total_charge", tgtB.c_str(),  kinB.c_str());
  
  hms_trk_eff_A    = get("hms_trk_eff",  tgtA.c_str(), kinA.c_str() );
  hms_trk_eff_B   = get("hms_trk_eff",  tgtB.c_str(), kinB.c_str() );    
  
  shms_trk_eff_A   = get("shms_trk_eff", tgtA.c_str(), kinA.c_str() );
  shms_trk_eff_B  = get("shms_trk_eff", tgtB.c_str(), kinB.c_str() );
  
  total_LT_A       = get("total_live_time", tgtA.c_str(),  kinA.c_str());
  total_LT_B      = get("total_live_time", tgtB.c_str(),  kinB.c_str());
  
  scale_factor_A   = 1. / ( Q_A * hms_trk_eff_A * shms_trk_eff_A * total_LT_A ) ;
  scale_factor_B   = 1. / ( Q_B * hms_trk_eff_B * shms_trk_eff_B * total_LT_B ) ;
  
  // set .root file names 
  string fname_A = Form("analyzed_files_combined_pass1/cafe_prod_%s_%s_combined.root", tgtA.c_str(), kinA.c_str());
  string fname_B = Form("analyzed_files_combined_pass1/cafe_prod_%s_%s_combined.root", tgtB.c_str(), kinB.c_str());

  // define TFile
  TFile *file_A  = NULL;
  TFile *file_B = NULL;
    
  // define histogram object arrays
  TH1F *H_hist_A  = 0;
  TH1F *H_hist_B = 0;
  //TH1F *H_hist_R = 0;
    
  // read TFile
  file_A = new TFile(fname_A.c_str());
  file_B = new TFile(fname_B.c_str());

   
  // get histogram objects
  file_A->cd(); 
  file_A->GetObject(hist_name.c_str(), H_hist_A);

  // scale histogram appropiately
  H_hist_A->Scale(scale_factor_A);
  
    
  file_B->cd();
  file_B->GetObject(hist_name.c_str(), H_hist_B);

  // scale histogram appropiately
  H_hist_B->Scale(scale_factor_B);
  

  // calculate the ratio (A * scale_factor / (B * scale_factor) )
  TH1F *H_hist_R = (TH1F*)H_hist_A->Clone("H_hist_R");
  H_hist_R->Divide(H_hist_A, H_hist_B);

  if(show_histos) {

    int font_type = 132;
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetTitleFont(font_type, "");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFont(font_type);
    gStyle->SetLegendTextSize(0.05);
  
    // set histos aethetics
  H_hist_A->SetLineColor(kRed);
  H_hist_A->SetFillColorAlpha(kRed, 0.40);
  H_hist_A->SetFillStyle(3004);

  H_hist_B->SetLineColor(kBlue);
  H_hist_B->SetFillColorAlpha(kBlue, 0.40);
  H_hist_B->SetFillStyle(3005);

  H_hist_R->SetLineColor(kBlack);
  H_hist_R->SetFillColorAlpha(kBlack, 0.40);
  H_hist_R->SetFillStyle(3005);
  
  // set histogram titles/labels/font
  H_hist_A->SetTitle(H_hist_A->GetTitle());
  
  H_hist_A->GetXaxis()->SetLabelSize(0.04);
  H_hist_A->GetYaxis()->SetLabelSize(0.04);

  H_hist_A->SetLabelFont(font_type, "XY");
  
  H_hist_R->SetLabelFont(font_type, "XY");
  H_hist_R->GetXaxis()->SetLabelSize(0.04);
  H_hist_R->GetYaxis()->SetLabelSize(0.04);
  H_hist_R->SetMarkerStyle(8);
  H_hist_R->SetMarkerSize(0.8);
  H_hist_R->SetMarkerColor(kBlack);
    
    TCanvas *c = new TCanvas("c", "", 800,1200);
    c->Divide(1,2);

    c->cd(1);
    H_hist_A->Draw("histE0");
    H_hist_B->Draw("samehistE0");

    // create legend ( displays hist legend label and integral counts)
    TLegend *leg1 = new TLegend(0.14,0.89,0.25,0.65);
    double h1_I, h2_I;
    double h1_Ierr, h2_Ierr;
    double nbins = H_hist_A->GetNbinsX();  //Get total number of bins (excluding overflow)
    h1_I = H_hist_A->IntegralAndError(1, nbins, h1_Ierr);
    h2_I = H_hist_B->IntegralAndError(1, nbins, h2_Ierr);
    
    leg1->AddEntry(H_hist_A,Form("#splitline{%s %s}{Integral: %.3f #pm %.2f}", tgtA.c_str(), kinA.c_str(), h1_I, h1_Ierr), "f");
    leg1->AddEntry(H_hist_B,Form("#splitline{%s %s}{Integral: %.3f #pm %.2f}", tgtB.c_str(), kinB.c_str(), h2_I, h2_Ierr), "f");
    // draw legend
    leg1->Draw();

    
    c->cd(2);
    H_hist_R->Draw();

    
  }
  
  //return H_hist_A;
  //return H_hist_B;
  return H_hist_R;

  
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
  //TString tgt[5] = {"LD2", "Be9", "B10", "B11", "C12"};
  //int clr[5]     = {2, 4, 6, 8, 9}; 
  //overlay_nuclei(5, tgt, clr, "MF", "kin_plots/H_Pm", "Missing Momentum [GeV/c]", "Normalized Counts", "Missing Momentum (light nuclei)");
  //overlay_nuclei(5, tgt, clr, "SRC", "kin_plots/H_Pm", "Missing Momentum [GeV/c]", "Normalized Counts", "Missing Momentum (light nuclei)");

  TH1F *h1_Ca48 = 0;
  
  h1_Ca48 = single_ratios("Ca48", "MF", "Ca40", "MF",  "kin_plots/H_Pm", true);

  
}
