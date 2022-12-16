
// generic function to compare (overlay) two histograms from different files (assumoing they have same binning)
void compare_histos(TString hist1="kin_plots/H_Pm", TString hist2="kin_plots/H_Pm", TString xlabel="X-label [units]",
		    TString ylabel="Y-label [units]", TString title="title", TString hist1_leg="hist1_legend_title", TString hist2_leg="hist2_legend_title", bool norm=true) {

  gStyle->SetOptStat(0);

  int font_type = 132;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetTitleFont(font_type, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(font_type);
  gStyle->SetLegendTextSize(0.03);
  
  // select the two files from which to compare the two histos
  TString file1_path="cafe_prod_LH2_heep_coin_16965_-1_histos_newOptics.root";
  TString file2_path="cafe_heep_coin_rad_analyzed.root ";
    
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

// cafe-specific plotting utility to overlay n nuclei histograms, assuming a generic filename that varies only with (target, kin
void overlay_nuclei(const int n, TString tgt[], int clr[], TString kin="", TString hist_name="",
		    TString xlabel="X-label [units]", TString ylabel="Y-label [units]", TString title="title"){


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
  TCanvas *c = new TCanvas("c", "c", 900, 700);
  TLegend *leg = new TLegend(0.14,0.89,0.25,0.78);
  
  for (int i=0; i<n; i++){

    // generic file name with specific target, kinematic
    fname[i] = Form("cafe_prod_%s_%s.root", tgt[i].Data(), kin.Data());

    
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
    H_hist->SetFillColorAlpha(clr[i], 0.40);
    H_hist->SetFillStyle(3002);
    
    // set y-range
    H_hist->GetYaxis()->SetRangeUser(0, H_hist->GetMaximum()+0.6*H_hist->GetMaximum());
    
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
    if(n==0) { H_hist->Draw("histE0") ;}
    else{
      H_hist->Draw("sameshistE0");
    }

    // add legend entry
    leg->AddEntry(H_hist,Form("%s %s", tgt[i].Data(), kin.Data()),"f");

  }
  leg->Draw();
  
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

  // set umber of nuclei to overlay, which target (tgt), color (clr), and kinematics ("MF", "SRC")
  TString tgt[3] = {"LD2", "Be9", "B10"};
  int clr[3]     = {1, 2, 3}; 
  overlay_nuclei(3, tgt, clr, "MF", "kin_plots/H_Pm");
    
}
