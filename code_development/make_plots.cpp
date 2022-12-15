#include "plot_utils.C"

void make_plots(TString hist1=""){

  TString file1_path="../../cafe_prod_LH2_heep_coin_16965_-1_histos_newOptics.root";
  TString file2_path="";
    
  //Open  ROOT files;
  TFile *file1 = NULL;
  TFile *file2 = NULL;

  file1 = new TFile(file1_path.Data());
  TH1F *H_f1 = 0;

  file1->cd();
  TString hist1_name=Form("%s", hist1.Data());
  file1->GetObject(hist1_name.Data(), H_f1);

  plot_hist(H_f1, "", "Counts", "", "");
  
  
}
