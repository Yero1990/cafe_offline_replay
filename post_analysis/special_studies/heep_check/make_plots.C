void make_plots(){

  TString data_fname="cafe_replay_sample_16975_30000.root";
  TString simc_fname="d2_pm120_jmlfsi_rad.root";

  
  TFile *fdata = new TFile(data_fname, "READ");
  TFile *fsimc = new TFile(simc_fname, "READ");

  fdata->cd();
  TTree *T = (TTree*)fdata->Get("T");

  fsimc->cd();
  TTree *SNT = (TTree*)fsimc->Get("SNT");
  
  TCut data_cuts = "P.gtr.dp>0&&P.gtr.dp<22&&P.cal.etottracknorm>0.8";
  TCut simc_cuts = "Weight*(e_delta>0&&e_delta<22)";
  
  const int nplots = 10;
  TH1F *H_data[nplots];
  TH1F *H_simc[nplots];

  TCanvas *c[nplots];

 
  for(int i=0; i<nplots; i++){

    if(i==0) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS e- momentum", 900, 900);
      H_data[i]=new TH1F("H_shms_kf", "SHMS e- momentum", 100,8,10);
      H_simc[i]=new TH1F("H_shms_kf_simc", "SHMS e- momentum", 100,8,10);     
      T->Draw("P.gtr.p>>H_shms_kf", data_cuts);
      SNT->Draw("e_pf/1000.>>H_shms_kf_simc", simc_cuts);
      H_data[i]->DrawNormalized();
      H_simc[i]->DrawNormalized("sames");
    }
    /*
    if(i==1) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS #delta", 900, 900);
      H_data[i]=new TH1F("H_shms_delta", "", 100,0,22);     
      T->Draw("P.gtr.dp>>H_shms_delta", data_cuts);  }

    if(i==2) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS xptar", 900, 900);
      H_data[i]=new TH1F("H_shms_xptar", "", 100,-0.1,0.1);     
      T->Draw("P.gtr.th>>H_shms_xptar", data_cuts);  }

    if(i==3) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS yptar", 900, 900);
      H_data[i]=new TH1F("H_shms_yptar", "", 100,-0.1,0.1);     
      T->Draw("P.gtr.ph>>H_shms_yptar", data_cuts);  }

    if(i==4) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS ytar", 900, 900);
      H_data[i]=new TH1F("H_shms_ytar", "", 100,-6,6);     
      T->Draw("P.gtr.y>>H_shms_ytar", data_cuts); }
      
    if(i==5) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 900);
      H_data[i]=new TH1F("H_shms_ecal", "SHMS Cal", 100, 0.5,2.6);
      T->Draw("P.cal.etottracknorm>>H_shms_ecal", data_cuts); }

    if(i==6) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 900);
      H_data[i]=new TH1F("H_W", "Invariant Mass", 100, 0.8,1.1);
      T->Draw("P.kin.primary.W>>H_W", data_cuts); }

    if(i==7) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 900);
      H_data[i]=new TH1F("H_xbj", "x-Bjorken", 100, 0.8,1.1);
      T->Draw("P.kin.primary.x_bj>>H_xbj", data_cuts); }

    if(i==8) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 900);
      H_data[i]=new TH1F("H_Q2", "Q2", 100, 0,4);
      T->Draw("P.kin.primary.Q2>>H_Q2", data_cuts); }

    if(i==9) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 900);
      H_data[i]=new TH1F("H_the", "e- angle", 100, 6,11);
      T->Draw("P.kin.primary.scat_ang_deg>>H_the", data_cuts); }

    if(i==10) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 900);
      H_data[i]=new TH1F("H_nu", "energy transfer, E-E'", 100, 0,4);
      T->Draw("P.kin.primary.nu>>H_nu", data_cuts); }
    */

    
  }


}
