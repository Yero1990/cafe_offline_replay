void make_plots(){

  // e- angle: 6.8 deg
  TString data_fname="initial/cafe_replay_optics_16036_500000.root ";
  TString simc_fname="initial/cafe_heep_singles_kin0_rad.root";

  
  TFile *fdata = new TFile(data_fname, "READ");
  TFile *fsimc = new TFile(simc_fname, "READ");

  fdata->cd();
  TTree *T = (TTree*)fdata->Get("T");

  fsimc->cd();
  TTree *SNT = (TTree*)fsimc->Get("SNT");
  
  TCut data_cuts = "P.gtr.dp>0&&P.gtr.dp<22&&P.cal.etottracknorm>0.8&&P.kin.primary.x_bj>0.9&&P.kin.primary.x_bj<1.1&&g.evtyp==1";
  TCut simc_cuts = "Weight*(e_delta>0&&e_delta<22&&(Q2/(2.*0.938*nu))>0.8&&(Q2/(2.*0.938*nu)<1.1))";

  //TCut data_cuts = "P.gtr.dp>0&&P.gtr.dp<22&&P.cal.etottracknorm>0.8&&g.evtyp==1";
  //TCut simc_cuts = "Weight*(e_delta>0&&e_delta<22)";
  
  const int nplots = 13;
  TH1F *H_data[nplots];
  TH1F *H_simc[nplots];

  TCanvas *c[nplots];

 
  for(int i=0; i<nplots; i++){

    if(i!=6) continue;
    
    if(i==0) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS e- momentum", 900, 900);
      H_data[i]=new TH1F("H_shms_kf", "SHMS e- momentum", 100,9,10.5);
      H_simc[i]=new TH1F("H_shms_kf_simc", "SHMS e- momentum", 100,9,10.5);     
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.gtr.p>>H_shms_kf(100,9,10.5)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_pf/1000.>>H_shms_kf_simc(100,9,10.5)", simc_cuts, "normhistEsames");
      

    }
    

    
    if(i==1) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS #delta", 900, 700);
      H_data[i]=new TH1F("H_shms_delta", "", 100, 0., 22);
      H_simc[i]=new TH1F("H_shms_delta_simc", "", 100, 0., 22);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.gtr.dp>>H_shms_delta(100,0,22)", data_cuts, "normhistE");  
      fsimc->cd();
      SNT->Draw("e_delta>>H_shms_delta_simc(100,0,22)", simc_cuts, "normhistEsames");  
    }
    
    if(i==2) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS xptar", 900, 700);
      H_data[i]=new TH1F("H_shms_xptar", "", 100,-0.06,0.06);
      H_simc[i]=new TH1F("H_shms_xptar_simc", "", 100,-0.06,0.06);
     
      fdata->cd();
      T->Draw("P.gtr.th>>H_shms_xptar(100,-0.06,0.06)", data_cuts, "normhistE");
       c[i]->Modified();c[i]->Update();
      fsimc->cd();
       H_simc[i]->SetLineColor(kRed);
      SNT->Draw("e_xptar>>H_shms_xptar_simc(100,-0.06,0.06)", simc_cuts, "normhistEsames");  
      
    }
    
    if(i==3) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS yptar", 900, 700);
      H_data[i]=new TH1F("H_shms_yptar", "", 100,-0.06,0.06);
      H_simc[i]=new TH1F("H_shms_yptar_simc", "", 100,-0.06,0.06);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.gtr.ph>>H_shms_yptar(100,-0.06,0.06)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_yptar>>H_shms_yptar_simc(100,-0.06,0.06)", simc_cuts, "normhistEsames");  

    }

    if(i==4) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS ytar", 900, 700);
      H_data[i]=new TH1F("H_shms_ytar", "", 100,-2,2);
      H_simc[i]=new TH1F("H_shms_ytar_simc", "", 100,-2,2);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.gtr.y>>H_shms_ytar(100,-2,2)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_ytar>>H_shms_ytar_simc(100,-2,2)", simc_cuts, "normhistEsames");  

    }
      
    if(i==5) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      H_data[i]=new TH1F("H_shms_ecal", "SHMS Cal", 100, 0.5,2.6);
      fdata->cd();
      T->Draw("P.cal.etottracknorm>>H_shms_ecal", data_cuts);
      
    }

    if(i==6) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      H_data[i]=new TH1F("H_W", "Invariant Mass", 100, 0.85,1.05);
      H_simc[i]=new TH1F("H_W_simc", "Invariant Mass", 100, 0.85,1.05);
      H_simc[i]->SetLineColor(kRed);
      fsimc->cd();
      SNT->Draw("W>>H_W_simc(100,0.85,1.05)", simc_cuts, "normhistE");
    fdata->cd();
      T->Draw("P.kin.primary.W>>H_W(100,0.85,1.05)", data_cuts, "normhistEsames");
    }

    if(i==7) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      H_data[i]=new TH1F("H_xbj", "x-Bjorken", 100, 0.8,1.1);
      H_simc[i]=new TH1F("H_xbj_simc", "x-Bjorken", 100, 0.8,1.1);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.kin.primary.x_bj>>H_xbj(100,0.8,1.1)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("Q2/(2.*0.938*nu)>>H_xbj_simc(100,0.8,1.1)", simc_cuts, "normhistEsames");

    }

    if(i==8) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      H_data[i]=new TH1F("H_Q2", "Q2", 100, 0,4);
      H_simc[i]=new TH1F("H_Q2_simc", "Q2", 100, 0,4);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.kin.primary.Q2>>H_Q2(100,0,4)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("Q2>>H_Q2_simc(100,0,4)", simc_cuts, "normhistEsames");
    }

    if(i==9) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      H_data[i]=new TH1F("H_the", "e- angle", 100, 4.5,11);
      H_simc[i]=new TH1F("H_the_simc", "e- angle", 100, 4.5,11);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.kin.primary.scat_ang_deg>>H_the(100,4.5,11)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("theta_e*180/3.14>>H_the_simc(100,4.5,11)", simc_cuts, "normhistEsames");

    }

    if(i==10) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      H_data[i]=new TH1F("H_nu", "energy transfer, E-E'", 100, 0.,2);
      H_simc[i]=new TH1F("H_nu_simc", "energy transfer, E-E'", 100, 0.,2);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.kin.primary.nu>>H_nu(100,0.,2)", data_cuts, "normhistE"); 
      fsimc->cd();
      SNT->Draw("nu>>H_nu_simc(100,0.,2)", simc_cuts, "normhistEsames");
    }

    if(i==11) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      H_data[i]=new TH1F("H_tarx", "x-target (lab)'", 100, 0., 0.4);
      H_simc[i]=new TH1F("H_tarx_simc", "x-target (lab)", 100, 0., 0.4);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.react.x>>H_tarx(100,0.,0.4)", data_cuts, "normhistE"); 
      fsimc->cd();
      SNT->Draw("tar_x>>H_tarx_simc(100,0.,0.4)", simc_cuts, "normhistEsames");
    }
    
    if(i==12) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      H_data[i]=new TH1F("H_tary", "y-target (lab)'", 100, -0.2, 0.2);
      H_simc[i]=new TH1F("H_tary_simc", "y-target (lab)", 100, -0.2, 0.2);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.react.y>>H_tary(100,-0.2,0.2)", data_cuts, "normhistE"); 
      fsimc->cd();
      SNT->Draw("tar_y>>H_tary_simc(100,-0.2,0.2)", simc_cuts, "normhistEsames");
    }	
    
  }

}
