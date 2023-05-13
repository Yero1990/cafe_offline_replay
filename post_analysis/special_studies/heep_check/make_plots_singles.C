void make_plots_singles(){

  // delta : data-simc
  /*             [rad]       [cm]        [GeV]       [deg]      [GeV]
    kin----------dyptar------dytar---------dE'--------dthe-------dW------
    kin0      
    kin1      
    kin2       
   */


  //------ cafe h(e,e') singles ------
  // e- angle: 6.8 deg
  TString data_fname="~/ROOTfiles/heep_singles/step1/cafe_replay_optics_16026_-1.root";
  TString simc_fname="~/ROOTfiles/heep_singles/step1/cafe_heep_singles_kin2_rad.root";

  // e- angle: 7.495 deg
  //TString data_fname="~/ROOTfiles/heep_singles/step1/cafe_replay_optics_16028_-1.root";
  //TString simc_fname="~/ROOTfiles/heep_singles/step1/cafe_heep_singles_kin1_rad.root";

  
  // e- angle: 8.295 deg
  //TString data_fname="~/ROOTfiles/heep_singles/step1/cafe_replay_optics_16036_-1.root";
  //TString simc_fname="~/ROOTfiles/heep_singles/step1/cafe_heep_singles_kin0_rad.root";

  
  //------ deuteron exp h(e,e') singles -----
  // e- angle: 7.704 deg (run 20867), delta_shms = +12 %
  //TString data_fname="~/ROOTfiles/cafe_replay_optics_20867_500000.root ";
  //TString simc_fname="~/ROOTfiles/d2_heep_scan_singles_rad_+12.root";
  
  // e- angle: 14.153 deg (run 20844), delta_shms = -8 %
  //TString data_fname="~/ROOTfiles/cafe_replay_optics_20844_500000.root ";
  //TString simc_fname="~/ROOTfiles/d2_heep_scan_singles_rad_-8.root";
  

  
  TFile *fdata = new TFile(data_fname, "READ");
  TFile *fsimc = new TFile(simc_fname, "READ");

  fdata->cd();
  TTree *T = (TTree*)fdata->Get("T");

  fsimc->cd();
  TTree *SNT = (TTree*)fsimc->Get("SNT");

  // for run 16962, (e,e'p) are mostly in SHMS angular range: xptar(P.gtr.th): (-0.015, 0.015) rad,  yptar(P.gtr.ph): (-0.01, 0.01) rad
  // therefore, if using singles run 16036, which was taken at the same kinematics, this range MUST be selected for W to line up between singles/coin data
  TCut data_cuts = "P.gtr.dp>0&&P.gtr.dp<22&&P.cal.etottracknorm>0.8&&P.kin.primary.x_bj>0.9&&P.kin.primary.x_bj<1.1&&g.evtyp==1&&abs(P.gtr.th)<0.01&&abs(P.gtr.ph)<0.01";
  TCut simc_cuts = "Weight*(e_delta>0&&e_delta<22&&(Q2/(2.*0.938*nu))>0.9&&(Q2/(2.*0.938*nu)<1.1)&&abs(e_xptar)<0.01&&abs(e_yptar)<0.01)";

  //TCut data_cuts = "P.gtr.dp>-10&&P.gtr.dp<22&&P.cal.etottracknorm>0.8&&P.kin.primary.x_bj>0.9&&P.kin.primary.x_bj<1.1&&g.evtyp==1";
  //TCut simc_cuts = "Weight*(e_delta>-10&&e_delta<22&&(Q2/(2.*0.938*nu))>0.9&&(Q2/(2.*0.938*nu)<1.1))";


  const int nplots = 17;
  TH1F *H_data[nplots];
  TH1F *H_simc[nplots];

  TCanvas *c[nplots];


  
  
  for(int i=0; i<nplots; i++){

    // only plot tarx,y,z
    //if( ((i!=11) && (i!=12) && (i!=15)) ) continue;

    // only plot xptar, yptar, ytar, delta
    if( (i!=1) && (i!=2) && (i!=3) && (i!=4) && (i!=16)) continue;

    // 2d correlations
    //if(i!=14) continue;

    // plot only kinematics (kf, th_e, Q2, xbj, nu, W)
    //if( (i!=0) &&  (i!=6) && (i!=9) && (i!=7) ) continue;
    
   if(i==0) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS e- momentum", 900, 900);
      fdata->cd();
      T->Draw("P.gtr.p>>H_shms_kf(100,9,10.5)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_pf/1000.>>H_shms_kf_simc(100,9,10.5)", simc_cuts, "normhistEsames");
      

    }
    

    
    if(i==1) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS #delta", 900, 700);
  
      fdata->cd();
      T->Draw("P.gtr.dp>>H_shms_delta(100,0,22)", data_cuts, "normhistE");  
      fsimc->cd();
      SNT->Draw("e_delta>>H_shms_delta_simc(100,0,22)", simc_cuts, "normhistEsames");  
    }
    
    if(i==2) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS xptar", 900, 700);
    
     
      fdata->cd();
      T->Draw("P.gtr.th>>H_shms_xptar(100,-0.06,0.06)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_xptar>>H_shms_xptar_simc(100,-0.06,0.06)", simc_cuts, "normhistEsames");  
      
    }
    
    if(i==3) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS yptar", 900, 700);
    
      fdata->cd();
      T->Draw("P.gtr.ph>>H_shms_yptar(200,-0.06,0.06)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_yptar>>H_shms_yptar_simc(200,-0.06,0.06)", simc_cuts, "normhistEsames");  

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
      fdata->cd();
      T->Draw("P.kin.primary.scat_ang_deg>>H_the(200,4.5,11)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("theta_e*180/3.14159265359>>H_the_simc(200,4.5,11)", simc_cuts, "normhistEsames");

    }

    if(i==10) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      fdata->cd();
      T->Draw("P.kin.primary.nu>>H_nu(200,0.4,1.4)", data_cuts, "normhistE"); 
      fsimc->cd();
      SNT->Draw("nu>>H_nu_simc(200,0.4,1.4)", simc_cuts, "normhistEsames");
    }

    if(i==11) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
     
      fdata->cd();
      T->Draw("P.react.x>>H_tarx(100,-0.4,0.5)", data_cuts, "normhistE");  // need negative sign applied in data to make directo comparisong to SIMC beam target X (different coord.)
      fsimc->cd();
      SNT->Draw("tar_x>>H_tarx_simc(100,-0.4,0.5)", simc_cuts, "normhistEsames");
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
    
      
      if(i==13) {

	// calculated energy transfer nu, as a function of the e- angle:
	// nu_calc = Eb**2 * (1 - cos (theta_e) ) / ( Mp + Eb * ( 1 - cos(theta_e) ) )

	// calculated e- scat angle as a function of nu
	// cos (the_calc) =  (Eb**2  - Eb * nu - nu * Mp) / ( Eb * (Eb-nu) )

	
	c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
	
	fdata->cd();

	c[i]->Divide(2,2);

	// overlay energy transfer, nu measured and calculated
	c[i]->cd(1);
	T->Draw("P.kin.primary.nu >> H_nu_meas_data(100, 0, 5)",   data_cuts); 
	T->Draw("10.549*10.549*(1. - cos(P.kin.primary.scat_ang_rad) ) / ( 0.938272 + 10.549*(1. - cos(P.kin.primary.scat_ang_rad) )) >> H_nu_calc_data(100,0,5)",   data_cuts, "sames"); 

	// plots nu_meas - nu_cal  vs. theta_e
	c[i]->cd(2);
	T->Draw("(P.kin.primary.nu - ( 10.549*10.549*(1. - cos(P.kin.primary.scat_ang_rad) ) / ( 0.938272 + 10.549*(1. - cos(P.kin.primary.scat_ang_rad) ) )  )):P.kin.primary.scat_ang_deg  >> H_dnu_vs_theta_e_data(100, 5, 12, 100,-0.2,0.2)",   data_cuts, "colz"); 

	// overlay e- scat. angle, the measured and calculated
	c[i]->cd(3);
	T->Draw("P.kin.primary.scat_ang_deg >> H_the_meas_data(100, 5, 12)",   data_cuts); 
	T->Draw("acos(( (10.549*10.549) - (10.549*P.kin.primary.nu) - (P.kin.primary.nu * 0.938272) ) / (10.549*(10.549 - P.kin.primary.nu)) )*180./3.14159265359 >> H_the_calc_data(100, 5, 12)",   data_cuts, "sames"); 

	c[i]->cd(4);
	T->Draw("(P.kin.primary.scat_ang_deg - acos(( (10.549*10.549) - (10.549*P.kin.primary.nu) - (P.kin.primary.nu * 0.938272) ) / (10.549*(10.549 - P.kin.primary.nu)) )*180./3.14159265359):P.kin.primary.nu >> H_dthe_vs_nu_data(100,0.6,2, 100, -0.3,0.3)", data_cuts, "colz");
	
	// ----------- SIMC-----
	TCanvas *c2_kin_corr = new TCanvas("c2_kin_corr", "", 900, 700);
	
	fsimc->cd();
	c2_kin_corr->Divide(2,2);

	// overlay energy transfer, nu measured and calculated
	c2_kin_corr->cd(1);
	SNT->Draw("nu >> H_nu_meas_simc(100,0,5)", simc_cuts, "histE");
	SNT->Draw("( 10.549*10.549*(1. - cos(theta_e) ) / ( 0.938272 + 10.549*(1. - cos(theta_e) )))  >> nu_calc_simc(100,0,5)", simc_cuts, "histEsames");

	// plots nu_meas - nu_cal  vs. theta_e
	c2_kin_corr->cd(2);
	SNT->Draw("(nu - ( 10.549*10.549*(1. - cos(theta_e) ) / ( 0.938272 + 10.549*(1. - cos(theta_e) ) ) )):theta_e*180/3.14159265359  >> H_dnu_simc(100, 5, 12, 100,-0.2,0.2)", simc_cuts, "colz");

	// overlay e- scat. angle, the measured and calculated
	c2_kin_corr->cd(3);
	SNT->Draw("theta_e * 180./TMath::Pi() >> H_the_meas_simc(100, 5, 12)",   simc_cuts); 
	SNT->Draw("acos(( (10.549*10.549) - (10.549 * nu) - (nu * 0.938272) ) / (10.549*(10.549 - nu)) )*180./3.14159265359 >> H_the_calc_simc(100, 5, 12)",   simc_cuts, "sames"); 

	c2_kin_corr->cd(4);
	SNT->Draw("(theta_e * 180./TMath::Pi() - acos(( (10.549*10.549) - (10.549*nu) - (nu * 0.938272) ) / (10.549*(10.549 - nu)) )*180./TMath::Pi()):nu >> H_dthe_vs_nu_simc(100,0.6,2, 100, -0.3,0.3)", simc_cuts, "colz");


	// NOTE: simc collimator calculatrion
	// SNT->Draw("(tar_x - e_xptar*e_zv*cos(6.8*3.14159265359/180.) + e_xptar*253):(e_ytar + e_yptar*253.-(0.019+40.*.01*0.052)*e_delta+(0.00019+40*.01*.00052)*e_delta*e_delta)>>(100,-15,15,100,-15,15)", "Weight*(e_delta>0)", "colz")

      }
    
      
    if(i==14) {


      TCanvas *c1_recon_corr = new TCanvas("c1_recon_corr", "recons correlations", 900, 700);      
      c1_recon_corr->Divide(3,2);

      // DATA
      fdata->cd();
      c1_recon_corr->cd(1);
      T->Draw("P.gtr.dp:P.gtr.th>>H_shms_delta_vs_xptar(100, -0.07, 0.07, 100,0,22)", data_cuts, "colz");  
      c1_recon_corr->cd(2);
      T->Draw("P.gtr.dp:P.gtr.ph>>H_shms_delta_vs_yptar(100, -0.07, 0.07, 100,0,22)", data_cuts, "colz");  
      c1_recon_corr->cd(3);
      T->Draw("P.gtr.dp:P.gtr.y>>H_shms_delta_vs_ytar(100, -2, 2, 100,0,22)", data_cuts, "colz");  
      
      

      // SIMC
      fsimc->cd();
      c1_recon_corr->cd(4);
      SNT->Draw("e_delta:e_xptar>>H_shms_delta_vs_xptar_simc(100, -0.07, 0.07, 100,0,22)", simc_cuts, "colz");  
      c1_recon_corr->cd(5);
      SNT->Draw("e_delta:e_yptar>>H_shms_delta_vs_yptar_simc(100, -0.07, 0.07, 100,0,22)", simc_cuts, "colz");  
      c1_recon_corr->cd(6);
      SNT->Draw("e_delta:e_ytar>>H_shms_delta_vs_ytar_simc(100, -2, 2, 100,0,22)", simc_cuts, "colz");  
      

    }
    

    if(i==15) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      
      fdata->cd();
      T->Draw("P.react.z>>H_tarz(100,-15,15)", data_cuts, "normhistE");  // need negative sign applied in data to make directo comparisong to SIMC beam target X (different coord.)
      fsimc->cd();
      //SNT->Draw("tar_z>>H_tarz_simc(100,-10,10)", simc_cuts, "normhistEsames");
      SNT->Draw("e_zv>>H_tarz_simc(100,-15,15)", simc_cuts, "normhistEsames");

    }


    if(i==16) {


      TCanvas *c1_recon_corr = new TCanvas("W_recon_corr", "W recons correlations", 900, 700);      
      c1_recon_corr->Divide(4,2);

      // DATA
      fdata->cd();
      c1_recon_corr->cd(1);
      T->Draw("P.kin.primary.W:P.gtr.th>>H_shms_W_vs_xptar(100, -0.07, 0.07, 100,0.9,1.)", data_cuts, "colz");  
      c1_recon_corr->cd(2);
      T->Draw("P.kin.primary.W:P.gtr.ph>>H_shms_W_vs_yptar(100, -0.07, 0.07, 100,0.9,1.)", data_cuts, "colz");  
      c1_recon_corr->cd(3);
      T->Draw("P.kin.primary.W:P.gtr.y>>H_shms_W_vs_ytar(100, -2, 2, 100,0.9,1.)", data_cuts, "colz");  
      c1_recon_corr->cd(4);
      T->Draw("P.kin.primary.W:P.gtr.dp>>H_shms_W_vs_delta(100, 0, 22, 100,0.9,1.)", data_cuts, "colz");  
      
      

      // SIMC
      fsimc->cd();
      c1_recon_corr->cd(5);
      SNT->Draw("W:e_xptar>>H_shms_W_vs_xptar_simc(100, -0.07, 0.07, 100,0.9,1.)", simc_cuts, "colz");  
      c1_recon_corr->cd(6);
      SNT->Draw("W:e_yptar>>H_shms_W_vs_yptar_simc(100, -0.07, 0.07, 100,0.9,1.)", simc_cuts, "colz");  
      c1_recon_corr->cd(7);
      SNT->Draw("W:e_ytar>>H_shms_W_vs_ytar_simc(100, -2, 2, 100,0.9,1.)", simc_cuts, "colz");  
      c1_recon_corr->cd(8);
      SNT->Draw("W:e_delta>>H_shms_W_vs_delta_simc(100, 0, 22, 100,0.9,1.)", simc_cuts, "colz");  
      

    }

    

    
  }
  
  

}
