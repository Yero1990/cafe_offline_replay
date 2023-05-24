void make_plots_coin(){



  // e- angle: 8.3 deg (run 16962)
  TString data_fname="~/ROOTfiles/heep_coin_optim/step4/cafe_replay_optics_16962_-1.root";
  //TString data_fname="~/ROOTfiles/heep_coin_optim/step1/cafe_replay_optics_16962_100000_eyptar_minus0p5mr.root";
  //TString data_fname="~/ROOTfiles/heep_coin_optim/step1/cafe_replay_optics_16962_100000_hyptar_plus1mr.root";
  TString simc_fname="~/ROOTfiles/heep_coin_optim/step4/cafe_heep_coin_kin0_rad.root";

  
  TFile *fdata = new TFile(data_fname, "READ");
  TFile *fsimc = new TFile(simc_fname, "READ");

  fdata->cd();
  TTree *T = (TTree*)fdata->Get("T");

  fsimc->cd();
  TTree *SNT = (TTree*)fsimc->Get("SNT");
  
  TCut data_cuts = "P.gtr.dp>0&&P.gtr.dp<22&&abs(H.gtr.dp)<10.&&P.cal.etottracknorm>0.8&&P.kin.primary.x_bj>0.9&&P.kin.primary.x_bj<1.1&&H.kin.secondary.emiss<0.1&&g.evtyp>=4&&abs(P.gtr.th)<0.01&&abs(P.gtr.ph)<0.01";
  TCut simc_cuts = "Weight*(e_delta>0&&e_delta<22&&abs(h_delta)<10.&&(Q2/(2.*0.938*nu))>0.9&&(Q2/(2.*0.938*nu)<1.1)&&Em<0.1&&abs(e_xptar)<0.01&&abs(e_yptar)<0.01)";

  //TCut data_cuts = "P.gtr.dp>0&&P.gtr.dp<22&&abs(H.gtr.dp)<10.&&P.cal.etottracknorm>0.8&&P.kin.primary.x_bj>0.9&&P.kin.primary.x_bj<1.1&&H.kin.secondary.emiss<0.1&&g.evtyp>=4&&abs(P.gtr.th)<0.01";
  //TCut simc_cuts = "Weight*(e_delta>0&&e_delta<22&&abs(h_delta)<10.&&(Q2/(2.*0.938*nu))>0.9&&(Q2/(2.*0.938*nu)<1.1)&&Em<0.1&&abs(e_xptar)<0.01)";
  
  
  const int nplots = 29;
  TH1F *H_data[nplots];
  TH1F *H_simc[nplots];

  TCanvas *c[nplots];


  
  
  for(int i=0; i<nplots; i++){


     // only plot tarx,y,z
    //if( ((i!=11) && (i!=12) && (i!=15)) ) continue;

    // only plot xptar,yptar,ytar,delta
    // if( ( (i!=1) && (i!=2) && (i!=3) && (i!=4) && (i!=16) && (i!=17) && (i!=18)) && (i!=19) && (i!=28) ) continue;

    //plot only kinematics (kf, th_e, Q2, xbj, nu, W, Em, Pmx,y,z, Pm,  pcal-pmeas, )
    if( (i!=0) &&  (i!=6) && (i!=9) && (i!=20) &&  (i!=21) && (i!=22) && (i!=23) && (i!=24) && (i!=25) && (i!=26) && (i!=27)  ) continue;

    //if((i!=26)) continue;
    
    
    if(i==0) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS e- momentum", 900, 900);
      fdata->cd();
      T->Draw("P.gtr.p>>H_shms_kf(300,9,10.5)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_pf/1000.>>H_shms_kf_simc(300,9,10.5)", simc_cuts, "normhistEsames");
      

    }
    
    
    if(i==1) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS #delta", 900, 700);
  
      fdata->cd();
      T->Draw("P.gtr.dp>>H_shms_delta(300,0,22)", data_cuts, "normhistE");  
      fsimc->cd();
      SNT->Draw("e_delta>>H_shms_delta_simc(300,0,22)", simc_cuts, "normhistEsames");  
    }
    
    if(i==2) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS xptar", 900, 700);
    
     
      fdata->cd();
      T->Draw("P.gtr.th>>H_shms_xptar(200,-0.06,0.06)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_xptar>>H_shms_xptar_simc(200,-0.06,0.06)", simc_cuts, "normhistEsames");  
      
    }
    
    if(i==3) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS yptar", 900, 700);
    
      fdata->cd();
      T->Draw("P.gtr.ph>>H_shms_yptar(200,-0.02,0.02)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_yptar>>H_shms_yptar_simc(200,-0.02,0.02)", simc_cuts, "normhistEsames");  

    }

    if(i==4) {
      c[i] = new TCanvas(Form("c%i",i), "SHMS ytar", 900, 700);
      H_data[i]=new TH1F("H_shms_ytar", "", 100,-2,2);
      H_simc[i]=new TH1F("H_shms_ytar_simc", "", 100,-2,2);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.gtr.y>>H_shms_ytar(200,-2,2)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("e_ytar>>H_shms_ytar_simc(200,-2,2)", simc_cuts, "normhistEsames");  

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
      T->Draw("P.kin.primary.scat_ang_deg>>H_the(100,7.5,9.5)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("theta_e*180/3.14>>H_the_simc(100,7.5,9.5)", simc_cuts, "normhistEsames");

    }

    if(i==10) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      H_data[i]=new TH1F("H_nu", "energy transfer, E-E'", 100, 0.,2);
      H_simc[i]=new TH1F("H_nu_simc", "energy transfer, E-E'", 100, 0.,2);
      H_simc[i]->SetLineColor(kRed);
      fdata->cd();
      T->Draw("P.kin.primary.nu>>H_nu(200,0.,2)", data_cuts, "normhistE"); 
      fsimc->cd();
      SNT->Draw("nu>>H_nu_simc(200,0.,2)", simc_cuts, "normhistEsames");
    }

    if(i==11) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
     
      fdata->cd();
      T->Draw("P.react.x>>H_tarx(100,-0.2,0.2)", data_cuts, "normhistE");  
      fsimc->cd();
      SNT->Draw("tar_x>>H_tarx_simc(100,-0.2,0.2)", simc_cuts, "normhistEsames");
    }
    
    if(i==12) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);

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
	T->Draw("acos(( (10.549*10.549) - (10.549*P.kin.primary.nu) - (P.kin.primary.nu * 0.938272) ) / (10.549*(10.549 - P.kin.primary.nu)) )*180./3.14 >> H_the_calc_data(100, 5, 12)",   data_cuts, "sames"); 

	c[i]->cd(4);
	T->Draw("(P.kin.primary.scat_ang_deg - acos(( (10.549*10.549) - (10.549*P.kin.primary.nu) - (P.kin.primary.nu * 0.938272) ) / (10.549*(10.549 - P.kin.primary.nu)) )*180./3.14):P.kin.primary.nu >> H_dthe_vs_nu_data(100,0.6,2, 100, -0.3,0.3)", data_cuts, "colz");
	
	// ----------- SIMC-----
	TCanvas *c2 = new TCanvas("c2", "", 900, 700);
	
	fsimc->cd();
	c2->Divide(2,2);

	// overlay energy transfer, nu measured and calculated
	c2->cd(1);
	SNT->Draw("nu >> H_nu_meas_simc(100,0,5)", simc_cuts, "histE");
	SNT->Draw("( 10.549*10.549*(1. - cos(theta_e) ) / ( 0.938272 + 10.549*(1. - cos(theta_e) )))  >> nu_calc_simc(100,0,5)", simc_cuts, "histEsames");

	// plots nu_meas - nu_cal  vs. theta_e
	c2->cd(2);
	SNT->Draw("(nu - ( 10.549*10.549*(1. - cos(theta_e) ) / ( 0.938272 + 10.549*(1. - cos(theta_e) ) ) )):theta_e*180/3.14  >> H_dnu_simc(100, 5, 12, 100,-0.2,0.2)", simc_cuts, "colz");

	// overlay e- scat. angle, the measured and calculated
	c2->cd(3);
	SNT->Draw("theta_e * 180./TMath::Pi() >> H_the_meas_simc(100, 5, 12)",   simc_cuts); 
	SNT->Draw("acos(( (10.549*10.549) - (10.549 * nu) - (nu * 0.938272) ) / (10.549*(10.549 - nu)) )*180./3.14 >> H_the_calc_simc(100, 5, 12)",   simc_cuts, "sames"); 

	c2->cd(4);
	SNT->Draw("(theta_e * 180./TMath::Pi() - acos(( (10.549*10.549) - (10.549*nu) - (nu * 0.938272) ) / (10.549*(10.549 - nu)) )*180./TMath::Pi()):nu >> H_dthe_vs_nu_simc(100,0.6,2, 100, -0.3,0.3)", simc_cuts, "colz");


	// NOTE: simc collimator calculatrion
	// SNT->Draw("(tar_x - e_xptar*e_zv*cos(6.8*3.14/180.) + e_xptar*253):(e_ytar + e_yptar*253.-(0.019+40.*.01*0.052)*e_delta+(0.00019+40*.01*.00052)*e_delta*e_delta)>>(100,-15,15,100,-15,15)", "Weight*(e_delta>0)", "colz")

      }
    
      
    if(i==14) {


      TCanvas *c1 = new TCanvas("c1", "recons correlations", 900, 700);      
      c1->Divide(3,2);

      // DATA
      fdata->cd();
      c1->cd(1);
      T->Draw("P.gtr.dp:P.gtr.th>>H_shms_delta_vs_xptar(100, -0.07, 0.07, 100,0,22)", data_cuts, "colz");  
      c1->cd(2);
      T->Draw("P.gtr.dp:P.gtr.ph>>H_shms_delta_vs_yptar(100, -0.07, 0.07, 100,0,22)", data_cuts, "colz");  
      c1->cd(3);
      T->Draw("P.gtr.dp:P.gtr.y>>H_shms_delta_vs_ytar(100, -2, 2, 100,0,22)", data_cuts, "colz");  
      
      

      // SIMC
      fsimc->cd();
      c1->cd(4);
      SNT->Draw("e_delta:e_xptar>>H_shms_delta_vs_xptar_simc(100, -0.07, 0.07, 100,0,22)", simc_cuts, "colz");  
      c1->cd(5);
      SNT->Draw("e_delta:e_yptar>>H_shms_delta_vs_yptar_simc(100, -0.07, 0.07, 100,0,22)", simc_cuts, "colz");  
      c1->cd(6);
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



  //--------------------------- proton side ---------------

  if(i==16) {
    c[i] = new TCanvas(Form("c%i",i), "HMS #delta", 900, 700);
    
    fdata->cd();
    T->Draw("H.gtr.dp>>H_hms_delta(200,-12,12)", data_cuts, "normhistE");  
    fsimc->cd();
    SNT->Draw("h_delta>>H_hms_delta_simc(200,-12,22)", simc_cuts, "normhistEsames");  
  }
  
  if(i==17) {
    c[i] = new TCanvas(Form("c%i",i), "HMS xptar", 900, 700);
    
    
    fdata->cd();
    T->Draw("H.gtr.th>>H_hms_xptar(200,-0.1,0.1)", data_cuts, "normhistE");
    fsimc->cd();
    SNT->Draw("h_xptar>>H_hms_xptar_simc(200,-0.1,0.1)", simc_cuts, "normhistEsames");  
    
  }
  
  if(i==18) {
    c[i] = new TCanvas(Form("c%i",i), "HMS yptar", 900, 700);
    
    fdata->cd();
    T->Draw("H.gtr.ph>>H_hms_yptar(200,-0.06,0.06)", data_cuts, "normhistE");
    fsimc->cd();
    SNT->Draw("h_yptar>>H_hms_yptar_simc(200,-0.06,0.06)", simc_cuts, "normhistEsames");  
    
  }
  
  if(i==19) {
    c[i] = new TCanvas(Form("c%i",i), "HMS ytar", 900, 700);
    fdata->cd();
    T->Draw("H.gtr.y>>H_hms_ytar(200,-7,7)", data_cuts, "normhistE");
    fsimc->cd();
    SNT->Draw("h_ytar>>H_hms_ytar_simc(200,-7,7)", simc_cuts, "normhistEsames");  
    
  }

  if(i==20) {
    c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
    fdata->cd();
    T->Draw("(H.kin.secondary.xangle-P.kin.primary.scat_ang_rad)*180/3.14>>H_thp(100,45,55)", data_cuts, "normhistE");
    fsimc->cd();
    SNT->Draw("theta_p*180/3.14>>H_thp_simc(100,45,55)", simc_cuts, "normhistEsames");
    
  }
  
  
  if(i==21) {
    c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
    fdata->cd();
    T->Draw("H.kin.secondary.emiss>>H_Em(100,-0.1,0.1)", data_cuts, "normhistE");
    fsimc->cd();
    SNT->Draw("Em>>H_Em_simc(100,-0.1,0.1)", simc_cuts, "normhistEsames");
    
  }
  
  
  if(i==22) {
    c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
    fdata->cd();
    T->Draw("H.kin.secondary.Prec_x>>H_Pmx(100,-0.1,0.1)", data_cuts, "normhistE");
    fsimc->cd();
    SNT->Draw("-1.*Pmx>>H_Pmx_simc(100,-0.1,0.1)", simc_cuts, "normhistEsames");
    
  }
  
  
  if(i==23) {
    c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
    fdata->cd();
    T->Draw("H.kin.secondary.Prec_y>>H_Pmy(100,-0.1,0.1)", data_cuts, "normhistE");
    fsimc->cd();
    SNT->Draw("Pmy>>H_Pmy_simc(100,-0.1,0.1)", simc_cuts, "normhistEsames");
    
  }
  
  if(i==24) {
    c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
    fdata->cd();
    T->Draw("H.kin.secondary.Prec_z>>H_Pmz(100,-0.1,0.1)", data_cuts, "normhistE");
    fsimc->cd();
    SNT->Draw("Pmz>>H_Pmz_simc(100,-0.1,0.1)", simc_cuts, "normhistEsames");
    
  }

  
  if(i==25) {
    c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
    fdata->cd();
    T->Draw("H.gtr.p>>H_Pf(200,1.5,2.1)", data_cuts, "normhistE");
    fsimc->cd();
    SNT->Draw("h_pf/1000.>>H_Pf_simc(200,1.5,2.1)", simc_cuts, "normhistEsames");
    
  }


    if(i==26) {

      // calculate dP = (Pcalc - Pmeas)   and dP / Pmeas : Pcalc (Eb, theta_p), Pmeas (delta),  Pcalc independent of delta
      // Pcalc formula only applies for the proton momentum of elastic hydrogen reactions h(e,e'p)

      // Pcalc =  2*Mp*Eb(Eb+Mp)*cos(theta_p) / ( Mp**2 + 2*Mp*Eb + Eb**2*sin^2(theta_p)  ) 

      
      c[i] = new TCanvas(Form("c%i",i), "", 1200, 800);
      
      fdata->cd();
      T->Draw("( 2*0.938272*10.549*(10.549+0.938272)*cos((H.kin.secondary.xangle-P.kin.primary.scat_ang_rad)) / ( pow(0.938272,2) + 2*0.938272*10.549 + pow(10.549,2)* pow(sin((H.kin.secondary.xangle-P.kin.primary.scat_ang_rad)),2) )  -  H.gtr.p ) >>H_dPf(100,-0.1,0.1)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("( 2*0.938272*10.549*(10.549+0.938272)*cos(theta_p) / ( pow(0.938272,2) + 2*0.938272*10.549 + pow(10.549,2)* pow(sin(theta_p),2) )  -  h_pf/1000. ) >>H_dPf_simc(100,-0.1,0.1)", simc_cuts, "normhistEsames");
      
      // ---------- dP vs. theta_p
      TCanvas *c_dPf = new TCanvas("c_dPf", "", 1200, 800);
      c_dPf->Divide(2,1);
      
      c_dPf->cd(1);
      fdata->cd();
      T->Draw("( 2*0.938272*10.549*(10.549+0.938272)*cos((H.kin.secondary.xangle-P.kin.primary.scat_ang_rad)) / ( pow(0.938272,2) + 2*0.938272*10.549 + pow(10.549,2)* pow(sin((H.kin.secondary.xangle-P.kin.primary.scat_ang_rad)),2) )  -  H.gtr.p ) : ((H.kin.secondary.xangle-P.kin.primary.scat_ang_rad)*180/3.14) >>H_dPf_vs_theta_p(100,45,55, 100,-0.2,0.2)", data_cuts, "colz");
      
      c_dPf->cd(2);
      fsimc->cd();
      SNT->Draw("( 2*0.938272*10.549*(10.549+0.938272)*cos(theta_p) / ( pow(0.938272,2) + 2*0.938272*10.549 + pow(10.549,2)* pow(sin(theta_p),2) )  -  h_pf/1000. ) : theta_p*180/3.14 >>H_dPf_vs_theta_p_simc(100,45,55,100,-0.2,0.2)", simc_cuts, "colz");
      
      
      //--------------kf_calc - fk_meas
      
      // calculate dkf = (kf_calc - kf_meas), calculated electron momentum    
      // kf_calc formula only applies for the electron of elastic hydrogen reactions h(e,e'p)
      // kf_calc = Mp * Eb / (Mp + 2.*Eb* sin^2(theta_e/2.) )
      
      TCanvas *c_dkf = new TCanvas("c_dkf", "", 1200, 800);
      fdata->cd();
      T->Draw("(( 0.938272*10.549 / (0.938272 + 2.*10.549 * pow(sin(P.kin.primary.scat_ang_rad/2.), 2) )) - P.gtr.p)>>H_dkf(100,-0.1,0.1)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("(( 0.938272*10.549 / (0.938272 + 2.*10.549 * pow(sin(theta_e/2.), 2) )) - e_pf/1000.) >>H_dkf_simc(100,-10,0.1)",simc_cuts, "normhistEsames" );
    
      TCanvas *c_2dkf = new TCanvas("c_2dkf", "", 1200, 800);
      c_2dkf->Divide(2,1);
      fdata->cd();
      c_2dkf->cd(1);
      T->Draw("(( 0.938272*10.549 / (0.938272 + 2.*10.549 * pow(sin(P.kin.primary.scat_ang_rad/2.), 2) )) - P.gtr.p):P.kin.primary.scat_ang_deg>>H_2dkf(100,6,10,100,-0.1,0.1)", data_cuts, "colz");
      fsimc->cd();
      c_2dkf->cd(2);
      SNT->Draw("(( 0.938272*10.549 / (0.938272 + 2.*10.549 * pow(sin(theta_e/2.), 2) )) - e_pf/1000.):(theta_e*180/3.14) >>H_2dkf_simc(100,6,10,100,-0.1,0.1)",simc_cuts, "colz" );
    
    }
    
    if(i==27) {
      c[i] = new TCanvas(Form("c%i",i), "", 900, 700);
      fdata->cd();
      T->Draw("H.kin.secondary.pmiss>>H_Pm(150,-0.01,0.05)", data_cuts, "normhistE");
      fsimc->cd();
      SNT->Draw("Pm>>H_Pm_simc(150,-0.01,0.05)", simc_cuts, "normhistEsames");
    
  }


    if(i==28) {


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
  
  //-------------------------------------------------------
  
}

}
