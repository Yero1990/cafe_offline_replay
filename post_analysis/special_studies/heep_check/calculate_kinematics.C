#include "TMath.h"

void calculate_kinematics(){

  /* Brief: this code defines the calculated (e,e'p) elastics kinematics
     using minimum required initial information, e.g.
     (Ebeam, e- momentum) or (Ebeam, e- angle) or (e- momentum, e- angle)

     The quantities are:
     kf_calc (Eb, the)
     Pf_calc (Eb. the)
     thp_calc (Eb, the)
     
     Then the calculated-measured quantities are filled into histograms for fitting, e.g.,
     dkf = (kf_calc - kf_meas),
     dPf =  (Pf_calc - Pf_meas),
     dthp = (thp_calc - thp_meas)

     keep in mind that the measured quantities do have energy loss whereas
     the calculated quantities do not, and eloss will have to be subtarcted
     from the mesured quantities. To achieve this, calculated the difference between
     the reconstructed and thrown variables from SIMC, then a fit to that difference
     should give the numerical value of eloss

   */


  // coincidence
  TString data_fname="~/ROOTfiles/heep_coin_optim/step4/cafe_replay_optics_16962_-1.root";
  TString simc_fname="~/ROOTfiles/heep_coin_optim/step4/cafe_heep_coin_kin0_rad.root";

  // e- angle: 6.8 deg
  //TString data_fname="~/ROOTfiles/heep_singles/step1/cafe_replay_optics_16026_-1.root";
  //TString simc_fname="~/ROOTfiles/heep_singles/step1/cafe_heep_singles_kin2_rad.root";

  // e- angle: 7.495 deg
  //TString data_fname="~/ROOTfiles/heep_singles/step1/cafe_replay_optics_16028_-1.root";
  //TString simc_fname="~/ROOTfiles/heep_singles/step1/cafe_heep_singles_kin1_rad.root";

  
  // e- angle: 8.295 deg
  //TString data_fname="~/ROOTfiles/heep_singles/step2/cafe_replay_optics_16036_-1.root";
  //TString simc_fname="~/ROOTfiles/heep_singles/step2/cafe_heep_singles_kin0_rad.root";

  //Read ROOTfile
  TFile *fdata = new TFile(data_fname.Data(), "READ");
  TFile *fsimc = new TFile(simc_fname.Data(), "READ");
  
  //Get the data tree
  fdata->cd();
  TTree *data_tree = (TTree*)fdata->Get("T");
  Long64_t data_nentries = data_tree->GetEntries();

  fsimc->cd();
  TTree *simc_tree = (TTree*)fsimc->Get("SNT");
  Long64_t simc_nentries = simc_tree->GetEntries();


  
  // define constants
  Double_t Mp = 0.938272;
  Double_t Eb = 10.549355;
  Double_t dtr = TMath::Pi() / 180.;
  
  // branch variables
  Double_t weight;
  Double_t th_e;
  Double_t th_p;
  Double_t th_q;
  Double_t th_pq;
  Double_t xangle;
  Double_t Pf;
  Double_t Pf_x;
  Double_t Pf_z;
  Double_t Pm_x;
  Double_t Pm_z;
  Double_t Pm_x_simc;
  Double_t Pm_z_simc;
  
  Double_t kf;
  Double_t xbj;
  Double_t Q2;
  Double_t nu;
  Double_t q;
  Double_t q_x;
  Double_t q_z;
  Double_t Em;
  Double_t hms_delta;
  Double_t shms_delta;
  Double_t e_xptar;
  Double_t e_yptar;
  Double_t h_xptar;
  Double_t h_yptar;
  
  // user-defined variables
  Double_t thp_calc;
  Double_t Pf_calc, Pfx_calc, Pfz_calc;
  Double_t kf_calc;

  Double_t dthp, dPf, dkf; // (calculated-measured) variables
  
  // data histograms
  TH1F *H_kf_calc  = new TH1F("H_kf_calc", "", 100, 8.5,10);
  TH1F *H_kf_meas  = new TH1F("H_kf_meas", "", 100, 8.5,10);
  TH1F *H_dkf      = new TH1F("H_dkf", "", 100, -0.05,0.05);
  TH2F *H2_dkf_v_the  = new TH2F("H2_dkf_v_the", "", 100,5,11, 100, -0.05,0.05);

	          
  TH1F *H_Pf_calc  = new TH1F("H_Pf_calc", "", 100, 1.3,3.3);
  TH1F *H_Pf_meas  = new TH1F("H_Pf_meas", "", 100, 1.3,2.3);
  TH1F *H_dPf      = new TH1F("H_dPf", "", 100, -0.05,0.05);
  TH2F *H2_dPf_v_thp  = new TH2F("H2_dPf_v_thp", "", 100, 45, 53, 100, -0.05,0.05);
     
  TH1F *H_thp_calc = new TH1F("H_thp_calc", "", 100, 45,53);
  TH1F *H_thp_meas = new TH1F("H_thp_meas", "", 100, 45,53);
  TH1F *H_dthp     = new TH1F("H_dthp", "#theta_{p,calc}-#theta_{p,meas}; d#theta_{p} [rad]", 100, -0.015,0.015);

  // simc histograms
  TH1F *H_kf_calc_simc  = new TH1F("H_kf_calc_simc", "", 100, 8.5,10);
  TH1F *H_kf_meas_simc  = new TH1F("H_kf_meas_simc", "", 100, 8.5,10);
  TH1F *H_dkf_simc      = new TH1F("H_dkf_simc", "", 100, -0.05,0.05);
  TH2F *H2_dkf_v_the_simc  = new TH2F("H2_dkf_v_the_simc", "", 100,5,11, 100, -0.05,0.05);
	          
  TH1F *H_Pf_calc_simc  = new TH1F("H_Pf_calc_simc", "", 100, 1.3,3.3);
  TH1F *H_Pf_meas_simc  = new TH1F("H_Pf_meas_simc", "", 100, 1.3,2.3);
  TH1F *H_dPf_simc      = new TH1F("H_dPf_simc", "", 100, -0.05,0.05);
  TH2F *H2_dPf_v_thp_simc  = new TH2F("H2_dPf_v_thp_simc", "", 100, 45, 53, 100, -0.05,0.05);
	          
  TH1F *H_thp_calc_simc = new TH1F("H_thp_calc_simc", "", 100, 45,53);
  TH1F *H_thp_meas_simc = new TH1F("H_thp_meas_simc", "", 100, 45,53);
  TH1F *H_dthp_simc     = new TH1F("H_dthp_simc", "#theta_{p,calc}-#theta_{p,meas}; d#theta_{p} [rad]", 100, -0.015,0.015);

  TH1F *H_Pmx_simc = new TH1F("H_Pmx_simc", "", 100,-0.02,0.1);
  TH1F *H_Pmz_simc = new TH1F("H_Pmz_simc", "", 100,-0.02,0.1);

  TH2F *H_Pmx_vs_exptar_simc = new TH2F("H_Pmx_vs_exptar_simc", "", 100, -0.1, 0.1, 100,-0.02,0.05);
  TH2F *H_Pmz_vs_exptar_simc = new TH2F("H_Pmz_vs_exptar_simc", "", 100, -0.1, 0.1, 100,-0.02,0.05);
  
  // MOMENTUM COMPONENTS (TESTING)
  TH1F *H_thq  = new TH1F("H_thq", "", 100, 44,54);

  TH1F *H_Pf_x  = new TH1F("H_Pf_x", "", 100, 0.8,2);
  TH1F *H_Pf_z  = new TH1F("H_Pf_z", "", 100, 0.8,2);
  TH1F *H_q_x  = new TH1F("H_q_x", "", 100, 0.8,2);
  TH1F *H_q_z  = new TH1F("H_q_z", "", 100, 0.8,2);

  TH1F *H_Pm_x  = new TH1F("H_Pm_x", "", 100, -0.1,0.1);
  TH1F *H_Pm_z  = new TH1F("H_Pm_z", "", 100, -0.1,0.1);
  TH2F *H_Pmx_vs_exptar = new TH2F("H_Pmx_vs_exptar", "", 100, -0.1, 0.1, 100,-0.02,0.05);
  TH2F *H_Pmz_vs_exptar = new TH2F("H_Pmz_vs_exptar", "", 100, -0.1, 0.1, 100,-0.02,0.05);

  
  fdata->cd();
  data_tree->SetBranchAddress("P.kin.primary.x_bj",         &xbj    );
  data_tree->SetBranchAddress("P.kin.primary.scat_ang_rad", &th_e   );
  data_tree->SetBranchAddress("H.kin.secondary.xangle",     &xangle );
  data_tree->SetBranchAddress("H.kin.secondary.emiss",      &Em     );
  data_tree->SetBranchAddress("H.gtr.p",                    &Pf     );
  data_tree->SetBranchAddress("P.gtr.p",                    &kf     );
  data_tree->SetBranchAddress("H.gtr.dp",                   &hms_delta     );
  data_tree->SetBranchAddress("P.gtr.dp",                   &shms_delta     );
  data_tree->SetBranchAddress("P.gtr.th",                   &e_xptar     );
  data_tree->SetBranchAddress("P.gtr.ph",                   &e_yptar     );
  data_tree->SetBranchAddress("H.gtr.th",                   &h_xptar     );
  data_tree->SetBranchAddress("H.gtr.ph",                   &h_yptar     );

  

  // loop over data
  for(int ientry=0; ientry<data_nentries; ientry++)
    {
      
      data_tree->GetEntry(ientry);
     
      // define calculated quantities

      th_p = xangle - th_e;  //detected hadron angle for each particle

      // all calculated quantities must be derived from (Eb, th_e)
      kf_calc = Mp*Eb / (Mp + 2.*Eb*pow(sin(th_e/2.),2) );  // final e- momentum kf_calc(Eb, th_e)
      Pfx_calc = -kf_calc * sin(th_e);                      // final proton momentum x-component 
      Pfz_calc = Eb - kf_calc * cos(th_e);                  // final proton momentum z-component
      Pf_calc = sqrt( pow(Pfx_calc,2) + pow(Pfz_calc,2) );      // final proton momentum (co-planar) Pf_calc (Eb, th_e)
      thp_calc = atan(abs(Pfx_calc)/abs(Pfz_calc));                       // final proton angle, thp_calc (Eb, th_e)

      // define (calculated-measured) quantities
      dthp = thp_calc - th_p;  // radians
      dPf  = Pf_calc - Pf;     // GeV
      dkf  = kf_calc - kf;     // GeV


      // define cuts
      Bool_t base_cuts;
      Bool_t base_cuts_noaccp;
      Bool_t base_cuts_singles;
      Bool_t base_cuts_singles_noaccp;

      base_cuts = (xbj > 0.9) && (xbj < 1.1) && Em<0.1 && shms_delta>0 && abs(hms_delta)<10. && abs(e_xptar)<0.01 && abs(e_yptar)<0.01  ;
      base_cuts_noaccp = (xbj > 0.9) && (xbj < 1.1) && Em<0.1 && shms_delta>0 && abs(hms_delta)<10.;

      base_cuts_singles = (xbj > 0.9) && (xbj < 1.1) && shms_delta>0 && abs(e_xptar)<0.01 && abs(e_yptar)<0.01 ;
      base_cuts_singles_noaccp = (xbj > 0.9) && (xbj < 1.1) && shms_delta>0 ;

      if(base_cuts){
	// fill histograms
	H_kf_calc -> Fill(kf_calc);
	H_kf_meas -> Fill(kf);
	H_dkf     -> Fill(dkf);
	H2_dkf_v_the -> Fill(th_e/dtr, dkf);
	  
	H_Pf_calc -> Fill(Pf_calc);
	H_Pf_meas -> Fill(Pf);
	H_dPf     -> Fill(dPf);
	H2_dPf_v_thp ->Fill(th_p/dtr, dPf);
	  
	H_thp_calc -> Fill(thp_calc/dtr);
	H_thp_meas -> Fill(th_p/dtr);
	H_dthp     -> Fill(dthp);
      }
      
      if(base_cuts_noaccp){
	H2_dkf_v_the -> Fill(th_e/dtr, dkf);
	H2_dPf_v_thp ->Fill(th_p/dtr, dPf);
      }
       
    }




  fsimc->cd();
  simc_tree->SetBranchAddress("Weight",     &weight    );
  simc_tree->SetBranchAddress("nu",         &nu    );
  simc_tree->SetBranchAddress("Q2",         &Q2    );
  simc_tree->SetBranchAddress("q",         &q    );
  simc_tree->SetBranchAddress("theta_e",    &th_e  ); // radians
  simc_tree->SetBranchAddress("theta_p",    &th_p );  // radians
  simc_tree->SetBranchAddress("theta_pq",    &th_pq );  // radians
  simc_tree->SetBranchAddress("Pmx",      &Pm_x_simc     );
  simc_tree->SetBranchAddress("Pmz",      &Pm_z_simc     );

  simc_tree->SetBranchAddress("Em",      &Em     );
  simc_tree->SetBranchAddress("h_pf",                    &Pf     ); // MeV
  simc_tree->SetBranchAddress("e_pf",                    &kf     ); // MeV
  simc_tree->SetBranchAddress("h_delta",                   &hms_delta     );
  simc_tree->SetBranchAddress("e_delta",                   &shms_delta     );
  simc_tree->SetBranchAddress("e_xptar",                  &e_xptar     );
  simc_tree->SetBranchAddress("e_yptar",                  &e_yptar     );
  simc_tree->SetBranchAddress("h_xptar",                  &h_xptar     );
  simc_tree->SetBranchAddress("h_yptar",                  &h_yptar     );


  // loop over simc
  for(int ientry=0; ientry<simc_nentries; ientry++)
    {
      
      simc_tree->GetEntry(ientry);
     

      xbj = Q2/(2.*Mp*nu);
      kf = kf/1000.;  // to GeV
      Pf = Pf/1000.; // to GeV

      //  th_q = acos( (Eb - kf*cos(th_e))/q );  // [rad]     
      th_q = th_p-th_pq; //rad

      q_x = q*sin(th_q);
      q_z = q*cos(th_q);

      Pf_x = Pf*sin(th_p);
      Pf_z = Pf*cos(th_p);

      Pm_x = -1.*(Pf_x - q_x);
      Pm_z = -1.*(Pf_z - q_z);

      
    
      H_thq->Fill(th_q, weight);
      H_q_x->Fill(q_x,  weight);
      H_q_z->Fill(q_z,  weight);
      H_Pm_x->Fill(Pm_x, weight);
      H_Pm_z->Fill(Pm_z, weight);
      H_Pmx_vs_exptar->Fill(e_xptar, Pm_x, weight);
      H_Pmz_vs_exptar->Fill(e_xptar, Pm_z, weight);
      
      H_Pf_x->Fill(Pf_x, weight);
      H_Pf_z->Fill(Pf_z, weight);

      H_Pmx_simc->Fill(Pm_x_simc, weight);
      H_Pmz_simc->Fill(Pm_z_simc, weight);

      H_Pmx_vs_exptar_simc->Fill(e_xptar, Pm_x_simc, weight);
      H_Pmz_vs_exptar_simc->Fill(e_xptar, Pm_z_simc, weight);
      
      // all calculated quantities must be derived from (Eb, th_e)
      kf_calc = Mp*Eb / (Mp + 2.*Eb*pow(sin(th_e/2.),2) );  // final e- momentum kf_calc(Eb, th_e)
      Pfx_calc = -kf_calc * sin(th_e);                      // final proton momentum x-component 
      Pfz_calc = Eb - kf_calc * cos(th_e);                  // final proton momentum z-component
      Pf_calc = sqrt( pow(Pfx_calc,2) + pow(Pfz_calc,2) );      // final proton momentum (co-planar) Pf_calc (Eb, th_e)
      thp_calc = atan(abs(Pfx_calc)/abs(Pfz_calc));                       // final proton angle, thp_calc (Eb, th_e)

      // define (calculated-measured) quantities
      dthp = thp_calc - th_p;  // radians
      dPf  = Pf_calc - Pf;     // GeV
      dkf  = kf_calc - kf;     // GeV


      // define cuts
      Bool_t base_cuts;
      Bool_t base_cuts_noaccp;
      Bool_t base_cuts_singles;
      Bool_t base_cuts_singles_noaccp;
	    
      base_cuts = (xbj > 0.9) && (xbj < 1.1) && Em<0.1 && shms_delta>0 && abs(hms_delta)<10. && abs(e_xptar)<0.01 && abs(e_yptar)<0.01;
      base_cuts_noaccp = (xbj > 0.9) && (xbj < 1.1) && Em<0.1 && shms_delta>0 && abs(hms_delta)<10.;

      base_cuts_singles = (xbj > 0.9) && (xbj < 1.1) && shms_delta>0 && abs(e_xptar)<0.01 && abs(e_yptar)<0.01;
      base_cuts_singles_noaccp = (xbj > 0.9) && (xbj < 1.1) && shms_delta>0;

      if(base_cuts){
	// fill histograms
	H_kf_calc_simc -> Fill(kf_calc);
	H_kf_meas_simc -> Fill(kf);
	H_dkf_simc     -> Fill(dkf);
	
	H_Pf_calc_simc -> Fill(Pf_calc);
	H_Pf_meas_simc -> Fill(Pf);
	H_dPf_simc     -> Fill(dPf);
			
	H_thp_calc_simc -> Fill(thp_calc/dtr);
	H_thp_meas_simc -> Fill(th_p/dtr);
	H_dthp_simc     -> Fill(dthp); // radians
      }

       if(base_cuts_noaccp){
	 H2_dkf_v_the_simc -> Fill(th_e/dtr, dkf);
	 H2_dPf_v_thp_simc ->Fill(th_p/dtr, dPf);
       }
      
    }



  /*  
  // write histograms to file
  TFile *myfile = new TFile("calculated_kinematics_output_heep_singles_16962_8p3deg.root", "RECREATE");
  H_kf_calc  ->Write();
  H_kf_meas  ->Write();
  H_dkf      ->Write();
  H2_dkf_v_the ->Write();
  H_Pf_calc  ->Write();
  H_Pf_meas  ->Write(); 
  H_dPf      ->Write(); 
  H2_dPf_v_thp->Write();	     
  H_thp_calc ->Write();
  H_thp_meas ->Write();
  H_dthp     ->Write();


  H_kf_calc_simc  ->Write();
  H_kf_meas_simc  ->Write();
  H_dkf_simc      ->Write();
  H2_dkf_v_the_simc ->Write();
  H_Pf_calc_simc  ->Write();
  H_Pf_meas_simc  ->Write(); 
  H_dPf_simc      ->Write(); 
  H2_dPf_v_thp_simc ->Write();	     
  H_thp_calc_simc ->Write();
  H_thp_meas_simc ->Write();
  H_dthp_simc     ->Write();
  
  myfile->Close();
  */
  
  TCanvas *c_dfk = new TCanvas("c_dkf", "", 800,800);
  H_dkf->SetLineColor(kBlack);
  H_dkf_simc->SetLineColor(kRed);
  H_dkf_simc     ->DrawNormalized("histE");
  H_dkf          ->DrawNormalized("histEsames");

  TCanvas *c_Pf = new TCanvas("c_Pf", "", 800,800);
  H_dPf->SetLineColor(kBlack);
  H_dPf_simc->SetLineColor(kRed);
  H_dPf_simc     ->DrawNormalized("histE");
  H_dPf          ->DrawNormalized("histEsames");

  TCanvas *c_dthp = new TCanvas("c_dthp", "", 800,800);
  H_dthp->SetLineColor(kBlack);
  H_dthp_simc->SetLineColor(kRed);
  H_dthp_simc     ->DrawNormalized("histE");
  H_dthp          ->DrawNormalized("histEsames");


  TCanvas *c_p = new TCanvas("c_p", "", 1000,1200);
  c_p->Divide(2,3);

  c_p->cd(1);
  //H_thq->Fill(th_q, weight);
  H_q_x->Draw("histE0");

  c_p->cd(2);
  H_q_z->Draw("histE0");

  c_p->cd(3);
  H_Pmx_simc->SetLineColor(kRed);
  H_Pm_x->Draw("normhistE0");
  H_Pmx_simc->Draw("normhistE0sames");

  c_p->cd(4);
  H_Pmz_simc->SetLineColor(kRed);
  H_Pm_z->Draw("normhistE0");
  H_Pmz_simc->Draw("normhistE0sames");
  
  c_p->cd(5);
  H_Pf_x->Draw("histE0");

  c_p->cd(6);
  H_Pf_z->Draw("histE0");

  TCanvas *c_2d = new TCanvas("c_2d", "", 1000,1000);
  c_2d->Divide(2,2);

  c_2d->cd(1);
  H_Pmx_vs_exptar->Draw("colz");
  c_2d->cd(2);
  H_Pmz_vs_exptar->Draw("colz");

  c_2d->cd(3);
  H_Pmx_vs_exptar_simc->Draw("colz");
  c_2d->cd(4);
  H_Pmz_vs_exptar_simc->Draw("colz");

  
}
