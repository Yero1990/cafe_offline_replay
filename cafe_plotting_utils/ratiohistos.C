
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
// This file when produces per-proton and per-neutron csv files that will later be used by python code to make nice histograms
// You need to manually comment and uncomment the normalization sections to switch between per-proton and per-neutron
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TImageDump.h"
#include <stdio.h>/////////////////////////
#include "Histo.h"
#include "Cuts.h"
#include "parse_utils.h"
#include "hist_utils.h"
using namespace std;

//-----------------------
// Constants
//------------------------
Double_t pi = 3.141592654;
Double_t dtr = pi/180.;
Double_t MP = 0.938272; //Proton Mass GeV
Double_t MD = 1.87561; //GeV
Double_t MN = 0.939566; //Neutron Mass GeV
Double_t me = 0.000510998; //GeV

void thingy(
	TH1F*, TH1F*, TH1F*,
	TH1F*, TH1F*, TH1F*,
	TH1F*, TH1F*,
	TH1F*, TH1F*,
	TFile*);

void thingpap(TH1F*, TH1F*, TH1F*, TFile*);
void thingpap2(TH1F*, TH1F*, TH1F*, const char*, const char*, TFile*);

void Draw1dE(TH1F*, const char*, const char*, TFile*);
void Draw2d(TH2F*, const char*, const char*, TFile*);
void Draw2de(TH2F*, const char*, const char*, TFile*);
void Draw2dh(TH2F*, const char*, const char*, TFile*);
void Draw2d2(TH2F*, const char*, const char*, TFile*);





//*****************************************************************************************************************************************************************************
//Begin Main
//*****************************************************************************************************************************************************************************
void ratiohistos() 
{
	
	//File paths
	TFile* outROOT = new TFile("outHist_ratios.root", "RECREATE");
	
	TString base = "Sec_Kin/H1_pmiss_ratio";
	TString base1 = "Sec_Kin/H1_pmiss_ratio_pid";
	TString base2 = "Sec_Kin/H1_pmiss_ratio_pid_acc";
	
	TString base211 = "Sec_Kin/H1_Ef_pid_acc_kin";
	TString base221 = "Sec_Kin/H1_Pf_pid_acc_kin";
	TString base212 = "Sec_Kin/H1_Ef_pid_acc_kin_mf";
	TString base222 = "Sec_Kin/H1_Pf_pid_acc_kin_mf";
	TString base23 = "Prim_Kin/H1_th_e_pid_acc_kin";
	TString base24 = "Sec_Kin/H1_th_p_pid_acc_kin";
	TString base25 = "Prim_Kin/H1_q_pid_acc_kin";
	TString base26 = "Prim_Kin/H1_nu_pid_acc_kin";


	TString base3 = "Sec_Kin/H1_pmiss_ratio_pid_acc_kin";
	TString base4 = "Prim_Kin/H1_xbj_ratio_pid_acc_kin_full";
	TString base5 = "Prim_Kin/H1_Q2_ratio_pid_acc_kin_full";
	
	TString base6 = "Sec_Kin/H1_pmiss_ratio_pid_acc_kin_full";
	TString base7 = "Prim_Kin/H1_thrq_ratio_pid_acc_kin_full";
	TString base8 = "Prim_Kin/H1_W_ratio_pid_acc_kin_full";
	TString base22 = "Prim_Kin/H1_W_ratio_pid_acc_kin_full_mf";
	TString base20 = "Prim_Kin/H1_Em_ratio_pid_acc_kin_full";

	TString base9 = "Sec_Kin/H1_pmiss_ratio_pid_acc_kin_full_mf";
	TString base10 = "Prim_Kin/H1_Q2_ratio_pid_acc_kin_full_mf";
	TString base11 = "Prim_Kin/H1_Em_ratio_pid_acc_kin_full_mf";

	TString base12 = "PID/H1_ep_ctime_ratio_pid_acc_kin_full";
	TString base21 = "PID/H1_ep_ctime_ratio_pid_acc_kin_full_mf";

	TString base13 = "PID/H1_pCalEtotTrkNorm_pid_acc_kin";//ratio_full
	TString base131 = "PID/H1_pCalEtotTrkNorm_ratio_pid_acc_kin_full_mf";
	TString base14 = "HMS_Accp/H1_hdelta_ratio_pid_acc_kin_full";
	TString base141 = "HMS_Accp/H1_hdelta_ratio_pid_acc_kin_full_mf";
	TString base15 = "SHMS_Accp/H1_edelta_ratio_pid_acc_kin_full";
	TString base151 = "SHMS_Accp/H1_edelta_ratio_pid_acc_kin_full_mf";

	TString base16 = "Prim_Kin/H2_Em_Pm_ratio_pid_acc_kin_full";
	TString base17 = "Prim_Kin/H2_xbj_Q2_ratio_pid_acc_kin_full";

	TString base18 = "HMS_Accp/H2_hXColl_hYColl_ratio_pid_acc_kin_full";
	TString base181 = "HMS_Accp/H2_hXColl_hYColl_ratio_pid_acc_kin_full_mf";
	TString base19 = "SHMS_Accp/H2_eXColl_eYColl_ratio_pid_acc_kin_full";
	TString base191 = "SHMS_Accp/H2_eXColl_eYColl_ratio_pid_acc_kin_full_mf";

	TString base192 = "Prim_Kin/H2_Em_Pm_ratio_pid_acc_kin_full";
	TString base193 = "Prim_Kin/H2_Em_Pm_ratio_pid_acc_kin_full2";
	TString base194 = "Prim_Kin/H2_Em_Pm_ratio_pid_acc_kin_full_mf";
	TString base195 = "Prim_Kin/H2_Em_Pm_ratio_pid_acc_kin_full_mf2";

	

	
	const Int_t NBins2 = 10;
	/*Double_t edges2[NBins2 + 1] = {
		0.0, 0.01, 0.02, 0.03, 0.04,
		0.05, 0.06, 0.07, 0.08, 0.09,
		0.10, 0.11, 0.12, 0.13, 0.14,
		0.15, 0.16, 0.17, 0.18, 0.19,
		0.20, 0.21, 0.22, 0.23, 0.24,
		0.25, 0.26, 0.27, 0.375, 0.395,
		0.415, 0.445, 0.475, 0.505, 0.535,
		0.565, 0.595, 0.625, 0.655, 0.685,
		0.700, 0.71//42
	};*/
	Double_t edges2[NBins2 + 1] = {
		0.3, 0.375, 0.415625, 0.45625, 0.496875,
		0.5375, 0.578125, 0.61875, 0.659375, 0.7,
		0.75//11
	};
	
	TH1F* H1_Ca40_sim = new TH1F("H1_Ca40_sim", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2);
	TH1F* H1_Ca48_sim = new TH1F("H1_Ca48_sim", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2);
	TH1F* H1_Fe54_sim = new TH1F("H1_Fe54_sim", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2);
	TH1F* H1_C12_sim = new TH1F("H1_C12_sim", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2);
	

	//TTree* t = new TTree("t", "tree from C12_Mom.csv");
	//t->ReadFile("C12_Mom.csv", "K/D:RHOKP:DRHOKP");
	//t->Draw("RHOKP : K", "", "LogY");









	float Be9_MF_Ef3 = 0;
	float B10_MF_Ef3 = 0;
	float B11_MF_Ef3 = 0;
	float C12_MF_Ef3 = 0;
	float Ca40_MF_Ef3 = 0;
	float Ca48_MF_Ef3 = 0;
	float Fe54_MF_Ef3 = 0;
	float Au197_MF_Ef3 = 0;
	float Be9_MF_Pf3 = 0;
	float B10_MF_Pf3 = 0;
	float B11_MF_Pf3 = 0;
	float C12_MF_Pf3 = 0;
	float Ca40_MF_Pf3 = 0;
	float Ca48_MF_Pf3 = 0;
	float Fe54_MF_Pf3 = 0;
	float Au197_MF_Pf3 = 0;
	float Be9_MF_Q23 = 0;
	float B10_MF_Q23 = 0;
	float B11_MF_Q23 = 0;
	float C12_MF_Q23 = 0;
	float Ca40_MF_Q23 = 0;
	float Ca48_MF_Q23 = 0;
	float Fe54_MF_Q23 = 0;
	float Au197_MF_Q23 = 0;
	float C12_MF_X3 = 0;
	float Ca40_MF_X3 = 0;
	float Ca48_MF_X3 = 0;
	float Fe54_MF_X3 = 0;
	

	//MF D2
	TString MF_D2 = "Mf/pass4/D2/Coin/MF_D2_Results.root";
	TFile* inROOT_MF_D2 = new TFile(MF_D2.Data(), "READ");
	TTree* inputtree_MF_D2 = (TTree*)inROOT_MF_D2->Get("T");
	TFile* file_MF_D2 = NULL;
	file_MF_D2 = new TFile(MF_D2.Data());
	TH1F* H1_pmiss_ratio_MF_D2 = 0;
	TH1F* H1_pmiss_ratio_MF_D2_pid = 0;
	TH1F* H1_pmiss_ratio_MF_D2_pid_acc = 0;
	TH1F* H1_pmiss_ratio_MF_D2_pid_acc_kin = 0;
	file_MF_D2->GetObject(base.Data(), H1_pmiss_ratio_MF_D2);
	file_MF_D2->GetObject(base1.Data(), H1_pmiss_ratio_MF_D2_pid);
	file_MF_D2->GetObject(base2.Data(), H1_pmiss_ratio_MF_D2_pid_acc);
	file_MF_D2->GetObject(base3.Data(), H1_pmiss_ratio_MF_D2_pid_acc_kin);
	TH1F* H1_W_ratio_MF_D2_pid_acc_kin_full = 0;
	file_MF_D2->GetObject(base8.Data(), H1_W_ratio_MF_D2_pid_acc_kin_full);
	TH1F* H1_xbj_ratio_MF_D2_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_MF_D2_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_MF_D2_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_MF_D2_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_D2_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_D2_pid_acc_kin_full_mf = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_MF_D2_pid_acc_kin_full_mf = 0;
	TH1F* H1_hdelta_ratio_MF_D2_pid_acc_kin_full_mf = 0;
	TH1F* H1_edelta_ratio_MF_D2_pid_acc_kin_full_mf = 0;
	file_MF_D2->GetObject(base12.Data(), H1_epctime_ratio_MF_D2_pid_acc_kin_full);
	file_MF_D2->GetObject(base21.Data(), H1_epctime_ratio_MF_D2_pid_acc_kin_full_mf);
	file_MF_D2->GetObject(base131.Data(), H1_pCalEtotTrkNorm_ratio_MF_D2_pid_acc_kin_full_mf);
	file_MF_D2->GetObject(base141.Data(), H1_hdelta_ratio_MF_D2_pid_acc_kin_full_mf);
	file_MF_D2->GetObject(base151.Data(), H1_edelta_ratio_MF_D2_pid_acc_kin_full_mf);
	file_MF_D2->GetObject(base4.Data(), H1_xbj_ratio_MF_D2_pid_acc_kin_full);
	file_MF_D2->GetObject(base5.Data(), H1_Q2_ratio_MF_D2_pid_acc_kin_full);
	file_MF_D2->GetObject(base6.Data(), H1_pmiss_ratio_MF_D2_pid_acc_kin_full);
	file_MF_D2->GetObject(base7.Data(), H1_thrq_ratio_MF_D2_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_MF_D2_pid_acc_kin_full_mf = 0;
	file_MF_D2->GetObject(base212.Data(), H1_Ef_ratio_MF_D2_pid_acc_kin_full_mf);
	TH1F* H1_Pf_ratio_MF_D2_pid_acc_kin_full_mf = 0;
	file_MF_D2->GetObject(base222.Data(), H1_Pf_ratio_MF_D2_pid_acc_kin_full_mf);
	TH1F* H1_the_ratio_MF_D2_pid_acc_kin_full = 0;
	file_MF_D2->GetObject(base23.Data(), H1_the_ratio_MF_D2_pid_acc_kin_full);
	TH1F* H1_thp_ratio_MF_D2_pid_acc_kin_full = 0;
	file_MF_D2->GetObject(base24.Data(), H1_thp_ratio_MF_D2_pid_acc_kin_full);
	TH1F* H1_q_ratio_MF_D2_pid_acc_kin_full = 0;
	file_MF_D2->GetObject(base25.Data(), H1_q_ratio_MF_D2_pid_acc_kin_full);
	TH1F* H1_nu_ratio_MF_D2_pid_acc_kin_full = 0;
	file_MF_D2->GetObject(base26.Data(), H1_nu_ratio_MF_D2_pid_acc_kin_full);


	//MF Be9
	TString MF_Be9 = "Mf/pass4/Be9/Coin/MF_Be9_Results.root";
	TFile* inROOT_MF_Be9 = new TFile(MF_Be9.Data(), "READ");
	TTree* inputtree_MF_Be9 = (TTree*)inROOT_MF_Be9->Get("T");
	TFile* file_MF_Be9 = NULL;
	file_MF_Be9 = new TFile(MF_Be9.Data());
	TH1F* H1_pmiss_ratio_MF_Be9 = 0;
	TH1F* H1_pmiss_ratio_MF_Be9_pid = 0;
	TH1F* H1_pmiss_ratio_MF_Be9_pid_acc = 0;
	TH1F* H1_pmiss_ratio_MF_Be9_pid_acc_kin = 0;
	file_MF_Be9->GetObject(base.Data(), H1_pmiss_ratio_MF_Be9);
	file_MF_Be9->GetObject(base1.Data(), H1_pmiss_ratio_MF_Be9_pid);
	file_MF_Be9->GetObject(base2.Data(), H1_pmiss_ratio_MF_Be9_pid_acc);
	file_MF_Be9->GetObject(base3.Data(), H1_pmiss_ratio_MF_Be9_pid_acc_kin);
	TH1F* H1_pmiss_ratio_MF_Be9_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_MF_Be9_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_MF_Be9_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Be9_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf = 0;
	TH1F* H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf = 0;
	TH1F* H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf = 0;
	TH1F* H1_W_ratio_MF_Be9_pid_acc_kin_full = 0;
	file_MF_Be9->GetObject(base12.Data(), H1_epctime_ratio_MF_Be9_pid_acc_kin_full);
	file_MF_Be9->GetObject(base21.Data(), H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf);
	file_MF_Be9->GetObject(base131.Data(), H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf);
	file_MF_Be9->GetObject(base141.Data(), H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf);
	file_MF_Be9->GetObject(base151.Data(), H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf);
	file_MF_Be9->GetObject(base9.Data(), H1_pmiss_ratio_MF_Be9_pid_acc_kin_full);
	file_MF_Be9->GetObject(base10.Data(), H1_Q2_ratio_MF_Be9_pid_acc_kin_full);
	file_MF_Be9->GetObject(base11.Data(), H1_Em_ratio_MF_Be9_pid_acc_kin_full);
	file_MF_Be9->GetObject(base22.Data(), H1_W_ratio_MF_Be9_pid_acc_kin_full);
	
	TH1F* H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf = 0;
	file_MF_Be9->GetObject(base212.Data(), H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf);
	TH1F* H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf = 0;
	file_MF_Be9->GetObject(base222.Data(), H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf);
	
	TH1F* H1_the_ratio_MF_Be9_pid_acc_kin_full = 0;
	file_MF_Be9->GetObject(base23.Data(), H1_the_ratio_MF_Be9_pid_acc_kin_full);
	TH1F* H1_thp_ratio_MF_Be9_pid_acc_kin_full = 0;
	file_MF_Be9->GetObject(base24.Data(), H1_thp_ratio_MF_Be9_pid_acc_kin_full);
	TH1F* H1_q_ratio_MF_Be9_pid_acc_kin_full = 0;
	file_MF_Be9->GetObject(base25.Data(), H1_q_ratio_MF_Be9_pid_acc_kin_full);
	TH1F* H1_nu_ratio_MF_Be9_pid_acc_kin_full = 0;
	file_MF_Be9->GetObject(base26.Data(), H1_nu_ratio_MF_Be9_pid_acc_kin_full);


	//MF B10
	TString MF_B10 = "Mf/pass4/B10/Coin/MF_B10_Results.root";
	TFile* inROOT_MF_B10 = new TFile(MF_B10.Data(), "READ");
	TTree* inputtree_MF_B10 = (TTree*)inROOT_MF_B10->Get("T");
	TFile* file_MF_B10 = NULL;
	file_MF_B10 = new TFile(MF_B10.Data());
	TH1F* H1_pmiss_ratio_MF_B10 = 0;
	TH1F* H1_pmiss_ratio_MF_B10_pid = 0;
	TH1F* H1_pmiss_ratio_MF_B10_pid_acc = 0;
	TH1F* H1_pmiss_ratio_MF_B10_pid_acc_kin = 0;
	file_MF_B10->GetObject(base.Data(), H1_pmiss_ratio_MF_B10);
	file_MF_B10->GetObject(base1.Data(), H1_pmiss_ratio_MF_B10_pid);
	file_MF_B10->GetObject(base2.Data(), H1_pmiss_ratio_MF_B10_pid_acc);
	file_MF_B10->GetObject(base3.Data(), H1_pmiss_ratio_MF_B10_pid_acc_kin);
	TH1F* H1_pmiss_ratio_MF_B10_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_MF_B10_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_MF_B10_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_B10_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf = 0;
	TH1F* H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf = 0;
	TH1F* H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf = 0;
	TH1F* H1_W_ratio_MF_B10_pid_acc_kin_full = 0;
	file_MF_B10->GetObject(base12.Data(), H1_epctime_ratio_MF_B10_pid_acc_kin_full);
	file_MF_B10->GetObject(base21.Data(), H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf);
	file_MF_B10->GetObject(base131.Data(), H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf);
	file_MF_B10->GetObject(base141.Data(), H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf);
	file_MF_B10->GetObject(base151.Data(), H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf);
	file_MF_B10->GetObject(base9.Data(), H1_pmiss_ratio_MF_B10_pid_acc_kin_full);
	file_MF_B10->GetObject(base10.Data(), H1_Q2_ratio_MF_B10_pid_acc_kin_full);
	file_MF_B10->GetObject(base11.Data(), H1_Em_ratio_MF_B10_pid_acc_kin_full);
	file_MF_B10->GetObject(base22.Data(), H1_W_ratio_MF_B10_pid_acc_kin_full);
	
	TH1F* H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf = 0;
	file_MF_B10->GetObject(base212.Data(), H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf);
	TH1F* H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf = 0;
	file_MF_B10->GetObject(base222.Data(), H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf);
	
	TH1F* H1_the_ratio_MF_B10_pid_acc_kin_full = 0;
	file_MF_B10->GetObject(base23.Data(), H1_the_ratio_MF_B10_pid_acc_kin_full);
	TH1F* H1_thp_ratio_MF_B10_pid_acc_kin_full = 0;
	file_MF_B10->GetObject(base24.Data(), H1_thp_ratio_MF_B10_pid_acc_kin_full);
	TH1F* H1_q_ratio_MF_B10_pid_acc_kin_full = 0;
	file_MF_B10->GetObject(base25.Data(), H1_q_ratio_MF_B10_pid_acc_kin_full);
	TH1F* H1_nu_ratio_MF_B10_pid_acc_kin_full = 0;
	file_MF_B10->GetObject(base26.Data(), H1_nu_ratio_MF_B10_pid_acc_kin_full);


	//MF B11
	TString MF_B11 = "Mf/pass4/B11/Coin/MF_B11_Results.root";
	TFile* inROOT_MF_B11 = new TFile(MF_B11.Data(), "READ");//.Data()
	TTree* inputtree_MF_B11 = (TTree*)inROOT_MF_B11->Get("T");
	TFile* file_MF_B11 = NULL;
	file_MF_B11 = new TFile(MF_B11.Data());
	TH1F* H1_pmiss_ratio_MF_B11 = 0;
	TH1F* H1_pmiss_ratio_MF_B11_pid = 0;
	TH1F* H1_pmiss_ratio_MF_B11_pid_acc = 0;
	TH1F* H1_pmiss_ratio_MF_B11_pid_acc_kin = 0;
	file_MF_B11->GetObject(base.Data(), H1_pmiss_ratio_MF_B11);
	file_MF_B11->GetObject(base1.Data(), H1_pmiss_ratio_MF_B11_pid);
	file_MF_B11->GetObject(base2.Data(), H1_pmiss_ratio_MF_B11_pid_acc);
	file_MF_B11->GetObject(base3.Data(), H1_pmiss_ratio_MF_B11_pid_acc_kin);
	TH1F* H1_pmiss_ratio_MF_B11_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_MF_B11_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_MF_B11_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_B11_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf = 0;
	TH1F* H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf = 0;
	TH1F* H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf = 0;
	TH1F* H1_W_ratio_MF_B11_pid_acc_kin_full = 0;
	file_MF_B11->GetObject(base12.Data(), H1_epctime_ratio_MF_B11_pid_acc_kin_full);
	file_MF_B11->GetObject(base21.Data(), H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf);
	file_MF_B11->GetObject(base131.Data(), H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf);
	file_MF_B11->GetObject(base141.Data(), H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf);
	file_MF_B11->GetObject(base151.Data(), H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf);
	file_MF_B11->GetObject(base9.Data(), H1_pmiss_ratio_MF_B11_pid_acc_kin_full);
	file_MF_B11->GetObject(base10.Data(), H1_Q2_ratio_MF_B11_pid_acc_kin_full);
	file_MF_B11->GetObject(base11.Data(), H1_Em_ratio_MF_B11_pid_acc_kin_full);
	file_MF_B11->GetObject(base22.Data(), H1_W_ratio_MF_B11_pid_acc_kin_full);
	
	TH1F* H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf = 0;
	file_MF_B11->GetObject(base212.Data(), H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf);
	TH1F* H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf = 0;
	file_MF_B11->GetObject(base222.Data(), H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf);
	
	TH1F* H1_the_ratio_MF_B11_pid_acc_kin_full = 0;
	file_MF_B11->GetObject(base23.Data(), H1_the_ratio_MF_B11_pid_acc_kin_full);
	TH1F* H1_thp_ratio_MF_B11_pid_acc_kin_full = 0;
	file_MF_B11->GetObject(base24.Data(), H1_thp_ratio_MF_B11_pid_acc_kin_full);
	TH1F* H1_q_ratio_MF_B11_pid_acc_kin_full = 0;
	file_MF_B11->GetObject(base25.Data(), H1_q_ratio_MF_B11_pid_acc_kin_full);
	TH1F* H1_nu_ratio_MF_B11_pid_acc_kin_full = 0;
	file_MF_B11->GetObject(base26.Data(), H1_nu_ratio_MF_B11_pid_acc_kin_full);


	//MF C12
	TString MF_C12 = "Mf/pass4/C12/Coin/MF_C12_Results.root";
	TFile* inROOT_MF_C12 = new TFile(MF_C12.Data(), "READ");//"READ"
	TTree* inputtree_MF_C12 = (TTree*)inROOT_MF_C12->Get("T");
	TFile* file_MF_C12 = NULL;
	file_MF_C12 = new TFile(MF_C12.Data()); 
	TH1F* H1_pmiss_ratio_MF_C12 = 0;
	TH1F* H1_pmiss_ratio_MF_C12_pid = 0;
	TH1F* H1_pmiss_ratio_MF_C12_pid_acc = 0;
	TH1F* H1_pmiss_ratio_MF_C12_pid_acc_kin = 0;
	file_MF_C12->GetObject(base.Data(), H1_pmiss_ratio_MF_C12);
	file_MF_C12->GetObject(base1.Data(), H1_pmiss_ratio_MF_C12_pid);
	file_MF_C12->GetObject(base2.Data(), H1_pmiss_ratio_MF_C12_pid_acc);
	file_MF_C12->GetObject(base3.Data(), H1_pmiss_ratio_MF_C12_pid_acc_kin);
	TH1F* H1_pmiss_ratio_MF_C12_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_MF_C12_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_MF_C12_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_C12_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	TH1F* H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	TH1F* H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	TH2F* H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	TH2F* H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	//TH2F* H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full = 0;
	//TH2F* H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full = 0;
	TH1F* H1_W_ratio_MF_C12_pid_acc_kin_full = 0;
	file_MF_C12->GetObject(base12.Data(), H1_epctime_ratio_MF_C12_pid_acc_kin_full);
	file_MF_C12->GetObject(base21.Data(), H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf);
	file_MF_C12->GetObject(base131.Data(), H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf);
	file_MF_C12->GetObject(base141.Data(), H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf);
	file_MF_C12->GetObject(base151.Data(), H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf);
	file_MF_C12->GetObject(base181.Data(), H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf);
	file_MF_C12->GetObject(base191.Data(), H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf);
	//file_MF_C12->GetObject(base16.Data(), H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full);
	//file_MF_C12->GetObject(base17.Data(), H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full);
	file_MF_C12->GetObject(base9.Data(), H1_pmiss_ratio_MF_C12_pid_acc_kin_full);
	file_MF_C12->GetObject(base10.Data(), H1_Q2_ratio_MF_C12_pid_acc_kin_full);
	file_MF_C12->GetObject(base11.Data(), H1_Em_ratio_MF_C12_pid_acc_kin_full);
	file_MF_C12->GetObject(base22.Data(), H1_W_ratio_MF_C12_pid_acc_kin_full);
	
	TH1F* H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	file_MF_C12->GetObject(base212.Data(), H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf);
	TH1F* H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	file_MF_C12->GetObject(base222.Data(), H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf);
	
	TH1F* H1_the_ratio_MF_C12_pid_acc_kin_full = 0;
	file_MF_C12->GetObject(base23.Data(), H1_the_ratio_MF_C12_pid_acc_kin_full);
	TH1F* H1_thp_ratio_MF_C12_pid_acc_kin_full = 0;
	file_MF_C12->GetObject(base24.Data(), H1_thp_ratio_MF_C12_pid_acc_kin_full);
	TH1F* H1_q_ratio_MF_C12_pid_acc_kin_full = 0;
	file_MF_C12->GetObject(base25.Data(), H1_q_ratio_MF_C12_pid_acc_kin_full);
	TH1F* H1_nu_ratio_MF_C12_pid_acc_kin_full = 0;
	file_MF_C12->GetObject(base26.Data(), H1_nu_ratio_MF_C12_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full = 0;
	file_MF_C12->GetObject(base192.Data(), H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full2 = 0;
	file_MF_C12->GetObject(base193.Data(), H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full2);
	TH2F* H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	file_MF_C12->GetObject(base194.Data(), H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf);
	TH2F* H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf2 = 0;
	file_MF_C12->GetObject(base195.Data(), H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf2);


	//MF Ca40
	TString MF_Ca40 = "Mf/pass4/Ca40/Coin/MF_Ca40_Results.root";
	TFile* inROOT_MF_Ca40 = new TFile(MF_Ca40.Data(), "READ");
	TTree* inputtree_MF_Ca40 = (TTree*)inROOT_MF_Ca40->Get("T");
	TFile* file_MF_Ca40 = NULL;
	file_MF_Ca40 = new TFile(MF_Ca40.Data());
	TH1F* H1_pmiss_ratio_MF_Ca40 = 0;
	TH1F* H1_pmiss_ratio_MF_Ca40_pid = 0;
	TH1F* H1_pmiss_ratio_MF_Ca40_pid_acc = 0;
	TH1F* H1_pmiss_ratio_MF_Ca40_pid_acc_kin = 0;
	file_MF_Ca40->GetObject(base.Data(), H1_pmiss_ratio_MF_Ca40);
	file_MF_Ca40->GetObject(base1.Data(), H1_pmiss_ratio_MF_Ca40_pid);
	file_MF_Ca40->GetObject(base2.Data(), H1_pmiss_ratio_MF_Ca40_pid_acc);
	file_MF_Ca40->GetObject(base3.Data(), H1_pmiss_ratio_MF_Ca40_pid_acc_kin);
	TH1F* H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_MF_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_MF_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf = 0;
	TH1F* H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf = 0;
	TH1F* H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf = 0;
	TH1F* H1_W_ratio_MF_Ca40_pid_acc_kin_full = 0;
	file_MF_Ca40->GetObject(base12.Data(), H1_epctime_ratio_MF_Ca40_pid_acc_kin_full);
	file_MF_Ca40->GetObject(base21.Data(), H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf);
	file_MF_Ca40->GetObject(base131.Data(), H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf);
	file_MF_Ca40->GetObject(base141.Data(), H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf);
	file_MF_Ca40->GetObject(base151.Data(), H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf);
	file_MF_Ca40->GetObject(base9.Data(), H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full);
	file_MF_Ca40->GetObject(base10.Data(), H1_Q2_ratio_MF_Ca40_pid_acc_kin_full);
	file_MF_Ca40->GetObject(base11.Data(), H1_Em_ratio_MF_Ca40_pid_acc_kin_full);
	file_MF_Ca40->GetObject(base22.Data(), H1_W_ratio_MF_Ca40_pid_acc_kin_full);
	
	TH1F* H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf = 0;
	file_MF_Ca40->GetObject(base212.Data(), H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf);
	TH1F* H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf = 0;
	file_MF_Ca40->GetObject(base222.Data(), H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf);
	
	TH1F* H1_the_ratio_MF_Ca40_pid_acc_kin_full = 0;
	file_MF_Ca40->GetObject(base23.Data(), H1_the_ratio_MF_Ca40_pid_acc_kin_full);
	TH1F* H1_thp_ratio_MF_Ca40_pid_acc_kin_full = 0;
	file_MF_Ca40->GetObject(base24.Data(), H1_thp_ratio_MF_Ca40_pid_acc_kin_full);
	TH1F* H1_q_ratio_MF_Ca40_pid_acc_kin_full = 0;
	file_MF_Ca40->GetObject(base25.Data(), H1_q_ratio_MF_Ca40_pid_acc_kin_full);
	TH1F* H1_nu_ratio_MF_Ca40_pid_acc_kin_full = 0;
	file_MF_Ca40->GetObject(base26.Data(), H1_nu_ratio_MF_Ca40_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full = 0;
	file_MF_Ca40->GetObject(base192.Data(), H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full2 = 0;
	file_MF_Ca40->GetObject(base193.Data(), H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full2);
	TH2F* H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full_mf = 0;
	file_MF_Ca40->GetObject(base194.Data(), H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full_mf);
	TH2F* H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full_mf2 = 0;
	file_MF_Ca40->GetObject(base195.Data(), H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full_mf2);


	//MF Ca48
	TString MF_Ca48 = "Mf/pass4/Ca48/Coin/MF_Ca48_Results.root";
	TFile* inROOT_MF_Ca48 = new TFile(MF_Ca48.Data(), "READ");
	TTree* inputtree_MF_Ca48 = (TTree*)inROOT_MF_Ca48->Get("T");
	TFile* file_MF_Ca48 = NULL;
	file_MF_Ca48 = new TFile(MF_Ca48.Data());
	TH1F* H1_pmiss_ratio_MF_Ca48 = 0;
	TH1F* H1_pmiss_ratio_MF_Ca48_pid = 0;
	TH1F* H1_pmiss_ratio_MF_Ca48_pid_acc = 0;
	TH1F* H1_pmiss_ratio_MF_Ca48_pid_acc_kin = 0;
	file_MF_Ca48->GetObject(base.Data(), H1_pmiss_ratio_MF_Ca48);
	file_MF_Ca48->GetObject(base1.Data(), H1_pmiss_ratio_MF_Ca48_pid);
	file_MF_Ca48->GetObject(base2.Data(), H1_pmiss_ratio_MF_Ca48_pid_acc);
	file_MF_Ca48->GetObject(base3.Data(), H1_pmiss_ratio_MF_Ca48_pid_acc_kin);
	TH1F* H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_MF_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_MF_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf = 0;
	TH1F* H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf = 0;
	TH1F* H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf = 0;
	TH1F* H1_W_ratio_MF_Ca48_pid_acc_kin_full = 0;
	file_MF_Ca48->GetObject(base12.Data(), H1_epctime_ratio_MF_Ca48_pid_acc_kin_full);
	file_MF_Ca48->GetObject(base21.Data(), H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf);
	file_MF_Ca48->GetObject(base131.Data(), H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf);
	file_MF_Ca48->GetObject(base141.Data(), H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf);
	file_MF_Ca48->GetObject(base151.Data(), H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf);
	file_MF_Ca48->GetObject(base9.Data(), H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full);
	file_MF_Ca48->GetObject(base10.Data(), H1_Q2_ratio_MF_Ca48_pid_acc_kin_full);
	file_MF_Ca48->GetObject(base11.Data(), H1_Em_ratio_MF_Ca48_pid_acc_kin_full);
	file_MF_Ca48->GetObject(base22.Data(), H1_W_ratio_MF_Ca48_pid_acc_kin_full);
	
	TH1F* H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf = 0;
	file_MF_Ca48->GetObject(base212.Data(), H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf);
	TH1F* H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf = 0;
	file_MF_Ca48->GetObject(base222.Data(), H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf);
	
	TH1F* H1_the_ratio_MF_Ca48_pid_acc_kin_full = 0;
	file_MF_Ca48->GetObject(base23.Data(), H1_the_ratio_MF_Ca48_pid_acc_kin_full);
	TH1F* H1_thp_ratio_MF_Ca48_pid_acc_kin_full = 0;
	file_MF_Ca48->GetObject(base24.Data(), H1_thp_ratio_MF_Ca48_pid_acc_kin_full);
	TH1F* H1_q_ratio_MF_Ca48_pid_acc_kin_full = 0;
	file_MF_Ca48->GetObject(base25.Data(), H1_q_ratio_MF_Ca48_pid_acc_kin_full);
	TH1F* H1_nu_ratio_MF_Ca48_pid_acc_kin_full = 0;
	file_MF_Ca48->GetObject(base26.Data(), H1_nu_ratio_MF_Ca48_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full = 0;
	file_MF_Ca48->GetObject(base192.Data(), H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full2 = 0;
	file_MF_Ca48->GetObject(base193.Data(), H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full2);
	TH2F* H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full_mf = 0;
	file_MF_Ca48->GetObject(base194.Data(), H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full_mf);
	TH2F* H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full_mf2 = 0;
	file_MF_Ca48->GetObject(base195.Data(), H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full_mf2);


	//MF Fe54
	TString MF_Fe54 = "Mf/pass4/Fe54/Coin/MF_Fe54_Results.root";
	TFile* inROOT_MF_Fe54 = new TFile(MF_Fe54.Data(), "READ");
	TTree* inputtree_MF_Fe54 = (TTree*)inROOT_MF_Fe54->Get("T");
	TFile* file_MF_Fe54 = NULL;
	file_MF_Fe54 = new TFile(MF_Fe54.Data());
	TH1F* H1_pmiss_ratio_MF_Fe54 = 0;
	TH1F* H1_pmiss_ratio_MF_Fe54_pid = 0;
	TH1F* H1_pmiss_ratio_MF_Fe54_pid_acc = 0;
	TH1F* H1_pmiss_ratio_MF_Fe54_pid_acc_kin = 0;
	file_MF_Fe54->GetObject(base.Data(), H1_pmiss_ratio_MF_Fe54);
	file_MF_Fe54->GetObject(base1.Data(), H1_pmiss_ratio_MF_Fe54_pid);
	file_MF_Fe54->GetObject(base2.Data(), H1_pmiss_ratio_MF_Fe54_pid_acc);
	file_MF_Fe54->GetObject(base3.Data(), H1_pmiss_ratio_MF_Fe54_pid_acc_kin);
	TH1F* H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_MF_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_MF_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Fe54_pid_acc_kin_full_mf = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_MF_Fe54_pid_acc_kin_full_mf = 0;
	TH1F* H1_hdelta_ratio_MF_Fe54_pid_acc_kin_full_mf = 0;
	TH1F* H1_edelta_ratio_MF_Fe54_pid_acc_kin_full_mf = 0;
	TH1F* H1_W_ratio_MF_Fe54_pid_acc_kin_full = 0;
	file_MF_Fe54->GetObject(base12.Data(), H1_epctime_ratio_MF_Fe54_pid_acc_kin_full);
	file_MF_Fe54->GetObject(base21.Data(), H1_epctime_ratio_MF_Fe54_pid_acc_kin_full_mf);
	file_MF_Fe54->GetObject(base131.Data(), H1_pCalEtotTrkNorm_ratio_MF_Fe54_pid_acc_kin_full_mf);
	file_MF_Fe54->GetObject(base141.Data(), H1_hdelta_ratio_MF_Fe54_pid_acc_kin_full_mf);
	file_MF_Fe54->GetObject(base151.Data(), H1_edelta_ratio_MF_Fe54_pid_acc_kin_full_mf);
	file_MF_Fe54->GetObject(base9.Data(), H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full);
	file_MF_Fe54->GetObject(base10.Data(), H1_Q2_ratio_MF_Fe54_pid_acc_kin_full);
	file_MF_Fe54->GetObject(base11.Data(), H1_Em_ratio_MF_Fe54_pid_acc_kin_full);
	file_MF_Fe54->GetObject(base22.Data(), H1_W_ratio_MF_Fe54_pid_acc_kin_full);
	
	TH1F* H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf = 0;
	file_MF_Fe54->GetObject(base212.Data(), H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf);
	TH1F* H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf = 0;
	file_MF_Fe54->GetObject(base222.Data(), H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf);
	
	TH1F* H1_the_ratio_MF_Fe54_pid_acc_kin_full = 0;
	file_MF_Fe54->GetObject(base23.Data(), H1_the_ratio_MF_Fe54_pid_acc_kin_full);
	TH1F* H1_thp_ratio_MF_Fe54_pid_acc_kin_full = 0;
	file_MF_Fe54->GetObject(base24.Data(), H1_thp_ratio_MF_Fe54_pid_acc_kin_full);
	TH1F* H1_q_ratio_MF_Fe54_pid_acc_kin_full = 0;
	file_MF_Fe54->GetObject(base25.Data(), H1_q_ratio_MF_Fe54_pid_acc_kin_full);
	TH1F* H1_nu_ratio_MF_Fe54_pid_acc_kin_full = 0;
	file_MF_Fe54->GetObject(base26.Data(), H1_nu_ratio_MF_Fe54_pid_acc_kin_full);


	//MF Au197
	TString MF_Au197 = "Mf/pass4/Au197/Coin/MF_Au197_Results.root";
	TFile* inROOT_MF_Au197 = new TFile(MF_Au197.Data(), "READ");
	TTree* inputtree_MF_Au197 = (TTree*)inROOT_MF_Au197->Get("T");
	TFile* file_MF_Au197 = NULL;
	file_MF_Au197 = new TFile(MF_Au197.Data());
	TH1F* H1_pmiss_ratio_MF_Au197 = 0;
	TH1F* H1_pmiss_ratio_MF_Au197_pid = 0;
	TH1F* H1_pmiss_ratio_MF_Au197_pid_acc = 0;
	TH1F* H1_pmiss_ratio_MF_Au197_pid_acc_kin = 0;
	file_MF_Au197->GetObject(base.Data(), H1_pmiss_ratio_MF_Au197);
	file_MF_Au197->GetObject(base1.Data(), H1_pmiss_ratio_MF_Au197_pid);
	file_MF_Au197->GetObject(base2.Data(), H1_pmiss_ratio_MF_Au197_pid_acc);
	file_MF_Au197->GetObject(base3.Data(), H1_pmiss_ratio_MF_Au197_pid_acc_kin);
	TH1F* H1_pmiss_ratio_MF_Au197_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_MF_Au197_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_MF_Au197_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Au197_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_MF_Au197_pid_acc_kin_full_mf = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_MF_Au197_pid_acc_kin_full_mf = 0;
	TH1F* H1_hdelta_ratio_MF_Au197_pid_acc_kin_full_mf = 0;
	TH1F* H1_edelta_ratio_MF_Au197_pid_acc_kin_full_mf = 0;
	TH1F* H1_W_ratio_MF_Au197_pid_acc_kin_full = 0;
	file_MF_Au197->GetObject(base12.Data(), H1_epctime_ratio_MF_Au197_pid_acc_kin_full);
	file_MF_Au197->GetObject(base21.Data(), H1_epctime_ratio_MF_Au197_pid_acc_kin_full_mf);
	file_MF_Au197->GetObject(base131.Data(), H1_pCalEtotTrkNorm_ratio_MF_Au197_pid_acc_kin_full_mf);
	file_MF_Au197->GetObject(base141.Data(), H1_hdelta_ratio_MF_Au197_pid_acc_kin_full_mf);
	file_MF_Au197->GetObject(base151.Data(), H1_edelta_ratio_MF_Au197_pid_acc_kin_full_mf);
	file_MF_Au197->GetObject(base9.Data(), H1_pmiss_ratio_MF_Au197_pid_acc_kin_full);
	file_MF_Au197->GetObject(base10.Data(), H1_Q2_ratio_MF_Au197_pid_acc_kin_full);
	file_MF_Au197->GetObject(base11.Data(), H1_Em_ratio_MF_Au197_pid_acc_kin_full);
	file_MF_Au197->GetObject(base22.Data(), H1_W_ratio_MF_Au197_pid_acc_kin_full);
	
	TH1F* H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf = 0;
	file_MF_Au197->GetObject(base212.Data(), H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf);
	TH1F* H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf = 0;
	file_MF_Au197->GetObject(base222.Data(), H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf);
	
	TH1F* H1_the_ratio_MF_Au197_pid_acc_kin_full = 0;
	file_MF_Au197->GetObject(base23.Data(), H1_the_ratio_MF_Au197_pid_acc_kin_full);
	TH1F* H1_thp_ratio_MF_Au197_pid_acc_kin_full = 0;
	file_MF_Au197->GetObject(base24.Data(), H1_thp_ratio_MF_Au197_pid_acc_kin_full);
	TH1F* H1_q_ratio_MF_Au197_pid_acc_kin_full = 0;
	file_MF_Au197->GetObject(base25.Data(), H1_q_ratio_MF_Au197_pid_acc_kin_full);
	TH1F* H1_nu_ratio_MF_Au197_pid_acc_kin_full = 0;
	file_MF_Au197->GetObject(base26.Data(), H1_nu_ratio_MF_Au197_pid_acc_kin_full);





	//SRC D2
	TString SRC_D2 = "Src/pass4/D2/Coin/SRC_D2_Results.root";
	TFile* inROOT_SRC_D2 = new TFile(SRC_D2.Data(), "READ");
	TTree* inputtree_SRC_D2 = (TTree*)inROOT_SRC_D2->Get("T");
	TFile* file_SRC_D2 = NULL;
	file_SRC_D2 = new TFile(SRC_D2.Data());
	TH1F* H1_pmiss_ratio_SRC_D2 = 0;
	TH1F* H1_pmiss_ratio_SRC_D2_pid = 0;
	TH1F* H1_pmiss_ratio_SRC_D2_pid_acc = 0;
	TH1F* H1_pmiss_ratio_SRC_D2_pid_acc_kin = 0;
	file_SRC_D2->GetObject(base.Data(), H1_pmiss_ratio_SRC_D2);
	file_SRC_D2->GetObject(base1.Data(), H1_pmiss_ratio_SRC_D2_pid);
	file_SRC_D2->GetObject(base2.Data(), H1_pmiss_ratio_SRC_D2_pid_acc);
	file_SRC_D2->GetObject(base3.Data(), H1_pmiss_ratio_SRC_D2_pid_acc_kin);
	TH1F* H1_W_ratio_SRC_D2_pid_acc_kin_full = 0;
	file_SRC_D2->GetObject(base8.Data(), H1_W_ratio_SRC_D2_pid_acc_kin_full);
	TH1F* H1_xbj_ratio_SRC_D2_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_SRC_D2_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_SRC_D2_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_SRC_D2_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_SRC_D2_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_SRC_D2_pid_acc_kin_full = 0;
	//TH1F* H1_epctime_ratio_SRC_D2_pid_acc_kin_full_alt = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_SRC_D2_pid_acc_kin_full = 0;
	TH1F* H1_hdelta_ratio_SRC_D2_pid_acc_kin_full = 0;
	TH1F* H1_edelta_ratio_SRC_D2_pid_acc_kin_full = 0;
	file_SRC_D2->GetObject(base12.Data(), H1_epctime_ratio_SRC_D2_pid_acc_kin_full);
	//file_SRC_D2->GetObject(base21.Data(), H1_epctime_ratio_SRC_D2_pid_acc_kin_full_alt);
	file_SRC_D2->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_D2_pid_acc_kin_full);
	file_SRC_D2->GetObject(base14.Data(), H1_hdelta_ratio_SRC_D2_pid_acc_kin_full);
	file_SRC_D2->GetObject(base15.Data(), H1_edelta_ratio_SRC_D2_pid_acc_kin_full);
	file_SRC_D2->GetObject(base4.Data(), H1_xbj_ratio_SRC_D2_pid_acc_kin_full);
	file_SRC_D2->GetObject(base5.Data(), H1_Q2_ratio_SRC_D2_pid_acc_kin_full);
	file_SRC_D2->GetObject(base6.Data(), H1_pmiss_ratio_SRC_D2_pid_acc_kin_full);
	file_SRC_D2->GetObject(base7.Data(), H1_thrq_ratio_SRC_D2_pid_acc_kin_full);
	file_SRC_D2->GetObject(base20.Data(), H1_Em_ratio_SRC_D2_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_SRC_D2_pid_acc_kin_full = 0;
	file_SRC_D2->GetObject(base211.Data(), H1_Ef_ratio_SRC_D2_pid_acc_kin_full);
	TH1F* H1_Pf_ratio_SRC_D2_pid_acc_kin_full = 0;
	file_SRC_D2->GetObject(base221.Data(), H1_Pf_ratio_SRC_D2_pid_acc_kin_full);
	TH1F* H1_the_ratio_SRC_D2_pid_acc_kin_full = 0;
	file_SRC_D2->GetObject(base23.Data(), H1_the_ratio_SRC_D2_pid_acc_kin_full);
	TH1F* H1_thp_ratio_SRC_D2_pid_acc_kin_full = 0;
	file_SRC_D2->GetObject(base24.Data(), H1_thp_ratio_SRC_D2_pid_acc_kin_full);
	TH1F* H1_q_ratio_SRC_D2_pid_acc_kin_full = 0;
	file_SRC_D2->GetObject(base25.Data(), H1_q_ratio_SRC_D2_pid_acc_kin_full);
	TH1F* H1_nu_ratio_SRC_D2_pid_acc_kin_full = 0;
	file_SRC_D2->GetObject(base26.Data(), H1_nu_ratio_SRC_D2_pid_acc_kin_full);


	//SRC Be9
	TString SRC_Be9 = "Src/pass4/Be9/Coin/SRC_Be9_Results.root";
	TFile* inROOT_SRC_Be9 = new TFile(SRC_Be9.Data(), "READ");
	TTree* inputtree_SRC_Be9 = (TTree*)inROOT_SRC_Be9->Get("T");
	TFile* file_SRC_Be9 = NULL;
	file_SRC_Be9 = new TFile(SRC_Be9.Data());
	TH1F* H1_pmiss_ratio_SRC_Be9 = 0;
	TH1F* H1_pmiss_ratio_SRC_Be9_pid = 0;
	TH1F* H1_pmiss_ratio_SRC_Be9_pid_acc = 0;
	TH1F* H1_pmiss_ratio_SRC_Be9_pid_acc_kin = 0;
	file_SRC_Be9->GetObject(base.Data(), H1_pmiss_ratio_SRC_Be9);
	file_SRC_Be9->GetObject(base1.Data(), H1_pmiss_ratio_SRC_Be9_pid);
	file_SRC_Be9->GetObject(base2.Data(), H1_pmiss_ratio_SRC_Be9_pid_acc);
	file_SRC_Be9->GetObject(base3.Data(), H1_pmiss_ratio_SRC_Be9_pid_acc_kin);
	TH1F* H1_xbj_ratio_SRC_Be9_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_SRC_Be9_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_SRC_Be9_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_SRC_Be9_pid_acc_kin_full = 0;
	TH1F* H1_W_ratio_SRC_Be9_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_SRC_Be9_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_SRC_Be9_pid_acc_kin_full = 0;
	//TH1F* H1_epctime_ratio_SRC_Be9_pid_acc_kin_full_alt = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_SRC_Be9_pid_acc_kin_full = 0;
	TH1F* H1_hdelta_ratio_SRC_Be9_pid_acc_kin_full = 0;
	TH1F* H1_edelta_ratio_SRC_Be9_pid_acc_kin_full = 0;
	file_SRC_Be9->GetObject(base12.Data(), H1_epctime_ratio_SRC_Be9_pid_acc_kin_full);
	//file_SRC_Be9->GetObject(base21.Data(), H1_epctime_ratio_SRC_Be9_pid_acc_kin_full_alt);
	file_SRC_Be9->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_Be9_pid_acc_kin_full);
	file_SRC_Be9->GetObject(base14.Data(), H1_hdelta_ratio_SRC_Be9_pid_acc_kin_full);
	file_SRC_Be9->GetObject(base15.Data(), H1_edelta_ratio_SRC_Be9_pid_acc_kin_full);
	file_SRC_Be9->GetObject(base4.Data(), H1_xbj_ratio_SRC_Be9_pid_acc_kin_full);
	file_SRC_Be9->GetObject(base5.Data(), H1_Q2_ratio_SRC_Be9_pid_acc_kin_full);
	file_SRC_Be9->GetObject(base6.Data(), H1_pmiss_ratio_SRC_Be9_pid_acc_kin_full);
	file_SRC_Be9->GetObject(base7.Data(), H1_thrq_ratio_SRC_Be9_pid_acc_kin_full);
	file_SRC_Be9->GetObject(base8.Data(), H1_W_ratio_SRC_Be9_pid_acc_kin_full);
	file_SRC_Be9->GetObject(base20.Data(), H1_Em_ratio_SRC_Be9_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_SRC_Be9_pid_acc_kin_full = 0;
	file_SRC_Be9->GetObject(base211.Data(), H1_Ef_ratio_SRC_Be9_pid_acc_kin_full);
	TH1F* H1_Pf_ratio_SRC_Be9_pid_acc_kin_full = 0;
	file_SRC_Be9->GetObject(base221.Data(), H1_Pf_ratio_SRC_Be9_pid_acc_kin_full);
	TH1F* H1_the_ratio_SRC_Be9_pid_acc_kin_full = 0;
	file_SRC_Be9->GetObject(base23.Data(), H1_the_ratio_SRC_Be9_pid_acc_kin_full);
	TH1F* H1_thp_ratio_SRC_Be9_pid_acc_kin_full = 0;
	file_SRC_Be9->GetObject(base24.Data(), H1_thp_ratio_SRC_Be9_pid_acc_kin_full);
	TH1F* H1_q_ratio_SRC_Be9_pid_acc_kin_full = 0;
	file_SRC_Be9->GetObject(base25.Data(), H1_q_ratio_SRC_Be9_pid_acc_kin_full);
	TH1F* H1_nu_ratio_SRC_Be9_pid_acc_kin_full = 0;
	file_SRC_Be9->GetObject(base26.Data(), H1_nu_ratio_SRC_Be9_pid_acc_kin_full);


	//SRC B10
	TString SRC_B10 = "Src/pass4/B10/Coin/SRC_B10_Results.root";
	TFile* inROOT_SRC_B10 = new TFile(SRC_B10.Data(), "READ");
	TTree* inputtree_SRC_B10 = (TTree*)inROOT_SRC_B10->Get("T");
	TFile* file_SRC_B10 = NULL;
	file_SRC_B10 = new TFile(SRC_B10.Data());
	TH1F* H1_pmiss_ratio_SRC_B10 = 0;
	TH1F* H1_pmiss_ratio_SRC_B10_pid = 0;
	TH1F* H1_pmiss_ratio_SRC_B10_pid_acc = 0;
	TH1F* H1_pmiss_ratio_SRC_B10_pid_acc_kin = 0;
	file_SRC_B10->GetObject(base.Data(), H1_pmiss_ratio_SRC_B10);
	file_SRC_B10->GetObject(base1.Data(), H1_pmiss_ratio_SRC_B10_pid);
	file_SRC_B10->GetObject(base2.Data(), H1_pmiss_ratio_SRC_B10_pid_acc);
	file_SRC_B10->GetObject(base3.Data(), H1_pmiss_ratio_SRC_B10_pid_acc_kin);
	TH1F* H1_xbj_ratio_SRC_B10_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_SRC_B10_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_SRC_B10_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_SRC_B10_pid_acc_kin_full = 0;
	TH1F* H1_W_ratio_SRC_B10_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_SRC_B10_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_SRC_B10_pid_acc_kin_full = 0;
	//TH1F* H1_epctime_ratio_SRC_B10_pid_acc_kin_full_alt = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full = 0;
	TH1F* H1_hdelta_ratio_SRC_B10_pid_acc_kin_full = 0;
	TH1F* H1_edelta_ratio_SRC_B10_pid_acc_kin_full = 0;
	file_SRC_B10->GetObject(base12.Data(), H1_epctime_ratio_SRC_B10_pid_acc_kin_full);
	//file_SRC_B10->GetObject(base21.Data(), H1_epctime_ratio_SRC_B10_pid_acc_kin_full_alt);
	file_SRC_B10->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full);
	file_SRC_B10->GetObject(base14.Data(), H1_hdelta_ratio_SRC_B10_pid_acc_kin_full);
	file_SRC_B10->GetObject(base15.Data(), H1_edelta_ratio_SRC_B10_pid_acc_kin_full);
	file_SRC_B10->GetObject(base4.Data(), H1_xbj_ratio_SRC_B10_pid_acc_kin_full);
	file_SRC_B10->GetObject(base5.Data(), H1_Q2_ratio_SRC_B10_pid_acc_kin_full);
	file_SRC_B10->GetObject(base6.Data(), H1_pmiss_ratio_SRC_B10_pid_acc_kin_full);
	file_SRC_B10->GetObject(base7.Data(), H1_thrq_ratio_SRC_B10_pid_acc_kin_full);
	file_SRC_B10->GetObject(base8.Data(), H1_W_ratio_SRC_B10_pid_acc_kin_full);
	file_SRC_B10->GetObject(base20.Data(), H1_Em_ratio_SRC_B10_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_SRC_B10_pid_acc_kin_full = 0;
	file_SRC_B10->GetObject(base211.Data(), H1_Ef_ratio_SRC_B10_pid_acc_kin_full);
	TH1F* H1_Pf_ratio_SRC_B10_pid_acc_kin_full = 0;
	file_SRC_B10->GetObject(base221.Data(), H1_Pf_ratio_SRC_B10_pid_acc_kin_full);
	TH1F* H1_the_ratio_SRC_B10_pid_acc_kin_full = 0;
	file_SRC_B10->GetObject(base23.Data(), H1_the_ratio_SRC_B10_pid_acc_kin_full);
	TH1F* H1_thp_ratio_SRC_B10_pid_acc_kin_full = 0;
	file_SRC_B10->GetObject(base24.Data(), H1_thp_ratio_SRC_B10_pid_acc_kin_full);
	TH1F* H1_q_ratio_SRC_B10_pid_acc_kin_full = 0;
	file_SRC_B10->GetObject(base25.Data(), H1_q_ratio_SRC_B10_pid_acc_kin_full);
	TH1F* H1_nu_ratio_SRC_B10_pid_acc_kin_full = 0;
	file_SRC_B10->GetObject(base26.Data(), H1_nu_ratio_SRC_B10_pid_acc_kin_full);


	//SRC B11
	TString SRC_B11 = "Src/pass4/B11/Coin/SRC_B11_Results.root";
	TFile* inROOT_SRC_B11 = new TFile(SRC_B11.Data(), "READ");
	TTree* inputtree_SRC_B11 = (TTree*)inROOT_SRC_B11->Get("T");
	TFile* file_SRC_B11 = NULL;
	file_SRC_B11 = new TFile(SRC_B11.Data());
	TH1F* H1_pmiss_ratio_SRC_B11 = 0;
	TH1F* H1_pmiss_ratio_SRC_B11_pid = 0;
	TH1F* H1_pmiss_ratio_SRC_B11_pid_acc = 0;
	TH1F* H1_pmiss_ratio_SRC_B11_pid_acc_kin = 0;
	file_SRC_B11->GetObject(base.Data(), H1_pmiss_ratio_SRC_B11);
	file_SRC_B11->GetObject(base1.Data(), H1_pmiss_ratio_SRC_B11_pid);
	file_SRC_B11->GetObject(base2.Data(), H1_pmiss_ratio_SRC_B11_pid_acc);
	file_SRC_B11->GetObject(base3.Data(), H1_pmiss_ratio_SRC_B11_pid_acc_kin);
	TH1F* H1_xbj_ratio_SRC_B11_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_SRC_B11_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_SRC_B11_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_SRC_B11_pid_acc_kin_full = 0;
	TH1F* H1_W_ratio_SRC_B11_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_SRC_B11_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_SRC_B11_pid_acc_kin_full = 0;
	//TH1F* H1_epctime_ratio_SRC_B11_pid_acc_kin_full_alt = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full = 0;
	TH1F* H1_hdelta_ratio_SRC_B11_pid_acc_kin_full = 0;
	TH1F* H1_edelta_ratio_SRC_B11_pid_acc_kin_full = 0;
	file_SRC_B11->GetObject(base12.Data(), H1_epctime_ratio_SRC_B11_pid_acc_kin_full);
	//file_SRC_B11->GetObject(base21.Data(), H1_epctime_ratio_SRC_B11_pid_acc_kin_full_alt);
	file_SRC_B11->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full);
	file_SRC_B11->GetObject(base14.Data(), H1_hdelta_ratio_SRC_B11_pid_acc_kin_full);
	file_SRC_B11->GetObject(base15.Data(), H1_edelta_ratio_SRC_B11_pid_acc_kin_full);
	file_SRC_B11->GetObject(base4.Data(), H1_xbj_ratio_SRC_B11_pid_acc_kin_full);
	file_SRC_B11->GetObject(base5.Data(), H1_Q2_ratio_SRC_B11_pid_acc_kin_full);
	file_SRC_B11->GetObject(base6.Data(), H1_pmiss_ratio_SRC_B11_pid_acc_kin_full);
	file_SRC_B11->GetObject(base7.Data(), H1_thrq_ratio_SRC_B11_pid_acc_kin_full);
	file_SRC_B11->GetObject(base8.Data(), H1_W_ratio_SRC_B11_pid_acc_kin_full);
	file_SRC_B11->GetObject(base20.Data(), H1_Em_ratio_SRC_B11_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_SRC_B11_pid_acc_kin_full = 0;
	file_SRC_B11->GetObject(base211.Data(), H1_Ef_ratio_SRC_B11_pid_acc_kin_full);
	TH1F* H1_Pf_ratio_SRC_B11_pid_acc_kin_full = 0;
	file_SRC_B11->GetObject(base221.Data(), H1_Pf_ratio_SRC_B11_pid_acc_kin_full);
	TH1F* H1_the_ratio_SRC_B11_pid_acc_kin_full = 0;
	file_SRC_B11->GetObject(base23.Data(), H1_the_ratio_SRC_B11_pid_acc_kin_full);
	TH1F* H1_thp_ratio_SRC_B11_pid_acc_kin_full = 0;
	file_SRC_B11->GetObject(base24.Data(), H1_thp_ratio_SRC_B11_pid_acc_kin_full);
	TH1F* H1_q_ratio_SRC_B11_pid_acc_kin_full = 0;
	file_SRC_B11->GetObject(base25.Data(), H1_q_ratio_SRC_B11_pid_acc_kin_full);
	TH1F* H1_nu_ratio_SRC_B11_pid_acc_kin_full = 0;
	file_SRC_B11->GetObject(base26.Data(), H1_nu_ratio_SRC_B11_pid_acc_kin_full);


	//SRC C12
	TString SRC_C12 = "Src/pass4/C12/Coin/SRC_C12_Results.root";
	TFile* inROOT_SRC_C12 = new TFile(SRC_C12.Data(), "READ");
	TTree* inputtree_SRC_C12 = (TTree*)inROOT_SRC_C12->Get("T");
	TFile* file_SRC_C12 = NULL;
	file_SRC_C12 = new TFile(SRC_C12.Data());
	TH1F* H1_pmiss_ratio_SRC_C12 = 0;
	TH1F* H1_pmiss_ratio_SRC_C12_pid = 0;
	TH1F* H1_pmiss_ratio_SRC_C12_pid_acc = 0;
	TH1F* H1_pmiss_ratio_SRC_C12_pid_acc_kin = 0;
	file_SRC_C12->GetObject(base.Data(), H1_pmiss_ratio_SRC_C12);
	file_SRC_C12->GetObject(base1.Data(), H1_pmiss_ratio_SRC_C12_pid);
	file_SRC_C12->GetObject(base2.Data(), H1_pmiss_ratio_SRC_C12_pid_acc);
	file_SRC_C12->GetObject(base3.Data(), H1_pmiss_ratio_SRC_C12_pid_acc_kin);
	TH1F* H1_xbj_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH1F* H1_W_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_SRC_C12_pid_acc_kin_full = 0;
	//TH1F* H1_epctime_ratio_SRC_C12_pid_acc_kin_full_alt = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH1F* H1_hdelta_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH1F* H1_edelta_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH2F* H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH2F* H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full = 0;
	TH2F* H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full = 0; 
	file_SRC_C12->GetObject(base12.Data(), H1_epctime_ratio_SRC_C12_pid_acc_kin_full);
	//file_SRC_C12->GetObject(base21.Data(), H1_epctime_ratio_SRC_C12_pid_acc_kin_full_alt);
	file_SRC_C12->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base14.Data(), H1_hdelta_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base15.Data(), H1_edelta_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base18.Data(), H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base19.Data(), H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base17.Data(), H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base4.Data(), H1_xbj_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base5.Data(), H1_Q2_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base6.Data(), H1_pmiss_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base7.Data(), H1_thrq_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base8.Data(), H1_W_ratio_SRC_C12_pid_acc_kin_full);
	file_SRC_C12->GetObject(base8.Data(), H1_Em_ratio_SRC_C12_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_SRC_C12_pid_acc_kin_full = 0;
	file_SRC_C12->GetObject(base211.Data(), H1_Ef_ratio_SRC_C12_pid_acc_kin_full);
	TH1F* H1_Pf_ratio_SRC_C12_pid_acc_kin_full = 0;
	file_SRC_C12->GetObject(base221.Data(), H1_Pf_ratio_SRC_C12_pid_acc_kin_full);
	TH1F* H1_the_ratio_SRC_C12_pid_acc_kin_full = 0;
	file_SRC_C12->GetObject(base23.Data(), H1_the_ratio_SRC_C12_pid_acc_kin_full);
	TH1F* H1_thp_ratio_SRC_C12_pid_acc_kin_full = 0;
	file_SRC_C12->GetObject(base24.Data(), H1_thp_ratio_SRC_C12_pid_acc_kin_full);
	TH1F* H1_q_ratio_SRC_C12_pid_acc_kin_full = 0;
	file_SRC_C12->GetObject(base25.Data(), H1_q_ratio_SRC_C12_pid_acc_kin_full);
	TH1F* H1_nu_ratio_SRC_C12_pid_acc_kin_full = 0;
	file_SRC_C12->GetObject(base26.Data(), H1_nu_ratio_SRC_C12_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full = 0;
	file_SRC_C12->GetObject(base192.Data(), H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full2 = 0;
	file_SRC_C12->GetObject(base193.Data(), H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full2);
	TH2F* H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full_mf = 0;
	file_SRC_C12->GetObject(base194.Data(), H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full_mf);
	TH2F* H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full_mf2 = 0;
	file_SRC_C12->GetObject(base195.Data(), H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full_mf2);


	//SRC Ca40
	TString SRC_Ca40 = "Src/pass4/Ca40/Coin/SRC_Ca40_Results.root";
	TFile* inROOT_SRC_Ca40 = new TFile(SRC_Ca40.Data(), "READ");
	TTree* inputtree_SRC_Ca40 = (TTree*)inROOT_SRC_Ca40->Get("T");
	TFile* file_SRC_Ca40 = NULL;
	file_SRC_Ca40 = new TFile(SRC_Ca40.Data());
	TH1F* H1_pmiss_ratio_SRC_Ca40 = 0;
	TH1F* H1_pmiss_ratio_SRC_Ca40_pid = 0;
	TH1F* H1_pmiss_ratio_SRC_Ca40_pid_acc = 0;
	TH1F* H1_pmiss_ratio_SRC_Ca40_pid_acc_kin = 0;
	file_SRC_Ca40->GetObject(base.Data(), H1_pmiss_ratio_SRC_Ca40);
	file_SRC_Ca40->GetObject(base1.Data(), H1_pmiss_ratio_SRC_Ca40_pid);
	file_SRC_Ca40->GetObject(base2.Data(), H1_pmiss_ratio_SRC_Ca40_pid_acc);
	file_SRC_Ca40->GetObject(base3.Data(), H1_pmiss_ratio_SRC_Ca40_pid_acc_kin);
	TH1F* H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_W_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	//TH1F* H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full_alt = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH1F* H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH2F* H2_xbj_Q2_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH2F* H2_hXColl_hYColl_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	TH2F* H2_eXColl_eYColl_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	file_SRC_Ca40->GetObject(base12.Data(), H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full);
	//file_SRC_Ca40->GetObject(base21.Data(), H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full_alt);
	file_SRC_Ca40->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base14.Data(), H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base15.Data(), H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base4.Data(), H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base5.Data(), H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base6.Data(), H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base7.Data(), H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base8.Data(), H1_W_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base8.Data(), H1_Em_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base17.Data(), H2_xbj_Q2_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base18.Data(), H2_hXColl_hYColl_ratio_SRC_Ca40_pid_acc_kin_full);
	file_SRC_Ca40->GetObject(base19.Data(), H2_eXColl_eYColl_ratio_SRC_Ca40_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	file_SRC_Ca40->GetObject(base211.Data(), H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full);
	TH1F* H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	file_SRC_Ca40->GetObject(base221.Data(), H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full);
	TH1F* H1_the_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	file_SRC_Ca40->GetObject(base23.Data(), H1_the_ratio_SRC_Ca40_pid_acc_kin_full);
	TH1F* H1_thp_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	file_SRC_Ca40->GetObject(base24.Data(), H1_thp_ratio_SRC_Ca40_pid_acc_kin_full);
	TH1F* H1_q_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	file_SRC_Ca40->GetObject(base25.Data(), H1_q_ratio_SRC_Ca40_pid_acc_kin_full);
	TH1F* H1_nu_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	file_SRC_Ca40->GetObject(base26.Data(), H1_nu_ratio_SRC_Ca40_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_SRC_Ca40_pid_acc_kin_full = 0;
	file_SRC_Ca40->GetObject(base192.Data(), H2_Em_Pm_ratio_SRC_Ca40_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_SRC_Ca40_pid_acc_kin_full2 = 0;
	file_SRC_Ca40->GetObject(base193.Data(), H2_Em_Pm_ratio_SRC_Ca40_pid_acc_kin_full2);
	TH2F* H2_Em_Pm_ratio_SRC_Ca40_pid_acc_kin_full_mf = 0;
	file_SRC_Ca40->GetObject(base194.Data(), H2_Em_Pm_ratio_SRC_Ca40_pid_acc_kin_full_mf);
	TH2F* H2_Em_Pm_ratio_SRC_Ca40_pid_acc_kin_full_mf2 = 0;
	file_SRC_Ca40->GetObject(base195.Data(), H2_Em_Pm_ratio_SRC_Ca40_pid_acc_kin_full_mf2);


	//SRC Ca48
	TString SRC_Ca48 = "Src/pass4/Ca48/Coin/SRC_Ca48_Results.root";
	TFile* inROOT_SRC_Ca48 = new TFile(SRC_Ca48.Data(), "READ");
	TTree* inputtree_SRC_Ca48 = (TTree*)inROOT_SRC_Ca48->Get("T");
	TFile* file_SRC_Ca48 = NULL;
	file_SRC_Ca48 = new TFile(SRC_Ca48.Data());
	TH1F* H1_pmiss_ratio_SRC_Ca48 = 0;
	TH1F* H1_pmiss_ratio_SRC_Ca48_pid = 0;
	TH1F* H1_pmiss_ratio_SRC_Ca48_pid_acc = 0;
	TH1F* H1_pmiss_ratio_SRC_Ca48_pid_acc_kin = 0;
	file_SRC_Ca48->GetObject(base.Data(), H1_pmiss_ratio_SRC_Ca48);
	file_SRC_Ca48->GetObject(base1.Data(), H1_pmiss_ratio_SRC_Ca48_pid);
	file_SRC_Ca48->GetObject(base2.Data(), H1_pmiss_ratio_SRC_Ca48_pid_acc);
	file_SRC_Ca48->GetObject(base3.Data(), H1_pmiss_ratio_SRC_Ca48_pid_acc_kin);
	TH1F* H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_W_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	//TH1F* H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full_alt = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	TH1F* H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	file_SRC_Ca48->GetObject(base12.Data(), H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full);
	//file_SRC_Ca48->GetObject(base21.Data(), H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full_alt);
	file_SRC_Ca48->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full);
	file_SRC_Ca48->GetObject(base14.Data(), H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full);
	file_SRC_Ca48->GetObject(base15.Data(), H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full);
	file_SRC_Ca48->GetObject(base4.Data(), H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full);
	file_SRC_Ca48->GetObject(base5.Data(), H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full);
	file_SRC_Ca48->GetObject(base6.Data(), H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full);
	file_SRC_Ca48->GetObject(base7.Data(), H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full);
	file_SRC_Ca48->GetObject(base8.Data(), H1_W_ratio_SRC_Ca48_pid_acc_kin_full);
	file_SRC_Ca48->GetObject(base8.Data(), H1_Em_ratio_SRC_Ca48_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	file_SRC_Ca48->GetObject(base211.Data(), H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full);
	TH1F* H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	file_SRC_Ca48->GetObject(base221.Data(), H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full);
	TH1F* H1_the_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	file_SRC_Ca48->GetObject(base23.Data(), H1_the_ratio_SRC_Ca48_pid_acc_kin_full);
	TH1F* H1_thp_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	file_SRC_Ca48->GetObject(base24.Data(), H1_thp_ratio_SRC_Ca48_pid_acc_kin_full);
	TH1F* H1_q_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	file_SRC_Ca48->GetObject(base25.Data(), H1_q_ratio_SRC_Ca48_pid_acc_kin_full);
	TH1F* H1_nu_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	file_SRC_Ca48->GetObject(base26.Data(), H1_nu_ratio_SRC_Ca48_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_SRC_Ca48_pid_acc_kin_full = 0;
	file_SRC_Ca48->GetObject(base192.Data(), H2_Em_Pm_ratio_SRC_Ca48_pid_acc_kin_full);
	TH2F* H2_Em_Pm_ratio_SRC_Ca48_pid_acc_kin_full2 = 0;
	file_SRC_Ca48->GetObject(base193.Data(), H2_Em_Pm_ratio_SRC_Ca48_pid_acc_kin_full2);
	TH2F* H2_Em_Pm_ratio_SRC_Ca48_pid_acc_kin_full_mf = 0;
	file_SRC_Ca48->GetObject(base194.Data(), H2_Em_Pm_ratio_SRC_Ca48_pid_acc_kin_full_mf);
	TH2F* H2_Em_Pm_ratio_SRC_Ca48_pid_acc_kin_full_mf2 = 0;
	file_SRC_Ca48->GetObject(base195.Data(), H2_Em_Pm_ratio_SRC_Ca48_pid_acc_kin_full_mf2);


	//SRC Fe54
	TString SRC_Fe54 = "Src/pass4/Fe54/Coin/SRC_Fe54_Results.root";
	TFile* inROOT_SRC_Fe54 = new TFile(SRC_Fe54.Data(), "READ");
	TTree* inputtree_SRC_Fe54 = (TTree*)inROOT_SRC_Fe54->Get("T");
	TFile* file_SRC_Fe54 = NULL;
	file_SRC_Fe54 = new TFile(SRC_Fe54.Data());
	TH1F* H1_pmiss_ratio_SRC_Fe54 = 0;
	TH1F* H1_pmiss_ratio_SRC_Fe54_pid = 0;
	TH1F* H1_pmiss_ratio_SRC_Fe54_pid_acc = 0;
	TH1F* H1_pmiss_ratio_SRC_Fe54_pid_acc_kin = 0;
	file_SRC_Fe54->GetObject(base.Data(), H1_pmiss_ratio_SRC_Fe54);
	file_SRC_Fe54->GetObject(base1.Data(), H1_pmiss_ratio_SRC_Fe54_pid);
	file_SRC_Fe54->GetObject(base2.Data(), H1_pmiss_ratio_SRC_Fe54_pid_acc);
	file_SRC_Fe54->GetObject(base3.Data(), H1_pmiss_ratio_SRC_Fe54_pid_acc_kin);
	TH1F* H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_W_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	//TH1F* H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full_alt = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	TH1F* H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	file_SRC_Fe54->GetObject(base12.Data(), H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full);
	//file_SRC_Fe54->GetObject(base21.Data(), H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full_alt);
	file_SRC_Fe54->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base14.Data(), H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base15.Data(), H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base4.Data(), H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base5.Data(), H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base6.Data(), H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base7.Data(), H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base8.Data(), H1_W_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base8.Data(), H1_Em_ratio_SRC_Fe54_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	file_SRC_Fe54->GetObject(base211.Data(), H1_Ef_ratio_SRC_Fe54_pid_acc_kin_full);
	TH1F* H1_Pf_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	file_SRC_Fe54->GetObject(base221.Data(), H1_Pf_ratio_SRC_Fe54_pid_acc_kin_full);
	TH1F* H1_the_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	file_SRC_Fe54->GetObject(base23.Data(), H1_the_ratio_SRC_Fe54_pid_acc_kin_full);
	TH1F* H1_thp_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	file_SRC_Fe54->GetObject(base24.Data(), H1_thp_ratio_SRC_Fe54_pid_acc_kin_full);
	TH1F* H1_q_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	file_SRC_Fe54->GetObject(base25.Data(), H1_q_ratio_SRC_Fe54_pid_acc_kin_full);
	TH1F* H1_nu_ratio_SRC_Fe54_pid_acc_kin_full = 0;
	file_SRC_Fe54->GetObject(base26.Data(), H1_nu_ratio_SRC_Fe54_pid_acc_kin_full);


	//SRC Au197
	TString SRC_Au197 = "Src/pass4/Au197/Coin/SRC_Au197_Results.root";
	TFile* inROOT_SRC_Au197 = new TFile(SRC_Au197.Data(), "READ");
	TTree* inputtree_SRC_Au197 = (TTree*)inROOT_SRC_Au197->Get("T");
	TFile* file_SRC_Au197 = NULL;
	file_SRC_Au197 = new TFile(SRC_Au197.Data());
	TH1F* H1_pmiss_ratio_SRC_Au197 = 0;
	TH1F* H1_pmiss_ratio_SRC_Au197_pid = 0;
	TH1F* H1_pmiss_ratio_SRC_Au197_pid_acc = 0;
	TH1F* H1_pmiss_ratio_SRC_Au197_pid_acc_kin = 0;
	file_SRC_Au197->GetObject(base.Data(), H1_pmiss_ratio_SRC_Au197);
	file_SRC_Au197->GetObject(base1.Data(), H1_pmiss_ratio_SRC_Au197_pid);
	file_SRC_Au197->GetObject(base2.Data(), H1_pmiss_ratio_SRC_Au197_pid_acc);
	file_SRC_Au197->GetObject(base3.Data(), H1_pmiss_ratio_SRC_Au197_pid_acc_kin);
	TH1F* H1_xbj_ratio_SRC_Au197_pid_acc_kin_full = 0;
	TH1F* H1_Q2_ratio_SRC_Au197_pid_acc_kin_full = 0;
	TH1F* H1_pmiss_ratio_SRC_Au197_pid_acc_kin_full = 0;
	TH1F* H1_thrq_ratio_SRC_Au197_pid_acc_kin_full = 0;
	TH1F* H1_W_ratio_SRC_Au197_pid_acc_kin_full = 0;
	TH1F* H1_Em_ratio_SRC_Au197_pid_acc_kin_full = 0;
	TH1F* H1_epctime_ratio_SRC_Au197_pid_acc_kin_full = 0;
	//TH1F* H1_epctime_ratio_SRC_Au197_pid_acc_kin_full_alt = 0;
	TH1F* H1_pCalEtotTrkNorm_ratio_SRC_Au197_pid_acc_kin_full = 0;
	TH1F* H1_hdelta_ratio_SRC_Au197_pid_acc_kin_full = 0;
	TH1F* H1_edelta_ratio_SRC_Au197_pid_acc_kin_full = 0;
	file_SRC_Au197->GetObject(base12.Data(), H1_epctime_ratio_SRC_Au197_pid_acc_kin_full);
	//file_SRC_Au197->GetObject(base21.Data(), H1_epctime_ratio_SRC_Au197_pid_acc_kin_full_alt);
	file_SRC_Au197->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_Au197_pid_acc_kin_full);
	file_SRC_Au197->GetObject(base14.Data(), H1_hdelta_ratio_SRC_Au197_pid_acc_kin_full);
	file_SRC_Au197->GetObject(base15.Data(), H1_edelta_ratio_SRC_Au197_pid_acc_kin_full);
	file_SRC_Au197->GetObject(base4.Data(), H1_xbj_ratio_SRC_Au197_pid_acc_kin_full);
	file_SRC_Au197->GetObject(base5.Data(), H1_Q2_ratio_SRC_Au197_pid_acc_kin_full);
	file_SRC_Au197->GetObject(base6.Data(), H1_pmiss_ratio_SRC_Au197_pid_acc_kin_full);
	file_SRC_Au197->GetObject(base7.Data(), H1_thrq_ratio_SRC_Au197_pid_acc_kin_full);
	file_SRC_Au197->GetObject(base8.Data(), H1_W_ratio_SRC_Au197_pid_acc_kin_full);
	file_SRC_Au197->GetObject(base8.Data(), H1_Em_ratio_SRC_Au197_pid_acc_kin_full);
	TH1F* H1_Ef_ratio_SRC_Au197_pid_acc_kin_full = 0;
	file_SRC_Au197->GetObject(base211.Data(), H1_Ef_ratio_SRC_Au197_pid_acc_kin_full);
	TH1F* H1_Pf_ratio_SRC_Au197_pid_acc_kin_full = 0;
	file_SRC_Au197->GetObject(base221.Data(), H1_Pf_ratio_SRC_Au197_pid_acc_kin_full);
	TH1F* H1_the_ratio_SRC_Au197_pid_acc_kin_full = 0;
	file_SRC_Au197->GetObject(base23.Data(), H1_the_ratio_SRC_Au197_pid_acc_kin_full);
	TH1F* H1_thp_ratio_SRC_Au197_pid_acc_kin_full = 0;
	file_SRC_Au197->GetObject(base24.Data(), H1_thp_ratio_SRC_Au197_pid_acc_kin_full);
	TH1F* H1_q_ratio_SRC_Au197_pid_acc_kin_full = 0;
	file_SRC_Au197->GetObject(base25.Data(), H1_q_ratio_SRC_Au197_pid_acc_kin_full);
	TH1F* H1_nu_ratio_SRC_Au197_pid_acc_kin_full = 0;
	file_SRC_Au197->GetObject(base26.Data(), H1_nu_ratio_SRC_Au197_pid_acc_kin_full);



	file_SRC_Ca40->GetObject(base3.Data(), H1_Ca40_sim);
	/*for (int i = 1; i <= H1_pmiss_ratio_SRC_Ca40_pid_acc_kin->GetNbinsX(); i++) {
		double x = H1_pmiss_ratio_SRC_Ca40_pid_acc_kin->GetXaxis()->GetBinCenter(i);
		double cont = H1_Ca40_sim->GetBinContent(i);
		// set now content in histogram b using contend of histogram A as weight
		H1_Ca40_sim->Fill(x, cont);
	}*/

	file_SRC_Ca48->GetObject(base3.Data(), H1_Ca48_sim);
	/*for (int i = 1; i <= H1_pmiss_ratio_SRC_Ca48_pid_acc_kin->GetNbinsX(); i++) {
		double x = H1_pmiss_ratio_SRC_Ca48_pid_acc_kin->GetXaxis()->GetBinCenter(i);
		double cont = H1_Ca48_sim->GetBinContent(i);
		// set now content in histogram b using contend of histogram A as weight
		H1_Ca48_sim->Fill(x, cont);
	}*/

	file_SRC_Fe54->GetObject(base3.Data(), H1_Fe54_sim);
	/*for (int i = 1; i <= H1_pmiss_ratio_SRC_Fe54_pid_acc_kin->GetNbinsX(); i++) {
		double x = H1_pmiss_ratio_SRC_Fe54_pid_acc_kin->GetXaxis()->GetBinCenter(i);
		double cont = H1_Fe54_sim->GetBinContent(i);
		// set now content in histogram b using contend of histogram A as weight
		H1_Fe54_sim->Fill(x, cont);
	}*/



	/*
	TString base12 = "PIDList/H1_epctime_ratio_pid_acc_kin_full";
	TString base13 = "PIDList/H1_pCalEtotTrkNorm_ratio_pid_acc_kin_full";
	TString base14 = "HMS_Accp/H1_hdelta_ratio_pid_acc_kin_full";
	TString base15 = "SHMS_Accp/H1_edelta_ratio_pid_acc_kin_full";

	TString base18 = "HMS_Accp/H2_hXColl_hYColl_ratio_pid_acc_kin_full";
	TString base19 = "SHMS_Accp/H2_eXColl_eYColl_ratio_pid_acc_kin_full";
	TString base16 = "Prim_Kin/H2_Em_Pm_ratio_pid_acc_kin_full";
	TString base17 = "Prim_Kin/H2_xbj_Q2_ratio_pid_acc_kin_full";
	*/
	
	
	///Subtract for B10
	//MF
	Double_t B10_MF_scale = 0;
	B10_MF_scale = ((0.5762 - 0.4432) / 0.5738) * (120.059 / 122.15);
	H1_pmiss_ratio_MF_B10_pid_acc_kin->Add(H1_pmiss_ratio_MF_C12_pid_acc_kin, -B10_MF_scale);
	H1_pmiss_ratio_MF_B10_pid_acc_kin_full->Add(H1_pmiss_ratio_MF_C12_pid_acc_kin_full, -B10_MF_scale);//Q2
	H1_Q2_ratio_MF_B10_pid_acc_kin_full->Add(H1_Q2_ratio_MF_C12_pid_acc_kin_full, -B10_MF_scale);
	H1_Em_ratio_MF_B10_pid_acc_kin_full->Add(H1_Em_ratio_MF_C12_pid_acc_kin_full, -B10_MF_scale);
	H1_W_ratio_MF_B10_pid_acc_kin_full->Add(H1_W_ratio_MF_C12_pid_acc_kin_full, -B10_MF_scale);
	H1_epctime_ratio_MF_B10_pid_acc_kin_full->Add(H1_epctime_ratio_MF_C12_pid_acc_kin_full, -B10_MF_scale);
	H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf->Add(H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf, -B10_MF_scale);
	H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf->Add(H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf, -B10_MF_scale);
	H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf->Add(H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf, -B10_MF_scale);
	H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf->Add(H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf, -B10_MF_scale);
	H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->Add(H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf, -B10_MF_scale);
	H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Add(H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf, -B10_MF_scale);
	H1_the_ratio_MF_B10_pid_acc_kin_full->Add(H1_the_ratio_MF_C12_pid_acc_kin_full, -B10_MF_scale);
	H1_thp_ratio_MF_B10_pid_acc_kin_full->Add(H1_thp_ratio_MF_C12_pid_acc_kin_full, -B10_MF_scale);
	H1_nu_ratio_MF_B10_pid_acc_kin_full->Add(H1_nu_ratio_MF_C12_pid_acc_kin_full, -B10_MF_scale);
	H1_q_ratio_MF_B10_pid_acc_kin_full->Add(H1_q_ratio_MF_C12_pid_acc_kin_full, -B10_MF_scale);
	//SRC
	Double_t B10_SRC_scale = 0;
	B10_SRC_scale = (0.5762 - 0.4432) / 0.5738 * (1099.007 / 1143.59);
	H1_pmiss_ratio_SRC_B10_pid_acc_kin->Add(H1_pmiss_ratio_SRC_C12_pid_acc_kin, -B10_SRC_scale);
	H1_pmiss_ratio_SRC_B10_pid_acc_kin_full->Add(H1_pmiss_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_thrq_ratio_SRC_B10_pid_acc_kin_full->Add(H1_thrq_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_Q2_ratio_SRC_B10_pid_acc_kin_full->Add(H1_Q2_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_xbj_ratio_SRC_B10_pid_acc_kin_full->Add(H1_xbj_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_Em_ratio_SRC_B10_pid_acc_kin_full->Add(H1_Em_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_W_ratio_SRC_B10_pid_acc_kin_full->Add(H1_W_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_epctime_ratio_SRC_B10_pid_acc_kin_full->Add(H1_epctime_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	//H1_epctime_ratio_SRC_B10_pid_acc_kin_full_alt->Add(H1_epctime_ratio_SRC_C12_pid_acc_kin_full_alt, -B10_SRC_scale);
	H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full->Add(H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_hdelta_ratio_SRC_B10_pid_acc_kin_full->Add(H1_hdelta_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_edelta_ratio_SRC_B10_pid_acc_kin_full->Add(H1_edelta_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_Ef_ratio_SRC_B10_pid_acc_kin_full->Add(H1_Ef_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_Pf_ratio_SRC_B10_pid_acc_kin_full->Add(H1_Pf_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_the_ratio_SRC_B10_pid_acc_kin_full->Add(H1_the_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_thp_ratio_SRC_B10_pid_acc_kin_full->Add(H1_thp_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_nu_ratio_SRC_B10_pid_acc_kin_full->Add(H1_nu_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);
	H1_q_ratio_SRC_B10_pid_acc_kin_full->Add(H1_q_ratio_SRC_C12_pid_acc_kin_full, -B10_SRC_scale);


	///Subtract for B11
	//MF
	Double_t B11_MF_scale = 0;
	//(b11 tgt thick - effective tgt thick) / c12 tgt thick * (b11 charge / c charge)
	B11_MF_scale = ((0.6328 - 0.4972) / 0.5738) * (262.547 / 122.15);
	H1_pmiss_ratio_MF_B11_pid_acc_kin->Add(H1_pmiss_ratio_MF_C12_pid_acc_kin, -B11_MF_scale);
	H1_pmiss_ratio_MF_B11_pid_acc_kin_full->Add(H1_pmiss_ratio_MF_C12_pid_acc_kin_full, -B11_MF_scale);
	H1_Q2_ratio_MF_B11_pid_acc_kin_full->Add(H1_Q2_ratio_MF_C12_pid_acc_kin_full, -B11_MF_scale);
	H1_Em_ratio_MF_B11_pid_acc_kin_full->Add(H1_Em_ratio_MF_C12_pid_acc_kin_full, -B11_MF_scale);
	H1_W_ratio_MF_B11_pid_acc_kin_full->Add(H1_W_ratio_MF_C12_pid_acc_kin_full, -B11_MF_scale);
	H1_epctime_ratio_MF_B11_pid_acc_kin_full->Add(H1_epctime_ratio_MF_C12_pid_acc_kin_full, -B11_MF_scale);
	H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf->Add(H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf, -B11_MF_scale);
	H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf->Add(H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf, -B11_MF_scale);
	H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf->Add(H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf, -B11_MF_scale);
	H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf->Add(H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf, -B11_MF_scale);
	H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->Add(H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf, -B11_MF_scale);
	H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Add(H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf, -B11_MF_scale);
	H1_the_ratio_MF_B11_pid_acc_kin_full->Add(H1_the_ratio_MF_C12_pid_acc_kin_full, -B11_MF_scale);
	H1_thp_ratio_MF_B11_pid_acc_kin_full->Add(H1_thp_ratio_MF_C12_pid_acc_kin_full, -B11_MF_scale);
	H1_nu_ratio_MF_B11_pid_acc_kin_full->Add(H1_nu_ratio_MF_C12_pid_acc_kin_full, -B11_MF_scale);
	H1_q_ratio_MF_B11_pid_acc_kin_full->Add(H1_q_ratio_MF_C12_pid_acc_kin_full, -B11_MF_scale);
	//SRC	
	Double_t B11_SRC_scale = 0;
	B11_SRC_scale = (0.6328 - 0.4972) / 0.5738 * (1108.109 / 1143.59);
	H1_pmiss_ratio_SRC_B11_pid_acc_kin->Add(H1_pmiss_ratio_SRC_C12_pid_acc_kin, -B11_SRC_scale);
	H1_pmiss_ratio_SRC_B11_pid_acc_kin_full->Add(H1_pmiss_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_thrq_ratio_SRC_B11_pid_acc_kin_full->Add(H1_thrq_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_Q2_ratio_SRC_B11_pid_acc_kin_full->Add(H1_Q2_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_xbj_ratio_SRC_B11_pid_acc_kin_full->Add(H1_xbj_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_Em_ratio_SRC_B11_pid_acc_kin_full->Add(H1_Em_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_W_ratio_SRC_B11_pid_acc_kin_full->Add(H1_W_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_epctime_ratio_SRC_B11_pid_acc_kin_full->Add(H1_epctime_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	//H1_epctime_ratio_SRC_B11_pid_acc_kin_full_alt->Add(H1_epctime_ratio_SRC_C12_pid_acc_kin_full_alt, -B11_SRC_scale);
	H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full->Add(H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_hdelta_ratio_SRC_B11_pid_acc_kin_full->Add(H1_hdelta_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_edelta_ratio_SRC_B11_pid_acc_kin_full->Add(H1_edelta_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_Ef_ratio_SRC_B11_pid_acc_kin_full->Add(H1_Ef_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_Pf_ratio_SRC_B11_pid_acc_kin_full->Add(H1_Pf_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_the_ratio_SRC_B11_pid_acc_kin_full->Add(H1_the_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_thp_ratio_SRC_B11_pid_acc_kin_full->Add(H1_thp_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_nu_ratio_SRC_B11_pid_acc_kin_full->Add(H1_nu_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);
	H1_q_ratio_SRC_B11_pid_acc_kin_full->Add(H1_q_ratio_SRC_C12_pid_acc_kin_full, -B11_SRC_scale);


	///Subtract for Ca40
	//MF
	Double_t Ca40_MF_oil_scale = 6.0 / 7.0 * (1.0 - 0.99484) * (148.587 / 122.15);
	H1_pmiss_ratio_MF_Ca40_pid_acc_kin->Add(H1_pmiss_ratio_MF_C12_pid_acc_kin, -Ca40_MF_oil_scale);
	H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->Add(H1_pmiss_ratio_MF_C12_pid_acc_kin_full, -Ca40_MF_oil_scale);
	H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->Add(H1_Q2_ratio_MF_C12_pid_acc_kin_full, -Ca40_MF_oil_scale);
	H1_Em_ratio_MF_Ca40_pid_acc_kin_full->Add(H1_Em_ratio_MF_C12_pid_acc_kin_full, -Ca40_MF_oil_scale);
	H1_epctime_ratio_MF_Ca40_pid_acc_kin_full->Add(H1_epctime_ratio_MF_C12_pid_acc_kin_full, -Ca40_MF_oil_scale);
	H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->Add(H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf , -Ca40_MF_oil_scale);
	H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->Add(H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf, -Ca40_MF_oil_scale);
	H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->Add(H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf, -Ca40_MF_oil_scale);
	H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->Add(H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf, -Ca40_MF_oil_scale);
	H1_W_ratio_MF_Ca40_pid_acc_kin_full->Add(H1_W_ratio_MF_C12_pid_acc_kin_full, -Ca40_MF_oil_scale);
	H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->Add(H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf, -Ca40_MF_oil_scale);
	H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Add(H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf, -Ca40_MF_oil_scale);
	H1_the_ratio_MF_Ca40_pid_acc_kin_full->Add(H1_the_ratio_MF_C12_pid_acc_kin_full, -Ca40_MF_oil_scale);
	H1_thp_ratio_MF_Ca40_pid_acc_kin_full->Add(H1_thp_ratio_MF_C12_pid_acc_kin_full, -Ca40_MF_oil_scale);
	H1_nu_ratio_MF_Ca40_pid_acc_kin_full->Add(H1_nu_ratio_MF_C12_pid_acc_kin_full, -Ca40_MF_oil_scale);
	H1_q_ratio_MF_Ca40_pid_acc_kin_full->Add(H1_q_ratio_MF_C12_pid_acc_kin_full, -Ca40_MF_oil_scale);
	//SRC
	Double_t Ca40_SRC_oil_scale = 6.0 / 7.0 * (1.0 - 0.99484) * (2354.885 / 1143.59);
	H1_pmiss_ratio_SRC_Ca40_pid_acc_kin->Add(H1_pmiss_ratio_SRC_C12_pid_acc_kin, -Ca40_SRC_oil_scale);
	H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_Q2_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_xbj_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_pmiss_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_thrq_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_Em_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_W_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_W_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_epctime_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	//H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full_alt->Add(H1_epctime_ratio_SRC_C12_pid_acc_kin_full_alt, -Ca40_SRC_oil_scale);
	H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_hdelta_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_edelta_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_Ef_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_Pf_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_the_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_the_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_thp_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_nu_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_nu_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);
	H1_q_ratio_SRC_Ca40_pid_acc_kin_full->Add(H1_q_ratio_SRC_C12_pid_acc_kin_full, -Ca40_SRC_oil_scale);

	///Subtract for Ca48
	//add oil
	//MF
	Double_t Ca48_MF_oil_scale = 6.0 / 7.0 * (1-0.994)*(161.123 / 112.15);
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin->Add(H1_pmiss_ratio_MF_C12_pid_acc_kin, -Ca48_MF_oil_scale);
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_pmiss_ratio_MF_C12_pid_acc_kin_full, -Ca48_MF_oil_scale);
	H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_Q2_ratio_MF_C12_pid_acc_kin_full, -Ca48_MF_oil_scale);
	H1_Em_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_Em_ratio_MF_C12_pid_acc_kin_full, -Ca48_MF_oil_scale);
	H1_epctime_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_epctime_ratio_MF_C12_pid_acc_kin_full, -Ca48_MF_oil_scale);
	H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf, -Ca48_MF_oil_scale);
	H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf, -Ca48_MF_oil_scale);
	H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf, -Ca48_MF_oil_scale);
	H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf, -Ca48_MF_oil_scale);
	H1_W_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_W_ratio_MF_C12_pid_acc_kin_full, -Ca48_MF_oil_scale);
	H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf, -Ca48_MF_oil_scale);
	H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf, -Ca48_MF_oil_scale);
	H1_the_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_the_ratio_MF_C12_pid_acc_kin_full, -Ca48_MF_oil_scale);
	H1_thp_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_thp_ratio_MF_C12_pid_acc_kin_full, -Ca48_MF_oil_scale);
	H1_nu_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_nu_ratio_MF_C12_pid_acc_kin_full, -Ca48_MF_oil_scale);
	H1_q_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_q_ratio_MF_C12_pid_acc_kin_full, -Ca48_MF_oil_scale);
	//SRC
	Double_t Ca48_SRC_oil_scale = 6.0 / 7.0 * (1.0 - 0.9907) * (2710.87 / 1143.59);
	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin->Add(H1_pmiss_ratio_SRC_C12_pid_acc_kin, -Ca48_SRC_oil_scale);
	H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_Q2_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_xbj_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_pmiss_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_thrq_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_Em_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_Em_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_W_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_W_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_epctime_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	//H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full_alt->Add(H1_epctime_ratio_SRC_C12_pid_acc_kin_full_alt, -Ca48_SRC_oil_scale);
	H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_hdelta_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_edelta_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_Ef_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_Pf_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_the_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_the_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_thp_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_thp_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_nu_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_nu_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);
	H1_q_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_q_ratio_SRC_C12_pid_acc_kin_full, -Ca48_SRC_oil_scale);


	//MF
	Double_t Ca48_MF_scale = 0;
	Ca48_MF_scale = (1.0509 - 0.9616) / 0.7851 * (161.123 / 148.587);
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin->Add(H1_pmiss_ratio_MF_Ca40_pid_acc_kin, -Ca48_MF_scale);
	H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_Q2_ratio_MF_Ca40_pid_acc_kin_full, -Ca48_MF_scale);
	H1_Em_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_Em_ratio_MF_Ca40_pid_acc_kin_full, -Ca48_MF_scale);
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full, -Ca48_MF_scale);
	H1_epctime_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_epctime_ratio_MF_Ca40_pid_acc_kin_full, -Ca48_MF_scale);
	H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf, -Ca48_MF_scale);
	H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf, -Ca48_MF_scale);
	H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf, -Ca48_MF_scale);
	H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf, -Ca48_MF_scale);
	H1_W_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_W_ratio_MF_Ca40_pid_acc_kin_full, -Ca48_MF_scale);
	H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf, -Ca48_MF_scale);
	H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Add(H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf, -Ca48_MF_scale);
	H1_the_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_the_ratio_MF_Ca40_pid_acc_kin_full, -Ca48_MF_scale);
	H1_thp_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_thp_ratio_MF_Ca40_pid_acc_kin_full, -Ca48_MF_scale);
	H1_nu_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_nu_ratio_MF_Ca40_pid_acc_kin_full, -Ca48_MF_scale);
	H1_q_ratio_MF_Ca48_pid_acc_kin_full->Add(H1_q_ratio_MF_Ca40_pid_acc_kin_full, -Ca48_MF_scale);
	//SRC
	Double_t Ca48_SRC_scale = 0;
	Ca48_SRC_scale = (1.0509 - 0.9616) / 0.7851 * (2710.868 / 2354.885);
	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin->Add(H1_pmiss_ratio_SRC_Ca40_pid_acc_kin, -Ca48_SRC_scale);
	H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_Em_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_Em_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_W_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_W_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	//H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full_alt->Add(H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full_alt, -Ca48_SRC_scale);
	H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_the_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_the_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_thp_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_thp_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_nu_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_nu_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);
	H1_q_ratio_SRC_Ca48_pid_acc_kin_full->Add(H1_q_ratio_SRC_Ca40_pid_acc_kin_full, -Ca48_SRC_scale);

	
	
	


	//final scaling per proton
	//rad corr * charge * trans * Z/A * areal_den
	/*H1_pmiss_ratio_MF_Be9_pid_acc_kin->Scale(1 / (0.618 * 128.669 * 0.4807 * 0.4444 * 0.9859));
	H1_pmiss_ratio_MF_B10_pid_acc_kin->Scale(1 / (0.618 * 120.059 * 0.4642 * 0.5 * 0.4432));
	H1_pmiss_ratio_MF_B11_pid_acc_kin->Scale(1 / (0.618 * 262.547 * 0.4496 * 0.45455 * 0.4972));
	H1_pmiss_ratio_MF_C12_pid_acc_kin->Scale(1 / (0.618 * 122.15 * 0.4368 * 0.5 * 0.5738));
	H1_pmiss_ratio_MF_Ca40_pid_acc_kin->Scale(1 / (0.577 * 148.587 * 0.2924 * 0.5 * 0.7851));
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin->Scale(1 / (0.577 * 161.123 * 0.2752 * 0.41667 * 0.9616));
	H1_pmiss_ratio_MF_Fe54_pid_acc_kin->Scale(1 / (0.577 * 285.417 * 0.2646 * 0.48148 * 0.367));
	H1_pmiss_ratio_MF_Au197_pid_acc_kin->Scale(1 / (0.451 * 305.285 * 0.1719 * 0.40102 * 0.4047));

	H1_pmiss_ratio_SRC_Be9_pid_acc_kin->Scale(1 / (0.742 * 982.403 * 0.4807 * 0.44444 * 0.9859));
	H1_pmiss_ratio_SRC_B10_pid_acc_kin->Scale(1 / (0.742 * 1099.007 * 0.4642 * 0.5 * 0.4432));
	H1_pmiss_ratio_SRC_B11_pid_acc_kin->Scale(1 / (0.742 * 1108.109 * 0.4496 * 0.45455 * 0.4972));
	H1_pmiss_ratio_SRC_C12_pid_acc_kin->Scale(1 / (0.742 * 1143.59 * 0.4368 * 0.5 * 0.5738));
	H1_pmiss_ratio_SRC_Ca40_pid_acc_kin->Scale(1 / (0.734 * 2354.885 * 0.2924 * 0.5 * 0.7851));
	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin->Scale(1 / (0.734 * 2710.868 * 0.2752 * 0.41667 * 0.9616));
	H1_pmiss_ratio_SRC_Fe54_pid_acc_kin->Scale(1 / (0.734 * 3112.736 * 0.2646 * 0.48148 * 0.367));
	H1_pmiss_ratio_SRC_Au197_pid_acc_kin->Scale(1 / (0.604 * 2408.938 * 0.1719 * 0.40102 * 0.4047));*/
	







	//Glauber transparency per proton
	// (1 / (rad corr * charge * trans * Z/A * areal_den) )
	Double_t Be9_MF_PP = (1.0 / (0.618 * 128.669 * 0.618 * 0.4444 * 0.9859));
	Double_t B10_MF_PP = (1.0 / (0.618 * 120.059 * 0.572 * 0.5 * 0.4432));
	Double_t B11_MF_PP = (1.0 / (0.618 * 262.547 * 0.550 * 0.45455 * 0.4972));
	Double_t C12_MF_PP = (1.0 / (0.618 * 122.15 * 0.527 * 0.5 * 0.5738));
	Double_t Ca40_MF_PP = (1.0 / (0.577 * 148.587 * 0.398 * 0.5 * 0.7851));
	Double_t Ca48_MF_PP = (1.0 / (0.577 * 161.123 * 0.359 * 0.41667 * 0.9616));
	Double_t Fe54_MF_PP = (1.0 / (0.577 * 285.417 * 0.347 * 0.48148 * 0.367));
	Double_t Au197_MF_PP = (1.0 / (0.451 * 305.285 * 0.232 * 0.40102 * 0.4047)); 

	H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PP);
	H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PP);
	H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PP);
	H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PP);
	H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PP);
	H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PP);
	H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PP);
	H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PP);
	H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PP);
	H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PP);
	H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PP);
	H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PP);
	H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PP);
	H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PP);
	H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PP);
	H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PP);



	H1_pmiss_ratio_MF_Be9_pid_acc_kin->Scale(Be9_MF_PP);
	H1_pmiss_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PP);
	H1_Q2_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PP);
	H1_Em_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PP);
	H1_epctime_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PP);
	H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PP);
	H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PP);
	H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PP);
	H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PP);
	H1_W_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PP);
	H1_the_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PP);
	H1_thp_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PP);
	H1_nu_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PP);
	H1_q_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PP);



	H1_pmiss_ratio_MF_B10_pid_acc_kin->Scale(B10_MF_PP);
	H1_pmiss_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PP);
	H1_Q2_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PP);
	H1_Em_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PP);
	H1_epctime_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PP);
	H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PP);
	H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PP);
	H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PP);
	H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PP);
	H1_W_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PP);
	H1_the_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PP);
	H1_thp_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PP);
	H1_nu_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PP);
	H1_q_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PP);



	H1_pmiss_ratio_MF_B11_pid_acc_kin->Scale(B11_MF_PP);
	H1_pmiss_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PP);
	H1_Q2_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PP);
	H1_Em_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PP);
	H1_epctime_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PP);
	H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PP);
	H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PP);
	H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PP);
	H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PP);
	H1_W_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PP);
	H1_the_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PP);
	H1_thp_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PP);
	H1_nu_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PP);
	H1_q_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PP);



	H1_pmiss_ratio_MF_C12_pid_acc_kin->Scale(C12_MF_PP);
	H1_pmiss_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PP);
	H1_Q2_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PP);
	H1_Em_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PP);
	H1_epctime_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PP);
	H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PP);
	H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PP);
	H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PP);
	H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PP);
	H1_W_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PP);
	H1_the_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PP);
	H1_thp_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PP);
	H1_nu_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PP);
	H1_q_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PP);



	H1_pmiss_ratio_MF_Ca40_pid_acc_kin->Scale(Ca40_MF_PP);
	H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PP);
	H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PP);
	H1_Em_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PP);
	H1_epctime_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PP);
	H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PP);
	H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PP);
	H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PP);
	H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PP);
	H1_W_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PP);
	H1_the_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PP);
	H1_thp_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PP);
	H1_nu_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PP);
	H1_q_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PP);



	H1_pmiss_ratio_MF_Ca48_pid_acc_kin->Scale(Ca48_MF_PP);
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PP);
	H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PP);
	H1_Em_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PP);
	H1_epctime_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PP);
	H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PP);
	H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PP);
	H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PP);
	H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PP);
	H1_W_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PP);	
	H1_the_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PP);
	H1_thp_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PP);
	H1_nu_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PP);
	H1_q_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PP);



	H1_pmiss_ratio_MF_Fe54_pid_acc_kin->Scale(Fe54_MF_PP);
	H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PP);
	H1_Q2_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PP);
	H1_Em_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PP);
	H1_epctime_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PP);
	H1_epctime_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PP);
	H1_pCalEtotTrkNorm_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PP);
	H1_hdelta_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PP);
	H1_edelta_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PP);
	H1_W_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PP);
	H1_the_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PP);
	H1_thp_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PP);
	H1_nu_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PP);
	H1_q_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PP);



	H1_pmiss_ratio_MF_Au197_pid_acc_kin->Scale(Au197_MF_PP);
	H1_pmiss_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PP);
	H1_Q2_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PP);
	H1_Em_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PP);
	H1_epctime_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PP);
	H1_epctime_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PP);
	H1_pCalEtotTrkNorm_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PP);
	H1_hdelta_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PP);
	H1_edelta_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PP);
	H1_W_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PP);
	H1_the_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PP);
	H1_thp_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PP);
	H1_nu_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PP);
	H1_q_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PP);



	//(1 / (rad corr * charge * trans * Z/A * areal_den) )
	Double_t Be9_SRC_PP = (1.0 / (0.742 * 982.403 * 0.618 * 0.4444 * 0.9859));
	Double_t B10_SRC_PP = (1.0 / (0.742 * 1099.007 * 0.572 * 0.5 * 0.4432));
	Double_t B11_SRC_PP = (1.0 / (0.742 * 1108.109 * 0.550 * 0.45455 * 0.4972));
	Double_t C12_SRC_PP = (1.0 / (0.742 * 1143.59 * 0.527 * 0.5 * 0.5738));
	Double_t Ca40_SRC_PP = (1.0 / (0.734 * 2354.885 * 0.398 * 0.5 * 0.7851));
	Double_t Ca48_SRC_PP = (1.0 / (0.734 * 2710.868 * 0.359 * 0.41667 * 0.9616));
	Double_t Fe54_SRC_PP = (1.0 / (0.734 * 3112.736 * 0.347 * 0.48148 * 0.367));
	Double_t Au197_SRC_PP = (1.0 / (0.604 * 2408.938 * 0.232 * 0.40102 * 0.4047));

	H1_pmiss_ratio_SRC_Be9_pid_acc_kin->Scale(Be9_SRC_PP);
	H1_pmiss_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_Q2_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_thrq_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_xbj_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_Em_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_W_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_epctime_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	//H1_epctime_ratio_SRC_Be9_pid_acc_kin_full_alt->Scale(Be9_SRC_PP);
	H1_pCalEtotTrkNorm_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_hdelta_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_edelta_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_Ef_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_Pf_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_the_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_thp_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_nu_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);
	H1_q_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PP);

	H1_pmiss_ratio_SRC_B10_pid_acc_kin->Scale(B10_SRC_PP);
	H1_pmiss_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_Q2_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_thrq_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_xbj_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_Em_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_W_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_epctime_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	//H1_epctime_ratio_SRC_B10_pid_acc_kin_full_alt->Scale(B10_SRC_PP);
	H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_hdelta_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_edelta_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_Ef_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_Pf_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_the_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_thp_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_nu_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);
	H1_q_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PP);

	H1_pmiss_ratio_SRC_B11_pid_acc_kin->Scale(B11_SRC_PP);
	H1_pmiss_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_Q2_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_thrq_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_xbj_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_Em_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_W_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_epctime_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	//H1_epctime_ratio_SRC_B11_pid_acc_kin_full_alt->Scale(B11_SRC_PP);
	H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_hdelta_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_edelta_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_Ef_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_Pf_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_the_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_thp_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_nu_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);
	H1_q_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PP);

	H1_pmiss_ratio_SRC_C12_pid_acc_kin->Scale(C12_SRC_PP);
	H1_pmiss_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_Q2_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_thrq_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_xbj_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_Em_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_W_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_epctime_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	//H1_epctime_ratio_SRC_C12_pid_acc_kin_full_alt->Scale(C12_SRC_PP);
	H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_hdelta_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_edelta_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_Ef_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_Pf_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_the_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_thp_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_nu_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);
	H1_q_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PP);

	H1_pmiss_ratio_SRC_Ca40_pid_acc_kin->Scale(Ca40_SRC_PP);
	H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_W_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	//H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full_alt->Scale(Ca40_SRC_PP);
	H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_the_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_nu_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);
	H1_q_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PP);

	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin->Scale(Ca48_SRC_PP);
	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_Em_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_W_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	//H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full_alt->Scale(Ca48_SRC_PP);
	H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_the_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_thp_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_nu_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);
	H1_q_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PP);

	H1_pmiss_ratio_SRC_Fe54_pid_acc_kin->Scale(Fe54_SRC_PP);
	H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_Em_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_W_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	//H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full_alt->Scale(Fe54_SRC_PP);
	H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_Ef_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_Pf_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_the_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_thp_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_nu_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);
	H1_q_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PP);

	H1_pmiss_ratio_SRC_Au197_pid_acc_kin->Scale(Au197_SRC_PP);
	H1_pmiss_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_Q2_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_thrq_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_xbj_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_Em_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_W_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_epctime_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	//H1_epctime_ratio_SRC_Au197_pid_acc_kin_full_alt->Scale(Au197_SRC_PP);
	H1_pCalEtotTrkNorm_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_hdelta_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_edelta_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_Ef_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_Pf_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_the_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_thp_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_nu_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	H1_q_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PP);
	



	//Glauber transparency per nucleus
	// (A / (rad corr * charge * trans * areal_den) )
	/*Double_t Be9_MF_PN = (9.0 / (0.618 * 128.669 * 0.618 * 0.9859));
	Double_t B10_MF_PN = (10.0 / (0.618 * 120.059 * 0.572 * 0.4432));
	Double_t B11_MF_PN = (11.0 / (0.618 * 262.547 * 0.550 * 0.4972));
	Double_t C12_MF_PN = (12.0 / (0.618 * 122.15 * 0.527 * 0.5738));
	Double_t Ca40_MF_PN = (40.0 / (0.577 * 148.587 * 0.398 * 0.7851));
	Double_t Ca48_MF_PN = (48.0 / (0.577 * 161.123 * 0.359 * 0.9616));
	Double_t Fe54_MF_PN = (54.0 / (0.577 * 285.417 * 0.347 * 0.367));
	Double_t Au197_MF_PN = (197.0 / (0.451 * 305.285 * 0.232 * 0.4047)); 
	
	H1_pmiss_ratio_MF_Be9_pid_acc_kin->Scale(Be9_MF_PN);
	H1_Q2_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PN);
	H1_Em_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PN);
	H1_pmiss_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PN);
	H1_epctime_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PN);
	H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PN);
	H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PN);
	H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PN);
	H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PN);
	H1_W_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PN);
	H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PN);
	H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_PN);
	H1_the_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PN);
	H1_thp_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PN);
	H1_nu_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PN);
	H1_q_ratio_MF_Be9_pid_acc_kin_full->Scale(Be9_MF_PN);

	H1_pmiss_ratio_MF_B10_pid_acc_kin->Scale(B10_MF_PN);
	H1_Q2_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PN);
	H1_Em_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PN);
	H1_pmiss_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PN);
	H1_epctime_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PN);
	H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PN);
	H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PN);
	H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PN);
	H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PN);
	H1_W_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PN);
	H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PN);
	H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_PN);
	H1_the_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PN);
	H1_thp_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PN);
	H1_nu_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PN);
	H1_q_ratio_MF_B10_pid_acc_kin_full->Scale(B10_MF_PN);

	H1_pmiss_ratio_MF_B11_pid_acc_kin->Scale(B11_MF_PN);
	H1_Q2_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PN);
	H1_Em_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PN);
	H1_pmiss_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PN);
	H1_epctime_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PN);
	H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PN);
	H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PN);
	H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PN);
	H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PN);
	H1_W_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PN);
	H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PN);
	H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_PN);
	H1_the_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PN);
	H1_thp_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PN);
	H1_nu_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PN);
	H1_q_ratio_MF_B11_pid_acc_kin_full->Scale(B11_MF_PN);

	H1_pmiss_ratio_MF_C12_pid_acc_kin->Scale(C12_MF_PN);
	H1_Q2_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PN);
	H1_Em_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PN);
	H1_pmiss_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PN);
	H1_epctime_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PN);
	H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PN);
	H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PN);
	H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PN);
	H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PN);
	H1_W_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PN);
	H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PN);
	H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_PN);
	H1_the_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PN);
	H1_thp_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PN);
	H1_nu_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PN);
	H1_q_ratio_MF_C12_pid_acc_kin_full->Scale(C12_MF_PN);

	H1_pmiss_ratio_MF_Ca40_pid_acc_kin->Scale(Ca40_MF_PN);
	H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PN);
	H1_Em_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PN);
	H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PN);
	H1_epctime_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PN);
	H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PN);
	H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PN);
	H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PN);
	H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PN);
	H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PN);
	H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_PN);
	H1_the_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PN);
	H1_thp_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PN);
	H1_nu_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PN);
	H1_q_ratio_MF_Ca40_pid_acc_kin_full->Scale(Ca40_MF_PN);
	
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin->Scale(Ca48_MF_PN);
	H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PN);
	H1_Em_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PN);
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PN);
	H1_epctime_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PN);
	H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PN);
	H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PN);
	H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PN);
	H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PN);
	H1_W_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PN);
	H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PN);
	H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_PN);
	H1_the_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PN);
	H1_thp_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PN);
	H1_nu_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PN);
	H1_q_ratio_MF_Ca48_pid_acc_kin_full->Scale(Ca48_MF_PN);

	H1_pmiss_ratio_MF_Fe54_pid_acc_kin->Scale(Fe54_MF_PN);
	H1_Q2_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PN);
	H1_Em_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PN);	
	H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PN);
	H1_epctime_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PN);
	H1_epctime_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PN);
	H1_pCalEtotTrkNorm_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PN);
	H1_hdelta_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PN);
	H1_edelta_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PN);
	H1_W_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PN);
	H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PN);
	H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_PN);
	H1_the_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PN);
	H1_thp_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PN);
	H1_nu_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PN);
	H1_q_ratio_MF_Fe54_pid_acc_kin_full->Scale(Fe54_MF_PN);

	H1_pmiss_ratio_MF_Au197_pid_acc_kin->Scale(Au197_MF_PN);
	H1_Q2_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PN);
	H1_Em_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PN);
	H1_pmiss_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PN);
	H1_epctime_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PN);
	H1_epctime_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PN);
	H1_pCalEtotTrkNorm_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PN);
	H1_hdelta_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PN);
	H1_edelta_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PN);
	H1_W_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PN);
	H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PN);
	H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_PN);
	H1_the_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PN);
	H1_thp_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PN);
	H1_nu_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PN);
	H1_q_ratio_MF_Au197_pid_acc_kin_full->Scale(Au197_MF_PN);


	//(A / (rad corr * charge * trans * areal_den) )
	Double_t Be9_SRC_PN = (9.0 / (0.742 * 982.403 * 0.618 * 0.9859));
	Double_t B10_SRC_PN = (10.0 / (0.742 * 1099.007 * 0.573 * 0.4432));
	Double_t B11_SRC_PN = (11.0 / (0.742 * 1108.109 * 0.550 * 0.4972));
	Double_t C12_SRC_PN = (12.0 / (0.742 * 1143.59 * 0.527 * 0.5738));
	Double_t Ca40_SRC_PN = (40.0 / (0.734 * 2354.885 * 0.396 * 0.7851));
	Double_t Ca48_SRC_PN = (48.0 / (0.734 * 2710.868 * 0.359 * 0.9616));
	Double_t Fe54_SRC_PN = (54.0 / (0.734 * 3112.736 * 0.347 * 0.367));
	Double_t Au197_SRC_PN = (197.0 / (0.604 * 2408.938 * 0.232 * 0.4047));

	H1_pmiss_ratio_SRC_Be9_pid_acc_kin->Scale(Be9_SRC_PN);
	H1_Q2_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_xbj_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_pmiss_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_thrq_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_Em_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_W_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_epctime_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	//H1_epctime_ratio_SRC_Be9_pid_acc_kin_full_alt->Scale(Be9_SRC_PN);
	H1_pCalEtotTrkNorm_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_hdelta_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_edelta_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_Ef_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_Pf_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_the_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_thp_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_nu_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);
	H1_q_ratio_SRC_Be9_pid_acc_kin_full->Scale(Be9_SRC_PN);

	H1_pmiss_ratio_SRC_B10_pid_acc_kin->Scale(B10_SRC_PN);
	H1_Q2_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_xbj_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_pmiss_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_thrq_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_Em_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_W_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_epctime_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	//H1_epctime_ratio_SRC_B10_pid_acc_kin_full_alt->Scale(B10_SRC_PN);
	H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_hdelta_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_edelta_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_Ef_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_Pf_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_the_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_thp_ratio_SRC_B10_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_nu_ratio_SRC_Be9_pid_acc_kin_full->Scale(B10_SRC_PN);
	H1_q_ratio_SRC_Be9_pid_acc_kin_full->Scale(B10_SRC_PN);

	H1_pmiss_ratio_SRC_B11_pid_acc_kin->Scale(B11_SRC_PN);
	H1_Q2_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_xbj_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_pmiss_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_thrq_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_Em_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_W_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_epctime_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	//H1_epctime_ratio_SRC_B11_pid_acc_kin_full_alt->Scale(B11_SRC_PN);
	H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_hdelta_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_edelta_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_Ef_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_Pf_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_the_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_thp_ratio_SRC_B11_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_nu_ratio_SRC_Be9_pid_acc_kin_full->Scale(B11_SRC_PN);
	H1_q_ratio_SRC_Be9_pid_acc_kin_full->Scale(B11_SRC_PN);

	H1_pmiss_ratio_SRC_C12_pid_acc_kin->Scale(C12_SRC_PN);
	H1_Q2_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_xbj_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_pmiss_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_thrq_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_Em_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_W_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_epctime_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	//H1_epctime_ratio_SRC_C12_pid_acc_kin_full_alt->Scale(C12_SRC_PN);
	H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_hdelta_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_edelta_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_Ef_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_Pf_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_the_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_thp_ratio_SRC_C12_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_nu_ratio_SRC_Be9_pid_acc_kin_full->Scale(C12_SRC_PN);
	H1_q_ratio_SRC_Be9_pid_acc_kin_full->Scale(C12_SRC_PN);

	H1_pmiss_ratio_SRC_Ca40_pid_acc_kin->Scale(Ca40_SRC_PN);
	H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_W_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	//H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full_alt->Scale(Ca40_SRC_PN);
	H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_the_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_nu_ratio_SRC_Be9_pid_acc_kin_full->Scale(Ca40_SRC_PN);
	H1_q_ratio_SRC_Be9_pid_acc_kin_full->Scale(Ca40_SRC_PN);

	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin->Scale(Ca48_SRC_PN);
	H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_Em_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_W_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	//H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full_alt->Scale(Ca48_SRC_PN);
	H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_the_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_thp_ratio_SRC_Ca48_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_nu_ratio_SRC_Be9_pid_acc_kin_full->Scale(Ca48_SRC_PN);
	H1_q_ratio_SRC_Be9_pid_acc_kin_full->Scale(Ca48_SRC_PN);

	H1_pmiss_ratio_SRC_Fe54_pid_acc_kin->Scale(Fe54_SRC_PN);
	H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_Em_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_W_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	//H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full_alt->Scale(Fe54_SRC_PN);
	H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_Ef_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_Pf_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_the_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_thp_ratio_SRC_Fe54_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_nu_ratio_SRC_Be9_pid_acc_kin_full->Scale(Fe54_SRC_PN);
	H1_q_ratio_SRC_Be9_pid_acc_kin_full->Scale(Fe54_SRC_PN);

	H1_pmiss_ratio_SRC_Au197_pid_acc_kin->Scale(Au197_SRC_PN);
	H1_Q2_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_xbj_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_pmiss_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_thrq_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_Em_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_W_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_epctime_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	//H1_epctime_ratio_SRC_Au197_pid_acc_kin_full_alt->Scale(Au197_SRC_PN);
	H1_pCalEtotTrkNorm_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_hdelta_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_edelta_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_Ef_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_Pf_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_the_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_thp_ratio_SRC_Au197_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_nu_ratio_SRC_Be9_pid_acc_kin_full->Scale(Au197_SRC_PN);
	H1_q_ratio_SRC_Be9_pid_acc_kin_full->Scale(Au197_SRC_PN);
	*/



	/*H1_pmiss_ratio_MF_Be9_pid_acc_kin->Scale(1.0 / (0.618 * 128.669 * 0.618 * 0.9859));
	H1_pmiss_ratio_MF_B10_pid_acc_kin->Scale(1.0 / (0.618 * 120.059 * 0.572 * 0.4432));
	H1_pmiss_ratio_MF_B11_pid_acc_kin->Scale(1.0 / (0.618 * 262.547 * 0.550 * 0.4972));
	H1_pmiss_ratio_MF_C12_pid_acc_kin->Scale(1.0 / (0.618 * 122.15 * 0.527 * 0.5738));
	H1_pmiss_ratio_MF_Ca40_pid_acc_kin->Scale(1.0 / (0.577 * 148.587 * 0.398 * 0.7851));
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin->Scale(1.0 / (0.577 * 161.123 * 0.359 * 0.9616));
	H1_pmiss_ratio_MF_Fe54_pid_acc_kin->Scale(1.0 / (0.577 * 285.417 * 0.347 * 0.367));
	H1_pmiss_ratio_MF_Au197_pid_acc_kin->Scale(1.0 / (0.451 * 305.285 * 0.232 * 0.4047));

	H1_pmiss_ratio_SRC_Be9_pid_acc_kin->Scale(1.0 / (0.742 * 982.403 * 0.618 * 0.9859));
	H1_pmiss_ratio_SRC_B10_pid_acc_kin->Scale(1.0 / (0.742 * 1099.007 * 0.572 * 0.4432));
	H1_pmiss_ratio_SRC_B11_pid_acc_kin->Scale(1.0 / (0.742 * 1108.109 * 0.550 * 0.4972));
	H1_pmiss_ratio_SRC_C12_pid_acc_kin->Scale(1.0 / (0.742 * 1143.59 * 0.527 * 0.5738));
	H1_pmiss_ratio_SRC_Ca40_pid_acc_kin->Scale(1.0 / (0.734 * 2354.885 * 0.398 * 0.7851));
	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin->Scale(1.0 / (0.734 * 2710.868 * 0.359 * 0.9616));
	H1_pmiss_ratio_SRC_Fe54_pid_acc_kin->Scale(1.0 / (0.734 * 3112.736 * 0.347 * 0.367));
	H1_pmiss_ratio_SRC_Au197_pid_acc_kin->Scale(1.0 / (0.604 * 2408.938 * 0.232 * 0.4047));*/

	
	/*H1_pmiss_ratio_MF_Be9_pid_acc_kin->Scale(9.0 / (0.618 * 128.669 * 0.4807 * 0.9859));
	H1_pmiss_ratio_MF_B10_pid_acc_kin->Scale(10.0 / (0.618 * 120.059 * 0.4642 * 0.4432));
	H1_pmiss_ratio_MF_B11_pid_acc_kin->Scale(11.0 / (0.618 * 262.547 * 0.4496 * 0.4972));
	H1_pmiss_ratio_MF_C12_pid_acc_kin->Scale(12.0 / (0.618 * 122.15 * 0.4368 * 0.5738));
	H1_pmiss_ratio_MF_Ca40_pid_acc_kin->Scale(40.0 / (0.577 * 148.587 * 0.2924 * 0.7851));
	H1_pmiss_ratio_MF_Ca48_pid_acc_kin->Scale(48.0 / (0.577 * 161.123 * 0.2752 * 0.9616));
	H1_pmiss_ratio_MF_Fe54_pid_acc_kin->Scale(54.0 / (0.577 * 285.417 * 0.2646 * 0.367));
	H1_pmiss_ratio_MF_Au197_pid_acc_kin->Scale(197.0 / (0.451 * 305.285 * 0.1719 * 0.4047));

	H1_pmiss_ratio_SRC_Be9_pid_acc_kin->Scale(9.0 / (0.742 * 982.403 * 0.4807 * 0.9859));
	H1_pmiss_ratio_SRC_B10_pid_acc_kin->Scale(10.0 / (0.742 * 1099.007 * 0.4642 * 0.4432));
	H1_pmiss_ratio_SRC_B11_pid_acc_kin->Scale(11.0 / (0.742 * 1108.109 * 0.4496 * 0.4972));
	H1_pmiss_ratio_SRC_C12_pid_acc_kin->Scale(12.0 / (0.742 * 1143.59 * 0.4368 * 0.5738));
	H1_pmiss_ratio_SRC_Ca40_pid_acc_kin->Scale(40.0 / (0.734 * 2354.885 * 0.2924 * 0.7851));
	H1_pmiss_ratio_SRC_Ca48_pid_acc_kin->Scale(48.0 / (0.734 * 2710.868 * 0.2752 * 0.9616));
	H1_pmiss_ratio_SRC_Fe54_pid_acc_kin->Scale(54.0 / (0.734 * 3112.736 * 0.2646 * 0.367));
	H1_pmiss_ratio_SRC_Au197_pid_acc_kin->Scale(197.0 / (0.604 * 2408.938 * 0.1719 * 0.4047));*/
	


	

	


	Be9_MF_Ef3 = H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->Integral();
	B10_MF_Ef3 = H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->Integral();
	B11_MF_Ef3 = H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->Integral();
	C12_MF_Ef3 = H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->Integral();
	Ca40_MF_Ef3 = H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->Integral();
	Ca48_MF_Ef3 = H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->Integral();
	Fe54_MF_Ef3 = H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->Integral();
	Au197_MF_Ef3 = H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->Integral();
	
	Be9_MF_Pf3 = H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->Integral();
	B10_MF_Pf3 = H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Integral();
	B11_MF_Pf3 = H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Integral();
	C12_MF_Pf3 = H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->Integral();
	Ca40_MF_Pf3 = H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Integral();
	Ca48_MF_Pf3 = H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Integral();
	Fe54_MF_Pf3 = H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->Integral();
	Au197_MF_Pf3 = H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->Integral();
	
	Be9_MF_Q23 = H1_q_ratio_MF_Be9_pid_acc_kin_full->Integral();
	B10_MF_Q23 = H1_q_ratio_MF_B10_pid_acc_kin_full->Integral();
	B11_MF_Q23 = H1_q_ratio_MF_B11_pid_acc_kin_full->Integral();
	C12_MF_Q23 = H1_q_ratio_MF_C12_pid_acc_kin_full->Integral();
	Ca40_MF_Q23 = H1_q_ratio_MF_Ca40_pid_acc_kin_full->Integral();
	Ca48_MF_Q23 = H1_q_ratio_MF_Ca48_pid_acc_kin_full->Integral();
	Fe54_MF_Q23 = H1_q_ratio_MF_Fe54_pid_acc_kin_full->Integral();
	Au197_MF_Q23 = H1_q_ratio_MF_Au197_pid_acc_kin_full->Integral();

	
	H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->Scale((Be9_MF_Q23 / Be9_MF_Ef3));
	H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->Scale((B10_MF_Q23 / B10_MF_Ef3));
	H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->Scale((B11_MF_Q23 / B11_MF_Ef3));
	H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->Scale((C12_MF_Q23 / C12_MF_Ef3));
	H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale((Ca40_MF_Q23 / Ca40_MF_Ef3));
	H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale((Ca48_MF_Q23 / Ca48_MF_Ef3));
	H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale((Fe54_MF_Q23 / Fe54_MF_Ef3));
	H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->Scale((Au197_MF_Q23 / Au197_MF_Ef3));
	
	/*
	H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->Scale((Be9_MF_Q23 / Be9_MF_Pf3));
	H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Scale((B10_MF_Q23 / B10_MF_Pf3));
	H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Scale((B11_MF_Q23 / B11_MF_Pf3));
	H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->Scale((C12_MF_Q23 / C12_MF_Pf3));
	H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale((Ca40_MF_Q23 / Ca40_MF_Pf3));
	H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale((Ca48_MF_Q23 / Ca48_MF_Pf3));
	H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale((Fe54_MF_Q23 / Fe54_MF_Pf3));
	H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->Scale((Au197_MF_Q23 / Au197_MF_Pf3));
	*/

	/*
	H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->Scale((1 / H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->Integral()));
	H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Scale((1 / H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Integral()));
	H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Scale((1 / H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Integral()));
	H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->Scale((1 / H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->Integral()));
	H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale((1 / H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Integral()));
	H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale((1 / H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Integral()));
	H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale((1 / H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->Integral()));
	H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->Scale((1 / H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->Integral()));
	*/

	H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->Scale((1 / Be9_MF_Pf3));
	H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Scale((1 / B10_MF_Pf3));
	H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Scale((1 / B11_MF_Pf3));
	H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->Scale((1 / C12_MF_Pf3));
	H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale((1 / Ca40_MF_Pf3));
	H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale((1 / Ca48_MF_Pf3));
	H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale((1 / Fe54_MF_Pf3));
	H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->Scale((1 / Au197_MF_Pf3));

	H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->Scale(Be9_MF_Q23);
	H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Scale(B10_MF_Q23);
	H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Scale(B11_MF_Q23);
	H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->Scale(C12_MF_Q23);
	H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Scale(Ca40_MF_Q23);
	H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Scale(Ca48_MF_Q23);
	H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->Scale(Fe54_MF_Q23);
	H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->Scale(Au197_MF_Q23);


	Be9_MF_Ef3 = H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->Integral();
	B10_MF_Ef3 = H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->Integral();
	B11_MF_Ef3 = H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->Integral();
	C12_MF_Ef3 = H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->Integral();
	Ca40_MF_Ef3 = H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->Integral();
	Ca48_MF_Ef3 = H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->Integral();
	Fe54_MF_Ef3 = H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->Integral();
	Au197_MF_Ef3 = H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->Integral();

	Be9_MF_Pf3 = H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->Integral();
	B10_MF_Pf3 = H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->Integral();
	B11_MF_Pf3 = H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->Integral();
	C12_MF_Pf3 = H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->Integral();
	Ca40_MF_Pf3 = H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->Integral();
	Ca48_MF_Pf3 = H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->Integral();
	Fe54_MF_Pf3 = H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->Integral();
	Au197_MF_Pf3 = H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->Integral();


	TH1F* H1_pmiss_MF_Be9_C12 = new TH1F(*H1_pmiss_ratio_MF_Be9_pid_acc_kin); H1_pmiss_MF_Be9_C12->Sumw2(); H1_pmiss_MF_Be9_C12->Divide(H1_pmiss_ratio_MF_C12_pid_acc_kin); H1_pmiss_MF_Be9_C12->SetTitle("MF Be9 to C12");
	//Draw1dE(H1_pmiss_MF_Be9_C12, "H1_pmiss_MF_Be9_C12.png", "H1_pmiss_MF_Be9_C12", outROOT);
	TH1F* H1_pmiss_MF_B10_C12 = new TH1F(*H1_pmiss_ratio_MF_B10_pid_acc_kin); H1_pmiss_MF_B10_C12->Sumw2(); H1_pmiss_MF_B10_C12->Divide(H1_pmiss_ratio_MF_C12_pid_acc_kin); H1_pmiss_MF_B10_C12->SetTitle("MF B10 to C12");
	//Draw1dE(H1_pmiss_MF_B10_C12, "H1_pmiss_MF_B10_C12.png", "H1_pmiss_MF_B10_C12", outROOT);
	TH1F* H1_pmiss_MF_B11_C12 = new TH1F(*H1_pmiss_ratio_MF_B11_pid_acc_kin); H1_pmiss_MF_B11_C12->Sumw2(); H1_pmiss_MF_B11_C12->Divide(H1_pmiss_ratio_MF_C12_pid_acc_kin); H1_pmiss_MF_B11_C12->SetTitle("MF B11 to C12");
	//Draw1dE(H1_pmiss_MF_B11_C12, "H1_pmiss_MF_B11_C12.png", "H1_pmiss_MF_B11_C12", outROOT);
	//TH1F* H1_pmiss_MF_Ca40_C12 = new TH1F(*H1_pmiss_ratio_MF_Ca40_pid_acc_kin); H1_pmiss_MF_Ca40_C12->Sumw2(); H1_pmiss_MF_Ca40_C12->Divide(H1_pmiss_ratio_MF_C12_pid_acc_kin); H1_pmiss_MF_Ca40_C12->SetTitle("MF Ca40 to C12");
	//Draw1dE(H1_pmiss_MF_Ca40_C12, "H1_pmiss_MF_Ca40_C12.png", "H1_pmiss_MF_Ca40_C12", outROOT);
	//TH1F* H1_pmiss_MF_Ca48_C12 = new TH1F(*H1_pmiss_ratio_MF_Ca48_pid_acc_kin); H1_pmiss_MF_Ca48_C12->Sumw2(); H1_pmiss_MF_Ca48_C12->Divide(H1_pmiss_ratio_MF_C12_pid_acc_kin); H1_pmiss_MF_Ca48_C12->SetTitle("MF Ca48 to C12");
	//Draw1dE(H1_pmiss_MF_Ca48_C12, "H1_pmiss_MF_Ca48_C12.png", "H1_pmiss_MF_Ca48_C12", outROOT);
	//TH1F* H1_pmiss_MF_Fe54_C12 = new TH1F(*H1_pmiss_ratio_MF_Fe54_pid_acc_kin); H1_pmiss_MF_Fe54_C12->Sumw2(); H1_pmiss_MF_Fe54_C12->Divide(H1_pmiss_ratio_MF_C12_pid_acc_kin); H1_pmiss_MF_Fe54_C12->SetTitle("MF Fe54 to C12");
	//Draw1dE(H1_pmiss_MF_Fe54_C12, "H1_pmiss_MF_Fe54_C12.png", "H1_pmiss_MF_Fe54_C12", outROOT);
	//TH1F* H1_pmiss_MF_Au197_C12 = new TH1F(*H1_pmiss_ratio_MF_Au197_pid_acc_kin); H1_pmiss_MF_Au197_C12->Sumw2(); H1_pmiss_MF_Au197_C12->Divide(H1_pmiss_ratio_MF_C12_pid_acc_kin); H1_pmiss_MF_Au197_C12->SetTitle("MF Au197 to C12");
	//Draw1dE(H1_pmiss_MF_Au197_C12, "H1_pmiss_MF_Au197_C12.png", "H1_pmiss_MF_Au197_C12", outROOT);

	TH1F* H1_pmiss_SRC_Be9_C12 = new TH1F(*H1_pmiss_ratio_SRC_Be9_pid_acc_kin); H1_pmiss_SRC_Be9_C12->Sumw2(); H1_pmiss_SRC_Be9_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin); H1_pmiss_SRC_Be9_C12->SetTitle("SRC Be9 to C12");
	//Draw1dE(H1_pmiss_SRC_Be9_C12, "H1_pmiss_SRC_Be9_C12.png", "H1_pmiss_SRC_Be9_C12", outROOT);
	TH1F* H1_pmiss_SRC_B10_C12 = new TH1F(*H1_pmiss_ratio_SRC_B10_pid_acc_kin); H1_pmiss_SRC_B10_C12->Sumw2(); H1_pmiss_SRC_B10_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin); H1_pmiss_SRC_B10_C12->SetTitle("SRC B10 to C12");
	//Draw1dE(H1_pmiss_SRC_B10_C12, "H1_pmiss_SRC_B10_C12.png", "H1_pmiss_SRC_B10_C12", outROOT);
	TH1F* H1_pmiss_SRC_B11_C12 = new TH1F(*H1_pmiss_ratio_SRC_B11_pid_acc_kin); H1_pmiss_SRC_B11_C12->Sumw2(); H1_pmiss_SRC_B11_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin); H1_pmiss_SRC_B11_C12->SetTitle("SRC B11 to C12");
	//Draw1dE(H1_pmiss_SRC_B11_C12, "H1_pmiss_SRC_B11_C12.png", "H1_pmiss_SRC_B11_C12", outROOT);
	TH1F* H1_pmiss_SRC_Ca40_C12 = new TH1F(*H1_pmiss_ratio_SRC_Ca40_pid_acc_kin); H1_pmiss_SRC_Ca40_C12->Sumw2(); H1_pmiss_SRC_Ca40_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin); H1_pmiss_SRC_Ca40_C12->SetTitle("SRC Ca40 to C12");
	//Draw1dE(H1_pmiss_SRC_Ca40_C12, "H1_pmiss_SRC_Ca40_C12.png", "H1_pmiss_SRC_Ca40_C12", outROOT);
	TH1F* H1_pmiss_SRC_Ca48_C12 = new TH1F(*H1_pmiss_ratio_SRC_Ca48_pid_acc_kin); H1_pmiss_SRC_Ca48_C12->Sumw2(); H1_pmiss_SRC_Ca48_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin); H1_pmiss_SRC_Ca48_C12->SetTitle("SRC Ca48 to C12");
	//Draw1dE(H1_pmiss_SRC_Ca48_C12, "H1_pmiss_SRC_Ca48_C12.png", "H1_pmiss_SRC_Ca48_C12", outROOT);
	TH1F* H1_pmiss_SRC_Fe54_C12 = new TH1F(*H1_pmiss_ratio_SRC_Fe54_pid_acc_kin); H1_pmiss_SRC_Fe54_C12->Sumw2(); H1_pmiss_SRC_Fe54_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin); H1_pmiss_SRC_Fe54_C12->SetTitle("SRC Fe54 to C12");
	//Draw1dE(H1_pmiss_SRC_Fe54_C12, "H1_pmiss_SRC_Fe54_C12.png", "H1_pmiss_SRC_Fe54_C12", outROOT);
	TH1F* H1_pmiss_SRC_Au197_C12 = new TH1F(*H1_pmiss_ratio_SRC_Au197_pid_acc_kin); H1_pmiss_SRC_Au197_C12->Sumw2(); H1_pmiss_SRC_Au197_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin); H1_pmiss_SRC_Au197_C12->SetTitle("SRC Au197 to C12");
	//Draw1dE(H1_pmiss_SRC_Au197_C12, "H1_pmiss_SRC_Au197_C12.png", "H1_pmiss_SRC_Au197_C12", outROOT);


	TH1F* H1_pmiss_MF_Ca48_Ca40 = new TH1F(*H1_pmiss_ratio_MF_Ca48_pid_acc_kin);
	H1_pmiss_MF_Ca48_Ca40->Sumw2();
	H1_pmiss_MF_Ca48_Ca40->Divide(H1_pmiss_ratio_MF_Ca40_pid_acc_kin);
	H1_pmiss_MF_Ca48_Ca40->SetTitle("MF Ca48 to Ca40");
	//Draw1dE(H1_pmiss_MF_Ca48_Ca40, "H1_pmiss_MF_Ca48_Ca40.png", "H1_pmiss_MF_Ca48_Ca40", outROOT);
	
	TH1F* H1_pmiss_MF_Fe54_Ca40 = new TH1F(*H1_pmiss_ratio_MF_Fe54_pid_acc_kin);
	H1_pmiss_MF_Fe54_Ca40->Sumw2();
	H1_pmiss_MF_Fe54_Ca40->Divide(H1_pmiss_ratio_MF_Ca40_pid_acc_kin);
	H1_pmiss_MF_Fe54_Ca40->SetTitle("MF Fe54 to Ca40");
	//Draw1dE(H1_pmiss_MF_Fe54_Ca40, "H1_pmiss_MF_Fe54_Ca40.png", "H1_pmiss_MF_Fe54_Ca40", outROOT);



	/*TH1F* H1_pmiss_SRC_Be9_C12 = new TH1F(*H1_pmiss_ratio_SRC_Be9_pid_acc_kin);
	H1_pmiss_SRC_Be9_C12->Sumw2();
	H1_pmiss_SRC_Be9_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin);
	H1_pmiss_SRC_Be9_C12->SetTitle("SRC Be9 to C12");
	Draw1dE(H1_pmiss_SRC_Be9_C12, "H1_pmiss_SRC_Be9_C12.png", "H1_pmiss_SRC_Be9_C12", outROOT);

	TH1F* H1_pmiss_SRC_B10_C12 = new TH1F(*H1_pmiss_ratio_SRC_B10_pid_acc_kin);
	H1_pmiss_SRC_B10_C12->Sumw2();
	H1_pmiss_SRC_B10_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin);
	H1_pmiss_SRC_B10_C12->SetTitle("SRC B10 to C12");
	Draw1dE(H1_pmiss_SRC_B10_C12, "H1_pmiss_SRC_B10_C12.png", "H1_pmiss_SRC_B10_C12", outROOT);

	TH1F* H1_pmiss_SRC_B11_C12 = new TH1F(*H1_pmiss_ratio_SRC_B11_pid_acc_kin);
	H1_pmiss_SRC_B11_C12->Sumw2();
	H1_pmiss_SRC_B11_C12->Divide(H1_pmiss_ratio_SRC_C12_pid_acc_kin);
	H1_pmiss_SRC_B11_C12->SetTitle("SRC B11 to C12");
	Draw1dE(H1_pmiss_SRC_B11_C12, "H1_pmiss_SRC_B11_C12.png", "H1_pmiss_SRC_B11_C12", outROOT);
	*/



	TH1F* H1_pmiss_SRC_Ca48_Ca40 = new TH1F(*H1_pmiss_ratio_SRC_Ca48_pid_acc_kin);
	H1_pmiss_SRC_Ca48_Ca40->Sumw2();
	H1_pmiss_SRC_Ca48_Ca40->Divide(H1_pmiss_ratio_SRC_Ca40_pid_acc_kin);
	H1_pmiss_SRC_Ca48_Ca40->SetTitle("SRC Ca48 to Ca40");
	Draw1dE(H1_pmiss_SRC_Ca48_Ca40, "H1_pmiss_SRC_Ca48_Ca40.png", "H1_pmiss_SRC_Ca48_Ca40", outROOT);

	TH1F* H1_pmiss_SRC_Fe54_Ca40 = new TH1F(*H1_pmiss_ratio_SRC_Fe54_pid_acc_kin);
	H1_pmiss_SRC_Fe54_Ca40->Sumw2();
	H1_pmiss_SRC_Fe54_Ca40->Divide(H1_pmiss_ratio_SRC_Ca40_pid_acc_kin);
	H1_pmiss_SRC_Fe54_Ca40->SetTitle("SRC Fe54 to Ca40");
	Draw1dE(H1_pmiss_SRC_Fe54_Ca40, "H1_pmiss_SRC_Fe54_Ca40.png", "H1_pmiss_SRC_Fe54_Ca40", outROOT);

	TH1F* H1_Ca48_Ca40_sim = new TH1F(*H1_Ca48_sim);
	H1_Ca48_Ca40_sim->Sumw2();
	H1_Ca48_Ca40_sim->Divide(H1_Ca40_sim);
	H1_Ca48_Ca40_sim->SetTitle("SRC Ca48 to Ca40");
	//Draw1dE(H1_Ca48_Ca40_sim, "H1_Ca48_Ca40_sim.png", "H1_Ca48_Ca40_sim", outROOT);

	TH1F* H1_Fe54_Ca40_sim = new TH1F(*H1_Fe54_sim);
	H1_Fe54_Ca40_sim->Sumw2();
	H1_Fe54_Ca40_sim->Divide(H1_Ca40_sim);
	H1_Fe54_Ca40_sim->SetTitle("SRC Fe54 to Ca40");
	//Draw1dE(H1_Fe54_Ca40_sim, "H1_Fe54_Ca40_sim.png", "H1_Fe54_Ca40_sim", outROOT);

	TH1F* H1_Fe54_Ca48_sim = new TH1F(*H1_Fe54_sim);
	H1_Fe54_Ca48_sim->Sumw2();
	H1_Fe54_Ca48_sim->Divide(H1_Ca48_sim);
	H1_Fe54_Ca48_sim->SetTitle("SRC Fe54 to Ca48");
	//Draw1dE(H1_Fe54_Ca48_sim, "H1_Fe54_Ca48_sim.png", "H1_Fe54_Ca48_sim", outROOT);







	/*
	TH1F* H1_pmiss_MF_Ca40_Ca48 = new TH1F(*H1_pmiss_ratio_MF_Ca40_pid_acc_kin);
	H1_pmiss_MF_Ca40_Ca48->Sumw2();
	H1_pmiss_MF_Ca40_Ca48->Divide(H1_pmiss_ratio_MF_Ca48_pid_acc_kin);
	H1_pmiss_MF_Ca40_Ca48->SetTitle("MF Ca40 to Ca48");
	//Draw1dE(H1_pmiss_MF_Ca40_Ca48, "H1_pmiss_MF_Ca40_Ca48.png", "H1_pmiss_MF_Ca40_Ca48", outROOT);

	TH1F* H1_pmiss_MF_Fe54_Ca48 = new TH1F(*H1_pmiss_ratio_MF_Fe54_pid_acc_kin);
	H1_pmiss_MF_Fe54_Ca48->Sumw2();
	H1_pmiss_MF_Fe54_Ca48->Divide(H1_pmiss_ratio_MF_Ca48_pid_acc_kin);
	H1_pmiss_MF_Fe54_Ca48->SetTitle("MF Fe54 to Ca48");
	//Draw1dE(H1_pmiss_MF_Fe54_Ca48, "H1_pmiss_MF_Fe54_Ca48.png", "H1_pmiss_MF_Fe54_Ca48", outROOT);

	TH1F* H1_pmiss_SRC_Ca40_Ca48 = new TH1F(*H1_pmiss_ratio_SRC_Ca40_pid_acc_kin);
	H1_pmiss_SRC_Ca40_Ca48->Sumw2();
	H1_pmiss_SRC_Ca40_Ca48->Divide(H1_pmiss_ratio_SRC_Ca48_pid_acc_kin);
	H1_pmiss_SRC_Ca40_Ca48->SetTitle("SRC Ca40 to Ca48");
	//Draw1dE(H1_pmiss_SRC_Ca40_Ca48, "H1_pmiss_SRC_Ca40_Ca48.png", "H1_pmiss_SRC_Ca40_Ca48", outROOT);
	*/
	TH1F* H1_pmiss_SRC_Fe54_Ca48 = new TH1F(*H1_pmiss_ratio_SRC_Fe54_pid_acc_kin);
	H1_pmiss_SRC_Fe54_Ca48->Sumw2();
	H1_pmiss_SRC_Fe54_Ca48->Divide(H1_pmiss_ratio_SRC_Ca48_pid_acc_kin);
	H1_pmiss_SRC_Fe54_Ca48->SetTitle("SRC Fe54 to Ca48");
	//Draw1dE(H1_pmiss_SRC_Fe54_Ca48, "H1_pmiss_SRC_Fe54_Ca48.png", "H1_pmiss_SRC_Fe54_Ca48", outROOT);
	

	//thingpap(H1_Ca48_Ca40_sim, H1_Fe54_Ca40_sim, H1_Fe54_Ca48_sim, outROOT);
	//thingpap(H1_pmiss_SRC_Ca48_Ca40, H1_pmiss_SRC_Fe54_Ca40, H1_pmiss_SRC_Fe54_Ca48, outROOT);
	
	//thingpap2(H1_xbj_ratio_SRC_Ca40_pid_acc_kin, H1_xbj_ratio_SRC_Ca48_pid_acc_kin, H1_xbj_ratio_SRC_Fe54_pid_acc_kin, "H1_pmiss_SRC.png", "H1_pmiss_SRC", outROOT);
	//thingpap2(H1_Q2_ratio_SRC_Ca40_pid_acc_kin, H1_Q2_ratio_SRC_Ca48_pid_acc_kin, H1_Q2_ratio_SRC_Fe54_pid_acc_kin, "H1_Q2_SRC.png", "H1_Q2_SRC", outROOT);
	//thingpap2(H1_pmiss_ratio_SRC_Ca40_pid_acc_kin, H1_pmiss_ratio_SRC_Ca48_pid_acc_kin, H1_pmiss_ratio_SRC_Fe54_pid_acc_kin, "H1_xbj_SRC.png", "H1_xbj_SRC", outROOT);
	/*Draw1dE(H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full, "H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full.png", "H1_xbj_ratio_SRC_Ca40_pid_acc_kin", outROOT);
	Draw1dE(H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full, "H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full.png", "H1_xbj_ratio_SRC_Ca48_pid_acc_kin", outROOT);
	Draw1dE(H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full, "H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full.png", "H1_xbj_ratio_SRC_Fe54_pid_acc_kin", outROOT);
	Draw1dE(H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full, "H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full.png", "H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full", outROOT);
	Draw1dE(H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full, "H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full.png", "H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full", outROOT);
	Draw1dE(H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full, "H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full.png", "H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full", outROOT);
	Draw1dE(H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full, "H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full.png", "H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full", outROOT);
	Draw1dE(H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full, "H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full.png", "H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full", outROOT);
	Draw1dE(H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full, "H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full.png", "H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full", outROOT);
	Draw1dE(H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full, "H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full.png", "H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full", outROOT);
	Draw1dE(H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full, "H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full.png", "H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full", outROOT);
	Draw1dE(H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full, "H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full.png", "H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full", outROOT);
	
	Draw1dE(H1_W_ratio_SRC_D2_pid_acc_kin_full, "H1_W_ratio_SRC_D2_pid_acc_kin_full.png", "H1_W_ratio_SRC_D2_pid_acc_kin_full", outROOT);
	
	//Draw1dE(H1_W_ratio_MF_Ca40_pid_acc_kin_full, "H1_W_ratio_MF_Ca40_pid_acc_kin_full.png", "H1_W_ratio_MF_Ca40_pid_acc_kin_full", outROOT);
	//Draw1dE(H1_W_ratio_MF_Ca48_pid_acc_kin_full, "H1_W_ratio_MF_Ca48_pid_acc_kin_full.png", "H1_W_ratio_MF_Ca48_pid_acc_kin_full", outROOT);
	//Draw1dE(H1_W_ratio_MF_Fe54_pid_acc_kin_full, "H1_W_ratio_MF_Fe54_pid_acc_kin_full.png", "H1_W_ratio_MF_Fe54_pid_acc_kin_full", outROOT);
	Draw1dE(H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full, "H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full.png", "H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full", outROOT);
	Draw1dE(H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full, "H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full.png", "H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full", outROOT);
	Draw1dE(H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full, "H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full.png", "H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full", outROOT);
	*/

	//file_SRC_Ca40->GetObject(base3.Data(), H1_Ca40_sim);
	//file_SRC_Ca48->GetObject(base3.Data(), H1_Ca48_sim);
	//file_SRC_Fe54->GetObject(base3.Data(), H1_Fe54_sim);


	/*thingy(
		H1_pmiss_MF_Be9_C12, H1_pmiss_MF_B10_C12, H1_pmiss_MF_B11_C12,
		H1_pmiss_SRC_Be9_C12, H1_pmiss_SRC_B10_C12, H1_pmiss_SRC_B11_C12,
		H1_pmiss_MF_Ca48_Ca40, H1_pmiss_MF_Fe54_Ca40,
		H1_pmiss_SRC_Ca48_Ca40, H1_pmiss_SRC_Fe54_Ca40,
		outROOT);*/

	/*thingy(
		H1_pmiss_MF_Be9_C12, H1_pmiss_MF_B10_C12, H1_pmiss_MF_B11_C12,
		H1_pmiss_SRC_Be9_C12, H1_pmiss_SRC_B10_C12, H1_pmiss_SRC_B11_C12,
		H1_pmiss_MF_Ca48_Ca40, H1_pmiss_MF_Fe54_Ca40,
		H1_pmiss_SRC_Ca48_Ca40, H1_pmiss_SRC_Fe54_Ca40,
		outROOT);*/

	/*
	thingy(
		H1_Be9_C12_sim, H1_B10_C12_sim, H1_B11_C12_sim,
		H1_pmiss_SRC_Be9_C12, H1_pmiss_SRC_B10_C12, H1_pmiss_SRC_B11_C12,
		H1_pmiss_MF_Ca48_Ca40, H1_pmiss_MF_Fe54_Ca40,
		H1_pmiss_SRC_Ca48_Ca40, H1_pmiss_SRC_Fe54_Ca40,
		outROOT);
	*/
	
	Draw2dh(H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full, "H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full.png", "title", outROOT);
	Draw2de(H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full, "H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full.png", "title", outROOT);
	//Draw2d(H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full, "H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full.png", "title", outROOT);
	Draw2d(H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full, "H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full.png", "title", outROOT);

	Draw2d2(H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf2, "H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full.png", "title", outROOT);
	Draw2d2(H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full_mf2, "H2_Em_Pm_ratio_MF_Ca40_pid_acc_kin_full.png", "title", outROOT);
	Draw2d2(H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full_mf2, "H2_Em_Pm_ratio_MF_Ca48_pid_acc_kin_full.png", "title", outROOT);

	/*

	cout << "Be9 MF 1: " << Be9_MF_1 << endl;
	cout << "B10 MF 1: " << B10_MF_1 << endl;
	cout << "B11 MF 1: " << B11_MF_1 << endl;
	cout << "C12 MF 1: " << C12_MF_1 << endl;
	cout << "Ca40 MF 1: " << Ca40_MF_1 << endl;
	cout << "Ca48 MF 1: " << Ca48_MF_1 << endl;
	cout << "Fe54 MF 1: " << Fe54_MF_1 << endl;
	cout << "Au197 MF 1: " << Au197_MF_1 << endl;
	
	cout << endl;

	cout << "Be9 MF 1: " << Be9_MF_Pm1 << endl;
	cout << "B10 MF 1: " << B10_MF_Pm1 << endl;
	cout << "B11 MF 1: " << B11_MF_Pm1 << endl;
	cout << "C12 MF 1: " << C12_MF_Pm1 << endl;
	cout << "Ca40 MF 1: " << Ca40_MF_Pm1 << endl;
	cout << "Ca48 MF 1: " << Ca48_MF_Pm1 << endl;
	cout << "Fe54 MF 1: " << Fe54_MF_Pm1 << endl;
	cout << "Au197 MF 1: " << Au197_MF_Pm1 << endl;

	cout << endl;

	cout << "Be9 MF 2: " << Be9_MF_2 << endl;
	cout << "B10 MF 2: " << B10_MF_2 << endl;
	cout << "B11 MF 2: " << B11_MF_2 << endl;
	cout << "C12 MF 2: " << C12_MF_2 << endl;
	cout << "Ca40 MF 2: " << Ca40_MF_2 << endl;
	cout << "Ca48 MF 2: " << Ca48_MF_2 << endl;
	cout << "Fe54 MF 2: " << Fe54_MF_2 << endl;
	cout << "Au197 MF 2: " << Au197_MF_2 << endl;

	cout << endl;

	cout << "Be9 MF 2: " << Be9_MF_Pm2 << endl;
	cout << "B10 MF 2: " << B10_MF_Pm2 << endl;
	cout << "B11 MF 2: " << B11_MF_Pm2 << endl;
	cout << "C12 MF 2: " << C12_MF_Pm2 << endl;
	cout << "Ca40 MF 2: " << Ca40_MF_Pm2 << endl;
	cout << "Ca48 MF 2: " << Ca48_MF_Pm2 << endl;
	cout << "Fe54 MF 2: " << Fe54_MF_Pm2 << endl;
	cout << "Au197 MF 2: " << Au197_MF_Pm2 << endl;

	cout << endl;
	*/
	

	cout << "Be9 MF Ef3: " << Be9_MF_Ef3 << endl;
	cout << "B10 MF Ef3: " << B10_MF_Ef3 << endl;
	cout << "B11 MF Ef3: " << B11_MF_Ef3 << endl;
	cout << "C12 MF Ef3: " << C12_MF_Ef3 << endl;
	cout << "Ca40 MF Ef3: " << Ca40_MF_Ef3 << endl;
	cout << "Ca48 MF Ef3: " << Ca48_MF_Ef3 << endl;
	cout << "Fe54 MF Ef3: " << Fe54_MF_Ef3 << endl;
	cout << "Au197 MF EF3: " << Au197_MF_Ef3 << endl;

	cout << endl;

	cout << "Be9 MF Pf3: " << Be9_MF_Pf3 << endl;
	cout << "B10 MF Pf3: " << B10_MF_Pf3 << endl;
	cout << "B11 MF Pf3: " << B11_MF_Pf3 << endl;
	cout << "C12 MF Pf3: " << C12_MF_Pf3 << endl;
	cout << "Ca40 MF Pf3: " << Ca40_MF_Pf3 << endl;
	cout << "Ca48 MF Pf3: " << Ca48_MF_Pf3 << endl;
	cout << "Fe54 MF Pf3: " << Fe54_MF_Pf3 << endl;
	cout << "Au197 MF Pf3: " << Au197_MF_Pf3 << endl;
	
	cout << endl;

	cout << "cent Be9 MF 3: " << Be9_MF_Q23 << endl;
	cout << "B10 MF 3: " << B10_MF_Q23 << endl;
	cout << "B11 MF 3: " << B11_MF_Q23 << endl;
	cout << "C12 MF 3: " << C12_MF_Q23 << endl;
	cout << "Ca40 MF 3: " << Ca40_MF_Q23 << endl;
	cout << "Ca48 MF 3: " << Ca48_MF_Q23 << endl;
	cout << "Fe54 MF 3: " << Fe54_MF_Q23 << endl;
	cout << "Au197 MF 3: " << Au197_MF_Q23 << endl;

	cout << endl;

	
	/*
	cout << "Be9 SRC 1: " << Be9_SRC_1 << endl;
	cout << "B10 SRC 1: " << B10_SRC_1 << endl;
	cout << "B11 SRC 1: " << B11_SRC_1 << endl;
	cout << "C12 SRC 1: " << C12_SRC_1 << endl;
	cout << "Ca40 SRC 1: " << Ca40_SRC_1 << endl;
	cout << "Ca48 SRC 1: " << Ca48_SRC_1 << endl;
	cout << "Fe54 SRC 1: " << Fe54_SRC_1 << endl;
	cout << "Au197 SRC 1: " << Au197_SRC_1 << endl;

	cout << endl;

	cout << "Be9 SRC 1: " << Be9_SRC_Pm1 << endl;
	cout << "B10 SRC 1: " << B10_SRC_Pm1 << endl;
	cout << "B11 SRC 1: " << B11_SRC_Pm1 << endl;
	cout << "C12 SRC 1: " << C12_SRC_Pm1 << endl;
	cout << "Ca40 SRC 1: " << Ca40_SRC_Pm1 << endl;
	cout << "Ca48 SRC 1: " << Ca48_SRC_Pm1 << endl;
	cout << "Fe54 SRC 1: " << Fe54_SRC_Pm1 << endl;
	cout << "Au197 SRC 1: " << Au197_SRC_Pm1 << endl;

	cout << endl;

	cout << "Be9 SRC 2: " << Be9_SRC_2 << endl;
	cout << "B10 SRC 2: " << B10_SRC_2 << endl;
	cout << "B11 SRC 2: " << B11_SRC_2 << endl;
	cout << "C12 SRC 2: " << C12_SRC_2 << endl;
	cout << "Ca40 SRC 2: " << Ca40_SRC_2 << endl;
	cout << "Ca48 SRC 2: " << Ca48_SRC_2 << endl;
	cout << "Fe54 SRC 2: " << Fe54_SRC_2 << endl;
	cout << "Au197 SRC 2: " << Au197_SRC_2 << endl;

	cout << endl;

	cout << "Be9 SRC 2: " << Be9_SRC_Pm2 << endl;
	cout << "B10 SRC 2: " << B10_SRC_Pm2 << endl;
	cout << "B11 SRC 2: " << B11_SRC_Pm2 << endl;
	cout << "C12 SRC 2: " << C12_SRC_Pm2 << endl;
	cout << "Ca40 SRC 2: " << Ca40_SRC_Pm2 << endl;
	cout << "Ca48 SRC 2: " << Ca48_SRC_Pm2 << endl;
	cout << "Fe54 SRC 2: " << Fe54_SRC_Pm2 << endl;
	cout << "Au197 SRC 2: " << Au197_SRC_Pm2 << endl;

	cout << endl;
	
	cout << "Be9 SRC 3: " << Be9_SRC_3 << endl;
	cout << "B10 SRC 3: " << B10_SRC_3 << endl;
	cout << "B11 SRC 3: " << B11_SRC_3 << endl;
	cout << "C12 SRC 3: " << C12_SRC_3 << endl;
	cout << "Ca40 SRC 3: " << Ca40_SRC_3 << endl;
	cout << "Ca48 SRC 3: " << Ca48_SRC_3 << endl;
	cout << "Fe54 SRC 3: " << Fe54_SRC_3 << endl;
	cout << "Au197 SRC 3: " << Au197_SRC_3 << endl;

	cout << endl;

	cout << "Be9 SRC 3: " << Be9_SRC_Pm3 << endl;
	cout << "B10 SRC 3: " << B10_SRC_Pm3 << endl;
	cout << "B11 SRC 3: " << B11_SRC_Pm3 << endl;
	cout << "C12 SRC 3: " << C12_SRC_Pm3 << endl;
	cout << "Ca40 SRC 3: " << Ca40_SRC_Pm3 << endl;
	cout << "Ca48 SRC 3: " << Ca48_SRC_Pm3 << endl;
	cout << "Fe54 SRC 3: " << Fe54_SRC_Pm3 << endl;
	cout << "Au197 SRC 3: " << Au197_SRC_Pm3 << endl;

	cout << endl;
	*/







	
	ofstream ofile;
	
	Double_t z_cont, z_cont_err;

	Int_t xb, yb, ib;

	Double_t x0, xlow, xup;
	Double_t y0, ylow, yup;
	
	Int_t Pm_b_mf;
	Int_t Nbins_Pm_mf = H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Pmbin_width_mf = H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	//Int_t Nbins_Pm_mf = H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetNbins();
	//Int_t Nbins_Pm_mf = H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetNbins();
	
	Int_t Q2_b_mf;
	Int_t Nbins_Q2_mf = H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Q2bin_width_mf = H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	//Int_t Nbins_Q2_mf = H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetNbins();
	//Int_t Nbins_Q2_mf = H1_Q2_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetNbins();
	
	Int_t Em_b_mf;
	Int_t Nbins_Em_mf = H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Embin_width_mf = H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	//Int_t Nbins_Em_mf = H1_Em_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetNbins();
	//Int_t Nbins_Em_mf = H1_Em_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetNbins();

	Int_t epctime_b_mf;
	Int_t Nbins_epctime_mf = H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetNbins();
	Double_t epctimebin_width_mf = H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);

	Int_t pCalEtotTrkNorm_b_mf;
	Int_t Nbins_pCalEtotTrkNorm_mf = H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetNbins();
	Double_t pCalEtotTrkNormbin_width_mf = H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);

	Int_t hdelta_b_mf;
	Int_t Nbins_hdelta_mf = H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetNbins();
	Double_t hdeltabin_width_mf = H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);

	Int_t edelta_b_mf;
	Int_t Nbins_edelta_mf = H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetNbins();
	Double_t edeltabin_width_mf = H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);

	Int_t Ef_b_mf;
	Int_t Nbins_Ef_mf = H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetNbins();
	Double_t Efbin_width_mf = H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);

	Int_t Pf_b_mf;
	Int_t Nbins_Pf_mf = H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetNbins();
	Double_t Pfbin_width_mf = H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);



	Int_t Pm_b;
	Int_t Nbins_Pm = H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Pmbin_width = H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	//Int_t Nbins_Pm = H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetNbins();

	Int_t Em_b;
	Int_t Nbins_Em = H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Embin_width = H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	//Int_t Nbins_Pm = H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetNbins();
	//Int_t Nbins_Pm = H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetNbins();
	//Int_t Nbins_Pm = H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetNbins();

	Int_t Q2_b;
	Int_t Nbins_Q2 = H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Q2bin_width = H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	//Int_t Nbins_Q2 = H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetNbins();
	//Int_t Nbins_Q2 = H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetNbins();

	Int_t xbj_b;
	Int_t Nbins_xbj = H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t xbjbin_width = H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	//Int_t Nbins_xbj = H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetNbins();
	//Int_t Nbins_xbj = H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetNbins();
	
	Int_t thrq_b;
	Int_t Nbins_thrq = H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t thrqbin_width = H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	//Int_t Nbins_thrq = H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetNbins();
	//Int_t Nbins_thrq = H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetNbins();

	//Int_t Em_b;
	//Int_t Nbins_Em = H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	//Double_t Embin_width = H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t epctime_b;
	Int_t Nbins_epctime = H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t epctimebin_width = H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t pCalEtotTrkNorm_b;
	Int_t Nbins_pCalEtotTrkNorm = H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t pCalEtotTrkNormbin_width = H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t hdelta_b;
	Int_t Nbins_hdelta = H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t hdeltabin_width = H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t edelta_b;
	Int_t Nbins_edelta = H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t edeltabin_width = H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t Ef_b;
	Int_t Nbins_Ef = H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Efbin_width = H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t Pf_b;
	Int_t Nbins_Pf = H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Pfbin_width = H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t the_b;
	Int_t Nbins_the = H1_the_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t thebin_width = H1_the_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t thp_b;
	Int_t Nbins_thp = H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t thpbin_width = H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t nu_b;
	Int_t Nbins_nu = H1_nu_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t nubin_width = H1_nu_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t q_b;
	Int_t Nbins_q = H1_q_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t qbin_width = H1_q_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	Int_t W_b;
	Int_t Nbins_W = H1_W_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Wbin_width = H1_W_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	Int_t Nbins_W_mf = H1_W_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Wbin_width_mf = H1_W_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	Int_t Nbins_Wee = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Wbin_widthee = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	Int_t Nbins_Weep = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetXaxis()->GetNbins();
	Double_t Wbin_widtheep = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);

	/*
	file_SRC_Fe54->GetObject(base12.Data(), H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base13.Data(), H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base14.Data(), H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base15.Data(), H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full);
	file_SRC_Fe54->GetObject(base8.Data(), H1_W_ratio_SRC_Fe54_pid_acc_kin_full);
	*/
	
	TString SRCBe9Pm = "SRCBe9Pm.csv";
	TString SRCBe9thrq = "SRCBe9thrq.csv";
	TString SRCBe9Q2 = "SRCBe9Q2.csv";
	TString SRCBe9xbj = "SRCBe9xbj.csv"; 
	TString SRCBe9Em = "SRCBe9Em.csv";
	TString SRCBe9epctime = "SRCBe9epctime.csv";
	TString SRCBe9epctimea = "SRCBe9epctimea.csv";
	TString SRCBe9pCalEtotTrkNorm = "SRCBe9pCalEtotTrkNorm.csv";
	TString SRCBe9hdelta = "SRCBe9hdelta.csv";
	TString SRCBe9edelta = "SRCBe9edelta.csv";
	TString SRCBe9W = "SRCBe9W.csv";
	TString SRCBe9Ef = "SRCBe9Ef.csv";
	TString SRCBe9Pf = "SRCBe9Pf.csv";
	TString SRCBe9the = "SRCBe9the.csv";
	TString SRCBe9thp = "SRCBe9thp.csv";
	TString SRCBe9nu = "SRCBe9nu.csv";
	TString SRCBe9q = "SRCBe9q.csv";

	TString SRCB10Pm = "SRCB10Pm.csv";
	TString SRCB10thrq = "SRCB10thrq.csv";
	TString SRCB10Q2 = "SRCB10Q2.csv";
	TString SRCB10xbj = "SRCB10xbj.csv"; 
	TString SRCB10Em = "SRCB10Em.csv";
	TString SRCB10epctime = "SRCB10epctime.csv";
	TString SRCB10epctimea = "SRCB10epctimea.csv";
	TString SRCB10pCalEtotTrkNorm = "SRCB10pCalEtotTrkNorm.csv";
	TString SRCB10hdelta = "SRCB10hdelta.csv";
	TString SRCB10edelta = "SRCB10edelta.csv";
	TString SRCB10W = "SRCB10W.csv";
	TString SRCB10Ef = "SRCB10Ef.csv";
	TString SRCB10Pf = "SRCB10Pf.csv";
	TString SRCB10the = "SRCB10the.csv";
	TString SRCB10thp = "SRCB10thp.csv";
	TString SRCB10nu = "SRCB10nu.csv";
	TString SRCB10q = "SRCB10q.csv";

	TString SRCB11Pm = "SRCB11Pm.csv";
	TString SRCB11thrq = "SRCB11thrq.csv";
	TString SRCB11Q2 = "SRCB11Q2.csv";
	TString SRCB11xbj = "SRCB11xbj.csv"; 
	TString SRCB11Em = "SRCB11Em.csv";
	TString SRCB11epctime = "SRCB11epctime.csv";
	TString SRCB11epctimea = "SRCB11epctimea.csv";
	TString SRCB11pCalEtotTrkNorm = "SRCB11pCalEtotTrkNorm.csv";
	TString SRCB11hdelta = "SRCB11hdelta.csv";
	TString SRCB11edelta = "SRCB11edelta.csv";
	TString SRCB11W = "SRCB11W.csv";
	TString SRCB11Ef = "SRCB11Ef.csv";
	TString SRCB11Pf = "SRCB11Pf.csv";
	TString SRCB11the = "SRCB11the.csv";
	TString SRCB11thp = "SRCB11thp.csv";
	TString SRCB11nu = "SRCB11nu.csv";
	TString SRCB11q = "SRCB11q.csv";

	TString SRCC12Pm = "SRCC12Pm.csv";
	TString SRCC12thrq = "SRCC12thrq.csv";
	TString SRCC12Q2 = "SRCC12Q2.csv";
	TString SRCC12xbj = "SRCC12xbj.csv"; 
	TString SRCC12Em = "SRCC12Em.csv";
	TString SRCC12epctime = "SRCC12epctime.csv";
	TString SRCC12epctimea = "SRCC12epctimea.csv";
	TString SRCC12pCalEtotTrkNorm = "SRCC12pCalEtotTrkNorm.csv";
	TString SRCC12hdelta = "SRCC12hdelta.csv";
	TString SRCC12edelta = "SRCC12edelta.csv";
	TString SRCC12W = "SRCC12W.csv";
	TString SRCC12Ef = "SRCC12Ef.csv";
	TString SRCC12Pf = "SRCC12Pf.csv";
	TString SRCC12the = "SRCC12the.csv";
	TString SRCC12thp = "SRCC12thp.csv";
	TString SRCC12nu = "SRCC12nu.csv";
	TString SRCC12q = "SRCC12q.csv";

	TString SRCCa40Pm = "SRCCa40Pm.csv";
	TString SRCCa40thrq = "SRCCa40thrq.csv";
	TString SRCCa40Q2 = "SRCCa40Q2.csv";
	TString SRCCa40xbj = "SRCCa40xbj.csv"; 
	TString SRCCa40Em = "SRCCa40Em.csv";
	TString SRCCa40epctime = "SRCCa40epctime.csv";
	TString SRCCa40epctimea = "SRCCa40epctimea.csv";
	TString SRCCa40pCalEtotTrkNorm = "SRCCa40pCalEtotTrkNorm.csv";
	TString SRCCa40hdelta = "SRCCa40hdelta.csv";
	TString SRCCa40edelta = "SRCCa40edelta.csv";
	TString SRCCa40W = "SRCCa40W.csv";
	TString SRCCa40Ef = "SRCCa40Ef.csv";
	TString SRCCa40Pf = "SRCCa40Pf.csv";
	TString SRCCa40the = "SRCCa40the.csv";
	TString SRCCa40thp = "SRCCa40thp.csv";
	TString SRCCa40nu = "SRCCa40nu.csv";
	TString SRCCa40q = "SRCCa40q.csv";

	TString SRCCa48Pm = "SRCCa48Pm.csv";
	TString SRCCa48thrq = "SRCCa48thrq.csv";
	TString SRCCa48Q2 = "SRCCa48Q2.csv";
	TString SRCCa48xbj = "SRCCa48xbj.csv";
	TString SRCCa48Em = "SRCCa48Em.csv";
	TString SRCCa48epctime = "SRCCa48epctime.csv";
	TString SRCCa48epctimea = "SRCCa48epctimea.csv";
	TString SRCCa48pCalEtotTrkNorm = "SRCCa48pCalEtotTrkNorm.csv";
	TString SRCCa48hdelta = "SRCCa48hdelta.csv";
	TString SRCCa48edelta = "SRCCa48edelta.csv";
	TString SRCCa48W = "SRCCa48W.csv";
	TString SRCCa48Ef = "SRCCa48Ef.csv";
	TString SRCCa48Pf = "SRCCa48Pf.csv";
	TString SRCCa48the = "SRCCa48the.csv";
	TString SRCCa48thp = "SRCCa48thp.csv";
	TString SRCCa48nu = "SRCCa48nu.csv";
	TString SRCCa48q = "SRCCa48q.csv";

	TString SRCFe54Pm = "SRCFe54Pm.csv";
	TString SRCFe54thrq = "SRCFe54thrq.csv";
	TString SRCFe54Q2 = "SRCFe54Q2.csv";
	TString SRCFe54xbj = "SRCFe54xbj.csv";
	TString SRCFe54Em = "SRCFe54Em.csv";
	TString SRCFe54epctime = "SRCFe54epctime.csv";
	TString SRCFe54epctimea = "SRCFe54epctimea.csv";
	TString SRCFe54pCalEtotTrkNorm = "SRCFe54pCalEtotTrkNorm.csv";
	TString SRCFe54hdelta = "SRCFe54hdelta.csv";
	TString SRCFe54edelta = "SRCFe54edelta.csv";
	TString SRCFe54W = "SRCFe54W.csv";
	TString SRCFe54Ef = "SRCFe54Ef.csv";
	TString SRCFe54Pf = "SRCFe54Pf.csv";
	TString SRCFe54the = "SRCFe54the.csv";
	TString SRCFe54thp = "SRCFe54thp.csv";
	TString SRCFe54nu = "SRCFe54nu.csv";
	TString SRCFe54q = "SRCFe54q.csv";

	TString SRCAu197Pm = "SRCAu197Pm.csv";
	TString SRCAu197thrq = "SRCAu197thrq.csv";
	TString SRCAu197Q2 = "SRCAu197Q2.csv";
	TString SRCAu197xbj = "SRCAu197xbj.csv";
	TString SRCAu197Em = "SRCAu197Em.csv";
	TString SRCAu197epctime = "SRCAu197epctime.csv";
	TString SRCAu197epctimea = "SRCAu197epctimea.csv";
	TString SRCAu197pCalEtotTrkNorm = "SRCAu197pCalEtotTrkNorm.csv";
	TString SRCAu197hdelta = "SRCAu197hdelta.csv";
	TString SRCAu197edelta = "SRCAu197edelta.csv";
	TString SRCAu197W = "SRCAu197W.csv";
	TString SRCAu197Ef = "SRCAu197Ef.csv";
	TString SRCAu197Pf = "SRCAu197Pf.csv";
	TString SRCAu197the = "SRCAu197the.csv";
	TString SRCAu197thp = "SRCAu197thp.csv";
	TString SRCAu197nu = "SRCAu197nu.csv";
	TString SRCAu197q = "SRCAu197q.csv";


	TString MFBe9Pm = "MFBe9Pm.csv";
	TString MFBe9W = "MFBe9W.csv";
	TString MFBe9Q2 = "MFBe9Q2.csv";
	TString MFBe9Em = "MFBe9Em.csv"; 
	TString MFBe9epctime = "MFBe9epctime.csv";
	TString MFBe9epctimea = "MFBe9epctimea.csv";
	TString MFBe9pCalEtotTrkNorm = "MFBe9pCalEtotTrkNorm.csv";
	TString MFBe9hdelta = "MFBe9hdelta.csv";
	TString MFBe9edelta = "MFBe9edelta.csv";
	TString MFBe9Ef = "MFBe9Ef.csv";
	TString MFBe9Pf = "MFBe9Pf.csv";
	TString MFBe9the = "MFBe9the.csv";
	TString MFBe9thp = "MFBe9thp.csv";
	TString MFBe9nu = "MFBe9nu.csv";
	TString MFBe9q = "MFBe9q.csv";

	TString MFB10Pm = "MFB10Pm.csv";
	TString MFB10W = "MFB10W.csv";
	TString MFB10Q2 = "MFB10Q2.csv";
	TString MFB10Em = "MFB10Em.csv"; 
	TString MFB10epctime = "MFB10epctime.csv";
	TString MFB10epctimea = "MFB10epctimea.csv";
	TString MFB10pCalEtotTrkNorm = "MFB10pCalEtotTrkNorm.csv";
	TString MFB10hdelta = "MFB10hdelta.csv";
	TString MFB10edelta = "MFB10edelta.csv";
	TString MFB10Ef = "MFB10Ef.csv";
	TString MFB10Pf = "MFB10Pf.csv";
	TString MFB10the = "MFB10the.csv";
	TString MFB10thp = "MFB10thp.csv";
	TString MFB10nu = "MFB10nu.csv";
	TString MFB10q = "MFB10q.csv";

	TString MFB11Pm = "MFB11Pm.csv";
	TString MFB11W = "MFB11W.csv";
	TString MFB11Q2 = "MFB11Q2.csv";
	TString MFB11Em = "MFB11Em.csv"; 
	TString MFB11epctime = "MFB11epctime.csv";
	TString MFB11epctimea = "MFB11epctimea.csv";
	TString MFB11pCalEtotTrkNorm = "MFB11pCalEtotTrkNorm.csv";
	TString MFB11hdelta = "MFB11hdelta.csv";
	TString MFB11edelta = "MFB11edelta.csv";
	TString MFB11Ef = "MFB11Ef.csv";
	TString MFB11Pf = "MFB11Pf.csv";
	TString MFB11the = "MFB11the.csv";
	TString MFB11thp = "MFB11thp.csv";
	TString MFB11nu = "MFB11nu.csv";
	TString MFB11q = "MFB11q.csv";

	TString MFC12Pm = "MFC12Pm.csv";
	TString MFC12W = "MFC12W.csv";
	TString MFC12Q2 = "MFC12Q2.csv";
	TString MFC12Em = "MFC12Em.csv"; 
	TString MFC12epctime = "MFC12epctime.csv";
	TString MFC12epctimea = "MFC12epctimea.csv";
	TString MFC12pCalEtotTrkNorm = "MFC12pCalEtotTrkNorm.csv";
	TString MFC12hdelta = "MFC12hdelta.csv";
	TString MFC12edelta = "MFC12edelta.csv";
	TString MFC12Ef = "MFC12Ef.csv";
	TString MFC12Pf = "MFC12Pf.csv";
	TString MFC12the = "MFC12the.csv";
	TString MFC12thp = "MFC12thp.csv";
	TString MFC12nu = "MFC12nu.csv";
	TString MFC12q = "MFC12q.csv";

	TString MFCa40Pm = "MFCa40Pm.csv";
	TString MFCa40W = "MFCa40W.csv";
	TString MFCa40Q2 = "MFCa40Q2.csv";
	TString MFCa40Em = "MFCa40Em.csv";
	TString MFCa40epctime = "MFCa40epctime.csv";
	TString MFCa40epctimea = "MFCa40epctimea.csv";
	TString MFCa40pCalEtotTrkNorm = "MFCa40pCalEtotTrkNorm.csv";
	TString MFCa40hdelta = "MFCa40hdelta.csv";
	TString MFCa40edelta = "MFCa40edelta.csv";
	TString MFCa40Ef = "MFCa40Ef.csv";
	TString MFCa40Pf = "MFCa40Pf.csv";
	TString MFCa40the = "MFCa40the.csv";
	TString MFCa40thp = "MFCa40thp.csv";
	TString MFCa40nu = "MFCa40nu.csv";
	TString MFCa40q = "MFCa40q.csv";

	TString MFCa48Pm = "MFCa48Pm.csv";
	TString MFCa48W = "MFCa48W.csv";
	TString MFCa48Q2 = "MFCa48Q2.csv";
	TString MFCa48Em = "MFCa48Em.csv";
	TString MFCa48epctime = "MFCa48epctime.csv";
	TString MFCa48epctimea = "MFCa48epctimea.csv";
	TString MFCa48pCalEtotTrkNorm = "MFCa48pCalEtotTrkNorm.csv";
	TString MFCa48hdelta = "MFCa48hdelta.csv";
	TString MFCa48edelta = "MFCa48edelta.csv";
	TString MFCa48Ef = "MFCa48Ef.csv";
	TString MFCa48Pf = "MFCa48Pf.csv";
	TString MFCa48the = "MFCa48the.csv";
	TString MFCa48thp = "MFCa48thp.csv";
	TString MFCa48nu = "MFCa48nu.csv";
	TString MFCa48q = "MFCa48q.csv";

	TString MFFe54Pm = "MFFe54Pm.csv";
	TString MFFe54W = "MFFe54W.csv";
	TString MFFe54Q2 = "MFFe54Q2.csv";
	TString MFFe54Em = "MFFe54Em.csv";
	TString MFFe54epctime = "MFFe54epctime.csv";
	TString MFFe54epctimea = "MFFe54epctimea.csv";
	TString MFFe54pCalEtotTrkNorm = "MFFe54pCalEtotTrkNorm.csv";
	TString MFFe54hdelta = "MFFe54hdelta.csv";
	TString MFFe54edelta = "MFFe54edelta.csv";
	TString MFFe54Ef = "MFFe54Ef.csv";
	TString MFFe54Pf = "MFFe54Pf.csv";
	TString MFFe54the = "MFFe54the.csv";
	TString MFFe54thp = "MFFe54thp.csv";
	TString MFFe54nu = "MFFe54nu.csv";
	TString MFFe54q = "MFFe54q.csv";

	TString MFAu197Pm = "MFAu197Pm.csv";
	TString MFAu197W = "MFAu197W.csv";
	TString MFAu197Q2 = "MFAu197Q2.csv";
	TString MFAu197Em = "MFAu197Em.csv";
	TString MFAu197epctime = "MFAu197epctime.csv";
	TString MFAu197epctimea = "MFAu197epctimea.csv";
	TString MFAu197pCalEtotTrkNorm = "MFAu197pCalEtotTrkNorm.csv";
	TString MFAu197hdelta = "MFAu197hdelta.csv";
	TString MFAu197edelta = "MFAu197edelta.csv";
	TString MFAu197Ef = "MFAu197Ef.csv";
	TString MFAu197Pf = "MFAu197Pf.csv";
	TString MFAu197the = "MFAu197the.csv";
	TString MFAu197thp = "MFAu197thp.csv";
	TString MFAu197nu = "MFAu197nu.csv";
	TString MFAu197q = "MFAu197q.csv";
	
	TString eeD2 = "eeD2.csv";
	TString eepD2 = "eepD2.csv";



	TString header = Form
	(
	"# x0:         x-axis central bin value \n"
	"# zcont:      bin content (z-axis) \n"
	"# zcont_err:  bin content error (z-axis) \n"
	"#                                        \n"
	"x0,zcont,zcont_err"
	);

	


	//MFBe9edelta
	ofile.open(MFBe9Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b_mf = 1; Pf_b_mf <= Nbins_Pf_mf; Pf_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->GetBin(Pf_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinContent(Pf_b_mf);
		z_cont_err = H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinError(Pf_b_mf);

		xlow = H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Pf_b_mf);
		xup = H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Pf_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();



	//MFB10edelta
	ofile.open(MFB10Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b_mf = 1; Pf_b_mf <= Nbins_Pf_mf; Pf_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->GetBin(Pf_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->GetBinContent(Pf_b_mf);
		z_cont_err = H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->GetBinError(Pf_b_mf);

		xlow = H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Pf_b_mf);
		xup = H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Pf_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB11edelta
	ofile.open(MFB11Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b_mf = 1; Pf_b_mf <= Nbins_Pf_mf; Pf_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->GetBin(Pf_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->GetBinContent(Pf_b_mf);
		z_cont_err = H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->GetBinError(Pf_b_mf);

		xlow = H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Pf_b_mf);
		xup = H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Pf_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFC12edelta
	ofile.open(MFC12Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b_mf = 1; Pf_b_mf <= Nbins_Pf_mf; Pf_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(Pf_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(Pf_b_mf);
		z_cont_err = H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(Pf_b_mf);

		xlow = H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Pf_b_mf);
		xup = H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Pf_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa40edelta
	ofile.open(MFCa40Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b_mf = 1; Pf_b_mf <= Nbins_Pf_mf; Pf_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBin(Pf_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinContent(Pf_b_mf);
		z_cont_err = H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinError(Pf_b_mf);

		xlow = H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Pf_b_mf);
		xup = H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Pf_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa48edelta
	ofile.open(MFCa48Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b_mf = 1; Pf_b_mf <= Nbins_Pf_mf; Pf_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBin(Pf_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinContent(Pf_b_mf);
		z_cont_err = H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinError(Pf_b_mf);

		xlow = H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Pf_b_mf);
		xup = H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Pf_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFFe54edelta
	ofile.open(MFFe54Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b_mf = 1; Pf_b_mf <= Nbins_Pf_mf; Pf_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBin(Pf_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinContent(Pf_b_mf);
		z_cont_err = H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinError(Pf_b_mf);

		xlow = H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Pf_b_mf);
		xup = H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Pf_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFAu197edelta
	ofile.open(MFAu197Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b_mf = 1; Pf_b_mf <= Nbins_Pf_mf; Pf_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->GetBin(Pf_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinContent(Pf_b_mf);
		z_cont_err = H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinError(Pf_b_mf);

		xlow = H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Pf_b_mf);
		xup = H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Pf_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();






	//SRCBe9Pm
	ofile.open(SRCBe9Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_SRC_Be9_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_pmiss_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10Pm
	ofile.open(SRCB10Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_SRC_B10_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_SRC_B10_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_pmiss_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11Pm
	ofile.open(SRCB11Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_SRC_B11_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_SRC_B11_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_pmiss_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	

	//SRCC12Pm
	ofile.open(SRCC12Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for(int Pm_b=1; Pm_b<=Nbins_Pm; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_SRC_C12_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_SRC_C12_pid_acc_kin_full->GetBinError(Pm_b);
		
		xlow = H1_pmiss_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup  = H1_pmiss_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0   = (xlow + xup)/2.; //H1_pmiss_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48Pm
	ofile.open(SRCCa40Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48Pm
	ofile.open(SRCCa48Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54Pm
	ofile.open(SRCFe54Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197Pm
	ofile.open(SRCAu197Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_SRC_Au197_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();












	//SRCBe9Pm
	ofile.open(SRCBe9Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_SRC_Be9_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Em_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10Em
	ofile.open(SRCB10Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_SRC_B10_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_SRC_B10_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Em_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11Em
	ofile.open(SRCB11Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_SRC_B11_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_SRC_B11_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Em_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12Em
	ofile.open(SRCC12Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_SRC_C12_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_SRC_C12_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Em_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48Em
	ofile.open(SRCCa40Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48Em
	ofile.open(SRCCa48Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54Em
	ofile.open(SRCFe54Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197Em
	ofile.open(SRCAu197Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_SRC_Au197_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();












	//SRCBe9Q2
	ofile.open(SRCBe9Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_SRC_Be9_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10Q2
	ofile.open(SRCB10Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_SRC_B10_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_SRC_B10_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11Q2
	ofile.open(SRCB11Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_SRC_B11_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_SRC_B11_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12Q2
	ofile.open(SRCC12Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_SRC_C12_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_SRC_C12_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40Q2
	ofile.open(SRCCa40Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48Q2
	ofile.open(SRCCa48Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54Q2
	ofile.open(SRCFe54Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thrq_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197Q2
	ofile.open(SRCAu197Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_SRC_Au197_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thrq_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	







	//SRCBe9thrq
	ofile.open(SRCBe9thrq.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thrq_b = 1; thrq_b <= Nbins_thrq; thrq_b++)
	{
		//Get 2d bin number
		ib = H1_thrq_ratio_SRC_Be9_pid_acc_kin_full->GetBin(thrq_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thrq_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(thrq_b);
		z_cont_err = H1_thrq_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(thrq_b);

		xlow = H1_thrq_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thrq_b);
		xup = H1_thrq_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thrq_b);
		x0 = (xlow + xup) / 2.; //H1_thrq_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(thrq_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10thrq
	ofile.open(SRCB10thrq.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thrq_b = 1; thrq_b <= Nbins_thrq; thrq_b++)
	{
		//Get 2d bin number
		ib = H1_thrq_ratio_SRC_B10_pid_acc_kin_full->GetBin(thrq_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thrq_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(thrq_b);
		z_cont_err = H1_thrq_ratio_SRC_B10_pid_acc_kin_full->GetBinError(thrq_b);

		xlow = H1_thrq_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thrq_b);
		xup = H1_thrq_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thrq_b);
		x0 = (xlow + xup) / 2.; //H1_thrq_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(thrq_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11thrq
	ofile.open(SRCB11thrq.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thrq_b = 1; thrq_b <= Nbins_thrq; thrq_b++)
	{
		//Get 2d bin number
		ib = H1_thrq_ratio_SRC_B11_pid_acc_kin_full->GetBin(thrq_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thrq_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(thrq_b);
		z_cont_err = H1_thrq_ratio_SRC_B11_pid_acc_kin_full->GetBinError(thrq_b);

		xlow = H1_thrq_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thrq_b);
		xup = H1_thrq_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thrq_b);
		x0 = (xlow + xup) / 2.; //H1_thrq_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(thrq_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12thrq
	ofile.open(SRCC12thrq.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thrq_b = 1; thrq_b <= Nbins_thrq; thrq_b++)
	{
		//Get 2d bin number
		ib = H1_thrq_ratio_SRC_C12_pid_acc_kin_full->GetBin(thrq_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thrq_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(thrq_b);
		z_cont_err = H1_thrq_ratio_SRC_C12_pid_acc_kin_full->GetBinError(thrq_b);

		xlow = H1_thrq_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thrq_b);
		xup = H1_thrq_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thrq_b);
		x0 = (xlow + xup) / 2.; //H1_thrq_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(thrq_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40thrq
	ofile.open(SRCCa40thrq.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thrq_b = 1; thrq_b <= Nbins_thrq; thrq_b++)
	{
		//Get 2d bin number
		ib = H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(thrq_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(thrq_b);
		z_cont_err = H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(thrq_b);

		xlow = H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thrq_b);
		xup = H1_thrq_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thrq_b);
		x0 = (xlow + xup) / 2.; //H1_thrqiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thrq_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48thrq
	ofile.open(SRCCa48thrq.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thrq_b = 1; thrq_b <= Nbins_thrq; thrq_b++)
	{
		//Get 2d bin number
		ib = H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(thrq_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(thrq_b);
		z_cont_err = H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(thrq_b);

		xlow = H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thrq_b);
		xup = H1_thrq_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thrq_b);
		x0 = (xlow + xup) / 2.; //H1_thrqiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thrq_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54thrq
	ofile.open(SRCFe54thrq.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thrq_b = 1; thrq_b <= Nbins_thrq; thrq_b++)
	{
		//Get 2d bin number
		ib = H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(thrq_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(thrq_b);
		z_cont_err = H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(thrq_b);

		xlow = H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thrq_b);
		xup = H1_thrq_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thrq_b);
		x0 = (xlow + xup) / 2.; //H1_thrqiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197thrq
	ofile.open(SRCAu197thrq.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thrq_b = 1; thrq_b <= Nbins_thrq; thrq_b++)
	{
		//Get 2d bin number
		ib = H1_thrq_ratio_SRC_Au197_pid_acc_kin_full->GetBin(thrq_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thrq_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(thrq_b);
		z_cont_err = H1_thrq_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(thrq_b);

		xlow = H1_thrq_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thrq_b);
		xup = H1_thrq_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thrq_b);
		x0 = (xlow + xup) / 2.; //H1_thrqiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();








	//SRCBe9xbj
	ofile.open(SRCBe9xbj.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int xbj_b = 1; xbj_b <= Nbins_xbj; xbj_b++)
	{
		//Get 2d bin number
		ib = H1_xbj_ratio_SRC_Be9_pid_acc_kin_full->GetBin(xbj_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_xbj_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(xbj_b);
		z_cont_err = H1_xbj_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(xbj_b);

		xlow = H1_xbj_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(xbj_b);
		xup = H1_xbj_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(xbj_b);
		x0 = (xlow + xup) / 2.; //H1_xbj_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10xbj
	ofile.open(SRCB10xbj.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int xbj_b = 1; xbj_b <= Nbins_xbj; xbj_b++)
	{
		//Get 2d bin number
		ib = H1_xbj_ratio_SRC_B10_pid_acc_kin_full->GetBin(xbj_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_xbj_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(xbj_b);
		z_cont_err = H1_xbj_ratio_SRC_B10_pid_acc_kin_full->GetBinError(xbj_b);

		xlow = H1_xbj_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(xbj_b);
		xup = H1_xbj_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(xbj_b);
		x0 = (xlow + xup) / 2.; //H1_xbj_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11xbj
	ofile.open(SRCB11xbj.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int xbj_b = 1; xbj_b <= Nbins_xbj; xbj_b++)
	{
		//Get 2d bin number
		ib = H1_xbj_ratio_SRC_B11_pid_acc_kin_full->GetBin(xbj_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_xbj_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(xbj_b);
		z_cont_err = H1_xbj_ratio_SRC_B11_pid_acc_kin_full->GetBinError(xbj_b);

		xlow = H1_xbj_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(xbj_b);
		xup = H1_xbj_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(xbj_b);
		x0 = (xlow + xup) / 2.; //H1_xbj_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12xbj
	ofile.open(SRCC12xbj.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int xbj_b = 1; xbj_b <= Nbins_xbj; xbj_b++)
	{
		//Get 2d bin number
		ib = H1_xbj_ratio_SRC_C12_pid_acc_kin_full->GetBin(xbj_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_xbj_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(xbj_b);
		z_cont_err = H1_xbj_ratio_SRC_C12_pid_acc_kin_full->GetBinError(xbj_b);

		xlow = H1_xbj_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(xbj_b);
		xup = H1_xbj_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(xbj_b);
		x0 = (xlow + xup) / 2.; //H1_xbj_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40xbj
	ofile.open(SRCCa40xbj.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int xbj_b = 1; xbj_b <= Nbins_xbj; xbj_b++)
	{
		//Get 2d bin number
		ib = H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(xbj_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(xbj_b);
		z_cont_err = H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(xbj_b);

		xlow = H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(xbj_b);
		xup = H1_xbj_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(xbj_b);
		x0 = (xlow + xup) / 2.; //H1_xbjiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48xbj
	ofile.open(SRCCa48xbj.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int xbj_b = 1; xbj_b <= Nbins_xbj; xbj_b++)
	{
		//Get 2d bin number
		ib = H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(xbj_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(xbj_b);
		z_cont_err = H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(xbj_b);

		xlow = H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(xbj_b);
		xup = H1_xbj_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(xbj_b);
		x0 = (xlow + xup) / 2.; //H1_xbjiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54xbj
	ofile.open(SRCFe54xbj.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int xbj_b = 1; xbj_b <= Nbins_xbj; xbj_b++)
	{
		//Get 2d bin number
		ib = H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(xbj_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(xbj_b);
		z_cont_err = H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(xbj_b);

		xlow = H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(xbj_b);
		xup = H1_xbj_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(xbj_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197xbj
	ofile.open(SRCAu197xbj.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int xbj_b = 1; xbj_b <= Nbins_xbj; xbj_b++)
	{
		//Get 2d bin number
		ib = H1_xbj_ratio_SRC_Au197_pid_acc_kin_full->GetBin(xbj_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_xbj_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(xbj_b);
		z_cont_err = H1_xbj_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(xbj_b);

		xlow = H1_xbj_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(xbj_b);
		xup = H1_xbj_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(xbj_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(xbj_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();





	//SRCBe9epctime
	ofile.open(SRCBe9epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_SRC_Be9_pid_acc_kin_full->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctime_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10epctime
	ofile.open(SRCB10epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_SRC_B10_pid_acc_kin_full->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_SRC_B10_pid_acc_kin_full->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctime_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11epctime
	ofile.open(SRCB11epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_SRC_B11_pid_acc_kin_full->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_SRC_B11_pid_acc_kin_full->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctime_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12epctime
	ofile.open(SRCC12epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_SRC_C12_pid_acc_kin_full->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_SRC_C12_pid_acc_kin_full->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctime_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40epctime
	ofile.open(SRCCa40epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctimeiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48epctime
	ofile.open(SRCCa48epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctimeiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54epctime
	ofile.open(SRCFe54epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197epctime
	ofile.open(SRCAu197epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_SRC_Au197_pid_acc_kin_full->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();





	




	//SRCBe9pCalEtotTrkNorm
	ofile.open(SRCBe9pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_SRC_Be9_pid_acc_kin_full->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10pCalEtotTrkNorm
	ofile.open(SRCB10pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11pCalEtotTrkNorm
	ofile.open(SRCB11pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12pCalEtotTrkNorm
	ofile.open(SRCC12pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40pCalEtotTrkNorm
	ofile.open(SRCCa40pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNormiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48pCalEtotTrkNorm
	ofile.open(SRCCa48pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNormiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54pCalEtotTrkNorm
	ofile.open(SRCFe54pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197pCalEtotTrkNorm
	ofile.open(SRCAu197pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_SRC_Au197_pid_acc_kin_full->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();









	//SRCBe9hdelta
	ofile.open(SRCBe9hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_SRC_Be9_pid_acc_kin_full->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdelta_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10hdelta
	ofile.open(SRCB10hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_SRC_B10_pid_acc_kin_full->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_SRC_B10_pid_acc_kin_full->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdelta_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11hdelta
	ofile.open(SRCB11hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_SRC_B11_pid_acc_kin_full->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_SRC_B11_pid_acc_kin_full->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdelta_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12hdelta
	ofile.open(SRCC12hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_SRC_C12_pid_acc_kin_full->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_SRC_C12_pid_acc_kin_full->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdelta_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40hdelta
	ofile.open(SRCCa40hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdeltaiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48hdelta
	ofile.open(SRCCa48hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdeltaiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54hdelta
	ofile.open(SRCFe54hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197hdelta
	ofile.open(SRCAu197hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_SRC_Au197_pid_acc_kin_full->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();






	


	//SRCBe9edelta
	ofile.open(SRCBe9edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_SRC_Be9_pid_acc_kin_full->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edelta_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10edelta
	ofile.open(SRCB10edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_SRC_B10_pid_acc_kin_full->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_SRC_B10_pid_acc_kin_full->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edelta_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11edelta
	ofile.open(SRCB11edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_SRC_B11_pid_acc_kin_full->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_SRC_B11_pid_acc_kin_full->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edelta_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12edelta
	ofile.open(SRCC12edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_SRC_C12_pid_acc_kin_full->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_SRC_C12_pid_acc_kin_full->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edelta_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40edelta
	ofile.open(SRCCa40edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edeltaiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48edelta
	ofile.open(SRCCa48edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edeltaiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54edelta
	ofile.open(SRCFe54edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197edelta
	ofile.open(SRCAu197edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_SRC_Au197_pid_acc_kin_full->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();














	//SRCBe9W
	ofile.open(SRCBe9W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_Be9_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_W_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10W
	ofile.open(SRCB10W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_B10_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_SRC_B10_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_W_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11W
	ofile.open(SRCB11W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_B11_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_SRC_B11_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_W_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12W
	ofile.open(SRCC12W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_C12_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_SRC_C12_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_W_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40W
	ofile.open(SRCCa40W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_Wiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48W
	ofile.open(SRCCa48W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_Wiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54W
	ofile.open(SRCFe54W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197W
	ofile.open(SRCAu197W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_Au197_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();

	













	//SRCBe9edelta
	ofile.open(SRCBe9nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_SRC_Be9_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nu_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10nu
	ofile.open(SRCB10nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_SRC_B10_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_SRC_B10_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nu_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11nu
	ofile.open(SRCB11nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_SRC_B11_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_SRC_B11_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nu_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12nu
	ofile.open(SRCC12nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_SRC_C12_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_SRC_C12_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nu_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40nu
	ofile.open(SRCCa40nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nuiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48nu
	ofile.open(SRCCa48nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nuiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54nu
	ofile.open(SRCFe54nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197nu
	ofile.open(SRCAu197nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_SRC_Au197_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();



	//SRCBe9edelta
	ofile.open(SRCBe9q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_SRC_Be9_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_q_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10q
	ofile.open(SRCB10q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_SRC_B10_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_SRC_B10_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_q_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11q
	ofile.open(SRCB11q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_SRC_B11_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_SRC_B11_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_q_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12q
	ofile.open(SRCC12q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_SRC_C12_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_SRC_C12_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_q_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40q
	ofile.open(SRCCa40q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_qiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48q
	ofile.open(SRCCa48q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_qiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54q
	ofile.open(SRCFe54q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197q
	ofile.open(SRCAu197q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_SRC_Au197_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();











	//SRCBe9edelta
	ofile.open(SRCBe9Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b = 1; Ef_b <= Nbins_Ef; Ef_b++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_SRC_Be9_pid_acc_kin_full->GetBin(Ef_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(Ef_b);
		z_cont_err = H1_Ef_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(Ef_b);

		xlow = H1_Ef_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Ef_b);
		xup = H1_Ef_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Ef_b);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10Ef
	ofile.open(SRCB10Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b = 1; Ef_b <= Nbins_Ef; Ef_b++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_SRC_B10_pid_acc_kin_full->GetBin(Ef_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(Ef_b);
		z_cont_err = H1_Ef_ratio_SRC_B10_pid_acc_kin_full->GetBinError(Ef_b);

		xlow = H1_Ef_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Ef_b);
		xup = H1_Ef_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Ef_b);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11Ef
	ofile.open(SRCB11Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b = 1; Ef_b <= Nbins_Ef; Ef_b++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_SRC_B11_pid_acc_kin_full->GetBin(Ef_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(Ef_b);
		z_cont_err = H1_Ef_ratio_SRC_B11_pid_acc_kin_full->GetBinError(Ef_b);

		xlow = H1_Ef_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Ef_b);
		xup = H1_Ef_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Ef_b);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12Ef
	ofile.open(SRCC12Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b = 1; Ef_b <= Nbins_Ef; Ef_b++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_SRC_C12_pid_acc_kin_full->GetBin(Ef_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(Ef_b);
		z_cont_err = H1_Ef_ratio_SRC_C12_pid_acc_kin_full->GetBinError(Ef_b);

		xlow = H1_Ef_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Ef_b);
		xup = H1_Ef_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Ef_b);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40Ef
	ofile.open(SRCCa40Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b = 1; Ef_b <= Nbins_Ef; Ef_b++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(Ef_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(Ef_b);
		z_cont_err = H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(Ef_b);

		xlow = H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Ef_b);
		xup = H1_Ef_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Ef_b);
		x0 = (xlow + xup) / 2.; //H1_Efiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48Ef
	ofile.open(SRCCa48Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b = 1; Ef_b <= Nbins_Ef; Ef_b++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(Ef_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(Ef_b);
		z_cont_err = H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(Ef_b);

		xlow = H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Ef_b);
		xup = H1_Ef_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Ef_b);
		x0 = (xlow + xup) / 2.; //H1_Efiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54Ef
	ofile.open(SRCFe54Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b = 1; Ef_b <= Nbins_Ef; Ef_b++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(Ef_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(Ef_b);
		z_cont_err = H1_Ef_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(Ef_b);

		xlow = H1_Ef_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Ef_b);
		xup = H1_Ef_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Ef_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197Ef
	ofile.open(SRCAu197Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b = 1; Ef_b <= Nbins_Ef; Ef_b++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_SRC_Au197_pid_acc_kin_full->GetBin(Ef_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(Ef_b);
		z_cont_err = H1_Ef_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(Ef_b);

		xlow = H1_Ef_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Ef_b);
		xup = H1_Ef_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Ef_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close(); 
	
	
	




	//SRCBe9Pf
	ofile.open(SRCBe9Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b = 1; Pf_b <= Nbins_Pf; Pf_b++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_SRC_Be9_pid_acc_kin_full->GetBin(Pf_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(Pf_b);
		z_cont_err = H1_Pf_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(Pf_b);

		xlow = H1_Pf_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pf_b);
		xup = H1_Pf_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pf_b);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10Pf
	ofile.open(SRCB10Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b = 1; Pf_b <= Nbins_Pf; Pf_b++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_SRC_B10_pid_acc_kin_full->GetBin(Pf_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(Pf_b);
		z_cont_err = H1_Pf_ratio_SRC_B10_pid_acc_kin_full->GetBinError(Pf_b);

		xlow = H1_Pf_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pf_b);
		xup = H1_Pf_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pf_b);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11Pf
	ofile.open(SRCB11Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b = 1; Pf_b <= Nbins_Pf; Pf_b++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_SRC_B11_pid_acc_kin_full->GetBin(Pf_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(Pf_b);
		z_cont_err = H1_Pf_ratio_SRC_B11_pid_acc_kin_full->GetBinError(Pf_b);

		xlow = H1_Pf_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pf_b);
		xup = H1_Pf_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pf_b);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12Pf
	ofile.open(SRCC12Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b = 1; Pf_b <= Nbins_Pf; Pf_b++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_SRC_C12_pid_acc_kin_full->GetBin(Pf_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(Pf_b);
		z_cont_err = H1_Pf_ratio_SRC_C12_pid_acc_kin_full->GetBinError(Pf_b);

		xlow = H1_Pf_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pf_b);
		xup = H1_Pf_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pf_b);
		x0 = (xlow + xup) / 2.; //H1_Pf_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40Pf
	ofile.open(SRCCa40Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b = 1; Pf_b <= Nbins_Pf; Pf_b++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(Pf_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(Pf_b);
		z_cont_err = H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(Pf_b);

		xlow = H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pf_b);
		xup = H1_Pf_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pf_b);
		x0 = (xlow + xup) / 2.; //H1_Pfiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48Pf
	ofile.open(SRCCa48Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b = 1; Pf_b <= Nbins_Pf; Pf_b++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(Pf_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(Pf_b);
		z_cont_err = H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(Pf_b);

		xlow = H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pf_b);
		xup = H1_Pf_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pf_b);
		x0 = (xlow + xup) / 2.; //H1_Pfiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54Pf
	ofile.open(SRCFe54Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b = 1; Pf_b <= Nbins_Pf; Pf_b++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(Pf_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(Pf_b);
		z_cont_err = H1_Pf_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(Pf_b);

		xlow = H1_Pf_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pf_b);
		xup = H1_Pf_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pf_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197Pf
	ofile.open(SRCAu197Pf.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pf_b = 1; Pf_b <= Nbins_Pf; Pf_b++)
	{
		//Get 2d bin number
		ib = H1_Pf_ratio_SRC_Au197_pid_acc_kin_full->GetBin(Pf_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Pf_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(Pf_b);
		z_cont_err = H1_Pf_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(Pf_b);

		xlow = H1_Pf_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pf_b);
		xup = H1_Pf_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pf_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pf_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close(); 
	
	
	






	//SRCBe9the
	ofile.open(SRCBe9the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_SRC_Be9_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_the_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10the
	ofile.open(SRCB10the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_SRC_B10_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_SRC_B10_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_the_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11the
	ofile.open(SRCB11the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_SRC_B11_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_SRC_B11_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_the_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12the
	ofile.open(SRCC12the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_SRC_C12_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_SRC_C12_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_the_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40the
	ofile.open(SRCCa40the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_theiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48the
	ofile.open(SRCCa48the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_theiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54the
	ofile.open(SRCFe54the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197the
	ofile.open(SRCAu197the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_SRC_Au197_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close(); 
	
	
	
	




	//SRCBe9thp
	ofile.open(SRCBe9thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_SRC_Be9_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_SRC_Be9_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_SRC_Be9_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_SRC_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thp_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB10thp
	ofile.open(SRCB10thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_SRC_B10_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_SRC_B10_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_SRC_B10_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_SRC_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thp_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCB11thp
	ofile.open(SRCB11thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_SRC_B11_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_SRC_B11_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_SRC_B11_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_SRC_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thp_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCC12thp
	ofile.open(SRCC12thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_SRC_C12_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_SRC_C12_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thp_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa40thp
	ofile.open(SRCCa40thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_SRC_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thpiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCCa48thp
	ofile.open(SRCCa48thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_SRC_Ca48_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_SRC_Ca48_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_SRC_Ca48_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_SRC_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thpiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCFe54thp
	ofile.open(SRCFe54thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_SRC_Fe54_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_SRC_Fe54_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_SRC_Fe54_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_SRC_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//SRCAu197thp
	ofile.open(SRCAu197thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_SRC_Au197_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_SRC_Au197_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_SRC_Au197_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_SRC_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();













	//MFBe9Q2
	ofile.open(MFBe9Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2_mf; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_MF_Be9_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_MF_Be9_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_MF_Be9_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB10Q2
	ofile.open(MFB10Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2_mf; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_MF_B10_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_MF_B10_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_MF_B10_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB11Q2
	ofile.open(MFB11Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2_mf; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_MF_B11_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_MF_B11_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_MF_B11_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFC12Q2
	ofile.open(MFC12Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2_mf; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_MF_C12_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_MF_C12_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_MF_C12_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa40Q2
	ofile.open(MFCa40Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2_mf; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa48Q2
	ofile.open(MFCa48Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2_mf; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Q2_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFFe54Q2
	ofile.open(MFFe54Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2_mf; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_MF_Fe54_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_MF_Fe54_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_MF_Fe54_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFAu197Q2
	ofile.open(MFAu197Q2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Q2_b = 1; Q2_b <= Nbins_Q2_mf; Q2_b++)
	{
		//Get 2d bin number
		ib = H1_Q2_ratio_MF_Au197_pid_acc_kin_full->GetBin(Q2_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Q2_ratio_MF_Au197_pid_acc_kin_full->GetBinContent(Q2_b);
		z_cont_err = H1_Q2_ratio_MF_Au197_pid_acc_kin_full->GetBinError(Q2_b);

		xlow = H1_Q2_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Q2_b);
		xup = H1_Q2_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Q2_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();









	//MFBe9Pm
	ofile.open(MFBe9Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm_mf; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_MF_Be9_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_MF_Be9_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_MF_Be9_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_Pmiss_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	
	
	//MFB10Pm
	ofile.open(MFB10Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm_mf; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_MF_B10_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_MF_B10_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_MF_B10_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_Pmiss_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	
	
	//MFB11Pm
	ofile.open(MFB11Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm_mf; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_MF_B11_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_MF_B11_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_MF_B11_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_Pmiss_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	
	
	//MFC12Pm
	ofile.open(MFC12Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm_mf; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_MF_C12_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_MF_C12_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_MF_C12_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_Pmiss_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	
	
	//MFCa40Pm
	ofile.open(MFCa40Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm_mf; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_Pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa48Pm
	ofile.open(MFCa48Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm_mf; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_Pmiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFFe54Pm
	ofile.open(MFFe54Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm_mf; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFAu197Pm
	ofile.open(MFAu197Pm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Pm_b = 1; Pm_b <= Nbins_Pm_mf; Pm_b++)
	{
		//Get 2d bin number
		ib = H1_pmiss_ratio_MF_Au197_pid_acc_kin_full->GetBin(Pm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pmiss_ratio_MF_Au197_pid_acc_kin_full->GetBinContent(Pm_b);
		z_cont_err = H1_pmiss_ratio_MF_Au197_pid_acc_kin_full->GetBinError(Pm_b);

		xlow = H1_pmiss_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Pm_b);
		xup = H1_pmiss_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Pm_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();















	//MFBe9Q2
	ofile.open(MFBe9Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em_mf; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_MF_Be9_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_MF_Be9_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_MF_Be9_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Emiss_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	
	
	//MFB10Q2
	ofile.open(MFB10Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em_mf; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_MF_B10_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_MF_B10_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_MF_B10_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Emiss_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	
	
	//MFB11Q2
	ofile.open(MFB11Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em_mf; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_MF_B11_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_MF_B11_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_MF_B11_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Emiss_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	
	
	//MFC12Q2
	ofile.open(MFC12Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em_mf; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_MF_C12_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_MF_C12_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_MF_C12_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Emiss_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
	
	
	//MFCa40Q2
	ofile.open(MFCa40Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em_mf; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Emiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa48Em
	ofile.open(MFCa48Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em_mf; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_MF_Ca48_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_MF_Ca48_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_MF_Ca48_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Emiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFFe54Em
	ofile.open(MFFe54Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em_mf; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_MF_Fe54_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_MF_Fe54_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_MF_Fe54_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFAu197Em
	ofile.open(MFAu197Em.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Em_mf; Em_b++)
	{
		//Get 2d bin number
		ib = H1_Em_ratio_MF_Au197_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_Em_ratio_MF_Au197_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_Em_ratio_MF_Au197_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_Em_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_Em_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Pm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();













	//SRCBe9epctime
	ofile.open( MFBe9epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime_mf; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctime_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB10epctime
	ofile.open( MFB10epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime_mf; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctime_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB11epctime
	ofile.open( MFB11epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime_mf; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctime_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFC12epctime
	ofile.open( MFC12epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime_mf; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctime_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa40epctime
	ofile.open( MFCa40epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime_mf; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctimeiss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa48epctime
	ofile.open( MFCa48epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime_mf; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_epctimeiss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFFe54epctime
	ofile.open( MFFe54epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime_mf; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFAu197epctime
	ofile.open( MFAu197epctime.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int epctime_b = 1; epctime_b <= Nbins_epctime_mf; epctime_b++)
	{
		//Get 2d bin number
		ib = H1_epctime_ratio_MF_Au197_pid_acc_kin_full_mf->GetBin(epctime_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_epctime_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinContent(epctime_b);
		z_cont_err = H1_epctime_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinError(epctime_b);

		xlow = H1_epctime_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(epctime_b);
		xup = H1_epctime_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(epctime_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(epctime_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();












	// MFBe9pCalEtotTrkNorm
	ofile.open( MFBe9pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm_mf; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNorm_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB10pCalEtotTrkNorm
	ofile.open( MFB10pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm_mf; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNorm_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB11pCalEtotTrkNorm
	ofile.open( MFB11pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm_mf; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNorm_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFC12pCalEtotTrkNorm
	ofile.open( MFC12pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm_mf; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNorm_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa40pCalEtotTrkNorm
	ofile.open( MFCa40pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm_mf; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNormiss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa48pCalEtotTrkNorm
	ofile.open( MFCa48pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm_mf; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_pCalEtotTrkNormiss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFFe54pCalEtotTrkNorm
	ofile.open( MFFe54pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm_mf; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFAu197pCalEtotTrkNorm
	ofile.open( MFAu197pCalEtotTrkNorm.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int pCalEtotTrkNorm_b = 1; pCalEtotTrkNorm_b <= Nbins_pCalEtotTrkNorm_mf; pCalEtotTrkNorm_b++)
	{
		//Get 2d bin number
		ib = H1_pCalEtotTrkNorm_ratio_MF_Au197_pid_acc_kin_full_mf->GetBin(pCalEtotTrkNorm_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_pCalEtotTrkNorm_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinContent(pCalEtotTrkNorm_b);
		z_cont_err = H1_pCalEtotTrkNorm_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinError(pCalEtotTrkNorm_b);

		xlow = H1_pCalEtotTrkNorm_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(pCalEtotTrkNorm_b);
		xup = H1_pCalEtotTrkNorm_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(pCalEtotTrkNorm_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(pCalEtotTrkNorm_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();









	// MFBe9hdelta
	ofile.open( MFBe9hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta_mf; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB10hdelta
	ofile.open( MFB10hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta_mf; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdelta_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB11hdelta
	ofile.open( MFB11hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta_mf; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdelta_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFC12hdelta
	ofile.open( MFC12hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta_mf; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdelta_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa40hdelta
	ofile.open( MFCa40hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta_mf; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdeltaiss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa48hdelta
	ofile.open( MFCa48hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta_mf; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_hdeltaiss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFFe54hdelta
	ofile.open( MFFe54hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta_mf; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFAu197hdelta
	ofile.open( MFAu197hdelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int hdelta_b = 1; hdelta_b <= Nbins_hdelta_mf; hdelta_b++)
	{
		//Get 2d bin number
		ib = H1_hdelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetBin(hdelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_hdelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinContent(hdelta_b);
		z_cont_err = H1_hdelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinError(hdelta_b);

		xlow = H1_hdelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(hdelta_b);
		xup = H1_hdelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(hdelta_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(hdelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();









	// MFBe9edelta
	ofile.open( MFBe9edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta_mf; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edelta_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB10edelta
	ofile.open( MFB10edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta_mf; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edelta_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB11edelta
	ofile.open( MFB11edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta_mf; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edelta_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFC12edelta
	ofile.open( MFC12edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta_mf; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edelta_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa40edelta
	ofile.open( MFCa40edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta_mf; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edeltaiss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa48edelta
	ofile.open( MFCa48edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta_mf; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_edeltaiss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFFe54edelta
	ofile.open( MFFe54edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta_mf; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFAu197edelta
	ofile.open( MFAu197edelta.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int edelta_b = 1; edelta_b <= Nbins_edelta_mf; edelta_b++)
	{
		//Get 2d bin number
		ib = H1_edelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetBin(edelta_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_edelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinContent(edelta_b);
		z_cont_err = H1_edelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinError(edelta_b);

		xlow = H1_edelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(edelta_b);
		xup = H1_edelta_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(edelta_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();














	// MFBe9W
	ofile.open( MFBe9W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W_mf; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_MF_Be9_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_MF_Be9_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_MF_Be9_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_W_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB10W
	ofile.open( MFB10W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W_mf; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_MF_B10_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_MF_B10_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_MF_B10_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_W_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFB11W
	ofile.open( MFB11W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W_mf; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_MF_B11_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_MF_B11_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_MF_B11_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_W_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFC12W
	ofile.open( MFC12W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W_mf; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_MF_C12_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_MF_C12_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_MF_C12_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_W_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa40W
	ofile.open( MFCa40W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W_mf; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_MF_Ca40_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_MF_Ca40_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_MF_Ca40_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_Wiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFCa48W
	ofile.open( MFCa48W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W_mf; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_MF_Ca48_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_MF_Ca48_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_MF_Ca48_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_Wiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFFe54W
	ofile.open( MFFe54W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W_mf; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_MF_Fe54_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_MF_Fe54_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_MF_Fe54_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	// MFAu197W
	ofile.open( MFAu197W.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_W_mf; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_MF_Au197_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_MF_Au197_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_MF_Au197_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(W_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();






























	//SRCBe9edelta
	ofile.open(MFBe9nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_MF_Be9_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_MF_Be9_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_MF_Be9_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nu_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB10nu
	ofile.open(MFB10nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_MF_B10_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_MF_B10_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_MF_B10_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nu_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB11nu
	ofile.open(MFB11nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_MF_B11_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_MF_B11_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_MF_B11_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nu_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFC12nu
	ofile.open(MFC12nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_MF_C12_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_MF_C12_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_MF_C12_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nu_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa40nu
	ofile.open(MFCa40nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_MF_Ca40_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_MF_Ca40_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_MF_Ca40_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nuiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa48nu
	ofile.open(MFCa48nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_MF_Ca48_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_MF_Ca48_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_MF_Ca48_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_nuiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFFe54nu
	ofile.open(MFFe54nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_MF_Fe54_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_MF_Fe54_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_MF_Fe54_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFAu197nu
	ofile.open(MFAu197nu.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int nu_b = 1; nu_b <= Nbins_nu; nu_b++)
	{
		//Get 2d bin number
		ib = H1_nu_ratio_MF_Au197_pid_acc_kin_full->GetBin(nu_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_nu_ratio_MF_Au197_pid_acc_kin_full->GetBinContent(nu_b);
		z_cont_err = H1_nu_ratio_MF_Au197_pid_acc_kin_full->GetBinError(nu_b);

		xlow = H1_nu_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(nu_b);
		xup = H1_nu_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(nu_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(nu_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();



	//MFBe9edelta
	ofile.open(MFBe9q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_MF_Be9_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_MF_Be9_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_MF_Be9_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_q_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB10q
	ofile.open(MFB10q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_MF_B10_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_MF_B10_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_MF_B10_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_q_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB11q
	ofile.open(MFB11q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_MF_B11_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_MF_B11_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_MF_B11_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_q_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFC12q
	ofile.open(MFC12q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_MF_C12_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_MF_C12_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_MF_C12_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_q_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa40q
	ofile.open(MFCa40q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_MF_Ca40_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_MF_Ca40_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_MF_Ca40_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_qiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa48q
	ofile.open(MFCa48q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_MF_Ca48_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_MF_Ca48_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_MF_Ca48_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_qiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFFe54q
	ofile.open(MFFe54q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_MF_Fe54_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_MF_Fe54_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_MF_Fe54_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFAu197q
	ofile.open(MFAu197q.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int q_b = 1; q_b <= Nbins_q; q_b++)
	{
		//Get 2d bin number
		ib = H1_q_ratio_MF_Au197_pid_acc_kin_full->GetBin(q_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_q_ratio_MF_Au197_pid_acc_kin_full->GetBinContent(q_b);
		z_cont_err = H1_q_ratio_MF_Au197_pid_acc_kin_full->GetBinError(q_b);

		xlow = H1_q_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(q_b);
		xup = H1_q_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(q_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(q_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();



















	//MFBe9edelta
	ofile.open(MFBe9Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b_mf = 1; Ef_b_mf <= Nbins_Ef_mf; Ef_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->GetBin(Ef_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinContent(Ef_b_mf);
		z_cont_err = H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->GetBinError(Ef_b_mf);

		xlow = H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Ef_b_mf);
		xup = H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Ef_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_Be9_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();



	//MFB10edelta
	ofile.open(MFB10Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b_mf = 1; Ef_b_mf <= Nbins_Ef_mf; Ef_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->GetBin(Ef_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->GetBinContent(Ef_b_mf);
		z_cont_err = H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->GetBinError(Ef_b_mf);

		xlow = H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Ef_b_mf);
		xup = H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Ef_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_B10_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB11edelta
	ofile.open(MFB11Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b_mf = 1; Ef_b_mf <= Nbins_Ef_mf; Ef_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->GetBin(Ef_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->GetBinContent(Ef_b_mf);
		z_cont_err = H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->GetBinError(Ef_b_mf);

		xlow = H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Ef_b_mf);
		xup = H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Ef_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_B11_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFC12edelta
	ofile.open(MFC12Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b_mf = 1; Ef_b_mf <= Nbins_Ef_mf; Ef_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(Ef_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(Ef_b_mf);
		z_cont_err = H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(Ef_b_mf);

		xlow = H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Ef_b_mf);
		xup = H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Ef_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa40edelta
	ofile.open(MFCa40Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b_mf = 1; Ef_b_mf <= Nbins_Ef_mf; Ef_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBin(Ef_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinContent(Ef_b_mf);
		z_cont_err = H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->GetBinError(Ef_b_mf);

		xlow = H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Ef_b_mf);
		xup = H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Ef_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_Ca40_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa48edelta
	ofile.open(MFCa48Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b_mf = 1; Ef_b_mf <= Nbins_Ef_mf; Ef_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBin(Ef_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinContent(Ef_b_mf);
		z_cont_err = H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->GetBinError(Ef_b_mf);

		xlow = H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Ef_b_mf);
		xup = H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Ef_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_Ca48_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFFe54edelta
	ofile.open(MFFe54Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b_mf = 1; Ef_b_mf <= Nbins_Ef_mf; Ef_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBin(Ef_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinContent(Ef_b_mf);
		z_cont_err = H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->GetBinError(Ef_b_mf);

		xlow = H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Ef_b_mf);
		xup = H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Ef_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_Fe54_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFAu197edelta
	ofile.open(MFAu197Ef.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Ef_b_mf = 1; Ef_b_mf <= Nbins_Ef_mf; Ef_b_mf++)
	{
		//Get 2d bin number
		ib = H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->GetBin(Ef_b_mf);

		//Get bin content, error, and center for each bin
		z_cont = H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinContent(Ef_b_mf);
		z_cont_err = H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->GetBinError(Ef_b_mf);

		xlow = H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(Ef_b_mf);
		xup = H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(Ef_b_mf);
		x0 = (xlow + xup) / 2.; //H1_Ef_ratio_MF_Au197_pid_acc_kin_full_mf->GetXaxis()->GetBinCenter(Ef_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();




	





























	//MFBe9the
	ofile.open(MFBe9the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_MF_Be9_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_MF_Be9_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_MF_Be9_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_the_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB10the
	ofile.open(MFB10the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_MF_B10_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_MF_B10_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_MF_B10_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_the_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB11the
	ofile.open(MFB11the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_MF_B11_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_MF_B11_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_MF_B11_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_the_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFC12the
	ofile.open(MFC12the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_MF_C12_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_MF_C12_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_MF_C12_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_the_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa40the
	ofile.open(MFCa40the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_MF_Ca40_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_MF_Ca40_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_MF_Ca40_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_theiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa48the
	ofile.open(MFCa48the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_MF_Ca48_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_MF_Ca48_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_MF_Ca48_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_theiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFFe54the
	ofile.open(MFFe54the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_MF_Fe54_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_MF_Fe54_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_MF_Fe54_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFAu197the
	ofile.open(MFAu197the.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int the_b = 1; the_b <= Nbins_the; the_b++)
	{
		//Get 2d bin number
		ib = H1_the_ratio_MF_Au197_pid_acc_kin_full->GetBin(the_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_the_ratio_MF_Au197_pid_acc_kin_full->GetBinContent(the_b);
		z_cont_err = H1_the_ratio_MF_Au197_pid_acc_kin_full->GetBinError(the_b);

		xlow = H1_the_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(the_b);
		xup = H1_the_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(the_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(the_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();








	//MFBe9thp
	ofile.open(MFBe9thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_MF_Be9_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_MF_Be9_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_MF_Be9_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thp_ratio_MF_Be9_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB10thp
	ofile.open(MFB10thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_MF_B10_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_MF_B10_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_MF_B10_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thp_ratio_MF_B10_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFB11thp
	ofile.open(MFB11thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_MF_B11_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_MF_B11_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_MF_B11_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thp_ratio_MF_B11_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFC12thp
	ofile.open(MFC12thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_MF_C12_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_MF_C12_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_MF_C12_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thp_ratio_MF_C12_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa40thp
	ofile.open(MFCa40thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_MF_Ca40_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_MF_Ca40_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_MF_Ca40_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thpiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFCa48thp
	ofile.open(MFCa48thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_MF_Ca48_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_MF_Ca48_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_MF_Ca48_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_MF_Ca48_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_thpiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFFe54thp
	ofile.open(MFFe54thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_MF_Fe54_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_MF_Fe54_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_MF_Fe54_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_MF_Fe54_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(thp_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


	//MFAu197thp
	ofile.open(MFAu197thp.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int thp_b = 1; thp_b <= Nbins_thp; thp_b++)
	{
		//Get 2d bin number
		ib = H1_thp_ratio_MF_Au197_pid_acc_kin_full->GetBin(thp_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_thp_ratio_MF_Au197_pid_acc_kin_full->GetBinContent(thp_b);
		z_cont_err = H1_thp_ratio_MF_Au197_pid_acc_kin_full->GetBinError(thp_b);

		xlow = H1_thp_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(thp_b);
		xup = H1_thp_ratio_MF_Au197_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(thp_b);
		x0 = (xlow + xup) / 2.; //H1_Q2iss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(edelta_b);

		ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();
























	//MFCa40Q2
	ofile.open(eeD2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int Em_b = 1; Em_b <= Nbins_Wee; Em_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetBin(Em_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetBinContent(Em_b);
		z_cont_err = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetBinError(Em_b);

		xlow = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(Em_b);
		xup = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(Em_b);
		x0 = (xlow + xup) / 2.; //H1_Emiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();



	//MFCa48Em
	ofile.open(eepD2.Data());
	//Write header to data file
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int W_b = 1; W_b <= Nbins_Weep; W_b++)
	{
		//Get 2d bin number
		ib = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetBin(W_b);

		//Get bin content, error, and center for each bin
		z_cont = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetBinContent(W_b);
		z_cont_err = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetBinError(W_b);

		xlow = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(W_b);
		xup = H1_W_ratio_SRC_D2_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(W_b);
		x0 = (xlow + xup) / 2.; //H1_Wiss_ratio_MF_Ca40_pid_acc_kin_full->GetXaxis()->GetBinCenter(Em_b);

		ofile  << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
	}// end x-bins loop
	ofile.close();


















	
	//TH2F* H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full = 0;
	TString MFC12exey = "MFC12exey.csv";
	Int_t xbins = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetNbinsY();
	Int_t ybins = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetNbinsX();
	Double_t xbin_width = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);
	Double_t ybin_width = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinWidth(1);
	ofile.open(MFC12exey.Data());
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int x = 1; x <= xbins; x++)
	{
		for (int y = 1; y <= ybins; y++)
		{
			//Get 2d bin number
			ib = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(x, y);

			//Get bin content, error, and center for each bin
			z_cont = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(x, y);
			z_cont_err = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(x, y);

			xlow = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(x);
			xup = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(x);
			x0 = (xlow + xup) / 2.; 

			ylow = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinLowEdge(y);
			yup = H2_eXColl_eYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinUpEdge(y);
			y0 = (ylow + yup) / 2.;

			ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
			ofile << std::setw(7) << ib << std::setw(10) << x << std::setw(10) << y << std::setw(14) << x0 << std::setw(12) << y0 << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;
		}
	}// end x-bins loop
	ofile.close();
	

	//TH2F* H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full = 0;
	TString SRCC12hxhy = "SRCC12hxhy.csv";
	xbins = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetNbinsY();
	ybins = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetNbinsX();
	xbin_width = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	ybin_width = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinWidth(1);
	ofile.open(SRCC12hxhy.Data());
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int x = 1; x <= xbins; x++)
	{
		for (int y = 1; y <= ybins; y++)
		{
			//Get 2d bin number
			ib = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetBin(x, y);

			//Get bin content, error, and center for each bin
			z_cont = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(x, y);
			z_cont_err = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetBinError(x, y);

			xlow = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(x);
			xup = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(x);
			x0 = (xlow + xup) / 2.;

			ylow = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinLowEdge(y);
			yup = H2_hXColl_hYColl_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinUpEdge(y);
			y0 = (ylow + yup) / 2.;

			ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
			ofile << std::setw(7) << ib << std::setw(10) << x << std::setw(10) << y << std::setw(14) << x0 << std::setw(12) << y0 << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;
		}
	}// end x-bins loop
	ofile.close();


	//TH2F* H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full = 0;
	TString MFC12hxhy = "MFC12hxhy.csv";
	xbins = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetNbinsY();
	ybins = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetNbinsX();
	xbin_width = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);
	ybin_width = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinWidth(1);
	ofile.open(MFC12hxhy.Data());
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int x = 1; x <= xbins; x++)
	{
		for (int y = 1; y <= ybins; y++)
		{
			//Get 2d bin number
			ib = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(x, y);

			//Get bin content, error, and center for each bin
			z_cont = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(x, y);
			z_cont_err = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(x, y);

			xlow = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(x);
			xup = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(x);
			x0 = (xlow + xup) / 2.;

			ylow = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinLowEdge(y);
			yup = H2_hXColl_hYColl_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinUpEdge(y);
			y0 = (ylow + yup) / 2.;

			ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
			ofile << std::setw(7) << ib << std::setw(10) << x << std::setw(10) << y << std::setw(14) << x0 << std::setw(12) << y0 << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;
		}
	}// end x-bins loop
	ofile.close(); 


	//TH2F* H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full = 0;
	TString SRCC12exey = "SRCC12exey.csv";
	xbins = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetNbinsY();
	ybins = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetNbinsX();
	xbin_width = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	ybin_width = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinWidth(1);
	ofile.open(SRCC12exey.Data());
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int x = 1; x <= xbins; x++)
	{
		for (int y = 1; y <= ybins; y++)
		{
			//Get 2d bin number
			ib = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetBin(x, y);

			//Get bin content, error, and center for each bin
			z_cont = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(x, y);
			z_cont_err = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetBinError(x, y);

			xlow = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(x);
			xup = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(x);
			x0 = (xlow + xup) / 2.;

			ylow = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinLowEdge(y);
			yup = H2_eXColl_eYColl_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinUpEdge(y);
			y0 = (ylow + yup) / 2.;

			ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
			ofile << std::setw(7) << ib << std::setw(10) << x << std::setw(10) << y << std::setw(14) << x0 << std::setw(12) << y0 << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;
		}
	}// end x-bins loop
	ofile.close();
	

	/*
	//TH2F* H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf = 0;
	TString MFC12EmPm = "MFC12EmPm.csv";
	xbins = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetNbinsY();
	ybins = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetNbinsX();
	xbin_width = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);
	ybin_width = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinWidth(1);
	ofile.open(MFC12EmPm.Data());
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int x = 1; x <= xbins; x++)
	{
		for (int y = 1; y <= ybins; y++)
		{
			//Get 2d bin number
			ib = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(x, y);

			//Get bin content, error, and center for each bin
			z_cont = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(x, y);
			z_cont_err = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(x, y);

			xlow = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(x);
			xup = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(x);
			x0 = (xlow + xup) / 2.;

			ylow = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinLowEdge(y);
			yup = H2_Em_Pm_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinUpEdge(y);
			y0 = (ylow + yup) / 2.;

			ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
			ofile << std::setw(7) << ib << std::setw(10) << x << std::setw(10) << y << std::setw(14) << x0 << std::setw(12) << y0 << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;
		}
	}// end x-bins loop
	ofile.close(); 
	*/
	
	
	//TH2F* H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full = 0;
	TString SRCC12EmPm = "SRCC12EmPm.csv";
	xbins = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetNbinsY();
	ybins = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetNbinsX();
	xbin_width = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	ybin_width = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinWidth(1);
	ofile.open(SRCC12EmPm.Data());
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int x = 1; x <= xbins; x++)
	{
		for (int y = 1; y <= ybins; y++)
		{
			//Get 2d bin number
			ib = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetBin(x, y);

			//Get bin content, error, and center for each bin
			z_cont = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(x, y);
			z_cont_err = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetBinError(x, y);

			xlow = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(x);
			xup = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(x);
			x0 = (xlow + xup) / 2.;

			ylow = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinLowEdge(y);
			yup = H2_Em_Pm_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinUpEdge(y);
			y0 = (ylow + yup) / 2.;

			ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
			ofile << std::setw(7) << ib << std::setw(10) << x << std::setw(10) << y << std::setw(14) << x0 << std::setw(12) << y0 << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;
		}
	}// end x-bins loop
	ofile.close();


	/*
	//TH2F* H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full = 0;
	TString MFC12XQ2 = "MFC12XQ2.csv";
	xbins = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetNbinsY();
	ybins = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetNbinsX();
	xbin_width = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinWidth(1);
	ybin_width = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinWidth(1);
	ofile.open(MFC12XQ2.Data());
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int x = 1; x <= xbins; x++)
	{
		for (int y = 1; y <= ybins; y++)
		{
			//Get 2d bin number
			ib = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetBin(x, y);

			//Get bin content, error, and center for each bin
			z_cont = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetBinContent(x, y);
			z_cont_err = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetBinError(x, y);

			xlow = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinLowEdge(x);
			xup = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetXaxis()->GetBinUpEdge(x);
			x0 = (xlow + xup) / 2.;

			ylow = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinLowEdge(y);
			yup = H2_xbj_Q2_ratio_MF_C12_pid_acc_kin_full_mf->GetYaxis()->GetBinUpEdge(y);
			y0 = (ylow + yup) / 2.;

			ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
			ofile << std::setw(7) << ib << std::setw(10) << x << std::setw(10) << y << std::setw(14) << x0 << std::setw(12) << y0 << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;
		}
	}// end x-bins loop
	ofile.close(); 
	*/
	
	
	//TH2F* H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full = 0;
	TString SRCC12XQ2 = "SRCC12XQ2.csv";
	xbins = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetNbinsY();
	ybins = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetNbinsX();
	xbin_width = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinWidth(1);
	ybin_width = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinWidth(1);
	ofile.open(SRCC12XQ2.Data());
	ofile << header.Data() << endl;
	//loop over x-bins
	for (int x = 1; x <= xbins; x++)
	{
		for (int y = 1; y <= ybins; y++)
		{
			//Get 2d bin number
			ib = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetBin(x, y);

			//Get bin content, error, and center for each bin
			z_cont = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetBinContent(x, y);
			z_cont_err = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetBinError(x, y);

			xlow = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinLowEdge(x);
			xup = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetXaxis()->GetBinUpEdge(x);
			x0 = (xlow + xup) / 2.;

			ylow = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinLowEdge(y);
			yup = H2_xbj_Q2_ratio_SRC_C12_pid_acc_kin_full->GetYaxis()->GetBinUpEdge(y);
			y0 = (ylow + yup) / 2.;

			ofile << std::setw(14) << x0 << "," << std::setw(12) << z_cont << "," << std::setw(14) << z_cont_err << endl;
			ofile << std::setw(7) << ib << std::setw(10) << x << std::setw(10) << y << std::setw(14) << x0 << std::setw(12) << y0 << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;
		}
	}// end x-bins loop
	ofile.close();
	


	

	







	/*cout << "Be9 SRC 3: " << Be9_SRC_3 << endl;
	cout << "B10 SRC 3: " << B10_SRC_3 << endl;
	cout << "B11 SRC 3: " << B11_SRC_3 << endl;
	cout << "C12 SRC 3: " << C12_SRC_3 << endl;
	cout << "Ca40 SRC 3: " << Ca40_SRC_3 << endl;
	cout << "Ca48 SRC 3: " << Ca48_SRC_3 << endl;
	cout << "Fe54 SRC 3: " << Fe54_SRC_3 << endl;
	cout << "Ca40 SRC Q2 3: " << Ca40_SRC_Q23 << endl;
	cout << "Ca48 SRC Q2 3: " << Ca48_SRC_Q23 << endl;
	cout << "Fe54 SRC Q2 3: " << Fe54_SRC_Q23 << endl;
	cout << "Ca40 SRC X 3: " << Ca40_SRC_X3 << endl;
	cout << "Ca48 SRC X 3: " << Ca48_SRC_X3 << endl;
	cout << "Fe54 SRC X 3: " << Fe54_SRC_X3 << endl;
	cout << "Au197 SRC 3: " << Au197_SRC_3 << endl;*/
	outROOT->Close();	//Close File
}
//*****************************************************************************************************************************************************************************
//End Main
//*****************************************************************************************************************************************************************************








void thingy(
	TH1F* hist1, TH1F* hist2, TH1F* hist3,
	TH1F* hist4, TH1F* hist5, TH1F* hist6,
	TH1F* hist7, TH1F* hist8,
	TH1F* hist9, TH1F* hist10,
	TFile* outROOT)
{
	int font_type = 132;
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);//0.12
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);


	hist1->SetLineWidth(2);
	hist2->SetLineWidth(2);
	hist3->SetLineWidth(2);
	hist4->SetLineWidth(2);
	hist5->SetLineWidth(2);
	hist6->SetLineWidth(2);
	hist7->SetLineWidth(2);
	hist8->SetLineWidth(2);
	hist9->SetLineWidth(2);
	hist10->SetLineWidth(2);


	// set histos aethetics
	hist1->SetLineColor(kRed);
	//hist1->SetFillColorAlpha(kRed, 0.40);
	//hist1->SetFillStyle(3004);

	hist2->SetLineColor(kGreen);
	//hist2->SetFillColorAlpha(kGreen, 0.40);
	//hist2->SetFillStyle(3005);

	hist3->SetLineColor(kBlue);
	//hist3->SetFillColorAlpha(kBlue, 0.40);
	//hist3->SetFillStyle(3006);

	//TLegend* leg = new TLegend(0.64, 0.89, 0.75, 0.78);
	//leg->AddEntry(hist1, "Be9/C12");
	//leg->AddEntry(hist2, "B10/C12");
	//leg->AddEntry(hist3, "B11/C12");
	//leg->Draw();

	hist4->SetLineColor(kRed);
	//hist4->SetFillColorAlpha(kRed, 0.40);
	//hist4->SetFillStyle(3004);

	hist5->SetLineColor(kGreen);
	//hist5->SetFillColorAlpha(kGreen, 0.40);
	//hist5->SetFillStyle(3005);

	hist6->SetLineColor(kBlue);
	//hist6->SetFillColorAlpha(kBlue, 0.40);
	//hist6->SetFillStyle(3006);

	hist7->SetLineColor(kRed);
	//hist7->SetFillColorAlpha(kRed, 0.40);
	//hist7->SetFillStyle(3004);

	hist8->SetLineColor(kBlue);
	//hist8->SetFillColorAlpha(kBlue, 0.40);
	//hist8->SetFillStyle(3006);

	hist9->SetLineColor(kRed);
	//hist9->SetFillColorAlpha(kRed, 0.40);
	//hist9->SetFillStyle(3004);

	hist10->SetLineColor(kBlue);
	//hist10->SetFillColorAlpha(kBlue, 0.40);
	//hist10->SetFillStyle(3006);


	hist1->GetYaxis()->SetRangeUser(0, 2.0);
	hist1->SetTitle("Light Nuclei");
	hist1->GetXaxis()->SetLabelSize(0.0);//0.04
	hist1->GetYaxis()->SetLabelSize(0.04);
	//hist1->GetYaxis()->SetTitle("Counts / mc");
	hist1->GetXaxis()->SetTitle("");//
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");
	hist2->GetYaxis()->SetRangeUser(0, 2.0);
	hist3->GetYaxis()->SetRangeUser(0, 2.0);
	hist4->GetYaxis()->SetRangeUser(0, 2.0);
	hist5->GetYaxis()->SetRangeUser(0, 2.0);
	hist6->GetYaxis()->SetRangeUser(0, 2.0);
	
	
	hist7->GetYaxis()->SetRangeUser(0, 2.0);
	hist7->SetTitle("Heavy Nuclei");
	hist7->GetXaxis()->SetLabelSize(0.04);
	hist7->GetYaxis()->SetLabelSize(0.04);
	//hist1->GetYaxis()->SetTitle("Counts / mc");
	//hist1->GetXaxis()->SetTitle(xlabel);
	hist7->GetYaxis()->CenterTitle();
	hist7->GetXaxis()->CenterTitle();
	hist7->SetLabelFont(font_type, "XY");
	hist7->SetTitleFont(font_type, "XY");
	hist7->SetTitleSize(0.05, "XY");
	hist7->SetTitleOffset(1., "XY"); 
	hist8->GetYaxis()->SetRangeUser(0, 2.0); 
	hist9->GetYaxis()->SetRangeUser(0, 2.0); 
	hist10->GetYaxis()->SetRangeUser(0, 2.0);

	


	TCanvas* c = new TCanvas("c", "c", 1024, 768); c->Divide(1, 2);

	c->cd(1);
	hist1->Draw("histE0");
	hist2->Draw("sameshistE0");
	hist3->Draw("sameshistE0");
	hist4->Draw("sameshistE0");
	hist5->Draw("sameshistE0");
	hist6->Draw("sameshistE0");

	//TLegend* leg = new TLegend(0.64, 0.89, 0.75, 0.78);
	//leg->AddEntry(hist1, "Be9/C12");
	//leg->AddEntry(hist2, "B10/C12");
	//leg->AddEntry(hist3, "B11/C12");
	//leg->Draw();


	c->cd(2);
	hist7->Draw("histE0");
	hist8->Draw("sameshistE0");
	hist9->Draw("sameshistE0");
	hist10->Draw("sameshistE0");

	//TLegend* leg2 = new TLegend(0.44, 0.69, 0.55, 0.58);
	//leg2->AddEntry(hist7, "Ca48/Ca40");
	//leg2->AddEntry(hist8, "Fe54/Ca40");
	//leg2->Draw();

	outROOT->cd(); c->Write(); c->Print("thing.png");

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}







void Draw1dE(TH1F* hist1, const char* name, const char* title, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(1);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");

	TCanvas* c = new TCanvas(name, title, 1366, 768); c->cd();

	hist1->Draw("histE");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}













void thingpap2(
	TH1F* hist1, TH1F* hist2, TH1F* hist3, const char* name, const char* title, TFile* outROOT)
{
	int font_type = 132;
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.12);//0.05
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);



	Double_t yaxisrange = 0;
	yaxisrange = hist1->GetMaximum();
	if (yaxisrange < hist2->GetMaximum()) { yaxisrange = hist2->GetMaximum(); }
	if (yaxisrange < hist3->GetMaximum()) { yaxisrange = hist3->GetMaximum(); }
	hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.3 * yaxisrange));

	hist1->SetLineWidth(2);
	hist1->SetLineColor(kRed);
	//hist1->SetFillColorAlpha(kRed, 0.40);
	//hist1->SetFillStyle(3004);

	hist2->SetLineWidth(2);
	hist2->SetLineColor(kGreen);
	//hist2->SetFillColorAlpha(kGreen, 0.40);
	//hist2->SetFillStyle(3005);
	
	hist3->SetLineWidth(2);
	hist3->SetLineColor(kBlue);
	//hist2->SetFillColorAlpha(kGreen, 0.40);
	//hist2->SetFillStyle(3005);


	//TGraphError
	//Try free floating data points
	//Try 1 parameter polynomial fit in root
	//

	hist1->GetXaxis()->SetRangeUser(0.0, 3.0);
	hist1->SetTitle("");
	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	//hist1->GetYaxis()->SetTitle("Counts / mc");
	hist1->GetXaxis()->SetTitle("Pm");
	//hist1->GetYaxis()->SetTitle("SRC Ratio per proton");
	hist1->GetYaxis()->SetTitle("SRC Ratio per nucleus");
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");
	//hist2->GetXaxis()->SetRangeUser(0.0, 2.0);
	//hist3->GetXaxis()->SetRangeUser(0.0, 2.0);

	TCanvas* c = new TCanvas("c", "c", 1024, 768);

	c->cd();
	hist1->Draw("histE0");
	hist2->Draw("sameshistE0");
	hist3->Draw("sameshistE0");



	//TLegend* leg2 = new TLegend(0.44, 0.69, 0.55, 0.58);
	TLegend* leg2 = new TLegend(0.64, 0.89, 0.75, 0.78);
	leg2->AddEntry(hist1, "Ca40");
	leg2->AddEntry(hist2, "Ca48");
	leg2->AddEntry(hist3, "Fe54");
	leg2->Draw();

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}























void thingpap(
	TH1F* hist1, TH1F* hist2, TH1F* hist3, TFile* outROOT)
{
	int font_type = 132;
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.12);//0.05
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);


	//hist1->SetLineWidth(2);
	//hist2->SetLineWidth(2);
	//hist2->SetLineWidth(2);

	// set histos aethetics
	hist1->SetLineColor(kRed);
	//hist1->SetFillColorAlpha(kRed, 0.40);
	//hist1->SetFillStyle(3004);

	//hist2->SetLineColor(kGreen);
	//hist2->SetFillColorAlpha(kGreen, 0.40);
	//hist2->SetFillStyle(3005);

	hist3->SetLineColor(kBlue);
	//hist2->SetFillColorAlpha(kGreen, 0.40);
	//hist2->SetFillStyle(3005);


	//TGraphError
	//Try free floating data points
	//Try 1 parameter polynomial fit in root
	//

	hist1->GetYaxis()->SetRangeUser(0.0, 1.6);
	//hist1->GetYaxis()->SetRangeUser(0.0, 1.7);
	hist1->GetXaxis()->SetRangeUser(0.35, 0.705);
	hist1->SetTitle("");
	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	//hist1->GetYaxis()->SetTitle("Counts / mc");
	hist1->GetXaxis()->SetTitle("Pm");
	//hist1->GetYaxis()->SetTitle("SRC Ratio per proton");
	hist1->GetYaxis()->SetTitle("SRC Ratio per nucleus");
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");
	//hist2->GetYaxis()->SetRangeUser(0.0, 1.6);
	//hist2->GetXaxis()->SetRangeUser(0.35, 0.705);
	hist3->GetYaxis()->SetRangeUser(0.0, 1.6);
	//hist3->GetYaxis()->SetRangeUser(0.0, 1.7);
	hist3->GetXaxis()->SetRangeUser(0.3, 0.705);

	TCanvas* c = new TCanvas("c", "c", 1024, 768);

	c->cd();
	hist1->Draw("histE0");
	//hist2->Draw("sameshistE0");
	hist3->Draw("sameshistE0");



	//TLegend* leg2 = new TLegend(0.44, 0.69, 0.55, 0.58);
	TLegend* leg2 = new TLegend(0.64, 0.89, 0.75, 0.78);
	leg2->AddEntry(hist1, "Ca48/Ca40");
	//leg2->AddEntry(hist2, "Fe54/Ca40");
	leg2->AddEntry(hist3, "Fe54/Ca48");
	leg2->Draw();

	outROOT->cd(); c->Write(); c->Print("thingpap.png");

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




/*
	TCutG* contam_gCut = new TCutG("contamCut", 5);
	contam_gCut->SetVarX("X");
	contam_gCut->SetVarY("Y");

	contam_gCut->SetPoint(0, -0.02, 0.0);
	contam_gCut->SetPoint(1, -0.02, 0.04);
	contam_gCut->SetPoint(3, 0.06, 0.04);
	contam_gCut->SetPoint(2, 0.07, 0.03);
	contam_gCut->SetPoint(4, 0.04, 0.0);
	contam_gCut->SetPoint(5, -0.02, 0.0);

	TLine* line1 = new TLine(x1i, y1i, x1f, y1f); line1->SetLineColor(2); line1->SetLineWidth(2); line1->Draw("SAME");
*/


void Draw2d(TH2F* hist1, const char* name, const char* title, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");


	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();

	hist1->Draw("Col");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}






















void Draw2d2(TH2F* hist1, const char* name, const char* title, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");


	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();

	hist1->Draw("Col");



	//contam_gCut->SetPoint(0, -0.02, 0.0);
	//contam_gCut->SetPoint(1, -0.02, 0.04);
	//contam_gCut->SetPoint(3, 0.06, 0.04);
	//contam_gCut->SetPoint(2, 0.07, 0.03);
	//contam_gCut->SetPoint(4, 0.04, 0.0);
	//contam_gCut->SetPoint(5, -0.02, 0.0);
	TLine* line1 = new TLine(0.0, -0.02, 0.04, -0.02); line1->SetLineColor(2); line1->SetLineWidth(5); line1->Draw("SAME");
	TLine* line2 = new TLine(0.04, -0.02, 0.04, 0.06); line2->SetLineColor(2); line2->SetLineWidth(5); line2->Draw("SAME");
	TLine* line3 = new TLine(0.04, 0.06, 0.03, 0.07); line3->SetLineColor(2); line3->SetLineWidth(5); line3->Draw("SAME");
	TLine* line4 = new TLine(0.03, 0.07, 0.0, 0.04); line4->SetLineColor(2); line4->SetLineWidth(5); line4->Draw("SAME");
	TLine* line5 = new TLine(0.0, 0.04, 0.0, -0.02); line5->SetLineColor(2); line5->SetLineWidth(5); line5->Draw("SAME");

	//P_{m} \le 0.275$ GeV/c and $-0.02 \le E_{m} \le0.09$
	TLine* line6 = new TLine(0.0, -0.02, 0.275, -0.02); line6->SetLineColor(2); line6->SetLineWidth(5); line6->Draw("SAME");
	TLine* line7 = new TLine(0.275, -0.02, 0.275, 0.09); line7->SetLineColor(2); line7->SetLineWidth(5); line7->Draw("SAME");
	TLine* line8 = new TLine(0.275, 0.09, 0.0, 0.09); line8->SetLineColor(2); line8->SetLineWidth(5); line8->Draw("SAME");
	TLine* line9 = new TLine(0.0, 0.09, 0.0, -0.02); line9->SetLineColor(2); line9->SetLineWidth(5); line9->Draw("SAME");


	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}











//Scaling the HMS/SHMS Collimator Cuts













void Draw2dh(TH2F* hist1, const char* name, const char* title, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");

	/*Double_t hms_scale = 1;
	Double_t hms_hsize = 4.575;
	Double_t hms_vsize = 11.646;
	hms_hsize = hms_scale * hms_hsize;
	hms_vsize = hms_scale * hms_vsize;

	hms_hsize, hms_vsize / 2.);
	hms_hsize / 2., hms_vsize);
	hms_hsize / 2., hms_vsize);
	hms_hsize, hms_vsize / 2.);
	hms_hsize, -hms_vsize / 2.);
	hms_hsize / 2., -hms_vsize);
	hms_hsize / 2., -hms_vsize);
	hms_hsize, -hms_vsize / 2.);
	hms_hsize, hms_vsize / 2.);*/

	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();

	hist1->Draw("Col");

	Double_t hms_scale_1 = 1.08;
	Double_t hms_scale_2 = 0.92;
	Double_t hms_hsize = 4.575;
	Double_t hms_vsize = 11.646;
	//Scaling the HMS/SHMS Collimator Cuts
	Double_t hms_hsize_1 = hms_scale_1 * hms_hsize;
	Double_t hms_vsize_1 = hms_scale_1 * hms_vsize;
	Double_t hms_hsize_2 = hms_scale_2 * hms_hsize;
	Double_t hms_vsize_2 = hms_scale_2 * hms_vsize;

	TLine* line1 = new TLine(hms_hsize, hms_vsize / 2, hms_hsize / 2, hms_vsize); line1->SetLineColor(2); line1->SetLineWidth(5); line1->Draw("SAME");
	TLine* line2 = new TLine(hms_hsize / 2, hms_vsize, -hms_hsize / 2, hms_vsize); line2->SetLineColor(2); line2->SetLineWidth(5); line2->Draw("SAME");
	TLine* line3 = new TLine(-hms_hsize / 2, hms_vsize, -hms_hsize, hms_vsize / 2); line3->SetLineColor(2); line3->SetLineWidth(5); line3->Draw("SAME");
	TLine* line4 = new TLine(-hms_hsize, hms_vsize / 2, -hms_hsize, -hms_vsize / 2); line4->SetLineColor(2); line4->SetLineWidth(5); line4->Draw("SAME");
	TLine* line5 = new TLine(-hms_hsize, -hms_vsize / 2, -hms_hsize / 2, -hms_vsize); line5->SetLineColor(2); line5->SetLineWidth(5); line5->Draw("SAME");
	TLine* line6 = new TLine(-hms_hsize / 2, -hms_vsize, hms_hsize / 2, -hms_vsize); line6->SetLineColor(2); line6->SetLineWidth(5); line6->Draw("SAME");
	TLine* line7 = new TLine(hms_hsize / 2, -hms_vsize, hms_hsize, -hms_vsize / 2); line7->SetLineColor(2); line7->SetLineWidth(5); line7->Draw("SAME");
	TLine* line8 = new TLine(hms_hsize, -hms_vsize / 2, hms_hsize, hms_vsize / 2); line8->SetLineColor(2); line8->SetLineWidth(5); line8->Draw("SAME");

	TLine* line9 = new TLine(hms_hsize_1, hms_vsize_1 / 2, hms_hsize_1 / 2, hms_vsize_1); line9->SetLineColor(2); line9->SetLineWidth(5); line9->SetLineStyle(10); line9->Draw("SAME");
	TLine* line10 = new TLine(hms_hsize_1 / 2, hms_vsize_1, -hms_hsize_1 / 2, hms_vsize_1); line10->SetLineColor(2); line10->SetLineWidth(5); line10->SetLineStyle(10); line10->Draw("SAME");
	TLine* line11 = new TLine(-hms_hsize_1 / 2, hms_vsize_1, -hms_hsize_1, hms_vsize_1 / 2); line11->SetLineColor(2); line11->SetLineWidth(5); line11->SetLineStyle(10); line11->Draw("SAME");
	TLine* line12 = new TLine(-hms_hsize_1, hms_vsize_1 / 2, -hms_hsize_1, -hms_vsize_1 / 2); line12->SetLineColor(2); line12->SetLineWidth(5); line12->SetLineStyle(10); line12->Draw("SAME");
	TLine* line13 = new TLine(-hms_hsize_1, -hms_vsize_1 / 2, -hms_hsize_1 / 2, -hms_vsize_1); line13->SetLineColor(2); line13->SetLineWidth(5); line13->SetLineStyle(10); line13->Draw("SAME");
	TLine* line14 = new TLine(-hms_hsize_1 / 2, -hms_vsize_1, hms_hsize_1 / 2, -hms_vsize_1); line14->SetLineColor(2); line14->SetLineWidth(5); line14->SetLineStyle(10); line14->Draw("SAME");
	TLine* line15 = new TLine(hms_hsize_1 / 2, -hms_vsize_1, hms_hsize_1, -hms_vsize_1 / 2); line15->SetLineColor(2); line15->SetLineWidth(5); line15->SetLineStyle(10); line15->Draw("SAME");
	TLine* line16 = new TLine(hms_hsize_1, -hms_vsize_1 / 2, hms_hsize_1, hms_vsize_1 / 2); line16->SetLineColor(2); line16->SetLineWidth(5); line16->SetLineStyle(10); line16->Draw("SAME");

	TLine* line17 = new TLine(hms_hsize_2, hms_vsize_2 / 2, hms_hsize_2 / 2, hms_vsize_2); line17->SetLineColor(2); line17->SetLineWidth(5); line17->SetLineStyle(10); line17->Draw("SAME");
	TLine* line18 = new TLine(hms_hsize_2 / 2, hms_vsize_2, -hms_hsize_2 / 2, hms_vsize_2); line18->SetLineColor(2); line18->SetLineWidth(5); line18->SetLineStyle(10); line18->Draw("SAME");
	TLine* line19 = new TLine(-hms_hsize_2 / 2, hms_vsize_2, -hms_hsize_2, hms_vsize_2 / 2); line19->SetLineColor(2); line19->SetLineWidth(5); line19->SetLineStyle(10); line19->Draw("SAME");
	TLine* line20 = new TLine(-hms_hsize_2, hms_vsize_2 / 2, -hms_hsize_2, -hms_vsize_2 / 2); line20->SetLineColor(2); line20->SetLineWidth(5); line20->SetLineStyle(10); line20->Draw("SAME");
	TLine* line21 = new TLine(-hms_hsize_2, -hms_vsize_2 / 2, -hms_hsize_2 / 2, -hms_vsize_2); line21->SetLineColor(2); line21->SetLineWidth(5); line21->SetLineStyle(10); line21->Draw("SAME");
	TLine* line22 = new TLine(-hms_hsize_2 / 2, -hms_vsize_2, hms_hsize_2 / 2, -hms_vsize_2); line22->SetLineColor(2); line22->SetLineWidth(5); line22->SetLineStyle(10); line22->Draw("SAME");
	TLine* line23 = new TLine(hms_hsize_2 / 2, -hms_vsize_2, hms_hsize_2, -hms_vsize_2 / 2); line23->SetLineColor(2); line23->SetLineWidth(5); line23->SetLineStyle(10); line23->Draw("SAME");
	TLine* line24 = new TLine(hms_hsize_2, -hms_vsize_2 / 2, hms_hsize_2, hms_vsize_2 / 2); line24->SetLineColor(2); line24->SetLineWidth(5); line24->SetLineStyle(10); line24->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}
















void Draw2de(TH2F* hist1, const char* name, const char* title, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");

	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();

	hist1->Draw("Col");

	Double_t shms_scale_1 = 1.08;
	Double_t shms_scale_2 = 0.92;
	Double_t shms_hsize = 8.5;
	Double_t shms_vsize = 12.5;
	Double_t shms_hsize_1 = shms_scale_1 * shms_hsize;
	Double_t shms_vsize_1 = shms_scale_1 * shms_vsize;
	Double_t shms_hsize_2 = shms_scale_2 * shms_hsize;
	Double_t shms_vsize_2 = shms_scale_2 * shms_vsize;

	TLine* line1 = new TLine(shms_hsize, shms_vsize / 2, shms_hsize / 2, shms_vsize); line1->SetLineColor(2); line1->SetLineWidth(5); line1->Draw("SAME");
	TLine* line2 = new TLine(shms_hsize / 2, shms_vsize, -shms_hsize / 2, shms_vsize); line2->SetLineColor(2); line2->SetLineWidth(5); line2->Draw("SAME");
	TLine* line3 = new TLine(-shms_hsize / 2, shms_vsize, -shms_hsize, shms_vsize / 2); line3->SetLineColor(2); line3->SetLineWidth(5); line3->Draw("SAME");
	TLine* line4 = new TLine(-shms_hsize, shms_vsize / 2, -shms_hsize, -shms_vsize / 2); line4->SetLineColor(2); line4->SetLineWidth(5); line4->Draw("SAME");
	TLine* line5 = new TLine(-shms_hsize, -shms_vsize / 2, -shms_hsize / 2, -shms_vsize); line5->SetLineColor(2); line5->SetLineWidth(5); line5->Draw("SAME");
	TLine* line6 = new TLine(-shms_hsize / 2, -shms_vsize, shms_hsize / 2, -shms_vsize); line6->SetLineColor(2); line6->SetLineWidth(5); line6->Draw("SAME");
	TLine* line7 = new TLine(shms_hsize / 2, -shms_vsize, shms_hsize, -shms_vsize / 2); line7->SetLineColor(2); line7->SetLineWidth(5); line7->Draw("SAME");
	TLine* line8 = new TLine(shms_hsize, -shms_vsize / 2, shms_hsize, shms_vsize / 2); line8->SetLineColor(2); line8->SetLineWidth(5); line8->Draw("SAME");

	TLine* line9 = new TLine(shms_hsize_1, shms_vsize_1 / 2, shms_hsize_1 / 2, shms_vsize_1); line9->SetLineColor(2); line9->SetLineWidth(5); line9->SetLineStyle(10); line9->Draw("SAME");
	TLine* line10 = new TLine(shms_hsize_1 / 2, shms_vsize_1, -shms_hsize_1 / 2, shms_vsize_1); line10->SetLineColor(2); line10->SetLineWidth(5); line10->SetLineStyle(10); line10->Draw("SAME");
	TLine* line11 = new TLine(-shms_hsize_1 / 2, shms_vsize_1, -shms_hsize_1, shms_vsize_1 / 2); line11->SetLineColor(2); line11->SetLineWidth(5); line11->SetLineStyle(10); line11->Draw("SAME");
	TLine* line12 = new TLine(-shms_hsize_1, shms_vsize_1 / 2, -shms_hsize_1, -shms_vsize_1 / 2); line12->SetLineColor(2); line12->SetLineWidth(5); line12->SetLineStyle(10); line12->Draw("SAME");
	TLine* line13 = new TLine(-shms_hsize_1, -shms_vsize_1 / 2, -shms_hsize_1 / 2, -shms_vsize_1); line13->SetLineColor(2); line13->SetLineWidth(5); line13->SetLineStyle(10); line13->Draw("SAME");
	TLine* line14 = new TLine(-shms_hsize_1 / 2, -shms_vsize_1, shms_hsize_1 / 2, -shms_vsize_1); line14->SetLineColor(2); line14->SetLineWidth(5); line14->SetLineStyle(10); line14->Draw("SAME");
	TLine* line15 = new TLine(shms_hsize_1 / 2, -shms_vsize_1, shms_hsize_1, -shms_vsize_1 / 2); line15->SetLineColor(2); line15->SetLineWidth(5); line15->SetLineStyle(10); line15->Draw("SAME");
	TLine* line16 = new TLine(shms_hsize_1, -shms_vsize_1 / 2, shms_hsize_1, shms_vsize_1 / 2); line16->SetLineColor(2); line16->SetLineWidth(5); line16->SetLineStyle(10); line16->Draw("SAME");

	TLine* line17 = new TLine(shms_hsize_2, shms_vsize_2 / 2, shms_hsize_2 / 2, shms_vsize_2); line17->SetLineColor(2); line17->SetLineWidth(5); line17->SetLineStyle(10); line17->Draw("SAME");
	TLine* line18 = new TLine(shms_hsize_2 / 2, shms_vsize_2, -shms_hsize_2 / 2, shms_vsize_2); line18->SetLineColor(2); line18->SetLineWidth(5); line18->SetLineStyle(10); line18->Draw("SAME");
	TLine* line19 = new TLine(-shms_hsize_2 / 2, shms_vsize_2, -shms_hsize_2, shms_vsize_2 / 2); line19->SetLineColor(2); line19->SetLineWidth(5); line19->SetLineStyle(10); line19->Draw("SAME");
	TLine* line20 = new TLine(-shms_hsize_2, shms_vsize_2 / 2, -shms_hsize_2, -shms_vsize_2 / 2); line20->SetLineColor(2); line20->SetLineWidth(5); line20->SetLineStyle(10); line20->Draw("SAME");
	TLine* line21 = new TLine(-shms_hsize_2, -shms_vsize_2 / 2, -shms_hsize_2 / 2, -shms_vsize_2); line21->SetLineColor(2); line21->SetLineWidth(5); line21->SetLineStyle(10); line21->Draw("SAME");
	TLine* line22 = new TLine(-shms_hsize_2 / 2, -shms_vsize_2, shms_hsize_2 / 2, -shms_vsize_2); line22->SetLineColor(2); line22->SetLineWidth(5); line22->SetLineStyle(10); line22->Draw("SAME");
	TLine* line23 = new TLine(shms_hsize_2 / 2, -shms_vsize_2, shms_hsize_2, -shms_vsize_2 / 2); line23->SetLineColor(2); line23->SetLineWidth(5); line23->SetLineStyle(10); line23->Draw("SAME");
	TLine* line24 = new TLine(shms_hsize_2, -shms_vsize_2 / 2, shms_hsize_2, shms_vsize_2 / 2); line24->SetLineColor(2); line24->SetLineWidth(5); line24->SetLineStyle(10); line24->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}
