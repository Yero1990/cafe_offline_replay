

//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
// This file when run asks for 4 inputs
// The kinematics: MF or SRC
// The target: LD2, Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
// Pass: pass1, pass2, pass3, pass4
// Singles or Coincidence: single, coin (note that the single option hasn't been used in forever and is not need and thus won't run without some work
// Mind the file paths and everything should run
// This file looks complex but is very simple it can be broken into 1.) Variable/Histogram/List declarations 2.) Reading in from report files
// 3.) 2 nested for loops where the outer one goes over files and the inner one goes over events and fills histograms
// 4.) function calls to print the made histograms as .png files
// Most of the apparent complexity comes from the sheer number of unique histograms I've been asked to make ~4 years
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
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
#include <stdio.h>

#include "../../header_files/Histo.h"
#include "../../header_files/Cuts.h"
#include "../../header_files/parse_utils.h"
#include "../../header_files/hist_utils.h"
using namespace std;


//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//-----------------------------------------------------------------------------Constants----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
Double_t pi = TMath::Pi();			//Pi
Double_t dtr = pi/180.;//pi				//Conversion between Rad & Deg
Double_t MP = 0.938272;				//Proton Mass GeV
Double_t MD = 1.87561;				//GeV
Double_t MN = 0.939566;				//Neutron Mass GeV
Double_t me = 0.000510998;			//GeV


//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//-----------------------------------------------------------------------------Functions----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
void plotstuff(TH1F*, const char*, const char*, TFile*);
void be_and_af(TH1F*, TH1F*, TH1F*, TH1F*, double, double, double, double, double, double, double, double, const char*, const char*, bool, TFile*);
void be_and_af2d(TH2F*, TH2F*, TH2F*, TH2F*, const char*, const char*, TFile*);


//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//----------------------------------------------------------------------------Begin Main----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//comments to indicate each arguement is and usage
void omniSkims(string user_runtype ="", string user_target="", string user_pass="", string user_s_or_c="")
{
	Double_t nentries11 = 0;
	Double_t nentries22 = 0;
	Double_t temp1 = 0;
	
	//--------------------------------------	
	//Histogram Settings
	//--------------------------------------
	gStyle->SetOptStat(1);				//(1) for a stats window. (0) for no stats window.
	gStyle->SetOptTitle(1);				//(1) for a title. (0) for no title.
	gStyle->SetLabelSize(0.05, "xyz");	//Text size for the axis labels.
	gStyle->SetTitleSize(0.05, "xyz");	//Text size for the title.
	gStyle->SetPadBottomMargin(0.15);	//Margin spacing on the bottom of the histogram.
	gStyle->SetPadTopMargin(0.08);		//Margin spacing on the top of the histogram.
	gStyle->SetPadLeftMargin(0.17);		//Margin spacing on the left of the histogram.
	gStyle->SetPadRightMargin(0.2);		//Margin spacing on the right of the histogram.
	gStyle->SetImageScaling(3.);		//????.


	//--------------------------------------	
	//LH2, H(e,e'p) Target File Count & Run Numbers
	//--------------------------------------
	Double_t Heep_File_Count[1] = { 6 };
	const char* Heep_File[1][6] = {
		{"16962" "16963", "16964", "16965", "16966", "16967"}
	};


	//--------------------------------------	
	//MF Target File Counts & Run Numbers
	//--------------------------------------
	Double_t MF_File_Count[9] = { 1, 3, 2, 3, 2, 2, 3/*5*/, 2, 4/*5*/ }; //LD2, Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
	const char* MF_File[9][5] = {
		{ "16973", /*"16975"*/"", "", "", "" },								//LD2
		{ "16983", "17099", "17100", "", "" },								//Be9
		{ "16984", "17101", "", "", "" },									//B10
		{ "16986", "16991", "17102", "", "" },								//B11
		{ "16977", "17098", "", "", "" },									//C12
		{ "16980", "17097", "", "", "" },									//Ca40
		{ /*"16978", "16979",*/ "17093", "17094", "17096", "", "" },		//Ca48: 16978 & 16979 suffer from contamination are are chosen to be exluded
		{ "16981", "16982", "", "", "" },									//Fe54
		{ /*"20789",*/ "20793", "20797", "20798", "20799", "" }				//Au197: Carlos said something is wrong with 20789
	};

	
	//--------------------------------------
	//SRC Target File Counts & Run Numbers
	//--------------------------------------
	Double_t SRC_File_Count[9] = { 1, 10, 8, 8, 13, 25, 21, 24, 30 }; //LD2, Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
	const char* SRC_File[9][31] = {
		{ "17134", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },																//LD2
		{ "17106", "17107", "17108", "17109", "17110", "17111", "17129", "17130", "17131", "17132", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },					//Be9
		{ "17112", "17113", "17114", "17115", "17125", "17126", "17127", "17128", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },								//B10
		{ "17116", "17117", "17118", "17119", "17120", "17121", "17122", "17123", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },								//B11
		{ "17076", "17078", "17079", "17080", "17081", "17082", "17083", "17084", "17085", "17086", "17087", "17088", "17089", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },	//C12
		{ "17008", "17009", "17011", "17012", "17013", "17014", "17015", "17016", "17017", "17018", "17020", "17021", "17022",
		"17023", "17024", "17025", "17027", "17028", "17029", "17030", "17031", "17032", "17033", "17034", "17035", "", "", "", "", "", "" },															//Ca40
		{ "17036", "17037", "17038", "17039", "17040", "17041", "17043", "17044", "17045", "17046", "17047", "17048", "17049",
		"17050", "17051", "17052", "17053", "17054", "17055", "17056", "17057", "", "", "", "", "", "", "", "", "", "" },																				//Ca48
		{ "17058", "17059", "17060", "17061", "17062", "17063", "17064", "17067", "17068", "17069", "17070", "17071", "17072",
		"17073", "17074", "17075", "17135", "17136", "17138", "17140", "17141", "17142", "17143", "17144", "", "", "", "", "", "", "" },																//Fe54
		{ "20800", "20801", "20802", "20803", "20804", "20805", "20806", "20807", "20808", "20809", "20810", "20811", "20812",
		"20815", "20816", "20817", "20819", "20821", "20822", "20823", "20824", "20825", "20826", "20827", "20828", "20829", "20830", "20831", "20832", "20833", "20834" }								//Au197
	};
	

	//SRC run by run corrections
	Double_t Ca48_Corr[21] = { 0.975, 0.9789, 0.9822, 0.985, 0.9869, 0.9883, 0.9888, 0.9902, 0.9912, 0.9914, 0.9922, 0.9929, 0.9932, 0.9934, 0.994, 0.9943, 0.9945, 0.9947, 0.9948, 0.9949, 0.995 };
	//MF oil contamination corrections is 0.995

	//--------------------------------------	
	//Run Kinematics
	//--------------------------------------
	if (user_runtype == "")
	{
		cout << "Enter run type: SRC, MF, or Heep." << endl;
		cin >> user_runtype;
	}	

	//do not touch
	int binning = -1;
	if (user_runtype == "SRC") { binning = 3; }
	else if (user_runtype == "MF") { binning = 0; }
	else if (user_runtype == "Heep") { binning = 0; }
	else { cout << "Invalid run type." << endl; return; }


	//--------------------------------------	
	//Target List & Initial Target
	//--------------------------------------
	const char* Target[10] = { "LD2", "Be9", "B10", "B11", "C12", "Ca40", "Ca48", "Fe54", "Au197", "LH2" };

	if (user_target == "")
	{
		cout << "Enter Target: LD2, Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197, or LH2." << endl;
		cin >> user_target;
	}

	int target_id = -1;
	Double_t Z_per_A = 1;
	Double_t transparency = 1; //transparencies are calculated by Glauber Calculations
	Double_t areal_den;
	Double_t correction = 1;//oil contamination
	Double_t Prot_Abs = 0.952;//Proton Aborption.
	Double_t ctime_offset;
	Double_t Eff_total = 1;
	if (user_target == "LD2") { target_id = 0; Z_per_A = 1.0/2.0;  transparency = 0.7937; areal_den = 1.67; correction = 1; } 
	else if (user_target == "Be9") { target_id = 1; Z_per_A = 4.0/9.0;  transparency = 0.4807; areal_den = 0.9859; correction = 1; }
	else if (user_target == "B10") { target_id = 2; Z_per_A = 0.5;  transparency = 0.4642; areal_den = 0.44323; correction = 1; }
	else if (user_target == "B11") { target_id = 3; Z_per_A = 5.0/11.0;  transparency = 0.4496; areal_den = 0.4972; correction = 1; }
	else if (user_target == "C12") { target_id = 4; Z_per_A = 0.5;  transparency = 0.4368; areal_den = 0.5738; correction = 1; }
	else if (user_target == "Ca40") { target_id = 5; Z_per_A = 0.5;  transparency = 0.2924; areal_den = 0.7851; correction = 1.0; }//MF and SRC oil contamination is 0.99484 for all runs
	else if (user_target == "Ca48") { target_id = 6; Z_per_A = 20.0/48.0;  transparency = 0.2752; areal_den = 0.96157; correction = 1.0; }
	else if (user_target == "Fe54") { target_id = 7; Z_per_A = 26.0/54.0;  transparency = 0.2646; areal_den = 0.367; correction = 1; }
	else if (user_target == "Au197") { target_id = 8; Z_per_A = 79.0/197.0;  transparency = 0.1719; areal_den = 0.4047; correction = 1; }
	else if (user_target == "LH2") { target_id = 9; Z_per_A = 1.0;  transparency = 1; areal_den = 0.9859; correction = 1; }
	else { cout << "Invalid Target." << endl; }

	Double_t file_count = 0; //How many files a target has for the given kinematics
	if (user_runtype == "SRC") { file_count = SRC_File_Count[target_id]; }
	else if (user_runtype == "MF") { file_count = MF_File_Count[target_id]; }
	else if (user_runtype == "Heep") { file_count = Heep_File_Count[0]; }

	//ctime offset
	if (user_runtype == "MF")
	{
		if (user_target == "LD2") { ctime_offset = 0; Eff_total = 1;}
		else if (user_target == "Be9") { ctime_offset = 0; Eff_total = 0.845;}//time offset is 94.75;}
		else if (user_target == "B10") { ctime_offset = 0; Eff_total = 0.858;}//time offset is 94.75;}
		else if (user_target == "B11") { ctime_offset = 0; Eff_total = 0.854;}//time offset is 94.75;}
		else if (user_target == "C12") { ctime_offset = 0; Eff_total = 0.860;}//time offset is 94.75;}
		else if (user_target == "Ca40") { ctime_offset = 0; Eff_total = 0.843;}//time offset is 94.75;}
		else if (user_target == "Ca48") { ctime_offset = 0; Eff_total = 0.823;}//time offset is 86.583;}
		else if (user_target == "Fe54") { ctime_offset = 0; Eff_total = 0.877;}//time offset is 94.75;}
		else if (user_target == "Au197") { ctime_offset = 0; Eff_total = 0.826;}//time offset is 79.650;}
		else if (user_target == "LH2") { ctime_offset = 0; Eff_total = 1;}
		else { cout << "Invalid Target." << endl; }
	}
	else if (user_runtype == "SRC")
	{
		if (user_target == "LD2") { ctime_offset = 0; Eff_total = 1;}
		else if (user_target == "Be9") { ctime_offset = 0; Eff_total = 0.808;}//time offset is 87.062; }
		else if (user_target == "B10") { ctime_offset = 0; Eff_total = 0.842;}//time offset is 87.062; }
		else if (user_target == "B11") { ctime_offset = 0; Eff_total = 0.835;}//time offset is 87.062; }
		else if (user_target == "C12") { ctime_offset = 0; Eff_total = 0.849;}//time offset is 87.062; }
		else if (user_target == "Ca40") { ctime_offset = 0; Eff_total = 0.834;}//time offset is 87.062; }
		else if (user_target == "Ca48") { ctime_offset = 0; Eff_total = 0.816; correction = Ca48_Corr[0]; }//time offset is 87.062; }
		else if (user_target == "Fe54") { ctime_offset = 0; Eff_total = 0.857;}//time offset is 87.062; }
		else if (user_target == "Au197") { ctime_offset = 0; Eff_total = 0.827;}//time offset is 80.154; }
		else if (user_target == "LH2") { ctime_offset = 0; Eff_total = 1;}
		else { cout << "Invalid Target." << endl; }
	}

	
	//--------------------------------------	
	//Calibration Pass
	//--------------------------------------
	if (user_pass == "")
	{
		cout << "Enter Pass: pass1, pass2, pass3, or pass4." << endl;
		cin >> user_pass;
	}


	//--------------------------------------	
	//Singles events or Coincidence Events
	//--------------------------------------
	if (user_s_or_c == "")
	{
		cout << "Enter Singles (sing) or Coincidence (coin)." << endl;
		cin >> user_s_or_c;
	}


	
	//******************************************************************************************************************************************************************
	//---------------------------------------------------------------------------Manage FILES---------------------------------------------------------------------------
	//******************************************************************************************************************************************************************

	//--------------------------------------	
	//Find and open the first file
	//--------------------------------------
	TFile* inROOT;
	TTree* inputtree;
	
	string path1_c;
	string path2_c;
	string path3_c;
	string path4_c;
	string path4_c1;


	if (user_s_or_c == "sing")
	{
		path1_c  = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS1/CAFE_OUTPUT";
		path2_c  = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS2/CAFE_OUTPUT";
		path3_c  = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS3/CAFE_OUTPUT";
		path4_c  = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS4/CAFE_OUTPUT";
		path4_c1 = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS4";
	}
	else if (user_s_or_c == "coin")
	{
		path1_c  = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS1/CAFE_OUTPUT";
		path2_c  = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS2/CAFE_OUTPUT";
		path3_c  = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS3/CAFE_OUTPUT";
		path4_c  = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS4/CAFE_OUTPUT";
		path4_c1 = "/cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS4/CAFE_OUTPUT";
	}


	//../../../../../../../../../OFFLINE/PASS1/CAFE_OUTPUT/ROOT/cafe_prod_LD2_SRC_17134_-1_skimmed.root
	//_skimmed vs _histos
	
	if (user_pass == "pass1")
	{
		if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/skimmed_pass1/cafe_prod_%s_SRC_%s_-1_skimmed.root", path1_c.c_str(), Target[target_id], SRC_File[target_id][0]), "READ"); }
		else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/skimmed_pass1/cafe_prod_%s_MF_%s_-1_skimmed.root", path1_c.c_str(), Target[target_id], MF_File[target_id][0]), "READ"); }
		else if (user_runtype == "Heep") { inROOT = new TFile(Form("%s/skimmed_pass1/cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", path1_c.c_str(), Heep_File[0][0]), "READ"); }
	}
	else if (user_pass == "pass2")
	{
		if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_SRC_%s_-1_skimmed.root", path2_c.c_str(), Target[target_id], SRC_File[target_id][0]), "READ"); }
		else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_MF_%s_-1_skimmed.root", path2_c.c_str(), Target[target_id], MF_File[target_id][0]), "READ"); }
		else if (user_runtype == "Heep") { inROOT = new TFile(Form("../../../../../../../../OFFLINE/PASS2/ROOT/1cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", Heep_File[0][0]), "READ"); }
	}
	else if (user_pass == "pass3")
	{
		if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_SRC_%s_-1_skimmed.root", path3_c.c_str(), Target[target_id], SRC_File[target_id][0]), "READ"); }
		else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_MF_%s_-1_skimmed.root", path3_c.c_str(), Target[target_id], MF_File[target_id][0]), "READ"); }
		else if (user_runtype == "Heep") { inROOT = new TFile(Form("../../../../../../../../OFFLINE/PASS3/ROOT/cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", Heep_File[0][0]), "READ"); }
	}
	else if (user_pass == "pass4")
	{
		if (user_s_or_c == "coin")
		{
			if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_SRC_%s_-1_skimmed.root", path4_c1.c_str(), Target[target_id], SRC_File[target_id][0]), "READ"); }
			else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_MF_%s_-1_skimmed.root", path4_c1.c_str(), Target[target_id], MF_File[target_id][0]), "READ"); }
			else if (user_runtype == "Heep") { inROOT = new TFile(Form("../../../../../../../../OFFLINE/PASS4/ROOT/cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", Heep_File[0][0]), "READ"); }
		}
		else if (user_s_or_c == "sing")
		{
			if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/ROOTfiles/cafe_replay_prod_%s_-1.root", path4_c1.c_str(), SRC_File[target_id][0]), "READ"); }
			else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/ROOTfiles/cafe_replay_prod_%s_-1.root", path4_c1.c_str(), MF_File[target_id][0]), "READ"); }
			else if (user_runtype == "Heep") { inROOT = new TFile(Form("../../../../../../../../OFFLINE/PASS4/ROOT/cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", Heep_File[0][0]), "READ"); }
		}
	}
	else{ cout << "Invalid pass." << endl; }
	inputtree = (TTree*)inROOT->Get("T");
	
	
	inputtree->SetBranchStatus("*", kFALSE);

	inputtree->SetBranchStatus("shms_collimator_cut_flag", kTRUE);

	inputtree->SetBranchStatus("H.dc.ntrack", kTRUE);
	inputtree->SetBranchStatus("H.hod.goodscinhit", kTRUE);
	inputtree->SetBranchStatus("H.cer.npeSum", kTRUE);
	inputtree->SetBranchStatus("H.cal.etotnorm", kTRUE);
	inputtree->SetBranchStatus("H.cal.etottracknorm", kTRUE);
	inputtree->SetBranchStatus("H.hod.betanotrack", kTRUE);

	inputtree->SetBranchStatus("CTime.epCoinTime_ROC2_center", kTRUE);////or not center

	inputtree->SetBranchStatus("P.dc.ntrack", kTRUE);
	inputtree->SetBranchStatus("P.hod.goodscinhit", kTRUE);
	inputtree->SetBranchStatus("P.ngcer.npeSum", kTRUE);
	inputtree->SetBranchStatus("P.hgcer.npesum", kTRUE);
	inputtree->SetBranchStatus("P.cal.etotnorm", kTRUE);
	inputtree->SetBranchStatus("P.cal.etottracknorm", kTRUE);
	inputtree->SetBranchStatus("P.hod.betanotrack", kTRUE);

	inputtree->SetBranchStatus("P.kin.primary.scat_ang_rad", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.W", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.Q2", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.x_bj", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.nu", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.q3m", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.q_x", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.q_y", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.q_z", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.th_q", kTRUE);
	inputtree->SetBranchStatus("P.kin.primary.ph_q", kTRUE);

	if (user_target == "LH2") { inputtree->SetBranchStatus("H.kin.secondary.emiss", kTRUE); }//Standard Missing Energy for H(e,e'p)
	else { inputtree->SetBranchStatus("H.kin.secondary.emiss_nuc", kTRUE); }

	inputtree->SetBranchStatus("H.kin.secondary.pmiss", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.Prec_x", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.Prec_y", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.Prec_z", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.pmiss_x", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.pmiss_y", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.pmiss_z", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.tx", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.tb", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.Mrecoil", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.th_xq", kTRUE);
	//inputtree->SetBranchStatus("H.kin.secondary.th_bq", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.th_bq", kTRUE);//recoil particle in-plane angle w.r.to q-vector
	inputtree->SetBranchStatus("H.kin.secondary.ph_xq", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.ph_bq", kTRUE);
	inputtree->SetBranchStatus("H.kin.secondary.xangle", kTRUE);

	inputtree->SetBranchStatus("H.dc.x_fp", kTRUE);
	inputtree->SetBranchStatus("H.dc.xp_fp", kTRUE);
	inputtree->SetBranchStatus("H.dc.y_fp", kTRUE);
	inputtree->SetBranchStatus("H.dc.yp_fp", kTRUE);
	inputtree->SetBranchStatus("H.gtr.y", kTRUE);
	inputtree->SetBranchStatus("H.gtr.ph", kTRUE);
	inputtree->SetBranchStatus("H.gtr.th", kTRUE);
	inputtree->SetBranchStatus("H.gtr.dp", kTRUE);
	inputtree->SetBranchStatus("H.react.x", kTRUE);
	inputtree->SetBranchStatus("H.react.y", kTRUE);
	inputtree->SetBranchStatus("H.react.z", kTRUE);
	inputtree->SetBranchStatus("H.extcor.xsieve", kTRUE);
	inputtree->SetBranchStatus("H.extcor.ysieve", kTRUE);

	inputtree->SetBranchStatus("P.dc.x_fp", kTRUE);
	inputtree->SetBranchStatus("P.dc.xp_fp", kTRUE);
	inputtree->SetBranchStatus("P.dc.y_fp", kTRUE);
	inputtree->SetBranchStatus("P.dc.yp_fp", kTRUE);
	inputtree->SetBranchStatus("P.gtr.y", kTRUE);
	inputtree->SetBranchStatus("P.gtr.ph", kTRUE);
	inputtree->SetBranchStatus("P.gtr.p", kTRUE);
	inputtree->SetBranchStatus("H.gtr.p", kTRUE);
	inputtree->SetBranchStatus("P.gtr.th", kTRUE);
	inputtree->SetBranchStatus("P.gtr.dp", kTRUE);
	inputtree->SetBranchStatus("P.react.x", kTRUE);
	inputtree->SetBranchStatus("P.react.y", kTRUE);
	inputtree->SetBranchStatus("P.react.z", kTRUE);
	inputtree->SetBranchStatus("P.extcor.xsieve", kTRUE);
	inputtree->SetBranchStatus("P.extcor.ysieve", kTRUE);

	inputtree->SetBranchStatus("g.evtyp", kTRUE);
	inputtree->SetBranchStatus("P.dc.TheRealGolden", kTRUE);
	
	//Double_t shms_coll_cut = 0;
	//inputtree->SetBranchAddress("shms_collimator_cut_flag", &shms_coll_cut);

	//--------------------------------------	
	//Recreate a file to output to
	//--------------------------------------
	TFile* outROOT;
	if (user_runtype == "MF")
	{
		if (user_target == "LD2")		{ outROOT = new TFile("MF_D2_Results.root", "RECREATE"); }
		else if (user_target == "Be9")	{ outROOT = new TFile("MF_Be9_Results.root", "RECREATE"); }
		else if (user_target == "B10")	{ outROOT = new TFile("MF_B10_Results.root", "RECREATE"); }
		else if (user_target == "B11")	{ outROOT = new TFile("MF_B11_Results.root", "RECREATE"); }
		else if (user_target == "C12")	{ outROOT = new TFile("MF_C12_Results.root", "RECREATE"); }
		else if (user_target == "Ca40") { outROOT = new TFile("MF_Ca40_Results.root", "RECREATE"); }
		else if (user_target == "Ca48") { outROOT = new TFile("MF_Ca48_Results.root", "RECREATE"); }
		else if (user_target == "Fe54") { outROOT = new TFile("MF_Fe54_Results.root", "RECREATE"); }
		else if (user_target == "Au197") { outROOT = new TFile("MF_Au197_Results.root", "RECREATE"); }
	}
	else if (user_runtype == "SRC")
	{
		if (user_target == "LD2")		{ outROOT = new TFile("SRC_D2_Results.root", "RECREATE"); }
		else if (user_target == "Be9")	{ outROOT = new TFile("SRC_Be9_Results.root", "RECREATE"); }
		else if (user_target == "B10")	{ outROOT = new TFile("SRC_B10_Results.root", "RECREATE"); }
		else if (user_target == "B11")	{ outROOT = new TFile("SRC_B11_Results.root", "RECREATE"); }
		else if (user_target == "C12")	{ outROOT = new TFile("SRC_C12_Results.root", "RECREATE"); }
		else if (user_target == "Ca40") { outROOT = new TFile("SRC_Ca40_Results.root", "RECREATE"); }
		else if (user_target == "Ca48") { outROOT = new TFile("SRC_Ca48_Results.root", "RECREATE"); }
		else if (user_target == "Fe54") { outROOT = new TFile("SRC_Fe54_Results.root", "RECREATE"); }
		else if (user_target == "Au197") { outROOT = new TFile("SRC_Au197_Results.root", "RECREATE"); }
	}
	else if (user_runtype == "Heep")
	{
		if (user_target == "LH2")		{ outROOT = new TFile("Heep_Results.root", "RECREATE"); }
	}
	else { cout << "Invalid run type" << endl; }


	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	//------------------------------------------------------------------------DECLARE HISTOGRAMS------------------------------------------------------------------------
	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	TList* PIDList = new TList();		//Create TList to store histograms
	TList* Prim_Kin = new TList();		//Create TList to store histograms
	TList* Sec_Kin = new TList();		//Create TList to store histograms
	TList* HMS_Accp = new TList();		//Create TList to store histograms
	TList* SHMS_Accp = new TList();		//Create TList to store histograms
	TList* HList2 = new TList();		//Create TList to store histograms

	//double hhod_GoodSinHit; inputtree->SetBranchAddress("H.hod.goodscinhit", &hhod_GoodScinHit);
	//double hdc_ntrack; inputtree->SetBranchAddress("H.dc.ntrack", &hdc_ntrack);
	//double phod_GoodSinHit; inputtree->SetBranchAddress("P.hod.goodscinhit", &phod_GoodScinHit);
	//double pdc_ntrack; inputtree->SetBranchAddress("P.dc.ntrack", &pdc_ntrack);

	Double_t evtyp; inputtree->SetBranchAddress("g.evtyp", &evtyp);
	Double_t Real_Golden; inputtree->SetBranchAddress("P.dc.TheRealGolden", &Real_Golden);

	//Find offset for ep_ctime
	/*
	Double_t deadTime;  inputtree->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw", &deadTime);
	Double_t ep_ctime;	inputtree->SetBranchAddress("CTime.epCoinTime_ROC2", &ep_ctime);
	TH1F* H1_ctime_peak = new TH1F("H1_ctime_peak", "", 2000, -50, 150);//600
	int entries = inputtree->GetEntries();
	for (int i = 0; i < entries; i++) //50,000
	{
		inputtree->GetEntry(i);
		H1_ctime_peak->Fill(ep_ctime);
	}
	int binmax_ctime = H1_ctime_peak->GetMaximumBin();
	double ctime_offset = H1_ctime_peak->GetXaxis()->GetBinCenter(binmax_ctime);
	cout << "offset" << ctime_offset << endl;
	*/
	Double_t Ef;
	Double_t Pf;
	Double_t ep_ctime;				
	TH1F* H1_ep_ctime = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime", PIDList, binning);						
	TH1F* H1_ep_ctime_pid = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc_1 = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_1", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc_2 = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_2", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc_kin = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_kin", PIDList, binning);						
									
	TH1F* H1_ep_ctime_pid_alt = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_alt", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc_alt = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_alt", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc_kin_alt = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_kin_alt", PIDList, binning);						


	TH1F* H1_ep_ctime_rand = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_rand", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_rand = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_rand", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc_rand = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_rand", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc_kin_rand = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_kin_rand", PIDList, binning);						

	TH1F* H1_ep_ctime_sub = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_sub", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_sub = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_sub", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc_sub = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_sub", PIDList, binning);						
	TH1F* H1_ep_ctime_pid_acc_kin_sub = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_kin_sub", PIDList, binning);						

	//----- HMS -----
	Double_t hdc_ntrk; TH1F* H1_hdc_ntrk = Histo::mk_hdc_ntrk(inputtree, hdc_ntrk, "H1_hdc_ntrk", PIDList, binning);
	Double_t hScinGood; TH1F* H1_hScinGood = Histo::mk_hScinGood(inputtree, hScinGood, "H1_hScinGood", PIDList, binning);
									
	Double_t hCerNpeSum;			TH1F* H1_hCerNpeSum =			Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum", PIDList, binning);					
	Double_t hCalEtotNorm;			TH1F* H1_hCalEtotNorm =			Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm", PIDList, binning);			
	Double_t hCalEtotTrkNorm;		TH1F* H1_hCalEtotTrkNorm =		Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm", PIDList, binning);	
	Double_t hHodBetaNtrk;			TH1F* H1_hHodBetaNtrk =			Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk", PIDList, binning);			
	Double_t hHodBetaTrk;			TH1F* H1_hHodBetaTrk =			Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk", PIDList, binning);				
	
	TH1F* H1_hdc_ntrk_pid = Histo::mk_hdc_ntrk(inputtree, hdc_ntrk, "H1_hdc_ntrk_pid", PIDList, binning);
	TH1F* H1_hScinGood_pid = Histo::mk_hScinGood(inputtree, hScinGood, "H1_hScinGood_pid", PIDList, binning);

	TH1F* H1_hCerNpeSum_pid =			Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum_pid", PIDList, binning);					
	TH1F* H1_hCalEtotNorm_pid =		Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm_pid", PIDList, binning);			
	TH1F* H1_hCalEtotTrkNorm_pid =	Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm_pid", PIDList, binning);	
	TH1F* H1_hHodBetaNtrk_pid =		Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk_pid", PIDList, binning);			
	TH1F* H1_hHodBetaTrk_pid =		Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk_pid", PIDList, binning);				
	
	TH1F* H1_hdc_ntrk_pid_acc = Histo::mk_hdc_ntrk(inputtree, hdc_ntrk, "H1_hdc_ntrk_pid_acc", PIDList, binning);
	TH1F* H1_hScinGood_pid_acc = Histo::mk_hScinGood(inputtree, hScinGood, "H1_hScinGood_pid_acc", PIDList, binning);

	TH1F* H1_hCerNpeSum_pid_acc = Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum_pid_acc", PIDList, binning);					
	TH1F* H1_hCalEtotNorm_pid_acc = Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm_pid_acc", PIDList, binning);			
	TH1F* H1_hCalEtotTrkNorm_pid_acc = Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm_pid_acc", PIDList, binning);	
	TH1F* H1_hHodBetaNtrk_pid_acc = Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk_pid_acc", PIDList, binning);			
	TH1F* H1_hHodBetaTrk_pid_acc = Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk_pid_acc", PIDList, binning);				
	
	TH1F* H1_hdc_ntrk_pid_acc_1 = Histo::mk_hdc_ntrk(inputtree, hdc_ntrk, "H1_hdc_ntrk_pid_acc_1", PIDList, binning);
	TH1F* H1_hScinGood_pid_acc_1 = Histo::mk_hScinGood(inputtree, hScinGood, "H1_hScinGood_pid_acc_1", PIDList, binning);

	TH1F* H1_hCerNpeSum_pid_acc_1 = Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum_pid_acc_1", PIDList, binning);					
	TH1F* H1_hCalEtotNorm_pid_acc_1 = Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm_pid_acc_1", PIDList, binning);			
	TH1F* H1_hCalEtotTrkNorm_pid_acc_1 = Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm_pid_acc_1", PIDList, binning);	
	TH1F* H1_hHodBetaNtrk_pid_acc_1 = Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk_pid_acc_1", PIDList, binning);			
	TH1F* H1_hHodBetaTrk_pid_acc_1 = Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk_pid_acc_1", PIDList, binning);				

	TH1F* H1_hdc_ntrk_pid_acc_2 = Histo::mk_hdc_ntrk(inputtree, hdc_ntrk, "H1_hdc_ntrk_pid_acc_2", PIDList, binning);
	TH1F* H1_hScinGood_pid_acc_2 = Histo::mk_hScinGood(inputtree, hScinGood, "H1_hScinGood_pid_acc_2", PIDList, binning);

	TH1F* H1_hCerNpeSum_pid_acc_2 = Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum_pid_acc_2", PIDList, binning);					
	TH1F* H1_hCalEtotNorm_pid_acc_2 = Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm_pid_acc_2", PIDList, binning);			
	TH1F* H1_hCalEtotTrkNorm_pid_acc_2 = Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm_pid_acc_2", PIDList, binning);	
	TH1F* H1_hHodBetaNtrk_pid_acc_2 = Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk_pid_acc_2", PIDList, binning);			
	TH1F* H1_hHodBetaTrk_pid_acc_2 = Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk_pid_acc_2", PIDList, binning);				

	TH1F* H1_hdc_ntrk_pid_acc_kin = Histo::mk_hdc_ntrk(inputtree, hdc_ntrk, "H1_hdc_ntrk_pid_acc_kin", PIDList, binning);
	TH1F* H1_hScinGood_pid_acc_kin = Histo::mk_hScinGood(inputtree, hScinGood, "H1_hScinGood_pid_acc_kin", PIDList, binning);

	TH1F* H1_hCerNpeSum_pid_acc_kin = Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum_pid_acc_kin", PIDList, binning);					
	TH1F* H1_hCalEtotNorm_pid_acc_kin = Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm_pid_acc_kin", PIDList, binning);			
	TH1F* H1_hCalEtotTrkNorm_pid_acc_kin = Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm_pid_acc_kin", PIDList, binning);	
	TH1F* H1_hHodBetaNtrk_pid_acc_kin = Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk_pid_acc_kin", PIDList, binning);			
	TH1F* H1_hHodBetaTrk_pid_acc_kin = Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk_pid_acc_kin", PIDList, binning);				

	//----- SHMS -----
	Double_t pdc_ntrk; TH1F* H1_pdc_ntrk = Histo::mk_pdc_ntrk(inputtree, pdc_ntrk, "H1_pdc_ntrk", PIDList, binning);
	Double_t pScinGood; TH1F* H1_pScinGood = Histo::mk_pScinGood(inputtree, pScinGood, "H1_pScinGood", PIDList, binning);

	Double_t pNGCerNpeSum;			TH1F* H1_pNGCerNpeSum =			Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum", PIDList, binning);			
	Double_t pHGCerNpeSum;			TH1F* H1_pHGCerNpeSum =			Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum", PIDList, binning);			
	Double_t pCalEtotNorm;			TH1F* H1_pCalEtotNorm =			Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm", PIDList, binning);			
	Double_t pCalEtotTrkNorm;		TH1F* H1_pCalEtotTrkNorm =		Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm", PIDList, binning);	
	Double_t pHodBetaNtrk;			TH1F* H1_pHodBetaNtrk =			Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk", PIDList, binning);			
	Double_t pHodBetaTrk;			TH1F* H1_pHodBetaTrk =			Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk", PIDList, binning);				
	

	TH1F* H1_pdc_ntrk_pid = Histo::mk_pdc_ntrk(inputtree, pdc_ntrk, "H1_pdc_ntrk_pid", PIDList, binning);
	TH1F* H1_pScinGood_pid = Histo::mk_pScinGood(inputtree, pScinGood, "H1_pScinGood_pid", PIDList, binning);

	TH1F* H1_pNGCerNpeSum_pid =		Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum_pid", PIDList, binning);			
	TH1F* H1_pHGCerNpeSum_pid =		Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum_pid", PIDList, binning);			
	TH1F* H1_pCalEtotNorm_pid =		Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm_pid", PIDList, binning);			
	TH1F* H1_pCalEtotTrkNorm_pid =	Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_pid", PIDList, binning);	
	TH1F* H1_pHodBetaNtrk_pid =		Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk_pid", PIDList, binning);			
	TH1F* H1_pHodBetaTrk_pid =		Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk_pid", PIDList, binning);				
	

	TH1F* H1_pdc_ntrk_pid_acc = Histo::mk_pdc_ntrk(inputtree, pdc_ntrk, "H1_pdc_ntrk_pid_acc", PIDList, binning);
	TH1F* H1_pScinGood_pid_acc = Histo::mk_pScinGood(inputtree, pScinGood, "H1_pScinGood_pid_acc", PIDList, binning);

	TH1F* H1_pNGCerNpeSum_pid_acc = Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum_pid_acc", PIDList, binning);			
	TH1F* H1_pHGCerNpeSum_pid_acc = Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum_pid_acc", PIDList, binning);			
	TH1F* H1_pCalEtotNorm_pid_acc = Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm_pid_acc", PIDList, binning);			
	TH1F* H1_pCalEtotTrkNorm_pid_acc = Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_pid_acc", PIDList, binning);	
	TH1F* H1_pHodBetaNtrk_pid_acc = Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk_pid_acc", PIDList, binning);			
	TH1F* H1_pHodBetaTrk_pid_acc = Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk_pid_acc", PIDList, binning);				
	

	TH1F* H1_pdc_ntrk_pid_acc_1 = Histo::mk_pdc_ntrk(inputtree, pdc_ntrk, "H1_pdc_ntrk_pid_acc_1", PIDList, binning);
	TH1F* H1_pScinGood_pid_acc_1 = Histo::mk_pScinGood(inputtree, pScinGood, "H1_pScinGood_pid_acc_1", PIDList, binning);

	TH1F* H1_pNGCerNpeSum_pid_acc_1 = Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum_pid_acc_1", PIDList, binning);			
	TH1F* H1_pHGCerNpeSum_pid_acc_1 = Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum_pid_acc_1", PIDList, binning);			
	TH1F* H1_pCalEtotNorm_pid_acc_1 = Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm_pid_acc_1", PIDList, binning);			
	TH1F* H1_pCalEtotTrkNorm_pid_acc_1 = Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_pid_acc_1", PIDList, binning);	
	TH1F* H1_pHodBetaNtrk_pid_acc_1 = Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk_pid_acc_1", PIDList, binning);			
	TH1F* H1_pHodBetaTrk_pid_acc_1 = Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk_pid_acc_1", PIDList, binning);				


	TH1F* H1_pdc_ntrk_pid_acc_2 = Histo::mk_pdc_ntrk(inputtree, pdc_ntrk, "H1_pdc_ntrk_pid_acc_2", PIDList, binning);
	TH1F* H1_pScinGood_pid_acc_2 = Histo::mk_pScinGood(inputtree, pScinGood, "H1_pScinGood_pid_acc_2", PIDList, binning);

	TH1F* H1_pNGCerNpeSum_pid_acc_2 = Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum_pid_acc_2", PIDList, binning);			
	TH1F* H1_pHGCerNpeSum_pid_acc_2 = Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum_pid_acc_2", PIDList, binning);			
	TH1F* H1_pCalEtotNorm_pid_acc_2 = Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm_pid_acc_2", PIDList, binning);			
	TH1F* H1_pCalEtotTrkNorm_pid_acc_2 = Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_pid_acc_2", PIDList, binning);	
	TH1F* H1_pHodBetaNtrk_pid_acc_2 = Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk_pid_acc_2", PIDList, binning);			
	TH1F* H1_pHodBetaTrk_pid_acc_2 = Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk_pid_acc_2", PIDList, binning);				


	TH1F* H1_pdc_ntrk_pid_acc_kin = Histo::mk_pdc_ntrk(inputtree, pdc_ntrk, "H1_pdc_ntrk_pid_acc_kin", PIDList, binning);
	TH1F* H1_pScinGood_pid_acc_kin = Histo::mk_pScinGood(inputtree, pScinGood, "H1_pScinGood_pid_acc_kin", PIDList, binning);

	TH1F* H1_pNGCerNpeSum_pid_acc_kin = Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum_pid_acc_kin", PIDList, binning);			
	TH1F* H1_pHGCerNpeSum_pid_acc_kin = Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum_pid_acc_kin", PIDList, binning);			
	TH1F* H1_pCalEtotNorm_pid_acc_kin = Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm_pid_acc_kin", PIDList, binning);			
	TH1F* H1_pCalEtotTrkNorm_pid_acc_kin = Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_pid_acc_kin", PIDList, binning);	
	TH1F* H1_pHodBetaNtrk_pid_acc_kin = Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk_pid_acc_kin", PIDList, binning);			
	TH1F* H1_pHodBetaTrk_pid_acc_kin = Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk_pid_acc_kin", PIDList, binning);

	//******************************************************************************************************************************************************************
	//-------------------------------------------------Primary Kinematics (electron kinematics) (USED BY DATA AND SIMC)-------------------------------------------------
	//******************************************************************************************************************************************************************
	
	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	Double_t th_e;				TH1F* H1_th_e =						Histo::mk_the(inputtree, th_e, "H1_th_e", Prim_Kin, binning);			//Electron scattering angle
	Double_t W;					TH1F* H1_W =						Histo::mk_W(inputtree, W, "H1_W", Prim_Kin, binning);					//Invariant mass
	Double_t Q2;				TH1F* H1_Q2 =						Histo::mk_Q2(inputtree, Q2, "H1_Q2", Prim_Kin, binning);				//Four-momentum trasfer
	Double_t x_bj;				TH1F* H1_xbj =						Histo::mk_xbj(inputtree, x_bj, "H1_xbj", Prim_Kin, binning);			//B-jorken X  scaling variable
	Double_t nu;				TH1F* H1_nu =						Histo::mk_nu(inputtree, nu, "H1_nu", Prim_Kin, binning);				//Energy Transfer
	Double_t q;					TH1F* H1_q =						Histo::mk_q(inputtree, q, "H1_q", Prim_Kin, binning);					//Magnitude of the 3-vector q
	Double_t q_x;				TH1F* H1_qx =						Histo::mk_qx(inputtree, q_x, "H1_qx", Prim_Kin, binning);				//x-component of the energy transfer
	Double_t q_y;				TH1F* H1_qy =						Histo::mk_qy(inputtree, q_y, "H1_qy", Prim_Kin, binning);				//y-component of the energy transfer
	Double_t q_z;				TH1F* H1_qz =						Histo::mk_qz(inputtree, q_z, "H1_qz", Prim_Kin, binning);				//z-component of the energy transfer
	Double_t th_q;				TH1F* H1_th_q =						Histo::mk_thq(inputtree, th_q, "H1_th_q", Prim_Kin, binning);			//Angle between q and +z (hall coord. system)
	Double_t ph_q;				TH1F* H1_ph_q =						Histo::mk_phq(inputtree, ph_q, "H1_ph_q", Prim_Kin, binning);			//Out of plane angle between beamline and q
	//Double_t omega;				TH1F* H1_omega =					Histo::mk_omega(inputtree, omega, "H1_omega", Prim_Kin, binning);		//
	
	TH1F* H1_th_e_pid =			Histo::mk_the(inputtree, th_e, "H1_th_e_pid", Prim_Kin, binning);			//Electron scattering angle
	TH1F* H1_W_pid =			Histo::mk_W(inputtree, W, "H1_W_pid", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid =			Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid =			Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid =			Histo::mk_nu(inputtree, nu, "H1_nu_pid", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid =			Histo::mk_q(inputtree, q, "H1_q_pid", Prim_Kin, binning);				//Magnitude of the 3-vector q
	TH1F* H1_qx_pid =			Histo::mk_qx(inputtree, q_x, "H1_qx_pid", Prim_Kin, binning);				//x-component of the energy transfer
	TH1F* H1_qy_pid =			Histo::mk_qy(inputtree, q_y, "H1_qy_pid", Prim_Kin, binning);				//y-component of the energy transfer
	TH1F* H1_qz_pid =			Histo::mk_qz(inputtree, q_z, "H1_qz_pid", Prim_Kin, binning);				//z-component of the energy transfer
	TH1F* H1_th_q_pid =			Histo::mk_thq(inputtree, th_q, "H1_th_q_pid", Prim_Kin, binning);			//Angle between q and +z (hall coord. system)
	TH1F* H1_ph_q_pid =			Histo::mk_phq(inputtree, ph_q, "H1_ph_q_pid", Prim_Kin, binning);			//Out of plane angle between beamline and q
	//TH1F* H1_omega_pid =		Histo::mk_omega(inputtree, omega, "H1_omega_pid", Prim_Kin, binning);		//

	TH1F* H1_th_e_pid_acc =		Histo::mk_the(inputtree, th_e, "H1_th_e_pid_acc", Prim_Kin, binning);			//Electron scattering angle
	TH1F* H1_W_pid_acc =		Histo::mk_W(inputtree, W, "H1_W_pid_acc", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc =		Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_acc =		Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_acc =		Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc =		Histo::mk_q(inputtree, q, "H1_q_pid_acc", Prim_Kin, binning);				//Magnitude of the 3-vector q
	TH1F* H1_qx_pid_acc =		Histo::mk_qx(inputtree, q_x, "H1_qx_pid_acc", Prim_Kin, binning);				//x-component of the energy transfer
	TH1F* H1_qy_pid_acc =		Histo::mk_qy(inputtree, q_y, "H1_qy_pid_acc", Prim_Kin, binning);				//y-component of the energy transfer
	TH1F* H1_qz_pid_acc =		Histo::mk_qz(inputtree, q_z, "H1_qz_pid_acc", Prim_Kin, binning);				//z-component of the energy transfer
	TH1F* H1_th_q_pid_acc =		Histo::mk_thq(inputtree, th_q, "H1_th_q_pid_acc", Prim_Kin, binning);			//Angle between q and +z (hall coord. system)
	TH1F* H1_ph_q_pid_acc =		Histo::mk_phq(inputtree, ph_q, "H1_ph_q_pid_acc", Prim_Kin, binning);			//Out of plane angle between beamline and q
	//TH1F* H1_omega_pid_acc =	Histo::mk_omega(inputtree, omega, "H1_omega_pid_acc", Prim_Kin, binning);		//

	TH1F* H1_th_e_pid_acc_1 = Histo::mk_the(inputtree, th_e, "H1_th_e_pid_acc_1", Prim_Kin, binning);			//Electron scattering angle
	TH1F* H1_W_pid_acc_1 = Histo::mk_W(inputtree, W, "H1_W_pid_acc_1", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc_1 = Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc_1", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_acc_1 = Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc_1", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_acc_1 = Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc_1", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc_1 = Histo::mk_q(inputtree, q, "H1_q_pid_acc_1", Prim_Kin, binning);				//Magnitude of the 3-vector q
	TH1F* H1_qx_pid_acc_1 = Histo::mk_qx(inputtree, q_x, "H1_qx_pid_acc_1", Prim_Kin, binning);				//x-component of the energy transfer
	TH1F* H1_qy_pid_acc_1 = Histo::mk_qy(inputtree, q_y, "H1_qy_pid_acc_1", Prim_Kin, binning);				//y-component of the energy transfer
	TH1F* H1_qz_pid_acc_1 = Histo::mk_qz(inputtree, q_z, "H1_qz_pid_acc_1", Prim_Kin, binning);				//z-component of the energy transfer
	TH1F* H1_th_q_pid_acc_1 = Histo::mk_thq(inputtree, th_q, "H1_th_q_pid_acc_1", Prim_Kin, binning);			//Angle between q and +z (hall coord. system)
	TH1F* H1_ph_q_pid_acc_1 = Histo::mk_phq(inputtree, ph_q, "H1_ph_q_pid_acc_1", Prim_Kin, binning);			//Out of plane angle between beamline and q
	//TH1F* H1_omega_pid_acc_1 =	Histo::mk_omega(inputtree, omega, "H1_omega_pid_acc_1", Prim_Kin, binning);		//

	TH1F* H1_th_e_pid_acc_2 = Histo::mk_the(inputtree, th_e, "H1_th_e_pid_acc_2", Prim_Kin, binning);			//Electron scattering angle
	TH1F* H1_W_pid_acc_2 = Histo::mk_W(inputtree, W, "H1_W_pid_acc_2", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc_2 = Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc_2", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_acc_2 = Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc_2", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_acc_2 = Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc_2", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc_2 = Histo::mk_q(inputtree, q, "H1_q_pid_acc_2", Prim_Kin, binning);				//Magnitude of the 3-vector q
	TH1F* H1_qx_pid_acc_2 = Histo::mk_qx(inputtree, q_x, "H1_qx_pid_acc_2", Prim_Kin, binning);				//x-component of the energy transfer
	TH1F* H1_qy_pid_acc_2 = Histo::mk_qy(inputtree, q_y, "H1_qy_pid_acc_2", Prim_Kin, binning);				//y-component of the energy transfer
	TH1F* H1_qz_pid_acc_2 = Histo::mk_qz(inputtree, q_z, "H1_qz_pid_acc_2", Prim_Kin, binning);				//z-component of the energy transfer
	TH1F* H1_th_q_pid_acc_2 = Histo::mk_thq(inputtree, th_q, "H1_th_q_pid_acc_2", Prim_Kin, binning);			//Angle between q and +z (hall coord. system)
	TH1F* H1_ph_q_pid_acc_2 = Histo::mk_phq(inputtree, ph_q, "H1_ph_q_pid_acc_2", Prim_Kin, binning);			//Out of plane angle between beamline and q
	//TH1F* H1_omega_pid_acc_2 =	Histo::mk_omega(inputtree, omega, "H1_omega_pid_acc_2", Prim_Kin, binning);		//

	TH1F* H1_th_e_pid_acc_kin =	Histo::mk_the(inputtree, th_e, "H1_th_e_pid_acc_kin", Prim_Kin, binning);			//Electron scattering angle
	TH1F* H1_W_pid_acc_kin =	Histo::mk_W(inputtree, W, "H1_W_pid_acc_kin", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc_kin =	Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc_kin", Prim_Kin, binning);				//Four-momentum trasfer
	
	TH1F* H1_Q2_pid_acc_kin_alt = Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc_kin_alt", Prim_Kin, binning);				//Four-momentum trasfer
	
	TH1F* H1_xbj_pid_acc_kin =	Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc_kin", Prim_Kin, binning);			//B-jorken X  scaling variable
	
	TH1F* H1_xbj_pid_acc_kin_alt = Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc_kin_alt", Prim_Kin, binning);			//B-jorken X  scaling variable

	TH1F* H1_nu_pid_acc_kin =	Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc_kin", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc_kin =	Histo::mk_q(inputtree, q, "H1_q_pid_acc_kin", Prim_Kin, binning);				//Magnitude of the 3-vector q
	TH1F* H1_qx_pid_acc_kin =	Histo::mk_qx(inputtree, q_x, "H1_qx_pid_acc_kin", Prim_Kin, binning);				//x-component of the energy transfer
	TH1F* H1_qy_pid_acc_kin =	Histo::mk_qy(inputtree, q_y, "H1_qy_pid_acc_kin", Prim_Kin, binning);				//y-component of the energy transfer
	TH1F* H1_qz_pid_acc_kin =	Histo::mk_qz(inputtree, q_z, "H1_qz_pid_acc_kin", Prim_Kin, binning);				//z-component of the energy transfer
	TH1F* H1_th_q_pid_acc_kin =	Histo::mk_thq(inputtree, th_q, "H1_th_q_pid_acc_kin", Prim_Kin, binning);			//Angle between q and +z (hall coord. system)
	TH1F* H1_ph_q_pid_acc_kin = Histo::mk_phq(inputtree, ph_q, "H1_ph_q_pid_acc_kin", Prim_Kin, binning);			//Out of plane angle between beamline and q
	//TH1F* H1_omega_pid_acc_kin = Histo::mk_omega(inputtree, omega, "H1_omega_pid_acc_kin", Prim_Kin, binning);	//


	//***** Accidentals ****
	TH1F* H1_W_rand =						Histo::mk_W(inputtree, W, "H1_W_rand", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_rand =						Histo::mk_Q2(inputtree, Q2, "H1_Q2_rand", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_rand =						Histo::mk_xbj(inputtree, x_bj, "H1_xbj_rand", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_rand =						Histo::mk_nu(inputtree, nu, "H1_nu_rand", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_rand =						Histo::mk_q(inputtree, q, "H1_q_rand", Prim_Kin, binning);					//Magnitude of the 3-vector q
	
	TH1F* H1_W_pid_rand =			Histo::mk_W(inputtree, W, "H1_W_pid_rand", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_rand =			Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_rand", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_rand =			Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_rand", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_rand =			Histo::mk_nu(inputtree, nu, "H1_nu_pid_rand", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_rand =			Histo::mk_q(inputtree, q, "H1_q_pid_rand", Prim_Kin, binning);				//Magnitude of the 3-vector q
	
	TH1F* H1_W_pid_acc_rand =		Histo::mk_W(inputtree, W, "H1_W_pid_acc_rand", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc_rand =		Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc_rand", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_acc_rand =		Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc_rand", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_acc_rand =		Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc_rand", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc_rand =		Histo::mk_q(inputtree, q, "H1_q_pid_acc_rand", Prim_Kin, binning);				//Magnitude of the 3-vector q
	
	TH1F* H1_W_pid_acc_kin_rand =	Histo::mk_W(inputtree, W, "H1_W_pid_acc_kin_rand", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc_kin_rand =	Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc_kin_rand", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_acc_kin_rand =	Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc_kin_rand", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_acc_kin_rand =	Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc_kin_rand", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc_kin_rand =	Histo::mk_q(inputtree, q, "H1_q_pid_acc_kin_rand", Prim_Kin, binning);				//Magnitude of the 3-vector q

	//***** Subtraction ****
	TH1F* H1_W_sub =						Histo::mk_W(inputtree, W, "H1_W_sub", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_sub =						Histo::mk_Q2(inputtree, Q2, "H1_Q2_sub", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_sub =						Histo::mk_xbj(inputtree, x_bj, "H1_xbj_sub", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_sub =						Histo::mk_nu(inputtree, nu, "H1_nu_sub", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_sub =						Histo::mk_q(inputtree, q, "H1_q_sub", Prim_Kin, binning);					//Magnitude of the 3-vector q

	TH1F* H1_W_pid_sub =			Histo::mk_W(inputtree, W, "H1_W_pid_sub", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_sub =			Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_sub", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_sub =			Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_sub", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_sub =			Histo::mk_nu(inputtree, nu, "H1_nu_pid_sub", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_sub =			Histo::mk_q(inputtree, q, "H1_q_pid_sub", Prim_Kin, binning);				//Magnitude of the 3-vector q

	TH1F* H1_W_pid_acc_sub =		Histo::mk_W(inputtree, W, "H1_W_pid_acc_sub", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc_sub =		Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc_sub", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_acc_sub =		Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc_sub", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_acc_sub =		Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc_sub", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc_sub =		Histo::mk_q(inputtree, q, "H1_q_pid_acc_sub", Prim_Kin, binning);				//Magnitude of the 3-vector q

	
	TH1F* H1_W_pid_acc_kin_sub =	Histo::mk_W(inputtree, W, "H1_W_pid_acc_kin_sub", Prim_Kin, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc_kin_sub =	Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc_kin_sub", Prim_Kin, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_acc_kin_sub =	Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc_kin_sub", Prim_Kin, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_acc_kin_sub =	Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc_kin_sub", Prim_Kin, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc_kin_sub =	Histo::mk_q(inputtree, q, "H1_q_pid_acc_kin_sub", Prim_Kin, binning);				//Magnitude of the 3-vector q




	/*const Int_t NBins_ct = 27;
	//MF
	Double_t edges_Em[NBins_ct + 1] = {
		-0.02, -0.0145, -0.009, -0.0035, 0.002,
		0.0075, 0.013, 0.0185, 0.024, 0.0295,
		0.035, 0.0405, 0.046, 0.0515, 0.057,
		0.0625, 0.068, 0.0735, 0.079, 0.0845,
		0.09, 0.0955, 0.101, 0.1065, 0.112,
		0.1175, 0.1230, .1285//28
	};//pmiss*/

	/*
	const Int_t NBins_eh = 27;
	//MF
	Double_t edges_eh[NBins_eh + 1] = {
		-0.02, -0.0145, -0.009, -0.0035, 0.002,
		0.0075, 0.013, 0.0185, 0.024, 0.0295,
		0.035, 0.0405, 0.046, 0.0515, 0.057,
		0.0625, 0.068, 0.0735, 0.079, 0.0845,
		0.09, 0.0955, 0.101, 0.1065, 0.112,
		0.1175, 0.1230, .1285//28
	};//hdelta

	const Int_t NBins_ed = 27;
	//MF
	Double_t edges_ed[NBins_ed + 1] = {
		-0.02, -0.0145, -0.009, -0.0035, 0.002,
		0.0075, 0.013, 0.0185, 0.024, 0.0295,
		0.035, 0.0405, 0.046, 0.0515, 0.057,
		0.0625, 0.068, 0.0735, 0.079, 0.0845,
		0.09, 0.0955, 0.101, 0.1065, 0.112,
		0.1175, 0.1230, .1285//28
	};//edelta


	const Int_t NBins_pcal = 27;
	//MF
	Double_t edges_pcal[NBins_pcal + 1] = {
		-0.02, -0.0145, -0.009, -0.0035, 0.002,
		0.0075, 0.013, 0.0185, 0.024, 0.0295,
		0.035, 0.0405, 0.046, 0.0515, 0.057,
		0.0625, 0.068, 0.0735, 0.079, 0.0845,
		0.09, 0.0955, 0.101, 0.1065, 0.112,
		0.1175, 0.1230, .1285//28
	};//edelta*/


	const Int_t NBins_Em_mf = 27;
	//MF
	Double_t edges_Em_mf[NBins_Em_mf + 1] = {
		-0.02, -0.0145, -0.009, -0.0035, 0.002,
		0.0075, 0.013, 0.0185, 0.024, 0.0295,
		0.035, 0.0405, 0.046, 0.0515, 0.057,
		0.0625, 0.068, 0.0735, 0.079, 0.0845,
		0.09, 0.0955, 0.101, 0.1065, 0.112,
		0.1175, 0.1230, .1285//28
	};//Em


	const Int_t NBins_Em = 25;
	//MF
	Double_t edges_Em[NBins_Em + 1] = {
		0.0, 0.02, 0.04, 0.06, 0.08,
		0.1, 0.12, 0.14, 0.16, 0.18, 
		0.2, 0.22, 0.24, 0.26, 0.28,
		0.30, 0.32, 0.34, 0.36, 0.38,
		0.4, 0.42, 0.44, 0.46, 0.48,
		0.5//26
	};//Em


	//MF
	const Int_t NBins_Pm_mf = 32;
	Double_t edges_Pm_mf[NBins_Pm_mf + 1] = {
		0.0, 0.01, 0.02, 0.03, 0.04,
		0.05, 0.06, 0.07, 0.08, 0.09,
		0.10, 0.11, 0.12, 0.13, 0.14,
		0.15, 0.16, 0.17, 0.18, 0.19,
		0.20, 0.21, 0.22, 0.23, 0.24,
		0.25, 0.26, 0.27, 0.28, 0.29,
		0.30, 0.31, 0.32//33
	};//pmiss
	//SRC
	const Int_t NBins_Pm = 26;
	Double_t edges_Pm[NBins_Pm + 1] = {
		0.0, 0.0175, 0.05, 0.0825, 0.115, 
		0.1475, 0.18, 0.2125, 0.245, 0.2775, 
		0.31, 0.3425, 0.375, 0.4075, 0.44, 
		0.4725, 0.505, 0.5375, 0.57, 0.6025, 
		0.635, 0.6675, 0.700, 0.7325, 0.765,
		0.7975, 0.83//27
	};//pmiss
	

	const Int_t NBins_Q2 = 30;
	Double_t edges_Q2[NBins_Q2 + 1] = {
		0.0, 0.1, 0.2, 0.3, 0.4,
		0.5, 0.6, 0.7, 0.8, 0.9,
		1.0, 1.1, 1.2, 1.3, 1.4,
		1.5, 1.6, 1.7, 1.8, 1.9,
		2.0, 2.1, 2.2, 2.3, 2.4,
		2.5, 2.6, 2.7, 2.8, 2.9,
		3.0//31
	};//Q2


	const Int_t NBins_xbj = 20;
	Double_t edges_xbj[NBins_xbj + 1] = {
		0.0, 0.1, 0.2, 0.3, 0.4,
		0.5, 0.6, 0.7, 0.8, 0.9,
		1.0, 1.1, 1.2, 1.3, 1.4,
		1.5, 1.6, 1.7, 1.8, 1.9,
		2.0//32
	};//xbj


	const Int_t NBins_W = 22;//40
	Double_t edges_W[NBins_W + 1] = {
		0.05, 0.15, 0.25, 0.35, 0.45,
		0.55, 0.65, 0.75, 0.85, 0.95,
		1.05, 1.15, 1.25, 1.35, 1.45,
		1.55, 1.65, 1.75, 1.85, 1.95,
		2.05, 2.15, 2.25//23
	};//W


	const Int_t NBins_thrq = 35;//
	Double_t edges_thrq[NBins_thrq + 1] = {
		0.0, 2.0, 4.0, 6.0, 8.0,
		10.0, 12.0, 14.0, 16.0, 18.0,
		20.0, 22.0, 24.0, 26.0, 28.0,
		30.0, 32.0, 34.0, 36.0, 38.0,
		40.0, 42.0, 44.0, 46.0, 48.0,
		50.0, 52.0, 54.0, 56.0, 58.0,
		60.0, 62.0, 64.0, 66.0, 68.0,
		70.0//36
	};//thrq

	TH1F* H1_xbj_ratio = new TH1F("H1_xbj_ratio", "x_{bj}; x_{bj}", NBins_xbj, edges_xbj); Prim_Kin->Add(H1_xbj_ratio);
	TH1F* H1_xbj_ratio_pid = new TH1F("H1_xbj_ratio_pid", "x_{bj}; x_{bj}", NBins_xbj, edges_xbj); Prim_Kin->Add(H1_xbj_ratio_pid);
	TH1F* H1_xbj_ratio_pid_acc = new TH1F("H1_xbj_ratio_pid_acc", "x_{bj}; x_{bj}", NBins_xbj, edges_xbj); Prim_Kin->Add(H1_xbj_ratio_pid_acc);
	
	TH1F* H1_ep_ctime_ratio_pid_acc_kin_full = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_ratio_pid_acc_kin_full", PIDList, binning);
	TH1F* H1_ep_ctime_ratio_pid_acc_kin_full_mf = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_ratio_pid_acc_kin_full_mf", PIDList, binning);
	TH1F* H1_pCalEtotTrkNorm_ratio_pid_acc_kin_full = Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_ratio_pid_acc_kin_full", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_pCalEtotTrkNorm_ratio_pid_acc_kin_full_mf = Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_ratio_pid_acc_kin_full_mf", PIDList, binning);			//Four-momentum trasfer
	Double_t hdelta;
	Double_t edelta;
	TH1F* H1_hdelta_ratio_pid_acc_kin_full = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_ratio_pid_acc_kin_full", HMS_Accp, binning);
	TH1F* H1_hdelta_ratio_pid_acc_kin_full_mf = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_ratio_pid_acc_kin_full_mf", HMS_Accp, binning);
	TH1F* H1_edelta_ratio_pid_acc_kin_full = Histo::mk_edelta(inputtree, edelta, "H1_edelta_ratio_pid_acc_kin_full", SHMS_Accp, binning);
	TH1F* H1_edelta_ratio_pid_acc_kin_full_mf = Histo::mk_edelta(inputtree, edelta, "H1_edelta_ratio_pid_acc_kin_full_mf", SHMS_Accp, binning);

	TH2F* H2_hXColl_hYColl_ratio_pid_acc_kin_full = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_ratio_pid_acc_kin_full", HMS_Accp, binning);
	TH2F* H2_hXColl_hYColl_ratio_pid_acc_kin_full_mf = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_ratio_pid_acc_kin_full_mf", HMS_Accp, binning);
	TH2F* H2_eXColl_eYColl_ratio_pid_acc_kin_full = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_ratio_pid_acc_kin_full", SHMS_Accp, binning);
	TH2F* H2_eXColl_eYColl_ratio_pid_acc_kin_full_mf = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_ratio_pid_acc_kin_full_mf", SHMS_Accp, binning);
	TH2F* H2_Em_Pm_ratio_pid_acc_kin_full = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_ratio_pid_acc_kin_full", Prim_Kin, binning);
	TH2F* H2_Em_Pm_ratio_pid_acc_kin_full2 = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_ratio_pid_acc_kin_full2", Prim_Kin, binning);
	TH2F* H2_Em_Pm_ratio_pid_acc_kin_full_mf = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_ratio_pid_acc_kin_full_mf", Prim_Kin, binning);
	TH2F* H2_Em_Pm_ratio_pid_acc_kin_full_mf2 = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_ratio_pid_acc_kin_full_mf2", Prim_Kin, binning);
	TH2F* H2_xbj_Q2_ratio_pid_acc_kin_full = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_ratio_pid_acc_kin_full", Prim_Kin, binning);

	TH1F* H1_Q2_ratio_pid_acc_kin_full = new TH1F("H1_Q2_ratio_pid_acc_kin_full", "Q^{2}; Q^{2} (GeV/c)^{2}", NBins_Q2, edges_Q2); Prim_Kin->Add(H1_Q2_ratio_pid_acc_kin_full);
	TH1F* H1_xbj_ratio_pid_acc_kin_full = new TH1F("H1_xbj_ratio_pid_acc_kin_full", "x_{bj}; x_{bj}", NBins_xbj, edges_xbj); Prim_Kin->Add(H1_xbj_ratio_pid_acc_kin_full);
	TH1F* H1_thrq_ratio_pid_acc_kin_full = new TH1F("H1_thrq_ratio_pid_acc_kin_full", "#theta_{rq}; #theta_{rq} (deg)", NBins_thrq, edges_thrq); Prim_Kin->Add(H1_thrq_ratio_pid_acc_kin_full);
	TH1F* H1_W_ratio_pid_acc_kin_full = new TH1F("H1_W_ratio_pid_acc_kin_full", "W_{bj}; W_{bj} (GeV)", NBins_W, edges_W); Prim_Kin->Add(H1_W_ratio_pid_acc_kin_full);
	TH1F* H1_W_ratio_pid_acc_kin_full_mf = new TH1F("H1_W_ratio_pid_acc_kin_full_mf", "W_{bj}; W_{bj} (GeV)", NBins_W, edges_W); Prim_Kin->Add(H1_W_ratio_pid_acc_kin_full_mf);

	TH1F* H1_Q2_ratio_pid_acc_kin_full_mf = new TH1F("H1_Q2_ratio_pid_acc_kin_full_mf", "Q^{2}; Q^{2} (GeV/c)^{2}", NBins_Q2, edges_Q2); Prim_Kin->Add(H1_Q2_ratio_pid_acc_kin_full_mf);
	TH1F* H1_Em_ratio_pid_acc_kin_full_mf = new TH1F("H1_Em_ratio_pid_acc_kin_full_mf", "E_{miss}; E_{miss} (GeV)", NBins_Em_mf, edges_Em_mf); Prim_Kin->Add(H1_Em_ratio_pid_acc_kin_full_mf);
	TH1F* H1_Em_ratio_pid_acc_kin_full = new TH1F("H1_Em_ratio_pid_acc_kin_full", "E_{miss}; E_{miss} (GeV)", NBins_Em, edges_Em); Prim_Kin->Add(H1_Em_ratio_pid_acc_kin_full);




	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	//Double_t Ki;			TH1F* H1_Ki = Histo::mk_ki(inputtree, Ki, "H1_Ki", Sec_Kin);
	//Double_t Kf;			TH1F* H1_Kf = Histo::mk_kf(inputtree, Kf, "H1_Kf", Sec_Kin);
	

	//******************************************************************************************************************************************************************
	//--------------------------------------------Secondary (Hadron) Kinematics (recoil and missing are used interchageably)--------------------------------------------
	//******************************************************************************************************************************************************************
	/*const Int_t NBins2 = 62;
	Double_t edges2[NBins2 + 1] = { 
		0.0, 0.01, 0.02, 0.03, 0.04, 
		0.05, 0.06, 0.07, 0.08, 0.09, 
		0.10, 0.11, 0.12, 0.13, 0.14, 
		0.15, 0.16, 0.17, 0.18, 0.19, 
		0.20, 0.21, 0.22, 0.23, 0.24, 
		0.25, 0.26, 0.27, 0.375, 0.385, 
		0.395, 0.405, 0.415, 0.425, 0.435, 
		0.445, 0.455, 0.465, 0.475, 0.485, 
		0.495, 0.505, 0.515, 0.525, 0.535, 
		0.545, 0.555, 0.565, 0.575, 0.585, 
		0.595, 0.605, 0.615, 0.625, 0.635, 
		0.645, 0.655, 0.665, 0.675, 0.685, 
		0.695, 0.700, 0.71//63
		};//momentum bin ranges*/


	/*const Int_t NBins2 = 41;
	Double_t edges2[NBins2 + 1] = {
		0.0, 0.01, 0.02, 0.03, 0.04,
		0.05, 0.06, 0.07, 0.08, 0.09,
		0.10, 0.11, 0.12, 0.13, 0.14,
		0.15, 0.16, 0.17, 0.18, 0.19,
		0.20, 0.21, 0.22, 0.23, 0.24,
		0.25, 0.26, 0.27, 0.375, 0.395,
		0.415, 0.445, 0.475, 0.505, 0.535,
		0.565, 0.595, 0.625, 0.655, 0.685,
		0.700, 0.71//42
	};//momentum bin ranges*/

	const Int_t NBins2 = 39;
	Double_t edges2[NBins2 + 1] = {
		0.0, 0.01, 0.02, 0.03, 0.04,
		0.05, 0.06, 0.07, 0.08, 0.09,
		0.10, 0.11, 0.12, 0.13, 0.14,
		0.15, 0.16, 0.17, 0.18, 0.19,
		0.20, 0.21, 0.22, 0.23, 0.24,
		0.25, 0.26, 0.27, 0.375, 0.4075,
		0.44, 0.4725, 0.505, 0.5375, 0.57,
		0.6025, 0.635, 0.6675, 0.700, 0.71//40
	};//momentum bin ranges

	/*const Int_t NBins2 = 41;
	Double_t edges2[NBins2 + 1] = {
		0.0, 0.01, 0.02, 0.03, 0.04,
		0.05, 0.06, 0.07, 0.08, 0.09,
		0.10, 0.11, 0.12, 0.13, 0.14,
		0.15, 0.16, 0.17, 0.18, 0.19,
		0.20, 0.21, 0.22, 0.23, 0.24,
		0.25, 0.26, 0.27, 0.375, 0.40208,
		0.42917, 0.45625, 0.48333, 0.51042, 0.5375,
		0.56458, 0.59167, 0.61875, 0.64583, 0.67292,
		0.7, 0.71//42
	};//momentum bin ranges*/

	/*const Int_t NBins2 = 10;
	Double_t edges2[NBins2 + 1] = {
		0.3, 0.375, 0.415625, 0.45625, 0.496875,
		0.5375, 0.578125, 0.61875, 0.659375, 0.7,
		0.75//11
	};//momentum bin ranges*/

	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	Double_t emiss;
	
	TH1F* H1_emiss;
	if (user_target == "LH2") { H1_emiss = Histo::mk_Em(inputtree, emiss, "H1_emiss", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else{ H1_emiss = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss", Sec_Kin, binning); }

	Double_t pmiss;			TH1F* H1_pmiss =			Histo::mk_Pm(inputtree, pmiss, "H1_pmiss", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw =		Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c", Sec_Kin, binning);
	TH1F* H1_pmiss_ratio = new TH1F("H1_pmiss_ratio", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio);
	TH1F* H1_pmiss_ratio_pid = new TH1F("H1_pmiss_ratio_pid", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid);
	TH1F* H1_pmiss_ratio_pid_acc = new TH1F("H1_pmiss_ratio_pid_acc", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc);
	TH1F* H1_pmiss_ratio_pid_acc_kin = new TH1F("H1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);
	TH1F* H1_pmiss_ratio_pid_acc_kin_full = new TH1F("H1_pmiss_ratio_pid_acc_kin_full", "p_{miss}; p_{miss} (GeV/c)", NBins_Pm, edges_Pm); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin_full);
	TH1F* H1_pmiss_ratio_pid_acc_kin_full_mf = new TH1F("H1_pmiss_ratio_pid_acc_kin_full_mf", "p_{miss}; p_{miss} (GeV/c)", NBins_Pm_mf, edges_Pm_mf); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin_full_mf);
	Double_t prec_x;		TH1F* H1_prec_x =			Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x", Sec_Kin, binning);
	Double_t prec_y;		TH1F* H1_prec_y =			Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y", Sec_Kin, binning); 
	Double_t prec_z;		TH1F* H1_prec_z =			Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z", Sec_Kin, binning); 
	Double_t pmiss_x;		TH1F* H1_pmiss_x =			Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x", Sec_Kin, binning); 
	Double_t pmiss_y;		TH1F* H1_pmiss_y =			Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y", Sec_Kin, binning); 
	Double_t pmiss_z;		TH1F* H1_pmiss_z =			Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z", Sec_Kin, binning); 
	Double_t Tx;			TH1F* H1_Tx =				Histo::mk_Tx(inputtree, Tx, "H1_Tx", Sec_Kin, binning);
	Double_t Tr;			TH1F* H1_Tr =				Histo::mk_Tr(inputtree, Tr, "H1_Tr", Sec_Kin, binning);
	Double_t mmiss;			TH1F* H1_mmiss =			Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss", Sec_Kin, binning);
	Double_t th_pq;			TH1F* H1_th_pq =			Histo::mk_thpq(inputtree, th_pq, "H1_th_pq", Sec_Kin, binning); 		//detected particle in-plane angle w.r.to q-vector
	Double_t cth_rq;			TH1F* H1_th_rq =			Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq", Sec_Kin, binning); 		//recoil particle in-plane angle w.r.to q-vector
	/*Double_t cth_rq;*/		TH1F* H1_cth_rq =			Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq", Sec_Kin, binning); 		//recoil particle in-plane angle w.r.to q-vector
	Double_t ph_pq;			TH1F* H1_ph_pq =			Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq", Sec_Kin, binning); 		//detected particle ???
	Double_t ph_rq;			TH1F* H1_ph_rq =			Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq", Sec_Kin, binning); 
	Double_t xangle;		TH1F* H1_xangle =			Histo::mk_xangle(inputtree, xangle, "H1_xangle", Sec_Kin, binning); 
	Double_t Pm_par;		TH1F* H1_Pm_par =			Histo::mk_Pm_par(inputtree, Pm_par, "H1_Pm_par", Sec_Kin, binning);
	Double_t Pm_per;		TH1F* H1_Pm_per =			Histo::mk_Pm_per(inputtree, Pm_per, "H1_Pm_per", Sec_Kin, binning);
	Double_t alpha_n_v;
	Double_t alpha_v;		TH1F* H1_alpha =			Histo::mk_alpha(inputtree, alpha_v, "H1_alpha", Sec_Kin, binning);


	TH1F* H1_emiss_pid;
	if (user_target == "LH2") { H1_emiss_pid = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid", Sec_Kin, binning); }

	TH1F* H1_pmiss_pid =			Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid", Sec_Kin, binning);							//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_pid =		Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_prec_x_pid =			Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x_pid", Sec_Kin, binning);
	TH1F* H1_prec_y_pid =			Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y_pid", Sec_Kin, binning);
	TH1F* H1_prec_z_pid =			Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z_pid", Sec_Kin, binning);
	TH1F* H1_pmiss_x_pid =			Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x_pid", Sec_Kin, binning);
	TH1F* H1_pmiss_y_pid =			Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y_pid", Sec_Kin, binning);
	TH1F* H1_pmiss_z_pid =			Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z_pid", Sec_Kin, binning);
	TH1F* H1_Tx_pid =				Histo::mk_Tx(inputtree, Tx, "H1_Tx_pid", Sec_Kin, binning);
	TH1F* H1_Tr_pid =				Histo::mk_Tr(inputtree, Tr, "H1_Tr_pid", Sec_Kin, binning);
	TH1F* H1_mmiss_pid =			Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss_pid", Sec_Kin, binning);
	TH1F* H1_th_pq_pid =			Histo::mk_thpq(inputtree, th_pq, "H1_th_pq_pid", Sec_Kin, binning); 						//detected particle in-plane angle w.r.to q-vector
	TH1F* H1_th_rq_pid =			Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid", Sec_Kin, binning); 						//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_cth_rq_pid =			Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq_pid", Sec_Kin, binning); 						//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_ph_pq_pid =			Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq_pid", Sec_Kin, binning); 						//detected particle ???
	TH1F* H1_ph_rq_pid =			Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq_pid", Sec_Kin, binning);
	TH1F* H1_xangle_pid =			Histo::mk_xangle(inputtree, xangle, "H1_xangle_pid", Sec_Kin, binning);
	TH1F* H1_Pm_par_pid =				Histo::mk_Pm_par(inputtree, Pm_par, "H1_Pm_par_pid", Sec_Kin, binning);
	TH1F* H1_Pm_per_pid =				Histo::mk_Pm_per(inputtree, Pm_per, "H1_Pm_per_pid", Sec_Kin, binning);
	TH1F* H1_alpha_pid =				Histo::mk_alpha(inputtree, alpha_v, "H1_alpha_pid", Sec_Kin, binning);


	TH1F* H1_emiss_pid_acc;
	if (user_target == "LH2") { H1_emiss_pid_acc = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_acc = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_acc", Sec_Kin, binning); }
	
	TH1F* H1_emiss_pid_acc_1;
	if (user_target == "LH2") { H1_emiss_pid_acc_1 = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc_1", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_acc_1 = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_acc_1", Sec_Kin, binning); }

	TH1F* H1_emiss_pid_acc_2;
	if (user_target == "LH2") { H1_emiss_pid_acc_2 = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc_2", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_acc_2 = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_acc_2", Sec_Kin, binning); }

	TH1F* H1_pmiss_pid_acc =		Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc", Sec_Kin, binning);						//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_pid_acc =	Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_acc", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_acc = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_acc", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_prec_x_pid_acc =		Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x_pid_acc", Sec_Kin, binning);
	TH1F* H1_prec_y_pid_acc =		Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y_pid_acc", Sec_Kin, binning);
	TH1F* H1_prec_z_pid_acc =		Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z_pid_acc", Sec_Kin, binning);
	TH1F* H1_pmiss_x_pid_acc =		Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x_pid_acc", Sec_Kin, binning);
	TH1F* H1_pmiss_y_pid_acc =		Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y_pid_acc", Sec_Kin, binning);
	TH1F* H1_pmiss_z_pid_acc =		Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z_pid_acc", Sec_Kin, binning);
	TH1F* H1_Tx_pid_acc =			Histo::mk_Tx(inputtree, Tx, "H1_Tx_pid_acc", Sec_Kin, binning);
	TH1F* H1_Tr_pid_acc =			Histo::mk_Tr(inputtree, Tr, "H1_Tr_pid_acc", Sec_Kin, binning);
	TH1F* H1_mmiss_pid_acc =		Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss_pid_acc", Sec_Kin, binning);
	TH1F* H1_th_pq_pid_acc =		Histo::mk_thpq(inputtree, th_pq, "H1_th_pq_pid_acc", Sec_Kin, binning); 					//detected particle in-plane angle w.r.to q-vector
	TH1F* H1_th_rq_pid_acc =		Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid_acc", Sec_Kin, binning); 					//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_cth_rq_pid_acc =		Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq_pid_acc", Sec_Kin, binning); 					//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_ph_pq_pid_acc =		Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq_pid_acc", Sec_Kin, binning); 					//detected particle ???
	TH1F* H1_ph_rq_pid_acc =		Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq_pid_acc", Sec_Kin, binning);
	TH1F* H1_xangle_pid_acc =		Histo::mk_xangle(inputtree, xangle, "H1_xangle_pid_acc", Sec_Kin, binning);
	TH1F* H1_Pm_par_pid_acc =		Histo::mk_Pm_par(inputtree, Pm_par, "H1_Pm_par_pid_acc", Sec_Kin, binning);
	TH1F* H1_Pm_per_pid_acc =		Histo::mk_Pm_per(inputtree, Pm_per, "H1_Pm_pe_pid_acc", Sec_Kin, binning);
	TH1F* H1_alpha_pid_acc =		Histo::mk_alpha(inputtree, alpha_v, "H1_alpha_pid_acc", Sec_Kin, binning);

	TH1F* H1_pmiss_pid_acc_1 =		Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc_1", Sec_Kin, binning);						//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_pid_acc_1 =	Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_acc_1", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_acc_1 = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_acc_1", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_prec_x_pid_acc_1 =		Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_prec_y_pid_acc_1 =		Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_prec_z_pid_acc_1 =		Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_pmiss_x_pid_acc_1 =	Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_pmiss_y_pid_acc_1 =	Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_pmiss_z_pid_acc_1 =	Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_Tx_pid_acc_1 =			Histo::mk_Tx(inputtree, Tx, "H1_Tx_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_Tr_pid_acc_1 =			Histo::mk_Tr(inputtree, Tr, "H1_Tr_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_mmiss_pid_acc_1 =		Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_th_pq_pid_acc_1 =		Histo::mk_thpq(inputtree, th_pq, "H1_th_pq_pid_acc_1", Sec_Kin, binning); 					//detected particle in-plane angle w.r.to q-vector
	TH1F* H1_th_rq_pid_acc_1 =		Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid_acc_1", Sec_Kin, binning); 					//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_cth_rq_pid_acc_1 =		Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq_pid_acc_1", Sec_Kin, binning); 					//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_ph_pq_pid_acc_1 =		Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq_pid_acc_1", Sec_Kin, binning); 					//detected particle ???
	TH1F* H1_ph_rq_pid_acc_1 =		Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_xangle_pid_acc_1 =		Histo::mk_xangle(inputtree, xangle, "H1_xangle_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_Pm_par_pid_acc_1 =		Histo::mk_Pm_par(inputtree, Pm_par, "H1_Pm_par_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_Pm_per_pid_acc_1 =		Histo::mk_Pm_per(inputtree, Pm_per, "H1_Pm_pe_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_alpha_pid_acc_1 =		Histo::mk_alpha(inputtree, alpha_v, "H1_alpha_pid_acc_1", Sec_Kin, binning);

	TH1F* H1_pmiss_pid_acc_2 =		Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc_2", Sec_Kin, binning);						//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_pid_acc_2 =	Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_acc_2", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_acc_2 = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_acc_2", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_prec_x_pid_acc_2 =		Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_prec_y_pid_acc_2 =		Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_prec_z_pid_acc_2 =		Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_pmiss_x_pid_acc_2 =	Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_pmiss_y_pid_acc_2 =	Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_pmiss_z_pid_acc_2 =	Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_Tx_pid_acc_2 =			Histo::mk_Tx(inputtree, Tx, "H1_Tx_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_Tr_pid_acc_2 =			Histo::mk_Tr(inputtree, Tr, "H1_Tr_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_mmiss_pid_acc_2 =		Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_th_pq_pid_acc_2 =		Histo::mk_thpq(inputtree, th_pq, "H1_th_pq_pid_acc_2", Sec_Kin, binning); 					//detected particle in-plane angle w.r.to q-vector
	TH1F* H1_th_rq_pid_acc_2 =		Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid_acc_2", Sec_Kin, binning); 					//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_cth_rq_pid_acc_2 =		Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq_pid_acc_2", Sec_Kin, binning); 					//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_ph_pq_pid_acc_2 =		Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq_pid_acc_2", Sec_Kin, binning); 					//detected particle ???
	TH1F* H1_ph_rq_pid_acc_2 =		Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_xangle_pid_acc_2 =		Histo::mk_xangle(inputtree, xangle, "H1_xangle_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_Pm_par_pid_acc_2 =		Histo::mk_Pm_par(inputtree, Pm_par, "H1_Pm_par_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_Pm_per_pid_acc_2 =		Histo::mk_Pm_per(inputtree, Pm_per, "H1_Pm_pe_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_alpha_pid_acc_2 =		Histo::mk_alpha(inputtree, alpha_v, "H1_alpha_pid_acc_2", Sec_Kin, binning);

	TH1F* H1_emiss_pid_acc_kin;
	if (user_target == "LH2") { H1_emiss_pid_acc_kin = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc_kin", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_acc_kin = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_acc_kin", Sec_Kin, binning); }
	
	TH1F* H1_emiss_pid_acc_kin_alt;
	if (user_target == "LH2") { H1_emiss_pid_acc_kin_alt = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc_kin_alt", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_acc_kin_alt = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_acc_kin_alt", Sec_Kin, binning); }


	TH1F* H1_pmiss_pid_acc_kin =	Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc_kin", Sec_Kin, binning);					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	
	TH1F* H1_Ef_pid_acc_kin = Histo::mk_kf(inputtree, Ef, "H1_Ef_pid_acc_kin", Sec_Kin, binning);//verify					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_Pf_pid_acc_kin = Histo::mk_Pf(inputtree, Pf, "H1_Pf_pid_acc_kin", Sec_Kin, binning);//verify					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_Ef_pid_acc_kin_mf = Histo::mk_kf(inputtree, Ef, "H1_Ef_pid_acc_kin_mf", Sec_Kin, binning);//verify					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_Pf_pid_acc_kin_mf = Histo::mk_Pf(inputtree, Pf, "H1_Pf_pid_acc_kin_mf", Sec_Kin, binning);//verify					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))

	TH1F* H1_pmiss_raw_pid_acc_kin = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_acc_kin", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_acc_kin = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_acc_kin", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))

	TH1F* H1_pmiss_pid_acc_kin_alt =	Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc_kin_alt", Sec_Kin, binning);					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	
	TH1F* H1_prec_x_pid_acc_kin =	Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_prec_y_pid_acc_kin =	Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_prec_z_pid_acc_kin =	Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_pmiss_x_pid_acc_kin =	Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_pmiss_y_pid_acc_kin =	Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_pmiss_z_pid_acc_kin =	Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_Tx_pid_acc_kin =		Histo::mk_Tx(inputtree, Tx, "H1_Tx_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_Tr_pid_acc_kin =		Histo::mk_Tr(inputtree, Tr, "H1_Tr_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_mmiss_pid_acc_kin =	Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_th_pq_pid_acc_kin =	Histo::mk_thpq(inputtree, th_pq, "H1_th_pq_pid_acc_kin", Sec_Kin, binning); 				//detected particle in-plane angle w.r.to q-vector
	TH1F* H1_th_rq_pid_acc_kin =	Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid_acc_kin", Sec_Kin, binning); 				//recoil particle in-plane angle w.r.to q-vector
	
	TH1F* H1_th_rq_pid_acc_kin_alt = Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid_acc_kin_alt", Sec_Kin, binning); 				//recoil particle in-plane angle w.r.to q-vector

	TH1F* H1_cth_rq_pid_acc_kin =	Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq_pid_acc_kin", Sec_Kin, binning); 				//recoil particle in-plane angle w.r.to q-vector
	
	//TH1F* H1_cth_rq_pid_acc_kin_alt = Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq_pid_acc_kin_alt", Sec_Kin, binning); 				//recoil particle in-plane angle w.r.to q-vector

	TH1F* H1_ph_pq_pid_acc_kin =	Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq_pid_acc_kin", Sec_Kin, binning); 				//detected particle ???
	TH1F* H1_ph_rq_pid_acc_kin =	Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_xangle_pid_acc_kin =	Histo::mk_xangle(inputtree, xangle, "H1_xangle_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_Pm_par_pid_acc_kin =	Histo::mk_Pm_par(inputtree, Pm_par, "H1_Pm_par_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_Pm_per_pid_acc_kin =	Histo::mk_Pm_per(inputtree, Pm_per, "H1_Pm_pe_pid_acc_kin", Sec_Kin, binning);
	TH1F* H1_alpha_pid_acc_kin =	Histo::mk_alpha(inputtree, alpha_v, "H1_alpha_pid_acc_kin", Sec_Kin, binning);
	//mk_alpha_Pm_per(TTree* input, const char* name, TList* list, int bin)






	//***** Accidentals *****
	TH1F* H1_pmiss_rand =			Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_rand", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_rand = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_rand", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_rand = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_rand", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_th_rq_rand =			Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_rand", Sec_Kin, binning); 		//recoil particle in-plane angle w.r.to q-vector
	
	TH1F* H1_pmiss_pid_rand =			Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_rand", Sec_Kin, binning);							//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_pid_rand = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_rand", Sec_Kin, binning);							//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_rand = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_rand", Sec_Kin, binning);							//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_th_rq_pid_rand =			Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid_rand", Sec_Kin, binning); 						//recoil particle in-plane angle w.r.to q-vector
	
	TH1F* H1_pmiss_pid_acc_rand =		Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc_rand", Sec_Kin, binning);						//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_pid_acc_rand = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_acc_rand", Sec_Kin, binning);						//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_acc_rand = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_acc_rand", Sec_Kin, binning);						//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_th_rq_pid_acc_rand =		Histo::mk_thrq(inputtree, cth_rq, "H1_th_pq_pid_acc_rand", Sec_Kin, binning); 					//detected particle in-plane angle w.r.to q-vector
	
	TH1F* H1_pmiss_pid_acc_kin_rand =	Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc_kin_rand", Sec_Kin, binning);					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_ratio_pid_acc_kin_rand = new TH1F("H1_pmiss_ratio_pid_acc_kin_rand", "p_{miss}; p_{miss} (GeV/c)", NBins2 , edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin_rand);
	TH1F* H1_pmiss_ratio_pid_acc_kin_rand_full = new TH1F("H1_pmiss_ratio_pid_acc_kin_rand_full", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin_rand_full);
	//TH1F* H1_pmiss_ratio_pid_acc_kin_rand_full_mf = new TH1F("H1_pmiss_ratio_pid_acc_kin_rand_full_mf", "p_{miss}; p_{miss} (GeV/c)", NBins2_Pm, edges2_Pm); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin_mf);
	TH1F* H1_Q2_ratio_pid_acc_kin_rand_full = new TH1F("H1_Q3_ratio_pid_acc_kin_rand_full", "Q^{2}; Q^{2} (GeV/c)^{2}", NBins2, edges2); Sec_Kin->Add(H1_Q2_ratio_pid_acc_kin_rand_full);
	TH1F* H1_xbj_ratio_pid_acc_kin_rand_full = new TH1F("H1_xbj_ratio_pid_acc_kin_rand_full", "x_{bj}; x_{bj}", NBins2, edges2); Sec_Kin->Add(H1_xbj_ratio_pid_acc_kin_rand_full);
	TH1F* H1_thrq_ratio_pid_acc_kin_rand_full = new TH1F("H1_thrq_ratio_pid_acc_kin_rand_full", "#theta_{rq}; #theta_{rq} (deg)", NBins2, edges2); Sec_Kin->Add(H1_thrq_ratio_pid_acc_kin_rand_full);
	TH1F* H1_W_ratio_pid_acc_kin_rand_full = new TH1F("H1_W_W_pid_acc_kin_rand_full", "W; W (GeV)", NBins2, edges2); Sec_Kin->Add(H1_W_ratio_pid_acc_kin_rand_full);
	TH1F* H1_pmiss_raw_pid_acc_kin_rand = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_acc_kin_rand", Sec_Kin, binning);					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_acc_kin_rand = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_acc_kin_rand", Sec_Kin, binning);					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_th_rq_pid_acc_kin_rand =	Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid_acc_kin_rand", Sec_Kin, binning); 				//recoil particle in-plane angle w.r.to q-vector
	
	TH1F* H1_emiss_rand;
	if (user_target == "LH2") { H1_emiss_rand = Histo::mk_Em(inputtree, emiss, "H1_emiss_rand", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_rand = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_rand", Sec_Kin, binning); }
	
	TH1F* H1_emiss_pid_rand;
	if (user_target == "LH2") { H1_emiss_pid_rand = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_rand", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_rand = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_rand", Sec_Kin, binning); }

	TH1F* H1_emiss_pid_acc_rand;
	if (user_target == "LH2") { H1_emiss_pid_acc_rand = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc_rand", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_acc_rand = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_acc_rand", Sec_Kin, binning); }

	TH1F* H1_emiss_pid_acc_kin_rand;
	if (user_target == "LH2") { H1_emiss_pid_acc_kin_rand = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc_kin_rand", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_acc_kin_rand = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_acc_kin_rand", Sec_Kin, binning); }



	//***** Subtraction *****
	TH1F* H1_pmiss_sub =			Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_sub", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_sub = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_sub", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_sub = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_sub", Sec_Kin, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_th_rq_sub =			Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_sub", Sec_Kin, binning); 		//recoil particle in-plane angle w.r.to q-vector

	TH1F* H1_pmiss_pid_sub =			Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_sub", Sec_Kin, binning);							//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_pid_sub = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_sub", Sec_Kin, binning);							//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_sub = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_sub", Sec_Kin, binning);							//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_th_rq_pid_sub =			Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid_sub", Sec_Kin, binning); 						//recoil particle in-plane angle w.r.to q-vector

	TH1F* H1_pmiss_pid_acc_sub =		Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc_sub", Sec_Kin, binning);						//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_pid_acc_sub = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_acc_sub", Sec_Kin, binning);						//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_acc_sub = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_acc_sub", Sec_Kin, binning);						//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_th_rq_pid_acc_sub =		Histo::mk_thrq(inputtree, cth_rq, "H1_th_pq_pid_acc_sub", Sec_Kin, binning); 					//detected particle in-plane angle w.r.to q-vector

	TH1F* H1_pmiss_pid_acc_kin_sub =	Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc_kin_sub", Sec_Kin, binning);					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_raw_pid_acc_kin_sub = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_raw_pid_acc_kin_sub", Sec_Kin, binning);					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_pmiss_c_pid_acc_kin_sub = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_c_pid_acc_kin_sub", Sec_Kin, binning);					//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_th_rq_pid_acc_kin_sub =	Histo::mk_thrq(inputtree, cth_rq, "H1_th_rq_pid_acc_kin_sub", Sec_Kin, binning); 				//recoil particle in-plane angle w.r.to q-vector
	
	TH1F* H1_emiss_sub;
	if (user_target == "LH2") { H1_emiss_sub = Histo::mk_Em(inputtree, emiss, "H1_emiss_sub", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_sub = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_sub", Sec_Kin, binning); }
	
	TH1F* H1_emiss_pid_sub;
	if (user_target == "LH2") { H1_emiss_pid_sub = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_sub", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_sub = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_sub", Sec_Kin, binning); }

	TH1F* H1_emiss_pid_acc_sub;
	if (user_target == "LH2") { H1_emiss_pid_acc_sub = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc_sub", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_acc_sub = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_acc_sub", Sec_Kin, binning); }

	TH1F* H1_emiss_pid_acc_kin_sub;
	if (user_target == "LH2") { H1_emiss_pid_acc_kin_sub = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc_kin_sub", Sec_Kin, binning); }		//Standard Missing Energy for H(e,e'p)
	else { H1_emiss_pid_acc_kin_sub = Histo::mk_Em_nuc(inputtree, emiss, "H1_emiss_pid_acc_kin_sub", Sec_Kin, binning); }















	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	//Double_t Pf;			TH1F* H1_Pf =				Histo::mk_Pf(inputtree, Pf, "H1_Pf", Sec_Kin, binning);
	
	Double_t th_p;			TH1F* H1_th_p = Histo::mk_thp(inputtree, th_p, "H1_th_p", Sec_Kin, binning);
	TH1F* H1_th_p_pid = Histo::mk_thp(inputtree, th_p, "H1_th_p_pid", Sec_Kin, binning);
	TH1F* H1_th_p_pid_acc = Histo::mk_thp(inputtree, th_p, "H1_th_p_pid_acc", Sec_Kin, binning);
	TH1F* H1_th_p_pid_acc_1 = Histo::mk_thp(inputtree, th_p, "H1_th_p_pid_acc_1", Sec_Kin, binning);
	TH1F* H1_th_p_pid_acc_2 = Histo::mk_thp(inputtree, th_p, "H1_th_p_pid_acc_2", Sec_Kin, binning);
	TH1F* H1_th_p_pid_acc_kin = Histo::mk_thp(inputtree, th_p, "H1_th_p_pid_acc_kin", Sec_Kin, binning);

	Double_t Em_src;		TH1F* H1_Em_src = Histo::mk_Em_src(inputtree, Em_src, "H1_Em_src", Sec_Kin, binning);

	

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//--------------------------------------------------------------------Acceptance Histogram Bins---------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//----------------------------------------
	//--------------- HMS Leafs --------------
	//----------------------------------------
	//----- Hadron Arm Focal Plane -----
	Double_t hxfp;			TH1F* H1_hxfp =						Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp", HMS_Accp, binning); 			
	Double_t hxpfp;			TH1F* H1_hxpfp =					Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp", HMS_Accp, binning); 		
	Double_t hyfp;			TH1F* H1_hyfp =						Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp", HMS_Accp, binning); 			
	Double_t hypfp;			TH1F* H1_hypfp =					Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp", HMS_Accp, binning); 		
	
	TH1F* H1_hxfp_pid = Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp_pid", HMS_Accp, binning); 			
	TH1F* H1_hxpfp_pid = Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp_pid", HMS_Accp, binning); 		
	TH1F* H1_hyfp_pid = Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp_pid", HMS_Accp, binning); 			
	TH1F* H1_hypfp_pid = Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp_pid", HMS_Accp, binning); 		
	
	TH1F* H1_hxfp_pid_acc = Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp_pid_acc", HMS_Accp, binning); 			
	TH1F* H1_hxpfp_pid_acc = Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp_pid_acc", HMS_Accp, binning); 		
	TH1F* H1_hyfp_pid_acc = Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp_pid_acc", HMS_Accp, binning); 			
	TH1F* H1_hypfp_pid_acc = Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp_pid_acc", HMS_Accp, binning); 		
	
	TH1F* H1_hxfp_pid_acc_1 = Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp_pid_acc_1", HMS_Accp, binning);
	TH1F* H1_hxpfp_pid_acc_1 = Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp_pid_acc_1", HMS_Accp, binning);
	TH1F* H1_hyfp_pid_acc_1 = Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp_pid_acc_1", HMS_Accp, binning);
	TH1F* H1_hypfp_pid_acc_1 = Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp_pid_acc_1", HMS_Accp, binning);

	TH1F* H1_hxfp_pid_acc_2 = Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp_pid_acc_2", HMS_Accp, binning);
	TH1F* H1_hxpfp_pid_acc_2 = Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp_pid_acc_2", HMS_Accp, binning);
	TH1F* H1_hyfp_pid_acc_2 = Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp_pid_acc_2", HMS_Accp, binning);
	TH1F* H1_hypfp_pid_acc_2 = Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp_pid_acc_2", HMS_Accp, binning);

	TH1F* H1_hxfp_pid_acc_kin = Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp_pid_acc_kin", HMS_Accp, binning); 			
	TH1F* H1_hxpfp_pid_acc_kin = Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp_pid_acc_kin", HMS_Accp, binning); 		
	TH1F* H1_hyfp_pid_acc_kin = Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp_pid_acc_kin", HMS_Accp, binning); 			
	TH1F* H1_hypfp_pid_acc_kin = Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp_pid_acc_kin", HMS_Accp, binning); 		

	//----- Hadron Arm Reconstructed Quantities -----
	Double_t hytar;			TH1F* H1_hytar =					Histo::mk_hytar(inputtree, hytar, "H1_hytar", HMS_Accp, binning); 		
	Double_t hyptar;		TH1F* H1_hyptar =					Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar", HMS_Accp, binning); 	
	Double_t hxptar;		TH1F* H1_hxptar =					Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar", HMS_Accp, binning); 	
	TH1F* H1_hdelta =					Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta", HMS_Accp, binning); 	
	
	TH1F* H1_hytar_pid = Histo::mk_hytar(inputtree, hytar, "H1_hytar_pid", HMS_Accp, binning); 		
	TH1F* H1_hyptar_pid = Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar_pid", HMS_Accp, binning); 	
	TH1F* H1_hxptar_pid = Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar_pid", HMS_Accp, binning); 	
	TH1F* H1_hdelta_pid = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_pid", HMS_Accp, binning); 	
	
	TH1F* H1_hytar_pid_acc = Histo::mk_hytar(inputtree, hytar, "H1_hytar_pid_acc", HMS_Accp, binning); 		
	TH1F* H1_hyptar_pid_acc = Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar_pid_acc", HMS_Accp, binning); 	
	TH1F* H1_hxptar_pid_acc = Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar_pid_acc", HMS_Accp, binning); 	
	TH1F* H1_hdelta_pid_acc = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_pid_acc", HMS_Accp, binning); 	
	
	TH1F* H1_hytar_pid_acc_1 = Histo::mk_hytar(inputtree, hytar, "H1_hytar_pid_acc_1", HMS_Accp, binning);
	TH1F* H1_hyptar_pid_acc_1 = Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar_pid_acc_1", HMS_Accp, binning);
	TH1F* H1_hxptar_pid_acc_1 = Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar_pid_acc_1", HMS_Accp, binning);
	TH1F* H1_hdelta_pid_acc_1 = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_pid_acc_1", HMS_Accp, binning);

	TH1F* H1_hytar_pid_acc_2 = Histo::mk_hytar(inputtree, hytar, "H1_hytar_pid_acc_2", HMS_Accp, binning);
	TH1F* H1_hyptar_pid_acc_2 = Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar_pid_acc_2", HMS_Accp, binning);
	TH1F* H1_hxptar_pid_acc_2 = Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar_pid_acc_2", HMS_Accp, binning);
	TH1F* H1_hdelta_pid_acc_2 = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_pid_acc_2", HMS_Accp, binning);

	TH1F* H1_hytar_pid_acc_kin = Histo::mk_hytar(inputtree, hytar, "H1_hytar_pid_acc_kin", HMS_Accp, binning); 		
	TH1F* H1_hyptar_pid_acc_kin = Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar_pid_acc_kin", HMS_Accp, binning); 	
	TH1F* H1_hxptar_pid_acc_kin = Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar_pid_acc_kin", HMS_Accp, binning); 	
	TH1F* H1_hdelta_pid_acc_kin = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_pid_acc_kin", HMS_Accp, binning); 	
	
	//----- Target Reconstruction (Hall Coord. System) -----
	Double_t htarx;			TH1F* H1_htarx =					Histo::mk_htarx(inputtree, htarx, "H1_htarx", HMS_Accp, binning); 		
	Double_t htary;			TH1F* H1_htary =					Histo::mk_htary(inputtree, htary, "H1_htary", HMS_Accp, binning); 		
	Double_t htarz;			TH1F* H1_htarz =					Histo::mk_htarz(inputtree, htarz, "H1_htarz", HMS_Accp, binning); 		
	
	TH1F* H1_htarx_pid = Histo::mk_htarx(inputtree, htarx, "H1_htarx_pid", HMS_Accp, binning); 		
	TH1F* H1_htary_pid = Histo::mk_htary(inputtree, htary, "H1_htary_pid", HMS_Accp, binning); 		
	TH1F* H1_htarz_pid = Histo::mk_htarz(inputtree, htarz, "H1_htarz_pid", HMS_Accp, binning); 		
	
	TH1F* H1_htarx_pid_acc = Histo::mk_htarx(inputtree, htarx, "H1_htarx_pid_acc", HMS_Accp, binning); 		
	TH1F* H1_htary_pid_acc = Histo::mk_htary(inputtree, htary, "H1_htary_pid_acc", HMS_Accp, binning); 		
	TH1F* H1_htarz_pid_acc = Histo::mk_htarz(inputtree, htarz, "H1_htarz_pid_acc", HMS_Accp, binning); 		
	
	TH1F* H1_htarx_pid_acc_1 = Histo::mk_htarx(inputtree, htarx, "H1_htarx_pid_acc_1", HMS_Accp, binning);
	TH1F* H1_htary_pid_acc_1 = Histo::mk_htary(inputtree, htary, "H1_htary_pid_acc_1", HMS_Accp, binning);
	TH1F* H1_htarz_pid_acc_1 = Histo::mk_htarz(inputtree, htarz, "H1_htarz_pid_acc_1", HMS_Accp, binning);

	TH1F* H1_htarx_pid_acc_2 = Histo::mk_htarx(inputtree, htarx, "H1_htarx_pid_acc_2", HMS_Accp, binning);
	TH1F* H1_htary_pid_acc_2 = Histo::mk_htary(inputtree, htary, "H1_htary_pid_acc_2", HMS_Accp, binning);
	TH1F* H1_htarz_pid_acc_2 = Histo::mk_htarz(inputtree, htarz, "H1_htarz_pid_acc_2", HMS_Accp, binning);

	TH1F* H1_htarx_pid_acc_kin = Histo::mk_htarx(inputtree, htarx, "H1_htarx_pid_acc_kin", HMS_Accp, binning); 		
	TH1F* H1_htary_pid_acc_kin = Histo::mk_htary(inputtree, htary, "H1_htary_pid_acc_kin", HMS_Accp, binning); 		
	TH1F* H1_htarz_pid_acc_kin = Histo::mk_htarz(inputtree, htarz, "H1_htarz_pid_acc_kin", HMS_Accp, binning); 		
	
	//----- HMS Collimator -----
	Double_t hXColl;		TH1F* H1_hXColl =					Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl", HMS_Accp, binning); 	
	Double_t hYColl;		TH1F* H1_hYColl =					Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl", HMS_Accp, binning); 	
	
	TH1F* H1_hXColl_pid = Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl_pid", HMS_Accp, binning); 	
	TH1F* H1_hYColl_pid = Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl_pid", HMS_Accp, binning); 	
	
	TH1F* H1_hXColl_pid_acc = Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl_pid_acc", HMS_Accp, binning); 	
	TH1F* H1_hYColl_pid_acc = Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl_pid_acc", HMS_Accp, binning); 	
	
	TH1F* H1_hXColl_pid_acc_1 = Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl_pid_acc_1", HMS_Accp, binning);
	TH1F* H1_hYColl_pid_acc_1 = Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl_pid_acc_1", HMS_Accp, binning);

	TH1F* H1_hXColl_pid_acc_2 = Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl_pid_acc_2", HMS_Accp, binning);
	TH1F* H1_hYColl_pid_acc_2 = Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl_pid_acc_2", HMS_Accp, binning);

	TH1F* H1_hXColl_pid_acc_kin = Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl_pid_acc_kin", HMS_Accp, binning); 	
	TH1F* H1_hYColl_pid_acc_kin = Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl_pid_acc_kin", HMS_Accp, binning); 	


	//----------------------------------------
	//-------------- SHMS Leafs --------------
	//----------------------------------------
	//----- Electron Arm Focal Plane -----
	Double_t exfp;			TH1F* H1_exfp =						Histo::mk_exfp(inputtree, exfp, "H1_exfp", SHMS_Accp, binning);			
	Double_t expfp;			TH1F* H1_expfp =					Histo::mk_expfp(inputtree, expfp, "H1_expfp", SHMS_Accp, binning);		
	Double_t eyfp;			TH1F* H1_eyfp =						Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp", SHMS_Accp, binning);			
	Double_t eypfp;			TH1F* H1_eypfp =					Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp", SHMS_Accp, binning);		
	
	TH1F* H1_exfp_pid =		Histo::mk_exfp(inputtree, exfp, "H1_exfp_pid", SHMS_Accp, binning);			
	TH1F* H1_expfp_pid =		Histo::mk_expfp(inputtree, expfp, "H1_expfp_pid", SHMS_Accp, binning);		
	TH1F* H1_eyfp_pid =		Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp_pid", SHMS_Accp, binning);			
	TH1F* H1_eypfp_pid =		Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp_pid", SHMS_Accp, binning);		
	
	TH1F* H1_exfp_pid_acc = Histo::mk_exfp(inputtree, exfp, "H1_exfp_pid_acc", SHMS_Accp, binning);			
	TH1F* H1_expfp_pid_acc = Histo::mk_expfp(inputtree, expfp, "H1_expfp_pid_acc", SHMS_Accp, binning);		
	TH1F* H1_eyfp_pid_acc = Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp_pid_acc", SHMS_Accp, binning);			
	TH1F* H1_eypfp_pid_acc = Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp_pid_acc", SHMS_Accp, binning);		
	
	TH1F* H1_exfp_pid_acc_1 = Histo::mk_exfp(inputtree, exfp, "H1_exfp_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_expfp_pid_acc_1 = Histo::mk_expfp(inputtree, expfp, "H1_expfp_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_eyfp_pid_acc_1 = Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_eypfp_pid_acc_1 = Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp_pid_acc_1", SHMS_Accp, binning);

	TH1F* H1_exfp_pid_acc_2 = Histo::mk_exfp(inputtree, exfp, "H1_exfp_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_expfp_pid_acc_2 = Histo::mk_expfp(inputtree, expfp, "H1_expfp_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_eyfp_pid_acc_2 = Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_eypfp_pid_acc_2 = Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp_pid_acc_2", SHMS_Accp, binning);

	TH1F* H1_exfp_pid_acc_kin = Histo::mk_exfp(inputtree, exfp, "H1_exfp_pid_acc_kin", SHMS_Accp, binning);			
	TH1F* H1_expfp_pid_acc_kin = Histo::mk_expfp(inputtree, expfp, "H1_expfp_pid_acc_kin", SHMS_Accp, binning);		
	TH1F* H1_eyfp_pid_acc_kin = Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp_pid_acc_kin", SHMS_Accp, binning);			
	TH1F* H1_eypfp_pid_acc_kin = Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp_pid_acc_kin", SHMS_Accp, binning);		

	//----- Electron Arm Reconstructed Quantities -----
	Double_t eytar;			TH1F* H1_eytar =					Histo::mk_eytar(inputtree, eytar, "H1_eytar", SHMS_Accp, binning);		
	Double_t eyptar;		TH1F* H1_eyptar =					Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar", SHMS_Accp, binning);		
	Double_t exptar;		TH1F* H1_exptar =					Histo::mk_exptar(inputtree, exptar, "H1_exptar", SHMS_Accp, binning);		
	TH1F* H1_edelta =					Histo::mk_edelta(inputtree, edelta, "H1_edelta", SHMS_Accp, binning);		
	
	TH1F* H1_eytar_pid =		Histo::mk_eytar(inputtree, eytar, "H1_eytar_pid", SHMS_Accp, binning);		
	TH1F* H1_eyptar_pid =		Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar_pid", SHMS_Accp, binning);		
	TH1F* H1_exptar_pid =		Histo::mk_exptar(inputtree, exptar, "H1_exptar_pid", SHMS_Accp, binning);		
	TH1F* H1_edelta_pid =		Histo::mk_edelta(inputtree, edelta, "H1_edelta_pid", SHMS_Accp, binning);		
	
	TH1F* H1_eytar_pid_acc = Histo::mk_eytar(inputtree, eytar, "H1_eytar_pid_acc", SHMS_Accp, binning);		
	TH1F* H1_eyptar_pid_acc = Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar_pid_acc", SHMS_Accp, binning);		
	TH1F* H1_exptar_pid_acc = Histo::mk_exptar(inputtree, exptar, "H1_exptar_pid_acc", SHMS_Accp, binning);		
	TH1F* H1_edelta_pid_acc = Histo::mk_edelta(inputtree, edelta, "H1_edelta_pid_acc", SHMS_Accp, binning);		
	
	TH1F* H1_eytar_pid_acc_1 = Histo::mk_eytar(inputtree, eytar, "H1_eytar_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_eyptar_pid_acc_1 = Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_exptar_pid_acc_1 = Histo::mk_exptar(inputtree, exptar, "H1_exptar_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_edelta_pid_acc_1 = Histo::mk_edelta(inputtree, edelta, "H1_edelta_pid_acc_1", SHMS_Accp, binning);

	TH1F* H1_eytar_pid_acc_2 = Histo::mk_eytar(inputtree, eytar, "H1_eytar_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_eyptar_pid_acc_2 = Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_exptar_pid_acc_2 = Histo::mk_exptar(inputtree, exptar, "H1_exptar_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_edelta_pid_acc_2 = Histo::mk_edelta(inputtree, edelta, "H1_edelta_pid_acc_2", SHMS_Accp, binning);

	TH1F* H1_eytar_pid_acc_kin = Histo::mk_eytar(inputtree, eytar, "H1_eytar_pid_acc_kin", SHMS_Accp, binning);		
	TH1F* H1_eyptar_pid_acc_kin = Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar_pid_acc_kin", SHMS_Accp, binning);		
	TH1F* H1_exptar_pid_acc_kin = Histo::mk_exptar(inputtree, exptar, "H1_exptar_pid_acc_kin", SHMS_Accp, binning);		
	TH1F* H1_edelta_pid_acc_kin = Histo::mk_edelta(inputtree, edelta, "H1_edelta_pid_acc_kin", SHMS_Accp, binning);		
	
	//----- Target Reconstruction (Hall Coord. System) -----
	Double_t etarx;			TH1F* H1_etarx =					Histo::mk_etarx(inputtree, etarx, "H1_etarx", SHMS_Accp, binning);		
	Double_t etary;			TH1F* H1_etary =					Histo::mk_etary(inputtree, etary, "H1_etary", SHMS_Accp, binning);		
	Double_t etarz;			TH1F* H1_etarz =					Histo::mk_etarz(inputtree, etarz, "H1_etarz", SHMS_Accp, binning);		
	
	TH1F* H1_etarx_pid =		Histo::mk_etarx(inputtree, etarx, "H1_etarx_pid", SHMS_Accp, binning);		
	TH1F* H1_etary_pid =		Histo::mk_etary(inputtree, etary, "H1_etary_pid", SHMS_Accp, binning);		
	TH1F* H1_etarz_pid =		Histo::mk_etarz(inputtree, etarz, "H1_etarz_pid", SHMS_Accp, binning);		
	
	TH1F* H1_etarx_pid_acc = Histo::mk_etarx(inputtree, etarx, "H1_etarx_pid_acc", SHMS_Accp, binning);		
	TH1F* H1_etary_pid_acc = Histo::mk_etary(inputtree, etary, "H1_etary_pid_acc", SHMS_Accp, binning);		
	TH1F* H1_etarz_pid_acc = Histo::mk_etarz(inputtree, etarz, "H1_etarz_pid_acc", SHMS_Accp, binning);		
	
	TH1F* H1_etarx_pid_acc_1 = Histo::mk_etarx(inputtree, etarx, "H1_etarx_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_etary_pid_acc_1 = Histo::mk_etary(inputtree, etary, "H1_etary_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_etarz_pid_acc_1 = Histo::mk_etarz(inputtree, etarz, "H1_etarz_pid_acc_1", SHMS_Accp, binning);

	TH1F* H1_etarx_pid_acc_2 = Histo::mk_etarx(inputtree, etarx, "H1_etarx_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_etary_pid_acc_2 = Histo::mk_etary(inputtree, etary, "H1_etary_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_etarz_pid_acc_2 = Histo::mk_etarz(inputtree, etarz, "H1_etarz_pid_acc_2", SHMS_Accp, binning);

	TH1F* H1_etarx_pid_acc_kin = Histo::mk_etarx(inputtree, etarx, "H1_etarx_pid_acc_kin", SHMS_Accp, binning);		
	TH1F* H1_etary_pid_acc_kin = Histo::mk_etary(inputtree, etary, "H1_etary_pid_acc_kin", SHMS_Accp, binning);		
	TH1F* H1_etarz_pid_acc_kin = Histo::mk_etarz(inputtree, etarz, "H1_etarz_pid_acc_kin", SHMS_Accp, binning);		
	
	//----- HMS Collimator -----
	Double_t eXColl;		TH1F* H1_eXColl =					Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl", SHMS_Accp, binning);		
	TH1F* H1_eXColl_pid =		Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl_pid", SHMS_Accp, binning);		
	TH1F* H1_eXColl_pid_acc = Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl_pid_acc", SHMS_Accp, binning);		
	TH1F* H1_eXColl_pid_acc_1 = Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_eXColl_pid_acc_2 = Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_eXColl_pid_acc_kin = Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl_pid_acc_kin", SHMS_Accp, binning);
	
	Double_t eYColl;		TH1F* H1_eYColl = Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl", SHMS_Accp, binning);
	TH1F* H1_eYColl_pid = Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl_pid", SHMS_Accp, binning);
	TH1F* H1_eYColl_pid_acc = Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl_pid_acc", SHMS_Accp, binning);
	TH1F* H1_eYColl_pid_acc_1 = Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_eYColl_pid_acc_2 = Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_eYColl_pid_acc_kin = Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl_pid_acc_kin", SHMS_Accp, binning);
	
	//----- Calculated -----
	Double_t ztar_diff;		TH1F* H1_ztar_diff =				Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff", SHMS_Accp, binning);		
	TH1F* H1_ztar_diff_pid =	Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff_pid", SHMS_Accp, binning);		
	TH1F* H1_ztar_diff_pid_acc = Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff_pid_acc", SHMS_Accp, binning);		
	TH1F* H1_ztar_diff_pid_acc_1 = Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff_pid_acc_1", SHMS_Accp, binning);
	TH1F* H1_ztar_diff_pid_acc_2 = Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff_pid_acc_2", SHMS_Accp, binning);
	TH1F* H1_ztar_diff_pid_acc_kin = Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff_pid_acc_kin", SHMS_Accp, binning);



	//--------------------------------------------------------
	//--------------------- 2D HISTOGRAMS --------------------
	//--------------------------------------------------------

	//----- Kinematic -----
	TH2F* H2_cthrq_Pm_pid_acc = Histo::mk_cthrq_Pm(inputtree, "H2_cthrq_pmiss_pid_acc", HList2, binning);
	TH2F* H2_cthrq_Pm_pid_acc_1 = Histo::mk_cthrq_Pm(inputtree, "H2_cthrq_pmiss_pid_acc_1", HList2, binning);
	TH2F* H2_cthrq_Pm_pid_acc_2 = Histo::mk_cthrq_Pm(inputtree, "H2_cthrq_pmiss_pid_acc_2", HList2, binning);
	TH2F* H2_cthrq_Pm_pid_acc_kin = Histo::mk_cthrq_Pm(inputtree, "H2_cthrq_pmiss_pid_acc_kin", HList2, binning);

	TH2F* H2_thrq_Pm_pid_acc = Histo::mk_thrq_Pm(inputtree, "H2_thrq_pmiss_pid_acc", HList2, binning);
	TH2F* H2_thrq_Pm_pid_acc_1 = Histo::mk_thrq_Pm(inputtree, "H2_thrq_pmiss_pid_acc_1", HList2, binning);
	TH2F* H2_thrq_Pm_pid_acc_2 = Histo::mk_thrq_Pm(inputtree, "H2_thrq_pmiss_pid_acc_2", HList2, binning);
	TH2F* H2_thrq_Pm_pid_acc_kin = Histo::mk_thrq_Pm(inputtree, "H2_thrq_pmiss_pid_acc_kin", HList2, binning);

	TH2F* H2_thrq_Em_pid_acc = Histo::mk_thrq_Em(inputtree, "H2_thrq_Em_pid_acc", HList2, binning);
	TH2F* H2_thrq_Em_pid_acc_1 = Histo::mk_thrq_Em(inputtree, "H2_thrq_Em_pid_acc_1", HList2, binning);
	TH2F* H2_thrq_Em_pid_acc_2 = Histo::mk_thrq_Em(inputtree, "H2_thrq_Em_pid_acc_2", HList2, binning);
	TH2F* H2_thrq_Em_pid_acc_kin = Histo::mk_thrq_Em(inputtree, "H2_thrq_Em_pid_acc_kin", HList2, binning);

	TH2F* H2_thrq_Q2_pid_acc = Histo::mk_thrq_Q2(inputtree, "H2_thrq_Q2_pid_acc", HList2, binning);
	TH2F* H2_thrq_Q2_pid_acc_1 = Histo::mk_thrq_Q2(inputtree, "H2_thrq_Q2_pid_acc_1", HList2, binning);
	TH2F* H2_thrq_Q2_pid_acc_2 = Histo::mk_thrq_Q2(inputtree, "H2_thrq_Q2_pid_acc_2", HList2, binning);
	TH2F* H2_thrq_Q2_pid_acc_kin = Histo::mk_thrq_Q2(inputtree, "H2_thrq_Q2_pid_acc_kin", HList2, binning);

	TH2F* H2_thrq_xbj_pid_acc = Histo::mk_thrq_xbj(inputtree, "H2_thrq_xbj_pid_acc", HList2, binning);
	TH2F* H2_thrq_xbj_pid_acc_1 = Histo::mk_thrq_xbj(inputtree, "H2_thrq_xbj_pid_acc_1", HList2, binning);
	TH2F* H2_thrq_xbj_pid_acc_2 = Histo::mk_thrq_xbj(inputtree, "H2_thrq_xbj_pid_acc_2", HList2, binning);
	TH2F* H2_thrq_xbj_pid_acc_kin = Histo::mk_thrq_xbj(inputtree, "H2_thrq_xbj_pid_acc_kin", HList2, binning);

	TH2F* H2_xbj_Q2_pid_acc = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_pid_acc", HList2, binning);
	TH2F* H2_xbj_Q2_pid_acc_1 = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_pid_acc_1", HList2, binning);
	TH2F* H2_xbj_Q2_pid_acc_2 = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_pid_acc_2", HList2, binning);
	TH2F* H2_xbj_Q2_pid_acc_kin = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_pid_acc_kin", HList2, binning);

	TH2F* H2_xbj_Em_pid_acc = Histo::mk_xbj_Em(inputtree, "H2_xbj_Em_pid_acc", HList2, binning);
	TH2F* H2_xbj_Em_pid_acc_1 = Histo::mk_xbj_Em(inputtree, "H2_xbj_Em_pid_acc_1", HList2, binning);
	TH2F* H2_xbj_Em_pid_acc_2 = Histo::mk_xbj_Em(inputtree, "H2_xbj_Em_pid_acc_2", HList2, binning);
	TH2F* H2_xbj_Em_pid_acc_kin = Histo::mk_xbj_Em(inputtree, "H2_xbj_Em_pid_acc_kin", HList2, binning);

	TH2F* H2_xbj_Pm_pid_acc = Histo::mk_xbj_Pm(inputtree, "H2_xbj_Pm_pid_acc", HList2, binning);
	TH2F* H2_xbj_Pm_pid_acc_1 = Histo::mk_xbj_Pm(inputtree, "H2_xbj_Pm_pid_acc_1", HList2, binning);
	TH2F* H2_xbj_Pm_pid_acc_2 = Histo::mk_xbj_Pm(inputtree, "H2_xbj_Pm_pid_acc_2", HList2, binning);
	TH2F* H2_xbj_Pm_pid_acc_kin = Histo::mk_xbj_Pm(inputtree, "H2_xbj_Pm_pid_acc_kin", HList2, binning);

	TH2F* H2_Em_Pm_pid_acc = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_pid_acc", HList2, binning);
	TH2F* H2_Em_Pm_pid_acc_1 = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_pid_acc_1", HList2, binning);
	TH2F* H2_Em_Pm_pid_acc_2 = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_pid_acc_2", HList2, binning);
	TH2F* H2_Em_Pm_pid_acc_kin = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_pid_acc_kin", HList2, binning);

	TH2F* H2_Em_Pm_c_pid_acc = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_c_pid_acc", HList2, binning);
	TH2F* H2_Em_Pm_c_pid_acc_1 = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_c_pid_acc_1", HList2, binning);
	TH2F* H2_Em_Pm_c_pid_acc_2 = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_c_pid_acc_2", HList2, binning);
	TH2F* H2_Em_Pm_c_pid_acc_kin = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_c_pid_acc_kin", HList2, binning);

	TH2F* H2_Em_Q2_pid_acc = Histo::mk_Em_Q2(inputtree, "H2_Em_Q2_pid_acc", HList2, binning);
	TH2F* H2_Em_Q2_pid_acc_1 = Histo::mk_Em_Q2(inputtree, "H2_Em_Q2_pid_acc_1", HList2, binning);
	TH2F* H2_Em_Q2_pid_acc_2 = Histo::mk_Em_Q2(inputtree, "H2_Em_Q2_pid_acc_2", HList2, binning);
	TH2F* H2_Em_Q2_pid_acc_kin = Histo::mk_Em_Q2(inputtree, "H2_Em_Q2_pid_acc_kin", HList2, binning);

	TH2F* H2_Q2_Pm_pid_acc = Histo::mk_Q2_Pm(inputtree, "H2_Q2_Pm_pid_acc", HList2, binning);
	TH2F* H2_Q2_Pm_pid_acc_1 = Histo::mk_Q2_Pm(inputtree, "H2_Q2_Pm_pid_acc_1", HList2, binning);
	TH2F* H2_Q2_Pm_pid_acc_2 = Histo::mk_Q2_Pm(inputtree, "H2_Q2_Pm_pid_acc_2", HList2, binning);
	TH2F* H2_Q2_Pm_pid_acc_kin = Histo::mk_Q2_Pm(inputtree, "H2_Q2_Pm_pid_acc_kin", HList2, binning);

	TH2F* H2_W_Pm_pid_acc = Histo::mk_W_Pm(inputtree, "H2_W_Pm_pid_acc", HList2, binning);
	TH2F* H2_W_Pm_pid_acc_1 = Histo::mk_W_Pm(inputtree, "H2_W_Pm_pid_acc_1", HList2, binning);
	TH2F* H2_W_Pm_pid_acc_2 = Histo::mk_W_Pm(inputtree, "H2_W_Pm_pid_acc_2", HList2, binning);
	TH2F* H2_W_Pm_pid_acc_kin = Histo::mk_W_Pm(inputtree, "H2_W_Pm_pid_acc_kin", HList2, binning);
	
	TH2F* H2_W_thrq_pid_acc = Histo::mk_W_thrq(inputtree, "H2_W_thrq_pid_acc", HList2, binning);
	TH2F* H2_W_thrq_pid_acc_1 = Histo::mk_W_thrq(inputtree, "H2_W_thrq_pid_acc_1", HList2, binning);
	TH2F* H2_W_thrq_pid_acc_2 = Histo::mk_W_thrq(inputtree, "H2_W_thrq_pid_acc_2", HList2, binning);
	TH2F* H2_W_thrq_pid_acc_kin = Histo::mk_W_thrq(inputtree, "H2_W_thrq_pid_acc_kin", HList2, binning);

	TH2F* H2_W_Em_pid_acc = Histo::mk_W_Em(inputtree, "H2_W_Em_pid_acc", HList2, binning);
	TH2F* H2_W_Em_pid_acc_1 = Histo::mk_W_Em(inputtree, "H2_W_Em_pid_acc_1", HList2, binning);
	TH2F* H2_W_Em_pid_acc_2 = Histo::mk_W_Em(inputtree, "H2_W_Em_pid_acc_2", HList2, binning);
	TH2F* H2_W_Em_pid_acc_kin = Histo::mk_W_Em(inputtree, "H2_W_Em_pid_acc_kin", HList2, binning);

	TH2F* H2_W_Q2_pid_acc = Histo::mk_W_Q2(inputtree, "H2_W_Q2_pid_acc", HList2, binning);
	TH2F* H2_W_Q2_pid_acc_1 = Histo::mk_W_Q2(inputtree, "H2_W_Q2_pid_acc_1", HList2, binning);
	TH2F* H2_W_Q2_pid_acc_2 = Histo::mk_W_Q2(inputtree, "H2_W_Q2_pid_acc_2", HList2, binning);
	TH2F* H2_W_Q2_pid_acc_kin = Histo::mk_W_Q2(inputtree, "H2_W_Q2_pid_acc_kin", HList2, binning);

	TH2F* H2_W_xbj_pid_acc = Histo::mk_W_xbj(inputtree, "H2_W_xbj_pid_acc", HList2, binning);
	TH2F* H2_W_xbj_pid_acc_1 = Histo::mk_W_xbj(inputtree, "H2_W_xbj_pid_acc_1", HList2, binning);
	TH2F* H2_W_xbj_pid_acc_2 = Histo::mk_W_xbj(inputtree, "H2_W_xbj_pid_acc_2", HList2, binning);
	TH2F* H2_W_xbj_pid_acc_kin = Histo::mk_W_xbj(inputtree, "H2_W_xbj_pid_acc_kin", HList2, binning);

	TH2F* H2_alpha_Pm_per_pid_acc = Histo::mk_alpha_Pm_per(inputtree, "H2_alpha_Pm_per_pid_acc", HList2, binning);
	TH2F* H2_alpha_Pm_per_pid_acc_1 = Histo::mk_alpha_Pm_per(inputtree, "H2_alpha_Pm_per_pid_acc_1", HList2, binning);
	TH2F* H2_alpha_Pm_per_pid_acc_2 = Histo::mk_alpha_Pm_per(inputtree, "H2_alpha_Pm_per_pid_acc_2", HList2, binning);
	TH2F* H2_alpha_Pm_per_pid_acc_kin = Histo::mk_alpha_Pm_per(inputtree, "H2_alpha_Pm_per_pid_acc_kin", HList2, binning);

	TH2F* H2_exptar_Pm_pid_acc = Histo::mk_exptar_Pm(inputtree, "H2_exptar_Pm_pid_acc", HList2, binning);
	TH2F* H2_exptar_Pm_pid_acc_1 = Histo::mk_exptar_Pm(inputtree, "H2_exptar_Pm_pid_acc_1", HList2, binning);
	TH2F* H2_exptar_Pm_pid_acc_2 = Histo::mk_exptar_Pm(inputtree, "H2_exptar_Pm_pid_acc_2", HList2, binning);
	TH2F* H2_exptar_Pm_pid_acc_kin = Histo::mk_exptar_Pm(inputtree, "H2_exptar_Pm_pid_acc_kin", HList2, binning);

	TH2F* H2_eyptar_Pm_pid_acc = Histo::mk_eyptar_Pm(inputtree, "H2_eyptar_Pm_pid_acc", HList2, binning);
	TH2F* H2_eyptar_Pm_pid_acc_1 = Histo::mk_eyptar_Pm(inputtree, "H2_eyptar_Pm_pid_acc_1", HList2, binning);
	TH2F* H2_eyptar_Pm_pid_acc_2 = Histo::mk_eyptar_Pm(inputtree, "H2_eyptar_Pm_pid_acc_2", HList2, binning);
	TH2F* H2_eyptar_Pm_pid_acc_kin = Histo::mk_eyptar_Pm(inputtree, "H2_eyptar_Pm_pid_acc_kin", HList2, binning);

	TH2F* H2_hxptar_Pm_pid_acc = Histo::mk_hxptar_Pm(inputtree, "H2_hxptar_Pm_pid_acc", HList2, binning);
	TH2F* H2_hxptar_Pm_pid_acc_1 = Histo::mk_hxptar_Pm(inputtree, "H2_hxptar_Pm_pid_acc_1", HList2, binning);
	TH2F* H2_hxptar_Pm_pid_acc_2 = Histo::mk_hxptar_Pm(inputtree, "H2_hxptar_Pm_pid_acc_2", HList2, binning);
	TH2F* H2_hxptar_Pm_pid_acc_kin = Histo::mk_hxptar_Pm(inputtree, "H2_hxptar_Pm_pid_acc_kin", HList2, binning);

	TH2F* H2_hyptar_Pm_pid_acc = Histo::mk_hyptar_Pm(inputtree, "H2_hyptar_Pm_pid_acc", HList2, binning);
	TH2F* H2_hyptar_Pm_pid_acc_1 = Histo::mk_hyptar_Pm(inputtree, "H2_hyptar_Pm_pid_acc_1", HList2, binning);
	TH2F* H2_hyptar_Pm_pid_acc_2 = Histo::mk_hyptar_Pm(inputtree, "H2_hyptar_Pm_pid_acc_2", HList2, binning);
	TH2F* H2_hyptar_Pm_pid_acc_kin = Histo::mk_hyptar_Pm(inputtree, "H2_hyptar_Pm_pid_acc_kin", HList2, binning);



	//----- Acceptance -----
	TH2F* H2_hXColl_hYColl = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl", HList2, binning);
	TH2F* H2_eXColl_eYColl = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl", HList2, binning);
	TH2F* H2_hxfp_hyfp = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp", HList2, binning);
	TH2F* H2_exfp_eyfp = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp", HList2, binning);
	TH2F* H2_exptar_eyptar = Histo::mk_exptar_eyptar(inputtree, "H2_exptar_eyptar", HList2, binning);
	TH2F* H2_hxptar_hyptar = Histo::mk_hxptar_hyptar(inputtree, "H2_hxptar_hyptar", HList2, binning);
	TH2F* H2_hxptar_exptar = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar", HList2, binning);
	TH2F* H2_hyptar_eyptar = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar", HList2, binning);
	TH2F* H2_hdelta_edelta = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta", HList2, binning);
	TH2F* H2_pCalEtotTrkNorm_edelta = Histo::mk_pCalEtotTrkNorm_edelta(inputtree, "H2_pCalEtotTrkNorm_edelta", HList2, binning);

	TH2F* H2_hXColl_hYColl_pid = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_pid", HList2, binning);
	TH2F* H2_eXColl_eYColl_pid = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_pid", HList2, binning);
	TH2F* H2_hxfp_hyfp_pid = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_pid", HList2, binning);
	TH2F* H2_exfp_eyfp_pid = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_pid", HList2, binning);
	TH2F* H2_exptar_eyptar_pid = Histo::mk_exptar_eyptar(inputtree, "H2_exptar_eyptar_pid", HList2, binning);
	TH2F* H2_hxptar_hyptar_pid = Histo::mk_hxptar_hyptar(inputtree, "H2_hxptar_hyptar_pid", HList2, binning);
	TH2F* H2_hxptar_exptar_pid = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_pid", HList2, binning);
	TH2F* H2_hyptar_eyptar_pid = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_pid", HList2, binning);
	TH2F* H2_hdelta_edelta_pid = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_pid", HList2, binning);
	TH2F* H2_pCalEtotTrkNorm_edelta_pid = Histo::mk_pCalEtotTrkNorm_edelta(inputtree, "H2_pCalEtotTrkNorm_edelta_pid", HList2, binning);

	TH2F* H2_hXColl_hYColl_pid_acc = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_pid_acc", HList2, binning);
	TH2F* H2_eXColl_eYColl_pid_acc = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_pid_acc", HList2, binning);
	TH2F* H2_hxfp_hyfp_pid_acc = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_pid_acc", HList2, binning);
	TH2F* H2_exfp_eyfp_pid_acc = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_pid_acc", HList2, binning);
	TH2F* H2_exptar_eyptar_pid_acc = Histo::mk_exptar_eyptar(inputtree, "H2_exptar_eyptar_pid_acc", HList2, binning);
	TH2F* H2_hxptar_hyptar_pid_acc = Histo::mk_hxptar_hyptar(inputtree, "H2_hxptar_hyptar_pid_acc", HList2, binning);
	TH2F* H2_hxptar_exptar_pid_acc = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_pid_acc", HList2, binning);
	TH2F* H2_hyptar_eyptar_pid_acc = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_pid_acc", HList2, binning);
	TH2F* H2_hdelta_edelta_pid_acc = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_pid_acc", HList2, binning);
	TH2F* H2_pCalEtotTrkNorm_edelta_pid_acc = Histo::mk_pCalEtotTrkNorm_edelta(inputtree, "H2_pCalEtotTrkNorm_edelta_pid_acc", HList2, binning);

	TH2F* H2_hXColl_hYColl_pid_acc_1 = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_pid_acc_1", HList2, binning);
	TH2F* H2_eXColl_eYColl_pid_acc_1 = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_pid_acc_1", HList2, binning);
	TH2F* H2_hxfp_hyfp_pid_acc_1 = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_pid_acc_1", HList2, binning);
	TH2F* H2_exfp_eyfp_pid_acc_1 = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_pid_acc_1", HList2, binning);
	TH2F* H2_exptar_eyptar_pid_acc_1 = Histo::mk_exptar_eyptar(inputtree, "H2_exptar_eyptar_pid_acc_1", HList2, binning);
	TH2F* H2_hxptar_hyptar_pid_acc_1 = Histo::mk_hxptar_hyptar(inputtree, "H2_hxptar_hyptar_pid_acc_1", HList2, binning);
	TH2F* H2_hxptar_exptar_pid_acc_1 = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_pid_acc_1", HList2, binning);
	TH2F* H2_hyptar_eyptar_pid_acc_1 = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_pid_acc_1", HList2, binning);
	TH2F* H2_hdelta_edelta_pid_acc_1 = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_pid_acc_1", HList2, binning);
	TH2F* H2_pCalEtotTrkNorm_edelta_pid_acc_1 = Histo::mk_pCalEtotTrkNorm_edelta(inputtree, "H2_pCalEtotTrkNorm_edelta_pid_acc_1", HList2, binning);

	TH2F* H2_hXColl_hYColl_pid_acc_2 = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_pid_acc_2", HList2, binning);
	TH2F* H2_eXColl_eYColl_pid_acc_2 = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_pid_acc_2", HList2, binning);
	TH2F* H2_hxfp_hyfp_pid_acc_2 = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_pid_acc_2", HList2, binning);
	TH2F* H2_exfp_eyfp_pid_acc_2 = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_pid_acc_2", HList2, binning);
	TH2F* H2_exptar_eyptar_pid_acc_2 = Histo::mk_exptar_eyptar(inputtree, "H2_exptar_eyptar_pid_acc_2", HList2, binning);
	TH2F* H2_hxptar_hyptar_pid_acc_2 = Histo::mk_hxptar_hyptar(inputtree, "H2_hxptar_hyptar_pid_acc_2", HList2, binning);
	TH2F* H2_hxptar_exptar_pid_acc_2 = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_pid_acc_2", HList2, binning);
	TH2F* H2_hyptar_eyptar_pid_acc_2 = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_pid_acc_2", HList2, binning);
	TH2F* H2_hdelta_edelta_pid_acc_2 = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_pid_acc_2", HList2, binning);
	TH2F* H2_pCalEtotTrkNorm_edelta_pid_acc_2 = Histo::mk_pCalEtotTrkNorm_edelta(inputtree, "H2_pCalEtotTrkNorm_edelta_pid_acc_2", HList2, binning);


	TH2F* H2_hXColl_hYColl_pid_acc_kin = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_pid_acc_kin", HList2, binning);
	TH2F* H2_eXColl_eYColl_pid_acc_kin = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_pid_acc_kin", HList2, binning);
	TH2F* H2_hxfp_hyfp_pid_acc_kin = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_pid_acc_kin", HList2, binning);
	TH2F* H2_exfp_eyfp_pid_acc_kin = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_pid_acc_kin", HList2, binning);
	TH2F* H2_exptar_eyptar_pid_acc_kin = Histo::mk_exptar_eyptar(inputtree, "H2_exptar_eyptar_pid_acc_kin", HList2, binning);
	TH2F* H2_hxptar_hyptar_pid_acc_kin = Histo::mk_hxptar_hyptar(inputtree, "H2_hxptar_hyptar_pid_acc_kin", HList2, binning);
	TH2F* H2_hxptar_exptar_pid_acc_kin = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_pid_acc_kin", HList2, binning);
	TH2F* H2_hyptar_eyptar_pid_acc_kin = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_pid_acc_kin", HList2, binning);
	TH2F* H2_hdelta_edelta_pid_acc_kin = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_pid_acc_kin", HList2, binning);
	TH2F* H2_pCalEtotTrkNorm_edelta_pid_acc_kin = Histo::mk_pCalEtotTrkNorm_edelta(inputtree, "H2_pCalEtotTrkNorm_edelta_pid_acc_kin", HList2, binning);






	//--------------------------------------------------------
	//Define HMS Collimator Shape
	//--------------------------------------------------------
	TCutG* contam_gCut = new TCutG("contamCut", 5);
	contam_gCut->SetVarX("X");
	contam_gCut->SetVarY("Y");

	contam_gCut->SetPoint(0, -0.02, 0.0);
	contam_gCut->SetPoint(1, -0.02, 0.04);
	contam_gCut->SetPoint(3, 0.06, 0.04);
	contam_gCut->SetPoint(2, 0.07, 0.03);
	contam_gCut->SetPoint(4, 0.04, 0.0);
	contam_gCut->SetPoint(5, -0.02, 0.0);


	Double_t hms_scale = 1;
	Double_t hms_hsize = 4.575;
	Double_t hms_vsize = 11.646;
	Double_t shms_scale = 1.0;
	Double_t shms_hsize = 8.5;// / 253.0;
	Double_t shms_vsize = 12.5;// / 253.0;

	//Scaling the HMS/SHMS Collimator Cuts
	hms_hsize = hms_scale * hms_hsize;
	hms_vsize = hms_scale * hms_vsize;

	shms_hsize = shms_scale * shms_hsize;
	shms_vsize = shms_scale * shms_vsize;

	//Define HMS Collimator Shape
	TCutG* hms_Coll_gCut = new TCutG("hmsCollCut", 8);
	hms_Coll_gCut->SetVarX("X");
	hms_Coll_gCut->SetVarY("Y");

	hms_Coll_gCut->SetPoint(0, hms_hsize, hms_vsize / 2.);
	hms_Coll_gCut->SetPoint(1, hms_hsize / 2., hms_vsize);
	hms_Coll_gCut->SetPoint(2, -hms_hsize / 2., hms_vsize);
	hms_Coll_gCut->SetPoint(3, -hms_hsize, hms_vsize / 2.);
	hms_Coll_gCut->SetPoint(4, -hms_hsize, -hms_vsize / 2.);
	hms_Coll_gCut->SetPoint(5, -hms_hsize / 2., -hms_vsize);
	hms_Coll_gCut->SetPoint(6, hms_hsize / 2., -hms_vsize);
	hms_Coll_gCut->SetPoint(7, hms_hsize, -hms_vsize / 2.);
	hms_Coll_gCut->SetPoint(8, hms_hsize, hms_vsize / 2.);

	//Define SHMS Collimator Shape
	TCutG* shms_Coll_gCut = new TCutG("shmsCollCut", 8);
	shms_Coll_gCut->SetVarX("X");
	shms_Coll_gCut->SetVarY("Y");

	shms_Coll_gCut->SetPoint(0, shms_hsize, shms_vsize / 2.);
	shms_Coll_gCut->SetPoint(1, shms_hsize / 2., shms_vsize);
	shms_Coll_gCut->SetPoint(2, -shms_hsize / 2., shms_vsize);
	shms_Coll_gCut->SetPoint(3, -shms_hsize, shms_vsize / 2.);
	shms_Coll_gCut->SetPoint(4, -shms_hsize, -shms_vsize / 2.);
	shms_Coll_gCut->SetPoint(5, -shms_hsize / 2., -shms_vsize);
	shms_Coll_gCut->SetPoint(6, shms_hsize / 2., -shms_vsize);
	shms_Coll_gCut->SetPoint(7, shms_hsize, -shms_vsize / 2.);
	shms_Coll_gCut->SetPoint(8, shms_hsize, shms_vsize / 2.);




	bool HMS_TRK = true;
	bool SHMS_TRK = true;
	bool ctime_L = false;
	bool ctime_R = false;
	bool ctime = false;
	bool pid_hms = false;
	bool pid_shms = false;
	bool accp_hms = false;
	bool accp_shms = false;
	bool contam_cut = false;
	bool kin = false;
	bool bevtyp = false;
	bool accp_z;


	Double_t etarx_corr;

	Double_t PSF;
	Double_t Ehtrk;
	Double_t Eetrk;
	Double_t Emtrk;
	Double_t Elt;
	Double_t charge = 0;
	Double_t weight;
	Double_t weight2;
	Double_t nentries;

	//--------------------------------------------------------
	//Get net charge
	//--------------------------------------------------------
	if (user_pass == "pass1")
	{
		if (user_runtype == "SRC")
		{
			for (int i = 0; i < file_count; i++)
			{
				charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_SRC_report_%s_-1.txt", path1_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);//charge is later
			}
		}
		else if (user_runtype == "MF")
		{
			for (int i = 0; i < file_count; i++)
			{
				charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_MF_report_%s_-1.txt", path1_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);//charge is later
			}
		}
		else if (user_runtype == "Heep")
		{
			for (int i = 0; i < file_count; i++)
			{
				charge = charge + stod(split(FindString("HMS BCM1  Charge", Form("../../../../../../../../OFFLINE/PASS1/pass1_heep_replayed/cafe_prod_%s_-1.report", Heep_File[0][i]))[0], ':')[1]);//charge is later
			}
		}
	}
	else if (user_pass == "pass2")
	{
		if (user_runtype == "SRC")
		{
			for (int i = 0; i < file_count; i++)
			{
				charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path2_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);//charge is later
			}
		}
		else if (user_runtype == "MF")
		{
			for (int i = 0; i < file_count; i++)
			{
				charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path2_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);//charge is later
			}
		}
		else if (user_runtype == "Heep")
		{
			for (int i = 0; i < file_count; i++)
			{
				charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("../../../../../../../../OFFLINE/PASS2/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);//charge is later
			}
		}
	}
	else if (user_pass == "pass3")
	{
		if (user_runtype == "SRC")
		{
			for (int i = 0; i < file_count; i++)
			{
				charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path3_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);//charge is later
			}
		}
		else if (user_runtype == "MF")
		{
			for (int i = 0; i < file_count; i++)
			{
				charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path3_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);//charge is later
			}
		}
		else if (user_runtype == "Heep")
		{
			for (int i = 0; i < file_count; i++)
			{
				charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("../../../../../../../../OFFLINE/PASS3/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);//charge is later
			}
		}
	}
	else if (user_pass == "pass4")
	{
		if (user_s_or_c == "coin")
		{
			if (user_runtype == "SRC")
			{
				for (int i = 0; i < file_count; i++)
				{
					charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);//charge is later
					cout << "charge per run: " << (charge - temp1) << endl;
				}
			}
			else if (user_runtype == "MF")
			{
				for (int i = 0; i < file_count; i++)
				{
					charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);//charge is later
					cout << "charge per run: " << (charge - temp1) << endl;
				}
			}
			else if (user_runtype == "Heep")
			{
				for (int i = 0; i < file_count; i++)
				{
					charge = charge + stod(split(FindString("BCM1_Charge [mC]", Form("../../../../../../../../OFFLINE/PASS4/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);//charge is later
				}
			}
		}
		else if (user_s_or_c == "sing")
		{
			if (user_runtype == "SRC")
			{
				for (int i = 0; i < file_count; i++)
				{
					charge = charge + stod(split(FindString("SHMS BCM1  Charge", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.report", path4_c1.c_str(), SRC_File[target_id][i]))[0], ':')[1]);//charge is later
					cout << "charge per run: " << (charge - temp1) << endl;
					temp1 = charge;
				}
			}
			else if (user_runtype == "MF")
			{
				for (int i = 0; i < file_count; i++)
				{
					charge = charge + stod(split(FindString("SHMS BCM1  Charge", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.report", path4_c1.c_str(), MF_File[target_id][i]))[0], ':')[1]);//charge is later
					cout << "charge per run: " << (charge - temp1) << endl;
					temp1 = charge; charge = charge + stod(split(FindString("SHMS BCM1  Charge", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.report", path4_c1.c_str(), MF_File[target_id][i]))[0], ':')[1]);//charge is later
				}
			}
			else if (user_runtype == "Heep")
			{
				for (int i = 0; i < file_count; i++)
				{
					charge = charge + stod(split(FindString("SHMS BCM1  Charge", Form("../../../../../../../../OFFLINE/PASS4/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);//charge is later
				}
			}
		}
	}
		
	

	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	//-------------------------------------------------------------------------- Loop ---------------------------------------------------------------------------
	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	int numLoops = 0;							//Initialize event count to 0
	int max_loops = 70000000;					//Define maximum number of events to process

	TVector3 v3_ebeam;
	TLorentzVector v4_ebeam;
	TVector3 v3_eprime;
	TLorentzVector v4_eprime;

	TVector3 v3_p;
	TLorentzVector v4_p;

	TVector3 v3_q;
	TLorentzVector v4_q;

	TVector3 v3_pr;
	TLorentzVector v4_pr;

	Double_t p4x = 0; inputtree->SetBranchAddress("P.kin.primary.p4x", &p4x);
	Double_t p4y = 0; inputtree->SetBranchAddress("P.kin.primary.p4y", &p4x);
	Double_t p4z = 0; inputtree->SetBranchAddress("P.kin.primary.p4z", &p4x);
	Double_t p4e = 0; inputtree->SetBranchAddress("P.kin.primary.p4e", &p4x);
	
	Double_t h4x = 0; inputtree->SetBranchAddress("H.kin.secondary.p4x", &h4x);
	Double_t h4y = 0; inputtree->SetBranchAddress("H.kin.secondary.p4y", &h4x);
	Double_t h4z = 0; inputtree->SetBranchAddress("H.kin.secondary.p4z", &h4x);
	Double_t h4e = 0; inputtree->SetBranchAddress("H.kin.secondary.p4e", &h4x);

	Double_t Er = 0;
	Double_t alpha_r = 0;
	
	//--------------------------------------------------------
	//Readdress all the branches
	//--------------------------------------------------------
	for (int i = 0; i < file_count; i++)
	{
		if (i > 0)
		{			
			if (user_pass == "pass1")
			{
				if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/skimmed_pass1/cafe_prod_%s_SRC_%s_-1_skimmed.root", path1_c.c_str(), Target[target_id], SRC_File[target_id][i]), "READ"); }
				else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/skimmed_pass1/cafe_prod_%s_MF_%s_-1_skimmed.root", path1_c.c_str(), Target[target_id], MF_File[target_id][i]), "READ"); }
				else if (user_runtype == "Heep") { inROOT = new TFile(Form("%s/skimmed_pass1/cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", path1_c.c_str(), Heep_File[0][i]), "READ"); }
			}
			else if (user_pass == "pass2")
			{
				if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_SRC_%s_-1_skimmed.root", path2_c.c_str(), Target[target_id], SRC_File[target_id][i]), "READ"); }
				else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_MF_%s_-1_skimmed.root", path2_c.c_str(), Target[target_id], MF_File[target_id][i]), "READ"); }
				else if (user_runtype == "Heep") { inROOT = new TFile(Form("../../../../../../../../OFFLINE/PASS2/ROOT/cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", Heep_File[0][i]), "READ"); }
			}
			else if (user_pass == "pass3")
			{
				if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_SRC_%s_-1_skimmed.root", path3_c.c_str(), Target[target_id], SRC_File[target_id][i]), "READ"); }
				else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_MF_%s_-1_skimmed.root", path3_c.c_str(), Target[target_id], MF_File[target_id][i]), "READ"); }
				else if (user_runtype == "Heep") { inROOT = new TFile(Form("../../../../../../../../OFFLINE/PASS3/ROOT/cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", Heep_File[0][i]), "READ"); }
			}
			else if (user_pass == "pass4")
			{
				if (user_s_or_c == "coin")
				{
					if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_SRC_%s_-1_skimmed.root", path4_c1.c_str(), Target[target_id], SRC_File[target_id][i]), "READ"); }
					else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/ROOT/cafe_prod_%s_MF_%s_-1_skimmed.root", path4_c1.c_str(), Target[target_id], MF_File[target_id][i]), "READ"); }
					else if (user_runtype == "Heep") { inROOT = new TFile(Form("../../../../../../../../OFFLINE/PASS4/ROOT/cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", Heep_File[0][i]), "READ"); }
				}
				else if (user_s_or_c == "sing")
				{

					if (user_runtype == "SRC") { inROOT = new TFile(Form("%s/ROOTfiles/cafe_replay_prod_%s_-1.root", path4_c1.c_str(), SRC_File[target_id][i]), "READ"); }
					else if (user_runtype == "MF") { inROOT = new TFile(Form("%s/ROOTfiles/cafe_replay_prod_%s_-1.root", path4_c1.c_str(), MF_File[target_id][i]), "READ"); }
					else if (user_runtype == "Heep") { inROOT = new TFile(Form("../../../../../../../../OFFLINE/PASS4/ROOT/cafe_prod_LH2_heep_coin_%s_-1_skimmed.root", Heep_File[0][i]), "READ"); }
				}
			}
			inputtree = (TTree*)inROOT->Get("T");

			inputtree->SetBranchStatus("*", kFALSE);

			inputtree->SetBranchStatus("shms_collimator_cut_flag", kTRUE);
			inputtree->SetBranchStatus("H.dc.ntrack", kTRUE);
			inputtree->SetBranchStatus("H.hod.goodscinhit", kTRUE);
			inputtree->SetBranchStatus("H.cer.npeSum", kTRUE);
			inputtree->SetBranchStatus("H.cal.etotnorm", kTRUE);
			inputtree->SetBranchStatus("H.cal.etottracknorm", kTRUE);
			inputtree->SetBranchStatus("H.hod.betanotrack", kTRUE);

			inputtree->SetBranchStatus("CTime.epCoinTime_ROC2_center", kTRUE);////or not center

			inputtree->SetBranchStatus("P.dc.ntrack", kTRUE);
			inputtree->SetBranchStatus("P.hod.goodscinhit", kTRUE);
			inputtree->SetBranchStatus("P.ngcer.npeSum", kTRUE);
			inputtree->SetBranchStatus("P.hgcer.npesum", kTRUE);
			inputtree->SetBranchStatus("P.cal.etotnorm", kTRUE);
			inputtree->SetBranchStatus("P.cal.etottracknorm", kTRUE);
			inputtree->SetBranchStatus("P.hod.betanotrack", kTRUE);

			inputtree->SetBranchStatus("P.kin.primary.scat_ang_rad", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.W", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.Q2", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.x_bj", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.nu", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.q3m", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.q_x", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.q_y", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.q_z", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.th_q", kTRUE);
			inputtree->SetBranchStatus("P.kin.primary.ph_q", kTRUE);

			if (user_target == "LH2") { inputtree->SetBranchStatus("H.kin.secondary.emiss", kTRUE); }		//Standard Missing Energy for H(e,e'p)
			else { inputtree->SetBranchStatus("H.kin.secondary.emiss_nuc", kTRUE); }

			inputtree->SetBranchStatus("H.kin.secondary.pmiss", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.Prec_x", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.Prec_y", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.Prec_z", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.pmiss_x", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.pmiss_y", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.pmiss_z", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.tx", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.tb", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.Mrecoil", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.th_xq", kTRUE);
			//inputtree->SetBranchStatus("H.kin.secondary.th_bq", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.th_bq", kTRUE);//recoil particle in-plane angle w.r.to q-vector
			inputtree->SetBranchStatus("H.kin.secondary.ph_xq", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.ph_bq", kTRUE);
			inputtree->SetBranchStatus("H.kin.secondary.xangle", kTRUE);

			inputtree->SetBranchStatus("H.dc.x_fp", kTRUE);
			inputtree->SetBranchStatus("H.dc.xp_fp", kTRUE);
			inputtree->SetBranchStatus("H.dc.y_fp", kTRUE);
			inputtree->SetBranchStatus("H.dc.yp_fp", kTRUE);
			inputtree->SetBranchStatus("H.gtr.y", kTRUE);
			inputtree->SetBranchStatus("H.gtr.ph", kTRUE);
			inputtree->SetBranchStatus("H.gtr.th", kTRUE);
			inputtree->SetBranchStatus("H.gtr.dp", kTRUE);
			inputtree->SetBranchStatus("H.react.x", kTRUE);
			inputtree->SetBranchStatus("H.react.y", kTRUE);
			inputtree->SetBranchStatus("H.react.z", kTRUE);
			inputtree->SetBranchStatus("H.extcor.xsieve", kTRUE);
			inputtree->SetBranchStatus("H.extcor.ysieve", kTRUE);

			inputtree->SetBranchStatus("P.dc.x_fp", kTRUE);
			inputtree->SetBranchStatus("P.dc.xp_fp", kTRUE);
			inputtree->SetBranchStatus("P.dc.y_fp", kTRUE);
			inputtree->SetBranchStatus("P.dc.yp_fp", kTRUE);
			inputtree->SetBranchStatus("P.gtr.y", kTRUE);
			inputtree->SetBranchStatus("P.gtr.ph", kTRUE);
			inputtree->SetBranchStatus("P.gtr.p", kTRUE);//verify
			inputtree->SetBranchStatus("H.gtr.p", kTRUE);//verify
			inputtree->SetBranchStatus("P.gtr.th", kTRUE);
			inputtree->SetBranchStatus("P.gtr.dp", kTRUE);
			inputtree->SetBranchStatus("P.react.x", kTRUE);
			inputtree->SetBranchStatus("P.react.y", kTRUE);
			inputtree->SetBranchStatus("P.react.z", kTRUE);
			inputtree->SetBranchStatus("P.extcor.xsieve", kTRUE);
			inputtree->SetBranchStatus("P.extcor.ysieve", kTRUE);

			inputtree->SetBranchStatus("g.evtyp", kTRUE);
			inputtree->SetBranchStatus("P.dc.TheRealGolden", kTRUE);


			//inputtree->SetBranchAddress("shms_collimator_cut_flag", &shms_coll_cut);

			
			inputtree->SetBranchAddress("H.dc.ntrack", &hdc_ntrk);
			inputtree->SetBranchAddress("H.hod.goodscinhit", &hScinGood);
			inputtree->SetBranchAddress("H.cer.npeSum", &hCerNpeSum);
			inputtree->SetBranchAddress("H.cal.etotnorm", &hCalEtotNorm);
			inputtree->SetBranchAddress("H.cal.etottracknorm", &hCalEtotTrkNorm);
			inputtree->SetBranchAddress("H.hod.betanotrack", &hHodBetaNtrk);
			
			inputtree->SetBranchAddress("CTime.epCoinTime_ROC2_center", &ep_ctime);////or not center

			inputtree->SetBranchAddress("P.dc.ntrack", &pdc_ntrk);
			inputtree->SetBranchAddress("P.hod.goodscinhit", &pScinGood);
			inputtree->SetBranchAddress("P.ngcer.npeSum", &pNGCerNpeSum);
			inputtree->SetBranchAddress("P.hgcer.npesum", &pHGCerNpeSum);
			inputtree->SetBranchAddress("P.cal.etotnorm", &pCalEtotNorm);
			inputtree->SetBranchAddress("P.cal.etottracknorm", &pCalEtotTrkNorm);
			inputtree->SetBranchAddress("P.hod.betanotrack", &pHodBetaNtrk);

			inputtree->SetBranchAddress("P.kin.primary.scat_ang_rad", &th_e);
			inputtree->SetBranchAddress("P.kin.primary.W", &W);
			inputtree->SetBranchAddress("P.kin.primary.Q2", &Q2);
			inputtree->SetBranchAddress("P.kin.primary.x_bj", &x_bj);
			inputtree->SetBranchAddress("P.kin.primary.nu", &nu);
			inputtree->SetBranchAddress("P.kin.primary.q3m", &q);
			inputtree->SetBranchAddress("P.kin.primary.q_x", &q_x);
			inputtree->SetBranchAddress("P.kin.primary.q_y", &q_y);
			inputtree->SetBranchAddress("P.kin.primary.q_z", &q_z);
			inputtree->SetBranchAddress("P.kin.primary.th_q", &th_q);
			inputtree->SetBranchAddress("P.kin.primary.ph_q", &ph_q);

			if (user_target == "LH2") { inputtree->SetBranchAddress("H.kin.secondary.emiss", &emiss); }		//Standard Missing Energy for H(e,e'p)
			else { inputtree->SetBranchAddress("H.kin.secondary.emiss_nuc", &emiss); }	

			inputtree->SetBranchAddress("H.kin.secondary.pmiss", &pmiss);
			inputtree->SetBranchAddress("H.kin.secondary.Prec_x", &prec_x);
			inputtree->SetBranchAddress("H.kin.secondary.Prec_y", &prec_y);
			inputtree->SetBranchAddress("H.kin.secondary.Prec_z", &prec_z);
			inputtree->SetBranchAddress("H.kin.secondary.pmiss_x", &pmiss_x);
			inputtree->SetBranchAddress("H.kin.secondary.pmiss_y", &pmiss_y);
			inputtree->SetBranchAddress("H.kin.secondary.pmiss_z", &pmiss_z);
			inputtree->SetBranchAddress("H.kin.secondary.tx", &Tx);
			inputtree->SetBranchAddress("H.kin.secondary.tb", &Tr);
			inputtree->SetBranchAddress("H.kin.secondary.Mrecoil", &mmiss);
			inputtree->SetBranchAddress("H.kin.secondary.th_xq", &th_pq);
			//inputtree->SetBranchAddress("H.kin.secondary.th_bq", &th_rq);
			inputtree->SetBranchAddress("H.kin.secondary.th_bq", &cth_rq);//recoil particle in-plane angle w.r.to q-vector
			inputtree->SetBranchAddress("H.kin.secondary.ph_xq", &ph_pq);
			inputtree->SetBranchAddress("H.kin.secondary.ph_bq", &ph_rq);
			inputtree->SetBranchAddress("H.kin.secondary.xangle", &xangle);

			inputtree->SetBranchAddress("H.dc.x_fp", &hxfp);
			inputtree->SetBranchAddress("H.dc.xp_fp", &hxpfp);
			inputtree->SetBranchAddress("H.dc.y_fp", &hyfp);
			inputtree->SetBranchAddress("H.dc.yp_fp", &hypfp);
			inputtree->SetBranchAddress("H.gtr.y", &hytar);
			inputtree->SetBranchAddress("H.gtr.ph", &hyptar);
			inputtree->SetBranchAddress("H.gtr.th", &hxptar);
			inputtree->SetBranchAddress("H.gtr.dp", &hdelta);
			inputtree->SetBranchAddress("H.react.x", &htarx);
			inputtree->SetBranchAddress("H.react.y", &htary);
			inputtree->SetBranchAddress("H.react.z", &htarz);
			inputtree->SetBranchAddress("H.extcor.xsieve", &hXColl);
			inputtree->SetBranchAddress("H.extcor.ysieve", &hYColl);

			inputtree->SetBranchAddress("P.dc.x_fp", &exfp);
			inputtree->SetBranchAddress("P.dc.xp_fp", &expfp);
			inputtree->SetBranchAddress("P.dc.y_fp", &eyfp);
			inputtree->SetBranchAddress("P.dc.yp_fp", &eypfp);
			inputtree->SetBranchAddress("P.gtr.y", &eytar);
			inputtree->SetBranchAddress("P.gtr.ph", &eyptar);
			inputtree->SetBranchAddress("P.gtr.p", &Ef);//verify
			inputtree->SetBranchAddress("H.gtr.p", &Pf);//verify
			inputtree->SetBranchAddress("P.gtr.th", &exptar);	
			inputtree->SetBranchAddress("P.gtr.dp", &edelta);
			inputtree->SetBranchAddress("P.react.x", &etarx);
			inputtree->SetBranchAddress("P.react.y", &etary);
			inputtree->SetBranchAddress("P.react.z", &etarz);
			inputtree->SetBranchAddress("P.extcor.xsieve", &eXColl);
			inputtree->SetBranchAddress("P.extcor.ysieve", &eYColl);

			inputtree->SetBranchAddress("g.evtyp", &evtyp);
			inputtree->SetBranchAddress("P.dc.TheRealGolden", &Real_Golden);

			/*
			H1_ctime_peak->Reset("ICESM");
			//H1_ctime_peak = new TH1F("H1_ctime_peak", "", 300, -50, 150);
			//Form("",)
			entries = inputtree->GetEntries();
			for (int i = 0; i < entries; i++) //50,000
			{
				inputtree->GetEntry(i);
				H1_ctime_peak->Fill(ep_ctime);
			}
			binmax_ctime = H1_ctime_peak->GetMaximumBin();
			ctime_offset = H1_ctime_peak->GetXaxis()->GetBinCenter(binmax_ctime);
			*/

			//--------------------------------------------------------
			//Update corrections
			//--------------------------------------------------------
			if (user_runtype == "MF")
			{
				if (user_target == "Be9") { ctime_offset = 0; }//86.583; }
				else if (user_target == "B10") { ctime_offset = 0; }//86.583; }
				else if (user_target == "B11") { ctime_offset = 0; }//94.75; }
				else if (user_target == "C12") { ctime_offset = 0; }//86.583; }
				else if (user_target == "Ca40") { ctime_offset = 0; }//86.583; }
			}

			if (user_runtype == "MF" && i == 2)
			{
				if (user_target == "B11") { ctime_offset = 0; }//86.583; }
			}

			if (user_runtype == "SRC")
			{
				if (user_target == "Ca48") { correction = Ca48_Corr[i]; }
			}
		}
		

		//--------------------------------------------------------
		//Update efficiencies
		//--------------------------------------------------------
		if (user_pass == "pass1")
		{
			if (user_runtype == "SRC")
			{
				Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_SRC_report_%s_-1.txt", path1_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_SRC_report_%s_-1.txt", path1_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_SRC_report_%s_-1.txt", path1_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_SRC_report_%s_-1.txt", Target[target_id], path1_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_SRC_report_%s_-1.txt", path1_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_SRC_report_%s_-1.txt", path1_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
			}
			else if (user_runtype == "MF")
			{
				Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_MF_report_%s_-1.txt", path1_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_MF_report_%s_-1.txt", path1_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_MF_report_%s_-1.txt", path1_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_MF_report_%s_-1.txt", path1_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_MF_report_%s_-1.txt", path1_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_report_pass1/cafe_prod_%s_MF_report_%s_-1.txt", path1_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
			}
			else if (user_runtype == "Heep")
			{
				Ehtrk = stod(split(FindString("SW_hms_hadron_trk_eff", Form("../../../../../../../../OFFLINE/PASS1/pass1_heep_replayed/cafe_prod_%s_-1.report", Heep_File[0][i]))[0], ':')[1]);
				Eetrk = stod(split(FindString("SW_shms_elec_trk_eff", Form("../../../../../../../../OFFLINE/PASS1/pass1_heep_replayed/cafe_prod_%s_-1.report", Heep_File[0][i]))[0], ':')[1]);
				Emtrk = stod(split(FindString("S_trk_eff", Form("../../../../../../../../OFFLINE/PASS1/pass1_heep_replayed/cafe_prod_%s_-1.report", Heep_File[0][i]))[0], ':')[1]);
				Elt = stod(split(FindString("SHMS Singles FID TRACK EFF"/*"T5_tLT"*/, Form("../../../../../../../../OFFLINE/PASS1/pass1_heep_replayed/cafe_prod_%s_-1.report", Heep_File[0][i]))[0], ':')[1]);
				PSF = stod(split(FindString("Ps2_factor", Form("../../../../../../../../OFFLINE/PASS1/pass1_heep_replayed/cafe_prod_%s_-1.report", Heep_File[0][i]))[0], ':')[1]);
				//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("../../../../../../../../OFFLINE/PASS1/pass1_heep_replayed/cafe_prod_%s_-1.report", Heep_File[target_id][i]))[0], ':')[1]);
			}
		}
		else if (user_pass == "pass2")
		{
			if (user_runtype == "SRC")
			{
				Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path2_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path2_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path2_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path2_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path2_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path2_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
			}
			else if (user_runtype == "MF")
			{
				Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path2_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path2_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path2_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path2_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path2_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path2_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
			}
			else if (user_runtype == "Heep")
			{
				Ehtrk = stod(split(FindString("hms_had_track_eff", Form("../../../../../../../../OFFLINE/PASS2/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
				Eetrk = stod(split(FindString("shms_elec_track_eff", Form("../../../../../../../../OFFLINE/PASS2/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
				Emtrk = stod(split(FindString("multi_track_eff", Form("../../../../../../../../OFFLINE/PASS2/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
				Elt = stod(split(FindString("T5_tLT", Form("../../../../../../../../OFFLINE/PASS2/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
				PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", path2_c.c_str(), Heep_File[target_id][i]))[0], ':')[1]);
				//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", path2_c.c_str(), Heep_File[target_id][i]))[0], ':')[1]);
			}
		}
		else if (user_pass == "pass3")
		{
			if (user_runtype == "SRC")
			{
				Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path3_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path3_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path3_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path3_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path3_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path3_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
			}
			else if (user_runtype == "MF")
			{
				Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path3_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path3_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path3_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path3_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path3_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path3_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
			}
			else if (user_runtype == "Heep")
			{
				Ehtrk = stod(split(FindString("hms_had_track_eff", Form("../../../../../../../../OFFLINE/PASS3/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
				Eetrk = stod(split(FindString("shms_elec_track_eff", Form("../../../../../../../../OFFLINE/PASS3/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
				//Emtrk = stod(split(FindString("shms_elec_track_eff", Form("../../../../../../../../OFFLINE/PASS3/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
				Elt = stod(split(FindString("T5_tLT", Form("../../../../../../../../OFFLINE/PASS3/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
				PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", path3_c.c_str(), Heep_File[target_id][i]))[0], ':')[1]);
				//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", path3_c.c_str(), Heep_File[target_id][i]))[0], ':')[1]);
			}
		}
		else if (user_pass == "pass4")
		{
			if (user_s_or_c == "coin")
			{
				if (user_runtype == "SRC")
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "MF")
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c1.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "Heep")
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("../../../../../../../../OFFLINE/PASS4/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("../../../../../../../../OFFLINE/PASS4/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("shms_elec_track_eff", Form("../../../../../../../../OFFLINE/PASS4/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("../../../../../../../../OFFLINE/PASS4/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", path4_c.c_str(), Heep_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", path4_c.c_str(), Heep_File[target_id][i]))[0], ':')[1]);
				}//target_thickness Count / (g/cm^2 * C *  
			}
			if (user_s_or_c == "sing")
			{
			/*Double_t SRC_File_Count[9] = { 1, 10, 8, 8, 13, 25, 21, 24, 30 }; //LD2, Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
			const char* SRC_File[9][31] = {
				{ "17134", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },																//LD2
				{ "17106", "17107", "17108", "17109", "17110", "17111", "17129", "17130", "17131", "17132", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },					//Be9
				{ "17112", "17113", "17114", "17115", "17125", "17126", "17127", "17128", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },								//B10
				{ "17116", "17117", "17118", "17119", "17120", "17121", "17122", "17123", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },								//B11
				{ "17076", "17078", "17079", "17080", "17081", "17082", "17083", "17084", "17085", "17086", "17087", "17088", "17089", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },	//C12
				{ "17008", "17009", "17011", "17012", "17013", "17014", "17015", "17016", "17017", "17018", "17020", "17021", "17022",
				"17023", "17024", "17025", "17027", "17028", "17029", "17030", "17031", "17032", "17033", "17034", "17035", "", "", "", "", "", "" },															//Ca40
				{ "17036", "17037", "17038", "17039", "17040", "17041", "17043", "17044", "17045", "17046", "17047", "17048", "17049",
				"17050", "17051", "17052", "17053", "17054", "17055", "17056", "17057", "", "", "", "", "", "", "", "", "", "" },																				//Ca48
				{ "17058", "17059", "17060", "17061", "17062", "17063", "17064", "17067", "17068", "17069", "17070", "17071", "17072",
				"17073", "17074", "17075", "17135", "17136", "17138", "17140", "17141", "17142", "17143", "17144", "", "", "", "", "", "", "" },																//Fe54
				{ "20800", "20801", "20802", "20803", "20804", "20805", "20806", "20807", "20808", "20809", "20810", "20811", "20812",
				"20815", "20816", "20817", "20819", "20821", "20822", "20823", "20824", "20825", "20826", "20827", "20828", "20829", "20830", "20831", "20832", "20833", "20834" }								//Au197
			};

			Double_t SRC_File_Count[9] = { 2, 13, 10, 11, 15, 27, 24, 26, 34 }; //LD2, Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
			const char* SRC_File[9][35] = {
				{ "17134", "16973", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },																//LD2
				{ "17106", "17107", "17108", "17109", "17110", "17111", "17129", "17130", "17131", "17132", "16983", "17099", "17100", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },					//Be9
				{ "17112", "17113", "17114", "17115", "17125", "17126", "17127", "17128", "16984", "17101", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },								//B10
				{ "17116", "17117", "17118", "17119", "17120", "17121", "17122", "17123", "16986", "16991", "17102", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },								//B11
				{ "17076", "17078", "17079", "17080", "17081", "17082", "17083", "17084", "17085", "17086", "17087", "17088", "17089", "16977", "17098", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },	//C12
				{ "17008", "17009", "17011", "17012", "17013", "17014", "17015", "17016", "17017", "17018", "17020", "17021", "17022",
				"17023", "17024", "17025", "17027", "17028", "17029", "17030", "17031", "17032", "17033", "17034", "17035", "16980", "17097", "", "", "", "" },															//Ca40
				{ "17036", "17037", "17038", "17039", "17040", "17041", "17043", "17044", "17045", "17046", "17047", "17048", "17049",
				"17050", "17051", "17052", "17053", "17054", "17055", "17056", "17057", "17093", "17094", "17096", "", "", "", "", "", "", "" },																				//Ca48
				{ "17058", "17059", "17060", "17061", "17062", "17063", "17064", "17067", "17068", "17069", "17070", "17071", "17072",
				"17073", "17074", "17075", "17135", "17136", "17138", "17140", "17141", "17142", "17143", "17144", "16981", "16982", "", "", "", "", "" },																//Fe54
				{ "20800", "20801", "20802", "20803", "20804", "20805", "20806", "20807", "20808", "20809", "20810", "20811", "20812",
				"20815", "20816", "20817", "20819", "20821", "20822", "20823", "20824", "20825", "20826", "20827", "20828", "20829", "20830", "20831", "20832", "20833", "20834", "20793", "20797", "20798", "20799" }								//Au197
			};
			*/
				if (user_runtype == "SRC" && user_target == "LD2" && i >= 1)
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
					//cout << Ehtrk << endl;
					//cout << Eetrk << endl;
					//cout << Emtrk << endl;
				}
				else if (user_runtype == "SRC" && user_target == "Be9" && i >= 10)
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "SRC" && user_target == "B10" && i >= 8)
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "SRC" && user_target == "B11" && i >= 8)
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "SRC" && user_target == "C12" && i >= 13)
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "SRC" && user_target == "Ca40" && i >= 25)
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "SRC" && user_target == "Ca48" && i >= 21)
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "SRC" && user_target == "Fe54" && i >= 24)
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "SRC" && user_target == "Au197" && i >= 30)
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "SRC")
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_SRC_report_%s_-1.txt", path4_c.c_str(), Target[target_id], SRC_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), SRC_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "MF")
				{
					Ehtrk = stod(split(FindString("hms_had_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("shms_elec_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					Emtrk = stod(split(FindString("multi_track_eff", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					Elt = stod(split(FindString("T5_tLT", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					PSF = stod(split(FindString("Ps2_factor", Form("%s/REPORT/cafe_prod_%s_MF_report_%s_-1.txt", path4_c.c_str(), Target[target_id], MF_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUTPUT/cafe_prod_%s_-1.txt", path4_c.c_str(), MF_File[target_id][i]))[0], ':')[1]);
				}
				else if (user_runtype == "Heep")
				{
					//Ehtrk = stod(split(FindString("hms_had_track_eff", Form("../../../../../../../../OFFLINE/PASS4/REPORT_OUTPUT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
					Eetrk = stod(split(FindString("SW_shms_elec_trk_eff", Form("../../../../../../../../OFFLINE/PASS4/REPORT_OUTPUT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
					//Emtrk = stod(split(FindString("shms_elec_track_eff", Form("../../../../../../../../OFFLINE/PASS4/REPORT_OUTPUT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
					//Elt = stod(split(FindString("T5_tLT", Form("../../../../../../../../OFFLINE/PASS4/REPORT_OUTPUT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", Heep_File[0][i]))[0], ':')[1]);
					PSF = stod(split(FindString("SW_Ps2_factor", Form("%s/REPORT_OUTPUT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", path4_c.c_str(), Heep_File[target_id][i]))[0], ':')[1]);
					//areal_den = stod(split(FindString("target_areal_density [g/cm2]", Form("%s/REPORT_OUPUT/cafe_prod_LH2_heep_coin_report_%s_-1.txt", path4_c.c_str(), Heep_File[target_id][i]))[0], ':')[1]);
				}//target_thickness Count / (g/cm^2 * C *  
			}
		}

		//--------------------------------------------------------
		//Set weight
		//--------------------------------------------------------
		if (user_s_or_c == "coin")
		{
			//weight = /*(1 - (6.0 / 7.0) * (1 - correction))*/ 1 / ((Ehtrk * Eetrk * Emtrk * Elt) * (/*charge * areal_den * transparency * Z_per_A **/ Prot_Abs));
			weight = 1 / (Ehtrk * Eetrk * Emtrk * Elt * Prot_Abs);
			//weight = 1 / (Eff_total * (charge * areal_den * transparency * Z_per_A));
			cout << "Ehtrk: " << Ehtrk << endl;
			cout << "Eetrk: " << Eetrk << endl;
			cout << "Emtrk: " << Emtrk << endl;
			cout << "Elt: " << Elt << endl;
			//cout << "Eff_total: " << Eff_total << endl;
			cout << "charge: " << charge << endl;
			//cout << "correction: " << correction << endl;
			//cout << "areal_den: " << areal_den << endl;
			//cout << "transparency: " << transparency << endl;
			cout << "Proton Absorbtion: " << Prot_Abs << endl;
			//cout << "z/a: " << Z_per_A << endl;
			//weight2 = 1;//
			//cout << "Weight2: " << weight2 << endl;
			//cout << "offset" << ctime_offset << endl;
		}
		else if (user_s_or_c == "sing")
		{
			//weight = 1 / (Eff_total * (charge * areal_den * transparency * Z_per_A)) * PSF;
			weight = PSF / ((Eetrk * Emtrk * Elt));///////////
			//weight = 1 / (Eff_total * (charge * areal_den * transparency * Z_per_A));
			cout << "Eetrk: " << Eetrk << endl;
			cout << "Emtrk: " << Emtrk << endl;
			cout << "Elt: " << Elt << endl;
			cout << "PSF: " << PSF << endl;
			cout << "charge: " << charge << endl;
			//cout << "Eff_total: " << Eff_total << endl;
			//cout << "charge: " << charge << endl;
			//cout << "correction: " << correction << endl;
			//cout << "areal_den: " << areal_den << endl;
			//cout << "transparency: " << transparency << endl;
			//cout << "z/a: " << Z_per_A << endl;
			cout << "Weight: " << weight << endl;
			//weight2 = 1;//
			//cout << "Weight2: " << weight2 << endl;
		}

		cout << "Weight: " << weight << endl;
		
		
		
		nentries = inputtree->GetEntries();		//Get the total number of entries
		//cout << "nentries: " << nentries << endl;



		for (int k = 0; k <= nentries; k++)//nentries
		{
			inputtree->GetEntry(k);				//Get the ith entry from the T TTree
			Em_src = nu - Tx - (sqrt(MN * MN + pmiss * pmiss) - MN);
			ztar_diff = htarz - etarz;
			th_p = xangle - th_e;
			//Pm_par = pmiss*cos(cth_rq);
			//cout << "pmiss:"<< pmiss <<endl;
			//cout << "thrq:" << th_rq << endl;
			//cout << "cthrq:" << cos(th_rq) << endl;
			//Pm_per = sqrt(pmiss*pmiss - Pm_par*Pm_par);
			//alpha_n_v = (emiss - Pm_par) / MP; //eerecoil
			//alpha_v = 2 - alpha_n_v;
			
			v3_ebeam.SetXYZ(0, 0, 10.5494);
			v4_ebeam.SetPxPyPzE(0, 0, 10.5494, 10.5494);
			v3_eprime.SetXYZ(p4x, p4y, p4z);
			v4_eprime.SetVect(v3_eprime);
			v4_eprime.SetE(p4e);
			//cout << "p4x: " << p4x << endl;
			//cout << "p4y: " << p4y << endl;
			//cout << "p4z: " << p4z << endl;
			//cout << "p4e: " << p4e << endl;
			//v3_q = v3_ebeam - v3_eprime;
			//v4_q = v4_ebeam - v4_eprime;

			v3_p.SetXYZ(h4x, h4y, h4z);
			v4_p.SetVect(v3_p);
			v4_p.SetE(h4e);
			//cout << "h4x: " << h4x << endl;
			//cout << "h4y: " << h4y << endl;
			//cout << "h4z: " << h4z << endl;
			//cout << "h4e: " << h4e << endl;

			//v3_q = v3_ebeam - v3_eprime;
			//v4_q = v4_ebeam - v4_eprime;
			
			v3_q.SetXYZ(q_x,q_y,q_z);//= v3_ebeam - v3_eprime;
			v4_q.SetVect(v3_q);// = v4_ebeam - v4_eprime;
			v4_q.SetE(nu);//q3m

			v3_pr.SetXYZ(prec_x, prec_y, prec_z);//v3_q - v3_p;
			v4_pr.SetVect(v3_pr);//v4_q - v4_p;

			//v3_q.SetXYZ(q_x,q_y,q_z);//= v3_ebeam - v3_eprime;
			//v4_q.SetVect(v3_q);// = v4_ebeam - v4_eprime;
			//v4_q.SetE(nu);

			//setmagq

			Er = sqrt(v3_pr.Mag() * v3_pr.Mag() + MP * MP);
			Pm_par = v3_pr.Mag() * cos(cth_rq);//fix
			Pm_per = sqrt(v3_pr.Mag() * v3_pr.Mag() - Pm_par * Pm_par);//fix

			alpha_v = (Er - Pm_par) / MP;//+?

			//TH1F* H1_Pm_par_pid_acc_kin = Histo::mk_Pm_par(inputtree, Pm_par, "H1_Pm_par_pid_acc_kin", Sec_Kin, binng);
			//TH1F* H1_Pm_per_pid_acc_kin = Histo::mk_Pm_per(inputtree, Pm_per, "H1_Pm_pe_pid_acc_kin", Sec_Kin, binning);
			//TH1F* H1_alpha_pid_acc_kin = Histo::mk_alpha(inputtree, alpha_v, "H1_alpha_pid_acc_kin", Sec_Kin, binning);
			//mk_alpha_Pm_per(TTree* input, const char* name, TList* list, int bin)
			

			//--------------------------------------------------------
			//Cuts
			//--------------------------------------------------------
			ctime_L = false;// Cuts::cTime_Rand_L(ep_ctime - ctime_offset);
			ctime_R = Cuts::cTime_Rand_R(ep_ctime - ctime_offset);
			HMS_TRK = true;// Cuts::HMS_Tr(hdc_ntrk, hScinGood, hCerNpeSum, hCalEtotNorm, hHodBetaNtrk);
			//problematic!!!!!!!!!!!
			SHMS_TRK = true;// Cuts::SHMS_Tr(pdc_ntrk, pScinGood, pNGCerNpeSum, pHGCerNpeSum, pCalEtotNorm, pHodBetaNtrk);
			//problematic!!!!!!!!!!!
			ctime = Cuts::cTime(ep_ctime - ctime_offset);
			pid_hms = Cuts::PID_HMS(hCalEtotTrkNorm, hCerNpeSum);
			pid_shms = Cuts::PID_SHMS(pCalEtotTrkNorm, pNGCerNpeSum, pHGCerNpeSum, 1);
			accp_hms = Cuts::Accp_HMS(hdelta, hxptar, hyptar, hms_Coll_gCut->IsInside(hYColl, hXColl));
			accp_shms = Cuts::Accp_SHMS(edelta, exptar, eyptar, shms_Coll_gCut->IsInside(eYColl, eXColl));
			contam_cut = contam_gCut->IsInside(emiss, pmiss);
			

			if (user_target != "LD2")
			{
				accp_z = true;// Cuts::Cut(-2, etarz, 2);// true;// Cuts::Accp_Z(ztar_diff);
			}
			else
			{
				accp_z = true;// Cuts::Cut(-5, etarz, 5);//true;
			}


			if (user_s_or_c == "sing")//sing
			{
				if (evtyp == 1 || evtyp == 3 || evtyp == 5 || evtyp == 7)
				{
					bevtyp = true;
				}
				else
				{
					bevtyp = false;
				}
			}
			else if (user_s_or_c == "coin")//cooin
			{
				if (evtyp >= 4)
				{
					bevtyp = true;
				}
				else
				{
					bevtyp = false;
				}
			}
			else{ cout << "evtyp error" << endl; }
				

			if (user_runtype == "SRC" && user_s_or_c != "sing")
			{
				kin = Cuts::kinSRC(Q2, pmiss, x_bj, ((cth_rq / dtr)), emiss, Em_src, user_target);
			}
			else if (user_runtype == "MF" && user_s_or_c != "sing")
			{
				kin = Cuts::kinMF(Q2, pmiss, emiss, ((cth_rq / dtr)), user_target);
			}
			else if (user_runtype == "Heep" && user_s_or_c != "sing")
			{ 
				kin = Cuts::kinHeepCoin(Q2, W, x_bj, emiss, mmiss); 
			}
			else if (user_s_or_c == "sing")
			{
				pid_hms = true;
				//pid_shms = true;
				Real_Golden = true;
				accp_hms = true;
				//accp_shms = true;
				accp_z = true;////
				kin = true;// Cuts::Cut(0.9, W, 1.05);
				//bevtyp
				ctime = true;
				HMS_TRK = true;
				SHMS_TRK = true;
			}
			else
			{
				cout << "Invalid user run type." << endl;
			}


			//--------------------------------------
			//PID Histograms Bins
			//--------------------------------------
			//----- HMS -----
			H1_hdc_ntrk->Fill(hdc_ntrk, weight);
			H1_hScinGood->Fill(hScinGood, weight);
			H1_pdc_ntrk->Fill(pdc_ntrk, weight);
			H1_pScinGood->Fill(pScinGood, weight);

			H1_hCerNpeSum->Fill(hCerNpeSum, weight);
			H1_hCalEtotNorm->Fill(hCalEtotNorm, weight);
			H1_hCalEtotTrkNorm->Fill(hCalEtotTrkNorm, weight);
			H1_hHodBetaNtrk->Fill(hHodBetaNtrk, weight);
			H1_hHodBetaTrk->Fill(hHodBetaTrk, weight);
			
			
			H1_pNGCerNpeSum->Fill(pNGCerNpeSum, weight);
			H1_pHGCerNpeSum->Fill(pHGCerNpeSum, weight);
			H1_pCalEtotNorm->Fill(pCalEtotNorm, weight);
			H1_pCalEtotTrkNorm->Fill(pCalEtotTrkNorm, weight);
			H1_pHodBetaNtrk->Fill(pHodBetaNtrk, weight);
			H1_pHodBetaTrk->Fill(pHodBetaTrk, weight);
			H1_ep_ctime->Fill(ep_ctime - ctime_offset, weight);

			//----- SHMS -----
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_hdc_ntrk_pid->Fill(hdc_ntrk, weight);
				H1_hScinGood_pid->Fill(hScinGood, weight);
				H1_pdc_ntrk_pid->Fill(pdc_ntrk, weight);
				H1_pScinGood_pid->Fill(pScinGood, weight);
				
				H1_hCerNpeSum_pid->Fill(hCerNpeSum, weight);
				H1_hCalEtotNorm_pid->Fill(hCalEtotNorm, weight);
				H1_hCalEtotTrkNorm_pid->Fill(hCalEtotTrkNorm, weight);
				H1_hHodBetaNtrk_pid->Fill(hHodBetaNtrk, weight);
				H1_hHodBetaTrk_pid->Fill(hHodBetaTrk, weight);
				H1_pNGCerNpeSum_pid->Fill(pNGCerNpeSum, weight);
				H1_pHGCerNpeSum_pid->Fill(pHGCerNpeSum, weight);
				H1_pCalEtotNorm_pid->Fill(pCalEtotNorm, weight);
				H1_pCalEtotTrkNorm_pid->Fill(pCalEtotTrkNorm, weight);
				H1_pHodBetaNtrk_pid->Fill(pHodBetaNtrk, weight);
				H1_pHodBetaTrk_pid->Fill(pHodBetaTrk, weight);
				H1_ep_ctime_pid->Fill(ep_ctime - ctime_offset, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && bevtyp && HMS_TRK && SHMS_TRK)
			{
				H1_ep_ctime_pid_alt->Fill(ep_ctime - ctime_offset, weight);
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_hdc_ntrk_pid_acc->Fill(hdc_ntrk, weight);
				H1_hScinGood_pid_acc->Fill(hScinGood, weight);
				H1_pdc_ntrk_pid_acc->Fill(pdc_ntrk, weight);
				H1_pScinGood_pid_acc->Fill(pScinGood, weight);
				
				H1_hCerNpeSum_pid_acc->Fill(hCerNpeSum, weight);
				H1_hCalEtotNorm_pid_acc->Fill(hCalEtotNorm, weight);
				H1_hCalEtotTrkNorm_pid_acc->Fill(hCalEtotTrkNorm, weight);
				H1_hHodBetaNtrk_pid_acc->Fill(hHodBetaNtrk, weight);
				H1_hHodBetaTrk_pid_acc->Fill(hHodBetaTrk, weight);
				H1_pNGCerNpeSum_pid_acc->Fill(pNGCerNpeSum, weight);
				H1_pHGCerNpeSum_pid_acc->Fill(pHGCerNpeSum, weight);
				H1_pCalEtotNorm_pid_acc->Fill(pCalEtotNorm, weight);
				H1_pCalEtotTrkNorm_pid_acc->Fill(pCalEtotTrkNorm, weight);
				H1_pHodBetaNtrk_pid_acc->Fill(pHodBetaNtrk, weight);
				H1_pHodBetaTrk_pid_acc->Fill(pHodBetaTrk, weight);
				H1_ep_ctime_pid_acc->Fill(ep_ctime - ctime_offset, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && HMS_TRK && SHMS_TRK)
			{
				H1_ep_ctime_pid_acc_alt->Fill(ep_ctime - ctime_offset, weight);
			}
			

			if (user_runtype == "SRC")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_hdc_ntrk_pid_acc_1->Fill(hdc_ntrk, weight);
					H1_hScinGood_pid_acc_1->Fill(hScinGood, weight);
					H1_pdc_ntrk_pid_acc_1->Fill(pdc_ntrk, weight);
					H1_pScinGood_pid_acc_1->Fill(pScinGood, weight); 
					//pid_shms && !contam_cut
					H1_hCerNpeSum_pid_acc_1->Fill(hCerNpeSum, weight);
					H1_hCalEtotNorm_pid_acc_1->Fill(hCalEtotNorm, weight);
					H1_hCalEtotTrkNorm_pid_acc_1->Fill(hCalEtotTrkNorm, weight);
					H1_hHodBetaNtrk_pid_acc_1->Fill(hHodBetaNtrk, weight);
					H1_hHodBetaTrk_pid_acc_1->Fill(hHodBetaTrk, weight);
					H1_pNGCerNpeSum_pid_acc_1->Fill(pNGCerNpeSum, weight);
					H1_pHGCerNpeSum_pid_acc_1->Fill(pHGCerNpeSum, weight);
					H1_pCalEtotNorm_pid_acc_1->Fill(pCalEtotNorm, weight);
					H1_pCalEtotTrkNorm_pid_acc_1->Fill(pCalEtotTrkNorm, weight);
					H1_pHodBetaNtrk_pid_acc_1->Fill(pHodBetaNtrk, weight);
					H1_pHodBetaTrk_pid_acc_1->Fill(pHodBetaTrk, weight);
					H1_ep_ctime_pid_acc_1->Fill(ep_ctime - ctime_offset, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_hdc_ntrk_pid_acc_2->Fill(hdc_ntrk, weight);
					H1_hScinGood_pid_acc_2->Fill(hScinGood, weight);
					H1_pdc_ntrk_pid_acc_2->Fill(pdc_ntrk, weight);
					H1_pScinGood_pid_acc_2->Fill(pScinGood, weight); 
					
					H1_hCerNpeSum_pid_acc_2->Fill(hCerNpeSum, weight);
					H1_hCalEtotNorm_pid_acc_2->Fill(hCalEtotNorm, weight);
					H1_hCalEtotTrkNorm_pid_acc_2->Fill(hCalEtotTrkNorm, weight);
					H1_hHodBetaNtrk_pid_acc_2->Fill(hHodBetaNtrk, weight);
					H1_hHodBetaTrk_pid_acc_2->Fill(hHodBetaTrk, weight);
					H1_pNGCerNpeSum_pid_acc_2->Fill(pNGCerNpeSum, weight);
					H1_pHGCerNpeSum_pid_acc_2->Fill(pHGCerNpeSum, weight);
					H1_pCalEtotNorm_pid_acc_2->Fill(pCalEtotNorm, weight);
					H1_pCalEtotTrkNorm_pid_acc_2->Fill(pCalEtotTrkNorm, weight);
					H1_pHodBetaNtrk_pid_acc_2->Fill(pHodBetaNtrk, weight);
					H1_pHodBetaTrk_pid_acc_2->Fill(pHodBetaTrk, weight);
					H1_ep_ctime_pid_acc_2->Fill(ep_ctime - ctime_offset, weight);
				}
			}
			else if (user_runtype == "MF")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_hdc_ntrk_pid_acc_1->Fill(hdc_ntrk, weight);
					H1_hScinGood_pid_acc_1->Fill(hScinGood, weight);
					H1_pdc_ntrk_pid_acc_1->Fill(pdc_ntrk, weight);
					H1_pScinGood_pid_acc_1->Fill(pScinGood, weight); 
					
					H1_hCerNpeSum_pid_acc_1->Fill(hCerNpeSum, weight);
					H1_hCalEtotNorm_pid_acc_1->Fill(hCalEtotNorm, weight);
					H1_hCalEtotTrkNorm_pid_acc_1->Fill(hCalEtotTrkNorm, weight);
					H1_hHodBetaNtrk_pid_acc_1->Fill(hHodBetaNtrk, weight);
					H1_hHodBetaTrk_pid_acc_1->Fill(hHodBetaTrk, weight);
					H1_pNGCerNpeSum_pid_acc_1->Fill(pNGCerNpeSum, weight);
					H1_pHGCerNpeSum_pid_acc_1->Fill(pHGCerNpeSum, weight);
					H1_pCalEtotNorm_pid_acc_1->Fill(pCalEtotNorm, weight);
					H1_pCalEtotTrkNorm_pid_acc_1->Fill(pCalEtotTrkNorm, weight);
					H1_pHodBetaNtrk_pid_acc_1->Fill(pHodBetaNtrk, weight);
					H1_pHodBetaTrk_pid_acc_1->Fill(pHodBetaTrk, weight);
					H1_ep_ctime_pid_acc_1->Fill(ep_ctime - ctime_offset, weight);

				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_hdc_ntrk_pid_acc_2->Fill(hdc_ntrk, weight);
					H1_hScinGood_pid_acc_2->Fill(hScinGood, weight);
					H1_pdc_ntrk_pid_acc_2->Fill(pdc_ntrk, weight);
					H1_pScinGood_pid_acc_2->Fill(pScinGood, weight); 
					
					H1_hCerNpeSum_pid_acc_2->Fill(hCerNpeSum, weight);
					H1_hCalEtotNorm_pid_acc_2->Fill(hCalEtotNorm, weight);
					H1_hCalEtotTrkNorm_pid_acc_2->Fill(hCalEtotTrkNorm, weight);
					H1_hHodBetaNtrk_pid_acc_2->Fill(hHodBetaNtrk, weight);
					H1_hHodBetaTrk_pid_acc_2->Fill(hHodBetaTrk, weight);
					H1_pNGCerNpeSum_pid_acc_2->Fill(pNGCerNpeSum, weight);
					H1_pHGCerNpeSum_pid_acc_2->Fill(pHGCerNpeSum, weight);
					H1_pCalEtotNorm_pid_acc_2->Fill(pCalEtotNorm, weight);
					H1_pCalEtotTrkNorm_pid_acc_2->Fill(pCalEtotTrkNorm, weight);
					H1_pHodBetaNtrk_pid_acc_2->Fill(pHodBetaNtrk, weight);
					H1_pHodBetaTrk_pid_acc_2->Fill(pHodBetaTrk, weight);
					H1_ep_ctime_pid_acc_2->Fill(ep_ctime - ctime_offset, weight);
				}
			}
			else if (user_runtype == "Heep")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_hdc_ntrk_pid_acc_1->Fill(hdc_ntrk, weight);
					H1_hScinGood_pid_acc_1->Fill(hScinGood, weight);
					H1_pdc_ntrk_pid_acc_1->Fill(pdc_ntrk, weight);
					H1_pScinGood_pid_acc_1->Fill(pScinGood, weight); 
					
					H1_hCerNpeSum_pid_acc_1->Fill(hCerNpeSum, weight);
					H1_hCalEtotNorm_pid_acc_1->Fill(hCalEtotNorm, weight);
					H1_hCalEtotTrkNorm_pid_acc_1->Fill(hCalEtotTrkNorm, weight);
					H1_hHodBetaNtrk_pid_acc_1->Fill(hHodBetaNtrk, weight);
					H1_hHodBetaTrk_pid_acc_1->Fill(hHodBetaTrk, weight);
					H1_pNGCerNpeSum_pid_acc_1->Fill(pNGCerNpeSum, weight);
					H1_pHGCerNpeSum_pid_acc_1->Fill(pHGCerNpeSum, weight);
					H1_pCalEtotNorm_pid_acc_1->Fill(pCalEtotNorm, weight);
					H1_pCalEtotTrkNorm_pid_acc_1->Fill(pCalEtotTrkNorm, weight);
					H1_pHodBetaNtrk_pid_acc_1->Fill(pHodBetaNtrk, weight);
					H1_pHodBetaTrk_pid_acc_1->Fill(pHodBetaTrk, weight);
					H1_ep_ctime_pid_acc_1->Fill(ep_ctime - ctime_offset, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1), Cuts::get_heep_xbj_cut(0), x_bj, Cuts::get_heep_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_hdc_ntrk_pid_acc_2->Fill(hdc_ntrk, weight);
					H1_hScinGood_pid_acc_2->Fill(hScinGood, weight);
					H1_pdc_ntrk_pid_acc_2->Fill(pdc_ntrk, weight);
					H1_pScinGood_pid_acc_2->Fill(pScinGood, weight); 
					
					H1_hCerNpeSum_pid_acc_2->Fill(hCerNpeSum, weight);
					H1_hCalEtotNorm_pid_acc_2->Fill(hCalEtotNorm, weight);
					H1_hCalEtotTrkNorm_pid_acc_2->Fill(hCalEtotTrkNorm, weight);
					H1_hHodBetaNtrk_pid_acc_2->Fill(hHodBetaNtrk, weight);
					H1_hHodBetaTrk_pid_acc_2->Fill(hHodBetaTrk, weight);
					H1_pNGCerNpeSum_pid_acc_2->Fill(pNGCerNpeSum, weight);
					H1_pHGCerNpeSum_pid_acc_2->Fill(pHGCerNpeSum, weight);
					H1_pCalEtotNorm_pid_acc_2->Fill(pCalEtotNorm, weight);
					H1_pCalEtotTrkNorm_pid_acc_2->Fill(pCalEtotTrkNorm, weight);
					H1_pHodBetaNtrk_pid_acc_2->Fill(pHodBetaNtrk, weight);
					H1_pHodBetaTrk_pid_acc_2->Fill(pHodBetaTrk, weight);
					H1_ep_ctime_pid_acc_2->Fill(ep_ctime - ctime_offset, weight);
				}
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			//if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && (evtyp == 1 || evtyp == 3 || evtyp == 5 || evtyp == 7) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_hdc_ntrk_pid_acc_kin->Fill(hdc_ntrk, weight);
				H1_hScinGood_pid_acc_kin->Fill(hScinGood, weight);
				H1_pdc_ntrk_pid_acc_kin->Fill(pdc_ntrk, weight);
				H1_pScinGood_pid_acc_kin->Fill(pScinGood, weight); 
				
				H1_hCerNpeSum_pid_acc_kin->Fill(hCerNpeSum, weight);
				H1_hCalEtotNorm_pid_acc_kin->Fill(hCalEtotNorm, weight);
				H1_hCalEtotTrkNorm_pid_acc_kin->Fill(hCalEtotTrkNorm, weight);
				H1_hHodBetaNtrk_pid_acc_kin->Fill(hHodBetaNtrk, weight);
				H1_hHodBetaTrk_pid_acc_kin->Fill(hHodBetaTrk, weight);
				H1_pNGCerNpeSum_pid_acc_kin->Fill(pNGCerNpeSum, weight);
				H1_pHGCerNpeSum_pid_acc_kin->Fill(pHGCerNpeSum, weight);
				//H1_pCalEtotNorm_pid_acc_kin->Fill(pCalEtotNorm, weight);
				if (evtyp == 4)
				{
					//H1_pCalEtotTrkNorm_pid_acc_kin->Fill(pCalEtotTrkNorm, weight);
				}
				H1_pHodBetaNtrk_pid_acc_kin->Fill(pHodBetaNtrk, weight);
				H1_pHodBetaTrk_pid_acc_kin->Fill(pHodBetaTrk, weight);
				H1_ep_ctime_pid_acc_kin->Fill(ep_ctime - ctime_offset, weight);
			}
			if (!contam_cut && pid_hms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_pCalEtotTrkNorm_pid_acc_kin->Fill(pCalEtotTrkNorm, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && HMS_TRK && SHMS_TRK)
			{
				H1_ep_ctime_pid_acc_kin_alt->Fill(ep_ctime - ctime_offset, weight);
				H1_ep_ctime_ratio_pid_acc_kin_full->Fill(ep_ctime - ctime_offset, weight);
				H1_ep_ctime_ratio_pid_acc_kin_full_mf->Fill(ep_ctime - ctime_offset, weight);
			}
			/*if (
				!contam_cut &&
				0.8 < pCalEtotNorm && pCalEtotNorm < 1.3 &&
				Real_Golden &&
				-10 < hdelta && hdelta < 10 &&
				hms_Coll_gCut->IsInside(hYColl, hXColl) &&
				0 < edelta && edelta < 22 &&
				shms_Coll_gCut->IsInside(eYColl, eXColl) &&
				Cuts::Cut(
					Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1),
					Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1),
					Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1),
					Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1)) &&
				bevtyp && HMS_TRK && SHMS_TRK
				)
			{

			}*/


			/*
			H1_ep_ctime->Fill(ep_ctime - ctime_offset, weight);

			if (!contam_cut && pid_hms && pid_shms && Real_Golden)
			{
				H1_ep_ctime_pid->Fill(ep_ctime - ctime_offset, weight);
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z)
			{
				H1_ep_ctime_pid_acc->Fill(ep_ctime - ctime_offset, weight);
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin)
			{
				H1_ep_ctime_pid_acc_kin->Fill(ep_ctime - ctime_offset, weight);
			}
			*/
			

			//--------------------------------------
			//Acceptance Histogram Bins
			//--------------------------------------

			//----- Hadron Arm Focal Plane -----
			H1_hxfp->Fill(hxfp, weight);
			H1_hxpfp->Fill(hxpfp / dtr, weight);
			H1_hyfp->Fill(hyfp, weight);
			H1_hypfp->Fill(hypfp / dtr, weight);
			//----- Hadron Arm Reconstructed Quantities -----
			H1_hytar->Fill(hytar, weight);
			H1_hyptar->Fill(hyptar / dtr, weight);
			H1_hxptar->Fill(hxptar / dtr, weight);
			H1_hdelta->Fill(hdelta, weight);
			//----- Target Reconstruction (Hall Coord. System) -----
			H1_htarx->Fill(htarx, weight);
			H1_htary->Fill(htary, weight);
			H1_htarz->Fill(htarz, weight);
			//----- HMS Collimator -----
			H1_hXColl->Fill(hXColl, weight);
			H1_hYColl->Fill(hYColl, weight);
			//----- Hadron Arm Focal Plane -----
			H1_exfp->Fill(exfp, weight);
			H1_expfp->Fill(expfp / dtr, weight);
			H1_eyfp->Fill(eyfp, weight);
			H1_eypfp->Fill(eypfp / dtr, weight);
			//----- Hadron Arm Reconstructed Quantities -----
			H1_eytar->Fill(eytar, weight);
			H1_eyptar->Fill(eyptar / dtr, weight);
			H1_exptar->Fill(exptar / dtr, weight);
			H1_edelta->Fill(edelta, weight);
			//----- Target Reconstruction (Hall Coord. System) -----
			H1_etarx->Fill(etarx, weight);
			H1_etary->Fill(etary, weight);
			H1_etarz->Fill(etarz, weight);
			//----- SHMS Collimator -----
			H1_eXColl->Fill(eXColl, weight);
			H1_eYColl->Fill(eYColl, weight);
			//----- -----
			H1_ztar_diff->Fill(ztar_diff, weight);

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				//----- Hadron Arm Focal Plane -----
				H1_hxfp_pid->Fill(hxfp, weight);
				H1_hxpfp_pid->Fill(hxpfp / dtr, weight);
				H1_hyfp_pid->Fill(hyfp, weight);
				H1_hypfp_pid->Fill(hypfp / dtr, weight);
				//----- Hadron Arm Reconstructed Quantities -----
				H1_hytar_pid->Fill(hytar, weight);
				H1_hyptar_pid->Fill(hyptar / dtr, weight);
				H1_hxptar_pid->Fill(hxptar / dtr, weight);
				H1_hdelta_pid->Fill(hdelta, weight);
				//----- Target Reconstruction (Hall Coord. System) -----
				H1_htarx_pid->Fill(htarx, weight);
				H1_htary_pid->Fill(htary, weight);
				H1_htarz_pid->Fill(htarz, weight);
				//----- HMS Collimator -----
				H1_hXColl_pid->Fill(hXColl, weight);
				H1_hYColl_pid->Fill(hYColl, weight);
				//----- Hadron Arm Focal Plane -----
				H1_exfp_pid->Fill(exfp, weight);
				H1_expfp_pid->Fill(expfp / dtr, weight);
				H1_eyfp_pid->Fill(eyfp, weight);
				H1_eypfp_pid->Fill(eypfp / dtr, weight);
				//----- Hadron Arm Reconstructed Quantities -----
				H1_eytar_pid->Fill(eytar, weight);
				H1_eyptar_pid->Fill(eyptar / dtr, weight);
				H1_exptar_pid->Fill(exptar / dtr, weight);
				H1_edelta_pid->Fill(edelta, weight);
				//----- Target Reconstruction (Hall Coord. System) -----
				H1_etarx_pid->Fill(etarx, weight);
				H1_etary_pid->Fill(etary, weight);
				H1_etarz_pid->Fill(etarz, weight);
				//----- SHMS Collimator -----
				H1_eXColl_pid->Fill(eXColl, weight);
				H1_eYColl_pid->Fill(eYColl, weight);
				//----- -----
				H1_ztar_diff_pid->Fill(ztar_diff, weight);
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				//----- Hadron Arm Focal Plane -----
				H1_hxfp_pid_acc->Fill(hxfp, weight);
				H1_hxpfp_pid_acc->Fill(hxpfp / dtr, weight);
				H1_hyfp_pid_acc->Fill(hyfp, weight);
				H1_hypfp_pid_acc->Fill(hypfp / dtr, weight);
				//----- Hadron Arm Reconstructed Quantities -----
				H1_hytar_pid_acc->Fill(hytar, weight);
				H1_hyptar_pid_acc->Fill(hyptar / dtr, weight);
				H1_hxptar_pid_acc->Fill(hxptar / dtr, weight);
				H1_hdelta_pid_acc->Fill(hdelta, weight);
				//----- Target Reconstruction (Hall Coord. System) -----
				H1_htarx_pid_acc->Fill(htarx, weight);
				H1_htary_pid_acc->Fill(htary, weight);
				H1_htarz_pid_acc->Fill(htarz, weight);
				//----- HMS Collimator -----
				H1_hXColl_pid_acc->Fill(hXColl, weight);
				H1_hYColl_pid_acc->Fill(hYColl, weight);
				//----- Hadron Arm Focal Plane -----
				H1_exfp_pid_acc->Fill(exfp, weight);
				H1_expfp_pid_acc->Fill(expfp / dtr, weight);
				H1_eyfp_pid_acc->Fill(eyfp, weight);
				H1_eypfp_pid_acc->Fill(eypfp / dtr, weight);
				//----- Hadron Arm Reconstructed Quantities -----
				H1_eytar_pid_acc->Fill(eytar, weight);
				H1_eyptar_pid_acc->Fill(eyptar / dtr, weight);
				H1_exptar_pid_acc->Fill(exptar / dtr, weight);
				H1_edelta_pid_acc->Fill(edelta, weight);
				//----- Target Reconstruction (Hall Coord. System) -----
				H1_etarx_pid_acc->Fill(etarx, weight);
				H1_etary_pid_acc->Fill(etary, weight);
				H1_etarz_pid_acc->Fill(etarz, weight);
				//----- SHMS Collimator -----
				H1_eXColl_pid_acc->Fill(eXColl, weight);
				H1_eYColl_pid_acc->Fill(eYColl, weight);
				//----- -----
				H1_ztar_diff_pid_acc->Fill(ztar_diff, weight);
			}


			if (user_runtype == "SRC")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					//----- Hadron Arm Focal Plane -----
					H1_hxfp_pid_acc_1->Fill(hxfp, weight);
					H1_hxpfp_pid_acc_1->Fill(hxpfp / dtr, weight);
					H1_hyfp_pid_acc_1->Fill(hyfp, weight);
					H1_hypfp_pid_acc_1->Fill(hypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_hytar_pid_acc_1->Fill(hytar, weight);
					H1_hyptar_pid_acc_1->Fill(hyptar / dtr, weight);
					H1_hxptar_pid_acc_1->Fill(hxptar / dtr, weight);
					H1_hdelta_pid_acc_1->Fill(hdelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_htarx_pid_acc_1->Fill(htarx, weight);
					H1_htary_pid_acc_1->Fill(htary, weight);
					H1_htarz_pid_acc_1->Fill(htarz, weight);
					//----- HMS Collimator -----
					H1_hXColl_pid_acc_1->Fill(hXColl, weight);
					H1_hYColl_pid_acc_1->Fill(hYColl, weight);
					//----- Hadron Arm Focal Plane -----
					H1_exfp_pid_acc_1->Fill(exfp, weight);
					H1_expfp_pid_acc_1->Fill(expfp / dtr, weight);
					H1_eyfp_pid_acc_1->Fill(eyfp, weight);
					H1_eypfp_pid_acc_1->Fill(eypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_eytar_pid_acc_1->Fill(eytar, weight);
					H1_eyptar_pid_acc_1->Fill(eyptar / dtr, weight);
					H1_exptar_pid_acc_1->Fill(exptar / dtr, weight);
					H1_edelta_pid_acc_1->Fill(edelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_etarx_pid_acc_1->Fill(etarx, weight);
					H1_etary_pid_acc_1->Fill(etary, weight);
					H1_etarz_pid_acc_1->Fill(etarz, weight);
					//----- SHMS Collimator -----
					H1_eXColl_pid_acc_1->Fill(eXColl, weight);
					H1_eYColl_pid_acc_1->Fill(eYColl, weight);
					//----- -----
					H1_ztar_diff_pid_acc_1->Fill(ztar_diff, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					//----- Hadron Arm Focal Plane -----
					H1_hxfp_pid_acc_2->Fill(hxfp, weight);
					H1_hxpfp_pid_acc_2->Fill(hxpfp / dtr, weight);
					H1_hyfp_pid_acc_2->Fill(hyfp, weight);
					H1_hypfp_pid_acc_2->Fill(hypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_hytar_pid_acc_2->Fill(hytar, weight);
					H1_hyptar_pid_acc_2->Fill(hyptar / dtr, weight);
					H1_hxptar_pid_acc_2->Fill(hxptar / dtr, weight);
					H1_hdelta_pid_acc_2->Fill(hdelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_htarx_pid_acc_2->Fill(htarx, weight);
					H1_htary_pid_acc_2->Fill(htary, weight);
					H1_htarz_pid_acc_2->Fill(htarz, weight);
					//----- HMS Collimator -----
					H1_hXColl_pid_acc_2->Fill(hXColl, weight);
					H1_hYColl_pid_acc_2->Fill(hYColl, weight);
					//----- Hadron Arm Focal Plane -----
					H1_exfp_pid_acc_2->Fill(exfp, weight);
					H1_expfp_pid_acc_2->Fill(expfp / dtr, weight);
					H1_eyfp_pid_acc_2->Fill(eyfp, weight);
					H1_eypfp_pid_acc_2->Fill(eypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_eytar_pid_acc_2->Fill(eytar, weight);
					H1_eyptar_pid_acc_2->Fill(eyptar / dtr, weight);
					H1_exptar_pid_acc_2->Fill(exptar / dtr, weight);
					H1_edelta_pid_acc_2->Fill(edelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_etarx_pid_acc_2->Fill(etarx, weight);
					H1_etary_pid_acc_2->Fill(etary, weight);
					H1_etarz_pid_acc_2->Fill(etarz, weight);
					//----- SHMS Collimator -----
					H1_eXColl_pid_acc_2->Fill(eXColl, weight);
					H1_eYColl_pid_acc_2->Fill(eYColl, weight);
					//----- -----
					H1_ztar_diff_pid_acc_2->Fill(ztar_diff, weight);
				}
			}
			else if (user_runtype == "MF")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					//----- Hadron Arm Focal Plane -----
					H1_hxfp_pid_acc_1->Fill(hxfp, weight);
					H1_hxpfp_pid_acc_1->Fill(hxpfp / dtr, weight);
					H1_hyfp_pid_acc_1->Fill(hyfp, weight);
					H1_hypfp_pid_acc_1->Fill(hypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_hytar_pid_acc_1->Fill(hytar, weight);
					H1_hyptar_pid_acc_1->Fill(hyptar / dtr, weight);
					H1_hxptar_pid_acc_1->Fill(hxptar / dtr, weight);
					H1_hdelta_pid_acc_1->Fill(hdelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_htarx_pid_acc_1->Fill(htarx, weight);
					H1_htary_pid_acc_1->Fill(htary, weight);
					H1_htarz_pid_acc_1->Fill(htarz, weight);
					//----- HMS Collimator -----
					H1_hXColl_pid_acc_1->Fill(hXColl, weight);
					H1_hYColl_pid_acc_1->Fill(hYColl, weight);
					//----- Hadron Arm Focal Plane -----
					H1_exfp_pid_acc_1->Fill(exfp, weight);
					H1_expfp_pid_acc_1->Fill(expfp / dtr, weight);
					H1_eyfp_pid_acc_1->Fill(eyfp, weight);
					H1_eypfp_pid_acc_1->Fill(eypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_eytar_pid_acc_1->Fill(eytar, weight);
					H1_eyptar_pid_acc_1->Fill(eyptar / dtr, weight);
					H1_exptar_pid_acc_1->Fill(exptar / dtr, weight);
					H1_edelta_pid_acc_1->Fill(edelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_etarx_pid_acc_1->Fill(etarx, weight);
					H1_etary_pid_acc_1->Fill(etary, weight);
					H1_etarz_pid_acc_1->Fill(etarz, weight);
					//----- SHMS Collimator -----
					H1_eXColl_pid_acc_1->Fill(eXColl, weight);
					H1_eYColl_pid_acc_1->Fill(eYColl, weight);
					//----- -----
					H1_ztar_diff_pid_acc_1->Fill(ztar_diff, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					//----- Hadron Arm Focal Plane -----
					H1_hxfp_pid_acc_2->Fill(hxfp, weight);
					H1_hxpfp_pid_acc_2->Fill(hxpfp / dtr, weight);
					H1_hyfp_pid_acc_2->Fill(hyfp, weight);
					H1_hypfp_pid_acc_2->Fill(hypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_hytar_pid_acc_2->Fill(hytar, weight);
					H1_hyptar_pid_acc_2->Fill(hyptar / dtr, weight);
					H1_hxptar_pid_acc_2->Fill(hxptar / dtr, weight);
					H1_hdelta_pid_acc_2->Fill(hdelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_htarx_pid_acc_2->Fill(htarx, weight);
					H1_htary_pid_acc_2->Fill(htary, weight);
					H1_htarz_pid_acc_2->Fill(htarz, weight);
					//----- HMS Collimator -----
					H1_hXColl_pid_acc_2->Fill(hXColl, weight);
					H1_hYColl_pid_acc_2->Fill(hYColl, weight);
					//----- Hadron Arm Focal Plane -----
					H1_exfp_pid_acc_2->Fill(exfp, weight);
					H1_expfp_pid_acc_2->Fill(expfp / dtr, weight);
					H1_eyfp_pid_acc_2->Fill(eyfp, weight);
					H1_eypfp_pid_acc_2->Fill(eypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_eytar_pid_acc_2->Fill(eytar, weight);
					H1_eyptar_pid_acc_2->Fill(eyptar / dtr, weight);
					H1_exptar_pid_acc_2->Fill(exptar / dtr, weight);
					H1_edelta_pid_acc_2->Fill(edelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_etarx_pid_acc_2->Fill(etarx, weight);
					H1_etary_pid_acc_2->Fill(etary, weight);
					H1_etarz_pid_acc_2->Fill(etarz, weight);
					//----- SHMS Collimator -----
					H1_eXColl_pid_acc_2->Fill(eXColl, weight);
					H1_eYColl_pid_acc_2->Fill(eYColl, weight);
					//----- -----
					H1_ztar_diff_pid_acc_2->Fill(ztar_diff, weight);
				}
			}
			else if (user_runtype == "Heep")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					//----- Hadron Arm Focal Plane -----
					H1_hxfp_pid_acc_1->Fill(hxfp, weight);
					H1_hxpfp_pid_acc_1->Fill(hxpfp / dtr, weight);
					H1_hyfp_pid_acc_1->Fill(hyfp, weight);
					H1_hypfp_pid_acc_1->Fill(hypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_hytar_pid_acc_1->Fill(hytar, weight);
					H1_hyptar_pid_acc_1->Fill(hyptar / dtr, weight);
					H1_hxptar_pid_acc_1->Fill(hxptar / dtr, weight);
					H1_hdelta_pid_acc_1->Fill(hdelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_htarx_pid_acc_1->Fill(htarx, weight);
					H1_htary_pid_acc_1->Fill(htary, weight);
					H1_htarz_pid_acc_1->Fill(htarz, weight);
					//----- HMS Collimator -----
					H1_hXColl_pid_acc_1->Fill(hXColl, weight);
					H1_hYColl_pid_acc_1->Fill(hYColl, weight);
					//----- Hadron Arm Focal Plane -----
					H1_exfp_pid_acc_1->Fill(exfp, weight);
					H1_expfp_pid_acc_1->Fill(expfp / dtr, weight);
					H1_eyfp_pid_acc_1->Fill(eyfp, weight);
					H1_eypfp_pid_acc_1->Fill(eypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_eytar_pid_acc_1->Fill(eytar, weight);
					H1_eyptar_pid_acc_1->Fill(eyptar / dtr, weight);
					H1_exptar_pid_acc_1->Fill(exptar / dtr, weight);
					H1_edelta_pid_acc_1->Fill(edelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_etarx_pid_acc_1->Fill(etarx, weight);
					H1_etary_pid_acc_1->Fill(etary, weight);
					H1_etarz_pid_acc_1->Fill(etarz, weight);
					//----- SHMS Collimator -----
					H1_eXColl_pid_acc_1->Fill(eXColl, weight);
					H1_eYColl_pid_acc_1->Fill(eYColl, weight);
					//----- -----
					H1_ztar_diff_pid_acc_1->Fill(ztar_diff, weight);
				}
				
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1), Cuts::get_heep_xbj_cut(0), x_bj, Cuts::get_heep_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					//----- Hadron Arm Focal Plane -----
					H1_hxfp_pid_acc_2->Fill(hxfp, weight);
					H1_hxpfp_pid_acc_2->Fill(hxpfp / dtr, weight);
					H1_hyfp_pid_acc_2->Fill(hyfp, weight);
					H1_hypfp_pid_acc_2->Fill(hypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_hytar_pid_acc_2->Fill(hytar, weight);
					H1_hyptar_pid_acc_2->Fill(hyptar / dtr, weight);
					H1_hxptar_pid_acc_2->Fill(hxptar / dtr, weight);
					H1_hdelta_pid_acc_2->Fill(hdelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_htarx_pid_acc_2->Fill(htarx, weight);
					H1_htary_pid_acc_2->Fill(htary, weight);
					H1_htarz_pid_acc_2->Fill(htarz, weight);
					//----- HMS Collimator -----
					H1_hXColl_pid_acc_2->Fill(hXColl, weight);
					H1_hYColl_pid_acc_2->Fill(hYColl, weight);
					//----- Hadron Arm Focal Plane -----
					H1_exfp_pid_acc_2->Fill(exfp, weight);
					H1_expfp_pid_acc_2->Fill(expfp / dtr, weight);
					H1_eyfp_pid_acc_2->Fill(eyfp, weight);
					H1_eypfp_pid_acc_2->Fill(eypfp / dtr, weight);
					//----- Hadron Arm Reconstructed Quantities -----
					H1_eytar_pid_acc_2->Fill(eytar, weight);
					H1_eyptar_pid_acc_2->Fill(eyptar / dtr, weight);
					H1_exptar_pid_acc_2->Fill(exptar / dtr, weight);
					H1_edelta_pid_acc_2->Fill(edelta, weight);
					//----- Target Reconstruction (Hall Coord. System) -----
					H1_etarx_pid_acc_2->Fill(etarx, weight);
					H1_etary_pid_acc_2->Fill(etary, weight);
					H1_etarz_pid_acc_2->Fill(etarz, weight);
					//----- SHMS Collimator -----
					H1_eXColl_pid_acc_2->Fill(eXColl, weight);
					H1_eYColl_pid_acc_2->Fill(eYColl, weight);
					//----- -----
					H1_ztar_diff_pid_acc_2->Fill(ztar_diff, weight);
				}
			}


			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				//----- Hadron Arm Focal Plane -----
				H1_hxfp_pid_acc_kin->Fill(hxfp, weight);
				H1_hxpfp_pid_acc_kin->Fill(hxpfp / dtr, weight);
				H1_hyfp_pid_acc_kin->Fill(hyfp, weight);
				H1_hypfp_pid_acc_kin->Fill(hypfp / dtr, weight);
				//----- Hadron Arm Reconstructed Quantities -----
				H1_hytar_pid_acc_kin->Fill(hytar, weight);
				H1_hyptar_pid_acc_kin->Fill(hyptar / dtr, weight);
				H1_hxptar_pid_acc_kin->Fill(hxptar / dtr, weight);
				H1_hdelta_pid_acc_kin->Fill(hdelta, weight);
				//----- Target Reconstruction (Hall Coord. System) -----
				H1_htarx_pid_acc_kin->Fill(htarx, weight);
				H1_htary_pid_acc_kin->Fill(htary, weight);
				H1_htarz_pid_acc_kin->Fill(htarz, weight);
				//----- HMS Collimator -----
				H1_hXColl_pid_acc_kin->Fill(hXColl, weight);
				H1_hYColl_pid_acc_kin->Fill(hYColl, weight);
				//----- Hadron Arm Focal Plane -----
				H1_exfp_pid_acc_kin->Fill(exfp, weight);
				H1_expfp_pid_acc_kin->Fill(expfp / dtr, weight);
				H1_eyfp_pid_acc_kin->Fill(eyfp, weight);
				H1_eypfp_pid_acc_kin->Fill(eypfp / dtr, weight);
				//----- Hadron Arm Reconstructed Quantities -----
				H1_eytar_pid_acc_kin->Fill(eytar, weight);
				H1_eyptar_pid_acc_kin->Fill(eyptar / dtr, weight);
				H1_exptar_pid_acc_kin->Fill(exptar / dtr, weight);
				//H1_edelta_pid_acc_kin->Fill(edelta, weight);
				//----- Target Reconstruction (Hall Coord. System) -----
				H1_etarx_pid_acc_kin->Fill(etarx, weight);
				H1_etary_pid_acc_kin->Fill(etary, weight);
				H1_etarz_pid_acc_kin->Fill(etarz, weight);
				//----- SHMS Collimator -----
				H1_eXColl_pid_acc_kin->Fill(eXColl, weight);
				H1_eYColl_pid_acc_kin->Fill(eYColl, weight);
				//----- -----
				H1_ztar_diff_pid_acc_kin->Fill(ztar_diff, weight);
			}
			

			//--------------------------------------
			//Primary Electron Kinematics
			//--------------------------------------
			H1_Q2->Fill(Q2, weight);
			H1_W->Fill(W, weight);
			H1_nu->Fill(nu, weight);
			H1_ph_q->Fill(ph_q / dtr, weight);
			H1_q->Fill(q, weight);
			H1_qx->Fill(q_x, weight);
			H1_qy->Fill(q_y, weight);
			H1_qz->Fill(q_z, weight);
			H1_th_e->Fill(th_e / dtr, weight);
			H1_th_q->Fill(th_q / dtr, weight);
			H1_xbj->Fill(x_bj, weight);
			H1_xbj_ratio->Fill(x_bj, weight);

			if ((ctime_L || ctime_R))
			{
				H1_Q2_rand->Fill(Q2, weight);
				H1_W_rand->Fill(W, weight);
				H1_nu_rand->Fill(nu, weight);
				H1_q_rand->Fill(q, weight);
				H1_xbj_rand->Fill(x_bj, weight);
				H1_ep_ctime_rand->Fill(ep_ctime - ctime_offset, weight);

				H1_emiss_rand->Fill(emiss, weight);
				H1_pmiss_rand->Fill(pmiss, weight);
				H1_pmiss_raw_rand->Fill(pmiss, weight2);
				H1_pmiss_c_rand->Fill(pmiss, weight2);
				H1_th_rq_rand->Fill(cth_rq / dtr, weight);
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_Q2_pid->Fill(Q2, weight);
				H1_W_pid->Fill(W, weight);
				H1_nu_pid->Fill(nu, weight);
				H1_ph_q_pid->Fill(ph_q / dtr, weight);
				H1_q_pid->Fill(q, weight);
				H1_qx_pid->Fill(q_x, weight);
				H1_qy_pid->Fill(q_y, weight);
				H1_qz_pid->Fill(q_z, weight);
				H1_th_e_pid->Fill(th_e / dtr, weight);
				H1_th_q_pid->Fill(th_q / dtr, weight);
				H1_xbj_pid->Fill(x_bj, weight);
				H1_xbj_ratio_pid->Fill(x_bj, weight);

			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && bevtyp && (ctime_L || ctime_R) && HMS_TRK && SHMS_TRK)
			{
				H1_Q2_pid_rand->Fill(Q2, weight);
				H1_W_pid_rand->Fill(W, weight);
				H1_nu_pid_rand->Fill(nu, weight);
				H1_q_pid_rand->Fill(q, weight);
				H1_xbj_pid_rand->Fill(x_bj, weight);
				H1_ep_ctime_pid_rand->Fill(ep_ctime - ctime_offset, weight);

				H1_emiss_pid_rand->Fill(emiss, weight);
				H1_pmiss_pid_rand->Fill(pmiss, weight);
				H1_pmiss_raw_pid_rand->Fill(pmiss, weight2);
				if (!contam_cut)
				{
					H1_pmiss_raw_pid_rand->Fill(pmiss, weight);
				}
				H1_th_rq_pid_rand->Fill(cth_rq / dtr, weight);
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_Q2_pid_acc->Fill(Q2, weight);
				H1_W_pid_acc->Fill(W, weight);
				H1_nu_pid_acc->Fill(nu, weight);
				H1_ph_q_pid_acc->Fill(ph_q / dtr, weight);
				H1_q_pid_acc->Fill(q, weight);
				H1_qx_pid_acc->Fill(q_x, weight);
				H1_qy_pid_acc->Fill(q_y, weight);
				H1_qz_pid_acc->Fill(q_z, weight);
				H1_th_e_pid_acc->Fill(th_e / dtr, weight);
				H1_th_q_pid_acc->Fill(th_q / dtr, weight);
				H1_xbj_pid_acc->Fill(x_bj, weight);
				H1_xbj_ratio_pid_acc->Fill(x_bj, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && (ctime_L || ctime_R) && HMS_TRK && SHMS_TRK)
			{
				H1_Q2_pid_acc_rand->Fill(Q2, weight);
				H1_W_pid_acc_rand->Fill(W, weight);
				H1_nu_pid_acc_rand->Fill(nu, weight);
				H1_q_pid_acc_rand->Fill(q, weight);
				H1_xbj_pid_acc_rand->Fill(x_bj, weight);
				H1_ep_ctime_pid_acc_rand->Fill(ep_ctime - ctime_offset, weight);

				H1_emiss_pid_acc_rand->Fill(emiss, weight);
				H1_pmiss_pid_acc_rand->Fill(pmiss, weight);
				H1_pmiss_raw_pid_acc_rand->Fill(pmiss, weight2);
				if (!contam_cut)
				{
					H1_pmiss_c_pid_acc_rand->Fill(pmiss, weight);
				}
				H1_th_rq_pid_acc_rand->Fill(cth_rq / dtr, weight);
			}
			
			if (user_runtype == "SRC")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_Q2_pid_acc_1->Fill(Q2, weight);
					H1_W_pid_acc_1->Fill(W, weight);
					H1_nu_pid_acc_1->Fill(nu, weight);
					H1_ph_q_pid_acc_1->Fill(ph_q / dtr, weight);
					H1_q_pid_acc_1->Fill(q, weight);
					H1_qx_pid_acc_1->Fill(q_x, weight);
					H1_qy_pid_acc_1->Fill(q_y, weight);
					H1_qz_pid_acc_1->Fill(q_z, weight);
					H1_th_e_pid_acc_1->Fill(th_e / dtr, weight);
					H1_th_q_pid_acc_1->Fill(th_q / dtr, weight);
					H1_xbj_pid_acc_1->Fill(x_bj, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_Q2_pid_acc_2->Fill(Q2, weight);
					H1_W_pid_acc_2->Fill(W, weight);
					H1_nu_pid_acc_2->Fill(nu, weight);
					H1_ph_q_pid_acc_2->Fill(ph_q / dtr, weight);
					H1_q_pid_acc_2->Fill(q, weight);
					H1_qx_pid_acc_2->Fill(q_x, weight);
					H1_qy_pid_acc_2->Fill(q_y, weight);
					H1_qz_pid_acc_2->Fill(q_z, weight);
					H1_th_e_pid_acc_2->Fill(th_e / dtr, weight);
					H1_th_q_pid_acc_2->Fill(th_q / dtr, weight);
					H1_xbj_pid_acc_2->Fill(x_bj, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_Q2_pid_acc_kin_alt->Fill(Q2, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_xbj_pid_acc_kin_alt->Fill(x_bj, weight);
				}
			}
			else if (user_runtype == "MF")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_Q2_pid_acc_1->Fill(Q2, weight);
					H1_W_pid_acc_1->Fill(W, weight);
					H1_nu_pid_acc_1->Fill(nu, weight);
					H1_ph_q_pid_acc_1->Fill(ph_q / dtr, weight);
					H1_q_pid_acc_1->Fill(q, weight);
					H1_qx_pid_acc_1->Fill(q_x, weight);
					H1_qy_pid_acc_1->Fill(q_y, weight);
					H1_qz_pid_acc_1->Fill(q_z, weight);
					H1_th_e_pid_acc_1->Fill(th_e / dtr, weight);
					H1_th_q_pid_acc_1->Fill(th_q / dtr, weight);
					H1_xbj_pid_acc_1->Fill(x_bj, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_Q2_pid_acc_2->Fill(Q2, weight);
					H1_W_pid_acc_2->Fill(W, weight);
					H1_nu_pid_acc_2->Fill(nu, weight);
					H1_ph_q_pid_acc_2->Fill(ph_q / dtr, weight);
					H1_q_pid_acc_2->Fill(q, weight);
					H1_qx_pid_acc_2->Fill(q_x, weight);
					H1_qy_pid_acc_2->Fill(q_y, weight);
					H1_qz_pid_acc_2->Fill(q_z, weight);
					H1_th_e_pid_acc_2->Fill(th_e / dtr, weight);
					H1_th_q_pid_acc_2->Fill(th_q / dtr, weight);
					H1_xbj_pid_acc_2->Fill(x_bj, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_Q2_pid_acc_kin_alt->Fill(Q2, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_xbj_pid_acc_kin_alt->Fill(x_bj, weight);
				}
			}
			else if (user_runtype == "Heep")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_Q2_pid_acc_1->Fill(Q2, weight);
					H1_W_pid_acc_1->Fill(W, weight);
					H1_nu_pid_acc_1->Fill(nu, weight);
					H1_ph_q_pid_acc_1->Fill(ph_q / dtr, weight);
					H1_q_pid_acc_1->Fill(q, weight);
					H1_qx_pid_acc_1->Fill(q_x, weight);
					H1_qy_pid_acc_1->Fill(q_y, weight);
					H1_qz_pid_acc_1->Fill(q_z, weight);
					H1_th_e_pid_acc_1->Fill(th_e / dtr, weight);
					H1_th_q_pid_acc_1->Fill(th_q / dtr, weight);
					H1_xbj_pid_acc_1->Fill(x_bj, weight);
				}
				
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1), Cuts::get_heep_xbj_cut(0), x_bj, Cuts::get_heep_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_Q2_pid_acc_2->Fill(Q2, weight);
					H1_W_pid_acc_2->Fill(W, weight);
					H1_nu_pid_acc_2->Fill(nu, weight);
					H1_ph_q_pid_acc_2->Fill(ph_q / dtr, weight);
					H1_q_pid_acc_2->Fill(q, weight);
					H1_qx_pid_acc_2->Fill(q_x, weight);
					H1_qy_pid_acc_2->Fill(q_y, weight);
					H1_qz_pid_acc_2->Fill(q_z, weight);
					H1_th_e_pid_acc_2->Fill(th_e / dtr, weight);
					H1_th_q_pid_acc_2->Fill(th_q / dtr, weight);
					H1_xbj_pid_acc_2->Fill(x_bj, weight);
				}
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_Q2_pid_acc_kin->Fill(Q2, weight);
				H1_W_pid_acc_kin->Fill(W, weight);
				H1_nu_pid_acc_kin->Fill(nu, weight);
				H1_ph_q_pid_acc_kin->Fill(ph_q / dtr, weight);
				H1_q_pid_acc_kin->Fill(q, weight);
				H1_qx_pid_acc_kin->Fill(q_x, weight);
				H1_qy_pid_acc_kin->Fill(q_y, weight);
				H1_qz_pid_acc_kin->Fill(q_z, weight);
				H1_th_e_pid_acc_kin->Fill(th_e / dtr, weight);
				H1_Pf_pid_acc_kin->Fill(Pf, weight);//verify
				H1_Ef_pid_acc_kin->Fill(Ef, weight);//verify
				H1_th_q_pid_acc_kin->Fill(th_q / dtr, weight);
				H1_xbj_pid_acc_kin->Fill(x_bj, weight);
				//H1_xbj_ratio_pid_acc_kin->Fill(x_bj, weight);
				
				
				//numLoops++;
			}
			
			etarx_corr = etarx - exptar * etarz * cos(8.3 * dtr);
			eXColl = etarx_corr + exptar * 253.;
			eYColl = etary + eyptar * 253. - (0.019 + 40. * .01 * 0.052) * edelta + (0.00019 + 40 * .01 * .00052) * edelta * edelta; //correct for HB horizontal bend

			if (0 <= edelta && edelta <= 22 && 1.5 <= Q2 && 0.85 <= pCalEtotTrkNorm && pCalEtotTrkNorm <= 1.3 && shms_Coll_gCut->IsInside(eYColl, eXColl) && ( (evtyp == 1) || (evtyp == 3) || (evtyp == 5) || (evtyp == 7) ) )
			{
				
				
				H1_edelta_pid_acc_kin->Fill(edelta, weight);
				//H1_Q2_pid_acc_kin->Fill(Q2, weight); 
				//H1_pCalEtotTrkNorm_pid_acc_kin->Fill(pCalEtotTrkNorm, weight); 
				H2_eXColl_eYColl_pid_acc_kin->Fill(eYColl, eXColl, weight); 
				H2_exptar_eyptar_pid_acc_kin->Fill(eyptar / dtr, exptar / dtr, weight);
				//H1_xbj_ratio_pid_acc_kin->Fill(x_bj, weight);
				
				
				//numLoops++;
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1)) && (ctime_L || ctime_R) && HMS_TRK && SHMS_TRK)
			{
				H1_Q2_ratio_pid_acc_kin_rand_full->Fill(Q2, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1)) && (ctime_L || ctime_R) && HMS_TRK && SHMS_TRK)
			{
				H1_xbj_ratio_pid_acc_kin_rand_full->Fill(x_bj, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1)) && (ctime_L || ctime_R) && HMS_TRK && SHMS_TRK)
			{
				H1_thrq_ratio_pid_acc_kin_rand_full->Fill(cth_rq / dtr, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1)) && (ctime_L || ctime_R) && HMS_TRK && SHMS_TRK)
			{
				H1_W_ratio_pid_acc_kin_rand_full->Fill(W, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1)) && (ctime_L || ctime_R) && HMS_TRK && SHMS_TRK)
			{
				H1_pmiss_ratio_pid_acc_kin_rand_full->Fill(pmiss, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && (ctime_L || ctime_R) && HMS_TRK && SHMS_TRK)
			{
				H1_Q2_pid_acc_kin_rand->Fill(Q2, weight);
				H1_W_pid_acc_kin_rand->Fill(W, weight);
				H1_nu_pid_acc_kin_rand->Fill(nu, weight);
				H1_q_pid_acc_kin_rand->Fill(q, weight);
				H1_xbj_pid_acc_kin_rand->Fill(x_bj, weight);
				H1_ep_ctime_pid_acc_kin_rand->Fill(ep_ctime - ctime_offset, weight);

				H1_emiss_pid_acc_kin_rand->Fill(emiss, weight);
				H1_pmiss_pid_acc_kin_rand->Fill(pmiss, weight);
				H1_pmiss_ratio_pid_acc_kin_rand->Fill(pmiss, weight);
				H1_pmiss_raw_pid_acc_kin_rand->Fill(pmiss, weight2);
				if (!contam_cut)
				{
					H1_pmiss_c_pid_acc_kin_rand->Fill(pmiss, weight);
				}
				H1_th_rq_pid_acc_kin_rand->Fill(cth_rq / dtr, weight);
			}

			//--------------------------------------
			//Secondary Hadron Kinematics
			//--------------------------------------
			H1_prec_x->Fill(prec_x, weight);
			H1_prec_y->Fill(prec_y, weight);
			H1_prec_z->Fill(prec_z, weight);
			H1_emiss->Fill(emiss, weight);
			//H1_emiss_nuc->Fill(emiss_nuc, weight);
			H1_ph_rq->Fill(ph_rq / dtr, weight);
			H1_ph_pq->Fill(ph_pq / dtr, weight);
			H1_pmiss->Fill(pmiss, weight);
			H1_pmiss_ratio->Fill(pmiss, weight);
			H1_pmiss_c->Fill(pmiss, weight2);
			H1_pmiss_raw->Fill(pmiss, weight2);
			H1_pmiss_x->Fill(pmiss_x, weight);
			H1_pmiss_y->Fill(pmiss_y, weight);
			H1_pmiss_z->Fill(pmiss_z, weight);
			H1_Tr->Fill(Tr, weight);
			H1_th_rq->Fill(cth_rq / dtr, weight);
			H1_cth_rq->Fill(cos(cth_rq), weight);
			H1_th_pq->Fill(th_pq / dtr, weight);
			H1_Tx->Fill(Tx, weight);
			H1_xangle->Fill(xangle / dtr, weight);
			H1_mmiss->Fill(mmiss, weight);
			H1_th_p->Fill(th_p / dtr, weight);
			//H1_omega->Fill(omega, weight);
			H1_Pm_par->Fill(Pm_par, weight);
			H1_Pm_per->Fill(Pm_per, weight);
			H1_alpha->Fill(alpha_v, weight);


			if (!contam_cut && pid_hms && pid_shms && Real_Golden && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_prec_x_pid->Fill(prec_x, weight);
				H1_prec_y_pid->Fill(prec_y, weight);
				H1_prec_z_pid->Fill(prec_z, weight);
				H1_emiss_pid->Fill(emiss, weight);
				//H1_emiss_nuc_pid->Fill(emiss_nuc, weight);
				H1_ph_rq_pid->Fill(ph_rq / dtr, weight);
				H1_ph_pq_pid->Fill(ph_pq / dtr, weight);
				H1_pmiss_pid->Fill(pmiss, weight);
				H1_pmiss_ratio_pid->Fill(pmiss, weight);
				H1_pmiss_raw_pid->Fill(pmiss, weight2);
				if (!contam_cut)
				{
					H1_pmiss_c_pid->Fill(pmiss, weight);
				}
				H1_pmiss_x_pid->Fill(pmiss_x, weight);
				H1_pmiss_y_pid->Fill(pmiss_y, weight);
				H1_pmiss_z_pid->Fill(pmiss_z, weight);
				H1_Tr_pid->Fill(Tr, weight);
				H1_th_rq_pid->Fill(cth_rq / dtr, weight);
				H1_cth_rq_pid->Fill(cos(cth_rq), weight);
				H1_th_pq_pid->Fill(th_pq / dtr, weight);
				H1_Tx_pid->Fill(Tx, weight);
				H1_xangle_pid->Fill(xangle / dtr, weight);
				H1_mmiss_pid->Fill(mmiss, weight);
				H1_th_p_pid->Fill(th_p / dtr, weight);
				//H1_omega_pid->Fill(omega, weight);
				H1_Pm_par_pid->Fill(Pm_par, weight);
				H1_Pm_per_pid->Fill(Pm_per, weight);
				H1_alpha_pid->Fill(alpha_v, weight);
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_prec_x_pid_acc->Fill(prec_x, weight);
				H1_prec_y_pid_acc->Fill(prec_y, weight);
				H1_prec_z_pid_acc->Fill(prec_z, weight);
				H1_emiss_pid_acc->Fill(emiss, weight);
				//H1_emiss_nuc_pid_acc->Fill(emiss_nuc, weight);
				H1_ph_rq_pid_acc->Fill(ph_rq / dtr, weight);
				H1_ph_pq_pid_acc->Fill(ph_pq / dtr, weight);
				H1_pmiss_pid_acc->Fill(pmiss, weight);
				H1_pmiss_ratio_pid_acc->Fill(pmiss, weight);
				H1_pmiss_raw_pid_acc->Fill(pmiss, weight2);
				if (!contam_cut)
				{
					H1_pmiss_c_pid_acc->Fill(pmiss, weight);
				}
				H1_pmiss_x_pid_acc->Fill(pmiss_x, weight);
				H1_pmiss_y_pid_acc->Fill(pmiss_y, weight);
				H1_pmiss_z_pid_acc->Fill(pmiss_z, weight);
				H1_Tr_pid_acc->Fill(Tr, weight);
				H1_th_rq_pid_acc->Fill(cth_rq / dtr, weight);
				H1_cth_rq_pid_acc->Fill(cos(cth_rq), weight);
				H1_th_pq_pid_acc->Fill(th_pq / dtr, weight);
				H1_Tx_pid_acc->Fill(Tx, weight);
				H1_xangle_pid_acc->Fill(xangle / dtr, weight);
				H1_mmiss_pid_acc->Fill(mmiss, weight);
				H1_th_p_pid_acc->Fill(th_p / dtr, weight);
				//H1_omega_pid_acc->Fill(omega, weight);
				H1_Pm_par_pid_acc->Fill(Pm_par, weight);
				H1_Pm_per_pid_acc->Fill(Pm_per, weight);
				H1_alpha_pid_acc->Fill(alpha_v, weight);
			}

			
			if (user_runtype == "SRC")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_prec_x_pid_acc_1->Fill(prec_x, weight);
					H1_prec_y_pid_acc_1->Fill(prec_y, weight);
					H1_prec_z_pid_acc_1->Fill(prec_z, weight);
					H1_emiss_pid_acc_1->Fill(emiss, weight);
					//H1_emiss_nuc_pid_acc_1->Fill(emiss_nuc, weight);
					H1_ph_rq_pid_acc_1->Fill(ph_rq / dtr, weight);
					H1_ph_pq_pid_acc_1->Fill(ph_pq / dtr, weight);
					H1_pmiss_pid_acc_1->Fill(pmiss, weight);
					H1_pmiss_raw_pid_acc_1->Fill(pmiss, weight2);
					if (!contam_cut)
					{
						H1_pmiss_c_pid_acc_1->Fill(pmiss, weight);
					}
					H1_pmiss_x_pid_acc_1->Fill(pmiss_x, weight);
					H1_pmiss_y_pid_acc_1->Fill(pmiss_y, weight);
					H1_pmiss_z_pid_acc_1->Fill(pmiss_z, weight);
					H1_Tr_pid_acc_1->Fill(Tr, weight);
					H1_th_rq_pid_acc_1->Fill(cth_rq / dtr, weight);
					H1_cth_rq_pid_acc_1->Fill(cos(cth_rq), weight);
					H1_th_pq_pid_acc_1->Fill(th_pq / dtr, weight);
					H1_Tx_pid_acc_1->Fill(Tx, weight);
					H1_xangle_pid_acc_1->Fill(xangle / dtr, weight);
					H1_mmiss_pid_acc_1->Fill(mmiss, weight);
					H1_th_p_pid_acc_1->Fill(th_p / dtr, weight);
					//H1_omega_pid_acc_1->Fill(omega, weight);
					H1_Pm_par_pid_acc_1->Fill(Pm_par, weight);
					H1_Pm_per_pid_acc_1->Fill(Pm_per, weight);
					H1_alpha_pid_acc_1->Fill(alpha_v, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_prec_x_pid_acc_2->Fill(prec_x, weight);
					H1_prec_y_pid_acc_2->Fill(prec_y, weight);
					H1_prec_z_pid_acc_2->Fill(prec_z, weight);
					H1_emiss_pid_acc_2->Fill(emiss, weight);
					//H1_emiss_nuc_pid_acc_2->Fill(emiss_nuc, weight);
					H1_ph_rq_pid_acc_2->Fill(ph_rq / dtr, weight);
					H1_ph_pq_pid_acc_2->Fill(ph_pq / dtr, weight);
					H1_pmiss_pid_acc_2->Fill(pmiss, weight);
					H1_pmiss_raw_pid_acc_2->Fill(pmiss, weight2);
					if (!contam_cut)
					{
						H1_pmiss_c_pid_acc_2->Fill(pmiss, weight);
					}
					H1_pmiss_x_pid_acc_2->Fill(pmiss_x, weight);
					H1_pmiss_y_pid_acc_2->Fill(pmiss_y, weight);
					H1_pmiss_z_pid_acc_2->Fill(pmiss_z, weight);
					H1_Tr_pid_acc_2->Fill(Tr, weight);
					H1_th_rq_pid_acc_2->Fill(cth_rq / dtr, weight);
					H1_cth_rq_pid_acc_2->Fill(cos(cth_rq), weight);
					H1_th_pq_pid_acc_2->Fill(th_pq / dtr, weight);
					H1_Tx_pid_acc_2->Fill(Tx, weight);
					H1_xangle_pid_acc_2->Fill(xangle / dtr, weight);
					H1_mmiss_pid_acc_2->Fill(mmiss, weight);
					H1_th_p_pid_acc_2->Fill(th_p / dtr, weight);
					//H1_omega_pid_acc_2->Fill(omega, weight);
					H1_Pm_par_pid_acc_2->Fill(Pm_par, weight);
					H1_Pm_per_pid_acc_2->Fill(Pm_per, weight);
					H1_alpha_pid_acc_2->Fill(alpha_v, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_pmiss_pid_acc_kin_alt->Fill(pmiss, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_emiss_pid_acc_kin_alt->Fill(emiss, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_th_rq_pid_acc_kin_alt->Fill(cth_rq / dtr, weight);
				}
			}
			else if (user_runtype == "MF")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_prec_x_pid_acc_1->Fill(prec_x, weight);
					H1_prec_y_pid_acc_1->Fill(prec_y, weight);
					H1_prec_z_pid_acc_1->Fill(prec_z, weight);
					H1_emiss_pid_acc_1->Fill(emiss, weight);
					//H1_emiss_nuc_pid_acc_1->Fill(emiss_nuc, weight);
					H1_ph_rq_pid_acc_1->Fill(ph_rq / dtr, weight);
					H1_ph_pq_pid_acc_1->Fill(ph_pq / dtr, weight);
					H1_pmiss_pid_acc_1->Fill(pmiss, weight);
					H1_pmiss_raw_pid_acc_1->Fill(pmiss, weight2);
					if (!contam_cut)
					{
						H1_pmiss_c_pid_acc_1->Fill(pmiss, weight);
					}
					H1_pmiss_x_pid_acc_1->Fill(pmiss_x, weight);
					H1_pmiss_y_pid_acc_1->Fill(pmiss_y, weight);
					H1_pmiss_z_pid_acc_1->Fill(pmiss_z, weight);
					H1_Tr_pid_acc_1->Fill(Tr, weight);
					H1_th_rq_pid_acc_1->Fill(cth_rq / dtr, weight);
					H1_cth_rq_pid_acc_1->Fill(cos(cth_rq), weight);
					H1_th_pq_pid_acc_1->Fill(th_pq / dtr, weight);
					H1_Tx_pid_acc_1->Fill(Tx, weight);
					H1_xangle_pid_acc_1->Fill(xangle / dtr, weight);
					H1_mmiss_pid_acc_1->Fill(mmiss, weight);
					H1_th_p_pid_acc_1->Fill(th_p / dtr, weight);
					//H1_omega_pid_acc_1->Fill(omega, weight);
					H1_Pm_par_pid_acc_1->Fill(Pm_par, weight);
					H1_Pm_per_pid_acc_1->Fill(Pm_per, weight);
					H1_alpha_pid_acc_1->Fill(alpha_v, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_prec_x_pid_acc_2->Fill(prec_x, weight);
					H1_prec_y_pid_acc_2->Fill(prec_y, weight);
					H1_prec_z_pid_acc_2->Fill(prec_z, weight);
					H1_emiss_pid_acc_2->Fill(emiss, weight);
					//H1_emiss_nuc_pid_acc_2->Fill(emiss_nuc, weight);
					H1_ph_rq_pid_acc_2->Fill(ph_rq / dtr, weight);
					H1_ph_pq_pid_acc_2->Fill(ph_pq / dtr, weight);
					H1_pmiss_pid_acc_2->Fill(pmiss, weight);
					H1_pmiss_raw_pid_acc_2->Fill(pmiss, weight2);
					if (!contam_cut)
					{
						H1_pmiss_c_pid_acc_2->Fill(pmiss, weight);
					}
					H1_pmiss_x_pid_acc_2->Fill(pmiss_x, weight);
					H1_pmiss_y_pid_acc_2->Fill(pmiss_y, weight);
					H1_pmiss_z_pid_acc_2->Fill(pmiss_z, weight);
					H1_Tr_pid_acc_2->Fill(Tr, weight);
					H1_th_rq_pid_acc_2->Fill(cth_rq / dtr, weight);
					H1_cth_rq_pid_acc_2->Fill(cos(cth_rq), weight);
					H1_th_pq_pid_acc_2->Fill(th_pq / dtr, weight);
					H1_Tx_pid_acc_2->Fill(Tx, weight);
					H1_xangle_pid_acc_2->Fill(xangle / dtr, weight);
					H1_mmiss_pid_acc_2->Fill(mmiss, weight);
					H1_th_p_pid_acc_2->Fill(th_p / dtr, weight);
					//H1_omega_pid_acc_2->Fill(omega, weight);
					H1_Pm_par_pid_acc_2->Fill(Pm_par, weight);
					H1_Pm_per_pid_acc_2->Fill(Pm_per, weight);
					H1_alpha_pid_acc_2->Fill(alpha_v, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_pmiss_pid_acc_kin_alt->Fill(pmiss, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_emiss_pid_acc_kin_alt->Fill(emiss, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_th_rq_pid_acc_kin_alt->Fill(cth_rq / dtr, weight);
				}
			}
			else if (user_runtype == "Heep")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_prec_x_pid_acc_1->Fill(prec_x, weight);
					H1_prec_y_pid_acc_1->Fill(prec_y, weight);
					H1_prec_z_pid_acc_1->Fill(prec_z, weight);
					H1_emiss_pid_acc_1->Fill(emiss, weight);
					//H1_emiss_nuc_pid_acc_1->Fill(emiss_nuc, weight);
					H1_ph_rq_pid_acc_1->Fill(ph_rq / dtr, weight);
					H1_ph_pq_pid_acc_1->Fill(ph_pq / dtr, weight);
					H1_pmiss_pid_acc_1->Fill(pmiss, weight);
					H1_pmiss_raw_pid_acc_1->Fill(pmiss, weight2);
					if (!contam_cut)
					{
						H1_pmiss_c_pid_acc_1->Fill(pmiss, weight);
					}
					H1_pmiss_x_pid_acc_1->Fill(pmiss_x, weight);
					H1_pmiss_y_pid_acc_1->Fill(pmiss_y, weight);
					H1_pmiss_z_pid_acc_1->Fill(pmiss_z, weight);
					H1_Tr_pid_acc_1->Fill(Tr, weight);
					H1_th_rq_pid_acc_1->Fill(cth_rq / dtr, weight);
					H1_cth_rq_pid_acc_1->Fill(cos(cth_rq), weight);
					H1_th_pq_pid_acc_1->Fill(th_pq / dtr, weight);
					H1_Tx_pid_acc_1->Fill(Tx, weight);
					H1_xangle_pid_acc_1->Fill(xangle / dtr, weight);
					H1_mmiss_pid_acc_1->Fill(mmiss, weight);
					H1_th_p_pid_acc_1->Fill(th_p / dtr, weight);
					//H1_omega_pid_acc_1->Fill(omega, weight);
					H1_Pm_par_pid_acc_1->Fill(Pm_par, weight);
					H1_Pm_per_pid_acc_1->Fill(Pm_per, weight);
					H1_alpha_pid_acc_1->Fill(alpha_v, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1), Cuts::get_heep_xbj_cut(0), x_bj, Cuts::get_heep_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H1_prec_x_pid_acc_2->Fill(prec_x, weight);
					H1_prec_y_pid_acc_2->Fill(prec_y, weight);
					H1_prec_z_pid_acc_2->Fill(prec_z, weight);
					H1_emiss_pid_acc_2->Fill(emiss, weight);
					//H1_emiss_nuc_pid_acc_2->Fill(emiss_nuc, weight);
					H1_ph_rq_pid_acc_2->Fill(ph_rq / dtr, weight);
					H1_ph_pq_pid_acc_2->Fill(ph_pq / dtr, weight);
					H1_pmiss_pid_acc_2->Fill(pmiss, weight);
					H1_pmiss_raw_pid_acc_2->Fill(pmiss, weight2);
					if (!contam_cut)
					{
						H1_pmiss_c_pid_acc_2->Fill(pmiss, weight);
					}
					H1_pmiss_x_pid_acc_2->Fill(pmiss_x, weight);
					H1_pmiss_y_pid_acc_2->Fill(pmiss_y, weight);
					H1_pmiss_z_pid_acc_2->Fill(pmiss_z, weight);
					H1_Tr_pid_acc_2->Fill(Tr, weight);
					H1_th_rq_pid_acc_2->Fill(cth_rq / dtr, weight);
					H1_cth_rq_pid_acc_2->Fill(cos(cth_rq), weight);
					H1_th_pq_pid_acc_2->Fill(th_pq / dtr, weight);
					H1_Tx_pid_acc_2->Fill(Tx, weight);
					H1_xangle_pid_acc_2->Fill(xangle / dtr, weight);
					H1_mmiss_pid_acc_2->Fill(mmiss, weight);
					H1_th_p_pid_acc_2->Fill(th_p / dtr, weight);
					//H1_omega_pid_acc_2->Fill(omega, weight);
					H1_Pm_par_pid_acc_2->Fill(Pm_par, weight);
					H1_Pm_per_pid_acc_2->Fill(Pm_per, weight);
					H1_alpha_pid_acc_2->Fill(alpha_v, weight);
				}
			}


			/*
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_hxfp_hyfp_pid_acc_1->Fill(hyfp, hxfp, weight);
					H2_exfp_eyfp_pid_acc_1->Fill(eyfp, exfp, weight);

					H2_exptar_eyptar_pid_acc_1->Fill(eyptar / dtr, exptar / dtr, weight);
					H2_hxptar_hyptar_pid_acc_1->Fill(hyptar / dtr, hxptar / dtr, weight);

					H2_hxptar_exptar_pid_acc_1->Fill(exptar / dtr, hxptar / dtr, weight);
					H2_hyptar_eyptar_pid_acc_1->Fill(eyptar / dtr, hyptar / dtr, weight);
					H2_hdelta_edelta_pid_acc_1->Fill(edelta, hdelta, weight);
					H2_hXColl_hYColl_pid_acc_1->Fill(hYColl, hXColl, weight);
					H2_eXColl_eYColl_pid_acc_1->Fill(eYColl, eXColl, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
			*/

			/*
			ctime_L = false;// Cuts::cTime_Rand_L(ep_ctime - ctime_offset);
			ctime_R = Cuts::cTime_Rand_R(ep_ctime - ctime_offset);
			HMS_TRK = true;// Cuts::HMS_Tr(hdc_ntrk, hScinGood, hCerNpeSum, hCalEtotNorm, hHodBetaNtrk);
			//problematic!!!!!!!!!!!
			SHMS_TRK = true;// Cuts::SHMS_Tr(pdc_ntrk, pScinGood, pNGCerNpeSum, pHGCerNpeSum, pCalEtotNorm, pHodBetaNtrk);
			ctime = Cuts::cTime(ep_ctime - ctime_offset);
			pid_hms = Cuts::PID_HMS(hCalEtotTrkNorm, hCerNpeSum);
			pid_shms = Cuts::PID_SHMS(pCalEtotTrkNorm, pNGCerNpeSum, pHGCerNpeSum, 1);
			accp_hms = Cuts::Accp_HMS(hdelta, hxptar, hyptar, hms_Coll_gCut->IsInside(hYColl, hXColl));
			accp_shms = Cuts::Accp_SHMS(edelta, exptar, eyptar, shms_Coll_gCut->IsInside(eYColl, eXColl));
			contam_cut = contam_gCut->IsInside(emiss, pmiss);
			*/

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_Q2_ratio_pid_acc_kin_full->Fill(Q2, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_xbj_ratio_pid_acc_kin_full->Fill(x_bj, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_thrq_ratio_pid_acc_kin_full->Fill((cth_rq / dtr), weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_W_ratio_pid_acc_kin_full->Fill(W, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_pmiss_ratio_pid_acc_kin_full->Fill(pmiss, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_Em_ratio_pid_acc_kin_full->Fill(emiss, weight);
			}
			if (!contam_cut && pid_hms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_pCalEtotTrkNorm_ratio_pid_acc_kin_full->Fill(pCalEtotTrkNorm, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && shms_Coll_gCut->IsInside(eYColl, eXColl) && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_edelta_ratio_pid_acc_kin_full->Fill(edelta, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && hms_Coll_gCut->IsInside(hYColl, hXColl) && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_hdelta_ratio_pid_acc_kin_full->Fill(hdelta, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && 0 < edelta && edelta < 22 && accp_hms && accp_z && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_eXColl_eYColl_ratio_pid_acc_kin_full->Fill(eYColl, eXColl, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_shms && -10 < hdelta && hdelta < 10 && accp_z && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_hXColl_hYColl_ratio_pid_acc_kin_full->Fill(hYColl, hXColl, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_Em_Pm_ratio_pid_acc_kin_full->Fill(pmiss, emiss, weight);
			}
			if (pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_Q2_cut(0), Q2, Cuts::get_SRC_Q2_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_Em_Pm_ratio_pid_acc_kin_full2->Fill(pmiss, emiss, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_Pm_cut(0), pmiss, Cuts::get_SRC_Pm_cut(1)) && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_xbj_Q2_ratio_pid_acc_kin_full->Fill(Q2, x_bj, weight);
			}
			if (!contam_cut && pid_hms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_pCalEtotTrkNorm_ratio_pid_acc_kin_full_mf->Fill(pCalEtotTrkNorm, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && shms_Coll_gCut->IsInside(eYColl, eXColl) && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_edelta_ratio_pid_acc_kin_full_mf->Fill(edelta, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && hms_Coll_gCut->IsInside(hYColl, hXColl) && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_hdelta_ratio_pid_acc_kin_full_mf->Fill(hdelta, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && 0 < edelta && edelta < 22 && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_eXColl_eYColl_ratio_pid_acc_kin_full_mf->Fill(eYColl, eXColl, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && -10 < hdelta && hdelta < 10 && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_hXColl_hYColl_ratio_pid_acc_kin_full_mf->Fill(hYColl, hXColl, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_Q2_ratio_pid_acc_kin_full_mf->Fill(Q2, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_pmiss_ratio_pid_acc_kin_full_mf->Fill(pmiss, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_W_ratio_pid_acc_kin_full_mf->Fill(W, weight);
				H1_Pf_pid_acc_kin_mf->Fill(Pf, weight);
				H1_Ef_pid_acc_kin_mf->Fill(Ef, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1), Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_Em_ratio_pid_acc_kin_full_mf->Fill(emiss, weight);
			}
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_Em_Pm_ratio_pid_acc_kin_full_mf->Fill(pmiss, emiss, weight);
			}
			if (pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Q2_cut(0), Q2, Cuts::get_MF_Q2_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_Em_Pm_ratio_pid_acc_kin_full_mf2->Fill(pmiss, emiss, weight);
			}


			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H1_prec_x_pid_acc_kin->Fill(prec_x, weight);
				H1_prec_y_pid_acc_kin->Fill(prec_y, weight);
				H1_prec_z_pid_acc_kin->Fill(prec_z, weight);
				H1_emiss_pid_acc_kin->Fill(emiss, weight);
				//H1_emiss_nuc_pid_acc_kin->Fill(emiss_nuc, weight);
				H1_ph_rq_pid_acc_kin->Fill(ph_rq / dtr, weight);
				H1_ph_pq_pid_acc_kin->Fill(ph_pq / dtr, weight);
				H1_pmiss_pid_acc_kin->Fill(pmiss, weight);

				


				H1_pmiss_ratio_pid_acc_kin->Fill(pmiss, weight);

				H1_pmiss_raw_pid_acc_kin->Fill(pmiss, weight2);
				if (!contam_cut)
				{
					H1_pmiss_c_pid_acc_kin->Fill(pmiss, weight);
				}
				H1_pmiss_x_pid_acc_kin->Fill(pmiss_x, weight);
				H1_pmiss_y_pid_acc_kin->Fill(pmiss_y, weight);
				H1_pmiss_z_pid_acc_kin->Fill(pmiss_z, weight);
				H1_Tr_pid_acc_kin->Fill(Tr, weight);
				H1_th_rq_pid_acc_kin->Fill(cth_rq / dtr, weight);
				H1_cth_rq_pid_acc_kin->Fill(cos(cth_rq), weight);
				H1_th_pq_pid_acc_kin->Fill(th_pq / dtr, weight);
				H1_Tx_pid_acc_kin->Fill(Tx, weight);
				H1_xangle_pid_acc_kin->Fill(xangle / dtr, weight);
				H1_mmiss_pid_acc_kin->Fill(mmiss, weight);
				H1_th_p_pid_acc_kin->Fill(th_p / dtr, weight);
				//H1_omega_pid_acc_kin->Fill(omega, weight);
				H1_Pm_par_pid_acc_kin->Fill(Pm_par, weight);
				H1_Pm_per_pid_acc_kin->Fill(Pm_per, weight);
				H1_alpha_pid_acc_kin->Fill(alpha_v, weight);
			}


			//--------------------------------------
			//2D Histograms
			//--------------------------------------
			H2_hxfp_hyfp->Fill(hyfp, hxfp, weight);
			H2_exfp_eyfp->Fill(eyfp, exfp, weight);
		
			H2_exptar_eyptar->Fill(eyptar / dtr, exptar / dtr, weight);
			H2_hxptar_hyptar->Fill(hyptar / dtr, hxptar / dtr, weight);

			H2_hxptar_exptar->Fill(exptar / dtr, hxptar / dtr, weight);
			H2_hyptar_eyptar->Fill(eyptar / dtr, hyptar / dtr, weight);
			H2_hdelta_edelta->Fill(edelta, hdelta, weight);
			H2_pCalEtotTrkNorm_edelta->Fill(edelta, pCalEtotTrkNorm, weight);

			H2_hXColl_hYColl->Fill(hYColl, hXColl, weight);
			H2_eXColl_eYColl->Fill(eYColl, eXColl, weight);

			
			if (!contam_cut && pid_hms && pid_shms && Real_Golden && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_hxfp_hyfp_pid->Fill(hyfp, hxfp, weight);
				H2_exfp_eyfp_pid->Fill(eyfp, exfp, weight);

				H2_exptar_eyptar_pid->Fill(eyptar / dtr, exptar / dtr, weight);
				H2_hxptar_hyptar_pid->Fill(hyptar / dtr, hxptar / dtr, weight);

				H2_hxptar_exptar_pid->Fill(exptar / dtr, hxptar / dtr, weight);
				H2_hyptar_eyptar_pid->Fill(eyptar / dtr, hyptar / dtr, weight);
				H2_hdelta_edelta_pid->Fill(edelta, hdelta, weight);
				H2_pCalEtotTrkNorm_edelta_pid->Fill(edelta, pCalEtotTrkNorm, weight);
				H2_hXColl_hYColl_pid->Fill(hYColl, hXColl, weight);
				H2_eXColl_eYColl_pid->Fill(eYColl, eXColl, weight);
			}

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_hxfp_hyfp_pid_acc->Fill(hyfp, hxfp, weight);
				H2_exfp_eyfp_pid_acc->Fill(eyfp, exfp, weight);

				H2_exptar_eyptar_pid_acc->Fill(eyptar / dtr, exptar / dtr, weight);
				H2_hxptar_hyptar_pid_acc->Fill(hyptar / dtr, hxptar / dtr, weight);

				H2_hxptar_exptar_pid_acc->Fill(exptar / dtr, hxptar / dtr, weight);
				H2_hyptar_eyptar_pid_acc->Fill(eyptar / dtr, hyptar / dtr, weight);
				H2_hdelta_edelta_pid_acc->Fill(edelta, hdelta, weight);
				H2_pCalEtotTrkNorm_edelta_pid_acc->Fill(edelta, pCalEtotTrkNorm, weight);
				H2_hXColl_hYColl_pid_acc->Fill(hYColl, hXColl, weight);
				H2_eXColl_eYColl_pid_acc->Fill(eYColl, eXColl, weight);
			}

			
			if (user_runtype == "SRC")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_hxfp_hyfp_pid_acc_1->Fill(hyfp, hxfp, weight);
					H2_exfp_eyfp_pid_acc_1->Fill(eyfp, exfp, weight);

					H2_exptar_eyptar_pid_acc_1->Fill(eyptar / dtr, exptar / dtr, weight);
					H2_hxptar_hyptar_pid_acc_1->Fill(hyptar / dtr, hxptar / dtr, weight);

					H2_hxptar_exptar_pid_acc_1->Fill(exptar / dtr, hxptar / dtr, weight);
					H2_hyptar_eyptar_pid_acc_1->Fill(eyptar / dtr, hyptar / dtr, weight);
					H2_hdelta_edelta_pid_acc_1->Fill(edelta, hdelta, weight);
					H2_pCalEtotTrkNorm_edelta_pid_acc_1->Fill(edelta, pCalEtotTrkNorm, weight);
					H2_hXColl_hYColl_pid_acc_1->Fill(hYColl, hXColl, weight);
					H2_eXColl_eYColl_pid_acc_1->Fill(eYColl, eXColl, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_hxfp_hyfp_pid_acc_2->Fill(hyfp, hxfp, weight);
					H2_exfp_eyfp_pid_acc_2->Fill(eyfp, exfp, weight);

					H2_exptar_eyptar_pid_acc_2->Fill(eyptar / dtr, exptar / dtr, weight);
					H2_hxptar_hyptar_pid_acc_2->Fill(hyptar / dtr, hxptar / dtr, weight);

					H2_hxptar_exptar_pid_acc_2->Fill(exptar / dtr, hxptar / dtr, weight);
					H2_hyptar_eyptar_pid_acc_2->Fill(eyptar / dtr, hyptar / dtr, weight);
					H2_hdelta_edelta_pid_acc_2->Fill(edelta, hdelta, weight);
					H2_pCalEtotTrkNorm_edelta_pid_acc_2->Fill(edelta, pCalEtotTrkNorm, weight);
					H2_hXColl_hYColl_pid_acc_2->Fill(hYColl, hXColl, weight);
					H2_eXColl_eYColl_pid_acc_2->Fill(eYColl, eXColl, weight);
				}
			}
			else if (user_runtype == "MF")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_hxfp_hyfp_pid_acc_1->Fill(hyfp, hxfp, weight);
					H2_exfp_eyfp_pid_acc_1->Fill(eyfp, exfp, weight);

					H2_exptar_eyptar_pid_acc_1->Fill(eyptar / dtr, exptar / dtr, weight);
					H2_hxptar_hyptar_pid_acc_1->Fill(hyptar / dtr, hxptar / dtr, weight);

					H2_hxptar_exptar_pid_acc_1->Fill(exptar / dtr, hxptar / dtr, weight);
					H2_hyptar_eyptar_pid_acc_1->Fill(eyptar / dtr, hyptar / dtr, weight);
					H2_hdelta_edelta_pid_acc_1->Fill(edelta, hdelta, weight);
					H2_pCalEtotTrkNorm_edelta_pid_acc_1->Fill(edelta, pCalEtotTrkNorm, weight);
					H2_hXColl_hYColl_pid_acc_1->Fill(hYColl, hXColl, weight);
					H2_eXColl_eYColl_pid_acc_1->Fill(eYColl, eXColl, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_hxfp_hyfp_pid_acc_2->Fill(hyfp, hxfp, weight);
					H2_exfp_eyfp_pid_acc_2->Fill(eyfp, exfp, weight);

					H2_exptar_eyptar_pid_acc_2->Fill(eyptar / dtr, exptar / dtr, weight);
					H2_hxptar_hyptar_pid_acc_2->Fill(hyptar / dtr, hxptar / dtr, weight);

					H2_hxptar_exptar_pid_acc_2->Fill(exptar / dtr, hxptar / dtr, weight);
					H2_hyptar_eyptar_pid_acc_2->Fill(eyptar / dtr, hyptar / dtr, weight);
					H2_hdelta_edelta_pid_acc_2->Fill(edelta, hdelta, weight);
					H2_pCalEtotTrkNorm_edelta_pid_acc_2->Fill(edelta, pCalEtotTrkNorm, weight);
					H2_hXColl_hYColl_pid_acc_2->Fill(hYColl, hXColl, weight);
					H2_eXColl_eYColl_pid_acc_2->Fill(eYColl, eXColl, weight);
				}
			}
			else if (user_runtype == "Heep")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_hxfp_hyfp_pid_acc_1->Fill(hyfp, hxfp, weight);
					H2_exfp_eyfp_pid_acc_1->Fill(eyfp, exfp, weight);

					H2_exptar_eyptar_pid_acc_1->Fill(eyptar / dtr, exptar / dtr, weight);
					H2_hxptar_hyptar_pid_acc_1->Fill(hyptar / dtr, hxptar / dtr, weight);

					H2_hxptar_exptar_pid_acc_1->Fill(exptar / dtr, hxptar / dtr, weight);
					H2_hyptar_eyptar_pid_acc_1->Fill(eyptar / dtr, hyptar / dtr, weight);
					H2_hdelta_edelta_pid_acc_1->Fill(edelta, hdelta, weight);
					H2_pCalEtotTrkNorm_edelta_pid_acc_1->Fill(edelta, pCalEtotTrkNorm, weight);
					H2_hXColl_hYColl_pid_acc_1->Fill(hYColl, hXColl, weight);
					H2_eXColl_eYColl_pid_acc_1->Fill(eYColl, eXColl, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1), Cuts::get_heep_xbj_cut(0), x_bj, Cuts::get_heep_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_hxfp_hyfp_pid_acc_2->Fill(hyfp, hxfp, weight);
					H2_exfp_eyfp_pid_acc_2->Fill(eyfp, exfp, weight);

					H2_exptar_eyptar_pid_acc_2->Fill(eyptar / dtr, exptar / dtr, weight);
					H2_hxptar_hyptar_pid_acc_2->Fill(hyptar / dtr, hxptar / dtr, weight);

					H2_hxptar_exptar_pid_acc_2->Fill(exptar / dtr, hxptar / dtr, weight);
					H2_hyptar_eyptar_pid_acc_2->Fill(eyptar / dtr, hyptar / dtr, weight);
					H2_hdelta_edelta_pid_acc_2->Fill(edelta, hdelta, weight);
					H2_pCalEtotTrkNorm_edelta_pid_acc_2->Fill(edelta, pCalEtotTrkNorm, weight);
					H2_hXColl_hYColl_pid_acc_2->Fill(hYColl, hXColl, weight);
					H2_eXColl_eYColl_pid_acc_2->Fill(eYColl, eXColl, weight);
				}
			}


			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_hxfp_hyfp_pid_acc_kin->Fill(hyfp, hxfp, weight);
				H2_exfp_eyfp_pid_acc_kin->Fill(eyfp, exfp, weight);

				//H2_exptar_eyptar_pid_acc_kin->Fill(eyptar / dtr, exptar / dtr, weight);
				H2_hxptar_hyptar_pid_acc_kin->Fill(hyptar / dtr, hxptar / dtr, weight);

				H2_hxptar_exptar_pid_acc_kin->Fill(exptar / dtr, hxptar / dtr, weight);
				H2_hyptar_eyptar_pid_acc_kin->Fill(eyptar / dtr, hyptar / dtr, weight);
				H2_hdelta_edelta_pid_acc_kin->Fill(edelta, hdelta, weight);
				H2_pCalEtotTrkNorm_edelta_pid_acc_kin->Fill(edelta, pCalEtotTrkNorm, weight);
				H2_hXColl_hYColl_pid_acc_kin->Fill(hYColl, hXColl, weight);
				//H2_eXColl_eYColl_pid_acc_kin->Fill(eYColl, eXColl, weight);
			}
			







			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_cthrq_Pm_pid_acc->Fill(pmiss, cos(cth_rq), weight);
				H2_thrq_Pm_pid_acc->Fill(pmiss, cth_rq / dtr, weight);
				H2_thrq_Em_pid_acc->Fill(emiss, cth_rq / dtr, weight);
				H2_thrq_Q2_pid_acc->Fill(Q2, cth_rq / dtr, weight);
				H2_thrq_xbj_pid_acc->Fill(x_bj, cth_rq / dtr, weight);
				H2_xbj_Q2_pid_acc->Fill(Q2, x_bj, weight);
				H2_xbj_Em_pid_acc->Fill(emiss, x_bj, weight);
				H2_xbj_Pm_pid_acc->Fill(pmiss, x_bj, weight);
				H2_Em_Pm_pid_acc->Fill(pmiss, emiss, weight);
				
				if (!contam_cut)
				{
					H2_Em_Pm_c_pid_acc->Fill(pmiss, emiss, weight);
				}
				
				H2_Em_Q2_pid_acc->Fill(Q2, emiss, weight);
				H2_Q2_Pm_pid_acc->Fill(pmiss, Q2, weight);
				H2_W_Pm_pid_acc->Fill(pmiss, W, weight);
				H2_W_thrq_pid_acc->Fill(cth_rq / dtr, W, weight);
				H2_W_Em_pid_acc->Fill(emiss, W, weight);
				H2_W_Q2_pid_acc->Fill(Q2, W, weight);
				H2_W_xbj_pid_acc->Fill(x_bj, W, weight); 
				H2_exptar_Pm_pid_acc->Fill(pmiss, exptar / dtr, weight);
				H2_eyptar_Pm_pid_acc->Fill(pmiss, eyptar / dtr, weight);
				H2_hxptar_Pm_pid_acc->Fill(pmiss, hxptar / dtr, weight);
				H2_hyptar_Pm_pid_acc->Fill(pmiss, hyptar / dtr, weight);
				H2_alpha_Pm_per_pid_acc->Fill(Pm_per, alpha_v, weight);
			}


			if (user_runtype == "SRC")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_cthrq_Pm_pid_acc_1->Fill(pmiss, cos(cth_rq), weight);
					H2_thrq_Pm_pid_acc_1->Fill(pmiss, cth_rq / dtr, weight);
					H2_thrq_Em_pid_acc_1->Fill(emiss, cth_rq / dtr, weight);
					H2_thrq_Q2_pid_acc_1->Fill(Q2, cth_rq / dtr, weight);
					H2_thrq_xbj_pid_acc_1->Fill(x_bj, cth_rq / dtr, weight);

					H2_xbj_Q2_pid_acc_1->Fill(Q2, x_bj, weight);
					H2_xbj_Em_pid_acc_1->Fill(emiss, x_bj, weight);
					H2_xbj_Pm_pid_acc_1->Fill(pmiss, x_bj, weight);

					H2_Em_Pm_pid_acc_1->Fill(pmiss, emiss, weight);
					if (!contam_cut)
					{
						H2_Em_Pm_c_pid_acc_1->Fill(pmiss, emiss, weight);
					}
					
					H2_Em_Q2_pid_acc_1->Fill(Q2, emiss, weight);

					H2_Q2_Pm_pid_acc_1->Fill(pmiss, Q2, weight);

					H2_W_Pm_pid_acc_1->Fill(pmiss, W, weight);
					H2_W_thrq_pid_acc_1->Fill(cth_rq / dtr, W, weight);
					H2_W_Em_pid_acc_1->Fill(emiss, W, weight);
					H2_W_Q2_pid_acc_1->Fill(Q2, W, weight);
					H2_W_xbj_pid_acc_1->Fill(x_bj, W, weight);

					H2_exptar_Pm_pid_acc_1->Fill(pmiss, exptar / dtr, weight);
					H2_eyptar_Pm_pid_acc_1->Fill(pmiss, eyptar / dtr, weight);
					H2_hxptar_Pm_pid_acc_1->Fill(pmiss, hxptar / dtr, weight);
					H2_hyptar_Pm_pid_acc_1->Fill(pmiss, hyptar / dtr, weight);

					H2_alpha_Pm_per_pid_acc_1->Fill(Pm_per, alpha_v, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_SRC_thrq_cut(0), (cth_rq / dtr), Cuts::get_SRC_thrq_cut(1), Cuts::get_SRC_xbj_cut(0), x_bj, Cuts::get_SRC_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_cthrq_Pm_pid_acc_2->Fill(pmiss, cos(cth_rq), weight);
					H2_thrq_Pm_pid_acc_2->Fill(pmiss, cth_rq / dtr, weight);
					H2_thrq_Em_pid_acc_2->Fill(emiss, cth_rq / dtr, weight);
					H2_thrq_Q2_pid_acc_2->Fill(Q2, cth_rq / dtr, weight);
					H2_thrq_xbj_pid_acc_2->Fill(x_bj, cth_rq / dtr, weight);

					H2_xbj_Q2_pid_acc_2->Fill(Q2, x_bj, weight);
					H2_xbj_Em_pid_acc_2->Fill(emiss, x_bj, weight);
					H2_xbj_Pm_pid_acc_2->Fill(pmiss, x_bj, weight);

					H2_Em_Pm_pid_acc_2->Fill(pmiss, emiss, weight);
					if (!contam_cut)
					{
						H2_Em_Pm_c_pid_acc_2->Fill(pmiss, emiss, weight);
					}
					
					H2_Em_Q2_pid_acc_2->Fill(Q2, emiss, weight);

					H2_Q2_Pm_pid_acc_2->Fill(pmiss, Q2, weight);

					H2_W_Pm_pid_acc_2->Fill(pmiss, W, weight);
					H2_W_thrq_pid_acc_2->Fill(cth_rq / dtr, W, weight);
					H2_W_Em_pid_acc_2->Fill(emiss, W, weight);
					H2_W_Q2_pid_acc_2->Fill(Q2, W, weight);
					H2_W_xbj_pid_acc_2->Fill(x_bj, W, weight);

					H2_exptar_Pm_pid_acc_2->Fill(pmiss, exptar / dtr, weight);
					H2_eyptar_Pm_pid_acc_2->Fill(pmiss, eyptar / dtr, weight);
					H2_hxptar_Pm_pid_acc_2->Fill(pmiss, hxptar / dtr, weight);
					H2_hyptar_Pm_pid_acc_2->Fill(pmiss, hyptar / dtr, weight);

					H2_alpha_Pm_per_pid_acc_2->Fill(Pm_per, alpha_v, weight);
				}
			}
			else if (user_runtype == "MF")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_cthrq_Pm_pid_acc_1->Fill(pmiss, cos(cth_rq), weight);
					H2_thrq_Pm_pid_acc_1->Fill(pmiss, cth_rq / dtr, weight);
					H2_thrq_Em_pid_acc_1->Fill(emiss, cth_rq / dtr, weight);
					H2_thrq_Q2_pid_acc_1->Fill(Q2, cth_rq / dtr, weight);
					H2_thrq_xbj_pid_acc_1->Fill(x_bj, cth_rq / dtr, weight);

					H2_xbj_Q2_pid_acc_1->Fill(Q2, x_bj, weight);
					H2_xbj_Em_pid_acc_1->Fill(emiss, x_bj, weight);
					H2_xbj_Pm_pid_acc_1->Fill(pmiss, x_bj, weight);

					H2_Em_Pm_pid_acc_1->Fill(pmiss, emiss, weight);
					if (!contam_cut)
					{
						H2_Em_Pm_c_pid_acc_1->Fill(pmiss, emiss, weight);
					}
					
					H2_Em_Q2_pid_acc_1->Fill(Q2, emiss, weight);

					H2_Q2_Pm_pid_acc_1->Fill(pmiss, Q2, weight);

					H2_W_Pm_pid_acc_1->Fill(pmiss, W, weight);
					H2_W_thrq_pid_acc_1->Fill(cth_rq / dtr, W, weight);
					H2_W_Em_pid_acc_1->Fill(emiss, W, weight);
					H2_W_Q2_pid_acc_1->Fill(Q2, W, weight);
					H2_W_xbj_pid_acc_1->Fill(x_bj, W, weight);

					H2_exptar_Pm_pid_acc_1->Fill(pmiss, exptar / dtr, weight);
					H2_eyptar_Pm_pid_acc_1->Fill(pmiss, eyptar / dtr, weight);
					H2_hxptar_Pm_pid_acc_1->Fill(pmiss, hxptar / dtr, weight);
					H2_hyptar_Pm_pid_acc_1->Fill(pmiss, hyptar / dtr, weight);

					H2_alpha_Pm_per_pid_acc_1->Fill(Pm_per, alpha_v, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_MF_Pm_cut(0), pmiss, Cuts::get_MF_Pm_cut(1), Cuts::get_MF_Em_cut(0), emiss, Cuts::get_MF_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_cthrq_Pm_pid_acc_2->Fill(pmiss, cos(cth_rq), weight);
					H2_thrq_Pm_pid_acc_2->Fill(pmiss, cth_rq / dtr, weight);
					H2_thrq_Em_pid_acc_2->Fill(emiss, cth_rq / dtr, weight);
					H2_thrq_Q2_pid_acc_2->Fill(Q2, cth_rq / dtr, weight);
					H2_thrq_xbj_pid_acc_2->Fill(x_bj, cth_rq / dtr, weight);

					H2_xbj_Q2_pid_acc_2->Fill(Q2, x_bj, weight);
					H2_xbj_Em_pid_acc_2->Fill(emiss, x_bj, weight);
					H2_xbj_Pm_pid_acc_2->Fill(pmiss, x_bj, weight);

					H2_Em_Pm_pid_acc_2->Fill(pmiss, emiss, weight);
					if (!contam_cut)
					{
						H2_Em_Pm_c_pid_acc_2->Fill(pmiss, emiss, weight);
					}
					
					H2_Em_Q2_pid_acc_2->Fill(Q2, emiss, weight);

					H2_Q2_Pm_pid_acc_2->Fill(pmiss, Q2, weight);

					H2_W_Pm_pid_acc_2->Fill(pmiss, W, weight);
					H2_W_thrq_pid_acc_2->Fill(cth_rq / dtr, W, weight);
					H2_W_Em_pid_acc_2->Fill(emiss, W, weight);
					H2_W_Q2_pid_acc_2->Fill(Q2, W, weight);
					H2_W_xbj_pid_acc_2->Fill(x_bj, W, weight);

					H2_exptar_Pm_pid_acc_2->Fill(pmiss, exptar / dtr, weight);
					H2_eyptar_Pm_pid_acc_2->Fill(pmiss, eyptar / dtr, weight);
					H2_hxptar_Pm_pid_acc_2->Fill(pmiss, hxptar / dtr, weight);
					H2_hyptar_Pm_pid_acc_2->Fill(pmiss, hyptar / dtr, weight);

					H2_alpha_Pm_per_pid_acc_2->Fill(Pm_per, alpha_v, weight);
				}
			}
			else if (user_runtype == "Heep")
			{
				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_cthrq_Pm_pid_acc_1->Fill(pmiss, cos(cth_rq), weight);
					H2_thrq_Pm_pid_acc_1->Fill(pmiss, cth_rq / dtr, weight);
					H2_thrq_Em_pid_acc_1->Fill(emiss, cth_rq / dtr, weight);
					H2_thrq_Q2_pid_acc_1->Fill(Q2, cth_rq / dtr, weight);
					H2_thrq_xbj_pid_acc_1->Fill(x_bj, cth_rq / dtr, weight);

					H2_xbj_Q2_pid_acc_1->Fill(Q2, x_bj, weight);
					H2_xbj_Em_pid_acc_1->Fill(emiss, x_bj, weight);
					H2_xbj_Pm_pid_acc_1->Fill(pmiss, x_bj, weight);

					H2_Em_Pm_pid_acc_1->Fill(pmiss, emiss, weight);
					if (!contam_cut)
					{
						H2_Em_Pm_c_pid_acc_1->Fill(pmiss, emiss, weight);
					}
					
					H2_Em_Q2_pid_acc_1->Fill(Q2, emiss, weight);

					H2_Q2_Pm_pid_acc_1->Fill(pmiss, Q2, weight);

					H2_W_Pm_pid_acc_1->Fill(pmiss, W, weight);
					H2_W_thrq_pid_acc_1->Fill(cth_rq / dtr, W, weight);
					H2_W_Em_pid_acc_1->Fill(emiss, W, weight);
					H2_W_Q2_pid_acc_1->Fill(Q2, W, weight);
					H2_W_xbj_pid_acc_1->Fill(x_bj, W, weight);

					H2_exptar_Pm_pid_acc_1->Fill(pmiss, exptar / dtr, weight);
					H2_eyptar_Pm_pid_acc_1->Fill(pmiss, eyptar / dtr, weight);
					H2_hxptar_Pm_pid_acc_1->Fill(pmiss, hxptar / dtr, weight);
					H2_hyptar_Pm_pid_acc_1->Fill(pmiss, hyptar / dtr, weight);

					H2_alpha_Pm_per_pid_acc_1->Fill(Pm_per, alpha_v, weight);
				}

				if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && bevtyp && Cuts::Cut(Cuts::get_heep_Em_cut(0), emiss, Cuts::get_heep_Em_cut(1), Cuts::get_heep_xbj_cut(0), x_bj, Cuts::get_heep_xbj_cut(1)) && ctime && HMS_TRK && SHMS_TRK)
				{
					H2_cthrq_Pm_pid_acc_2->Fill(pmiss, cos(cth_rq), weight);
					H2_thrq_Pm_pid_acc_2->Fill(pmiss, cth_rq / dtr, weight);
					H2_thrq_Em_pid_acc_2->Fill(emiss, cth_rq / dtr, weight);
					H2_thrq_Q2_pid_acc_2->Fill(Q2, cth_rq / dtr, weight);
					H2_thrq_xbj_pid_acc_2->Fill(x_bj, cth_rq / dtr, weight);

					H2_xbj_Q2_pid_acc_2->Fill(Q2, x_bj, weight);
					H2_xbj_Em_pid_acc_2->Fill(emiss, x_bj, weight);
					H2_xbj_Pm_pid_acc_2->Fill(pmiss, x_bj, weight);

					H2_Em_Pm_pid_acc_2->Fill(pmiss, emiss, weight);
					if (!contam_cut)
					{
						H2_Em_Pm_c_pid_acc_2->Fill(pmiss, emiss, weight);
					}
					
					H2_Em_Q2_pid_acc_2->Fill(Q2, emiss, weight);

					H2_Q2_Pm_pid_acc_2->Fill(pmiss, Q2, weight);

					H2_W_Pm_pid_acc_2->Fill(pmiss, W, weight);
					H2_W_thrq_pid_acc_2->Fill(cth_rq / dtr, W, weight);
					H2_W_Em_pid_acc_2->Fill(emiss, W, weight);
					H2_W_Q2_pid_acc_2->Fill(Q2, W, weight);
					H2_W_xbj_pid_acc_2->Fill(x_bj, W, weight);

					H2_exptar_Pm_pid_acc_2->Fill(pmiss, exptar / dtr, weight);
					H2_eyptar_Pm_pid_acc_2->Fill(pmiss, eyptar / dtr, weight);
					H2_hxptar_Pm_pid_acc_2->Fill(pmiss, hxptar / dtr, weight);
					H2_hyptar_Pm_pid_acc_2->Fill(pmiss, hyptar / dtr, weight);

					H2_alpha_Pm_per_pid_acc_2->Fill(Pm_per, alpha_v, weight);
				}
			}
			

			if (!contam_cut && pid_hms && pid_shms && Real_Golden && accp_hms && accp_shms && accp_z && kin && bevtyp && ctime && HMS_TRK && SHMS_TRK)
			{
				H2_cthrq_Pm_pid_acc_kin->Fill(pmiss, cos(cth_rq), weight);
				H2_thrq_Pm_pid_acc_kin->Fill(pmiss, cth_rq / dtr, weight);
				H2_thrq_Em_pid_acc_kin->Fill(emiss, cth_rq / dtr, weight);
				H2_thrq_Q2_pid_acc_kin->Fill(Q2, cth_rq / dtr, weight);
				H2_thrq_xbj_pid_acc_kin->Fill(x_bj, cth_rq / dtr, weight);

				H2_xbj_Q2_pid_acc_kin->Fill(Q2, x_bj, weight);
				H2_xbj_Em_pid_acc_kin->Fill(emiss, x_bj, weight);
				H2_xbj_Pm_pid_acc_kin->Fill(pmiss, x_bj, weight);

				H2_Em_Pm_pid_acc_kin->Fill(pmiss, emiss, weight);
				if (!contam_cut)
				{
					H2_Em_Pm_c_pid_acc_kin->Fill(pmiss, emiss, weight);
				}
				
				H2_Em_Q2_pid_acc_kin->Fill(Q2, emiss, weight);

				H2_Q2_Pm_pid_acc_kin->Fill(pmiss, Q2, weight);

				H2_W_Pm_pid_acc_kin->Fill(pmiss, W, weight);
				H2_W_thrq_pid_acc_kin->Fill(cth_rq / dtr, W, weight);
				H2_W_Em_pid_acc_kin->Fill(emiss, W, weight);
				H2_W_Q2_pid_acc_kin->Fill(Q2, W, weight);
				H2_W_xbj_pid_acc_kin->Fill(x_bj, W, weight);

				H2_exptar_Pm_pid_acc_kin->Fill(pmiss, exptar / dtr, weight);
				H2_eyptar_Pm_pid_acc_kin->Fill(pmiss, eyptar / dtr, weight);
				H2_hxptar_Pm_pid_acc_kin->Fill(pmiss, hxptar / dtr, weight);
				H2_hyptar_Pm_pid_acc_kin->Fill(pmiss, hyptar / dtr, weight);

				H2_alpha_Pm_per_pid_acc_kin->Fill(Pm_per, alpha_v, weight);
			}



			if (numLoops % 100000 == 0 ) { cout << numLoops << endl; } //leave the loop after completing the loop max_loops times
			//if (numLoops == 200000 ) { break; } //leave the loop after completing the loop max_loops times

			numLoops++;
		}

		inROOT->Close();
		sleep(1);
	}//End Loop

	

	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	//-----------------------------------------------------------------------Write Out Histograms-----------------------------------------------------------------------
	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	
	
	
	//--------------------------------------
	//PID Histograms Bins
	//--------------------------------------
	//----- HMS -----
	double y2f = 0;
	double y1f = 0;

	bool not_a_flag = false;

	//--------------------------------------
	//Primary Electron Kinematics
	//--------------------------------------
	Double_t Q2_min = -0.1, Q2_max = 0.1; 	bool Q2_flag = false;
	Double_t xbj_min = -0.1, xbj_max = 0.1;	bool xbj_flag = false;
	Double_t Em_min = -0.1, Em_max = 0.1;	bool Em_flag = false;
	Double_t W_min = -0.1, W_max = 0.1;	bool W_flag = false;
	Double_t MM_min = -0.1, MM_max = 0.1;	bool MM_flag = false;
	Double_t Pm_min = -0.1, Pm_max = 0.1;	bool Pm_flag = false;
	Double_t thrq_min = -0.1, thrq_max = 0.1;	bool thrq_flag = false;


	if (user_runtype == "Heep")
	{
		Q2_flag = Cuts::get_heep_Q2_flag();			Q2_min = Cuts::get_heep_Q2_cut(0);		Q2_max = Cuts::get_heep_Q2_cut(1);
		xbj_flag = Cuts::get_heep_xbj_flag();		xbj_min = Cuts::get_heep_xbj_cut(0);	xbj_max = Cuts::get_heep_xbj_cut(1);
		Em_flag = Cuts::get_heep_Em_flag();			Em_min = Cuts::get_heep_Em_cut(0);		Em_max = Cuts::get_heep_Em_cut(1);
		W_flag = Cuts::get_heep_W_flag();			W_min = Cuts::get_heep_W_cut(0);		W_max = Cuts::get_heep_W_cut(1);
		MM_flag = Cuts::get_heep_MM_flag();			MM_min = Cuts::get_heep_MM_cut(0);		MM_max = Cuts::get_heep_MM_cut(1);
	}
	else if (user_runtype == "MF")
	{
		Q2_flag = Cuts::get_MF_Q2_flag();			Q2_min = Cuts::get_MF_Q2_cut(0);		Q2_max = Cuts::get_MF_Q2_cut(1);
		Pm_flag = Cuts::get_MF_Pm_flag();			Pm_min = Cuts::get_MF_Pm_cut(0);		Pm_max = Cuts::get_MF_Pm_cut(1);
		thrq_flag = Cuts::get_MF_thrq_flag();		thrq_min = Cuts::get_MF_thrq_cut(0);	thrq_max = Cuts::get_MF_thrq_cut(1);

		if (user_target == "D2")
		{
			Em_flag = Cuts::get_d2MF_Em_flag();		Em_min = Cuts::get_d2MF_Em_cut(0);		Em_max = Cuts::get_d2MF_Em_cut(1);
		}
		else
		{
			Em_flag = Cuts::get_MF_Em_flag();		Em_min = Cuts::get_MF_Em_cut(0);		Em_max = Cuts::get_MF_Em_cut(1);
		}
	}
	else if (user_runtype == "SRC")
	{
		Q2_flag = Cuts::get_SRC_Q2_flag();			Q2_min = Cuts::get_SRC_Q2_cut(0);		Q2_max = Cuts::get_SRC_Q2_cut(1);
		Pm_flag = Cuts::get_SRC_Pm_flag();			Pm_min = Cuts::get_SRC_Pm_cut(0);		Pm_max = Cuts::get_SRC_Pm_cut(1);
		xbj_flag = Cuts::get_SRC_xbj_flag();		xbj_min = Cuts::get_SRC_xbj_cut(0);		xbj_max = Cuts::get_SRC_xbj_cut(1);
		thrq_flag = Cuts::get_SRC_thrq_flag();		thrq_min = Cuts::get_SRC_thrq_cut(0);	thrq_max = Cuts::get_SRC_thrq_cut(1);

		if (user_target == "D2")
		{
			Em_flag = Cuts::get_d2SRC_Em_flag();	Em_min = Cuts::get_d2SRC_Em_cut(0);		Em_max = Cuts::get_d2SRC_Em_cut(1);
		}
		else
		{
			Em_flag = Cuts::get_SRC_Em_flag();		Em_min = Cuts::get_SRC_Em_cut(0);		Em_max = Cuts::get_SRC_Em_cut(1);
		}
	}





	//----- Coin -----
	bool ePctime_flag = Cuts::get_ePctime_flag();
	Double_t ePctime_min = Cuts::get_ePctime_cut(0); Double_t ePctime_max = Cuts::get_ePctime_cut(1);
	//H1_ep_ctime->Scale(scale); H1_ep_ctime_pid->Scale(scale); H1_ep_ctime_pid_acc->Scale(scale); H1_ep_ctime_pid_acc_kin->Scale(scale);
	be_and_af(H1_ep_ctime, H1_ep_ctime_pid, H1_ep_ctime_pid_acc, H1_ep_ctime_pid_acc_kin, ePctime_min, 0, ePctime_min, y1f, ePctime_max, 0, ePctime_max, y2f, "ep_ctime.png", "ep_ctime", ePctime_flag, outROOT);
	
	be_and_af(H1_ep_ctime, H1_ep_ctime_pid_alt, H1_ep_ctime_pid_acc_alt, H1_ep_ctime_pid_acc_kin_alt, ePctime_min, 0, ePctime_min, y1f, ePctime_max, 0, ePctime_max, y2f, "ep_ctime_alt.png", "ep_ctime_alt", ePctime_flag, outROOT);


	

	//----- SHMS -----
	bool pdc_ntrk_flag = Cuts::get_pdc_ntrk_flag();
	Double_t pdc_ntrk_min = Cuts::get_pdc_ntrk_cut(0);
	be_and_af(H1_pdc_ntrk, H1_pdc_ntrk_pid, H1_pdc_ntrk_pid_acc, H1_pdc_ntrk_pid_acc_kin, pdc_ntrk_min, 0, pdc_ntrk_min, y1f, pdc_ntrk_min, 0, pdc_ntrk_min, y2f, "pdc_ntrk.png", "pdc_ntrk", pdc_ntrk_flag, outROOT);
	
	bool pScinGood_flag = Cuts::get_pScinGood_flag();
	be_and_af(H1_pScinGood, H1_pScinGood_pid, H1_pScinGood_pid_acc, H1_pScinGood_pid_acc_kin, 0, 0, 0, y1f, 0, 0, 0, y2f, "pScinGood.png", "pScinGood", pScinGood_flag, outROOT);

	bool petot_trkNorm_pid_flag = Cuts::get_petot_trkNorm_pid_flag();
	Double_t petot_trkNorm_pid_min = Cuts::get_petot_trkNorm_pid_cut(0), petot_trkNorm_pid_max = Cuts::get_petot_trkNorm_pid_cut(1);
	be_and_af(H1_pCalEtotTrkNorm, H1_pCalEtotTrkNorm_pid, H1_pCalEtotTrkNorm_pid_acc, H1_pCalEtotTrkNorm_pid_acc_kin, petot_trkNorm_pid_min, 0, petot_trkNorm_pid_min, y1f, petot_trkNorm_pid_max, 0, petot_trkNorm_pid_max, y2f, "pCalEtotTrkNorm.png", "pCalEtotTrkNorm", petot_trkNorm_pid_flag, outROOT);

	bool pngcer_pid_flag = Cuts::get_pngcer_pid_flag();
	Double_t pngcer_pid_min = Cuts::get_pngcer_pid_cut(0), pngcer_pid_max = Cuts::get_pngcer_pid_cut(1);
	be_and_af(H1_pNGCerNpeSum, H1_pNGCerNpeSum_pid, H1_pNGCerNpeSum_pid_acc, H1_pNGCerNpeSum_pid_acc_kin, pngcer_pid_min, 0, pngcer_pid_min, y1f, pngcer_pid_max, 0, pngcer_pid_max, y2f, "pNGCerNpeSum.png", "pNGCerNpeSum", pngcer_pid_flag, outROOT);

	bool phgcer_pid_flag = Cuts::get_phgcer_pid_flag();
	Double_t phgcer_pid_min = Cuts::get_phgcer_pid_cut(0), phgcer_pid_max = Cuts::get_phgcer_pid_cut(1);
	be_and_af(H1_pHGCerNpeSum, H1_pHGCerNpeSum_pid, H1_pHGCerNpeSum_pid_acc, H1_pHGCerNpeSum_pid_acc_kin, phgcer_pid_min, 0, phgcer_pid_min, y1f, phgcer_pid_max, 0, phgcer_pid_max, y2f, "pHGCerNpeSum.png", "pHGCerNpeSum", phgcer_pid_flag, outROOT);
		
	be_and_af(H1_pCalEtotNorm, H1_pCalEtotNorm_pid, H1_pCalEtotNorm_pid_acc, H1_pCalEtotNorm_pid_acc_kin, 0, 0, 0, y1f, 0, 0, 0, y2f, "pCalEtotNorm.png", "pCalEtotNorm", not_a_flag, outROOT);
	be_and_af(H1_pHodBetaNtrk, H1_pHodBetaNtrk_pid, H1_pHodBetaNtrk_pid_acc, H1_pHodBetaNtrk_pid_acc_kin, 0, 0, 0, y1f, 0, 0, 0, y2f, "pHodBetaNtrk.png", "pHodBetaNtrk", not_a_flag, outROOT);
	
	//	H1_pHodBetaTrk_pid_acc_kin->Fill(pHodBetaTrk, weight);




	
	

	//----- HMS -----
	bool hdc_ntrk_flag = Cuts::get_hdc_ntrk_flag();
	Double_t hdc_ntrk_min = Cuts::get_hdc_ntrk_cut(0);
	be_and_af(H1_hdc_ntrk, H1_hdc_ntrk_pid, H1_hdc_ntrk_pid_acc, H1_hdc_ntrk_pid_acc_kin, petot_trkNorm_pid_min, 0, petot_trkNorm_pid_min, y1f, petot_trkNorm_pid_min, 0, petot_trkNorm_pid_min, y2f, "hdc_ntrk.png", "hdc_ntrk", hdc_ntrk_flag, outROOT);

	bool hScinGood_flag = Cuts::get_hScinGood_flag();
	be_and_af(H1_hScinGood, H1_hScinGood_pid, H1_hScinGood_pid_acc, H1_hScinGood_pid_acc_kin, 0, 0, 0, y1f, 0, 0, 0, y2f, "hScinGood.png", "hScinGood", hScinGood_flag, outROOT);
	
	bool hetot_trkNorm_pid_flag = Cuts::get_hetot_trkNorm_pid_flag();
	Double_t hetot_trkNorm_pid_min = Cuts::get_hetot_trkNorm_pid_cut(0), hetot_trkNorm_pid_max = Cuts::get_hetot_trkNorm_pid_cut(1);
	be_and_af(H1_hCalEtotTrkNorm, H1_hCalEtotTrkNorm_pid, H1_hCalEtotTrkNorm_pid_acc, H1_hCalEtotTrkNorm_pid_acc_kin, hetot_trkNorm_pid_min, 0, hetot_trkNorm_pid_min, y1f, hetot_trkNorm_pid_max, 0, hetot_trkNorm_pid_max, y2f, "hCalEtotTrkNorm.png", "hCalEtotTrkNorm", hetot_trkNorm_pid_flag, outROOT);

	bool hcer_pid_flag = Cuts::get_hcer_pid_flag();
	Double_t hcer_pid_min = Cuts::get_hcer_pid_cut(0), hcer_pid_max = Cuts::get_hcer_pid_cut(1);
	be_and_af(H1_hCerNpeSum, H1_hCerNpeSum_pid, H1_hCerNpeSum_pid_acc, H1_hCerNpeSum_pid_acc_kin, hcer_pid_min, 0, hcer_pid_min, y1f, hcer_pid_min, 0, hcer_pid_min, y2f, "hCerNpeSum.png", "hCerNpeSum", hcer_pid_flag, outROOT);

	be_and_af(H1_hCalEtotNorm, H1_hCalEtotNorm_pid, H1_hCalEtotNorm_pid_acc, H1_hCalEtotNorm_pid_acc_kin, 0, 0, 0, y1f, 0, 0, 0, y2f, "hCalEtotNorm.png", "hCalEtotNorm", not_a_flag, outROOT);
	be_and_af(H1_hHodBetaNtrk, H1_hHodBetaNtrk_pid, H1_hHodBetaNtrk_pid_acc, H1_hHodBetaNtrk_pid_acc_kin, 0, 0, 0, y1f, 0, 0, 0, y2f, "hHodBetaNtrk.png", "hHodBetaNtrk", not_a_flag, outROOT);
	
	//H1_hHodBetaTrk_pid_acc_kin->Fill(hHodBetaTrk, weight);


	//----- Hadron Arm Focal Plane -----
	be_and_af(H1_hxfp_pid_acc, H1_hxfp_pid_acc_1, H1_hxfp_pid_acc_2, H1_hxfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hxfp.png", "hxfp", not_a_flag, outROOT);
	be_and_af(H1_hxpfp_pid_acc, H1_hxpfp_pid_acc_1, H1_hxpfp_pid_acc_2, H1_hxpfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hxpfp.png", "hxpfp", not_a_flag, outROOT);
	be_and_af(H1_hyfp_pid_acc, H1_hyfp_pid_acc_1, H1_hyfp_pid_acc_2, H1_hyfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hyfp.png", "hyfp", not_a_flag, outROOT);
	be_and_af(H1_hypfp_pid_acc, H1_hypfp_pid_acc_1, H1_hypfp_pid_acc_2, H1_hypfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hypfp.png", "hypfp", not_a_flag, outROOT);
	
	//----- Hadron Arm Reconstructed Quantities -----
	be_and_af(H1_hytar_pid_acc, H1_hytar_pid_acc_1, H1_hytar_pid_acc_2, H1_hytar_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hytar.png", "hytar", not_a_flag, outROOT);
	
	bool hyptar_flag = Cuts::get_hyptar_flag();
	Double_t hyptar_min = Cuts::get_hyptar_cut(0), hyptar_max = Cuts::get_hyptar_cut(1);
	be_and_af(H1_hyptar_pid_acc, H1_hyptar_pid_acc_1, H1_hyptar_pid_acc_2, H1_hyptar_pid_acc_kin, hyptar_min, 0, hyptar_min, y1f, hyptar_max, 0, hyptar_max, y2f, "hyptar.png", "hyptar", hyptar_flag, outROOT);
	
	bool hxptar_flag = Cuts::get_hxptar_flag();
	Double_t hxptar_min = Cuts::get_hxptar_cut(0), hxptar_max = Cuts::get_hxptar_cut(1);
	be_and_af(H1_hxptar_pid_acc, H1_hxptar_pid_acc_1, H1_hxptar_pid_acc_2, H1_hxptar_pid_acc_kin, hxptar_min, 0, hxptar_min, y1f, hxptar_max, 0, hxptar_max, y2f, "hxptar.png", "hxptar", hxptar_flag, outROOT);
	
	bool hdelta_flag = Cuts::get_hdelta_flag();
	Double_t hdelta_min = Cuts::get_hdelta_cut(0), hdelta_max = Cuts::get_hdelta_cut(1);
	be_and_af(H1_hdelta, H1_hdelta_pid, H1_hdelta_pid_acc, H1_hdelta_pid_acc_kin, hdelta_min, 0, hdelta_min, y1f, hdelta_max, 0, hdelta_max, y2f, "hdelta.png", "hdelta", hdelta_flag, outROOT);
	
	//----- Target Reconstruction (Hall Coord. System) -----
	be_and_af(H1_htarx_pid_acc, H1_htarx_pid_acc_1, H1_htarx_pid_acc_2, H1_htarx_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "htarx.png", "htarx", not_a_flag, outROOT);
	be_and_af(H1_htary_pid_acc, H1_htary_pid_acc_1, H1_htary_pid_acc_2, H1_htary_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "htary.png", "htary", not_a_flag, outROOT);
	be_and_af(H1_htarz_pid_acc, H1_htarz_pid_acc_1, H1_htarz_pid_acc_2, H1_htarz_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "htarz.png", "htarz", not_a_flag, outROOT);
	
	//----- HMS Collimator -----
	be_and_af(H1_hXColl_pid_acc, H1_hXColl_pid_acc_1, H1_hXColl_pid_acc_2, H1_hXColl_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hXColl.png", "hXColl", not_a_flag, outROOT);
	be_and_af(H1_hYColl_pid_acc, H1_hYColl_pid_acc_1, H1_hYColl_pid_acc_2, H1_hYColl_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hYColl.png", "hYColl", not_a_flag, outROOT);


	//----- Hadron Arm Focal Plane -----
	be_and_af(H1_exfp_pid_acc, H1_exfp_pid_acc_1, H1_exfp_pid_acc_2, H1_exfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "exfp.png", "exfp", not_a_flag, outROOT);
	be_and_af(H1_expfp_pid_acc, H1_expfp_pid_acc_1, H1_expfp_pid_acc_2, H1_expfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "expfp.png", "expfp", not_a_flag, outROOT);
	be_and_af(H1_eyfp_pid_acc, H1_eyfp_pid_acc_1, H1_eyfp_pid_acc_2, H1_eyfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eyfp.png", "eyfp", not_a_flag, outROOT);
	be_and_af(H1_eypfp_pid_acc, H1_eypfp_pid_acc_1, H1_eypfp_pid_acc_2, H1_eypfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eypfp.png", "eypfp", not_a_flag, outROOT);
	
	//----- Hadron Arm Reconstructed Quantities -----
	be_and_af(H1_eytar_pid_acc, H1_eytar_pid_acc_1, H1_eytar_pid_acc_2, H1_eytar_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eytar.png", "eytar", not_a_flag, outROOT);
	
	bool eyptar_flag = Cuts::get_eyptar_flag();
	Double_t eyptar_min = Cuts::get_eyptar_cut(0), eyptar_max = Cuts::get_eyptar_cut(1);
	be_and_af(H1_eyptar_pid_acc, H1_eyptar_pid_acc_1, H1_eyptar_pid_acc_2, H1_eyptar_pid_acc_kin, eyptar_min, 0, eyptar_min, y1f, eyptar_max, 0, eyptar_max, y2f, "eyptar.png", "eyptar", eyptar_flag, outROOT);
	
	bool exptar_flag = Cuts::get_exptar_flag();
	Double_t exptar_min = Cuts::get_exptar_cut(0), exptar_max = Cuts::get_exptar_cut(1);
	be_and_af(H1_exptar_pid_acc, H1_exptar_pid_acc_1, H1_exptar_pid_acc_2, H1_exptar_pid_acc_kin, exptar_min, 0, exptar_min, y1f, exptar_max, 0, exptar_max, y2f, "exptar.png", "exptar", exptar_flag, outROOT);
	
	bool edelta_flag = Cuts::get_edelta_flag();
	Double_t edelta_min = Cuts::get_edelta_cut(0), edelta_max = Cuts::get_edelta_cut(1);
	be_and_af(H1_edelta, H1_edelta_pid, H1_edelta_pid_acc, H1_edelta_pid_acc_kin, edelta_min, 0, edelta_min, y1f, edelta_max, 0, edelta_max, y2f, "edelta.png", "edelta", edelta_flag, outROOT);
	//----- Target Reconstruction (Hall Coord. System) -----
	be_and_af(H1_etarx_pid_acc, H1_etarx_pid_acc_1, H1_etarx_pid_acc_2, H1_etarx_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "etarx.png", "etarx", not_a_flag, outROOT);
	be_and_af(H1_etary_pid_acc, H1_etary_pid_acc_1, H1_etary_pid_acc_2, H1_etary_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "etary.png", "etary", not_a_flag, outROOT);
	be_and_af(H1_etarz_pid_acc, H1_etarz_pid_acc_1, H1_etarz_pid_acc_2, H1_etarz_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "etarz.png", "etarz", not_a_flag, outROOT);
	//----- SHMS Collimator -----
	be_and_af(H1_eXColl_pid_acc, H1_eXColl_pid_acc_1, H1_eXColl_pid_acc_2, H1_eXColl_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eXColl.png", "eXColl", not_a_flag, outROOT);
	be_and_af(H1_eYColl_pid_acc, H1_eYColl_pid_acc_1, H1_eYColl_pid_acc_2, H1_eYColl_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eYColl.png", "eYColl", not_a_flag, outROOT);


	//----- ztarDiff -----
	bool ztarDiff_flag = Cuts::get_ztarDiff_flag();
	Double_t ztarDiff_min = Cuts::get_ztarDiff_cut(0), ztarDiff_max = Cuts::get_ztarDiff_cut(1);
	be_and_af(H1_ztar_diff_pid_acc, H1_ztar_diff_pid_acc_1, H1_ztar_diff_pid_acc_2, H1_ztar_diff_pid_acc_kin, ztarDiff_min, 0, ztarDiff_min, y1f, ztarDiff_max, 0, ztarDiff_max, y2f, "ztar_diff.png", "ztar_diff", ztarDiff_flag, outROOT);


	be_and_af(H1_Q2_pid_acc, H1_Q2_pid_acc_1, H1_Q2_pid_acc_2, H1_Q2_pid_acc_kin, Q2_min, 0, Q2_min, y1f, Q2_max, 0, Q2_max, y2f, "Q2.png", "Q2", Q2_flag, outROOT);
	be_and_af(H1_xbj_pid_acc, H1_xbj_pid_acc_1, H1_xbj_pid_acc_2, H1_xbj_pid_acc_kin, xbj_min, 0, xbj_min, y1f, xbj_max, 0, xbj_max, y2f, "xbj.png", "xbj", xbj_flag, outROOT);
	be_and_af(H1_W_pid_acc, H1_W_pid_acc_1, H1_W_pid_acc_2, H1_W_pid_acc_kin, W_min, 0, W_min, y1f, W_max, 0, W_max, y2f, "W.png", "W", W_flag, outROOT);
	be_and_af(H1_nu_pid_acc, H1_nu_pid_acc_1, H1_nu_pid_acc_2, H1_nu_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "nu.png", "nu", not_a_flag, outROOT);
	be_and_af(H1_ph_q_pid_acc, H1_ph_q_pid_acc_1, H1_ph_q_pid_acc_2, H1_ph_q_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "ph_q.png", "ph_q", not_a_flag, outROOT);
	be_and_af(H1_q_pid_acc, H1_q_pid_acc_1, H1_q_pid_acc_2, H1_q_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "q.png", "q", not_a_flag, outROOT);
	be_and_af(H1_qx_pid_acc, H1_qx_pid_acc_1, H1_qx_pid_acc_2, H1_qx_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "qx.png", "qx", not_a_flag, outROOT);
	be_and_af(H1_qy_pid_acc, H1_qy_pid_acc_1, H1_qy_pid_acc_2, H1_qy_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "qy.png", "qy", not_a_flag, outROOT);
	be_and_af(H1_qz_pid_acc, H1_qz_pid_acc_1, H1_qz_pid_acc_2, H1_qz_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "qz.png", "qz", not_a_flag, outROOT);
	be_and_af(H1_th_e_pid_acc, H1_th_e_pid_acc_1, H1_th_e_pid_acc_2, H1_th_e_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "the.png", "the", not_a_flag, outROOT);
	be_and_af(H1_th_q_pid_acc, H1_th_q_pid_acc_1, H1_th_q_pid_acc_2, H1_th_q_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "thq.png", "thq", not_a_flag, outROOT);



















	be_and_af(H1_ep_ctime_rand, H1_ep_ctime_pid_rand, H1_ep_ctime_pid_acc_rand, H1_ep_ctime_pid_acc_kin_rand, ePctime_min, 0, ePctime_min, y1f, ePctime_max, 0, ePctime_max, y2f, "ep_ctime_rand.png", "ep_ctime_rand", ePctime_flag, outROOT);

	Double_t P_scale_fac = 0.2;//(4/20)
	H1_ep_ctime_rand->Scale(P_scale_fac);
	H1_ep_ctime_pid_rand->Scale(P_scale_fac);
	H1_ep_ctime_pid_acc_rand->Scale(P_scale_fac);
	H1_ep_ctime_pid_acc_kin_rand->Scale(P_scale_fac);

	H1_W_rand->Scale(P_scale_fac);
	H1_W_pid_rand->Scale(P_scale_fac);
	H1_W_pid_acc_rand->Scale(P_scale_fac);
	H1_W_pid_acc_kin_rand->Scale(P_scale_fac);

	H1_Q2_rand->Scale(P_scale_fac);
	H1_Q2_pid_rand->Scale(P_scale_fac);
	H1_Q2_pid_acc_rand->Scale(P_scale_fac);
	H1_Q2_pid_acc_kin_rand->Scale(P_scale_fac);

	H1_xbj_rand->Scale(P_scale_fac);
	H1_xbj_pid_rand->Scale(P_scale_fac);
	H1_xbj_pid_acc_rand->Scale(P_scale_fac);
	H1_xbj_pid_acc_kin_rand->Scale(P_scale_fac);

	H1_nu_rand->Scale(P_scale_fac);
	H1_nu_pid_rand->Scale(P_scale_fac);
	H1_nu_pid_acc_rand->Scale(P_scale_fac);
	H1_nu_pid_acc_kin_rand->Scale(P_scale_fac);

	H1_q_rand->Scale(P_scale_fac);
	H1_q_pid_rand->Scale(P_scale_fac);
	H1_q_pid_acc_rand->Scale(P_scale_fac);
	H1_q_pid_acc_kin_rand->Scale(P_scale_fac);

	H1_emiss_rand->Scale(P_scale_fac);
	H1_emiss_pid_rand->Scale(P_scale_fac);
	H1_emiss_pid_acc_rand->Scale(P_scale_fac);
	H1_emiss_pid_acc_kin_rand->Scale(P_scale_fac);

	H1_pmiss_rand->Scale(P_scale_fac);
	H1_pmiss_pid_rand->Scale(P_scale_fac);
	H1_pmiss_pid_acc_rand->Scale(P_scale_fac);
	H1_pmiss_pid_acc_kin_rand->Scale(P_scale_fac);
	H1_pmiss_ratio_pid_acc_kin_rand->Scale(P_scale_fac);
	H1_pmiss_ratio_pid_acc_kin_rand_full->Scale(P_scale_fac);
	H1_Q2_ratio_pid_acc_kin_rand_full->Scale(P_scale_fac);
	H1_xbj_ratio_pid_acc_kin_rand_full->Scale(P_scale_fac);
	H1_thrq_ratio_pid_acc_kin_rand_full->Scale(P_scale_fac);
	H1_W_ratio_pid_acc_kin_rand_full->Scale(P_scale_fac);

	H1_pmiss_raw_rand->Scale(P_scale_fac);
	H1_pmiss_raw_pid_rand->Scale(P_scale_fac);
	H1_pmiss_raw_pid_acc_rand->Scale(P_scale_fac);
	H1_pmiss_raw_pid_acc_kin_rand->Scale(P_scale_fac);

	H1_pmiss_c_rand->Scale(P_scale_fac);
	H1_pmiss_c_pid_rand->Scale(P_scale_fac);
	H1_pmiss_c_pid_acc_rand->Scale(P_scale_fac);
	H1_pmiss_c_pid_acc_kin_rand->Scale(P_scale_fac);

	H1_th_rq_rand->Scale(P_scale_fac);
	H1_th_rq_pid_rand->Scale(P_scale_fac);
	H1_th_rq_pid_acc_rand->Scale(P_scale_fac);
	H1_th_rq_pid_acc_kin_rand->Scale(P_scale_fac);
	

	H1_ep_ctime_sub->Add(H1_ep_ctime, H1_ep_ctime_rand, 1, -1);
	H1_ep_ctime_pid_sub->Add(H1_ep_ctime_pid, H1_ep_ctime_pid_rand, 1, -1);
	H1_ep_ctime_pid_acc_sub->Add(H1_ep_ctime_pid_acc, H1_ep_ctime_pid_acc_rand, 1, -1);
	H1_ep_ctime_pid_acc_kin_sub->Add(H1_ep_ctime_pid_acc_kin, H1_ep_ctime_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_ep_ctime_sub, H1_ep_ctime_pid_sub, H1_ep_ctime_pid_acc_sub, H1_ep_ctime_pid_acc_kin_sub, ePctime_min, 0, ePctime_min, y1f, ePctime_max, 0, ePctime_max, y2f, "ep_ctime_sub.png", "ep_ctime_sub", ePctime_flag, outROOT);

	H1_W_sub->Add(H1_W, H1_W_rand, 1, -1);
	H1_W_pid_sub->Add(H1_W_pid, H1_W_pid_rand, 1, -1);
	H1_W_pid_acc_sub->Add(H1_W_pid_acc, H1_W_pid_acc_rand, 1, -1);
	H1_W_pid_acc_kin_sub->Add(H1_W_pid_acc_kin, H1_W_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_W_rand, H1_W_pid_rand, H1_W_pid_acc_rand, H1_W_pid_acc_kin_rand, W_min, 0, W_min, y1f, W_max, 0, W_max, y2f, "W_rand.png", "W_rand", W_flag, outROOT);
	be_and_af(H1_W_sub, H1_W_pid_sub, H1_W_pid_acc_sub, H1_W_pid_acc_kin_sub, W_min, 0, W_min, y1f, W_max, 0, W_max, y2f, "W_sub.png", "W_sub", W_flag, outROOT);

	H1_Q2_sub->Add(H1_Q2, H1_Q2_rand, 1, -1);
	H1_Q2_pid_sub->Add(H1_Q2_pid, H1_Q2_pid_rand, 1, -1);
	H1_Q2_pid_acc_sub->Add(H1_Q2_pid_acc, H1_Q2_pid_acc_rand, 1, -1);
	H1_Q2_pid_acc_kin_sub->Add(H1_Q2_pid_acc_kin, H1_Q2_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_Q2_rand, H1_Q2_pid_rand, H1_Q2_pid_acc_rand, H1_Q2_pid_acc_kin_rand, Q2_min, 0, Q2_min, y1f, Q2_max, 0, Q2_max, y2f, "Q2_rand.png", "Q2_rand", Q2_flag, outROOT);
	be_and_af(H1_Q2_sub, H1_Q2_pid_sub, H1_Q2_pid_acc_sub, H1_Q2_pid_acc_kin_sub, Q2_min, 0, Q2_min, y1f, Q2_max, 0, Q2_max, y2f, "Q2_sub.png", "Q2_sub", Q2_flag, outROOT);

	H1_xbj_sub->Add(H1_xbj, H1_xbj_rand, 1, -1);
	H1_xbj_pid_sub->Add(H1_xbj_pid, H1_xbj_pid_rand, 1, -1);
	H1_xbj_pid_acc_sub->Add(H1_xbj_pid_acc, H1_xbj_pid_acc_rand, 1, -1);
	H1_xbj_pid_acc_kin_sub->Add(H1_xbj_pid_acc_kin, H1_xbj_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_xbj_rand, H1_xbj_pid_rand, H1_xbj_pid_acc_rand, H1_xbj_pid_acc_kin_rand, xbj_min, 0, xbj_min, y1f, xbj_max, 0, xbj_max, y2f, "xbj_rand.png", "xbj_rand", xbj_flag, outROOT);
	be_and_af(H1_xbj_sub, H1_xbj_pid_sub, H1_xbj_pid_acc_sub, H1_xbj_pid_acc_kin_sub, xbj_min, 0, xbj_min, y1f, xbj_max, 0, xbj_max, y2f, "xbj_sub.png", "xbj_sub", xbj_flag, outROOT);

	H1_nu_sub->Add(H1_nu, H1_nu_rand, 1, -1);
	H1_nu_pid_sub->Add(H1_nu_pid, H1_nu_pid_rand, 1, -1);
	H1_nu_pid_acc_sub->Add(H1_nu_pid_acc, H1_nu_pid_acc_rand, 1, -1);
	H1_nu_pid_acc_kin_sub->Add(H1_nu_pid_acc_kin, H1_nu_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_nu_rand, H1_nu_pid_rand, H1_nu_pid_acc_rand, H1_nu_pid_acc_kin_rand, 0, 0, 0, 0, 0, 0, 0, 0, "nu_rand.png", "nu_rand", not_a_flag, outROOT);
	be_and_af(H1_nu_sub, H1_nu_pid_sub, H1_nu_pid_acc_sub, H1_nu_pid_acc_kin_sub, 0, 0, 0, 0, 0, 0, 0, 0, "nu_sub.png", "nu_sub", not_a_flag, outROOT);

	H1_q_sub->Add(H1_q, H1_q_rand, 1, -1);
	H1_q_pid_sub->Add(H1_q_pid, H1_q_pid_rand, 1, -1);
	H1_q_pid_acc_sub->Add(H1_q_pid_acc, H1_q_pid_acc_rand, 1, -1);
	H1_q_pid_acc_kin_sub->Add(H1_q_pid_acc_kin, H1_q_pid_acc_kin_rand, 1, -1);
	




	H1_emiss_sub->Add(H1_emiss, H1_emiss_rand, 1, -1);
	H1_emiss_pid_sub->Add(H1_emiss_pid, H1_emiss_pid_rand, 1, -1);
	H1_emiss_pid_acc_sub->Add(H1_emiss_pid_acc, H1_emiss_pid_acc_rand, 1, -1);
	H1_emiss_pid_acc_kin_sub->Add(H1_emiss_pid_acc_kin, H1_emiss_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_emiss_rand, H1_emiss_pid_rand, H1_emiss_pid_acc_rand, H1_emiss_pid_acc_kin_rand, 0, 0, 0, 0, 0, 0, 0, 0, "emiss_rand.png", "emiss_rand", not_a_flag, outROOT);
	be_and_af(H1_emiss_sub, H1_emiss_pid_sub, H1_emiss_pid_acc_sub, H1_emiss_pid_acc_kin_sub, 0, 0, 0, 0, 0, 0, 0, 0, "emiss_sub.png", "emiss_sub", not_a_flag, outROOT);
	
	H1_pmiss_sub->Add(H1_pmiss, H1_pmiss_rand, 1, -1);
	H1_pmiss_pid_sub->Add(H1_pmiss_pid, H1_pmiss_pid_rand, 1, -1);
	H1_pmiss_pid_acc_sub->Add(H1_pmiss_pid_acc, H1_pmiss_pid_acc_rand, 1, -1);
	H1_pmiss_pid_acc_kin_sub->Add(H1_pmiss_pid_acc_kin, H1_pmiss_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_pmiss_rand, H1_pmiss_pid_rand, H1_pmiss_pid_acc_rand, H1_pmiss_pid_acc_kin_rand, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_rand.png", "pmiss_rand", not_a_flag, outROOT);
	be_and_af(H1_pmiss_sub, H1_pmiss_pid_sub, H1_pmiss_pid_acc_sub, H1_pmiss_pid_acc_kin_sub, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_sub.png", "pmiss_sub", not_a_flag, outROOT);

	H1_pmiss_raw_sub->Add(H1_pmiss_raw, H1_pmiss_raw_rand, 1, -1);
	H1_pmiss_raw_pid_sub->Add(H1_pmiss_raw_pid, H1_pmiss_raw_pid_rand, 1, -1);
	H1_pmiss_raw_pid_acc_sub->Add(H1_pmiss_raw_pid_acc, H1_pmiss_raw_pid_acc_rand, 1, -1);
	H1_pmiss_raw_pid_acc_kin_sub->Add(H1_pmiss_raw_pid_acc_kin, H1_pmiss_raw_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_pmiss_raw_rand, H1_pmiss_raw_pid_rand, H1_pmiss_raw_pid_acc_rand, H1_pmiss_raw_pid_acc_kin_rand, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_raw_rand.png", "pmiss_raw_rand", not_a_flag, outROOT);
	be_and_af(H1_pmiss_raw_sub, H1_pmiss_raw_pid_sub, H1_pmiss_raw_pid_acc_sub, H1_pmiss_raw_pid_acc_kin_sub, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_raw_sub.png", "pmiss_raw_sub", not_a_flag, outROOT);

	H1_pmiss_c_sub->Add(H1_pmiss_c, H1_pmiss_c_rand, 1, -1);
	H1_pmiss_c_pid_sub->Add(H1_pmiss_c_pid, H1_pmiss_c_pid_rand, 1, -1);
	H1_pmiss_c_pid_acc_sub->Add(H1_pmiss_c_pid_acc, H1_pmiss_c_pid_acc_rand, 1, -1);
	H1_pmiss_c_pid_acc_kin_sub->Add(H1_pmiss_c_pid_acc_kin, H1_pmiss_c_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_pmiss_c_rand, H1_pmiss_c_pid_rand, H1_pmiss_c_pid_acc_rand, H1_pmiss_c_pid_acc_kin_rand, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_c_rand.png", "pmiss_c_rand", not_a_flag, outROOT);
	be_and_af(H1_pmiss_c_sub, H1_pmiss_c_pid_sub, H1_pmiss_c_pid_acc_sub, H1_pmiss_c_pid_acc_kin_sub, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_c_sub.png", "pmiss_c_sub", not_a_flag, outROOT);


	H1_th_rq_sub->Add(H1_th_rq, H1_th_rq_rand, 1, -1);
	H1_th_rq_pid_sub->Add(H1_th_rq_pid, H1_th_rq_pid_rand, 1, -1);	
	H1_th_rq_pid_acc_sub->Add(H1_th_rq_pid_acc, H1_th_rq_pid_acc_rand, 1, -1); 
	H1_th_rq_pid_acc_kin_sub->Add(H1_th_rq_pid_acc_kin, H1_th_rq_pid_acc_kin_rand, 1, -1);
	be_and_af(H1_th_rq_rand, H1_th_rq_pid_rand, H1_th_rq_pid_acc_rand, H1_th_rq_pid_acc_kin_rand, 0, 0, 0, 0, 0, 0, 0, 0, "th_rq_rand.png", "th_rq_rand", not_a_flag, outROOT);
	be_and_af(H1_th_rq_sub, H1_th_rq_pid_sub, H1_th_rq_pid_acc_sub, H1_th_rq_pid_acc_kin_sub, 0, 0, 0, 0, 0, 0, 0, 0, "th_rq_sub.png", "th_rq_sub", not_a_flag, outROOT);

	

	/////////////////////////////////////////////////////////////
	//H1_pmiss_pid_acc_kin_sub->Add(H1_pmiss_pid_acc_kin, H1_pmiss_pid_acc_kin_rand, 1, -1);

	
	H1_pmiss_ratio_pid_acc_kin->Add(H1_pmiss_ratio_pid_acc_kin, H1_pmiss_ratio_pid_acc_kin_rand, 1, -1);// = new TH1F("h1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);
	H1_pmiss_ratio_pid_acc_kin_full->Add(H1_pmiss_ratio_pid_acc_kin_full, H1_pmiss_ratio_pid_acc_kin_rand_full, 1, -1);// = new TH1F("h1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);
	H1_Q2_ratio_pid_acc_kin_full->Add(H1_Q2_ratio_pid_acc_kin_full, H1_Q2_ratio_pid_acc_kin_rand_full, 1, -1);// = new TH1F("h1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);
	H1_xbj_ratio_pid_acc_kin_full->Add(H1_xbj_ratio_pid_acc_kin_full, H1_xbj_ratio_pid_acc_kin_rand_full, 1, -1);// = new TH1F("h1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);
	H1_thrq_ratio_pid_acc_kin_full->Add(H1_thrq_ratio_pid_acc_kin_full, H1_thrq_ratio_pid_acc_kin_rand_full, 1, -1);// = new TH1F("h1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);
	H1_W_ratio_pid_acc_kin_full->Add(H1_W_ratio_pid_acc_kin_full, H1_W_ratio_pid_acc_kin_rand_full, 1, -1);// = new TH1F("h1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);

	H1_pmiss_ratio_pid_acc_kin_full_mf->Add(H1_pmiss_ratio_pid_acc_kin_full_mf, H1_pmiss_ratio_pid_acc_kin_rand_full, 1, -1);// = new TH1F("h1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);
	H1_Em_ratio_pid_acc_kin_full_mf->Add(H1_Em_ratio_pid_acc_kin_full_mf, H1_pmiss_ratio_pid_acc_kin_rand_full, 1, -1);// = new TH1F("h1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);
	H1_Q2_ratio_pid_acc_kin_full_mf->Add(H1_Q2_ratio_pid_acc_kin_full_mf, H1_Q2_ratio_pid_acc_kin_rand_full, 1, -1);// = new TH1F("h1_pmiss_ratio_pid_acc_kin", "p_{miss}; p_{miss} (GeV/c)", NBins2, edges2); Sec_Kin->Add(H1_pmiss_ratio_pid_acc_kin);



	Double_t total_bins_scale;
	Double_t total_bins_raw;
	
	total_bins_scale = H1_pmiss_pid_acc_kin->GetNbinsX();

	total_bins_raw = H1_pmiss_raw_pid_acc_kin->GetNbinsX();

	total_bins_raw = H1_pmiss_c_pid_acc_kin->GetNbinsX();


	Double_t Pm_total = 0.0;
	Double_t Pm_real = 0.0;
	Double_t Pm_rand = 0.0;
	Double_t Pm_total_err = 0.0;
	Double_t Pm_real_err = 0.0;
	Double_t Pm_rand_err = 0.0;
	Double_t Pm_total_rate = 0.0;
	Double_t Pm_real_rate = 0.0;

	Pm_total = H1_pmiss_pid_acc_kin->IntegralAndError(1, total_bins_scale, Pm_total_err);
	Pm_real = H1_pmiss_pid_acc_kin_sub->IntegralAndError(1, total_bins_scale, Pm_real_err);
	Pm_rand = H1_pmiss_pid_acc_kin_rand->IntegralAndError(1, total_bins_scale, Pm_rand_err);

	//Pm_total_rate = Pm_total / total_time_bcm_cut;
	//Pm_real_rate = Pm_real / total_time_bcm_cut;

	cout << "Pm_scale_rate: " << Pm_total << endl;
	cout << "Pm_scale_real_rate: " << Pm_real << endl;

	Pm_total = H1_pmiss_raw_pid_acc_kin->IntegralAndError(1, total_bins_raw, Pm_total_err);
	Pm_real = H1_pmiss_raw_pid_acc_kin_sub->IntegralAndError(1, total_bins_raw, Pm_real_err);
	Pm_rand = H1_pmiss_raw_pid_acc_kin_rand->IntegralAndError(1, total_bins_raw, Pm_rand_err);
	
	cout << "Pm_raw_rate: " << Pm_total << endl;
	cout << "Pm_raw_real_rate: " << Pm_real << endl;


	Pm_total = H1_pmiss_c_pid_acc_kin->IntegralAndError(1, total_bins_raw, Pm_total_err);
	Pm_real = H1_pmiss_c_pid_acc_kin_sub->IntegralAndError(1, total_bins_raw, Pm_real_err);
	Pm_rand = H1_pmiss_c_pid_acc_kin_rand->IntegralAndError(1, total_bins_raw, Pm_rand_err);

	cout << "Pm_c_rate: " << Pm_total << endl;
	cout << "Pm_c_real_rate: " << Pm_real << endl;
	









	//--------------------------------------
	//Secondary Hadron Kinematics
	//--------------------------------------
	be_and_af(H1_pmiss_pid_acc, H1_pmiss_pid_acc_1, H1_pmiss_pid_acc_2, H1_pmiss_pid_acc_kin, Pm_min, 0, Pm_min, y1f, Pm_max, 0, Pm_max, y2f, "pmiss.png", "pmiss", Pm_flag, outROOT);
	
	
	be_and_af(H1_th_rq_pid_acc, H1_th_rq_pid_acc_1, H1_th_rq_pid_acc_2, H1_th_rq_pid_acc_kin, thrq_min, 0, thrq_min, y1f, thrq_max, 0, thrq_max, y2f, "th_rq.png", "th_rq", thrq_flag, outROOT);
	be_and_af(H1_cth_rq_pid_acc, H1_cth_rq_pid_acc_1, H1_cth_rq_pid_acc_2, H1_cth_rq_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "cth_rq.png", "cth_rq", false, outROOT);//true
	
	
	be_and_af(H1_emiss_pid_acc, H1_emiss_pid_acc_1, H1_emiss_pid_acc_2, H1_emiss_pid_acc_kin, Em_min, 0, Em_min, y1f, Em_max, 0, Em_max, y2f, "emiss.png", "emiss", Em_flag, outROOT);
	//be_and_af(H1_emiss_nuc_pid_acc, H1_emiss_nuc_pid_acc_1, H1_emiss_nuc_pid_acc_2, H1_emiss_nuc_pid_acc_kin, Em_min, 0, Em_min, y1f, Em_max, 0, Em_max, y2f, "emiss_nuc.png", "emiss_nuc", Em_flag, outROOT);

	
	be_and_af(H1_prec_x_pid_acc, H1_prec_x_pid_acc_1, H1_prec_x_pid_acc_2, H1_prec_x_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Prec_x.png", "Prec_x", not_a_flag, outROOT);
	be_and_af(H1_prec_y_pid_acc, H1_prec_y_pid_acc_1, H1_prec_y_pid_acc_2, H1_prec_y_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Prec_y.png", "Prec_y", not_a_flag, outROOT);
	be_and_af(H1_prec_z_pid_acc, H1_prec_z_pid_acc_1, H1_prec_z_pid_acc_2, H1_prec_z_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Prec_z.png", "Prec_z", not_a_flag, outROOT);
	be_and_af(H1_ph_rq_pid_acc, H1_ph_rq_pid_acc_1, H1_ph_rq_pid_acc_2, H1_ph_rq_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "ph_rq.png", "ph_rq", not_a_flag, outROOT);
	be_and_af(H1_ph_pq_pid_acc, H1_ph_pq_pid_acc_1, H1_ph_pq_pid_acc_2, H1_ph_pq_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "ph_pq.png", "ph_pq", not_a_flag, outROOT);
	be_and_af(H1_pmiss_x_pid_acc, H1_pmiss_x_pid_acc_1, H1_pmiss_x_pid_acc_2, H1_pmiss_x_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_x.png", "pmiss_x", not_a_flag, outROOT);
	be_and_af(H1_pmiss_y_pid_acc, H1_pmiss_y_pid_acc_1, H1_pmiss_y_pid_acc_2, H1_pmiss_y_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_y.png", "pmiss_y", not_a_flag, outROOT);
	be_and_af(H1_pmiss_z_pid_acc, H1_pmiss_z_pid_acc_1, H1_pmiss_z_pid_acc_2, H1_pmiss_z_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_z.png", "pmiss_z", not_a_flag, outROOT);
	be_and_af(H1_Tr_pid_acc, H1_Tr_pid_acc_1, H1_Tr_pid_acc_2, H1_Tr_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Tr.png", "Tr", not_a_flag, outROOT);
	be_and_af(H1_th_pq_pid_acc, H1_th_pq_pid_acc_1, H1_th_pq_pid_acc_2, H1_th_pq_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "th_pq.png", "th_pq", not_a_flag, outROOT);
	be_and_af(H1_xangle_pid_acc, H1_xangle_pid_acc_1, H1_xangle_pid_acc_2, H1_xangle_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "xangle.png", "xangle", not_a_flag, outROOT);
	be_and_af(H1_Tx_pid_acc, H1_Tx_pid_acc_1, H1_Tx_pid_acc_2, H1_Tx_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Tx.png", "Tx", not_a_flag, outROOT);
	be_and_af(H1_mmiss_pid_acc, H1_mmiss_pid_acc_1, H1_mmiss_pid_acc_2, H1_mmiss_pid_acc_kin, MM_min, 0, MM_min, y1f, MM_max, 0, MM_max, y2f, "mmiss.png", "mmiss", MM_flag, outROOT);
	be_and_af(H1_th_p_pid_acc, H1_th_p_pid_acc_1, H1_th_p_pid_acc_2, H1_th_p_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "th_p.png", "th_p", not_a_flag, outROOT);
	//be_and_af(H1_omega, H1_omega_pid_acc_1, H1_omega_pid_acc_2, H1_omega_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "omega.png", "omega", not_a_flag, outROOT);
	be_and_af(H1_Pm_par_pid_acc, H1_Pm_par_pid_acc_1, H1_Pm_par_pid_acc_2, H1_Pm_par_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Pm_par.png", "Pm_par", not_a_flag, outROOT);
	be_and_af(H1_Pm_per_pid_acc, H1_Pm_per_pid_acc_1, H1_Pm_per_pid_acc_2, H1_Pm_per_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Pm_per.png", "Pm_per", not_a_flag, outROOT);
	be_and_af(H1_alpha_pid_acc, H1_alpha_pid_acc_1, H1_alpha_pid_acc_2, H1_alpha_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "alpha.png", "alpha", not_a_flag, outROOT);


	be_and_af2d(H2_hxfp_hyfp_pid_acc, H2_hxfp_hyfp_pid_acc_1, H2_hxfp_hyfp_pid_acc_2, H2_hxfp_hyfp_pid_acc_kin, "H2_hxfp_hyfp.png", "H2_hxfp_hyfp", outROOT);
	be_and_af2d(H2_exfp_eyfp_pid_acc, H2_exfp_eyfp_pid_acc_1, H2_exfp_eyfp_pid_acc_2, H2_exfp_eyfp_pid_acc_kin, "H2_exfp_eyfp.png", "H2_exfp_eyfp", outROOT);
	
	be_and_af2d(H2_exptar_eyptar_pid_acc, H2_exptar_eyptar_pid_acc_1, H2_exptar_eyptar_pid_acc_2, H2_exptar_eyptar_pid_acc_kin, "H2_exptar_eyptar.png", "H2_exptar_eyptar", outROOT);
	be_and_af2d(H2_hxptar_hyptar_pid_acc, H2_hxptar_hyptar_pid_acc_1, H2_hxptar_hyptar_pid_acc_2, H2_hxptar_hyptar_pid_acc_kin, "H2_hxptar_hyptar.png", "H2_hxptar_hyptar", outROOT);

	be_and_af2d(H2_hxptar_exptar_pid_acc, H2_hxptar_exptar_pid_acc_1, H2_hxptar_exptar_pid_acc_2, H2_hxptar_exptar_pid_acc_kin, "H2_hxptar_exptar.png", "H2_hxptar_exptar", outROOT);
	be_and_af2d(H2_hyptar_eyptar_pid_acc, H2_hyptar_eyptar_pid_acc_1, H2_hyptar_eyptar_pid_acc_2, H2_hyptar_eyptar_pid_acc_kin, "H2_hyptar_eyptar.png", "H2_hyptar_eyptar", outROOT);
	be_and_af2d(H2_hdelta_edelta_pid_acc, H2_hdelta_edelta_pid_acc_1, H2_hdelta_edelta_pid_acc_2, H2_hdelta_edelta_pid_acc_kin, "H2_hdelta_edelta.png", "H2_hdelta_edelta", outROOT);
	be_and_af2d(H2_pCalEtotTrkNorm_edelta_pid_acc, H2_pCalEtotTrkNorm_edelta_pid_acc_1, H2_pCalEtotTrkNorm_edelta_pid_acc_2, H2_pCalEtotTrkNorm_edelta_pid_acc_kin, "H2_pCalEtotTrkNorm_edelta.png", "H2_pCalEtotTrkNorm_edelta", outROOT);
	be_and_af2d(H2_hXColl_hYColl_pid_acc, H2_hXColl_hYColl_pid_acc_1, H2_hXColl_hYColl_pid_acc_2, H2_hXColl_hYColl_pid_acc_kin, "H2_hXColl_hYColl.png", "H2_hXColl_hYColl", outROOT);
	be_and_af2d(H2_eXColl_eYColl_pid_acc, H2_eXColl_eYColl_pid_acc_1, H2_eXColl_eYColl_pid_acc_2, H2_eXColl_eYColl_pid_acc_kin, "H2_eXColl_eYColl.png", "H2_eXColl_eYColl", outROOT);
	
	be_and_af2d(H2_alpha_Pm_per_pid_acc, H2_alpha_Pm_per_pid_acc_1, H2_alpha_Pm_per_pid_acc_2, H2_alpha_Pm_per_pid_acc_kin, "H2_alpha_Pm_per.png", "H2_alpha_Pm_per", outROOT);
	
	
	
	be_and_af2d(H2_cthrq_Pm_pid_acc, H2_cthrq_Pm_pid_acc_1, H2_cthrq_Pm_pid_acc_2, H2_cthrq_Pm_pid_acc_kin, "H2_cthrq_Pm.png", "H2_cthrq_Pm", outROOT);
	be_and_af2d(H2_thrq_Pm_pid_acc, H2_thrq_Pm_pid_acc_1, H2_thrq_Pm_pid_acc_2, H2_thrq_Pm_pid_acc_kin, "H2_thrq_Pm.png", "H2_thrq_Pm", outROOT);
	be_and_af2d(H2_thrq_Em_pid_acc, H2_thrq_Em_pid_acc_1, H2_thrq_Em_pid_acc_2, H2_thrq_Em_pid_acc_kin, "H2_thrq_Em.png", "H2_thrq_Em", outROOT);
	be_and_af2d(H2_thrq_Q2_pid_acc, H2_thrq_Q2_pid_acc_1, H2_thrq_Q2_pid_acc_2, H2_thrq_Q2_pid_acc_kin, "H2_thrq_Q2.png", "H2_thrq_Q2", outROOT);
	be_and_af2d(H2_thrq_xbj_pid_acc, H2_thrq_xbj_pid_acc_1, H2_thrq_xbj_pid_acc_2, H2_thrq_xbj_pid_acc_kin, "H2_thrq_xbj.png", "H2_thrq_xbj", outROOT);

	be_and_af2d(H2_xbj_Q2_pid_acc, H2_xbj_Q2_pid_acc_1, H2_xbj_Q2_pid_acc_2, H2_xbj_Q2_pid_acc_kin, "H2_xbj_Q2.png", "H2_xbj_Q2", outROOT);
	be_and_af2d(H2_xbj_Em_pid_acc, H2_xbj_Em_pid_acc_1, H2_xbj_Em_pid_acc_2, H2_xbj_Em_pid_acc_kin, "H2_xbj_Em.png", "H2_xbj_Em", outROOT);
	be_and_af2d(H2_xbj_Pm_pid_acc, H2_xbj_Pm_pid_acc_1, H2_xbj_Pm_pid_acc_2, H2_xbj_Pm_pid_acc_kin, "H2_xbj_Pm.png", "H2_xbj_Pm", outROOT);
	
	be_and_af2d(H2_Em_Pm_pid_acc, H2_Em_Pm_pid_acc_1, H2_Em_Pm_pid_acc_2, H2_Em_Pm_pid_acc_kin, "H2_Em_Pm.png", "H2_Em_Pm", outROOT);
	be_and_af2d(H2_Em_Pm_c_pid_acc, H2_Em_Pm_c_pid_acc_1, H2_Em_Pm_c_pid_acc_2, H2_Em_Pm_c_pid_acc_kin, "H2_Em_Pm_c.png", "H2_Em_Pm_c", outROOT);
	be_and_af2d(H2_Em_Q2_pid_acc, H2_Em_Q2_pid_acc_1, H2_Em_Q2_pid_acc_2, H2_Em_Q2_pid_acc_kin, "H2_Em_Q2.png", "H2_Em_Q2", outROOT);

	be_and_af2d(H2_Q2_Pm_pid_acc, H2_Q2_Pm_pid_acc_1, H2_Q2_Pm_pid_acc_2, H2_Q2_Pm_pid_acc_kin, "H2_Q2_Pm.png", "H2_Q2_Pm", outROOT);

	be_and_af2d(H2_W_Pm_pid_acc, H2_W_Pm_pid_acc_1, H2_W_Pm_pid_acc_2, H2_W_Pm_pid_acc_kin, "H2_W_Pm.png", "H2_W_Pm", outROOT);
	be_and_af2d(H2_W_thrq_pid_acc, H2_W_thrq_pid_acc_1, H2_W_thrq_pid_acc_2, H2_W_thrq_pid_acc_kin, "H2_W_thrq.png", "H2_W_thrq", outROOT);
	be_and_af2d(H2_W_Em_pid_acc, H2_W_Em_pid_acc_1, H2_W_Em_pid_acc_2, H2_W_Em_pid_acc_kin, "H2_W_Em.png", "H2_W_Em", outROOT);
	be_and_af2d(H2_W_Q2_pid_acc, H2_W_Q2_pid_acc_1, H2_W_Q2_pid_acc_2, H2_W_Q2_pid_acc_kin, "H2_W_Q2.png", "H2_W_Q2", outROOT);
	be_and_af2d(H2_W_xbj_pid_acc, H2_W_xbj_pid_acc_1, H2_W_xbj_pid_acc_2, H2_W_xbj_pid_acc_kin, "H2_W_xbj.png", "H2_W_xbj", outROOT);

	be_and_af2d(H2_exptar_Pm_pid_acc, H2_exptar_Pm_pid_acc_1, H2_exptar_Pm_pid_acc_2, H2_exptar_Pm_pid_acc_kin, "H2_exptar_Pm.png", "H2_exptar_Pm", outROOT);
	be_and_af2d(H2_eyptar_Pm_pid_acc, H2_eyptar_Pm_pid_acc_1, H2_eyptar_Pm_pid_acc_2, H2_eyptar_Pm_pid_acc_kin, "H2_eyptar_Pm.png", "H2_eyptar_Pm", outROOT);
	be_and_af2d(H2_hxptar_Pm_pid_acc, H2_hxptar_Pm_pid_acc_1, H2_hxptar_Pm_pid_acc_2, H2_hxptar_Pm_pid_acc_kin, "H2_hxptar_Pm.png", "H2_hxptar_Pm", outROOT);
	be_and_af2d(H2_hyptar_Pm_pid_acc, H2_hyptar_Pm_pid_acc_1, H2_hyptar_Pm_pid_acc_2, H2_hyptar_Pm_pid_acc_kin, "H2_hyptar_Pm.png", "H2_hyptar_Pm", outROOT);

	cout << "int1: " << H1_pmiss_ratio_pid_acc_kin->Integral() << endl;
	cout << "int2: " << H1_pmiss_c_pid_acc_kin_sub->Integral() << endl;
	
	Double_t sum22 = 0;;
	Double_t mi = H1_pmiss_c_pid_acc_kin_sub->GetXaxis()->GetXmin();
	cout << "mi: " << mi << endl;
	Double_t ma = H1_pmiss_c_pid_acc_kin_sub->GetXaxis()->GetXmax();
	cout << "ma: " << ma << endl;
	Double_t biins = H1_pmiss_c_pid_acc_kin_sub->GetXaxis()->GetNbins();
	Double_t binw = 1.3 / 100.0;
	Double_t i = mi;;
	//while(i<ma)
	//{
	for(int i = 0; i < 100; i++)
	{
		sum22 = sum22 + H1_pmiss_ratio_pid_acc_kin->GetBinContent(i);
		cout << "content: " << H1_pmiss_ratio_pid_acc_kin->GetBinContent(i) << endl;
		//i = i + binw;
		cout << "i: " << i << endl;
	}
	cout << "sum: " << sum22 << endl;

	outROOT->mkdir("PID");										//Make directories to store histograms based on Kinematic
	outROOT->cd("PID");											//Write Kinematics histos to kin_plots directory
	PIDList->Write();

	outROOT->mkdir("Prim_Kin");									//Make directories to store histograms based on Kinematic
	outROOT->cd("Prim_Kin");									//Write Kinematics histos to kin_plots directory
	Prim_Kin->Write();

	outROOT->mkdir("Sec_Kin");									//Make directories to store histograms based on Kinematic
	outROOT->cd("Sec_Kin");										//Write Kinematics histos to kin_plots directory
	Sec_Kin->Write();

	outROOT->mkdir("HMS_Accp");									//Make directories to store histograms based on Kinematic
	outROOT->cd("HMS_Accp");									//Write Kinematics histos to kin_plots directory
	HMS_Accp->Write();

	outROOT->mkdir("SHMS_Accp");								//Make directories to store histograms based on Kinematic
	outROOT->cd("SHMS_Accp");									//Write Kinematics histos to kin_plots directory
	SHMS_Accp->Write();

	outROOT->mkdir("2D_Hist");									//Make directories to store histograms based on Kinematic
	outROOT->cd("2D_Hist");										//Write Kinematics histos to kin_plots directory
	HList2->Write();

	outROOT->Close();											//Close File
	//good->Print();											//Print reduced cloned tree
	//good.Write();												//Write reduced cloned tree
}
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//-----------------------------------------------------------------------------End Main-----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=



void plotstuff(TH1F* hist1, const char* name, const char* title, TFile* outRoot)
{
	TCanvas* c = new TCanvas(name, title, 1024, 768); c->cd();
	
	hist1->Draw();

	outRoot->cd(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}



void be_and_af(TH1F* hist1, TH1F* hist2, TH1F* hist3, TH1F* hist4, double x1i, double y1i, double x1f, double y1f, double x2i, double y2i, double x2f, double y2f, const char* name, const char* title, bool line, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1024, 768); c->Divide(2,2);

	c->cd(1); hist1->Draw("HIST");

	if (line)
	{
		y1f = hist1->GetMaximum();
		y2f = hist1->GetMaximum();
		TLine* line1 = new TLine(x1i, y1i, x1f, y1f); line1->SetLineColor(2); line1->SetLineWidth(2); line1->Draw("SAME");
		TLine* line2 = new TLine(x2i, y2i, x2f, y2f); line2->SetLineColor(2); line2->SetLineWidth(2); line2->Draw("SAME");
	}

	c->cd(2); hist2->Draw("HIST SAME");

	if (line)
	{
		y1f = hist2->GetMaximum();
		y2f = hist2->GetMaximum();
		TLine* line3 = new TLine(x1i, y1i, x1f, y1f); line3->SetLineColor(2); line3->SetLineWidth(2); line3->Draw("SAME");
		TLine* line4 = new TLine(x2i, y2i, x2f, y2f); line4->SetLineColor(2); line4->SetLineWidth(2); line4->Draw("SAME");
	}

	c->cd(3); hist3->Draw("HIST SAME");

	if (line)
	{
		y1f = hist3->GetMaximum();
		y2f = hist3->GetMaximum();
		TLine* line4 = new TLine(x1i, y1i, x1f, y1f); line4->SetLineColor(2); line4->SetLineWidth(2); line4->Draw("SAME");
		TLine* line5 = new TLine(x2i, y2i, x2f, y2f); line5->SetLineColor(2); line5->SetLineWidth(2); line5->Draw("SAME");
	}

	c->cd(4); hist4->Draw("HIST SAME");

	if (line)
	{
		y1f = hist4->GetMaximum();
		y2f = hist4->GetMaximum();
		TLine* line6 = new TLine(x1i, y1i, x1f, y1f); line6->SetLineColor(2); line6->SetLineWidth(2); line6->Draw("SAME");
		TLine* line7 = new TLine(x2i, y2i, x2f, y2f); line7->SetLineColor(2); line7->SetLineWidth(2); line7->Draw("SAME");
	}

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}



void be_and_af2d(TH2F* hist1, TH2F* hist2, TH2F* hist3, TH2F* hist4, const char* name, const char* title, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1024, 768); c->Divide(2, 2);

	c->cd(1); hist1->Draw("Col");
	c->cd(2); hist2->Draw("Col SAME");
	c->cd(3); hist3->Draw("Col SAME");
	c->cd(4); hist4->Draw("Col SAME");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}
