
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
#include <stdio.h>
#include "../../header_files/parse_utils.h"
#include "../../header_files/hist_utils.h"
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

void Draw1d(TH1F*, const char*, const char*, TFile*);
void Draw1dE(TH1F*, const char*, const char*, TFile*);
void Draw2d(TH2F*, const char*, const char*, TFile*);
void HotWire(TH2F*, const char*, const char*, double, double, double, double, const char*, const char*, TFile*);
void OverLap2(TH1F*, TH1F*, double, double, double, double, const char*, const char*, bool, TFile*);
void Residuals(TH1F*, TH1F*, double, double, double, double, const char*, const char*, bool, TFile*);
void OverLap3(TH1F*, TH1F*, TH1F*, double, double, double, double, const char*, const char*, const char*, bool, TFile*);
void OverLap5(TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, double, double, double, double, const char*, const char*, const char*, bool, TFile*);
void OverLay3(TGraph*, TGraph*, TGraph*, const char*, const char*, const char*, bool, TFile*);
void OverLay5(TGraph*, TGraph*, TGraph*, TGraph*, TGraph*, const char*, const char*, const char*, bool, TFile*);
void Compare2(TH2F*, TH2F*, double, double, double, double, const char*, const char*, bool, TFile*);

void Draw1d_W(TH1F*, const char*, const char*, TFile*);

void compare_histos_light(TString, TString, TString, TString, TString, TString, TString, TString, TString, TString,
	TString, TString, TString,
	TString, TString, TString, TString, TString,
	bool, TFile*);

void compare_histos_heavy(TString, TString, TString, TString, TString, TString,
	TString, TString, TString,
	TString, TString, TString,
	bool, TFile*);



//*****************************************************************************************************************************************************************************
//Begin Main
//*****************************************************************************************************************************************************************************
void abslhs2()
{
	TFile* outROOT = new TFile("Hresult.root", "RECREATE");
	
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
	
	Double_t pcal_min = 0.8; Double_t pcal_max = 1.3;
	Double_t edelta_min = 0.0; Double_t edelta_max = 22;
	Double_t hdelta_min = -10; Double_t hdelta_max = 10;
	Double_t preactz_min1 = -7; Double_t preactz_max1 = 7;
	Double_t preactz_min2 = -2; Double_t preactz_max2 = 2;
	Double_t W_min = 0.9; Double_t W_max = 1.05;
	Double_t Em_min = -0.02; Double_t Em_max = 0.09;
	Double_t ctime_min = 85; Double_t ctime_max = 150;
	
	int nentries;

	Double_t Z_per_A = 1;
	Double_t transparency = 1;
	Double_t areal_den = 1;// 0.7231;
	Double_t Prot_Abs = 0.952;
	


	/////*************************************************//////
	Double_t hms_scale = 1.0;
	Double_t hms_hsize = 4.575;// 0.0275;
	Double_t hms_vsize = 11.646;// 0.07;
	Double_t shms_scale = 1.0;
	Double_t shms_hsize = 8.5;// 0.0511;
	Double_t shms_vsize = 12.5;// 0.0751;

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


	
	/////*************************************************//////
	Double_t hms_scale2 = 0.5;
	Double_t hms_hsize2 = 4.575;// 0.0275;
	Double_t hms_vsize2 = 11.646;// 0.07;
	Double_t shms_scale2 = 0.5;
	Double_t shms_hsize2 = 8.5;// 0.0511;
	Double_t shms_vsize2 = 12.5;// 0.0751;

	//Scaling the HMS/SHMS Collimator Cuts
	hms_hsize2 = hms_scale2 * hms_hsize2;
	hms_vsize2 = hms_scale2 * hms_vsize2;

	shms_hsize2 = shms_scale2 * shms_hsize2;
	shms_vsize2 = shms_scale2 * shms_vsize2;

	//Define HMS Collimator Shape
	TCutG* hms_Coll_gCut2 = new TCutG("hmsCollCut2", 8);
	hms_Coll_gCut2->SetVarX("X");
	hms_Coll_gCut2->SetVarY("Y");

	hms_Coll_gCut2->SetPoint(0, hms_hsize2, hms_vsize2 / 2.);
	hms_Coll_gCut2->SetPoint(1, hms_hsize2 / 2., hms_vsize2);
	hms_Coll_gCut2->SetPoint(2, -hms_hsize2 / 2., hms_vsize2);
	hms_Coll_gCut2->SetPoint(3, -hms_hsize2, hms_vsize2 / 2.);
	hms_Coll_gCut2->SetPoint(4, -hms_hsize2, -hms_vsize2 / 2.);
	hms_Coll_gCut2->SetPoint(5, -hms_hsize2 / 2., -hms_vsize2);
	hms_Coll_gCut2->SetPoint(6, hms_hsize2 / 2., -hms_vsize2);
	hms_Coll_gCut2->SetPoint(7, hms_hsize2, -hms_vsize2 / 2.);
	hms_Coll_gCut2->SetPoint(8, hms_hsize2, hms_vsize2 / 2.);

	//Define SHMS Collimator Shape
	TCutG* shms_Coll_gCut2 = new TCutG("shmsCollCut2", 8);
	shms_Coll_gCut2->SetVarX("X");
	shms_Coll_gCut2->SetVarY("Y");

	shms_Coll_gCut2->SetPoint(0, shms_hsize2, shms_vsize2 / 2.);
	shms_Coll_gCut2->SetPoint(1, shms_hsize2 / 2., shms_vsize2);
	shms_Coll_gCut2->SetPoint(2, -shms_hsize2 / 2., shms_vsize2);
	shms_Coll_gCut2->SetPoint(3, -shms_hsize2, shms_vsize2 / 2.);
	shms_Coll_gCut2->SetPoint(4, -shms_hsize2, -shms_vsize2 / 2.);
	shms_Coll_gCut2->SetPoint(5, -shms_hsize2 / 2., -shms_vsize2);
	shms_Coll_gCut2->SetPoint(6, shms_hsize2 / 2., -shms_vsize2);
	shms_Coll_gCut2->SetPoint(7, shms_hsize2, -shms_vsize2 / 2.);
	shms_Coll_gCut2->SetPoint(8, shms_hsize2, shms_vsize2 / 2.);







	//../../../../../../c-deuteron/cyero/worksim


	/////*************************************************//////
	//string path62_eep = "../../../../../../../../../cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS4/ROOTfiles/cafe_replay_prod_16962_-1.root";
	string path62_eep = "../../../../../../../../../../../../cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS4/ROOTfiles/cafe_replay_prod_16962_-1.root";
	//string path62_eep = "../../../../../../cafe_prod_LH2_heep_coin_16962_-1_skimmed.root";
	TFile* inROOT_62 = new TFile(path62_eep.c_str(), "READ");
	TTree* inputtree_62 = (TTree*)inROOT_62->Get("T");
	
	inputtree_62->SetBranchStatus("*", kFALSE);
	inputtree_62->SetBranchStatus("P.cal.etotnorm", kTRUE);
	inputtree_62->SetBranchStatus("CTime.epCoinTime_ROC2", kTRUE);
	inputtree_62->SetBranchStatus("P.react.z", kTRUE);
	inputtree_62->SetBranchStatus("P.gtr.th", kTRUE);
	inputtree_62->SetBranchStatus("P.gtr.ph", kTRUE);
	inputtree_62->SetBranchStatus("H.gtr.th", kTRUE);
	inputtree_62->SetBranchStatus("H.gtr.ph", kTRUE);
	inputtree_62->SetBranchStatus("P.gtr.y", kTRUE);
	inputtree_62->SetBranchStatus("P.gtr.dp", kTRUE);
	inputtree_62->SetBranchStatus("H.gtr.dp", kTRUE);
	inputtree_62->SetBranchStatus("H.kin.secondary.emiss", kTRUE);
	inputtree_62->SetBranchStatus("H.kin.secondary.pmiss", kTRUE);
	inputtree_62->SetBranchStatus("H.kin.secondary.pmiss_x", kTRUE);
	inputtree_62->SetBranchStatus("H.kin.secondary.pmiss_y", kTRUE);
	inputtree_62->SetBranchStatus("H.kin.secondary.pmiss_z", kTRUE);
	inputtree_62->SetBranchStatus("P.kin.primary.W", kTRUE);
	inputtree_62->SetBranchStatus("P.kin.primary.Q2", kTRUE);
	inputtree_62->SetBranchStatus("g.evtyp", kTRUE);
	inputtree_62->SetBranchStatus("P.extcor.xsieve", kTRUE);
	inputtree_62->SetBranchStatus("P.extcor.ysieve", kTRUE);
	inputtree_62->SetBranchStatus("H.extcor.xsieve", kTRUE);
	inputtree_62->SetBranchStatus("H.extcor.ysieve", kTRUE);
	
	Double_t Ehtrk_62 = 0.996;
	Double_t Eetrk_62 = 0.984;
	Double_t Emtrk_62 = 0.9934;
	Double_t Elt_62 = 0.941;
	Double_t charge_62 = 1;// 11.146;
	Double_t ctime_offset_62 = 1;
	Double_t weight_62;
	Double_t pcal_62; inputtree_62->SetBranchAddress("P.cal.etotnorm", &pcal_62);
	Double_t evtyp_62; inputtree_62->SetBranchAddress("g.evtyp", &evtyp_62);
	Double_t ctime_62; inputtree_62->SetBranchAddress("CTime.epCoinTime_ROC2", &ctime_62);//_center
	Double_t preactz_62; inputtree_62->SetBranchAddress("P.react.z", &preactz_62);
	Double_t exptar_62; inputtree_62->SetBranchAddress("P.gtr.th", &exptar_62);
	Double_t eyptar_62; inputtree_62->SetBranchAddress("P.gtr.ph", &eyptar_62);
	Double_t hxptar_62; inputtree_62->SetBranchAddress("H.gtr.th", &hxptar_62);
	Double_t hyptar_62; inputtree_62->SetBranchAddress("H.gtr.ph", &hyptar_62);
	Double_t eytar_62; inputtree_62->SetBranchAddress("P.gtr.y", &eytar_62);
	Double_t eXColl_62; inputtree_62->SetBranchAddress("P.extcor.xsieve", &eXColl_62);
	Double_t eYColl_62; inputtree_62->SetBranchAddress("P.extcor.ysieve", &eYColl_62);
	Double_t hXColl_62; inputtree_62->SetBranchAddress("H.extcor.xsieve", &hXColl_62);
	Double_t hYColl_62; inputtree_62->SetBranchAddress("H.extcor.ysieve", &hYColl_62);
	Double_t edelta_62; inputtree_62->SetBranchAddress("P.gtr.dp", &edelta_62);
	Double_t hdelta_62; inputtree_62->SetBranchAddress("H.gtr.dp", &hdelta_62);
	Double_t emiss_62; inputtree_62->SetBranchAddress("H.kin.secondary.emiss", &emiss_62);
	Double_t pmiss_62; inputtree_62->SetBranchAddress("H.kin.secondary.pmiss", &pmiss_62);
	Double_t pmiss_x_62; inputtree_62->SetBranchAddress("H.kin.secondary.pmiss_x", &pmiss_x_62);
	Double_t pmiss_y_62; inputtree_62->SetBranchAddress("H.kin.secondary.pmiss_y", &pmiss_y_62);
	Double_t pmiss_z_62; inputtree_62->SetBranchAddress("H.kin.secondary.pmiss_z", &pmiss_z_62);
	Double_t W_62; inputtree_62->SetBranchAddress("P.kin.primary.W", &W_62);
	Double_t Q2_62; inputtree_62->SetBranchAddress("P.kin.primary.Q2", &Q2_62);
	//Double_t xbj_62; inputtree_62->SetBranchAddress("P.kin.primary.x_bj", &xbj_62);
	//Double_t theta_e_62; inputtree_62->SetBranchAddress("exptar_62", &theta_e_62);
	//Double_t theta_p_62; inputtree_62->SetBranchAddress("exptar_62", &theta_p_62);
	


	/////*************************************************//////
	string path68_ee = "../../../../../../../../../../../../cache/hallc/c-cafe-2022/analysis/OFFLINE/PASS4/ROOTfiles/cafe_replay_prod_16968_5000000.root";
	//string path68_ee = "../../../../../../cafe_sample_LH2_heep_singles_16968_5000000_skimmed.root";
	TFile* inROOT_68 = new TFile(path68_ee.c_str(), "READ");
	TTree* inputtree_68 = (TTree*)inROOT_68->Get("T"); 
	inputtree_68->SetBranchStatus("*", kFALSE);
	inputtree_68->SetBranchStatus("P.cal.etotnorm", kTRUE);
	inputtree_68->SetBranchStatus("CTime.epCoinTime_ROC2", kTRUE);
	inputtree_68->SetBranchStatus("P.react.z", kTRUE);
	inputtree_68->SetBranchStatus("P.gtr.th", kTRUE);
	inputtree_68->SetBranchStatus("P.gtr.ph", kTRUE);
	inputtree_68->SetBranchStatus("H.gtr.th", kTRUE);
	inputtree_68->SetBranchStatus("H.gtr.ph", kTRUE);
	inputtree_68->SetBranchStatus("P.gtr.y", kTRUE);
	inputtree_68->SetBranchStatus("P.gtr.dp", kTRUE);
	inputtree_68->SetBranchStatus("H.gtr.dp", kTRUE);
	inputtree_68->SetBranchStatus("H.kin.secondary.emiss", kTRUE);
	inputtree_68->SetBranchStatus("H.kin.secondary.pmiss", kTRUE);
	inputtree_68->SetBranchStatus("H.kin.secondary.pmiss_x", kTRUE);
	inputtree_68->SetBranchStatus("H.kin.secondary.pmiss_y", kTRUE);
	inputtree_68->SetBranchStatus("H.kin.secondary.pmiss_z", kTRUE);
	inputtree_68->SetBranchStatus("P.kin.primary.W", kTRUE);
	inputtree_68->SetBranchStatus("P.kin.primary.Q2", kTRUE);
	inputtree_68->SetBranchStatus("g.evtyp", kTRUE);
	inputtree_68->SetBranchStatus("P.extcor.xsieve", kTRUE);
	inputtree_68->SetBranchStatus("P.extcor.ysieve", kTRUE);
	inputtree_68->SetBranchStatus("H.extcor.xsieve", kTRUE);
	inputtree_68->SetBranchStatus("H.extcor.ysieve", kTRUE);
	
	Double_t PS2_scalefac_68 = 9.0;
	Double_t Eetrk_68 = 0.983;
	Double_t Emtrk_68 = 0.9927;
	Double_t Elt_68 = 1.003;
	Double_t charge_68 = 1;// 15.006;
	Double_t ctime_offset_68 = 1;
	Double_t weight_68;

	Double_t evtyp_68; inputtree_68->SetBranchAddress("g.evtyp", &evtyp_68);
	Double_t pcal_68; inputtree_68->SetBranchAddress("P.cal.etotnorm", &pcal_68);
	Double_t ctime_68; inputtree_68->SetBranchAddress("CTime.epCoinTime_ROC2", &ctime_68);
	Double_t preactz_68; inputtree_68->SetBranchAddress("P.react.z", &preactz_68);
	Double_t exptar_68; inputtree_68->SetBranchAddress("P.gtr.th", &exptar_68);
	Double_t eyptar_68; inputtree_68->SetBranchAddress("P.gtr.ph", &eyptar_68);
	Double_t hxptar_68; inputtree_68->SetBranchAddress("H.gtr.th", &hxptar_68);
	Double_t hyptar_68; inputtree_68->SetBranchAddress("H.gtr.ph", &hyptar_68);
	Double_t eytar_68; inputtree_68->SetBranchAddress("P.gtr.y", &eytar_68);
	Double_t eXColl_68; inputtree_68->SetBranchAddress("P.extcor.xsieve", &eXColl_68);
	Double_t eYColl_68; inputtree_68->SetBranchAddress("P.extcor.ysieve", &eYColl_68);
	Double_t hXColl_68; inputtree_68->SetBranchAddress("H.extcor.xsieve", &hXColl_68);
	Double_t hYColl_68; inputtree_68->SetBranchAddress("H.extcor.ysieve", &hYColl_68);
	Double_t edelta_68; inputtree_68->SetBranchAddress("P.gtr.dp", &edelta_68);
	Double_t hdelta_68; inputtree_68->SetBranchAddress("H.gtr.dp", &hdelta_68);
	Double_t emiss_68; inputtree_68->SetBranchAddress("H.kin.secondary.emiss", &emiss_68);
	Double_t pmiss_68; inputtree_68->SetBranchAddress("H.kin.secondary.pmiss", &pmiss_68);
	Double_t pmiss_x_68; inputtree_68->SetBranchAddress("H.kin.secondary.pmiss_x", &pmiss_x_68);
	Double_t pmiss_y_68; inputtree_68->SetBranchAddress("H.kin.secondary.pmiss_y", &pmiss_y_68);
	Double_t pmiss_z_68; inputtree_68->SetBranchAddress("H.kin.secondary.pmiss_z", &pmiss_z_68);
	Double_t W_68; inputtree_68->SetBranchAddress("P.kin.primary.W", &W_68);
	Double_t Q2_68; inputtree_68->SetBranchAddress("P.kin.primary.Q2", &Q2_68);
	//Double_t xbj_68; inputtree_68->SetBranchAddress("P.kin.primary.x_bj", &xbj_68);
	


	/////*************************************************//////
	string pathsim_eep = "cafe_heep_coin_kin0_rad.root";
	TFile* inROOT_eep = new TFile(pathsim_eep.c_str(), "READ");
	TTree* inputtree_eep = (TTree*)inROOT_eep->Get("SNT");
	Double_t weight_eep; inputtree_eep->SetBranchAddress("Weight",  &weight_eep);
	Double_t Normfac_eep; inputtree_eep->SetBranchAddress("Normfac",  &Normfac_eep);
	Double_t nentries_eep = 1;
	Double_t FullWeight_eep;// = (Normfac_eep * weight_eep) / ((nentries_eep * areal_den * transparency * Z_per_A));

	//Double_t pcal_eep; inputtree_eep->SetBranchAddress("P.cal.etotnorm", &pcal_eep);
	//Double_t ctime_eep; inputtree_eep->SetBranchAddress("CTime.epCoinTime_ROC2_center", &ctime_eep);
	Double_t exptar_eep; inputtree_eep->SetBranchAddress("e_xptar", &exptar_eep);
	Double_t eyptar_eep; inputtree_eep->SetBranchAddress("e_yptar", &eyptar_eep);
	Double_t hxptar_eep; inputtree_eep->SetBranchAddress("h_xptar", &hxptar_eep);
	Double_t hyptar_eep; inputtree_eep->SetBranchAddress("h_yptar", &hyptar_eep);
	Double_t eytar_eep; inputtree_eep->SetBranchAddress("e_ytar", &eytar_eep);
	Double_t hytar_eep; inputtree_eep->SetBranchAddress("h_ytar", &hytar_eep);
	Double_t edelta_eep; inputtree_eep->SetBranchAddress("e_delta", &edelta_eep);
	Double_t hdelta_eep; inputtree_eep->SetBranchAddress("h_delta", &hdelta_eep);
	Double_t emiss_eep; inputtree_eep->SetBranchAddress("Em", &emiss_eep);
	Double_t pmiss_eep; inputtree_eep->SetBranchAddress("Pm", &pmiss_eep);
	Double_t pmiss_x_eep; inputtree_eep->SetBranchAddress("Pmx", &pmiss_x_eep);
	Double_t pmiss_y_eep; inputtree_eep->SetBranchAddress("Pmy", &pmiss_y_eep);
	Double_t pmiss_z_eep; inputtree_eep->SetBranchAddress("Pmz", &pmiss_z_eep);
	Double_t preactz_eep; inputtree_eep->SetBranchAddress("e_zv", &preactz_eep);
	Double_t W_eep; inputtree_eep->SetBranchAddress("W", &W_eep);
	Double_t Q2_eep; inputtree_eep->SetBranchAddress("Q2", &Q2_eep);
	//Double_t xbj_eep; inputtree_eep->SetBranchAddress("P.kin.primary.x_bj", &xbj_eep);
	Double_t eXColl_eep;// inputtree_eep->SetBranchAddress("P.extcor.xsieve", &eXColl_eep);
	Double_t eYColl_eep;// inputtree_eep->SetBranchAddress("P.extcor.ysieve", &eYColl_eep);
	Double_t hXColl_eep;// inputtree_eep->SetBranchAddress("H.extcor.xsieve", &hXColl_eep);
	Double_t hYColl_eep;// inputtree_eep->SetBranchAddress("H.extcor.ysieve", &hYColl_eep);
	Double_t h_angle_eep = 48.3;//(stod(split(split(FindString("spec%p%theta", infile.Data())[0], '!')[0], '=')[1]));
	Double_t e_angle_eep = 8.3;//(stod(split(split(FindString("spec%e%theta", infile.Data())[0], '!')[0], '=')[1]));
	Double_t hreactz_eep; inputtree_eep->SetBranchAddress("h_zv", &hreactz_eep);
	Double_t tarx_eep; inputtree_eep->SetBranchAddress("tar_x", &tarx_eep);
	Double_t htarx_corr_eep;// = tarx_eep - hxptar_eep * hreactz_eep * cos(h_angle_eep * dtr);
	Double_t etarx_corr_eep;// = tarx_eep - exptar_eep * preactz_eep * cos(e_angle_eep * dtr);
	
	

	/////*************************************************//////
	string pathsim_ee = "cafe_heep_singles_kin0_rad.root";
	TFile* inROOT_ee = new TFile(pathsim_ee.c_str(), "READ");
	TTree* inputtree_ee = (TTree*)inROOT_ee->Get("SNT");
	Double_t weight_ee; inputtree_ee->SetBranchAddress("Weight", &weight_ee);
	Double_t Normfac_ee; inputtree_ee->SetBranchAddress("Normfac", &Normfac_ee);
	Double_t nentries_ee = 1;
	Double_t FullWeight_ee;// = (Normfac_ee * weight_ee) / ((nentries_ee * areal_den * transparency * Z_per_A));
	//cout << "FullWeight ee: " << FullWeight_ee << endl;

	//Double_t pcal_ee; inputtree_ee->SetBranchAddress("P.cal.etotnorm", &pcal_ee);
	//Double_t ctime_ee; inputtree_ee->SetBranchAddress("CTime.epCoinTime_ROC2_center", &ctime_ee);
	Double_t exptar_ee; inputtree_ee->SetBranchAddress("e_xptar", &exptar_ee);
	Double_t eyptar_ee; inputtree_ee->SetBranchAddress("e_yptar", &eyptar_ee);
	Double_t hxptar_ee; inputtree_ee->SetBranchAddress("h_xptar", &hxptar_ee);
	Double_t hyptar_ee; inputtree_ee->SetBranchAddress("h_yptar", &hyptar_ee);
	Double_t eXColl_ee;// inputtree_ee->SetBranchAddress("P.extcor.xsieve", &eXColl_ee);
	Double_t eYColl_ee;// inputtree_ee->SetBranchAddress("P.extcor.ysieve", &eYColl_ee);
	Double_t hXColl_ee;// inputtree_ee->SetBranchAddress("H.extcor.xsieve", &hXColl_ee);
	Double_t hYColl_ee;// inputtree_ee->SetBranchAddress("H.extcor.ysieve", &hYColl_ee);
	Double_t eytar_ee; inputtree_ee->SetBranchAddress("e_ytar", &eytar_ee);
	Double_t hytar_ee; inputtree_ee->SetBranchAddress("h_ytar", &hytar_ee);
	Double_t edelta_ee; inputtree_ee->SetBranchAddress("e_delta", &edelta_ee);
	//Dobuel_t hdelta_ee; inputtree_ee->SetBranchAddress("h_delta", &hdelta_ee);
	//Double_t emiss_ee; inputtree_ee->SetBranchAddress("Em", &emiss_ee);
	//Double_t pmiss_ee; inputtree_ee->SetBranchAddress("Pm", &pmiss_ee);
	//Double_t pmiss_x_ee; inputtree_ee->SetBranchAddress("Pmx", &pmiss_x_ee);
	//Double_t pmiss_y_ee; inputtree_ee->SetBranchAddress("Pmy", &pmiss_y_ee);
	//Double_t pmiss_z_ee; inputtree_ee->SetBranchAddress("Pmz", &pmiss_z_ee);
	Double_t preactz_ee; inputtree_ee->SetBranchAddress("e_zv", &preactz_ee);
	Double_t W_ee; inputtree_ee->SetBranchAddress("W", &W_ee);
	Double_t Q2_ee; inputtree_ee->SetBranchAddress("Q2", &Q2_ee);
	//Double_t xbj_ee; inputtree_ee->SetBranchAddress("P.kin.primary.x_bj", &xbj_ee);
	Double_t h_angle_ee = 48.3;//(stod(split(split(FindString("spec%p%theta", infile.Data())[0], '!')[0], '=')[1]));
	Double_t hreactz_ee; inputtree_ee->SetBranchAddress("h_zv", &hreactz_ee);
	Double_t tarx_ee; inputtree_ee->SetBranchAddress("tar_x", &tarx_ee);
	Double_t e_angle_ee = 8.295;//(stod(split(split(FindString("spec%e%theta", infile.Data())[0], '!')[0], '=')[1]));
	Double_t htarx_corr_ee;// = tarx_ee - hxptar_ee * hreactz_ee * cos(h_angle_ee * dtr);
	Double_t etarx_corr_ee;// = tarx_ee - exptar_ee * preactz_ee * cos(e_angle_ee * dtr);
	


	


	TList* List = new TList();		//Create TList to store histograms
	TList* List1 = new TList();		//Create TList to store histograms


	//**** pcaletotnorm ****//
	TH1F* H1_pcal_62 = new TH1F("H1_pcal_62", "SHMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}", 50, 0.0, 2.0);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pcal_62_1 = new TH1F("H1_pcal_62_1", "SHMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}", 50, 0.0, 2.0);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pcal_62);
	List1->Add(H1_pcal_62_1);

	TH1F* H1_pcal_68 = new TH1F("H1_pcal_68", "SHMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}", 50, 0.0, 2.0);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pcal_68_1 = new TH1F("H1_pcal_68_1", "SHMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}", 50, 0.0, 2.0);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pcal_68);
	List1->Add(H1_pcal_68_1);
	
	//**** ctime ****//
	TH1F* H1_coin_62 = new TH1F("H1_coin_62", "Coincidence Time; coin (ns)", 100, 80, 150);							//Create a 1D Histogram of Coincidence Time from run 16968
	TH1F* H1_coin_62_1 = new TH1F("H1_coin_62_1", "Coincidence Time; coin (ns)", 100, 80, 150);			//Create a 1D Histogram of Coincidence Time from run 16968
	List->Add(H1_coin_62);
	List1->Add(H1_coin_62_1);

	TH1F* H1_coin_68 = new TH1F("H1_coin_68", "Coincidence Time; coin (ns)", 100, -5, 200);							//Create a 1D Histogram of Coincidence Time from run 16968
	TH1F* H1_coin_68_1 = new TH1F("H1_coin_68_1", "Coincidence Time; coin (ns)", 100, -5, 200);			//Create a 1D Histogram of Coincidence Time from run 16968
	List->Add(H1_coin_68);
	List1->Add(H1_coin_68_1);

		

	//**** preactz ****//
	TH1F* H1_preactz_62 = new TH1F("H1_preactz_62", "SHMS z-Target (Lab); z-Target (cm)", 50, -10, 10);					//Create a 1D Histogram of SHMS z-Target (Lab) from run 16968
	TH1F* H1_preactz_62_1 = new TH1F("H1_preactz_62_1", "SHMS z-Target (Lab); z-Target (cm)", 50, -10, 10);	//Create a 1D Histogram of SHMS z-Target (Lab) from run 16968
	List->Add(H1_preactz_62);
	List1->Add(H1_preactz_62_1);

	TH1F* H1_preactz_68 = new TH1F("H1_preactz_68", "SHMS z-Target (Lab); z-Target (cm)", 50, -10, 10);					//Create a 1D Histogram of SHMS z-Target (Lab) from run 16968
	TH1F* H1_preactz_68_1 = new TH1F("H1_preactz_68_1", "SHMS z-Target (Lab); z-Target (cm)", 50, -10, 10);	//Create a 1D Histogram of SHMS z-Target (Lab) from run 16968
	List->Add(H1_preactz_68);
	List1->Add(H1_preactz_68_1);

	TH1F* H1_preactz_ee = new TH1F("H1_preactz_ee", "SHMS z-Target (Lab); z-Target (cm)", 50, -10, 10);					//Create a 1D Histogram of SHMS z-Target (Lab) from run 16968
	TH1F* H1_preactz_ee_1 = new TH1F("H1_preactz_ee_1", "SHMS z-Target (Lab); z-Target (cm)", 50, -10, 10);	//Create a 1D Histogram of SHMS z-Target (Lab) from run 16968
	List->Add(H1_preactz_ee);
	List1->Add(H1_preactz_ee_1);
	
	TH1F* H1_preactz_eep = new TH1F("H1_preactz_eep", "SHMS z-Target (Lab); z-Target (cm)", 50, -10, 10);					//Create a 1D Histogram of SHMS z-Target (Lab) from run 16968
	TH1F* H1_preactz_eep_1 = new TH1F("H1_preactz_eep_1", "SHMS z-Target (Lab); z-Target (cm)", 50, -10, 10);	//Create a 1D Histogram of SHMS z-Target (Lab) from run 16968
	List->Add(H1_preactz_eep);
	List1->Add(H1_preactz_eep_1);
	


	//**** SHMS X'_{tar}  ****//
	TH1F* H1_exptar_62 = new TH1F("H1_exptar_62", "SHMS X'_{tar}; X'_{tar} (rad)", 200, -0.05, 0.05);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_exptar_62_1 = new TH1F("H1_exptar_62_1", "SHMS X'_{tar}; X'_{tar} (rad)", 200, -0.05, 0.05);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_exptar_62);
	List1->Add(H1_exptar_62_1);

	TH1F* H1_exptar_68 = new TH1F("H1_exptar_68", "SHMS X'_{tar}; X'_{tar} (rad)", 200, -0.05, 0.05);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_exptar_68_1 = new TH1F("H1_exptar_68_1", "SHMS X'_{tar}; X'_{tar} (rad)", 200, -0.05, 0.05);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_exptar_68);
	List1->Add(H1_exptar_68_1);

	TH1F* H1_exptar_eep = new TH1F("H1_exptar_eep", "SHMS X'_{tar}; X'_{tar} (rad)", 200, -0.05, 0.05);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_exptar_eep_1 = new TH1F("H1_exptar_eep_1", "SHMS X'_{tar}; X'_{tar} (rad)", 200, -0.05, 0.05);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_exptar_eep);
	List1->Add(H1_exptar_eep_1);

	TH1F* H1_exptar_ee = new TH1F("H1_exptar_ee", "SHMS X'_{tar}; X'_{tar} (rad)", 200, -0.05, 0.05);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_exptar_ee_1 = new TH1F("H1_exptar_ee_1", "SHMS X'_{tar}; X'_{tar} (rad)", 200, -0.05, 0.05);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_exptar_ee);
	List1->Add(H1_exptar_ee_1);

	

	//**** SHMS Y'_{tar}  ****//
	TH1F* H1_eyptar_62 = new TH1F("H1_eyptar_62", "SHMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.05, 0.05);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_eyptar_62_1 = new TH1F("H1_eyptar_62_1", "SHMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.05, 0.05);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_eyptar_62);
	List1->Add(H1_eyptar_62_1);

	TH1F* H1_eyptar_68 = new TH1F("H1_eyptar_68", "SHMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.05, 0.05);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_eyptar_68_1 = new TH1F("H1_eyptar_68_1", "SHMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.05, 0.05);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_eyptar_68);
	List1->Add(H1_eyptar_68_1);

	TH1F* H1_eyptar_eep = new TH1F("H1_eyptar_eep", "SHMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.05, 0.05);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_eyptar_eep_1 = new TH1F("H1_eyptar_eep_1", "SHMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.05, 0.05);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_eyptar_eep);
	List1->Add(H1_eyptar_eep_1);

	TH1F* H1_eyptar_ee = new TH1F("H1_eyptar_ee", "SHMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.05, 0.05);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_eyptar_ee_1 = new TH1F("H1_eyptar_ee_1", "SHMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.05, 0.05);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_eyptar_ee);
	List1->Add(H1_eyptar_ee_1);


	


	//**** HMS X'_{tar}  ****//
	TH1F* H1_hxptar_62 = new TH1F("H1_hxptar_62", "HMS X'_{tar}; X'_{tar} (rad)", 200, -0.1, 0.1);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_hxptar_62_1 = new TH1F("H1_hxptar_62_1", "HMS X'_{tar}; X'_{tar} (rad)", 200, -0.1, 0.1);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_hxptar_62);
	List1->Add(H1_hxptar_62_1);

	TH1F* H1_hxptar_68 = new TH1F("H1_hxptar_68", "HMS X'_{tar}; X'_{tar} (rad)", 200, -0.1, 0.1);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_hxptar_68_1 = new TH1F("H1_hxptar_68_1", "HMS X'_{tar}; X'_{tar} (rad)", 200, -0.1, 0.1);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_hxptar_68);
	List1->Add(H1_hxptar_68_1);

	TH1F* H1_hxptar_eep = new TH1F("H1_hxptar_eep", "HMS X'_{tar}; X'_{tar} (rad)", 200, -0.1, 0.1);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_hxptar_eep_1 = new TH1F("H1_hxptar_eep_1", "HMS X'_{tar}; X'_{tar} (rad)", 200, -0.1, 0.1);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_hxptar_eep);
	List1->Add(H1_hxptar_eep_1);
	

	//**** HMS Y'_{tar}  ****//
	TH1F* H1_hyptar_62 = new TH1F("H1_hyptar_62", "HMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.1, 0.1);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_hyptar_62_1 = new TH1F("H1_hyptar_62_1", "HMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.1, 0.1);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_hyptar_62);
	List1->Add(H1_hyptar_62_1);

	TH1F* H1_hyptar_68 = new TH1F("H1_hyptar_68", "HMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.1, 0.1);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_hyptar_68_1 = new TH1F("H1_hyptar_68_1", "HMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.1, 0.1);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_hyptar_68);
	List1->Add(H1_hyptar_68_1);

	TH1F* H1_hyptar_eep = new TH1F("H1_hyptar_eep", "HMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.1, 0.1);					//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	TH1F* H1_hyptar_eep_1 = new TH1F("H1_hyptar_eep_1", "HMS Y'_{tar}; Y'_{tar} (rad)", 200, -0.1, 0.1);		//Create a 1D Histogram of SHMS X'_{tar} from run 16968
	List->Add(H1_hyptar_eep);
	List1->Add(H1_hyptar_eep_1);
	

	//**** SHMS Momentum Acceptance  ****//
	TH1F* H1_edelta_62 = new TH1F("H1_edelta_62", "SHMS Momentum Acceptance; #sigma (%)", 50, -1, 23);				//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	TH1F* H1_edelta_62_1 = new TH1F("H1_edelta_62_1", "SHMS Momentum Acceptance; #sigma (%)", 50, -1, 23);	//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	List->Add(H1_edelta_62);
	List1->Add(H1_edelta_62_1);

	TH1F* H1_edelta_68 = new TH1F("H1_edelta_68", "SHMS Momentum Acceptance; #sigma (%)", 50, -1, 23);				//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	TH1F* H1_edelta_68_1 = new TH1F("H1_edelta_68_1", "SHMS Momentum Acceptance; #sigma (%)", 50, -1, 23);	//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	List->Add(H1_edelta_68);
	List1->Add(H1_edelta_68_1);

	TH1F* H1_edelta_eep = new TH1F("H1_edelta_eep", "SHMS Momentum Acceptance; #sigma (%)", 50, -1, 23);				//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	TH1F* H1_edelta_eep_1 = new TH1F("H1_edelta_eep_1", "SHMS Momentum Acceptance; #sigma (%)", 50, -1, 23);	//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	List->Add(H1_edelta_eep);
	List1->Add(H1_edelta_eep_1);

	TH1F* H1_edelta_ee = new TH1F("H1_edelta_ee", "SHMS Momentum Acceptance; #sigma (%)", 50, -1, 23);				//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	TH1F* H1_edelta_ee_1 = new TH1F("H1_edelta_ee_1", "SHMS Momentum Acceptance; #sigma (%)", 50, -1, 23);	//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	List->Add(H1_edelta_ee);
	List1->Add(H1_edelta_ee_1);



	//**** HMS Momentum Acceptance  ****//
	TH1F* H1_hdelta_62 = new TH1F("H1_hdelta_62", "HMS Momentum Acceptance; #sigma (%)", 50, -12, 12);				//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	TH1F* H1_hdelta_62_1 = new TH1F("H1_hdelta_62_1", "HMS Momentum Acceptance; #sigma (%)", 50, -12, 12);	//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	List->Add(H1_hdelta_62);
	List1->Add(H1_hdelta_62_1);

	TH1F* H1_hdelta_68 = new TH1F("H1_hdelta_68", "HMS Momentum Acceptance; #sigma (%)", 50, -12, 12);				//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	TH1F* H1_hdelta_68_1 = new TH1F("H1_hdelta_68_1", "HMS Momentum Acceptance; #sigma (%)", 50, -12, 12);	//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	List->Add(H1_hdelta_68);
	List1->Add(H1_hdelta_68_1);

	TH1F* H1_hdelta_eep = new TH1F("H1_hdelta_eep", "HMS Momentum Acceptance; #sigma (%)", 50, -12, 12);				//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	TH1F* H1_hdelta_eep_1 = new TH1F("H1_hdelta_eep_1", "HMS Momentum Acceptance; #sigma (%)", 50, -12, 12);	//Create a 1D Histogram of SHMS Momentum Acceptance from run 16968
	List->Add(H1_hdelta_eep);
	List1->Add(H1_hdelta_eep_1);
	

	//Missing Energy (GeV)
	TH1F* H1_emiss_62 = new TH1F("H1_emiss_62", "Missing Energy; emiss (GeV)", 50, -0.03, 0.1);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_emiss_62_1 = new TH1F("H1_emiss_62_1", "Missing Energy; emiss (GeV)", 50, -0.03, 0.1);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_emiss_62);
	List1->Add(H1_emiss_62_1);

	TH1F* H1_emiss_68 = new TH1F("H1_emiss_68", "Missing Energy; emiss (GeV)", 50, -0.03, 0.1);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_emiss_68_1 = new TH1F("H1_emiss_68_1", "Missing Energy; emiss (GeV)", 50, -0.03, 0.1);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_emiss_68);
	List1->Add(H1_emiss_68_1);

	TH1F* H1_emiss_eep = new TH1F("H1_emiss_eep", "Missing Energy; emiss (GeV)", 50, -0.03, 0.1);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_emiss_eep_1 = new TH1F("H1_emiss_eep_1", "Missing Energy; emiss (GeV)", 50, -0.03, 0.1);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_emiss_eep);
	List1->Add(H1_emiss_eep_1);
	

	//Missing Energy (GeV)
	TH1F* H1_pmiss_62 = new TH1F("H1_pmiss_62", "Missing Momentum; P_{miss} (GeV)", 50, -0.05, 0.2);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_62_1 = new TH1F("H1_pmiss_62_1", "Missing Momentum; P_{miss} (GeV)", 50, -0.05, 0.2);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_62);
	List1->Add(H1_pmiss_62_1);

	TH1F* H1_pmiss_68 = new TH1F("H1_pmiss_68", "Missing Momentum; P_{miss} (GeV)", 50, -0.05, 0.2);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_68_1 = new TH1F("H1_pmiss_68_1", "Missing Momentum; P_{miss} (GeV)", 50, -0.05, 0.2);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_68);
	List1->Add(H1_pmiss_68_1);

	TH1F* H1_pmiss_eep = new TH1F("H1_pmiss_eep", "Missing Momentum; P_{miss} (GeV)", 50, -0.05, 0.2);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_eep_1 = new TH1F("H1_pmiss_eep_1", "Missing Momentum; P_{miss} (GeV)", 50, -0.05, 0.2);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_eep);
	List1->Add(H1_pmiss_eep_1);
	

	//Missing Energy (GeV)
	TH1F* H1_pmiss_x_62 = new TH1F("H1_pmiss_x_62", "P_{miss, x}; P_{miss, x} (GeV)", 50, -0.5, 0.5);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_x_62_1 = new TH1F("H1_pmiss_x_62_1", "P_{miss, x}; P_{miss, x} (GeV)", 50, -0.5, 0.5);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_x_62);
	List1->Add(H1_pmiss_x_62_1);

	TH1F* H1_pmiss_x_68 = new TH1F("H1_pmiss_x_68", "P_{miss, x}; P_{miss, x} (GeV)", 50, -0.5, 0.5);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_x_68_1 = new TH1F("H1_pmiss_x_68_1", "P_{miss, x}; P_{miss, x} (GeV)", 50, -0.5, 0.5);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_x_68);
	List1->Add(H1_pmiss_x_68_1);

	TH1F* H1_pmiss_x_eep = new TH1F("H1_pmiss_x_eep", "P_{miss, x}; P_{miss, x} (GeV)", 50, -0.5, 0.5);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_x_eep_1 = new TH1F("H1_pmiss_x_eep_1", "P_{miss, x}; P_{miss, x} (GeV)", 50, -0.5, 0.5);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_x_eep);
	List1->Add(H1_pmiss_x_eep_1);
	

	//Missing Energy (GeV)
	TH1F* H1_pmiss_y_62 = new TH1F("H1_pmiss_y_62", "P_{miss, y}; P_{miss, y} (GeV)", 50, -0.5, 0.5);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_y_62_1 = new TH1F("H1_pmiss_y_62_1", "P_{miss, y}; P_{miss, y} (GeV)", 50, -0.5, 0.5);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_y_62);
	List1->Add(H1_pmiss_y_62_1);

	TH1F* H1_pmiss_y_68 = new TH1F("H1_pmiss_y_68", "P_{miss, y}; P_{miss, y} (GeV)", 50, -0.5, 0.5);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_y_68_1 = new TH1F("H1_pmiss_y_68_1", "P_{miss, y}; P_{miss, y} (GeV)", 50, -0.5, 0.5);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_y_68);
	List1->Add(H1_pmiss_y_68_1);

	TH1F* H1_pmiss_y_eep = new TH1F("H1_pmiss_y_eep", "P_{miss, y}; P_{miss, y} (GeV)", 50, -0.5, 0.5);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_y_eep_1 = new TH1F("H1_pmiss_y_eep_1", "P_{miss, y}; P_{miss, y} (GeV)", 50, -0.5, 0.5);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_y_eep);
	List1->Add(H1_pmiss_y_eep_1);
	





	//Missing Energy (GeV)
	TH1F* H1_pmiss_z_62 = new TH1F("H1_pmiss_z_62", "P_{miss, z}; P_{miss, z} (GeV)", 50, -0.5, 0.5);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_z_62_1 = new TH1F("H1_pmiss_z_62_1", "P_{miss, z}; P_{miss, z} (GeV)", 50, -0.5, 0.5);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_z_62);
	List1->Add(H1_pmiss_z_62_1);


	TH1F* H1_pmiss_z_68 = new TH1F("H1_pmiss_z_68", "P_{miss, z}; P_{miss, z} (GeV)", 50, -0.5, 0.5);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_z_68_1 = new TH1F("H1_pmiss_z_68_1", "P_{miss, z}; P_{miss, z} (GeV)", 50, -0.5, 0.5);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_z_68);
	List1->Add(H1_pmiss_z_68_1);


	TH1F* H1_pmiss_z_eep = new TH1F("H1_pmiss_z_eep", "P_{miss, z}; P_{miss, z} (GeV)", 50, -0.5, 0.5);						//Create a 1D Histogram of Missing Energy from run 16968
	TH1F* H1_pmiss_z_eep_1 = new TH1F("H1_pmiss_z_eep_1", "P_{miss, z}; P_{miss, z} (GeV)", 50, -0.5, 0.5);		//Create a 1D Histogram of Missing Energy from run 16968
	List->Add(H1_pmiss_z_eep);
	List1->Add(H1_pmiss_z_eep_1);
	



	//Invariant Mass (GeV)
	TH1F* H1_W_62 = new TH1F("H1_W_62", "Invariant Mass; W (GeV)", 50, 0.85, 1.1);										//Create a 1D Histogram of Invariant Mass from run 16968
	TH1F* H1_W_62_1 = new TH1F("H1_W_62_1", "Invariant Mass; W (GeV)", 50, 0.85, 1.1);						//Create a 1D Histogram of Invariant Mass from run 16968
	List->Add(H1_W_62);
	List1->Add(H1_W_62_1);

	TH1F* H1_W_68 = new TH1F("H1_W_68", "Invariant Mass; W (GeV)", 50, 0.85, 1.1);										//Create a 1D Histogram of Invariant Mass from run 16968
	TH1F* H1_W_68_1 = new TH1F("H1_W_68_1", "Invariant Mass; W (GeV)", 50, 0.85, 1.1);						//Create a 1D Histogram of Invariant Mass from run 16968
	List->Add(H1_W_68);
	List1->Add(H1_W_68_1);

	const Int_t NBins_W = 25;//40
	Double_t edges_W[NBins_W + 1] = {
		0.900, 0.906, 0.912, 0.918, 0.924,
		0.93, 0.936, 0.942, 0.948, 0.954, 
		0.96, 0.966, 0.972, 0.978, 0.984,
		0.99, 0.996, 1.002, 1.008, 1.014,
		1.02, 1.026, 1.032, 1.038, 1.044,
		1.05//26
	};//W
	TH1F* H1_W_ratio_pid_acc_kin_data_eep = new TH1F("H1_W_ratio_pid_acc_kin_data_eep", "W_{bj}; W_{bj} (GeV)", NBins_W, edges_W); List->Add(H1_W_ratio_pid_acc_kin_data_eep);
	TH1F* H1_W_ratio_pid_acc_kin_sim_eep = new TH1F("H1_W_ratio_pid_acc_kin_sim_eep", "W_{bj}; W_{bj} (GeV)", NBins_W, edges_W); List->Add(H1_W_ratio_pid_acc_kin_sim_eep);
	TH1F* H1_W_ratio_pid_acc_kin_data_ee = new TH1F("H1_W_ratio_pid_acc_kin_data_ee", "W_{bj}; W_{bj} (GeV)", NBins_W, edges_W); List->Add(H1_W_ratio_pid_acc_kin_data_ee);
	TH1F* H1_W_ratio_pid_acc_kin_sim_ee = new TH1F("H1_W_ratio_pid_acc_kin_sim_ee", "W_{bj}; W_{bj} (GeV)", NBins_W, edges_W); List->Add(H1_W_ratio_pid_acc_kin_sim_ee);

	TH1F* H1_W_ee = new TH1F("H1_W_ee", "Invariant Mass; W (GeV)", 50, 0.85, 1.1);										//Create a 1D Histogram of Invariant Mass from run 16968
	TH1F* H1_W_ee_1 = new TH1F("H1_W_ee_1", "Invariant Mass; W (GeV)", 50, 0.85, 1.1);						//Create a 1D Histogram of Invariant Mass from run 16968
	List->Add(H1_W_ee);
	List1->Add(H1_W_ee_1);

	TH1F* H1_W_eep = new TH1F("H1_W_eep", "Invariant Mass; W (GeV)", 50, 0.85, 1.1);										//Create a 1D Histogram of Invariant Mass from run 16968
	TH1F* H1_W_eep_1 = new TH1F("H1_W_eep_1", "Invariant Mass; W (GeV)", 50, 0.85, 1.1);						//Create a 1D Histogram of Invariant Mass from run 16968
	List->Add(H1_W_eep);
	List1->Add(H1_W_eep_1);
	

	
	//4-Momentum Transfer (GeV/c)^2
	TH1F* H1_Q2_62 = new TH1F("H1_Q2_62", "4-Momentum Transfer; Q^{2} (Gev/c)^{2}", 50, 1, 3);										//Create a 1D Histogram of Invariant Mass from run 16968
	TH1F* H1_Q2_62_1 = new TH1F("H1_Q2_62_1", "4-Momentum Transfer; Q^{2} (Gev/c)^{2}", 50, 1, 3);						//Create a 1D Histogram of Invariant Mass from run 16968
	List->Add(H1_Q2_62);
	List1->Add(H1_Q2_62_1);

	TH1F* H1_Q2_68 = new TH1F("H1_Q2_68", "4-Momentum Transfer; Q^{2} (Gev/c)^{2}", 50, 1, 3);										//Create a 1D Histogram of Invariant Mass from run 16968
	TH1F* H1_Q2_68_1 = new TH1F("H1_Q2_68_1", "4-Momentum Transfer; Q^{2} (Gev/c)^{2}", 50, 1, 3);						//Create a 1D Histogram of Invariant Mass from run 16968
	List->Add(H1_Q2_68);
	List1->Add(H1_Q2_68_1);

	TH1F* H1_Q2_ee = new TH1F("H1_Q2_ee", "4-Momentum Transfer; Q^{2} (Gev/c)^{2}", 50, 1, 3);										//Create a 1D Histogram of Invariant Mass from run 16968
	TH1F* H1_Q2_ee_1 = new TH1F("H1_Q2_ee_1", "4-Momentum Transfer; Q^{2} (Gev/c)^{2}", 50, 1, 3);						//Create a 1D Histogram of Invariant Mass from run 16968
	List->Add(H1_Q2_ee);
	List1->Add(H1_Q2_ee_1);

	TH1F* H1_Q2_eep = new TH1F("H1_Q2_eep", "4-Momentum Transfer; Q^{2} (Gev/c)^{2}", 50, 1, 3);										//Create a 1D Histogram of Invariant Mass from run 16968
	TH1F* H1_Q2_eep_1 = new TH1F("H1_Q2_eep_1", "4-Momentum Transfer; Q^{2} (Gev/c)^{2}", 50, 1, 3);						//Create a 1D Histogram of Invariant Mass from run 16968
	List->Add(H1_Q2_eep);
	List1->Add(H1_Q2_eep_1);
	

	//SHMS Momentum Acceptance vs HMS Momentum Acceptance
	TH2F* H2_hdelta_edelta_62 = new TH2F("H2_hdelta_edelta_62", "SHMS vs HMS; #sigma_{h} (%); #sigma_{e} (%)", 50, -12, 12, 50, 7, 13);					//Create a 2D Histogram of HMS Momentum Acceptance vs SHMS Momentum Acceptance from run 16968
	TH2F* H2_hdelta_edelta_62_1 = new TH2F("H2_hdelta_edelta_62_1", "SHMS vs HMS; #sigma_{h} (%); #sigma_{e} (%)", 50, -12, 12, 50, 7, 13);					//Create a 2D Histogram of HMS Momentum Acceptance vs SHMS Momentum Acceptance from run 16968
	List->Add(H2_hdelta_edelta_62);
	List1->Add(H2_hdelta_edelta_62_1);

	TH2F* H2_hdelta_edelta_68 = new TH2F("H2_hdelta_edelta_68", "SHMS vs HMS; #sigma_{h} (%); #sigma_{e} (%)", 50, -12, 12, 50, 7, 13);					//Create a 2D Histogram of HMS Momentum Acceptance vs SHMS Momentum Acceptance from run 16968
	TH2F* H2_hdelta_edelta_68_1 = new TH2F("H2_hdelta_edelta_68_1", "SHMS vs HMS; #sigma_{h} (%); #sigma_{e} (%)", 50, -12, 12, 50, 7, 13);					//Create a 2D Histogram of HMS Momentum Acceptance vs SHMS Momentum Acceptance from run 16968
	List->Add(H2_hdelta_edelta_68);
	List1->Add(H2_hdelta_edelta_68_1);

	TH2F* H2_hdelta_edelta_eep = new TH2F("H2_hdelta_edelta_eep", "SHMS vs HMS; #sigma_{h} (%); #sigma_{e} (%)", 50, -12, 12, 50, 7, 13);					//Create a 2D Histogram of HMS Momentum Acceptance vs SHMS Momentum Acceptance from run 16968
	TH2F* H2_hdelta_edelta_eep_1 = new TH2F("H2_hdelta_edelta_eep_1", "SHMS vs HMS; #sigma_{h} (%); #sigma_{e} (%)", 50, -12, 12, 50, 7, 13);					//Create a 2D Histogram of HMS Momentum Acceptance vs SHMS Momentum Acceptance from run 16968
	List->Add(H2_hdelta_edelta_eep);
	List1->Add(H2_hdelta_edelta_eep_1);
	



	//SHMS X'_{tar} vs SHMS Y'_{tar}
	TH2F* H2_eyptar_exptar_62 = new TH2F("H2_eyptar_exptar_62", "SHMS X'_{tar} vs SHMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.05, 0.05, 200, -0.05, 0.05);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	TH2F* H2_eyptar_exptar_62_1 = new TH2F("H2_eyptar_exptar_62_1", "SHMS X'_{tar} vs SHMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.05, 0.05, 200, -0.05, 0.05);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	List->Add(H2_eyptar_exptar_62);
	List1->Add(H2_eyptar_exptar_62_1);

	TH2F* H2_eyptar_exptar_68 = new TH2F("H2_eyptar_exptar_68", "SHMS X'_{tar} vs SHMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.05, 0.05, 200, -0.05, 0.05);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	TH2F* H2_eyptar_exptar_68_1 = new TH2F("H2_eyptar_exptar_68_1", "SHMS X'_{tar} vs SHMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.05, 0.05, 200, -0.05, 0.05);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	List->Add(H2_eyptar_exptar_68);
	List1->Add(H2_eyptar_exptar_68_1);

	TH2F* H2_eyptar_exptar_eep = new TH2F("H2_eyptar_exptar_eep", "SHMS X'_{tar} vs SHMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.05, 0.05, 200, -0.05, 0.05);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	TH2F* H2_eyptar_exptar_eep_1 = new TH2F("H2_eyptar_exptar_eep_1", "SHMS X'_{tar} vs SHMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.05, 0.05, 200, -0.05, 0.05);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	List->Add(H2_eyptar_exptar_eep);
	List1->Add(H2_eyptar_exptar_eep_1);

	TH2F* H2_eyptar_exptar_ee = new TH2F("H2_eyptar_exptar_ee", "SHMS X'_{tar} vs SHMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.05, 0.05, 200, -0.05, 0.05);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	TH2F* H2_eyptar_exptar_ee_1 = new TH2F("H2_eyptar_exptar_ee_1", "SHMS X'_{tar} vs SHMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.05, 0.05, 200, -0.05, 0.05);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	List->Add(H2_eyptar_exptar_ee);
	List1->Add(H2_eyptar_exptar_ee_1);





	TH2F* H2_eYColl_eXColl_62 = new TH2F("H2_eYColl_eXColl_62", "SHMS Collimator; SHMS Y-Collimator; SHMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	TH2F* H2_eYColl_eXColl_62_1 = new TH2F("H2_eYColl_eXColl_62_1", "SHMS Collimator; SHMS Y-Collimator; SHMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	List->Add(H2_eYColl_eXColl_62);
	List1->Add(H2_eYColl_eXColl_62_1);

	TH2F* H2_eYColl_eXColl_68 = new TH2F("H2_eYColl_eXColl_68", "SHMS Collimator; SHMS Y-Collimator; SHMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	TH2F* H2_eYColl_eXColl_68_1 = new TH2F("H2_eYColl_eXColl_68_1", "SHMS Collimator; SHMS Y-Collimator; SHMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	List->Add(H2_eYColl_eXColl_68);
	List1->Add(H2_eYColl_eXColl_68_1);

	TH2F* H2_eYColl_eXColl_eep = new TH2F("H2_eYColl_eXColl_eep", "SHMS Collimator; SHMS Y-Collimator; SHMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	TH2F* H2_eYColl_eXColl_eep_1 = new TH2F("H2_eYColl_eXColl_eep_1", "SHMS Collimator; SHMS Y-Collimator; SHMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	List->Add(H2_eYColl_eXColl_eep);
	List1->Add(H2_eYColl_eXColl_eep_1);

	TH2F* H2_eYColl_eXColl_ee = new TH2F("H2_eYColl_eXColl_ee", "SHMS Collimator; SHMS Y-Collimator; SHMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	TH2F* H2_eYColl_eXColl_ee_1 = new TH2F("H2_eYColl_eXColl_ee_1", "SHMS Collimator; SHMS Y-Collimator; SHMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of SHMS X'_{tar} vs SHMS Y'_{tar} from run 16968
	List->Add(H2_eYColl_eXColl_ee);
	List1->Add(H2_eYColl_eXColl_ee_1);



	//HMS Momentum Acceptance
	TH2F* H2_hyptar_hxptar_62 = new TH2F("H2_hyptar_hxptar_62", "HMS X'_{tar} vs HMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.1, 0.1, 200, -0.1, 0.1);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	TH2F* H2_hyptar_hxptar_62_1 = new TH2F("H2_hyptar_hxptar_62_1", "HMS X'_{tar} vs HMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.1, 0.1, 200, -0.1, 0.1);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	List->Add(H2_hyptar_hxptar_62);
	List1->Add(H2_hyptar_hxptar_62_1);

	TH2F* H2_hyptar_hxptar_68 = new TH2F("H2_hyptar_hxptar_68", "HMS X'_{tar} vs HMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.1, 0.1, 200, -0.1, 0.1);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	TH2F* H2_hyptar_hxptar_68_1 = new TH2F("H2_hyptar_hxptar_68_1", "HMS X'_{tar} vs HMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.1, 0.1, 200, -0.1, 0.1);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	List->Add(H2_hyptar_hxptar_68);
	List1->Add(H2_hyptar_hxptar_68_1);


	TH2F* H2_hyptar_hxptar_eep = new TH2F("H2_hyptar_hxptar_eep", "HMS X'_{tar} vs HMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.1, 0.1, 200, -0.1, 0.1);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	TH2F* H2_hyptar_hxptar_eep_1 = new TH2F("H2_hyptar_hxptar_eep_1", "HMS X'_{tar} vs HMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.1, 0.1, 200, -0.1, 0.1);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	List->Add(H2_hyptar_hxptar_eep);
	List1->Add(H2_hyptar_hxptar_eep_1);

	TH2F* H2_hyptar_hxptar_ee = new TH2F("H2_hyptar_hxptar_ee", "HMS X'_{tar} vs HMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.1, 0.1, 200, -0.1, 0.1);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	TH2F* H2_hyptar_hxptar_ee_1 = new TH2F("H2_hyptar_hxptar_ee_1", "HMS X'_{tar} vs HMS Y'_{tar}; X'_{etar} (rad); Y'_{etar} (rad)", 200, -0.1, 0.1, 200, -0.1, 0.1);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	List->Add(H2_hyptar_hxptar_ee);
	List1->Add(H2_hyptar_hxptar_ee_1);




	TH2F* H2_hYColl_hXColl_62 = new TH2F("H2_hYColl_hXColl_62", "HMS Collimator; HMS Y-Collimator; HMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	TH2F* H2_hYColl_hXColl_62_1 = new TH2F("H2_hYColl_hXColl_62_1", "HMS Collimator; HMS Y-Collimator; HMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	List->Add(H2_hYColl_hXColl_62);
	List1->Add(H2_hYColl_hXColl_62_1);

	TH2F* H2_hYColl_hXColl_68 = new TH2F("H2_hYColl_hXColl_68", "HMS Collimator; HMS Y-Collimator; HMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	TH2F* H2_hYColl_hXColl_68_1 = new TH2F("H2_hYColl_hXColl_68_1", "HMS Collimator; HMS Y-Collimator; HMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	List->Add(H2_hYColl_hXColl_68);
	List1->Add(H2_hYColl_hXColl_68_1);

	TH2F* H2_hYColl_hXColl_eep = new TH2F("H2_hYColl_hXColl_eep", "HMS Collimator; HMS Y-Collimator; HMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	TH2F* H2_hYColl_hXColl_eep_1 = new TH2F("H2_hYColl_hXColl_eep_1", "HMS Collimator; HMS Y-Collimator; HMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	List->Add(H2_hYColl_hXColl_eep);
	List1->Add(H2_hYColl_hXColl_eep_1);

	TH2F* H2_hYColl_hXColl_ee = new TH2F("H2_hYColl_hXColl_ee", "HMS Collimator; HMS Y-Collimator; HMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	TH2F* H2_hYColl_hXColl_ee_1 = new TH2F("H2_hYColl_hXColl_ee_1", "HMS Collimator; HMS Y-Collimator; HMS X-Collimator (cm)", 100, -15, 15, 100, -15, 15);					//Create a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	List->Add(H2_hYColl_hXColl_ee);
	List1->Add(H2_hYColl_hXColl_ee_1);


	
	


	int numLoops = 0;
	int max_loops = 70000000;

	nentries = inputtree_62->GetEntries();		//Get the total number of entries

	for (int k = 0; k <= nentries; k++)//nentries
	{
		inputtree_62->GetEntry(k);				//Get the ith entry from the T TTree

		weight_62 = 1 / ((Ehtrk_62 * Eetrk_62 * Emtrk_62 * Elt_62) * (charge_62 * areal_den * transparency * Z_per_A * Prot_Abs));

		if (
			evtyp_62 >= 4 &&
			pcal_min <= pcal_62 && pcal_62 <= pcal_max &&
			edelta_min <= edelta_62 && edelta_62 <= edelta_max &&
			preactz_min1 <= preactz_62 && preactz_62 <= preactz_max1 &&
			shms_Coll_gCut->IsInside(eYColl_62, eXColl_62) &&
			ctime_min <= ctime_62 && ctime_62 <= ctime_max &&//150
			hdelta_min <= hdelta_62 && hdelta_62 <= hdelta_max &&
			Em_min <= emiss_62 && emiss_62 <= Em_max
			)
		{
			H1_W_ratio_pid_acc_kin_data_eep->Fill(W_62, weight_62);
		}

		//0
		if (
			evtyp_62 >= 4 &&
			pcal_min <= pcal_62 && pcal_62 <= pcal_max &&
			edelta_min <= edelta_62 && edelta_62 <= edelta_max &&
			preactz_min1 <= preactz_62 && preactz_62 <= preactz_max1 &&
			W_min <= W_62 && W_62 <= W_max &&//0.9,1.05
			shms_Coll_gCut->IsInside(eYColl_62, eXColl_62) &&
			ctime_min <= ctime_62 && ctime_62 <= ctime_max &&//150
			hdelta_min <= hdelta_62 && hdelta_62 <= hdelta_max &&
			Em_min <= emiss_62 && emiss_62 <= Em_max
			)
		{
			H1_pcal_62->Fill(pcal_62, weight_62);
			H1_coin_62->Fill(ctime_62, weight_62);
			H1_preactz_62->Fill(preactz_62, weight_62);
			H1_exptar_62->Fill(exptar_62, weight_62);
			H1_eyptar_62->Fill(eyptar_62, weight_62);
			H1_hxptar_62->Fill(hxptar_62, weight_62);
			H1_hyptar_62->Fill(hyptar_62, weight_62);

			

			H1_edelta_62->Fill(edelta_62, weight_62);
			H1_hdelta_62->Fill(hdelta_62, weight_62);
			H1_emiss_62->Fill(emiss_62, weight_62);
			H1_pmiss_62->Fill(pmiss_62, weight_62);
			H1_pmiss_x_62->Fill(pmiss_x_62, weight_62);
			H1_pmiss_y_62->Fill(pmiss_y_62, weight_62);
			H1_pmiss_z_62->Fill(pmiss_z_62, weight_62);
			H1_W_62->Fill(W_62, weight_62);
			H1_Q2_62->Fill(Q2_62, weight_62);
			H2_hdelta_edelta_62->Fill(hdelta_62, edelta_62, weight_62);
			H2_eyptar_exptar_62->Fill(exptar_62, eyptar_62, weight_62);
			H2_hyptar_hxptar_62->Fill(hyptar_62, hxptar_62, weight_62);

			/*
			H1_eXColl_62->Fill(eXColl_62, weight_62);
			H1_eYColl_62->Fill(eYColl_62, weight_62);
			H1_hXColl_62->Fill(hXColl_62, weight_62);
			H1_hYColl_62->Fill(hYColl_62, weight_62);*/
			H2_eYColl_eXColl_62->Fill(eXColl_62, eYColl_62, weight_62);
			H2_hYColl_hXColl_62->Fill(hXColl_62, hYColl_62, weight_62);
		}


		//1
		if (
			evtyp_62 >= 4 &&
			pcal_min <= pcal_62 && pcal_62 <= pcal_max &&
			edelta_min <= edelta_62 && edelta_62 <= edelta_max &&
			preactz_min2 <= preactz_62 && preactz_62 <= preactz_max2 &&
			W_min <= W_62 && W_62 <= W_max &&
			shms_Coll_gCut->IsInside(eYColl_62, eXColl_62) &&
			ctime_min <= ctime_62 && ctime_62 <= ctime_max &&//150
			hdelta_min <= hdelta_62 && hdelta_62 <= hdelta_max &&
			Em_min <= emiss_62 && emiss_62 <= Em_max
			)
		{
			H1_pcal_62_1->Fill(pcal_62, weight_62);
			H1_coin_62_1->Fill(ctime_62, weight_62);
			H1_preactz_62_1->Fill(preactz_62, weight_62);
			H1_exptar_62_1->Fill(exptar_62, weight_62);
			H1_eyptar_62_1->Fill(eyptar_62, weight_62);
			H1_hxptar_62_1->Fill(hxptar_62, weight_62);
			H1_hyptar_62_1->Fill(hyptar_62, weight_62);



			H1_edelta_62_1->Fill(edelta_62, weight_62);
			H1_hdelta_62_1->Fill(hdelta_62, weight_62);
			H1_emiss_62_1->Fill(emiss_62, weight_62);
			H1_pmiss_62_1->Fill(pmiss_62, weight_62);
			H1_pmiss_x_62_1->Fill(pmiss_x_62, weight_62);
			H1_pmiss_y_62_1->Fill(pmiss_y_62, weight_62);
			H1_pmiss_z_62_1->Fill(pmiss_z_62, weight_62);
			H1_W_62_1->Fill(W_62, weight_62);
			H1_Q2_62_1->Fill(Q2_62, weight_62);
			H2_hdelta_edelta_62_1->Fill(hdelta_62, edelta_62, weight_62);
			H2_eyptar_exptar_62_1->Fill(exptar_62, eyptar_62, weight_62);
			H2_hyptar_hxptar_62_1->Fill(hyptar_62, hxptar_62, weight_62);

			/*
			H1_eXColl_62_1->Fill(eXColl_62, weight_62);
			H1_eYColl_62_1->Fill(eYColl_62, weight_62);
			H1_hXColl_62_1->Fill(hXColl_62, weight_62);
			H1_hYColl_62_1->Fill(hYColl_62, weight_62);*/
			H2_eYColl_eXColl_62_1->Fill(eXColl_62, eYColl_62, weight_62);
			H2_hYColl_hXColl_62_1->Fill(hXColl_62, hYColl_62, weight_62);
		}

		if (numLoops % 10000 == 0) { cout << numLoops << endl; } //leave the loop after completing the loop max_loops times
		//if (numLoops >= 10000000) { break; }
		numLoops++;
	}

	Draw1d(H1_pcal_62, "H1_pcal_62.png", "H1_pcal_62", outROOT);
	Draw1d(H1_pcal_62_1, "H1_pcal_62_1.png", "H1_pcal_62_1", outROOT);

	Draw1d(H1_coin_62, "H1_coin_62.png", "H1_coin_62", outROOT);
	Draw1d(H1_coin_62_1, "H1_coin_62_1.png", "H1_coin_62_1", outROOT);

	Draw1d(H1_preactz_62, "H1_preactz_62.png", "H1_preactz_62", outROOT);
	Draw1d(H1_preactz_62_1, "H1_preactz_62_1.png", "H1_preactz_62_1", outROOT);

	Draw1d(H1_exptar_62, "H1_exptar_62.png", "H1_exptar_62", outROOT);
	Draw1d(H1_exptar_62_1, "H1_exptar_62_1.png", "H1_exptar_62_1", outROOT);

	Draw1d(H1_eyptar_62, "H1_eyptar_62.png", "H1_eyptar_62", outROOT);
	Draw1d(H1_eyptar_62_1, "H1_eyptar_62_1.png", "H1_eyptar_62_1", outROOT);

	Draw1d(H1_hxptar_62, "H1_hxptar_62.png", "H1_hxptar_62", outROOT);
	Draw1d(H1_hxptar_62_1, "H1_hxptar_62_1.png", "H1_hxptar_62_1", outROOT);

	Draw1d(H1_hyptar_62, "H1_hyptar_62.png", "H1_hyptar_62", outROOT);
	Draw1d(H1_hyptar_62_1, "H1_hyptar_62_1.png", "H1_hyptar_62_1", outROOT);

	Draw1d(H1_edelta_62, "H1_edelta_62.png", "H1_edelta_62", outROOT);
	Draw1d(H1_edelta_62_1, "H1_edelta_62_1.png", "H1_edelta_62_1", outROOT);

	Draw1d(H1_hdelta_62, "H1_hdelta_62.png", "H1_hdelta_62", outROOT);
	Draw1d(H1_hdelta_62_1, "H1_hdelta_62_1.png", "H1_hdelta_62_1", outROOT);

	Draw1d(H1_emiss_62, "H1_emiss_62.png", "H1_emiss_62", outROOT);
	Draw1d(H1_emiss_62_1, "H1_emiss_62_1.png", "H1_emiss_62_1", outROOT);

	Draw1d(H1_pmiss_62, "H1_pmiss_62.png", "H1_pmiss_62", outROOT);
	Draw1d(H1_pmiss_62_1, "H1_pmiss_62_1.png", "H1_pmiss_62_1", outROOT);

	Draw1d(H1_pmiss_x_62, "H1_pmiss_x_62.png", "H1_pmiss_x_62", outROOT);
	Draw1d(H1_pmiss_x_62_1, "H1_pmiss_x_62_1.png", "H1_pmiss_x_62_1", outROOT);

	Draw1d(H1_pmiss_y_62, "H1_pmiss_y_62.png", "H1_pmiss_y_62", outROOT);
	Draw1d(H1_pmiss_y_62_1, "H1_pmiss_y_62_1.png", "H1_pmiss_y_62_1", outROOT);

	Draw1d(H1_pmiss_z_62, "H1_pmiss_z_62.png", "H1_pmiss_z_62", outROOT);
	Draw1d(H1_pmiss_z_62_1, "H1_pmiss_z_62_1.png", "H1_pmiss_z_62_1", outROOT);

	Draw1d(H1_W_62, "H1_W_62.png", "H1_W_62", outROOT);
	Draw1d(H1_W_62_1, "H1_W_62_1.png", "H1_W_62_1", outROOT);

	Draw1d(H1_Q2_62, "H1_Q2_62.png", "H1_Q2_62", outROOT);
	Draw1d(H1_Q2_62_1, "H1_Q2_62_1.png", "H1_Q2_62_1", outROOT);

	Draw2d(H2_eyptar_exptar_62, "H2_eyptar_exptar_62.png", "H2_eyptar_exptar_62", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_eyptar_exptar_62_1, "H2_eyptar_exptar_62_1.png", "H2_eyptar_exptar_62_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_hdelta_edelta_62, "H2_hdelta_edelta_62.png", "H2_hdelta_edelta_62", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_hdelta_edelta_62_1, "H2_hdelta_edelta_62_1.png", "H2_hdelta_edelta_62_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_hyptar_hxptar_62, "H2_hyptar_hxptar_62.png", "H2_hyptar_hxptar_62", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_hyptar_hxptar_62_1, "H2_hyptar_hxptar_62_1.png", "H2_hyptar_hxptar_62_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_eYColl_eXColl_62, "H2_eYColl_eXColl_62.png", "H2_eYColl_eXColl_62", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_eYColl_eXColl_62_1, "H2_eYColl_eXColl_62_1.png", "H2_eYColl_eXColl_62_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_hYColl_hXColl_62, "H2_hYColl_hXColl_62.png", "H2_hYColl_hXColl_62", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_hYColl_hXColl_62_1, "H2_hYColl_hXColl_62_1.png", "H2_hYColl_hXColl_62_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968


	


















	numLoops = 0;
	max_loops = 70000000;

	nentries = inputtree_68->GetEntries();		//Get the total number of entries

	for (int k = 0; k <= nentries; k++)//nentries
	{
		inputtree_68->GetEntry(k);				//Get the ith entry from the T TTree

		weight_68 = PS2_scalefac_68 / ((Eetrk_68 * Emtrk_68 * Elt_68) * (charge_68 * areal_den * transparency * Z_per_A));

		if (
			(evtyp_68 == 1 || evtyp_68 == 3 || evtyp_68 == 5 || evtyp_68 == 7) &&
			pcal_min <= pcal_68 && pcal_68 <= pcal_max &&
			edelta_min <= edelta_68 && edelta_68 <= edelta_max &&
			preactz_min1 <= preactz_68 && preactz_68 <= preactz_max1 &&
			shms_Coll_gCut->IsInside(eYColl_68, eXColl_68)
			)
		{
			H1_W_ratio_pid_acc_kin_data_ee->Fill(W_68, weight_68);
		}

		//0
		if (
			(evtyp_68 == 1 || evtyp_68 == 3 || evtyp_68 == 5 || evtyp_68 == 7) &&
			pcal_min <= pcal_68 && pcal_68 <= pcal_max &&
			edelta_min <= edelta_68 && edelta_68 <= edelta_max &&
			preactz_min1 <= preactz_68 && preactz_68 <= preactz_max1 &&
			W_min <= W_68 && W_68 <= W_max &&
			shms_Coll_gCut->IsInside(eYColl_68, eXColl_68)
			)
		{
			H1_pcal_68->Fill(pcal_68, weight_68);
			H1_coin_68->Fill(ctime_68, weight_68);
			H1_preactz_68->Fill(preactz_68, weight_68);
			H1_exptar_68->Fill(exptar_68, weight_68);
			H1_eyptar_68->Fill(eyptar_68, weight_68);
			H1_hxptar_68->Fill(hxptar_68, weight_68);
			H1_hyptar_68->Fill(hyptar_68, weight_68);
			H1_edelta_68->Fill(edelta_68, weight_68);
			H1_hdelta_68->Fill(hdelta_68, weight_68);
			H1_emiss_68->Fill(emiss_68, weight_68);
			H1_pmiss_68->Fill(pmiss_68, weight_68);
			H1_pmiss_x_68->Fill(pmiss_x_68, weight_68);
			H1_pmiss_y_68->Fill(pmiss_y_68, weight_68);
			H1_pmiss_z_68->Fill(pmiss_z_68, weight_68);
			H1_W_68->Fill(W_68, weight_68);
			H1_Q2_68->Fill(Q2_68, weight_68);
			H2_hdelta_edelta_68->Fill(hdelta_68, edelta_68, weight_68);
			H2_eyptar_exptar_68->Fill(exptar_68, eyptar_68, weight_68);
			H2_hyptar_hxptar_68->Fill(hyptar_68, hxptar_68, weight_68);

			/*H1_eXColl_68->Fill(eXColl_68, weight_68);
			H1_eYColl_68->Fill(eYColl_68, weight_68);
			H1_hXColl_68->Fill(hXColl_68, weight_68);
			H1_hYColl_68->Fill(hYColl_68, weight_68);*/
			H2_eYColl_eXColl_68->Fill(eXColl_68, eYColl_68, weight_68);
			H2_hYColl_hXColl_68->Fill(hXColl_68, hYColl_68, weight_68);
		}

		//1
		if (
			(evtyp_68 == 1 || evtyp_68 == 3 || evtyp_68 == 5 || evtyp_68 == 7) &&
			pcal_min <= pcal_68 && pcal_68 <= pcal_max &&
			edelta_min <= edelta_68 && edelta_68 <= edelta_max &&
			preactz_min2 <= preactz_68 && preactz_68 <= preactz_max2 &&
			W_min <= W_68 && W_68 <= W_max &&
			shms_Coll_gCut->IsInside(eYColl_68, eXColl_68)
			)
		{
			H1_pcal_68_1->Fill(pcal_68, weight_68);
			H1_coin_68_1->Fill(ctime_68, weight_68);
			H1_preactz_68_1->Fill(preactz_68, weight_68);
			H1_exptar_68_1->Fill(exptar_68, weight_68);
			H1_eyptar_68_1->Fill(eyptar_68, weight_68);
			H1_hxptar_68_1->Fill(hxptar_68, weight_68);
			H1_hyptar_68_1->Fill(hyptar_68, weight_68);
			H1_edelta_68_1->Fill(edelta_68, weight_68);
			H1_hdelta_68_1->Fill(hdelta_68, weight_68);
			H1_emiss_68_1->Fill(emiss_68, weight_68);
			H1_pmiss_68_1->Fill(pmiss_68, weight_68);
			H1_pmiss_x_68_1->Fill(pmiss_x_68, weight_68);
			H1_pmiss_y_68_1->Fill(pmiss_y_68, weight_68);
			H1_pmiss_z_68_1->Fill(pmiss_z_68, weight_68);
			H1_W_68_1->Fill(W_68, weight_68);
			H1_Q2_68_1->Fill(Q2_68, weight_68);
			H2_hdelta_edelta_68_1->Fill(hdelta_68, edelta_68, weight_68);
			H2_eyptar_exptar_68_1->Fill(exptar_68, eyptar_68, weight_68);
			H2_hyptar_hxptar_68_1->Fill(hyptar_68, hxptar_68, weight_68);

			/*H1_eXColl_68_1->Fill(eXColl_68, weight_68);
			H1_eYColl_68_1->Fill(eYColl_68, weight_68);
			H1_hXColl_68_1->Fill(hXColl_68, weight_68);
			H1_hYColl_68_1->Fill(hYColl_68, weight_68);*/
			H2_eYColl_eXColl_68_1->Fill(eXColl_68, eYColl_68, weight_68);
			H2_hYColl_hXColl_68_1->Fill(hXColl_68, hYColl_68, weight_68);
		}

		if (numLoops % 10000 == 0) { cout << numLoops << endl; } //leave the loop after completing the loop max_loops times
		//if (numLoops >= 100000000) { break; }
		numLoops++;
	}

	Draw1d(H1_pcal_68, "H1_pcal_68.png", "H1_pcal_68", outROOT);
	Draw1d(H1_pcal_68_1, "H1_pcal_68_1.png", "H1_pcal_68_1", outROOT);

	Draw1d(H1_preactz_68, "H1_preactz_68.png", "H1_preactz_68", outROOT);
	Draw1d(H1_preactz_68_1, "H1_preactz_68_1.png", "H1_preactz_68_1", outROOT);

	Draw1d(H1_exptar_68, "H1_exptar_68.png", "H1_exptar_68", outROOT);
	Draw1d(H1_exptar_68_1, "H1_exptar_68_1.png", "H1_exptar_68_1", outROOT);

	Draw1d(H1_eyptar_68, "H1_eyptar_68.png", "H1_eyptar_68", outROOT);
	Draw1d(H1_eyptar_68_1, "H1_eyptar_68_1.png", "H1_eyptar_68_1", outROOT);

	Draw1d(H1_hxptar_68, "H1_hxptar_68.png", "H1_hxptar_68", outROOT);
	Draw1d(H1_hxptar_68_1, "H1_hxptar_68_1.png", "H1_hxptar_68_1", outROOT);

	Draw1d(H1_hyptar_68, "H1_hyptar_68.png", "H1_hyptar_68", outROOT);
	Draw1d(H1_hyptar_68_1, "H1_hyptar_68_1.png", "H1_hyptar_68_1", outROOT);

	Draw1d(H1_edelta_68, "H1_edelta_68.png", "H1_edelta_68", outROOT);
	Draw1d(H1_edelta_68_1, "H1_edelta_68_1.png", "H1_edelta_68_1", outROOT);

	Draw1d(H1_hdelta_68, "H1_hdelta_68.png", "H1_hdelta_68", outROOT);
	Draw1d(H1_hdelta_68_1, "H1_hdelta_68_1.png", "H1_hdelta_68_1", outROOT);

	Draw1d(H1_emiss_68, "H1_emiss_68.png", "H1_emiss_68", outROOT);
	Draw1d(H1_emiss_68_1, "H1_emiss_68_1.png", "H1_emiss_68_1", outROOT);

	Draw1d(H1_pmiss_68, "H1_pmiss_68.png", "H1_pmiss_68", outROOT);
	Draw1d(H1_pmiss_68_1, "H1_pmiss_68_1.png", "H1_pmiss_68_1", outROOT);

	Draw1d(H1_pmiss_x_68, "H1_pmiss_x_68.png", "H1_pmiss_x_68", outROOT);
	Draw1d(H1_pmiss_x_68_1, "H1_pmiss_x_68_1.png", "H1_pmiss_x_68_1", outROOT);

	Draw1d(H1_pmiss_y_68, "H1_pmiss_y_68.png", "H1_pmiss_y_68", outROOT);
	Draw1d(H1_pmiss_y_68_1, "H1_pmiss_y_68_1.png", "H1_pmiss_y_68_1", outROOT);

	Draw1d(H1_pmiss_z_68, "H1_pmiss_z_68.png", "H1_pmiss_z_68", outROOT);
	Draw1d(H1_pmiss_z_68_1, "H1_pmiss_z_68_1.png", "H1_pmiss_z_68_1", outROOT);

	Draw1d(H1_W_68, "H1_W_68.png", "H1_W_68", outROOT);
	Draw1d(H1_W_68_1, "H1_W_68_1.png", "H1_W_68_1", outROOT);
	
	Draw1d(H1_Q2_68, "H1_Q2_68.png", "H1_Q2_68", outROOT);
	Draw1d(H1_Q2_68_1, "H1_Q2_68_1.png", "H1_Q2_68_1", outROOT);

	Draw2d(H2_eyptar_exptar_68, "H2_eyptar_exptar_68.png", "H2_eyptar_exptar_68", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_eyptar_exptar_68_1, "H2_eyptar_exptar_68_1.png", "H2_eyptar_exptar_68_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_hdelta_edelta_68, "H2_hdelta_edelta_68.png", "H2_hdelta_edelta_68", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_hdelta_edelta_68_1, "H2_hdelta_edelta_68_1.png", "H2_hdelta_edelta_68_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_hyptar_hxptar_68, "H2_hyptar_hxptar_68.png", "H2_hyptar_hxptar_68", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_hyptar_hxptar_68_1, "H2_hyptar_hxptar_68_1.png", "H2_hyptar_hxptar_68_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_eYColl_eXColl_68, "H2_eYColl_eXColl_68.png", "H2_eYColl_eXColl_68", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_eYColl_eXColl_68_1, "H2_eYColl_eXColl_68_1.png", "H2_eYColl_eXColl_68_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_hYColl_hXColl_68, "H2_hYColl_hXColl_68.png", "H2_hYColl_hXColl_68", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_hYColl_hXColl_68_1, "H2_hYColl_hXColl_68_1.png", "H2_hYColl_hXColl_68_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968









	
	



	



	




	
	numLoops = 0;
	max_loops = 70000000;

	nentries = inputtree_eep->GetEntries();		//Get the total number of entries

	for (int k = 0; k <= nentries; k++)//nentries
	{
		inputtree_eep->GetEntry(k);				//Get the ith entry from the T TTree

		FullWeight_eep = (Normfac_eep * weight_eep * 11.146) / ((nentries * areal_den * transparency * Z_per_A));

		//Define Collimator (same as in HCANA)
		htarx_corr_eep = tarx_eep - hxptar_eep * hreactz_eep * cos(h_angle_eep * dtr);
		etarx_corr_eep = tarx_eep - exptar_eep * preactz_eep * cos(e_angle_eep * dtr);
		hXColl_eep = htarx_corr_eep + hxptar_eep*168.;   //in cm
		hYColl_eep = hytar_eep + hyptar_eep*168.;
		eXColl_eep = etarx_corr_eep + exptar_eep*253.;
		eYColl_eep = eytar_eep + eyptar_eep*253.-(0.019+40.*.01*0.052)*edelta_eep+(0.00019+40*.01*.00052)*edelta_eep*edelta_eep; //correct for HB horizontal bend

		if (
			edelta_min <= edelta_eep && edelta_eep <= edelta_max &&
			hdelta_min <= hdelta_eep && hdelta_eep <= hdelta_max &&
			preactz_min1 <= preactz_eep && preactz_eep <= preactz_max1 &&
			shms_Coll_gCut->IsInside(eYColl_eep, eXColl_eep) &&
			Em_min <= emiss_eep && emiss_eep <= Em_max
			)
		{
			H1_W_ratio_pid_acc_kin_sim_eep->Fill(W_eep, FullWeight_eep);
		}

		//0
		if (
			//pcal_min <= pcal_eep && pcal_eep <= pcal_max &&
			//ctime
			edelta_min <= edelta_eep && edelta_eep <= edelta_max &&
			hdelta_min <= hdelta_eep && hdelta_eep <= hdelta_max &&
			preactz_min1 <= preactz_eep && preactz_eep <= preactz_max1 &&
			W_min <= W_eep && W_eep <= W_max &&
			shms_Coll_gCut->IsInside(eYColl_eep, eXColl_eep) &&
			Em_min <= emiss_eep && emiss_eep <= Em_max
			)
		{
			//H1_pcal_eep->Fill(pcal_eep, FullWeight_eep);
			//H1_coin_eep->Fill(coin_eep, FullWeight_eep);
			H1_preactz_eep->Fill(preactz_eep, FullWeight_eep);
			H1_exptar_eep->Fill(exptar_eep, FullWeight_eep);
			H1_eyptar_eep->Fill(eyptar_eep, FullWeight_eep);
			H1_hxptar_eep->Fill(hxptar_eep, FullWeight_eep);
			H1_hyptar_eep->Fill(hyptar_eep, FullWeight_eep);
			H1_edelta_eep->Fill(edelta_eep, FullWeight_eep);
			H1_hdelta_eep->Fill(hdelta_eep, FullWeight_eep);
			H1_emiss_eep->Fill(emiss_eep, FullWeight_eep);
			H1_pmiss_eep->Fill(pmiss_eep, FullWeight_eep);
			H1_pmiss_x_eep->Fill(pmiss_x_eep, FullWeight_eep);
			H1_pmiss_y_eep->Fill(pmiss_y_eep, FullWeight_eep);
			H1_pmiss_z_eep->Fill(pmiss_z_eep, FullWeight_eep);
			H1_W_eep->Fill(W_eep, FullWeight_eep);
			H1_Q2_eep->Fill(Q2_eep, FullWeight_eep);
			H2_hdelta_edelta_eep->Fill(hdelta_eep, edelta_eep, FullWeight_eep);
			H2_eyptar_exptar_eep->Fill(exptar_eep, eyptar_eep, FullWeight_eep);
			H2_hyptar_hxptar_eep->Fill(hyptar_eep, hxptar_eep, FullWeight_eep);

			//H1_eXColl_eep->Fill(eXColl_eep, weight_eep);
			//H1_eYColl_eep->Fill(eYColl_eep, weight_eep);
			//H1_hXColl_eep->Fill(hXColl_eep, weight_eep);
			//H1_hYColl_eep->Fill(hYColl_eep, weight_eep);
			H2_eYColl_eXColl_eep->Fill(eXColl_eep, eYColl_eep, weight_eep);
			H2_hYColl_hXColl_eep->Fill(hXColl_eep, hYColl_eep, weight_eep);
		}
		
		//1
		if (
			//pcal_min <= pcal_eep && pcal_eep <= pcal_max &&
			//ctime
			edelta_min <= edelta_eep && edelta_eep <= edelta_max &&
			hdelta_min <= hdelta_eep && hdelta_eep <= hdelta_max &&
			preactz_min2 <= preactz_eep && preactz_eep <= preactz_max2 &&
			W_min <= W_eep && W_eep <= W_max &&
			shms_Coll_gCut->IsInside(eYColl_eep, eXColl_eep) &&
			Em_min <= emiss_eep && emiss_eep <= Em_max
			)
		{
			//H1_pcal_eep_1->Fill(pcal_eep, FullWeight_eep);
			//H1_coin_eep_1->Fill(coin_eep, FullWeight_eep);
			H1_preactz_eep_1->Fill(preactz_eep, FullWeight_eep);
			H1_exptar_eep_1->Fill(exptar_eep, FullWeight_eep);
			H1_eyptar_eep_1->Fill(eyptar_eep, FullWeight_eep);
			H1_hxptar_eep_1->Fill(hxptar_eep, FullWeight_eep);
			H1_hyptar_eep_1->Fill(hyptar_eep, FullWeight_eep);
			H1_edelta_eep_1->Fill(edelta_eep, FullWeight_eep);
			H1_hdelta_eep_1->Fill(hdelta_eep, FullWeight_eep);
			H1_emiss_eep_1->Fill(emiss_eep, FullWeight_eep);
			H1_pmiss_eep_1->Fill(pmiss_eep, FullWeight_eep);
			H1_pmiss_x_eep_1->Fill(pmiss_x_eep, FullWeight_eep);
			H1_pmiss_y_eep_1->Fill(pmiss_y_eep, FullWeight_eep);
			H1_pmiss_z_eep_1->Fill(pmiss_z_eep, FullWeight_eep);
			H1_W_eep_1->Fill(W_eep, FullWeight_eep);
			H1_Q2_eep_1->Fill(Q2_eep, FullWeight_eep);
			H2_hdelta_edelta_eep_1->Fill(hdelta_eep, edelta_eep, FullWeight_eep);
			H2_eyptar_exptar_eep_1->Fill(exptar_eep, eyptar_eep, FullWeight_eep);
			H2_hyptar_hxptar_eep_1->Fill(hyptar_eep, hxptar_eep, FullWeight_eep);

			//H1_eXColl_eep_1->Fill(eXColl_eep, weight_eep);
			//H1_eYColl_eep_1->Fill(eYColl_eep, weight_eep);
			//H1_hXColl_eep_1->Fill(hXColl_eep, weight_eep);
			//H1_hYColl_eep_1->Fill(hYColl_eep, weight_eep);
			H2_eYColl_eXColl_eep_1->Fill(eXColl_eep, eYColl_eep, weight_eep);
			H2_hYColl_hXColl_eep_1->Fill(hXColl_eep, hYColl_eep, weight_eep);
		}

		if (numLoops % 10000 == 0) { cout << numLoops << endl; } //leave the loop after completing the loop max_loops times
		numLoops++;
	}
	
	Draw1d(H1_exptar_eep, "H1_exptar_eep.png", "H1_exptar_eep", outROOT);
	Draw1d(H1_exptar_eep_1, "H1_exptar_eep_1.png", "H1_exptar_eep_1", outROOT);

	Draw1d(H1_eyptar_eep, "H1_eyptar_eep.png", "H1_eyptar_eep", outROOT);
	Draw1d(H1_eyptar_eep_1, "H1_eyptar_eep_1.png", "H1_eyptar_eep_1", outROOT);

	Draw1d(H1_hxptar_eep, "H1_hxptar_eep.png", "H1_hxptar_eep", outROOT);
	Draw1d(H1_hxptar_eep_1, "H1_hxptar_eep_1.png", "H1_hxptar_eep_1", outROOT);

	Draw1d(H1_hyptar_eep, "H1_hyptar_eep.png", "H1_hyptar_eep", outROOT);
	Draw1d(H1_hyptar_eep_1, "H1_hyptar_eep_1.png", "H1_hyptar_eep_1", outROOT);

	Draw1d(H1_edelta_eep, "H1_edelta_eep.png", "H1_edelta_eep", outROOT);
	Draw1d(H1_edelta_eep_1, "H1_edelta_eep_1.png", "H1_edelta_eep_1", outROOT);

	Draw1d(H1_hdelta_eep, "H1_hdelta_eep.png", "H1_hdelta_eep", outROOT);
	Draw1d(H1_hdelta_eep_1, "H1_hdelta_eep_1.png", "H1_hdelta_eep_1", outROOT);

	Draw1d(H1_emiss_eep, "H1_emiss_eep.png", "H1_emiss_eep", outROOT);
	Draw1d(H1_emiss_eep_1, "H1_emiss_eep_1.png", "H1_emiss_eep_1", outROOT);

	Draw1d(H1_pmiss_eep, "H1_pmiss_eep.png", "H1_pmiss_eep", outROOT);
	Draw1d(H1_pmiss_eep_1, "H1_pmiss_eep_1.png", "H1_pmiss_eep_1", outROOT);

	Draw1d(H1_pmiss_x_eep, "H1_pmiss_x_eep.png", "H1_pmiss_x_eep", outROOT);
	Draw1d(H1_pmiss_x_eep_1, "H1_pmiss_x_eep_1.png", "H1_pmiss_x_eep_1", outROOT);

	Draw1d(H1_pmiss_y_eep, "H1_pmiss_y_eep.png", "H1_pmiss_y_eep", outROOT);
	Draw1d(H1_pmiss_y_eep_1, "H1_pmiss_y_eep_1.png", "H1_pmiss_y_eep_1", outROOT);

	Draw1d(H1_pmiss_z_eep, "H1_pmiss_z_eep.png", "H1_pmiss_z_eep", outROOT);
	Draw1d(H1_pmiss_z_eep_1, "H1_pmiss_z_eep_1.png", "H1_pmiss_z_eep_1", outROOT);

	Draw1d(H1_W_eep, "H1_W_eep.png", "H1_W_eep", outROOT);
	Draw1d(H1_W_eep_1, "H1_W_eep_1.png", "H1_W_eep_1", outROOT);

	Draw1d(H1_Q2_eep, "H1_Q2_eep.png", "H1_Q2_eep", outROOT);
	Draw1d(H1_Q2_eep_1, "H1_Q2_eep_1.png", "H1_Q2_eep_1", outROOT);


	Draw2d(H2_eyptar_exptar_eep, "H2_eyptar_exptar_eep.png", "H2_eyptar_exptar_eep", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_eyptar_exptar_eep_1, "H2_eyptar_exptar_eep_1.png", "H2_eyptar_exptar_eep_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_hdelta_edelta_eep, "H2_hdelta_edelta_eep.png", "H2_hdelta_edelta_eep", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_hdelta_edelta_eep_1, "H2_hdelta_edelta_eep_1.png", "H2_hdelta_edelta_eep_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_hyptar_hxptar_eep, "H2_hyptar_hxptar_eep.png", "H2_hyptar_hxptar_eep", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_hyptar_hxptar_eep_1, "H2_hyptar_hxptar_eep_1.png", "H2_hyptar_hxptar_eep_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_eYColl_eXColl_eep, "H2_eYColl_eXColl_eep.png", "H2_eYColl_eXColl_eep", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_eYColl_eXColl_eep_1, "H2_eYColl_eXColl_eep_1.png", "H2_eYColl_eXColl_eep_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_hYColl_hXColl_eep, "H2_hYColl_hXColl_eep.png", "H2_hYColl_hXColl_eep", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_hYColl_hXColl_eep_1, "H2_hYColl_hXColl_eep_1.png", "H2_hYColl_hXColl_eep_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	

	


























	
	numLoops = 0;
	max_loops = 70000000;

	nentries = inputtree_ee->GetEntries();		//Get the total number of entries

	for (int k = 0; k <= nentries; k++)//nentries
	{
		inputtree_ee->GetEntry(k);				//Get the ith entry from the T TTree

		FullWeight_ee = (Normfac_ee * weight_ee * 15.006) / ((nentries * areal_den * transparency * Z_per_A));

		//Define Collimator (same as in HCANA)
		htarx_corr_ee = tarx_ee - hxptar_ee * hreactz_ee * cos(h_angle_ee * dtr);
		etarx_corr_ee = tarx_ee - exptar_ee * preactz_ee * cos(e_angle_ee * dtr);
		hXColl_ee = htarx_corr_ee + hxptar_ee * 168.;   //in cm
		hYColl_ee = hytar_ee + hyptar_ee * 168.;
		eXColl_ee = etarx_corr_ee + exptar_ee * 253.;
		eYColl_ee = eytar_ee + eyptar_ee * 253. - (0.019 + 40. * .01 * 0.052) * edelta_ee + (0.00019 + 40 * .01 * .00052) * edelta_ee * edelta_ee; //correct for HB horizontal bend

		

		if (
			edelta_min <= edelta_ee && edelta_ee <= edelta_max &&
			preactz_min1 <= preactz_ee && preactz_ee <= preactz_max1 &&
			shms_Coll_gCut->IsInside(eYColl_ee, eXColl_ee)
			)
		{
			H1_W_ratio_pid_acc_kin_sim_ee->Fill(W_ee, FullWeight_ee);
		}


		//0
		if (
			edelta_min <= edelta_ee && edelta_ee <= edelta_max &&
			W_min <= W_ee && W_ee <= W_max &&
			preactz_min1 <= preactz_ee && preactz_ee <= preactz_max1 &&
			shms_Coll_gCut->IsInside(eYColl_ee, eXColl_ee)
			)
		{
			//H1_pcal_ee->Fill(pcal_ee, FullWeight_ee);
			H1_preactz_ee->Fill(preactz_ee, FullWeight_ee);
			H1_exptar_ee->Fill(exptar_ee, FullWeight_ee);
			H1_eyptar_ee->Fill(eyptar_ee, FullWeight_ee);
			H1_edelta_ee->Fill(edelta_ee, FullWeight_ee);
			H1_W_ee->Fill(W_ee, FullWeight_ee);
			H1_Q2_ee->Fill(Q2_ee, FullWeight_ee);
			H2_eyptar_exptar_ee->Fill(exptar_ee, eyptar_ee, FullWeight_ee);

			//H1_eXColl_ee->Fill(eXColl_ee, weight_ee);
			//H1_eYColl_ee->Fill(eYColl_ee, weight_ee);
			//H1_hXColl_ee->Fill(hXColl_ee, weight_ee);
			//H1_hYColl_ee->Fill(hYColl_ee, weight_ee);
			H2_eYColl_eXColl_ee->Fill(eXColl_ee, eYColl_ee, weight_ee);
			H2_hYColl_hXColl_ee->Fill(hXColl_ee, hYColl_ee, weight_ee);
		}

		//1
		if (
			edelta_min <= edelta_ee && edelta_ee <= edelta_max &&
			W_min <= W_ee && W_ee <= W_max &&
			preactz_min2 <= preactz_ee && preactz_ee <= preactz_max2 &&
			shms_Coll_gCut->IsInside(eYColl_ee, eXColl_ee)
			)
		{
			//H1_pcal_ee_1->Fill(pcal_ee, FullWeight_ee);
			H1_preactz_ee_1->Fill(preactz_ee, FullWeight_ee);
			H1_exptar_ee_1->Fill(exptar_ee, FullWeight_ee);
			H1_eyptar_ee_1->Fill(eyptar_ee, FullWeight_ee);
			H1_edelta_ee_1->Fill(edelta_ee, FullWeight_ee);
			H1_W_ee_1->Fill(W_ee, FullWeight_ee);
			H1_Q2_ee_1->Fill(Q2_ee, FullWeight_ee);
			H2_eyptar_exptar_ee_1->Fill(exptar_ee, eyptar_ee, FullWeight_ee);

			//H1_eXColl_ee_1->Fill(eXColl_ee, weight_ee);
			//H1_eYColl_ee_1->Fill(eYColl_ee, weight_ee);
			//H1_hXColl_ee_1->Fill(hXColl_ee, weight_ee);
			//H1_hYColl_ee_1->Fill(hYColl_ee, weight_ee);
			H2_eYColl_eXColl_ee_1->Fill(eXColl_ee, eYColl_ee, weight_ee);
			H2_hYColl_hXColl_ee_1->Fill(hXColl_ee, hYColl_ee, weight_ee);
		}

		if (numLoops % 10000 == 0) { cout << numLoops << endl; } //leave the loop after completing the loop max_loops times
		numLoops++;
	}
	
	Draw1d(H1_exptar_ee, "H1_exptar_ee.png", "H1_exptar_ee", outROOT);
	Draw1d(H1_exptar_ee_1, "H1_exptar_ee_1.png", "H1_exptar_ee_1", outROOT);

	Draw1d(H1_eyptar_ee, "H1_eyptar_ee.png", "H1_eyptar_ee", outROOT);
	Draw1d(H1_eyptar_ee_1, "H1_eyptar_ee_1.png", "H1_eyptar_ee_1", outROOT);

	Draw1d(H1_edelta_ee, "H1_edelta_ee.png", "H1_edelta_ee", outROOT);
	Draw1d(H1_edelta_ee_1, "H1_edelta_ee_1.png", "H1_edelta_ee_1", outROOT);

	Draw1d(H1_preactz_ee, "H1_preactz_ee.png", "H1_preactz_ee", outROOT);
	Draw1d(H1_preactz_ee_1, "H1_preactz_ee_1.png", "H1_preactz_ee_1", outROOT);

	Draw1d(H1_preactz_ee, "H1_preactz_ee.png", "H1_preactz_ee", outROOT);
	Draw1d(H1_preactz_ee_1, "H1_preactz_ee_1.png", "H1_preactz_ee_1", outROOT);

	Draw1d(H1_W_ee, "H1_W_ee.png", "H1_W_ee", outROOT);
	Draw1d(H1_W_ee_1, "H1_W_ee_1.png", "H1_W_ee_1", outROOT);

	Draw1d(H1_Q2_ee, "H1_Q2_ee.png", "H1_Q2_ee", outROOT);
	Draw1d(H1_Q2_ee_1, "H1_Q2_ee_1.png", "H1_Q2_ee_1", outROOT);

	Draw2d(H2_eyptar_exptar_ee, "H2_eyptar_exptar_ee.png", "H2_eyptar_exptar_ee", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_eyptar_exptar_ee_1, "H2_eyptar_exptar_ee_1.png", "H2_eyptar_exptar_ee_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968

	Draw2d(H2_eYColl_eXColl_ee, "H2_eYColl_eXColl_ee.png", "H2_eYColl_eXColl_ee", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968
	Draw2d(H2_eYColl_eXColl_ee_1, "H2_eYColl_eXColl_ee_1.png", "H2_eYColl_eXColl_ee_1", outROOT);																					//Draw a 2D Histogram of HMS X'_{tar} vs HMS Y'_{tar} from run 16968









	






	Double_t yf=1;
	
	OverLap2(H1_W_eep, H1_W_62, 1, 0, 1, yf, "H1_W_eep_over.png", "Invariant Mass", false, outROOT);						//Draw a 1D Histogram of Invariant Mass from run 16968
	OverLap2(H1_pmiss_eep, H1_pmiss_62, 1, 0, 1, yf, "H1_pmiss_eep_over.png", "Missing Momentum", false, outROOT);						//Draw a 1D Histogram of Invariant Mass from run 16968
	OverLap2(H1_emiss_eep, H1_emiss_62, 1, 0, 1, yf, "H1_emiss_eep_over.png", "Missing Energy", false, outROOT);						//Draw a 1D Histogram of Invariant Mass from run 16968
	OverLap2(H1_W_ee, H1_W_68, 1, 0, 1, yf, "H1_W_ee_over.png", "Invariant Mass", false, outROOT);						//Draw a 1D Histogram of Invariant Mass from run 16968


	OverLap2(H1_W_eep_1, H1_W_62_1, 1, 0, 1, yf, "H1_W_eep_over_1.png", "Invariant Mass", false, outROOT);						//Draw a 1D Histogram of Invariant Mass from run 16968
	OverLap2(H1_pmiss_eep_1, H1_pmiss_62_1, 1, 0, 1, yf, "H1_pmiss_eep_over_1.png", "Missing Momentum", false, outROOT);						//Draw a 1D Histogram of Invariant Mass from run 16968
	OverLap2(H1_emiss_eep_1, H1_emiss_62_1, 1, 0, 1, yf, "H1_emiss_eep_over_1.png", "Missing Energy", false, outROOT);						//Draw a 1D Histogram of Invariant Mass from run 16968
	OverLap2(H1_W_ee_1, H1_W_68_1, 1, 0, 1, yf, "H1_W_ee_over_1.png", "Invariant Mass", false, outROOT);						//Draw a 1D Histogram of Invariant Mass from run 16968


	TH1F* H1_W_eep_ratio = new TH1F("H1_W_eep_ratio", "Invariant Mass; W [GeV]; Counts / mC", 50, 0.85, 1.1);
	H1_W_eep_ratio->Sumw2();
	TH1F* H1_W_eep_ratio_1 = new TH1F("H1_W_eep_ratio_1", "Invariant Mass; W [GeV]; Counts / mC", 50, 0.85, 1.1);
	H1_W_eep_ratio_1->Sumw2();
	List->Add(H1_W_eep_ratio);
	List1->Add(H1_W_eep_ratio_1);



	TH1F* H1_W_ee_ratio = new TH1F("H1_W_ee_ratio", "Invariant Mass; W [GeV]; Counts / mC", 50, 0.85, 1.1);
	H1_W_ee_ratio->Sumw2();
	TH1F* H1_W_ee_ratio_1 = new TH1F("H1_W_ee_ratio_1", "Invariant Mass; W [GeV]; Counts / mC", 50, 0.85, 1.1);
	H1_W_ee_ratio_1->Sumw2();
	List->Add(H1_W_ee_ratio);
	List1->Add(H1_W_ee_ratio_1);

	H1_W_eep_ratio->Divide(H1_W_62, H1_W_eep);
	H1_W_eep_ratio_1->Divide(H1_W_62_1, H1_W_eep_1);

	H1_W_ee_ratio->Divide(H1_W_68, H1_W_ee);
	H1_W_ee_ratio_1->Divide(H1_W_68_1, H1_W_ee_1);


	Draw1d_W(H1_W_eep_ratio, "H1_W_eep_ratio.png", "H1_W_eep_ratio", outROOT);
	Draw1d_W(H1_W_eep_ratio_1, "H1_W_eep_ratio_1.png", "H1_W_eep_ratio_1", outROOT);


	Draw1d_W(H1_W_ee_ratio, "H1_W_ee_ratio.png", "H1_W_ee_ratio", outROOT);
	Draw1d_W(H1_W_ee_ratio_1, "H1_W_ee_ratio_1.png", "H1_W_ee_ratio_1", outROOT);
	


	//H1_W_ratio_pid_acc_kin_data_eep
	//H1_W_ratio_pid_acc_kin_sim_eep	
	Double_t num = H1_W_ratio_pid_acc_kin_data_ee->GetXaxis()->GetNbins();
	for (int i = 0; i <= num; i++)
	{
		cout << "data content: " << H1_W_ratio_pid_acc_kin_data_ee->GetBinContent(i) << endl;
		cout << "data error: " << H1_W_ratio_pid_acc_kin_data_ee->GetBinError(i) << endl;
	}
	cout << endl;
	cout << endl;
	

	num = H1_W_ratio_pid_acc_kin_sim_ee->GetXaxis()->GetNbins();
	for (int i = 0; i <= num; i++)
	{
		cout << "content: " << H1_W_ratio_pid_acc_kin_sim_ee->GetBinContent(i) << endl;
		cout << "error: " << H1_W_ratio_pid_acc_kin_sim_ee->GetBinError(i) << endl;
	}
	cout << endl;
	cout << endl;



	outROOT->mkdir("List");										//Make directories to store histograms based on Kinematic
	outROOT->cd("List");											//Write Kinematics histos to kin_plots directory
	List->Write();

	outROOT->mkdir("List1");									//Make directories to store histograms based on Kinematic
	outROOT->cd("List1");									//Write Kinematics histos to kin_plots directory
	List1->Write();

	outROOT->Close();	//Close File
	inROOT_62->Close();	//Close File
	inROOT_68->Close();	//Close File
	inROOT_eep->Close();//Close File
	inROOT_ee->Close();	//Close File
}
//*****************************************************************************************************************************************************************************
//End Main
//*****************************************************************************************************************************************************************************



void Draw1d_W(TH1F* hist1, const char* name, const char* title, TFile* outROOT)
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
	hist1->GetYaxis()->CenterTitle();
	hist1->GetYaxis()->SetRangeUser(0.0, 2.0);
	hist1->GetXaxis()->CenterTitle();
	hist1->GetXaxis()->SetRangeUser(0.85, 1.1);
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");

	hist1->GetYaxis()->SetTitleSize(0.0);//0.04, 0 for removing the label altoghether

	TCanvas* c = new TCanvas(name, title, 1366, 768); c->cd();

	hist1->Draw("histE");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void Draw1d(TH1F* hist1, const char* name, const char* title, TFile* outROOT)
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

	hist1->Draw("hist");

	outROOT->cd(); c->Write(); c->Print(name);

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



void Draw2d(TH2F* hist1, const char* name, const char* title, TFile* outROOT)
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
	
	hist1->Draw("Col");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}


void OverLap2(TH1F* hist1, TH1F* hist2, double xi, double yi, double xf, double yf, const char* name, const char* title, bool line, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1366, 768); c->cd();

	int font_type = 132;
	gStyle->SetOptStat(1);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
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

	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = hist1->GetMaximum();
	if (yaxisrange < hist2->GetMaximum()) { yaxisrange = hist2->GetMaximum(); }
	hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.3 * yaxisrange));

	//hist1->SetStats(0);
	hist1->SetLineWidth(2);
	hist1->SetLineColor(kRed);
	hist1->SetFillColorAlpha(kRed, 0.40);
	hist1->SetFillStyle(3004);
	hist1->Draw("histE0");

	//hist2->SetStats(0);
	hist2->SetLineWidth(2);
	hist2->SetLineColor(kBlue);
	hist2->SetFillColorAlpha(kBlue, 0.40);
	hist2->SetFillStyle(3005);
	hist2->Draw("samehistE0");

	float h1_I, h2_I;
	float h1_Ierr, h2_Ierr;
	float nbins = hist1->GetNbinsX();  //Get total number of bins (excluding overflow)
	h1_I = hist1->Integral();// AndError(0, nbins, h1_Ierr);//1, nbins, h1_Ierr
	h2_I = hist2->Integral();


	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box
	legend->AddEntry(hist1, Form("Sim | Integral: %.3f", h1_I), "f");
	legend->AddEntry(hist2, Form("Data | Integral: %.3f", h2_I));
	legend->Draw("SAME");

	if (line)
	{
		yf = hist1->GetMaximum();
		if (hist2->GetMaximum() > yf) { yf = hist2->GetMaximum(); }
		TLine* line = new TLine(xi, yi, xf, yf); line->SetLineColor(1); line->SetLineWidth(2); line->Draw("SAME");
	}
	
	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void OverLap3(TH1F* hist1, TH1F* hist2, TH1F* hist3, double xi, double yi, double xf, double yf, const char* name, const char* title, const char* xaxis, bool line, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1366, 768); c->cd();

	int font_type = 132;

	gStyle->SetOptStat(1);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->Scale(1. / hist1->Integral(), "width");
	hist2->Scale(1. / hist2->Integral(), "width");
	hist3->Scale(1. / hist3->Integral(), "width");

	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");

	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = hist1->GetMaximum();
	if (yaxisrange < hist2->GetMaximum()) { yaxisrange = hist2->GetMaximum(); }
	if (yaxisrange < hist3->GetMaximum()) { yaxisrange = hist3->GetMaximum(); }
	hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.25 * yaxisrange));

	
	hist1->SetLineWidth(2);
	hist1->SetLineColor(kRed);
	hist1->SetFillColorAlpha(kRed, 0.40);
	hist1->SetFillStyle(3004);
	hist1->SetTitle(title);
	hist1->GetXaxis()->SetTitle(xaxis);
	hist1->Draw("histE0");

	hist2->SetLineWidth(2);
	hist2->SetLineColor(kGreen);
	hist2->SetFillColorAlpha(kGreen, 0.40);
	hist2->SetFillStyle(3005);
	//hist2->SetTitle(title);
	//hist2->GetXaxis()->SetTitle(xaxis);
	hist2->Draw("samehistE0");

	hist3->SetLineWidth(2);
	hist3->SetLineColor(kBlue);
	hist3->SetFillColorAlpha(kBlue, 0.40);
	hist3->SetFillStyle(3006);
	//hist3->SetTitle(title);
	//hist3->GetXaxis()->SetTitle(xaxis);
	hist3->Draw("samehistE0");
	


	

	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box
	legend->AddEntry(hist1, "Ca40", "lp"); legend->AddEntry(hist2, "Ca48", "lp"); legend->AddEntry(hist3, "Fe54", "lp"); legend->Draw("SAME");


	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void OverLap5(TH1F* hist1, TH1F* hist2, TH1F* hist3, TH1F* hist4, TH1F* hist5, double xi, double yi, double xf, double yf, const char* name, const char* title, const char* xaxis, bool mf, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(1);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->SetStats(0);
	hist1->SetLineWidth(2);
	hist1->SetLineColor(kRed);
	hist1->SetFillColorAlpha(kRed, 0.40);
	hist1->SetFillStyle(3004);
	hist1->SetTitle(title);
	hist1->GetXaxis()->SetTitle(xaxis);
	hist1->Scale(1. / hist1->Integral(), "width");

	hist2->SetStats(0);
	hist2->SetLineWidth(2);
	hist2->SetLineColor(kMagenta);
	hist2->SetFillColorAlpha(kMagenta, 0.40);
	hist2->SetFillStyle(3005);
	hist2->SetTitle(title);
	hist2->GetXaxis()->SetTitle(xaxis);
	hist2->Scale(1. / hist2->Integral(), "width");

	hist3->SetStats(0);
	hist3->SetLineWidth(2);
	hist3->SetLineColor(kGreen);
	hist3->SetFillColorAlpha(kGreen, 0.40);
	hist3->SetFillStyle(3006);
	hist3->SetTitle(title);
	hist3->GetXaxis()->SetTitle(xaxis);
	hist3->Scale(1. / hist3->Integral(), "width");

	hist4->SetStats(0);
	hist4->SetLineWidth(2);
	hist4->SetLineColor(kCyan);
	hist4->SetFillColorAlpha(kCyan, 0.40);
	hist4->SetFillStyle(3007);
	hist4->SetTitle(title);
	hist4->GetXaxis()->SetTitle(xaxis);
	hist4->Scale(1. / hist4->Integral(), "width");

	hist5->SetStats(0);
	hist5->SetLineWidth(2);
	hist5->SetLineColor(kBlue);
	hist5->SetFillColorAlpha(kBlue, 0.40);
	hist5->SetFillStyle(3008);
	hist5->SetTitle(title);
	hist5->GetXaxis()->SetTitle(xaxis);
	hist5->Scale(1. / hist5->Integral(), "width");

	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = hist1->GetMaximum();
	if (yaxisrange < hist2->GetMaximum()) { yaxisrange = hist2->GetMaximum(); }
	if (yaxisrange < hist3->GetMaximum()) { yaxisrange = hist3->GetMaximum(); }
	if (yaxisrange < hist4->GetMaximum()) { yaxisrange = hist4->GetMaximum(); }
	if (yaxisrange < hist5->GetMaximum()) { yaxisrange = hist5->GetMaximum(); }
	hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.25 * yaxisrange));
	
	// set histogram titles/labels/font
	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");

	TCanvas* c = new TCanvas(name, title, 1366, 768); c->cd();
	
	if (mf)	{hist1->Draw("HIST"); hist2->Draw("HIST SAME"); hist3->Draw("HIST SAME"); hist4->Draw("HIST SAME"); hist5->Draw("HIST SAME");}
	else	{hist3->Draw("HIST"); hist1->Draw("HIST SAME"); hist2->Draw("HIST SAME"); hist4->Draw("HIST SAME"); hist5->Draw("HIST SAME");}

	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box
	legend->AddEntry(hist1, "D2", "lp"); legend->AddEntry(hist2, "Be9", "lp"); legend->AddEntry(hist3, "B10", "lp"); legend->AddEntry(hist4, "B11", "lp"); legend->AddEntry(hist5, "C12", "lp"); legend->Draw("SAME");

	TLine* line = new TLine(xi, yi, xf, yf); line->SetLineColor(1); line->SetLineWidth(2); line->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);
	
	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void Compare2(TH2F* hist1, TH2F* hist2, double xi, double yi, double xf, double yf, const char* name, const char* title, bool line, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1366, 768); c->Divide(2, 1);

	int font_type = 132;

	gStyle->SetOptStat(1);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.12);
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

	hist2->GetXaxis()->SetLabelSize(0.04);
	hist2->GetYaxis()->SetLabelSize(0.04);
	hist2->GetYaxis()->CenterTitle();
	hist2->GetXaxis()->CenterTitle();
	hist2->SetLabelFont(font_type, "XY");
	hist2->SetTitleFont(font_type, "XY");
	hist2->SetTitleSize(0.05, "XY");
	hist2->SetTitleOffset(1., "XY");

	if (line) {
		c->cd(2); hist1->Draw("Col"); TLine* line1 = new TLine(xi, yi, xf, yf); line1->SetLineColor(2); line1->SetLineWidth(2); line1->Draw("SAME");
		c->cd(1); hist2->Draw("Col"); TLine* line2 = new TLine(xi, yi, xf, yf); line2->SetLineColor(2); line2->SetLineWidth(2); line2->Draw("SAME");
	}
	else {
		c->cd(2); hist1->Draw("Col");
		c->cd(1); hist2->Draw("Col");
	}
	
	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void OverLay3(TGraph* hist1, TGraph* hist2, TGraph* hist3, const char* name, const char* title, const char* xaxis, bool hms, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1366, 768); c->cd();

	if (hms) { hist1->GetYaxis()->SetRangeUser(250, 450); hist2->GetYaxis()->SetRangeUser(250, 450); hist3->GetYaxis()->SetRangeUser(250, 450); }
	else	 {hist1->GetYaxis()->SetRangeUser(200, 350); hist2->GetYaxis()->SetRangeUser(200, 350); hist3->GetYaxis()->SetRangeUser(200, 350);}
	
	//hist1->GetXaxis()->SetBinLabel(7, "1u1"); hist1->GetXaxis()->SetBinLabel(14, "1u2"); hist1->GetXaxis()->SetBinLabel(21, "2u1"); hist1->GetXaxis()->SetBinLabel(28, "2u2");
	//hist1->GetXaxis()->SetBinLabel(35, "1x1"); hist1->GetXaxis()->SetBinLabel(42, "1x2"); hist1->GetXaxis()->SetBinLabel(49, "2x1"); hist1->GetXaxis()->SetBinLabel(56, "2x2");
	//hist1->GetXaxis()->SetBinLabel(63, "1v1"); hist1->GetXaxis()->SetBinLabel(70, "1v2"); hist1->GetXaxis()->SetBinLabel(77, "2v1"); hist1->GetXaxis()->SetBinLabel(84, "2v2");

	hist1->SetMarkerColor(2); hist1->SetMarkerStyle(20); hist1->SetMarkerSize(2); hist1->SetLineColor(2); hist1->SetLineWidth(2); hist1->SetTitle(title); hist1->GetXaxis()->SetTitle(xaxis); hist1->Draw("");
	hist2->SetMarkerColor(8); hist2->SetMarkerStyle(21); hist2->SetMarkerSize(2); hist2->SetLineColor(8); hist2->SetLineWidth(2); hist2->SetTitle(title); hist2->GetXaxis()->SetTitle(xaxis); hist2->Draw("PL SAME");
	hist3->SetMarkerColor(4); hist3->SetMarkerStyle(22); hist3->SetMarkerSize(2); hist3->SetLineColor(4); hist3->SetLineWidth(2); hist3->SetTitle(title); hist3->GetXaxis()->SetTitle(xaxis); hist3->Draw("PL SAME");

	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box
	legend->AddEntry(hist1, "Ca40", "lp"); legend->AddEntry(hist2, "Ca48", "lp"); legend->AddEntry(hist3, "Fe54", "lp"); legend->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void OverLay5(TGraph* hist1, TGraph* hist2, TGraph* hist3, TGraph* hist4, TGraph* hist5, const char* name, const char* title, const char* xaxis, bool mf, TFile * outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1366, 768); c->cd();

	if (mf) {hist1->GetYaxis()->SetRangeUser(250, 450); hist2->GetYaxis()->SetRangeUser(250, 450); hist3->GetYaxis()->SetRangeUser(250, 450); hist4->GetYaxis()->SetRangeUser(250, 450); hist5->GetYaxis()->SetRangeUser(250, 450);}
	else	{hist1->GetYaxis()->SetRangeUser(200, 350); hist2->GetYaxis()->SetRangeUser(200, 350); hist3->GetYaxis()->SetRangeUser(200, 350); hist4->GetYaxis()->SetRangeUser(200, 350); hist5->GetYaxis()->SetRangeUser(200, 350);}

	//hist1->GetXaxis()->SetBinLabel(10, "1u1"); hist1->GetXaxis()->SetBinLabel(20, "1u2"); hist1->GetXaxis()->SetBinLabel(30, "2u1"); hist1->GetXaxis()->SetBinLabel(40, "2u2");
	//hist1->GetXaxis()->SetBinLabel(50, "1x1"); hist1->GetXaxis()->SetBinLabel(60, "1x2"); hist1->GetXaxis()->SetBinLabel(70, "2x1"); hist1->GetXaxis()->SetBinLabel(80, "2x2");
	//hist1->GetXaxis()->SetBinLabel(90, "1v1"); hist1->GetXaxis()->SetBinLabel(100, "1v2"); hist1->GetXaxis()->SetBinLabel(110, "2v1"); hist1->GetXaxis()->SetBinLabel(120, "2v2");

	hist1->SetMarkerColor(2); hist1->SetMarkerStyle(20); hist1->SetMarkerSize(2); hist1->SetLineColor(2); hist1->SetLineWidth(2); hist1->SetTitle(title); hist1->GetXaxis()->SetTitle(xaxis); hist1->Draw("");
	hist2->SetMarkerColor(6); hist2->SetMarkerStyle(21); hist2->SetMarkerSize(2); hist2->SetLineColor(6); hist2->SetLineWidth(2); hist2->SetTitle(title); hist2->GetXaxis()->SetTitle(xaxis); hist2->Draw("PL SAME");
	hist3->SetMarkerColor(8); hist3->SetMarkerStyle(22); hist3->SetMarkerSize(2); hist3->SetLineColor(8); hist3->SetLineWidth(2); hist3->SetTitle(title); hist3->GetXaxis()->SetTitle(xaxis); hist3->Draw("PL SAME");
	hist4->SetMarkerColor(7); hist4->SetMarkerStyle(23); hist4->SetMarkerSize(2); hist4->SetLineColor(7); hist4->SetLineWidth(2); hist4->SetTitle(title); hist4->GetXaxis()->SetTitle(xaxis); hist4->Draw("PL SAME");
	hist5->SetMarkerColor(4); hist5->SetMarkerStyle(29); hist5->SetMarkerSize(2); hist5->SetLineColor(4); hist5->SetLineWidth(2); hist5->SetTitle(title); hist5->GetXaxis()->SetTitle(xaxis); hist5->Draw("PL SAME");
	
	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box,
	legend->AddEntry(hist1, "D2", "lp"); legend->AddEntry(hist2, "Be9", "lp"); legend->AddEntry(hist3, "B10", "lp"); legend->AddEntry(hist4, "B11", "lp"); legend->AddEntry(hist5, "C12", "lp"); legend->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);
	
	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}

























void compare_histos_light(
	TString file1_path, TString hist1,
	TString file2_path, TString hist2,
	TString file3_path, TString hist3,
	TString file4_path, TString hist4,
	TString file5_path, TString hist5,
	TString xlabel, TString ylabel, TString title,
	TString hist1_leg, TString hist2_leg, TString hist3_leg, TString hist4_leg, TString hist5_leg,
	bool norm,
	TFile* outROOT)
{

	int font_type = 132;

	gStyle->SetOptStat(1);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);


	//Open  ROOT files;
	TFile* file1 = NULL;
	TFile* file2 = NULL;
	TFile* file3 = NULL;
	TFile* file4 = NULL;
	TFile* file5 = NULL;

	file1 = new TFile(file1_path.Data());
	file2 = new TFile(file2_path.Data());
	file3 = new TFile(file3_path.Data());
	file4 = new TFile(file4_path.Data());
	file5 = new TFile(file5_path.Data());

	// declare 1D histos
	TH1F* H_hist1 = 0;
	TH1F* H_hist2 = 0;
	TH1F* H_hist3 = 0;
	TH1F* H_hist4 = 0;
	TH1F* H_hist5 = 0;

	// get histogram objects
	file1->cd();
	file1->GetObject(hist1.Data(), H_hist1);
	file2->cd();
	file2->GetObject(hist2.Data(), H_hist2);
	file3->cd();
	file3->GetObject(hist3.Data(), H_hist3);
	file4->cd();
	file4->GetObject(hist4.Data(), H_hist4);
	file5->cd();
	file5->GetObject(hist5.Data(), H_hist5);

	double h1_I, h2_I, h3_I, h4_I, h5_I;
	double h1_Ierr, h2_Ierr, h3_Ierr, h4_Ierr, h5_Ierr;
	double nbins = H_hist1->GetNbinsX();  //Get total number of bins (excluding overflow)
	h1_I = H_hist1->IntegralAndError(1, nbins, h1_Ierr);
	h2_I = H_hist2->IntegralAndError(1, nbins, h2_Ierr);
	h3_I = H_hist3->IntegralAndError(1, nbins, h3_Ierr);
	h4_I = H_hist4->IntegralAndError(1, nbins, h4_Ierr);
	h5_I = H_hist5->IntegralAndError(1, nbins, h5_Ierr);

	H_hist1->SetLineWidth(2);
	H_hist2->SetLineWidth(2);
	H_hist3->SetLineWidth(2);
	H_hist4->SetLineWidth(2);
	H_hist5->SetLineWidth(2);

	// set histos aethetics
	H_hist1->SetLineColor(kRed);
	H_hist1->SetFillColorAlpha(kRed, 0.40);
	H_hist1->SetFillStyle(3004);
	H_hist1->Scale(1. / H_hist1->Integral(), "width");

	H_hist2->SetLineColor(kMagenta);
	H_hist2->SetFillColorAlpha(kMagenta, 0.40);
	H_hist2->SetFillStyle(3005);
	H_hist2->Scale(1. / H_hist2->Integral(), "width");

	H_hist3->SetLineColor(kGreen);
	H_hist3->SetFillColorAlpha(kGreen, 0.40);
	H_hist3->SetFillStyle(3006);
	H_hist3->Scale(1. / H_hist3->Integral(), "width");

	H_hist4->SetLineColor(kCyan);
	H_hist4->SetFillColorAlpha(kCyan, 0.40);
	H_hist4->SetFillStyle(3007);
	H_hist4->Scale(1. / H_hist4->Integral(), "width");

	H_hist5->SetLineColor(kBlue);
	H_hist5->SetFillColorAlpha(kBlue, 0.40);
	H_hist5->SetFillStyle(3008);
	H_hist5->Scale(1. / H_hist5->Integral(), "width");

	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = H_hist1->GetMaximum();
	if (yaxisrange < H_hist2->GetMaximum()) { yaxisrange = H_hist2->GetMaximum(); }
	if (yaxisrange < H_hist3->GetMaximum()) { yaxisrange = H_hist3->GetMaximum(); }
	if (yaxisrange < H_hist4->GetMaximum()) { yaxisrange = H_hist4->GetMaximum(); }
	if (yaxisrange < H_hist5->GetMaximum()) { yaxisrange = H_hist5->GetMaximum(); }
	H_hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.25 * yaxisrange));


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


	TCanvas* c = new TCanvas("c", "c", 1366, 768);

	H_hist1->Draw("histE0");
	H_hist2->Draw("sameshistE0");
	H_hist3->Draw("sameshistE0");
	H_hist4->Draw("sameshistE0");
	H_hist5->Draw("sameshistE0");

	// create legend ( displays hist legend label and integral counts)
	TLegend* leg = new TLegend(0.14, 0.89, 0.25, 0.78);
	leg->AddEntry(H_hist1, Form("%s | Integral: %.3f", hist1_leg.Data(), h1_I), "f");
	leg->AddEntry(H_hist2, Form("%s | Integral: %.3f", hist2_leg.Data(), h2_I));
	leg->AddEntry(H_hist3, Form("%s | Integral: %.3f", hist3_leg.Data(), h3_I));
	leg->AddEntry(H_hist4, Form("%s | Integral: %.3f", hist4_leg.Data(), h4_I));
	leg->AddEntry(H_hist5, Form("%s | Integral: %.3f", hist5_leg.Data(), h5_I));
	// draw legend
	leg->Draw();

	outROOT->cd(); c->Write(); c->Print(Form("%s.png", title.Data()));

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}

















void compare_histos_heavy(
	TString file1_path, TString hist1,
	TString file2_path, TString hist2,
	TString file3_path, TString hist3,
	TString xlabel, TString ylabel, TString title,
	TString hist1_leg, TString hist2_leg, TString hist3_leg,
	bool norm,
	TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(1);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);


	//Open  ROOT files;
	TFile* file1 = NULL;
	TFile* file2 = NULL;
	TFile* file3 = NULL;

	file1 = new TFile(file1_path.Data());
	file2 = new TFile(file2_path.Data());
	file3 = new TFile(file3_path.Data());

	// declare 1D histos
	TH1F* H_hist1 = 0;
	TH1F* H_hist2 = 0;
	TH1F* H_hist3 = 0;

	// get histogram objects
	file1->cd();
	file1->GetObject(hist1.Data(), H_hist1);
	file2->cd();
	file2->GetObject(hist2.Data(), H_hist2);
	file3->cd();
	file3->GetObject(hist3.Data(), H_hist3);

	double h1_I, h2_I, h3_I;
	double h1_Ierr, h2_Ierr, h3_Ierr;
	double nbins = H_hist1->GetNbinsX();  //Get total number of bins (excluding overflow)
	h1_I = H_hist1->IntegralAndError(1, nbins, h1_Ierr);
	h2_I = H_hist2->IntegralAndError(1, nbins, h2_Ierr);
	h3_I = H_hist3->IntegralAndError(1, nbins, h3_Ierr);

	H_hist1->SetLineWidth(2);
	H_hist2->SetLineWidth(2);
	H_hist3->SetLineWidth(2);

	// set histos aethetics
	H_hist1->SetLineColor(kRed);
	H_hist1->SetFillColorAlpha(kRed, 0.40);
	H_hist1->SetFillStyle(3004);
	H_hist1->Scale(1. / H_hist1->Integral(), "width");

	H_hist2->SetLineColor(kGreen);
	H_hist2->SetFillColorAlpha(kGreen, 0.40);
	H_hist2->SetFillStyle(3005);
	H_hist2->Scale(1. / H_hist2->Integral(), "width");

	H_hist3->SetLineColor(kBlue);
	H_hist3->SetFillColorAlpha(kBlue, 0.40);
	H_hist3->SetFillStyle(3006);
	H_hist3->Scale(1. / H_hist3->Integral(), "width");

	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = H_hist1->GetMaximum();
	if (yaxisrange < H_hist2->GetMaximum()) { yaxisrange = H_hist2->GetMaximum(); }
	if (yaxisrange < H_hist3->GetMaximum()) { yaxisrange = H_hist3->GetMaximum(); }
	H_hist1->GetYaxis()->SetRangeUser(0, yaxisrange + 0.25 * yaxisrange);

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


	TCanvas* c = new TCanvas("c", "c", 1366, 768);

	H_hist1->Draw("histE0");
	H_hist2->Draw("sameshistE0");
	H_hist3->Draw("sameshistE0");

	// create legend ( displays hist legend label and integral counts)
	TLegend* leg = new TLegend(0.14, 0.89, 0.25, 0.78);
	leg->AddEntry(H_hist1, Form("%s | Integral: %.3f", hist1_leg.Data(), h1_I), "f");
	leg->AddEntry(H_hist2, Form("%s | Integral: %.3f", hist2_leg.Data(), h2_I));
	leg->AddEntry(H_hist3, Form("%s | Integral: %.3f", hist3_leg.Data(), h3_I));
	leg->Draw();

	outROOT->cd(); c->Write(); c->Print(Form("%s.png", title.Data()));

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}