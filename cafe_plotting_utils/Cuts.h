#ifndef Cuts_H
#define Cuts_H
#include <string>
using namespace std;

class Cuts
{
	public:
	//----- DATA ANALYSIS CUTS -----
	//lumi			baseCuts = e_delta>=-10		&& e_delta<=22		&& PID_SHMS
	//optics		baseCuts = PID_SHMS
	//heep_singles	baseCuts = Accp_SHMS	&& PID_SHMS		&& kinHeepSing	&& cTime
	//heep_coin		baseCuts = Accp_HMS		&& Accp_SHMS	&& Accp_Z		&& PID_HMS	&& PID_SHMS		&& kinHeepCoin
	//MF			baseCuts = Accp_HMS		&& Accp_SHMS	&& Accp_Z		&& PID_HMS	&& PID_SHMS		&& kinMF
	//SRC			baseCuts = Accp_HMS		&& Accp_SHMS	&& Accp_Z		&& PID_HMS	&& PID_SHMS		&& kinSRC

	//----- SIMC ANALYSIS CUTS-----
	//heep_coin		baseCuts = Accp_HMS		&& Accp_SHMS	&& Accp_Z		&& kinHeepCoin
	//MF			baseCuts = Accp_HMS		&& Accp_SHMS	&& Accp_Z		&& kinMF
	//SRC			baseCuts = Accp_HMS		&& Accp_SHMS	&& Accp_Z		&& kinSRC

	//lumi_baseCuts = e_delta>=-10 && e_delta<=22 && Cuts:PID_SHMS(petot_trkNorm, pngcer_npeSum_pid, phgcer_npeSum_pid)

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//---------------------------------------------------------------------------Tracking Eff---------------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//----- HMS Tracking Efficiency-----
	static bool HMS_Tr(double hdc_ntrk, double hScinGood, double hcer_npeSum, double hetotnorm, double hBetaNtrk) 
	{
		Double_t hdc_ntrk_min = 1;							//minimum number of tracks required
		Double_t hnpeSum_min = 0, hnpeSum_max = 0.5;		//HMS Cherenkov Total Number of Photo - Electrons(search "npeSum" in / hcana / src / THcCherenkov.cxx)
		Double_t hetotnorm_min = 0, hetotnorm_max = 0.7;	//HMS Calorimeter Total Energy / Spec.Central Momentum(search "etotnorm" in / hcana / src / THcShower.cxx)
		Double_t hBetaNtrk_min = 0.5, hBetaNtrk_max = 1.5;	//HMS Beta(p / E) Calculated from Scintillator Hits(without using Tracking Info) (search "betanotrack" in / hcana / src / THcHodoscope.cxx)

		if(hdc_ntrk_min > hdc_ntrk && hdc_ntrk != -100) { return false; }
		if(hScinGood != 1 && hScinGood != -100) { return false; }
		if( (hnpeSum_min > hcer_npeSum			|| hcer_npeSum > hnpeSum_max)		&& hcer_npeSum != -100)		{ return false; }
		if( (hetotnorm_min > hetotnorm			|| hetotnorm > hetotnorm_max)		&& hetotnorm != -100)		{ return false; }
		if( (hBetaNtrk_min > hBetaNtrk			|| hBetaNtrk > hBetaNtrk_max)		&& hBetaNtrk != -100)		{ return false; }
		return true;
	}

	//---- SHMS Tracking Efficiency-----
	static bool SHMS_Tr(double pdc_ntrk, double pScinGood, double pngcer_npeSum, double phgcer_npeSum, double petotnorm, double pBetaNtrk) 
	{
		Double_t pdc_ntrk_min = 1; //minimum number of tracks required
		Double_t pngcer_npeSum_min = 0.5, pngcer_npeSum_max = 100;	//SHMS Noble Gas Cherenkov Total Number of Photo - Electrons(search "npeSum" in / hcana / src / THcCherenkov.cxx)
		Double_t phgcer_npeSum_min = 0, phgcer_npeSum_max = 1.5;	//SHMS Heavy Gas Cherenkov(will NOT be used by CaFe experiment)
		Double_t petotnorm_min = 0.7, petotnorm_max = 1.1;			//SHMS Calorimeter Total Energy / Spec.Central Momentum(search "etotnorm" in / hcana / src / THcShower.cxx)
		Double_t pBetaNtrk_min = 0.5, pBetaNtrk_max = 1.5;			//SHMS Beta(p / E) Calculated from Scintillator Hits(without using Tracking Info) (search "betanotrack" in / hcana / src / THcHodoscope.cxx)

		if (pdc_ntrk_min > pdc_ntrk) { return false; }
		if(pScinGood != 1 && pScinGood != -100) { return false; }
		if( (pngcer_npeSum_min > pngcer_npeSum		|| pngcer_npeSum > pngcer_npeSum_max)	&& pngcer_npeSum != -100)	{ return false; }
		if( (phgcer_npeSum_min > phgcer_npeSum		|| phgcer_npeSum > phgcer_npeSum_max)	&& phgcer_npeSum != -100)	{ return false; }
		if( (petotnorm_min > petotnorm				|| petotnorm > petotnorm_max)			&& petotnorm != -100)		{ return false; }
		if( (pBetaNtrk_min > pBetaNtrk				|| pBetaNtrk > pBetaNtrk_max)			&& pBetaNtrk != -100)		{ return false; }
		return true;
	}

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//----------------------------------------------------------------------Standard Generic Cuts-----------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//---- Coincidence Time Cuts -----
	static bool cTime(double ePctime) 
	{
		Double_t ePctime_min = -3.0, ePctime_max = 3.0; //electron - Proton main coincidence time cut [ns]

		if( (ePctime_min > ePctime || ePctime > ePctime_max) && ePctime != -100) { return false; }
		return true;
	}

	static bool cTime_Rand(double ePctime_R, double ePctime_L) 
	{
		Double_t ePctime_R_min = 5.0, ePctime_R_max = 10.0; //accidental window selection(right of main coin.peak)
		Double_t ePctime_L_min = -5, ePctime_L_max = -10.0; //accidental window selection(left of main coin.peak)

		if(ePctime_R_min > ePctime_R	|| ePctime_R > ePctime_R_max)	{ return true; }
		if(ePctime_L_min > ePctime_L	|| ePctime_L > ePctime_L_max)	{ return true; }
		if (ePctime_L == -100			|| ePctime_R == -100)			{ return true; }
		return false;
	}

	//------ HMS PID Cuts ------
	static bool PID_HMS(double hetot_trkNorm, double hcer_npeSum) 
	{
		Double_t hetot_trkNorm_min = 0, hetot_trkNorm_max = 0.6;	//(HMS PID) Calorimeter Total Energy Normalized By Track Momentum
		Double_t hcer_npeSum_min = 0, hcer_npeSum_max = 0.5;		//(HMS PID) Gas Cherenkov

		if( (hetot_trkNorm_min > hetot_trkNorm		|| hetot_trkNorm > hetot_trkNorm_max)	&& hetot_trkNorm != -100)	{ return false; }
		if( (hcer_npeSum_min > hcer_npeSum			|| hcer_npeSum > hcer_npeSum_max)		&& hcer_npeSum != -100)		{ return false; }
		return true;
	}

	//------ SHMS PID Cuts ------
	static bool PID_SHMS(double petot_trkNorm, double pngcer_npeSum_pid, double phgcer_npeSum_pid) 
	{
		Double_t petot_trkNorm_min = 0.6, petot_trkNorm_max = 1.2;			//(SHMS PID) Calorimeter Total Energy Normalized By Track Momentum
		Double_t pngcer_npeSum_pid_min = 1.5, pngcer_npeSum_pid_max = 100;	//(SHMS PID) Noble Gas Cherenkov Total Number of Photo-Electrons
		Double_t phgcer_npeSum_pid_min = 0.0, phgcer_npeSum_pid_max = 1.5;	//(SHMS PID) Heavy Gas Cherenkov (Pi / K Separation) --not needed for CaFe

		if( (petot_trkNorm_min > petot_trkNorm			|| petot_trkNorm > petot_trkNorm_max)			&& petot_trkNorm != -100)		{ return false; }
		if( (pngcer_npeSum_pid_min > pngcer_npeSum_pid	|| pngcer_npeSum_pid > pngcer_npeSum_pid_max)	&& pngcer_npeSum_pid != -100)	{ return false; }
		if( (phgcer_npeSum_pid_min > phgcer_npeSum_pid	|| phgcer_npeSum_pid > phgcer_npeSum_pid_max)	&& phgcer_npeSum_pid != -100)	{ return false; } //disabled
		return true;
	}

	
	//------ HMS Acceptance Cuts ------
	static bool Accp_HMS(double hdelta, double hxptar, double hyptar, TCutG* hmsColl, double hYColl, double hXColl) 
	{
		Double_t hdelta_min = -10, hdelta_max = 10;			//hadron arm momentum acceptance, Delta			[%]
		Double_t hxptar_min = -0.06, hxptar_max = 0.06;		//hadron out - of - plane angular range, xptar	[Rad]
		Double_t hyptar_min = -0.035, hyptar_max = 0.035;	//hadron in - plane angular range, yptar		[Rad]

		if( (hdelta_min > hdelta || hdelta > hdelta_max)	&& hdelta != -100)		{ return false; }
		if( (hxptar_min > hxptar || hxptar > hxptar_max)	&& hxptar != -100)		{ return false; } //disabled
		if( (hyptar_min > hyptar || hyptar > hyptar_max)	&& hyptar != -100)		{ return false; } //disabled
		if( !hmsColl->IsInside(hYColl, hXColl)				&& hxptar == -100 && hyptar == -100)		{ return false; }
		return true;
	}

	//------ SHMS Acceptance Cuts ------
	static bool Accp_SHMS(double edelta, double exptar, double eyptar, TCutG* shmsColl, double eYColl, double eXColl) 
	{
		Double_t edelta_min = 0, edelta_max = 22;			//electron arm momentum acceptance, Delta			[%]
		Double_t exptar_min = -0.04, exptar_max = 0.04;		//electron out - of - plane angular range, xptar	[Rad]
		Double_t eyptar_min = -0.024, eyptar_max = 0.024;	//electron in - plane angular range, yptar			[Rad]

		if( (edelta_min > edelta || edelta > edelta_max)	&& edelta != -100)		{ return false; }
		if( (exptar_min > exptar || exptar > exptar_max)	&& exptar != -100)		{ return false; } //disabled
		if( (eyptar_min > eyptar || eyptar > eyptar_max)	&& eyptar != -100)		{ return false; } //disabled
		if( !shmsColl->IsInside(eYColl, eXColl)				&& exptar == -100 && eyptar == -100)	{ return false; }
		return true;
	}
	

	//------ Z-Reaction Acceptance Cuts ------
	static bool Accp_Z(double ztarDiff) 
	{
		Double_t ztarDiff_min = -2, ztarDiff_max = 2; //Z - Reaction Vertex Difference,  [cm]

		if( (ztarDiff_min > ztarDiff || ztarDiff > ztarDiff_max) && ztarDiff == -100) { return false; } //disabled
		return true;
	}

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//-------------------------------------------------------------------Cafe Specific Kinematic Cuts-------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//---- CaFe H(e,e'p) Elastics Kinematic Cuts -----
	static constexpr Double_t heep_Q2_min = 4, heep_Q2_max = 5;			//4 - Momentum Transfers	[(GeV)^2]
	static constexpr Double_t heep_W_min = 0.85, heep_W_max = 1.05;		//Invariant Mass, W			[GeV]
	static constexpr Double_t heep_xbj_min = 0.8, heep_xbj_max = 1.2;	//x - Bjorken

	static bool kinHeepSing(double heep_Q2, double heep_W, double heep_xbj)	
	{
		if ((heep_Q2_min > heep_Q2		|| heep_Q2 > heep_Q2_max)	&& heep_Q2 != -100)		{ return false; } //disabled
		if ((heep_W_min > heep_W		|| heep_W > heep_W_max)		&& heep_W != -100)		{ return false; }
		if ((heep_xbj_min > heep_xbj	|| heep_xbj > heep_xbj_max) && heep_xbj != -100)	{ return false; } //disabled
		return true;
	}
	
	static bool kinHeepCoin(double heep_Q2, double heep_W, double heep_xbj, double heep_Em, double heep_MM)	
	{
		Double_t heep_Em_min = -0.05, heep_Em_max = 0.05; //Missing Energy		[GeV]
		Double_t heep_MM_min = -0.15, heep_MM_max = 0.15; //Missing Mass, MM	[GeV]

		if( (heep_Q2_min>= heep_Q2		|| heep_Q2 > heep_Q2_max)	&& heep_Q2 != -100)		{ return false; } //disabled
		if( (heep_W_min > heep_W		|| heep_W > heep_W_max)		&& heep_W != -100)		{ return false; }
		if( (heep_xbj_min > heep_xbj	|| heep_xbj > heep_xbj_max)	&& heep_xbj != -100)	{ return false; } //disabled
		if( (heep_Em_min > heep_Em		|| heep_Em > heep_Em_max)	&& heep_Em != -100)		{ return false; } //disabled
		if( (heep_MM_min > heep_MM		|| heep_MM > heep_MM_max)	&& heep_MM != -100)		{ return false; } //disabled
		return true;
	}

	//---- CaFe-Specific A(e,e'p) Mean-Field (MF) Kinematic Cuts -----
	static bool kinMF(double MF_Q2, double MF_Pm, double Em_nuc, double MF_Em, string tgt_type) 
	{
		Double_t MF_Q2_min = 1.8, MF_Q2_max = 100;			//4 - Momentum Transfers	[(GeV)^2]
		Double_t MF_Pm_min = 0.0, MF_Pm_max = 0.250;		//Missing Momentum			[GeV]
		Double_t d2MF_Em_min = -0.02, d2MF_Em_max = 0.1;	//Missing Energy			[GeV] (applied IF - AND - ONLY - IF the deuterium target is used)
		Double_t MF_Em_min = -0.02, MF_Em_max = 0.180;		//Missing Energy			[GeV] (applied ONLY for A > 2 nuclei)

		if( (MF_Q2_min > MF_Q2		|| MF_Q2 > MF_Q2_max)		&& MF_Q2 != -100)	{ return false; }
		if( (MF_Pm_min > MF_Pm		|| MF_Pm > MF_Pm_max)		&& MF_Pm != -100)	{ return false; }
		
		if(tgt_type == "LD2"){
			if( (d2MF_Em_min > Em_nuc || Em_nuc > d2MF_Em_max) && Em_nuc != -100) { return false; }
		}
		else{
			if( (MF_Em_min > MF_Em || MF_Em > MF_Em_max) && MF_Em != -100) { return false; }
		}
		
		return true;
	}

	//---- CaFe-Specific A(e,e'p) SRC Kinematics Cuts -----
	static bool kinSRC(double SRC_Q2, double SRC_Pm, double SRC_xbj, double SRC_thrq, double Em_nuc, double Em_src, string tgt_type)	
	{
		Double_t SRC_Q2_min = 1.8, SRC_Q2_max = 100;		//4 - Momentum Transfers							[(GeV)^2]
		Double_t SRC_Pm_min = 0.350, SRC_Pm_max = 0.7;		//Missing Momentum									[GeV]
		Double_t SRC_xbj_min = 1.2, SRC_xbj_max = 100;		//x - Bjorken
		Double_t SRC_thrq_min = 0, SRC_thrq_max = 40;		//in - plane recoil(undetected) angle, theta_rq		[Deg]
		Double_t d2SRC_Em_min = -0.02, d2SRC_Em_max = 0.1;	//Missing Energy									[GeV] (applied IF - AND - ONLY - IF the deuterium target is used)

		if( (SRC_Q2_min > SRC_Q2		|| SRC_Q2 > SRC_Q2_max)		&& SRC_Q2 != -100)		{ return false; }
		if( (SRC_Pm_min > SRC_Pm		|| SRC_Pm > SRC_Pm_max)		&& SRC_Pm != -100)		{ return false; }
		if( (SRC_xbj_min > SRC_xbj		|| SRC_xbj > SRC_xbj_max)	&& SRC_xbj != -100)		{ return false; }
		if( (SRC_thrq_min > SRC_thrq	|| SRC_thrq > SRC_thrq_max)	&& SRC_thrq != -100)	{ return false; }
		
		if(tgt_type == "LD2"){
			if( (d2SRC_Em_min > Em_nuc || Em_nuc > d2SRC_Em_max) && Em_nuc != -100) { return false; }
		}
		else{
			if( (0 > Em_src || Em_nuc > Em_src) && Em_src != -100) { return false; }
		}

		return true;
	}

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//-----------------------------------------------------------------------------Any Cuts-----------------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	static bool Cut(double min1, double val1, double max1) 
	{
		if(min1 > val1 || val1 > max1) { return false; }
		return true;
	}

	static bool Cut(double min1, double val1, double max1, double min2, double val2, double max2) 
	{
		if (min1 > val1 || val1 > max1 || min2 > val2 || val2 > max2) { return false; }
		return true;
	}

	static bool Cut(double min1, double val1, double max1, double min2, double val2, double max2, double min3, double val3, double max3) 
	{
		if (min1 > val1 || val1 > max1 || min2 > val2 || val2 > max2 || min3 > val3 || val3 > max3) { return false; }
		return true;
	}

	static bool Cut(double min1, double val1, double max1, double min2, double val2, double max2, double min3, double val3, double max3, double min4, double val4, double max4) 
	{
		if (min1 > val1 || val1 > max1 || min2 > val2 || val2 > max2 || min3 > val3 || val3 > max3 || min4 > val4 || val4 > max4) { return false; }
		return true;
	}

	static bool Cut(double min1, double val1, double max1, double min2, double val2, double max2, double min3, double val3, double max3, double min4, double val4, double max4, double min5, double val5, double max5) 
	{
		if (min1 > val1 || val1 > max1 || min2 > val2 || val2 > max2 || min3 > val3 || val3 > max3 || min4 > val4 || val4 > max4 || min5 > val5 || val5 > max5) { return false; }
		return true;
	}

	static bool Cut(double min1, double val1, double max1, double min2, double val2, double max2, double min3, double val3, double max3, double min4, double val4, double max4, double min5, double val5, double max5, double min6, double val6, double max6) 
	{
		if (min1 > val1 || val1 > max1 || min2 > val2 || val2 > max2 || min3 > val3 || val3 > max3 || min4 > val4 || val4 > max4 || min5 > val5 || val5 > max5 || min6 > val6 || val6 > max6) { return false; }
		return true;
	}

	private:
};

#endif