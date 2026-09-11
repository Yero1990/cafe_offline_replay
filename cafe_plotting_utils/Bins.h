#ifndef Bins_H
#define Bins_H
using namespace std;

class Bins
{
	public:
	//static Double_t get_Variable(int i)	{Double_t Variable[] = {nbins_MF, xmin_MF, xmax_MF, nbins_SRC, xmin_SRC, xmax_SRC}; return Variable[i];}
	
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//------------------------------------------------------------------------PID Histograms Bins-----------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	static Double_t get_coin(int i)			{ Double_t coin[] = { 200, -45, -45, 400, -80, -80 };		return coin[i];}		//coincidence time [ns]
	//----HMS DETECTORS (FOR PID / TRACKING EFF. CHECKS)----
	static Double_t get_hcer(int i)			{ Double_t hcer[] = { 100, 0.001, 15, 100, 0.001, 15 };		return hcer[i];}		//cherenkov NPE sum
	static Double_t get_hcal(int i)			{ Double_t hcal[] = { 100, 0.001, 1.5, 100, 0.001, 1.5 };	return hcal[i];}		//calorimeter etotnorm / etottracknorm
	static Double_t get_hbeta(int i)		{ Double_t hbeta[] = { 100, 0.5, 1.5, 100, 0.5, 1.5 };		return hbeta[i];}		//calculated beta
	//----SHMS DETECTORS (FOR PID / TRACKING EFF. CHECKS)----
	static Double_t get_pngcer(int i)		{ Double_t pngcer[] = { 100, 0.001, 25, 100, 0.001, 25 };	return pngcer[i];}		//noble gas cherenkov npe sum
	static Double_t get_phgcer(int i)		{ Double_t phgcer[] = { 100, 0.001, 25, 100, 0.001, 25 };	return phgcer[i];}		//heavy gas cherenkob npe sum
	static Double_t get_pcal(int i)			{ Double_t pcal[] = { 100, 0.001, 1.5, 100, 0.001, 1.5 };	return pcal[i];}		//calorimeter etotnorm / etottracknorm
	static Double_t get_pbeta(int i)		{ Double_t pbeta[] = { 100, 0.5, 1.5, 100, 0.5, 1.5 };		return pbeta[i];}		//calculated beta

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//-------------------------------------------------------------Primary Kinematics (electron kinematics)-------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	static Double_t get_Q2(int i)			{ Double_t Q2[] = { 100, 1, 3, 100, 1, 3 };					return Q2[i];}			//4-momentum transfer [GeV]
	static Double_t get_W(int i)			{ Double_t W[] = { 100, 0.5, 1.5, 100, 0, 3 };				return W[i];}			//invariant mass [GeV]
	static Double_t get_W2(int i)			{ Double_t W2[] = { 100, -0.1, 1.5, 100, -0.1, 1.5 };		return W2[i];}			//invariant mass squared [(GeV)^2]
	static Double_t get_epsilon(int i)		{ Double_t epsilon[] = { 100, 0, 100, 100, 0, 100 };		return epsilon[i];}		//virtual photon polarization transfer
	static Double_t get_nu(int i)			{ Double_t nu[] = { 100, 0.5, 2, 100, 0.5, 2 };				return nu[i];}			//energy transfer [GeV]
	static Double_t get_omega(int i)		{ Double_t omega[] = { 100, 0, 4.5, 100, 0, 4.5 };			return omega[i];}		//????????????
	static Double_t get_phq(int i)			{ Double_t phq[] = { 100, 0, 360, 100, 0, 360 };			return phq[i];}			//out-of-plane angle between q-vector and +z (beam) [Deg]
	static Double_t get_q(int i)			{ Double_t q[] = { 100, 1, 3, 100, 1.2, 2.2 };				return q[i];}			//3-momentum transfer [GeV]
	static Double_t get_qx(int i)			{ Double_t qx[] = { 100, -5, 5, 100, -5, 5 };				return qx[i];}			//3-momentum transfer (x-component) [GeV]
	static Double_t get_qy(int i)			{ Double_t qy[] = { 100, -5, 5, 100, -5, 5 };				return qy[i];}			//3-momentum transfer (y-component) [GeV]
	static Double_t get_qz(int i)			{ Double_t qz[] = { 100, 0, 7, 100, 0, 7 };					return qz[i];}			//3-momentum transfer (z-component) [GeV]
	static Double_t get_the(int i)			{ Double_t the[] = { 100, 5, 10, 100, 5, 15 };				return the[i];}			//electron arm central angle [Deg]
	static Double_t get_thq(int i)			{ Double_t thq[] = { 100, 30, 70, 100, 30, 70 };			return thq[i];}			//in-plane angle between q-vector and +z (beam) [Deg]
	static Double_t get_X(int i)			{ Double_t X[] = { 100, 0.7, 1.2, 100, 0.1, 2.5 };			return X[i];}			//Bjorken-X
	static Double_t get_kf(int i)			{ Double_t kf[] = { 100, 6, 11, 100, 8, 11 };				return kf[i];}			//final electron arm momentum [GeV]

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//-------------------------------------------------------------------Secondary Kinematics (Hadron)------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	static Double_t get_Erecoil(int i)		{ Double_t Erecoil[] = { 100, -1, 1, 100, -1, 1 };			return Erecoil[i];}		//Total Energy of Recoil System [GeV/c]
	static Double_t get_MMk(int i)			{ Double_t MMk[] = { 100, -1, 1, 100, -1, 1 };				return MMk[i];}			//Missing Mass (kaon) [GeV]
	static Double_t get_MM(int i)			{ Double_t MM[] = { 200, -0.2, 1.2, 200, -0.2, 1.2 };		return MM[i];}			//Missing Mass (Undetected Recoil System Mass) [GeV]
	static Double_t get_MM2(int i)			{ Double_t MM2[] = { 100, -1.5, 1.5, 100, -1.5, 1.5 };		return MM2[i];}			//Missing Mass Squared [(GeV)^2]
	static Double_t get_MMpi(int i)			{ Double_t MMpi[] = { 100, -1, 1, 100, -1, 1 };				return MMpi[i];}		//Missing mass (pion) [GeV]
	static Double_t get_get_MandelS(int i)	{ Double_t MandelS[] = { 100, 0, 10, 100, 0, 10 };			return MandelS[i];}		//Mandelstam s for secondary vertex
	static Double_t get_get_MandelT(int i)	{ Double_t MandelT[] = { 100, 0, 10, 100, 0, 10 };			return MandelT[i];}		//Mandelstam t for secondary vertex
	static Double_t get_get_MandelU(int i)	{ Double_t MandelU[] = { 100, 0, 10, 100, 0, 10 };			return MandelU[i];}		//Mandelstam u for secondary vertex
	static Double_t get_Mrecoil(int i)		{ Double_t Mrecoil[] = { 100, -1, 1, 100, -1, 1 };			return Mrecoil[i];}		//Invariant Mass of Recoil System [GeV]
	static Double_t get_get_Em(int i)		{ Double_t Em[] = { 100, -0.1, 0.25, 100, -0.1, 0.25 };		return Em[i];}			//Missing Energy [GeV]
	static Double_t get_get_Em_nuc(int i)	{ Double_t Em_nuc[] = { 100, -0.1, 0.25, 100, -0.1, 0.25 }; return Em_nuc[i];}		//Nuclear Missing Energy [GeV]
	static Double_t get_Pm(int i)			{ Double_t Pm[] = { 100, -0.05, 0.5, 60, -0.05, 0.75 };		return Pm[i];}			//Missing (Recoil) Momentum [GeV/c]
	static Double_t get_get_Pmx_lab(int i)	{ Double_t Pmx_lab[] = { 80, -0.2, 0.2, 80, -0.5, -0.5 };	return Pmx_lab[i];}		//Missing Momentum (x-component, lab) [GeV]
	static Double_t get_get_Pmy_lab(int i)	{ Double_t Pmy_lab[] = { 80, -0.2, 0.2, 80, -0.5, 0.5 };	return Pmy_lab[i];}		//Missing Momentum (y-component, lab) [GeV]
	static Double_t get_get_Pmz_lab(int i)	{ Double_t Pmz_lab[] = { 80, -0.2, 0.4, 80, -0.2, 0.6 };	return Pmz_lab[i];}		//Missing Momentum (z-component, lab) [GeV]
	static Double_t get_get_Pmx_q(int i)	{ Double_t Pmx_q[] = { 80, -2, 2, 80, -2, 2 };				return Pmx_q[i];}		//Missing Momentum (x-component, q-frame: +z_lab rotated to +q) [GeV]
	static Double_t get_get_Pmy_q(int i)	{ Double_t Pmy_q[] = { 80, -2, 2, 80, -2, 2 };				return Pmy_q[i];}		//Missing Momentum (y-component, q-frame: +z_lab rotated to +q) [GeV]
	static Double_t get_get_Pmz_q(int i)	{ Double_t Pmz_q[] = { 80, -2, 3, 80, -2, 3 };				return Pmz_q[i];}		//Missing Momentum (z-component, q-frame: +z_lab rotated to +q) [GeV]
	static Double_t get_phxq(int i)			{ Double_t phxq[] = { 120, -10, 370, 100, -10, 360 };		return phxq[i];}		//Out-of-plane angle between detected particle and q [Deg]
	static Double_t get_phrq(int i)			{ Double_t phrq[] = { 120, -10, 370, 100, -10, 360 };		return phrq[i];}		//Out-of-plane angle between recoil system and q [Deg]
	static Double_t get_phxq_cm(int i)		{ Double_t phxq_cm[] = { 100, 0, 360, 120, -10, 370 };		return phxq_cm[i];}		//Out-of-plane angle between detected particle and q in CM frame [Deg]
	static Double_t get_phrq_cm(int i)		{ Double_t phrq_cm[] = { 100, 0, 360, 120, -10, 370 };		return phrq_cm[i];}		//Out-of-plane angle between recoil system and q in CM frame [Deg]
	static Double_t get_Ttot_cm(int i)		{ Double_t Ttot_cm[] = { 100, 0, 10, 100, 0, 10 };			return Ttot_cm[i];}		//total kinetic energy in CM frame [GeV]
	static Double_t get_thrq(int i)			{ Double_t thrq[] = { 120, -20, 200, 120, 0, 70 };			return thrq[i];}		//In-plane angle between the recoil system and q [Deg]
	static Double_t get_thrq_cm(int i)		{ Double_t thrq_cm[] = { 100, 0, 360, 100, 0, 360 };		return thrq_cm[i];}		//In-plane angle between the recoil system and q in CM frame [Deg]
	static Double_t get_thxq(int i)			{ Double_t thxq[] = { 120, -10, 25, 100, 0, 25 };			return thxq[i];}		//in-plane angle between detected particle and q [Deg]
	static Double_t get_thxq_cm(int i)		{ Double_t thxq_cm[] = { 100, 0, 360, 100, 0, 360 };		return thxq_cm[i];}		//in-plane angle between detected particle and q in CM frame [Deg]
	static Double_t get_thx(int i)			{ Double_t thx[] = { 100, 40, 55, 100, 60, 70 };			return thx[i];}			//hadron arm central angle [Deg]
	static Double_t get_thx_cm(int i)		{ Double_t thx_cm[] = { 100, -1, 1, 100, -1, 1 };			return thx_cm[i];}		//hadrom arm central angle [Deg]
	static Double_t get_thb_cm(int i)		{ Double_t thb_cm[] = { 100, -0.03, 10, 100, -0.03, 10 };	return thb_cm[i];}
	static Double_t get_tb(int i)			{ Double_t tb[] = { 100, 0, 10, 100, 0, 10 };				return tb[i];}
	static Double_t get_tb_cm(int i)		{ Double_t tb_cm[] = { 100, -0.01, 0.2, 100, -0.01, 0.2 };	return tb_cm[i];}
	static Double_t get_px_cm(int i)		{ Double_t px_cm[] = { 100, -1, 100, 100, -1, 100 };		return px_cm[i];}
	static Double_t get_ph_bq(int i)		{ Double_t ph_bq[] = { 100, -30, 90, 100, -30, 90 };		return ph_bq[i];}
	static Double_t get_phi_pq(int i)		{ Double_t phi_pq[] = { 100, -45, 180, 100, -45, 180 };		return phi_pq[i];}
	static Double_t get_phb_cm(int i)		{ Double_t phb_cm[] = { 100, -1, 1, 100, -1, 1 };			return phb_cm[i];}
	static Double_t get_phx_cm(int i)		{ Double_t phx_cm[] = { 100, -10, 10, 100, -10, 10 };		return phx_cm[i];}
	static Double_t get_tx(int i)			{ Double_t tx[] = { 100, 0, 5, 100, -5, 5 };				return tx[i];}			//Kinetic Energy (Detected Hadron, X) [GeV]
	static Double_t get_tx_cm(int i)		{ Double_t tx_cm[] = { 100, -5, 5, 100, -5, 10 };			return tx_cm[i];}		//Kinetic Energy (Detected Hadron, in CM frame) [GeV]
	static Double_t get_Tr(int i)			{ Double_t Tr[] = { 100, -5, 5, 100, -5, 5 };				return Tr[i];}			//Kinetic Energy (Undetected Recoil System, R) [GeV]
	static Double_t get_Tr_cm(int i)		{ Double_t Tr_cm[] = { 100, 0, 10, 100, 0, 10 };			return Tr_cm[i];}		//Kinetic Energy (Recoil System, in CM frame) [GeV]
	static Double_t get_Pf(int i)			{ Double_t Pf[] = { 100, 1.5, 2.2, 100, 1, 1.8 };			return Pf[i];}			//final hadron arm momentum [GeV/c]
	static Double_t get_xangle(int i)		{ Double_t xangle[] = { 100, -1, 1, 100, -1, 1 };			return xangle[i];}
	
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//--------------------------------------------------------------------Acceptance Histogram Bins---------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//----Electron Arm Focal Plane-----
	static Double_t get_exfp(int i)			{ Double_t exfp[] = { 100, -15, 40, 100, -15, 40 };			return exfp[i];}		//[cm]
	static Double_t get_expfp(int i)		{ Double_t expfp[] = { 100, -0.15, 0.15, 100, -0.15, 0.15 };return expfp[i];}		//[Rad]
	static Double_t get_eyfp(int i)			{ Double_t eyfp[] = { 100, -25, 15, 100, -25, 15 };			return eyfp[i];}		//[cm]
	static Double_t get_eypfp(int i)		{ Double_t eypfp[] = { 100, -0.06, 0.06, 100, -0.06, 0.06 };return eypfp[i];}		//[Rad]
	//----Electron Arm Reconstructed-----
	static Double_t get_eytar(int i)		{ Double_t eytar[] = { 100, -6, 6, 100, -6, 6 };			return eytar[i];}		//[cm]
	static Double_t get_eyptar(int i)		{ Double_t eyptar[] = { 100, -0.15, 0.15, 100, -0.15, 0.15 };return eyptar[i];}		//[Rad]
	static Double_t get_exptar(int i)		{ Double_t exptar[] = { 100, -0.15, 0.15, 100, -0.15, 0.15 };return exptar[i];}		//[Rad]
	static Double_t get_edelta(int i)		{ Double_t edelta[] = { 100, -30, 30, 100, -30, 30 };		return edelta[i];}		//[percent]
	//----Hadron Arm Focal Plane-----
	static Double_t get_hxfp(int i)			{ Double_t hxfp[] = { 100, -50, 50, 100, -50, 50 };			return hxfp[i];}		//[cm]
	static Double_t get_hxpfp(int i)		{ Double_t hxpfp[] = { 100, -0.2, 0.2, 100, -0.2, 0.2 };	return hxpfp[i];}		//[Rad]
	static Double_t get_hyfp(int i)			{ Double_t hyfp[] = { 100, -35, 35, 100, -35, 35 };			return hyfp[i];}		//[cm]
	static Double_t get_hypfp(int i)		{ Double_t hypfp[] = { 100, -0.2, 0.2, 100, -0.2, 0.2 };	return hypfp[i];}		//[Rad]
	//----Hadron Arm Reconstructed-----
	static Double_t get_hytar(int i)		{ Double_t hytar[] = { 100, -6, 6, 100, -6, 6 };			return hytar[i];}		//[cm]
	static Double_t get_hyptar(int i)		{ Double_t hyptar[] = { 100, -0.2, 0.2, 100, -0.2, 0.2 };	return hyptar[i];}		//[Rad]
	static Double_t get_hxptar(int i)		{ Double_t hxptar[] = { 100, -0.2, 0.2, 100, -0.2, 0.2 };	return hxptar[i];}		//[Rad]
	static Double_t get_hdelta(int i)		{ Double_t hdelta[] = { 100, -30, 30, 100, -30, 30 };		return hdelta[i];}		//[percent]
	//----Target Quantities----
	static Double_t get_tarx(int i)			{ Double_t tarx[] = { 100, -0.5, 0.5, 100, -0.5, 0.5 };		return tarx[i];}		//[cm]
	static Double_t get_tary(int i)			{ Double_t tary[] = { 100, -0.5, 0.5, 100, -0.5, 0.5 };		return tary[i];}		//[cm]
	static Double_t get_tarz(int i)			{ Double_t tarz[] = { 100, -15, 15, 100, -15, 15 };			return tarz[i];}		//[cm]
	static Double_t get_ztar_diff(int i)	{ Double_t ztar_diff[] = { 100, -10, 10, 100, -10, 10 };	return ztar_diff[i];}	//[percent]
	//----Collimator Quantities----
	static Double_t get_hXColl(int i)		{ Double_t hXColl[] = { 100, -15, 15, 100, -15, 15 };		return hXColl[i];}		//[cm]
	static Double_t get_hYColl(int i)		{ Double_t hYColl[] = { 100, -15, 15, 100, -15, 15 };		return hYColl[i];}		//[cm]
	static Double_t get_eXColl(int i)		{ Double_t eXColl[] = { 100, -15, 15, 100, -15, 15 };		return eXColl[i];}		//[cm]
	static Double_t get_eYColl(int i)		{ Double_t eYColl[] = { 100, -15, 15, 100, -15, 15 };		return eYColl[i];}		//[cm]

	private:
};

#endif