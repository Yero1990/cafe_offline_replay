//Main Calibration Code
#include "DC_calib.h"
#include "DC_calib.C"
#include <iostream>
#include <ctime>
using namespace std;

int main_calib()
{

  //prevent root from displaying graphs while executing
  gROOT->SetBatch(1);


  //measure execution time
  clock_t cl;
  cl = clock();
  
  //pid_prot -> applies cut on hodoscope beta (no track) 
  
  //template arguments
  //DC_calib obj("spec", "path/to/rootfile.root", runNUM, eventNUm, "pid_flag", "calib_mode"); pid_flag: "pid_elec" or "pid_kFALSE", calib_mode: "wire" or "card"
               
  DC_calib obj("HMS", "../../../ROOTfiles/dccalib/cafe_replay_dccalib_17134_-1.root ", 17134, -1, "pid_prot", "card");
  //DC_calib obj("SHMS", "../../../ROOTfiles/dccalib/cafe_replay_dccalib_16975_-1_dcUnCalib.root", 16975, -1, "pid_elec", "card");
  
  obj.setup_Directory();
  obj.SetPlaneNames();
  obj.GetDCLeafs();
  obj.AllocateDynamicArrays();
  obj.SetTdcOffset();
  obj.CreateHistoNames();
  obj.EventLoop("FillUncorrectedTimes");
  obj.Calculate_tZero();
  obj.EventLoop("ApplyT0Correction");
  obj.WriteTZeroParam();
  obj.WriteLookUpTable();
  obj.WriteToFile(1);  //set argument to (1) for debugging
 

  //stop clock
 cl = clock() - cl;
 cout << "execution time: " << cl/(double)CLOCKS_PER_SEC << " sec" << endl;

  return 0;
}
