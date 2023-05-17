void steering(){

  //string cmd=Form("root -l -q -b \"UTILS_CAFE/online_scripts/make_online_plots.cpp(%d, %d, %i, \\\"%s\\\", \\\"%s\\\", \\\"%s\\\", \\\"%s\\\", \\\"%s\\\", 1)\" ", run, evtNum, simc_exist, tgt_type.Data(), replay_type.Data(), analysis_cut.Data(), data_OutputFileName.Data(), simc_OutputFileName_rad.Data());
  //cout << cmd.c_str() << endl;

  //gSystem->Exec(cmd.c_str());

  Double_t Q2_min = 10.5;
  string cmd=Form("python ./post_analysis/special_studies/systematic_cuts_study/genRandCutsv2.py %.3f", Q2_min);
  gSystem->Exec(cmd.c_str());


}
