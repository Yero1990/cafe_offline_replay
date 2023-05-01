#include "../../../UTILS_CAFE/UTILS/parse_utils.h"
#include <iostream>
#include <string>
#include <vector>


// prototype code for implementing systematics studies to cafe analysis code

using namespace std;
void read_syst_csv()
{

  string csv_file = "cafe_systematics_cuts_file.csv";

  ifstream myFileStream(csv_file.c_str());

  if(!myFileStream.is_open()){
    cout << Form("File %s failed to open",csv_file.c_str()) << endl; 
  }


  // define the cut min/max values (the cut variables would already
  // exist in the baseAnalyzer code, but in this testing code, they
  // must be defined
  int ientry;
  
  Double_t c_MF_Q2_min;
  Double_t c_MF_Q2_max;
  
  Double_t c_d2MF_Em_min;
  Double_t c_d2MF_Em_max;

  Double_t c_MF_Em_min;
  Double_t c_MF_Em_max;

  Double_t c_MF_Pm_min;
  Double_t c_MF_Pm_max;
  
  Double_t c_SRC_Q2_min;
  Double_t c_SRC_Q2_max;


  Double_t c_SRC_Xbj_min;
  Double_t c_SRC_Xbj_max;

  Double_t c_SRC_thrq_min;
  Double_t c_SRC_thrq_max;

  Double_t c_SRC_Pm_min;
  Double_t c_SRC_Pm_max;

  Double_t hms_scale;
  Double_t shms_scale;
  
  string line;
  vector<string> parsed_header;
 
  
  int row_cnt = 0;
  string col_str;
  double row_value;
  bool debug = true;


  // this while loop should be outside the data event loop, such that for each instance of the
  // random cut combination, a loop over all events is done, and the corresponding integrated Pmiss
  // and binned Pmiss are extracted for each systematic cut combination 

  // NOTE: when implementing in baseAnalyzer, should consider creating a separate event loop method,
  // EventLoopSystematics(), or add a new requirement flag: e.g., analysis_type=="data_systematics" within the
  // EventLoop() methods, which would be for precisely these studies, and would be de-cluttered of the
  // other unnecessary variable / processes done inside the normal event loop
  
  
  // read the file: line by line
  while(getline(myFileStream, line)) {

    stringstream ss(line);

    // ignore comments
    if(line[0]=='#') continue;


    
    //parsed vector, whose elemetnts are the column headers
    parsed_header = parse_line(line, ','); // returns a vector of strings separated by a comma

    // read each of the cuts for each entry
    ientry = atoi(parsed_header[0].c_str()) ;

    c_MF_Q2_min = atof(parsed_header[1].c_str());
    c_MF_Q2_max = atof(parsed_header[2].c_str());

    c_d2MF_Em_min = atof(parsed_header[3].c_str());
    c_d2MF_Em_max = atof(parsed_header[4].c_str());
    
    c_MF_Em_min   = atof(parsed_header[3].c_str());
    c_MF_Em_max   = atof(parsed_header[4].c_str());

    c_MF_Pm_min  = atof(parsed_header[5].c_str());
    c_MF_Pm_max  = atof(parsed_header[6].c_str());

    c_SRC_Q2_min = atof(parsed_header[1].c_str());
    c_SRC_Q2_max = atof(parsed_header[2].c_str());

    c_SRC_Xbj_min = atof(parsed_header[7].c_str());
    c_SRC_Xbj_max = atof(parsed_header[8].c_str());

    c_SRC_thrq_min  = atof(parsed_header[9].c_str());
    c_SRC_thrq_max  = atof(parsed_header[10].c_str());

    c_SRC_Pm_min =  atof(parsed_header[11].c_str());
    c_SRC_Pm_max =  atof(parsed_header[12].c_str());

    hms_scale  =  atof(parsed_header[13].c_str());
    shms_scale =  atof(parsed_header[14].c_str());
    
    
    if(debug){
      cout << Form("ientry: %i ", ientry) << endl;
      cout << Form("Q2_min,max_mf: %.3f, %.3f ",    c_MF_Q2_min,    c_MF_Q2_max    ) << endl;
      cout << Form("Em_min,max_mf: %.3f, %.3f ",    c_MF_Em_min,    c_MF_Em_max    ) << endl;
      cout << Form("Pm_min,max_mf: %.3f, %.3f ",    c_MF_Pm_min,    c_MF_Pm_max    ) << endl;
      
      cout << Form("Q2_min,max_src: %.3f, %.3f ",   c_SRC_Q2_min,   c_SRC_Q2_max   ) << endl;
      cout << Form("Xbj_min,max_src: %.3f, %.3f ",  c_SRC_Xbj_min,  c_SRC_Xbj_max  ) << endl;
      cout << Form("thrq_min,max_src: %.3f, %.3f ", c_SRC_thrq_min, c_SRC_thrq_max ) << endl;
      cout << Form("Pm_min,max_src: %.3f, %.3f ",   c_SRC_Pm_min,   c_SRC_Pm_max   ) << endl;

      cout << Form("hms_scale, shms_scale: %.3f, %.3f ",  hms_scale, shms_scale    ) << endl;
    }
    
    row_cnt++;

    
  } // end loop over line


  
  
}
