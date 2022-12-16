void read_csv(){

  ifstream fin;
  string name;
  int rollnum,roll;
  roll=1;
  fin.open("summary_files_pass1/EmissCut_100MeV/cafe_prod_Ca40_SRC_report_summary.csv ");

  while(!fin.eof())
    {
      fin>>rollnum;

      cout << rollnum << endl;
      // if(roll==rollnum)
      //	{
      //fin>>name;
      //  cout<<name;
      //  exit(0);
	  //}
    }

}
