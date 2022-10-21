#ifndef constants_H
#define constants_H

void removeDupWord(std::string str, double *xarray)
{
  int i = 0;
  istringstream ss(str);

  string word; // for storing each word
  while (ss >> word)
  {
    // print the read word
    double number = std::stod(word);
    xarray[i] = number;
    //cout << xarray[i] << endl;
    i = i + 1;
  }
}




struct parameters{

  double xstart[6];
  double u0[3];
  double LU;
  double TU;
  double tof;
  double Tmax;
  int N;
  double mstart;
  double mu_earth;
  double Re;
  double dt;
  double J2  = 1.08262668e-3; 
  int method;
  double Isp;
  double g0;
  int meanel;
  double Rs; 
  int eclipses;
  double nu;

};

parameters  createParams(){


  parameters params;
  fstream myfile;

  myfile.open("params.txt", ios::in); // open a file to perform read operation using file object
  if (myfile.is_open())
  { // checking whether the file is open
    string tp;
   // double x0[6];

    int i = 0;
    while (getline(myfile, tp))
    { // read data from file object and put it into string
    if (i==0){
       removeDupWord(tp, params.xstart);
    }
    if (i==1){
      params.LU = stod(tp);
    }
    if (i==2){
      params.TU = stod(tp);
    }
    if (i==3){
      params.tof = stod(tp);
    }
    if (i==4){
      params.Tmax = stod(tp); //1.0/params.LU/100/1000*params.TU*params.TU;      
    }
    if (i==5){
      params.N = stoi(tp);
    }
    if (i==6){
      params.mstart = stod(tp);
    }
    if (i==7){
      params.mu_earth = stod(tp); // 3.986e5/(params.LU*params.LU*params.LU)*params.TU*params.TU;      
    }
    if (i==8){
      params.Re =stod(tp); //6378.1363/params.LU;
    }
    if (i==9){
      params.dt = stod(tp);
    }
    if (i==10){
      params.method = stoi(tp);
    }
    if (i==11){
      params.Isp = stod(tp);
    }
      if (i==12){
      params.g0 = stod(tp);
    }
      if (i==13){
      params.meanel = stoi(tp);
    }
    if (i==14){
      params.Rs = stod(tp);
    }
     if (i==15){
      params.eclipses = stoi(tp);
    }
      if (i==16){
      params.nu = stod(tp);
    }
      i++;
    }
    
    //cout <<  params.Tmax<< endl;
    myfile.close(); // close the file object.
  }


  //  cout <<  params.TU<< endl;
  //  cout <<  params.LU<< endl;
   

  return params;
}

#endif