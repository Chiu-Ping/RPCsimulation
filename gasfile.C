#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TFile.h>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {

  // TApplication app("app", &argc, argv);
 
  const double pressure = AtmosphericPressure;
  const double temperature = ZeroCelsius + 25;
 
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("c2f4h2", 75.8, "ic4h10", 4., "sf6", 0.31);
 
  // Set the field range to be covered by the gas table. 
  const int nFields = 18;
  const double emin = 0.;
  const double emax = 72000.;
  double egap = (emax - emin)/nFields;
  
  //calculate begin
  const int ncoll = 2;

  for(int i = 0; i < nFields; i++){
		  
	  double Ef = 0.;
	  double alpha =0.;
	  double beta = 0.;
	  double eta = 0.;
	  double vdrift = 0.;
	  double vdriftY = 0.;
	  double vdriftZ = 0.;
	  double dl = 0.;
	  double dt = 0.;
	  double vxerr = 0.;
	  double vyerr = 0.;
	  double vzerr = 0.;
	  double dlerr = 0.;
	  double dterr = 0.;
	  double alphaerr = 0.;
	  double etaerr = 0.;
	  double alphatof = 0.;
	  double lor=0.;
	  double lorerr=0.;
	  
	  Ef = emin + (i + 1)*egap;
	  
	  // Run Magboltz to generate the gas table.
	  gas->RunMagboltz(Ef, 0., 0., 
	                   ncoll, false, vdriftZ, vdriftY, 
					   vdrift, dl, dt, 
					   alpha, eta, lor,
					   vzerr, vyerr, vxerr, 
					   dlerr, dterr, 
					   alphaerr, etaerr, lorerr, alphatof);
		
	  beta = alpha - eta;
	  
      //alpha readout		
	  if(i == 0){
		  ofstream outfile_gasfile_alpha("gasfile_alpha.txt",ios::out);
		  if(!outfile_gasfile_alpha){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_gasfile_alpha << Ef << " " << alpha << " " << alpha*alphaerr*0.01 <<"\n";
		  outfile_gasfile_alpha.close();  
	  }
	  else{
		  ofstream outfile_gasfile_alpha("gasfile_alpha.txt",ios::app);
		  if(!outfile_gasfile_alpha){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_gasfile_alpha << Ef << " " << alpha << " " << alpha*alphaerr*0.01 << "\n";
		  outfile_gasfile_alpha.close();
	  }
	  
	  
	  //eta readout		
	  if(i == 0){
		  ofstream outfile_gasfile_eta("gasfile_eta.txt",ios::out);
		  if(!outfile_gasfile_eta){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_gasfile_eta << Ef << " " << eta << " " << eta*etaerr*0.01 << "\n";
		  outfile_gasfile_eta.close(); 
	  }
	  else{
		  ofstream outfile_gasfile_eta("gasfile_eta.txt",ios::app);
		  if(!outfile_gasfile_eta){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_gasfile_eta << Ef << " " << eta << " " << eta*etaerr*0.01 << "\n";
		  outfile_gasfile_eta.close();
	  }
	  
	  //beta readout		
	  if(i == 0){
		  ofstream outfile_gasfile_beta("gasfile_beta.txt",ios::out);
		  if(!outfile_gasfile_beta){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_gasfile_beta << Ef << " " << beta <<"\n";
		  outfile_gasfile_beta.close();  
	  }
	  else{
		  ofstream outfile_gasfile_beta("gasfile_beta.txt",ios::app);
		  if(!outfile_gasfile_beta){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_gasfile_beta << Ef << " " << beta << "\n";
		  outfile_gasfile_beta.close();
	  }
	  
	  //vdrift readout		
	  if(i == 0){
		  ofstream outfile_gasfile_vdrift("gasfile_vdrift.txt",ios::out);
		  if(!outfile_gasfile_vdrift){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_gasfile_vdrift << Ef << " " << vdrift << " " << vdrift*vxerr*0.01 << "\n";
		  outfile_gasfile_vdrift.close();
	  }
	  else{
		   ofstream outfile_gasfile_vdrift("gasfile_vdrift.txt",ios::app);
		  if(!outfile_gasfile_vdrift){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_gasfile_vdrift << Ef << " " << vdrift << " " << vdrift*vxerr*0.01 <<"\n";
		  outfile_gasfile_vdrift.close();
	  }
	  
	  //dl,dt readout
	  if(i == 0){
		  ofstream outfile_dl_dt("dl_dt.txt",ios::out);
		  if(!outfile_dl_dt){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_dl_dt << Ef << " " << dl << " " << dt << "\n";
		  outfile_dl_dt.close();
	  }
	  else{
		   ofstream outfile_dl_dt("dl_dt.txt",ios::app);
		  if(!outfile_dl_dt){
			  cerr << "open error!" << endl;
			  exit(1);
		  }
		  outfile_dl_dt << Ef << " " << dl << " " << dt <<"\n";
		  outfile_dl_dt.close();
	  }
  }	  
std::cout<<"Work Done"<<endl;		  
  // app.Run(kTRUE);

}
