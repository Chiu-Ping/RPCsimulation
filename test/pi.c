#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TFile.h>

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"
#include "Plotting.hh"
#include "FundamentalConstants.hh"

using namespace Garfield;

int main(int argc,char * argv[]){
	
	TApplication app("app",&argc,argv);
	plottingEngine.SetDefaultStyle();
	
	TH1::StatOverflows(kTRUE);
	TH1F* hClusterSize=new TH1F("hClusterSize","Cluster Size",310,-10,300);
	TH1F* hClusterNum=new TH1F("hClusterNum","Cluster Number",50,0,50);
	
	const double temperature=ZeroCelsius+25.;
	const double pressure=AtmosphericPressure;
	MediumMagboltz *gas=new MediumMagboltz();
	gas->SetTemperature(temperature);
	gas->SetPressure(pressure);
	gas->SetComposition("c2f4h2",75.8,"ic4h10",4,"sf6",0.31);
	
	const double width=.1;
	SolidBox * box=new SolidBox(width/2.,0.,0.,width/2.,48.7/2.,48.7/2.);
	GeometrySimple* geo=new GeometrySimple();
	geo->AddSolid(box,gas);
	
	ComponentConstant* comp=new ComponentConstant();
	comp->SetGeometry(geo);
	comp->SetElectricField(0.,0.,0.);
	
	Sensor* sensor=new Sensor();
	sensor->AddComponent(comp);
	
	TrackHeed* track=new TrackHeed();
	track->SetSensor(sensor);
	track->SetParticle("muon");
	track->SetMomentum(120.e9);
	
	const int nEvents=100000;
	track->EnableDebugging();
	for(int i=0;i<nEvents;i++){
		if(i==1)track->DisableDebugging();
		if(i%1000==0)std::cout<<i<<"/"<<nEvents<<"\n";
		
		double x0=0.,y0=0.,z0=0.,t0=0.;
		double dx0=1.,dy0=0.,dz0=0.;
		track->NewTrack(x0,y0,z0,t0,dx0,dy0,dz0);
		double xc=0.,yc=0.,zc=0.,tc=0.;
		int nc=0;
		int ncluster=0;
		double ec=0.;
		double extra=0.;
		double esum=0.;
		double nsum=0.;
		
		while(track->GetCluster(xc,yc,zc,tc,nc,ec,extra)){
			ncluster++;
			esum+=ec;
			nsum+=nc;
			hClusterSize->Fill(nc);
		}
		hClusterNum->Fill(ncluster);
	}
	
	TCanvas* c1=new TCanvas();
	hClusterSize->GetXaxis()->SetTitle("Cluster Size");
	c1->SetLogy();
	hClusterSize->Draw();
	hClusterSize->SaveAs("clustersize.root");
	
	TCanvas* c2=new TCanvas();
	hClusterNum->GetXaxis()->SetTitle("Number of Clusters");
	hClusterNum->Draw();
	hClusterNum->SaveAs("clusternum.root");
	
	std::cout<<"Work Done"<<std::endl;
	app.Run(kTRUE);
	
}
