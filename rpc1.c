#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>

//#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "Plotting.hh"

//unit: cm ns ohm pF
//simulation of rpc with positive electrode in the above

using namespace Garfield;
using namespace std;
using namespace TMath;

struct aval_cluster{
	double z_in_gap;//distance
	int ioni_size;//initial size
	int aval_size;//number of electrons in a proceeding avalanche
	double aval_time;//when avalanche begins
	int isaval;//whether avalanche begins
	struct aval_cluster * next;
};
//create struct aval_cluster
struct aval_cluster*make_aval_cluster(double z,int clustersize,double t){
	struct aval_cluster*p;
	p=(struct aval_cluster*)malloc(sizeof(struct aval_cluster));
	p->z_in_gap=z;
	p->ioni_size=clustersize;
	p->aval_size=clustersize;
	p->aval_time=t;
	p->isaval=0;
	p->next=NULL;
	return p;
}
//write avalanche process
struct aval_cluster*write_aval_cluster(struct aval_cluster *p,double z,double clustersize,double t){
	p->z_in_gap=z;
	p->aval_time+=t;
	p->aval_size=clustersize;
	return p;
}
//judge if avalanche begins
struct aval_cluster*aval_begins(struct aval_cluster * p,int on_off){
	p->isaval=on_off;
	return p;
}
//free space
struct aval_cluster*free_aval_cluster(struct aval_cluster *pcluster){
	struct aval_cluster *p,*q=NULL;
	if(pcluster!=NULL&&pcluster->next!=NULL){
		p=pcluster->next;
		while(p!=NULL){
			q=p->next;free(p);p=q;
		}
	}
	pcluster->next=NULL;
	return(pcluster);
}

int main(int argc,char *argv[]){
	
	TApplication app("app",&argc,argv);
	plottingEngine.SetDefaultStyle();
	
    //geometry parameter unit [cm]
    const double gap =.1;
    const double thickness_bakelite = 0.149;//mylar included
    //const double thickness_insulatingfilm = 0.029;
  
    //electromagnetic parameter 
    const double permittivity_bakelite = 10.;
    //const double permittivity_insulatingfilm =  3.3;  //mylar film
    const double R = 1000.; // Unit [ohm]
    const double C = 10.; // Unit [pF]
  
    //voltage [v]
    const double voltage = 6500.;
	
	const int gasfile_nSteps = 20;
    double E_eta[gasfile_nSteps] = {0.}, gasfile_eta[gasfile_nSteps] = {0.};
    double E_vdrift[gasfile_nSteps] = {0.}, gasfile_vdrift[gasfile_nSteps] = {0.};
    double Error = 0.;
    //draw the gasfile
    TGraphErrors* gasfile_alpha_graph = new TGraphErrors("gasfile_alpha.txt","%lg %lg %lg");
    TGraphErrors* gasfile_eta_graph = new TGraphErrors("gasfile_eta.txt","%lg %lg %lg");
    TGraphErrors* gasfile_vdrift_graph = new TGraphErrors("gasfile_vdrift.txt","%lg %lg %lg");
    //read in beta gasfile
    ifstream infile_eta("gasfile_eta.txt",ios::in);
    if(!infile_eta){
	    cerr << "open error!" << endl;
	    exit(1);
    }
    for(int i_gasfile = 0; i_gasfile < 3*gasfile_nSteps; i_gasfile++){
	    if(i_gasfile%3 == 0){
		    infile_eta >> E_eta[i_gasfile/3];
	    }
	    else{
		    if(i_gasfile%3 == 1){
		    infile_eta >> gasfile_eta[(i_gasfile-1)/3];
	    }
	        else{
			    infile_eta >> Error;
		    }
	    }
    }
    infile_eta.close();
    //read in vdrift gasfile
    ifstream infile_vdrift("gasfile_vdrift.txt",ios::in);
    if(!infile_vdrift){
	    cerr << "open error!" << endl;
	    exit(1);
    }
    for(int i_gasfile = 0; i_gasfile < 3*gasfile_nSteps; i_gasfile++){
	    if(i_gasfile%3 == 0){
		    infile_vdrift >> E_vdrift[i_gasfile/3];
	    }
	    else{
		    if(i_gasfile%3 ==1){
		        infile_vdrift >> gasfile_vdrift[(i_gasfile-1)/3];
	        }
	        else{
			    infile_vdrift >> Error;
	        }
        }
    }
    infile_vdrift.close();
	TFile *clusternum = new TFile("clusternum.root");
    TH1F *den = (TH1F*) clusternum->Get("hClusterNum");
    TFile *clustersize = new TFile("clustersize.root");
    TH1F *size = (TH1F*) clustersize->Get("hClusterSize");
    
    
    //simualtion parameter
    const double tStep = 0.02;  //Unit [ns]
    const int nSteps = 4000;
    const int nEvents = 1000;
    const int N_clt = 500;
    const double nesatu = 1.6e7;
    
    //preparation
    //Unit [V/cm]
    double Ef = 0.;
    double Ew = 0.;
    Ef = voltage/gap ;
    Ew = 1/(gap + 2.*thickness_bakelite/permittivity_bakelite);
    
    //polynomial fitting of the alpha~E(3 order) 
    //linear interpolation of the eta~E vdrift~E
    //Unit [cm^-1] [cm/ns]
    double alpha =0.;
    double alpha0=0.;
    double eta = 0.;
    double beta = 0.;
    double vdrift = 0.;
    //alpha fitting
    double par_alpha[4] = {0.};
    TF1* fun_alpha = new TF1("fun_alpha", "pol3", 20000., 70000.);
    gasfile_alpha_graph->Fit("fun_alpha", "R");
    fun_alpha->GetParameters(par_alpha);
    alpha0 = par_alpha[0] + par_alpha[1]*Ef + par_alpha[2]*Ef*Ef + par_alpha[3]*Ef*Ef*Ef;
    //beta interpolation
    double E_eta_interpo[2] = {0.}, eta_interpo[2] = {0.};
    for(int i_gasfile = 0;;i_gasfile++){
	    if(E_eta[i_gasfile] < Ef){
		    E_eta_interpo[0] = E_eta[i_gasfile];
		    eta_interpo[0] = gasfile_eta[i_gasfile];
	    }
	    else{
		    E_eta_interpo[1] = E_eta[i_gasfile];
		    eta_interpo[1] = gasfile_eta[i_gasfile];
		    break;
	    }
    }
    eta = (eta_interpo[1]-eta_interpo[0])/(E_eta_interpo[1]-E_eta_interpo[0])*(Ef - E_eta_interpo[0]) + eta_interpo[0];
    double E_vdrift_interpo[2] = {0.}, vdrift_interpo[2] = {0.};
    for(int i_gasfile = 0;;i_gasfile++){
	    if(E_vdrift[i_gasfile] < Ef){
		    E_vdrift_interpo[0] = E_vdrift[i_gasfile];
		    vdrift_interpo[0] = gasfile_vdrift[i_gasfile];
	    }
	    else{
		    E_vdrift_interpo[1] = E_vdrift[i_gasfile];
		    vdrift_interpo[1] = gasfile_vdrift[i_gasfile];
		    break;
	    }
    }
    vdrift = -((vdrift_interpo[1]-vdrift_interpo[0])/(E_vdrift_interpo[1]-E_vdrift_interpo[0])*(Ef - E_vdrift_interpo[0]) + vdrift_interpo[0]);
   	cout << "Ef = " << Ef << "\n";
    cout << "Ew = " << Ew << "\n";
    cout << "alpha = " << alpha0 << "\n";
    cout << "eta = " << eta << "\n";
    cout << "vdrift = " << vdrift << "\n";
   	
   	TCanvas* c1 = new TCanvas();
    gasfile_alpha_graph->GetXaxis()->SetTitle("E [V/cm]");
    gasfile_alpha_graph->GetYaxis()->SetTitle("alpha [cm^-1]");
    gasfile_alpha_graph->Draw();
    TCanvas* c2 = new TCanvas();
    gasfile_eta_graph->GetXaxis()->SetTitle("E [V/cm]");
    gasfile_eta_graph->GetYaxis()->SetTitle("eta [cm^-1]");
    c2->SetLogy();
    gasfile_eta_graph->Draw();
    TCanvas* c3 = new TCanvas();
    gasfile_vdrift_graph->GetXaxis()->SetTitle("E [V/cm]");
    gasfile_vdrift_graph->GetYaxis()->SetTitle("vdrift [cm/ns]");
    gasfile_vdrift_graph->Draw();
   	
   	
   	
   	//array to hold results
   	int N_t[nSteps]={0};
   	double i_t[nSteps]={0.};
   	double v_t[nSteps]={0.};
   	double Q_t[nSteps]={0.};
   	
   	//avalanche parameters
   	double dStep=0.;
   	double kratio=0.;
   	double Astep=0.;
   	double stepsigma=0.;
   	double rnddis=0.;
   	double s_rnd=0.;
   	double cltmean=0.;
   	double cltsigma=0.;
   	
   	double z_0=0.;
   	double z=0.;
   	double z_drift=0.;
   	double start_0=0.;
   	double t_0=0.;
   	int t_chn=0;
   	int N_ava=0;
   	int N_this=0;
   	int N_amp=0;
   	int N_pri_electr=0;
   	int N_aval_electr1=0;
   	int N_aval_electr=0;
   	double gain=0.;
   	int issatu=0;
   	
   	int ncluster = 0, ncluster0 = 0;
    double ncluster1 = 0., ncluster2 = 0., ncluster_rnd = 0.;
    int cl_size = 0, siz0 = 0;
    double siz1 = 0., siz2 = 0., siz_rnd = 0.;
    //[ns]
    start_0 += thickness_bakelite*Sqrt(permittivity_bakelite)/SpeedOfLight;
    
    struct aval_cluster *head,*p_temp,*q_temp;
    head=(struct aval_cluster*)malloc(sizeof(struct aval_cluster));
    head->next=NULL;
    //events
    for(int i_evt=0;i_evt<nEvents;i_evt++){
    	for(int ii=0;ii<nSteps;ii++){
    		N_t[ii]=0;
    		i_t[ii]=0.;
    		v_t[ii]=0.;
    		Q_t[ii]=0.;
    	}
    	
		N_aval_electr=0;
		N_aval_electr1=0;
		N_pri_electr=0;
		
		//get cluster number
    	ncluster1 = den->GetRandom();
	    //ncluster0 = int(ncluster1);
	    //ncluster2 = ncluster1 - ncluster0;
	    //ncluster_rnd = gRandom->Rndm();
	    //if(ncluster_rnd < ncluster2) { ncluster = ncluster0 + 1;}
	    //else {ncluster = ncluster0;}
	    ncluster=ncluster1;
	    
    	q_temp=head;
        
    	for(int i_cluster=0;i_cluster<ncluster;i_cluster++){
    		z_0=(gRandom->Rndm())*gap;
    		t_0=start_0+z_0/SpeedOfLight;
    		
    		
			//get cluster size
    		siz1 = size->GetRandom();
		    //siz0 = int(siz1);
		    //siz2 = siz1 - siz0;
		    //siz_rnd = gRandom->Rndm();
		    //if(siz_rnd < siz2) { cl_size = siz0 + 1;}
		    //else {cl_size = siz0;}
		    cl_size=siz1;
		    N_pri_electr += cl_size;
		    
            //create clusters
    		p_temp=make_aval_cluster(z_0,cl_size,t_0);
    		q_temp->next=p_temp;
			q_temp=p_temp;
    		p_temp=p_temp->next;
    	}
    	
    	z=gap;
    	
    	issatu=1;//if number of electrons reach saturation (obsolete)
    	N_aval_electr1=N_pri_electr;
    	
    	//time steps
    	for(int i_tSteps=0;i_tSteps<nSteps;i_tSteps++){
    		alpha=alpha0*nesatu/(nesatu+N_aval_electr);
			beta=(alpha-eta);
	        kratio=eta/alpha;
	        dStep=tStep*vdrift;
	        Astep=Exp(-beta*dStep);
	        stepsigma=Sqrt((1+kratio)*Astep*(Astep-1)/(1-kratio));
	        rnddis=kratio*(Astep-1)/(Astep-kratio);
	        
	        z+=dStep;
	        
	        p_temp=head;q_temp=head->next;
	        //loop over clusters
	        for(int i_cluster=0;i_cluster<ncluster;i_cluster++){
	        	N_ava=q_temp->aval_size;
	        	N_this=0;
	        	//avalanche of clusters
	        	if((q_temp->aval_time)<i_tSteps*tStep&&((q_temp->z_in_gap)+dStep)>0){
	        		q_temp=aval_begins(q_temp,1);
	        		if(N_ava<N_clt&&issatu){
	        			for(int i_electron=0;i_electron<N_ava;i_electron++){
	        				s_rnd=gRandom->Rndm();
	        				if(s_rnd>=rnddis){
	        					N_amp = int(Log((Astep-kratio)*(1-s_rnd)/(1-kratio)/Astep)/
										  Log(1-(1-kratio)/(Astep-kratio)));
							    N_this += (N_amp + 1);
	        				}
	        			}
	        		}
	        		//clt
	        		else if (issatu){
	        			cltmean=N_ava*Astep;
	      			    cltsigma=Sqrt(N_ava*1.)*stepsigma;
	      			    N_amp=int(gRandom->Gaus(cltmean,cltsigma));
	      			    N_this+=N_amp;
	        		}
	        		else {
	        			issatu=0;
	        			N_this=N_ava;
	        		}
	        		N_ava=N_this;
	        		//write cluster information
	        		z_drift=dStep+q_temp->z_in_gap;
					q_temp=write_aval_cluster(q_temp,z_drift,N_ava,0);
	      	        q_temp=q_temp->next;
	        	}
	        	//clusters not involved in avalanche
	        	else{
	        		q_temp=aval_begins(q_temp,0);
	        	    q_temp=q_temp->next;
	        	}
	        }
	        p_temp=head->next;
	        N_aval_electr=0;
	        //read total electrons in gap
	        for(int i_cluster=0;i_cluster<ncluster;i_cluster++){
	        	if(p_temp!=NULL&&p_temp->isaval!=0){
	        		N_aval_electr+=p_temp->aval_size;
	        		
	        	}
	        	p_temp=p_temp->next;
	        }
	        N_t[i_tSteps]=N_aval_electr;
	        if(N_aval_electr>N_aval_electr1)N_aval_electr1=N_aval_electr;
    	}
    	
    	//save result
    	for(int i_result=0;i_result<nSteps;i_result++){
    		
    		i_t[i_result]=-N_t[i_result]*ElementaryCharge*Ew*(0-vdrift)*1.0e-6;//A
    		Q_t[i_result]=-i_t[i_result]*tStep*1.0e3;//pC
    		if(i_result>0){
    			Q_t[i_result]+=Q_t[i_result-1];
    		}
    		v_t[i_result]=i_t[i_result]*tStep*Exp((i_result*tStep*1.0e3)/(R*C));
    		if(i_result>0){
    			v_t[i_result]+=v_t[i_result-1];
    		}
    	}
    	for(int i_result=0;i_result<nSteps;i_result++){
    		v_t[i_result]*=Exp((-i_result*tStep*1.0e3)/(R*C))*1.0e3/C;//V
    	}
    	
    	gain=double(N_aval_electr1)/N_pri_electr;
    	
    	stringstream stream;
	    string filename_electron = "";
	    string filename_current = "";
	    string filename_voltage = "";
	    string filename_charge = "";
	    string i_str = "";
	    stream << i_evt+1;
	    stream >> i_str;
	    filename_electron = "./result/electron/electronnum_" + i_str + ".txt";
	    filename_current = "./result/current/signal_current_" + i_str + ".txt";
	    filename_voltage = "./result/voltage/signal_voltage_" + i_str + ".txt";
	    filename_charge = "./result/charge/signal_charge_" + i_str + ".txt";
	    
	    //electron readout out
	    ofstream outfile_electron(filename_electron.c_str(),ios::out);
	    if(!outfile_electron){
		    cerr << "open error!(electron)" << endl;
		    exit(1);
	    }
	    for(int signalbin = 0; signalbin < nSteps; signalbin++){
		    outfile_electron << start_0 + (signalbin + 0.5)*(tStep) << " " << N_t[signalbin] << "\n";
	    }
	    outfile_electron.close();
        
	    //current readout[mA]
	    ofstream outfile_current(filename_current.c_str(),ios::out);
	    if(!outfile_current){
		    cerr << "open error!(current)" << endl;
		    exit(1);
	    }
	    for(int signalbin = 0; signalbin < nSteps; signalbin++){
		    outfile_current << start_0 + (signalbin + 0.5)*(tStep) << " " << i_t[signalbin]/2*1e3 << "\n";
	    }
	    outfile_current.close();
	    
	    //voltage readout[V]
	    ofstream outfile_voltage(filename_voltage.c_str(),ios::out);
	    if(!outfile_voltage){
		    cerr << "open error!(voltage)" << endl;
		    exit(1);
	    }
	    for(int signalbin = 0; signalbin < nSteps; signalbin++){
		    outfile_voltage << start_0 + (signalbin + 0.5)*(tStep) << " " << v_t[signalbin]/2 << "\n";
	    }
	    outfile_voltage.close();
	    
	    //charge readout[pC]
	    ofstream outfile_charge(filename_charge.c_str(),ios::out);
	    if(!outfile_charge){
		    cerr << "open error!(charge)" << endl;
		    exit(1);
	    }
	    for(int signalbin = 0; signalbin < nSteps; signalbin++){
		    outfile_charge << start_0 + (signalbin + 0.5)*(tStep) << " " << Q_t[signalbin]/2 << "\n";
	    }
	    outfile_charge.close();
	
	    //gain readout
	    if(i_evt == 0){
		    ofstream outfile_gain("./result/gain/gain.txt",ios::out);
		    if(!outfile_gain){
		        cerr << "open error!(gain,evt==0)" << endl;
		        exit(1);
		    }
		    outfile_gain << gain << "\n";
		    outfile_gain.close();
	    }
	    else{
		    ofstream outfile_gain("./result/gain/gain.txt",ios::app);
		    if(!outfile_gain){
		        cerr << "open error!(gain)" << endl;
		        exit(1);
		    }
		    outfile_gain << gain << "\n";
		    outfile_gain.close();
	    }
    
    //mark of progress
	if(i_evt%100==0)cout<<i_evt<<"/"<<nEvents<<"\n";	
	//free nodes
	head=free_aval_cluster(head);	
    }
    free(head);
	cout<<"Work Done"<<endl;
	app.Run(kTRUE); 
}
