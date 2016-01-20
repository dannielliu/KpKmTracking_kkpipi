#define KaonTrack_cxx
#include "KaonTrack.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h>
#include <TVector3.h>
#include <iostream>
#include <fstream>
using namespace std;

double FitSpectrum(TTree *dataraw, double beame, const char* namesfx=0, double *nsig=0, double *esig=0);
double FitSpectrum2(TTree *dataraw, double beame, const char* namesfx=0, double *nsig=0, double *esig=0);
double SimFit(TTree *dataraw1, TTree *dataraw2, double beame, const char* namesfx=0, double *nsig=0, double *esig=0);
double SimFit(TH1 *dataraw1, TH1 *dataraw2, double beame, const char* namesfx=0, double *nsig=0, double *esig=0);
double SimFit2(TTree *dataraw1, TTree *dataraw2, double beame, const char* namesfx=0, double *nsig=0, double *esig=0);
double SimFit3(TTree *dataraw1, TTree *dataraw2, double beame, const char* namesfx=0, double *nsig=0, double *esig=0);

class TrackingAlg
{
public:
   double mass;
   int nran;
   vector<double> pnode; // nran+1

   vector<double> kapnpran; // nran
   vector<double> kapepran;
   vector<double> kapneff;
   vector<double> kapeeff;
   
   vector<double> kamnpran; // nran
   vector<double> kamepran;
   vector<double> kamneff;
   vector<double> kameeff;
   
   vector<int> ValidRange;

   TH1D *hMmiss34 ;// = new TH1D("hMmiss34" ,"M(K)",100,0.2,0.8);
   TH1D *hMmiss4  ;// = new TH1D("hMmiss4" ,"M(K)",100,0.2,0.8);
   vector<TTree*> data3trkp; // nran
   vector<TTree*> data3trkm; // nran
   vector<TTree*> data4trkp;
   vector<TTree*> data4trkm;
   vector<TTree*> dataTtrkp;
   vector<TTree*> dataTtrkm;
   TTree* data3trkpAll;
   TTree* data3trkmAll;
   TTree* data4trkpAll;
   TTree* data4trkmAll;
   TTree* dataTtrkpAll;
   TTree* dataTtrkmAll;

   vector<TH1D*> datah3trkp; // nran
   vector<TH1D*> datah3trkm; // nran
   vector<TH1D*> datah4trkp;
   vector<TH1D*> datah4trkm;
   vector<TH1D*> datahTtrkp;
   vector<TH1D*> datahTtrkm;
   TH1D* datah3trkpAll;
   TH1D* datah3trkmAll;
   TH1D* datah4trkpAll;
   TH1D* datah4trkmAll;
   TH1D* datahTtrkpAll;
   TH1D* datahTtrkmAll;

public:
   TrackingAlg()
   {
     init();
   }
  
   ~TrackingAlg()
   {
     cout<<"Destruct tracking alg"<<endl;
     delete hMmiss34;
     delete hMmiss4;
     for (int i=0; i<nran; i++){
       delete data3trkp.at(i);
       delete data3trkm.at(i);
       delete data4trkp.at(i);
       delete data4trkm.at(i);
       delete dataTtrkp.at(i);
       delete dataTtrkm.at(i);
     }
   }
   
   void init()
   {
     cout<<"Construct tracking alg"<<endl;
     // n pt range
     nran=30;
     cout<<nran<<endl;
     for (int i=0;i<25;i++){
       pnode.push_back(0+i*0.05);
     }
     for (int i=25;i<nran+1;i++){
       pnode.push_back(1.2+(i-24)*0.1);
     }

     hMmiss34  = new TH1D("hMmiss34" ,"M(K)",100,0.2,0.8);
     hMmiss4  = new TH1D("hMmiss4" ,"M(K)",100,0.2,0.8);
     TTree* tree;
     char namei[200];
     for (int i=0; i<nran; i++){
       cout<<i<<"th tree constructing"<<endl;
       
       sprintf(namei,"data3trkp_%d",i);
       tree = new TTree(namei,namei);
       tree->Branch("mass", &mass, "mass/D");
       data3trkp.push_back(tree);
       sprintf(namei,"data3trkm_%d",i);
       tree = new TTree(namei,namei);
       tree->Branch("mass", &mass, "mass/D");
       data3trkm.push_back(tree);
       
       sprintf(namei,"data4trkp_%d",i);
       tree = new TTree(namei,namei);
       tree->Branch("mass", &mass, "mass/D");
       data4trkp.push_back(tree);
       sprintf(namei,"data4trkm_%d",i);
       tree = new TTree(namei,namei);
       tree->Branch("mass", &mass, "mass/D");
       data4trkm.push_back(tree);
       
       sprintf(namei,"dataTtrkp_%d",i); // both 3 and 4 tracks
       tree = new TTree(namei,namei);
       tree->Branch("mass", &mass, "mass/D");
       dataTtrkp.push_back(tree);
       sprintf(namei,"dataTtrkm_%d",i); // both 3 and 4 tracks
       tree = new TTree(namei,namei);
       tree->Branch("mass", &mass, "mass/D");
       dataTtrkm.push_back(tree);
     }
     sprintf(namei,"data3trkp_All");
     data3trkpAll = new TTree(namei,namei);
     data3trkpAll->Branch("mass", &mass, "mass/D");
     sprintf(namei,"data3trkm_All");
     data3trkmAll = new TTree(namei,namei);
     data3trkmAll->Branch("mass", &mass, "mass/D");
     
     sprintf(namei,"data4trkp_All");
     data4trkpAll = new TTree(namei,namei);
     data4trkpAll->Branch("mass", &mass, "mass/D");
     sprintf(namei,"data4trkm_All");
     data4trkmAll = new TTree(namei,namei);
     data4trkmAll->Branch("mass", &mass, "mass/D");
     
     sprintf(namei,"dataTtrkp"); // both 3 and 4 tracks
     dataTtrkpAll = new TTree(namei,namei);
     dataTtrkpAll->Branch("mass", &mass, "mass/D");
     sprintf(namei,"dataTtrkm"); // both 3 and 4 tracks
     dataTtrkmAll = new TTree(namei,namei);
     dataTtrkmAll->Branch("mass", &mass, "mass/D");

     // ini hist
     TH1D* hist;
     for (int i=0; i<nran; i++){
       cout<<i<<"th histogram constructing"<<endl;
       
       sprintf(namei,"datah3trkp_%d",i);
       hist = new TH1D(namei,namei,100,0.3,0.8);
       datah3trkp.push_back(hist);
       sprintf(namei,"datah3trkm_%d",i);
       hist = new TH1D(namei,namei,100,0.3,0.8);
       datah3trkm.push_back(hist);
       
       sprintf(namei,"datah4trkp_%d",i);
       hist = new TH1D(namei,namei,100,0.3,0.8);
       datah4trkp.push_back(hist);
       sprintf(namei,"datah4trkm_%d",i);
       hist = new TH1D(namei,namei,100,0.3,0.8);
       datah4trkm.push_back(hist);
       
       sprintf(namei,"datahTtrkp_%d",i); // both 3 and 4 tracks
       hist = new TH1D(namei,namei,100,0.3,0.8);
       datahTtrkp.push_back(hist);
       sprintf(namei,"datahTtrkm_%d",i); // both 3 and 4 tracks
       hist = new TH1D(namei,namei,100,0.3,0.8);
       datahTtrkm.push_back(hist);
     }
     sprintf(namei,"datah3trkp_All");
     datah3trkpAll = new TH1D(namei,namei,100,0.3,0.8);
     sprintf(namei,"datah3trkm_All");
     datah3trkmAll = new TH1D(namei,namei,100,0.3,0.8);
     
     sprintf(namei,"datah4trkp_All");
     datah4trkpAll = new TH1D(namei,namei,100,0.3,0.8);
     sprintf(namei,"datah4trkm_All");
     datah4trkmAll = new TH1D(namei,namei,100,0.3,0.8);
     
     sprintf(namei,"datahTtrkp"); // both 3 and 4 tracks
     datahTtrkpAll = new TH1D(namei,namei,100,0.3,0.8);
     sprintf(namei,"datahTtrkm"); // both 3 and 4 tracks
     datahTtrkmAll = new TH1D(namei,namei,100,0.3,0.8);

   }
 
   bool FitDataSet()
   {
     cout<<"fitting data set"<<endl;
     if (nran<1) return false;
     double nsig;
     char suffix[1000];
     double n4trk;
     double e4trk;
     double nTtrk;
     double eTtrk;
     for (int i=0; i<nran; i++)
     //int i=1;
     {
       double pave = (pnode[i+1]+pnode[i])/2;
       cout<<"\n\n\n\n#################\n"<< i << endl;
       cout<<"size of 4 trk is "<< data4trkp.at(i)->GetBranch("mass")->GetEntries()<<endl;
       cout<<"size of 3/4 trk is "<< dataTtrkp.at(i)->GetEntries()<<endl;
       sprintf(suffix,"Trk4_%.3f", pave);
       nsig = FitSpectrum2(data4trkp.at(i), pave, suffix, &n4trk, &e4trk);
       sprintf(suffix,"All_%.3f", pave);
       nsig = FitSpectrum2(dataTtrkp.at(i), pave, suffix, &nTtrk, &eTtrk);
       double TrkEff = n4trk/nTtrk;
       // uncertainty, considering correlation, 原文龙
       double num = n4trk;
       double Num = nTtrk;
     //double covnn = pow(e4trk,2);
     //double covnN = pow(e4trk,2);
     //double covNn = pow(e4trk,2);
     //double covNN = pow(eTtrk,2);
       double EffErr = 1./Num*sqrt((1-2*TrkEff)*pow(e4trk,2)+pow(TrkEff,2)*pow(eTtrk,2));
       //double EffErr = TrkEff*sqrt(pow(e4trk/n4trk,2)+pow(eTtrk/nTtrk,2));
 
       kapnpran.push_back( (pnode[i+1]+pnode[i])/2);
       kapepran.push_back( (pnode[i+1]-pnode[i])/2);
       kapneff.push_back( TrkEff);
       kapeeff.push_back( EffErr);
     }
     return true;
   }
   bool FitDataSet(int i)
   {
     cout<<"fitting data set"<<endl;
     if (nran<1) return false;
     if (nran<i) return false;
     double nsig;
     char suffix[1000];
     double n4trk;
     double e4trk;
     double nTtrk;
     double eTtrk;
     //int i=1;
     {
       double pave = (pnode[i+1]+pnode[i])/2;
       cout<<"\n\n\n\n#################\n"<< i << "\t average p "<< pave << endl;
       cout<<"size of 4 trk is "<< data4trkp.at(i)->GetBranch("mass")->GetEntries()<<endl;
       cout<<"size of 3/4 trk is "<< dataTtrkp.at(i)->GetEntries()<<endl;
       sprintf(suffix,"Trk4_%.3f", pave);
       nsig = FitSpectrum(data4trkp.at(i), pave, suffix, &n4trk, &e4trk);
       sprintf(suffix,"All_%.3f", pave);
       nsig = FitSpectrum(dataTtrkp.at(i), pave, suffix, &nTtrk, &eTtrk);
       double TrkEff = n4trk/nTtrk;
       // uncertainty, considering correlation, 原文龙
       double num = n4trk;
       double Num = nTtrk;
     //double covnn = pow(e4trk,2);
     //double covnN = pow(e4trk,2);
     //double covNn = pow(e4trk,2);
     //double covNN = pow(eTtrk,2);
       double EffErr = 1./Num*sqrt((1-2*TrkEff)*pow(e4trk,2)+pow(TrkEff,2)*pow(eTtrk,2));
       //double EffErr = TrkEff*sqrt(pow(e4trk/n4trk,2)+pow(eTtrk/nTtrk,2));
 
       kapnpran.push_back( (pnode[i+1]+pnode[i])/2);
       kapepran.push_back( (pnode[i+1]-pnode[i])/2);
       kapneff.push_back( TrkEff);
       kapeeff.push_back( EffErr);
     }
     return true;
   }
   bool SimFitDataSet()
   {
     cout<<"fitting data set"<<endl;
     if (nran<1) return false;
     double nsig;
     char suffix[1000];
     double n2[2], nerr2[2];
     double &n4trk=n2[0];
     double &e4trk=nerr2[0];
     double &nTtrk=n2[1];
     double &eTtrk=nerr2[1];
     for (int i=0; i<nran; i++)
     {
       double pave = (pnode[i+1]+pnode[i])/2;
       cout<<"\n\n\n\n#################\n"<< i << endl;
       cout<<"size of 4 trk is "<< data4trkp.at(i)->GetBranch("mass")->GetEntries()<<endl;
       cout<<"size of 3/4 trk is "<< dataTtrkp.at(i)->GetEntries()<<endl;
       sprintf(suffix,"Trkp_%.3f", pave);
       nsig = SimFit(data4trkp.at(i),dataTtrkp.at(i), pave, suffix, n2, nerr2);
       double TrkEff = n4trk/nTtrk;
       // uncertainty, considering correlation, 原文龙
       double num = n4trk;
       double Num = nTtrk;
     //double covnn = pow(e4trk,2);
     //double covnN = pow(e4trk,2);
     //double covNn = pow(e4trk,2);
     //double covNN = pow(eTtrk,2);
       double EffErr = 1./Num*sqrt((1-2*TrkEff)*pow(e4trk,2)+pow(TrkEff,2)*pow(eTtrk,2));
       //double EffErr = TrkEff*sqrt(pow(e4trk/n4trk,2)+pow(eTtrk/nTtrk,2));
 
       kapnpran.push_back( (pnode[i+1]+pnode[i])/2);
       kapepran.push_back( (pnode[i+1]-pnode[i])/2);
       kapneff.push_back( TrkEff);
       kapeeff.push_back( EffErr);
     }
     for (int i=0; i<nran; i++)
     {
       double pave = (pnode[i+1]+pnode[i])/2;
       cout<<"\n\n\n\n#################\n"<< i << endl;
       cout<<"size of 4 trk is "<< data4trkm.at(i)->GetBranch("mass")->GetEntries()<<endl;
       cout<<"size of 3/4 trk is "<< dataTtrkm.at(i)->GetEntries()<<endl;
       sprintf(suffix,"Trkm_%.3f", pave);
       nsig = SimFit(data4trkm.at(i),dataTtrkm.at(i), pave, suffix, n2, nerr2);
       double TrkEff = n4trk/nTtrk;
       // uncertainty, considering correlation, 原文龙
       double num = n4trk;
       double Num = nTtrk;
     //double covnn = pow(e4trk,2);
     //double covnN = pow(e4trk,2);
     //double covNn = pow(e4trk,2);
     //double covNN = pow(eTtrk,2);
       double EffErr = 1./Num*sqrt((1-2*TrkEff)*pow(e4trk,2)+pow(TrkEff,2)*pow(eTtrk,2));
       //double EffErr = TrkEff*sqrt(pow(e4trk/n4trk,2)+pow(eTtrk/nTtrk,2));
 
       kamnpran.push_back( (pnode[i+1]+pnode[i])/2);
       kamepran.push_back( (pnode[i+1]-pnode[i])/2);
       kamneff.push_back( TrkEff);
       kameeff.push_back( EffErr);
     }
 
     return true;
   }
   bool SimFitDataHist()
   {
     cout<<"fitting data set"<<endl;
     if (nran<1) return false;
     double nsig;
     char suffix[1000];
     double n2[2], nerr2[2];
     double &n4trk=n2[0];
     double &e4trk=nerr2[0];
     double &nTtrk=n2[1];
     double &eTtrk=nerr2[1];
     for (int i=0; i<nran; i++)
     {
       double pave = (pnode[i+1]+pnode[i])/2;
       cout<<"\n\n\n\n#################\n"<< i << endl;
       cout<<"size of 4 trk is "<< datah4trkp.at(i)->GetEntries()<<endl;
       cout<<"size of 3/4 trk is "<< datahTtrkp.at(i)->GetEntries()<<endl;
       sprintf(suffix,"Trkp_%.3f", pave);
       nsig = SimFit(datah4trkp.at(i),datahTtrkp.at(i), pave, suffix, n2, nerr2);
       double TrkEff = n4trk/nTtrk;
       // uncertainty, considering correlation, 原文龙
       double num = n4trk;
       double Num = nTtrk;
     //double covnn = pow(e4trk,2);
     //double covnN = pow(e4trk,2);
     //double covNn = pow(e4trk,2);
     //double covNN = pow(eTtrk,2);
       double EffErr = 1./Num*sqrt((1-2*TrkEff)*pow(e4trk,2)+pow(TrkEff,2)*pow(eTtrk,2));
       //double EffErr = TrkEff*sqrt(pow(e4trk/n4trk,2)+pow(eTtrk/nTtrk,2));
 
       kapnpran.push_back( (pnode[i+1]+pnode[i])/2);
       kapepran.push_back( (pnode[i+1]-pnode[i])/2);
       kapneff.push_back( TrkEff);
       kapeeff.push_back( EffErr);
     }
     for (int i=0; i<nran; i++)
     {
       double pave = (pnode[i+1]+pnode[i])/2;
       cout<<"\n\n\n\n#################\n"<< i << endl;
       cout<<"size of 4 trk is "<< datah4trkm.at(i)->GetEntries()<<endl;
       cout<<"size of 3/4 trk is "<< datahTtrkm.at(i)->GetEntries()<<endl;
       sprintf(suffix,"Trkm_%.3f", pave);
       nsig = SimFit(datah4trkm.at(i),datahTtrkm.at(i), pave, suffix, n2, nerr2);
       double TrkEff = n4trk/nTtrk;
       // uncertainty, considering correlation, 原文龙
       double num = n4trk;
       double Num = nTtrk;
     //double covnn = pow(e4trk,2);
     //double covnN = pow(e4trk,2);
     //double covNn = pow(e4trk,2);
     //double covNN = pow(eTtrk,2);
       double EffErr = 1./Num*sqrt((1-2*TrkEff)*pow(e4trk,2)+pow(TrkEff,2)*pow(eTtrk,2));
       //double EffErr = TrkEff*sqrt(pow(e4trk/n4trk,2)+pow(eTtrk/nTtrk,2));
 
       kamnpran.push_back( (pnode[i+1]+pnode[i])/2);
       kamepran.push_back( (pnode[i+1]-pnode[i])/2);
       kamneff.push_back( TrkEff);
       kameeff.push_back( EffErr);
     }
 
     return true;
   }
   bool SimFitDataSet(int i)
   {
     cout<<"fitting data set"<<endl;
     if (nran<1) return false;
     if (nran<i) return false;
     double nsig;
     char suffix[1000];
     double n2[2], nerr2[2];
     double &n4trk=n2[0];
     double &e4trk=nerr2[0];
     double &nTtrk=n2[1];
     double &eTtrk=nerr2[1];
     //int i=1;
     {
       double pave = (pnode[i+1]+pnode[i])/2;
       cout<<"\n\n\n\n#################\n"<< i << "\t average p "<< pave << endl;
       cout<<"size of 4 trk is "<< data4trkp.at(i)->GetBranch("mass")->GetEntries()<<endl;
       cout<<"size of 3/4 trk is "<< dataTtrkp.at(i)->GetEntries()<<endl;
       sprintf(suffix,"Trkp_%.3f", pave);
       nsig = SimFit(data4trkp.at(i),dataTtrkp.at(i), pave,suffix,n2,nerr2);
     }
     {
       double pave = (pnode[i+1]+pnode[i])/2;
       cout<<"\n\n\n\n#################\n"<< i << "\t average p "<< pave << endl;
       cout<<"size of 4 trk is "<< data4trkm.at(i)->GetBranch("mass")->GetEntries()<<endl;
       cout<<"size of 3/4 trk is "<< dataTtrkm.at(i)->GetEntries()<<endl;
       sprintf(suffix,"Trkm_%.3f", pave);
       nsig = SimFit(data4trkm.at(i),dataTtrkm.at(i), pave,suffix,n2,nerr2);
     }
     return true;
   }

   void ShowInHist()
   {
     char name[100];
     char drawopt[100];
     TH1D *hmKm;
     TCanvas *c1 = new TCanvas();
     for (int i=0; i<nran; i++){
       sprintf(name,"mass_%d",i);
       hmKm = new TH1D(name,name,100,0.1,0.9);
       sprintf(drawopt,"mass>>%s",name);
       data4trkp.at(i)->Draw(drawopt);
       //hmKm->Draw();
       //hmKm->Write(name);
       sprintf(name,"mass_%d.pdf",i);
       c1->Print(name);
       delete hmKm;
     }
     delete c1;
   }

   void CreateEffVsPt(TFile *file=0)
   {
     cout<<"Creating Canvas"<<endl;
     if (file!=0) file->cd();
     TCanvas *c1 = new TCanvas();
     TGraphErrors *graph = new TGraphErrors(nran,&kapnpran.at(0),&kapneff.at(0),&kapepran.at(0),&kapeeff.at(0));
     graph->GetXaxis()->SetTitle("p_{t} (GeV/c)");
     graph->GetYaxis()->SetTitle("eff");
     graph->Draw();
     c1->Write("Eff_p");
     
     TCanvas *c2 = new TCanvas();
     TGraphErrors *graphm = new TGraphErrors(nran,&kamnpran.at(0),&kamneff.at(0),&kamepran.at(0),&kameeff.at(0));
     graphm->GetXaxis()->SetTitle("p_{t} (GeV/c)");
     graphm->GetYaxis()->SetTitle("eff");
     graphm->Draw();
     c2->Write("Eff_m");

   }

};

double GetEnergy(int run)
{
  run = abs(run);
  if (run>=39335 && run<=39618) return 3.08;
  if (run>=39711 && run<=39738) return 3.02;
  if (run>=39680 && run<=39710) return 3.00;
  if (run>=39651 && run<=39679) return 2.981;
  if (run>=39619 && run<=39650) return 2.95;
  
  if (run>=39775 && run<=40069) return 2.90;
  if (run>=40128 && run<=40296) return 2.6444;
  if (run>=40300 && run<=40435) return 2.6464;
  if (run>=40436 && run<=40439) return 2.70;
  if (run>=40440 && run<=40443) return 2.80;
  if (run>=40459 && run<=40769) return 2.396;
  if (run>=40771 && run<=40776) return 2.5;
  if (run>=40777 && run<=40804) return 2.6444;//separated beam
  if (run>=40806 && run<=40951) return 2.3864;
  if (run>=40989 && run<=41121) return 2.2;
  if (run>=41122 && run<=41239) return 2.2324;
  if (run>=41240 && run<=41411) return 2.3094;
  if (run>=41416 && run<=41532) return 2.175;
  if (run>=41533 && run<=41570) return 2.15;
  if (run>=41588 && run<=41727) return 2.1;

  if (run>=41729 && run<=41909) return 2.0;
  if (run>=41911 && run<=41958) return 2.05;
  if (run>=41959 && run<=41999) return 2.2324; // separated beam
  
  if (run>=29677 && run<=30367) return 4.258;
  if (run>=31561 && run<=31981) return 4.257;
  if (run>-36773 && run<=38140) return 4.415;
  return -1;
}

void KaonTrack::LoopA(TrackingAlg* trkalg)
{
   if (fChain == 0) return;
   cout<<"Loop..."<<endl;
   Long64_t nentries = fChain->GetEntriesFast();


   fChain->GetEntry(0);
   double Ecm = GetEnergy(run);
   //TLorentzVector EPcms(Ecm*sin(0.011),0,0,Ecm);
   cout<<"Beam energy is "<< Ecm << endl;

   TLorentzVector trk[4],kapmiss,kammiss;
   double mpi = 0.13957;
   double mka = 0.493677;
   double mk0 = 0.497614;
   TLorentzVector tot4p(0.011*Ecm,0,0,Ecm);

   TH1D *hpprec = new TH1D("hpprec" ,"hpprec",200,0,2);
   TH1D *hpptru = new TH1D("hpptru" ,"hpptru",200,0,2);
   TH1D *hpmrec = new TH1D("hpmrec" ,"hpmrec",200,0,2);
   TH1D *hpmtru = new TH1D("hpmtru" ,"hpmtru",200,0,2);
   TH1D *hppdif = new TH1D("hppdif" ,"hppdif",200,-1,1);
   TH1D *hpmdif = new TH1D("hpmdif" ,"hpmdif",200,-1,1);
   
   TH1D *hpang = new TH1D("hpang","#theta_{+}",180,0,180);
   TH1D *hmang = new TH1D("hmang","#theta_{-}",180,0,180);
   
   double &mass = trkalg->mass;
   Long64_t nbytes = 0, nb = 0;
   int counter[30]={0};
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%100000 == 0) cout<<"current entry is "<< jentry << std::endl;
      
      if (indexmc>5) continue; 

      if (Ktag!=0xf) continue;
      trk[0].SetVectMag(TVector3(pippx,pippy,pippz),mpi);
      trk[1].SetVectMag(TVector3(pimpx,pimpy,pimpz),mpi);
      trk[2].SetVectMag(TVector3(kampx,kampy,kampz),mka);
         kapmiss = tot4p - (trk[0]+trk[1]+trk[2]);
      trk[3].SetVectMag(TVector3(kappx,kappy,kappz),mka);
         kammiss = tot4p - (trk[0]+trk[1]+trk[3]);
      
      double pprec = trk[3].Rho();
      double pptru = kapmiss.Rho();
      double pmrec = trk[2].Rho();
      double pmtru = kammiss.Rho();
      double ppdif = pprec - pptru;
      double pmdif = pmrec - pmtru;

      double pang = trk[3].Angle(kapmiss.Vect())*180/TMath::Pi();
      double mang = trk[2].Angle(kammiss.Vect())*180/TMath::Pi();
      double pmass = kapmiss.M();
      double mmass = kammiss.M();
      
      hpprec->Fill(pprec);
      hpptru->Fill(pptru);
      hpmrec->Fill(pmrec);
      hpmtru->Fill(pmtru);
      if (pang<5 && fabs(pmass-mka)<0.05) hppdif->Fill(ppdif);
      if (mang<5 && fabs(pmass-mka)<0.05) hpmdif->Fill(pmdif);
      
      if (fabs(ppdif)<0.05 && fabs(pmass-mka)<0.05) hpang->Fill(pang);
      if (fabs(pmdif)<0.05 && fabs(pmass-mka)<0.05) hmang->Fill(mang);
   }
   char suffix[1000];

   TCanvas *c1 = new TCanvas();
   hpprec->Draw();
   c1->Write("pprec");
   c1->Print("pprec.pdf");
   hpptru->Draw();
   c1->Write("pptru");
   c1->Print("pptru.pdf");
   hpmrec->Draw();
   c1->Write("pmrec");
   c1->Print("pmrec.pdf");
   hpmtru->Draw();
   c1->Write("pmtru");
   c1->Print("pmtru.pdf");
   hppdif->Draw();
   c1->Write("ppdif");
   c1->Print("ppdif.pdf");
   hpmdif->Draw();
   c1->Write("pmdif");
   c1->Print("pmdif.pdf");
   
   hpang->Draw();
   c1->Write("pang");
   c1->Print("pang.pdf");
   hmang->Draw();
   c1->Write("mang");
   c1->Print("mang.pdf");

   

 //for (int i=0;i<20;i++){
 //  cout<<"Part "<<i<<": "<< counter[i]<<endl;
 //}

   return;
}

void KaonTrack::LoopH(TrackingAlg* trkalg)
{
//   In a ROOT session, you can do:
//      Root > .L KaonTrack.C
//      Root > KaonTrack t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   cout<<"Loop..."<<endl;
   Long64_t nentries = fChain->GetEntriesFast();


   fChain->GetEntry(0);
   double Ecm = GetEnergy(run);
   //TLorentzVector EPcms(Ecm*sin(0.011),0,0,Ecm);
   cout<<"Beam energy is "<< Ecm << endl;

   TLorentzVector trk[4],kapmiss,kammiss;
   double mpi = 0.13957;
   double mka = 0.493677;
   double mk0 = 0.497614;
   TLorentzVector tot4p(0.011*Ecm,0,0,Ecm);

   TH1D *hm1 = new TH1D("hm1","M(#pi #pi K_{miss} #pi)",200,1,Ecm*1.2);
   TH1D *hm2 = new TH1D("hm2","M(#pi #pi K #pi)",200,1,Ecm*1.2);
   TH1D *hm3 = new TH1D("hm3","M(#pi #pi K_{miss}/K #pi)",200,1,Ecm*1.2);
   TH1D *hm4 = new TH1D("hm4","M(#pi #pi)",200,0,Ecm*1.2);
   TH1D *hth = new TH1D("hth","#theta",100,-1,1);
   TH1D *hp  = new TH1D("hp" ,"hp",200,0,2);
   TH1D *hppt  = new TH1D("hppt" ,"hppt",200,0,2);
   TH1D *hmpt  = new TH1D("hmpt" ,"hmpt",200,0,2);
   
   TH1D *hMmiss34  = trkalg->hMmiss34;
   TH1D *hMmiss4  = trkalg->hMmiss4;
   TH1D *hMrec4   = new TH1D("hMrec4" ,"M(K)",100,0.2,0.8);
   //TH1D *hPmiss34  = new TH1D("hPmiss34" ,"hp",200,0,Ecm);
   //TH1D *hPmiss4  = new TH1D("hPmiss4" ,"hp",200,0,Ecm);


   //const int nran = 8;
   //double prange[nran+1] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
 //const int nran = 7;
 //double prange[nran+1] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4};
   //double prange[nran+1] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4};
   int &nran = trkalg->nran;
   double *prange = &(trkalg->pnode.at(0));

   double &mass = trkalg->mass;
   TH1D *datah3trkp[nran];
   TH1D *datah3trkm[nran];
   TH1D *datah4trkp[nran];
   TH1D *datah4trkm[nran];
   TH1D *datahTtrkp[nran];
   TH1D *datahTtrkm[nran];
   char namei[100];
   for (int i=0; i<nran; i++){
     datah3trkp[i] = trkalg->datah3trkp.at(i);
     datah3trkm[i] = trkalg->datah3trkm.at(i);
     datah4trkp[i] = trkalg->datah4trkp.at(i);
     datah4trkm[i] = trkalg->datah4trkm.at(i);
     datahTtrkp[i] = trkalg->datahTtrkp.at(i);
     datahTtrkm[i] = trkalg->datahTtrkm.at(i);
   }
   TH1D* datah3trkpAll = trkalg->datah3trkpAll; //= new TH1D(namei,namei);
   TH1D* datah3trkmAll = trkalg->datah3trkmAll; //= new TH1D(namei,namei);
   TH1D* datah4trkpAll = trkalg->datah4trkpAll; //= new TH1D(namei,namei);
   TH1D* datah4trkmAll = trkalg->datah4trkmAll; //= new TH1D(namei,namei);
   TH1D* datahTtrkpAll = trkalg->datahTtrkpAll; //= new TH1D(namei,namei);
   TH1D* datahTtrkmAll = trkalg->datahTtrkmAll; //= new TH1D(namei,namei);

   Long64_t nbytes = 0, nb = 0;
   int counter[30]={0};
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%100000 == 0) cout<<"current entry is "<< jentry << std::endl;
      
      if (Ktag<=3) continue;
      trk[0].SetVectMag(TVector3(pippx,pippy,pippz),mpi);
      trk[1].SetVectMag(TVector3(pimpx,pimpy,pimpz),mpi);
      
    //if (Ktag&0x4) { trk[2].SetVectMag(TVector3(kappx,kappy,kappz),mka);
    //   kapmiss = tot4p - (trk[0]+trk[1]+trk[2]);}
    //if (Ktag&0x8) { trk[3].SetVectMag(TVector3(kampx,kampy,kampz),mka);
    //   kammiss = tot4p - (trk[0]+trk[1]+trk[3]);}
      
      if (Ktag&0x4) { trk[2].SetVectMag(TVector3(kampx,kampy,kampz),mka);
         kapmiss = tot4p - (trk[0]+trk[1]+trk[2]);}
      if (Ktag&0x8) { trk[3].SetVectMag(TVector3(kappx,kappy,kappz),mka);
         kammiss = tot4p - (trk[0]+trk[1]+trk[3]);}
      
      //else trk[3].SetVectMag((tot4p - (trk[0]+trk[1]+trk[2])).Vect(),mka);
      TLorentzVector trktot;// = trk[0]+trk[1]+trk[2]+trk[3];

      ///////////////////////////////////
      ////    when kap miss       ///////
      ///////////////////////////////////
      if (Ktag&0x4){
        double pkm = kapmiss.Rho();
        // transverse momentum
        double ptkm = kapmiss.Pt();
        double Coskm = kapmiss.CosTheta();
        mass = kapmiss.M();
        
        if (fabs(mass-0.493677)<0.1) hppt->Fill(ptkm);
        int idx = -1;
        for (int i=0;i<nran;i++){
          if (ptkm>prange[i] && ptkm<prange[i+1]){idx=i; break;}
        }
        if (idx>=nran || idx<0) goto kapend;
        if (fabs(Coskm)>0.93) goto kapend;
        
        hMmiss34->Fill(mass);
        //data3trkp[idx]->Fill();
        datahTtrkp[idx]->Fill(mass);
        //data3trkpAll->Fill();
        datahTtrkpAll->Fill(mass);
        if ((Ktag&0x8) || (Ktag&0x20)) 
        {
              double pprec = trk[3].Rho();
              double pptru = kapmiss.Rho();
              double ppdif = pprec - pptru;
              double pang = trk[3].Angle(kapmiss.Vect())*180/TMath::Pi();
	  if (pang>5) continue;
	  if (ppdif>0.05) continue;

              datah4trkp[idx]->Fill(mass);
              datah4trkpAll->Fill(mass);
              hMmiss4->Fill(mass);
              hMrec4->Fill(trk[3].M());
        }
      }
      kapend:
   
      ///////////////////////////////////
      ////    when kam miss       ///////
      ///////////////////////////////////
      if (Ktag&0x8){
        double pkm = kammiss.Rho();
        // transverse momentum
        double ptkm = kammiss.Pt();
        double Coskm = kammiss.CosTheta();
        mass = kammiss.M();
        if (fabs(mass-0.493677)<0.1) hmpt->Fill(ptkm);
        int idx = -1;
        for (int i=0;i<nran;i++){
          if (ptkm>prange[i] && ptkm<prange[i+1]){idx=i; break;}
        }
        if (idx>=nran || idx<0) continue;
        if (fabs(Coskm)>0.93) continue;

        //if (jentry>200) break;
        counter[idx]++;
        hMmiss34->Fill(mass);
        //data3trkp[idx]->Fill();
        datahTtrkm[idx]->Fill(mass);
        //data3trkpAll->Fill();
        datahTtrkmAll->Fill(mass);
        if ((Ktag&0x4)||(Ktag&0x10) )
        {
              double pmrec = trk[2].Rho();
              double pmtru = kammiss.Rho();
              double pmdif = pmrec - pmtru;
              double mang = trk[2].Angle(kammiss.Vect())*180/TMath::Pi();
	  if (mang>5) continue;
	  if (pmdif>0.05) continue;

              datah4trkm[idx]->Fill(mass);
              datah4trkmAll->Fill(mass);
              hMmiss4->Fill(mass);
              hMrec4->Fill(trk[2].M());
        }
      } 

   }
   char suffix[1000];
   TCanvas *c1 = new TCanvas();
 //sprintf(suffix,"Trk4_All");
 //FitSpectrum(data4trkAll, 3, suffix);
 //sprintf(suffix,"All_All");
 //FitSpectrum(dataTtrkAll, 3, suffix);
   hMmiss4->SetLineColor(2);
   hMmiss4->Draw();
   hMrec4->SetLineColor(3);
   hMrec4->Draw("same");
   c1->Write("CompreMissAndRecM");
   
   hppt->Draw();
   c1->Print("hppt.pdf");
   hmpt->Draw();
   c1->Print("hmpt.pdf");
 //hm1->Write();
 //hm2->Write();
 //hm3->Write();
 //hm4->Write();
 //hth->Write();
 //hp->Write();
   
   hMmiss34->Write();
   hMmiss4->Write();
   //hPmiss34->Write();
   //hPmiss4->Write();

   for (int i=0;i<20;i++){
     cout<<"Part "<<i<<": "<< counter[i]<<endl;
   }

   delete hm1;
   delete hm2;
   delete hm3;
   delete hm4;
   delete hth;
   delete hp;
   delete hMrec4;
   return;
}

void KaonTrack::Loop(TrackingAlg* trkalg)
{
//   In a ROOT session, you can do:
//      Root > .L KaonTrack.C
//      Root > KaonTrack t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   cout<<"Loop..."<<endl;
   Long64_t nentries = fChain->GetEntriesFast();


   fChain->GetEntry(0);
   double Ecm = GetEnergy(run);
   //TLorentzVector EPcms(Ecm*sin(0.011),0,0,Ecm);
   cout<<"Beam energy is "<< Ecm << endl;

   TLorentzVector trk[4],kapmiss,kammiss;
   double mpi = 0.13957;
   double mka = 0.493677;
   double mk0 = 0.497614;
   TLorentzVector tot4p(0.011*Ecm,0,0,Ecm);

   TH1D *hm1 = new TH1D("hm1","M(#pi #pi K_{miss} #pi)",200,1,Ecm*1.2);
   TH1D *hm2 = new TH1D("hm2","M(#pi #pi K #pi)",200,1,Ecm*1.2);
   TH1D *hm3 = new TH1D("hm3","M(#pi #pi K_{miss}/K #pi)",200,1,Ecm*1.2);
   TH1D *hm4 = new TH1D("hm4","M(#pi #pi)",200,0,Ecm*1.2);
   TH1D *hth = new TH1D("hth","#theta",100,-1,1);
   TH1D *hp  = new TH1D("hp" ,"hp",200,0,2);
   TH1D *hppt  = new TH1D("hppt" ,"hppt",200,0,2);
   TH1D *hmpt  = new TH1D("hmpt" ,"hmpt",200,0,2);
   
   TH1D *hMmiss34  = trkalg->hMmiss34;
   TH1D *hMmiss4  = trkalg->hMmiss4;
   TH1D *hMrec4   = new TH1D("hMrec4" ,"M(K)",100,0.2,0.8);
   //TH1D *hPmiss34  = new TH1D("hPmiss34" ,"hp",200,0,Ecm);
   //TH1D *hPmiss4  = new TH1D("hPmiss4" ,"hp",200,0,Ecm);


   //const int nran = 8;
   //double prange[nran+1] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
 //const int nran = 7;
 //double prange[nran+1] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4};
   //double prange[nran+1] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4};
   int &nran = trkalg->nran;
   double *prange = &(trkalg->pnode.at(0));

   double &mass = trkalg->mass;
   TTree *data3trkp[nran];
   TTree *data3trkm[nran];
   TTree *data4trkp[nran];
   TTree *data4trkm[nran];
   TTree *dataTtrkp[nran];
   TTree *dataTtrkm[nran];
   char namei[100];
   for (int i=0; i<nran; i++){
     data3trkp[i] = trkalg->data3trkp.at(i);
     data3trkm[i] = trkalg->data3trkm.at(i);
     data4trkp[i] = trkalg->data4trkp.at(i);
     data4trkm[i] = trkalg->data4trkm.at(i);
     dataTtrkp[i] = trkalg->dataTtrkp.at(i);
     dataTtrkm[i] = trkalg->dataTtrkm.at(i);
   }
   TTree* data3trkpAll = trkalg->data3trkpAll; //= new TTree(namei,namei);
   TTree* data3trkmAll = trkalg->data3trkmAll; //= new TTree(namei,namei);
   TTree* data4trkpAll = trkalg->data4trkpAll; //= new TTree(namei,namei);
   TTree* data4trkmAll = trkalg->data4trkmAll; //= new TTree(namei,namei);
   TTree* dataTtrkpAll = trkalg->dataTtrkpAll; //= new TTree(namei,namei);
   TTree* dataTtrkmAll = trkalg->dataTtrkmAll; //= new TTree(namei,namei);

   Long64_t nbytes = 0, nb = 0;
   int counter[30]={0};
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%100000 == 0) cout<<"current entry is "<< jentry << std::endl;
      
      if (Ktag<=3) continue;
      trk[0].SetVectMag(TVector3(pippx,pippy,pippz),mpi);
      trk[1].SetVectMag(TVector3(pimpx,pimpy,pimpz),mpi);
      
    //if (Ktag&0x4) { trk[2].SetVectMag(TVector3(kappx,kappy,kappz),mka);
    //   kapmiss = tot4p - (trk[0]+trk[1]+trk[2]);}
    //if (Ktag&0x8) { trk[3].SetVectMag(TVector3(kampx,kampy,kampz),mka);
    //   kammiss = tot4p - (trk[0]+trk[1]+trk[3]);}
      
      if (Ktag&0x4) { trk[2].SetVectMag(TVector3(kampx,kampy,kampz),mka);
         kapmiss = tot4p - (trk[0]+trk[1]+trk[2]);}
      if (Ktag&0x8) { trk[3].SetVectMag(TVector3(kappx,kappy,kappz),mka);
         kammiss = tot4p - (trk[0]+trk[1]+trk[3]);}
      
      //else trk[3].SetVectMag((tot4p - (trk[0]+trk[1]+trk[2])).Vect(),mka);
      TLorentzVector trktot;// = trk[0]+trk[1]+trk[2]+trk[3];

      ///////////////////////////////////
      ////    when kap miss       ///////
      ///////////////////////////////////
      if (Ktag&0x4){
        double pkm = kapmiss.Rho();
        // transverse momentum
        double ptkm = kapmiss.Pt();
        double Coskm = kapmiss.CosTheta();
        mass = kapmiss.M();
        
        if (fabs(mass-0.493677)<0.1) hppt->Fill(ptkm);
        int idx = -1;
        for (int i=0;i<nran;i++){
          if (ptkm>prange[i] && ptkm<prange[i+1]){idx=i; break;}
        }
        if (idx>=nran || idx<0) goto kapend;
        if (fabs(Coskm)>0.93) goto kapend;
       
        hMmiss34->Fill(mass);
        //data3trkp[idx]->Fill();
        dataTtrkp[idx]->Fill();
        //data3trkpAll->Fill();
        dataTtrkpAll->Fill();
        if ((Ktag&0x8) || (Ktag&0x20)) 
        {
              double pprec = trk[3].Rho();
              double pptru = kapmiss.Rho();
              double ppdif = pprec - pptru;
              double pang = trk[3].Angle(kapmiss.Vect())*180/TMath::Pi();
	  if (pang>5) continue;
	  if (ppdif>0.05) continue;

              data4trkp[idx]->Fill();
              data4trkpAll->Fill();
              hMmiss4->Fill(mass);
              hMrec4->Fill(trk[3].M());
        }
      }
      kapend:
   
      ///////////////////////////////////
      ////    when kam miss       ///////
      ///////////////////////////////////
      if (Ktag&0x8){
        double pkm = kammiss.Rho();
        // transverse momentum
        double ptkm = kammiss.Pt();
        double Coskm = kammiss.CosTheta();
        mass = kammiss.M();
        if (fabs(mass-0.493677)<0.1) hmpt->Fill(ptkm);
        int idx = -1;
        for (int i=0;i<nran;i++){
          if (ptkm>prange[i] && ptkm<prange[i+1]){idx=i; break;}
        }
        if (idx>=nran || idx<0) continue;
        if (fabs(Coskm)>0.93) continue;
             
        counter[idx]++;
        hMmiss34->Fill(mass);
        dataTtrkm[idx]->Fill();
        dataTtrkmAll->Fill();
        if ((Ktag&0x4)||(Ktag&0x10) )
        {
              double pmrec = trk[2].Rho();
              double pmtru = kammiss.Rho();
              double pmdif = pmrec - pmtru;
              double mang = trk[2].Angle(kammiss.Vect())*180/TMath::Pi();
	  if (mang>5) continue;
	  if (pmdif>0.05) continue;

              data4trkm[idx]->Fill();
              data4trkmAll->Fill();
              hMmiss4->Fill(mass);
              hMrec4->Fill(trk[2].M());
        }
      }

   }
   char suffix[1000];
   TCanvas *c1 = new TCanvas();
 //sprintf(suffix,"Trk4_All");
 //FitSpectrum(data4trkAll, 3, suffix);
 //sprintf(suffix,"All_All");
 //FitSpectrum(dataTtrkAll, 3, suffix);
   hMmiss4->SetLineColor(2);
   hMmiss4->Draw();
   hMrec4->SetLineColor(3);
   hMrec4->Draw("same");
   c1->Write("CompreMissAndRecM");
   
   hppt->Draw();
   c1->Print("hppt.pdf");
   hmpt->Draw();
   c1->Print("hmpt.pdf");
 //hm1->Write();
 //hm2->Write();
 //hm3->Write();
 //hm4->Write();
 //hth->Write();
 //hp->Write();
   
   hMmiss34->Write();
   hMmiss4->Write();
   //hPmiss34->Write();
   //hPmiss4->Write();

   for (int i=0;i<20;i++){
     cout<<"Part "<<i<<": "<< counter[i]<<endl;
   }
   
   delete hm1;
   delete hm2;
   delete hm3;
   delete hm4;
   delete hth;
   delete hp;
   delete hMrec4;
   return;
}

void KaonTrack::Loop(TrackingAlg* trkalg,int parti)
{
//   In a ROOT session, you can do:
//      Root > .L KaonTrack.C
//      Root > KaonTrack t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   cout<<"Loop..."<<endl;
   Long64_t nentries = fChain->GetEntriesFast();


   fChain->GetEntry(0);
   double Ecm = GetEnergy(run);
   //TLorentzVector EPcms(Ecm*sin(0.011),0,0,Ecm);
   cout<<"Beam energy is "<< Ecm << endl;

   TLorentzVector trk[4],kapmiss,kammiss;
   double mpi = 0.13957;
   double mka = 0.493677;
   double mk0 = 0.497614;
   TLorentzVector tot4p(0.011*Ecm,0,0,Ecm);

   TH1D *hm1 = new TH1D("hm1","M(#pi #pi K_{miss} #pi)",200,1,Ecm*1.2);
   TH1D *hm2 = new TH1D("hm2","M(#pi #pi K #pi)",200,1,Ecm*1.2);
   TH1D *hm3 = new TH1D("hm3","M(#pi #pi K_{miss}/K #pi)",200,1,Ecm*1.2);
   TH1D *hm4 = new TH1D("hm4","M(#pi #pi)",200,0,Ecm*1.2);
   TH1D *hth = new TH1D("hth","#theta",100,-1,1);
   TH1D *hp  = new TH1D("hp" ,"hp" ,200,0,2);
   TH1D *hpt = new TH1D("hpt","hpt",200,0,2);
   TH1D *hMmiss = new TH1D("hMmiss","hMmiss",100,0.1,0.9);
   
   TH1D *hMmiss34  = trkalg->hMmiss34;
   TH1D *hMmiss4  = trkalg->hMmiss4;
   TH1D *hMrec4   = new TH1D("hMrec4" ,"M(K)",100,0.2,0.8);
   
   TH1D *hprec   = new TH1D("hprec" ,"p(K+ rec)",100,0,1);
   TH1D *hpmiss   = new TH1D("hpmiss" ,"p(K+ miss)",100,0,1);
   TH1D *hpdiff   = new TH1D("hpdiff" ,"prec-pmiss",100,-1,1);
   TH1D *hang = new TH1D("hang","#theta_{-}",180,0,180);
   TH1D *hptdif   = new TH1D("hptdif" ,"ptrec-ptmiss",100,-1,1);

   //TH1D *hPmiss34  = new TH1D("hPmiss34" ,"hp",200,0,Ecm);
   //TH1D *hPmiss4  = new TH1D("hPmiss4" ,"hp",200,0,Ecm);


   //const int nran = 8;
   //double prange[nran+1] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
 //const int nran = 7;
 //double prange[nran+1] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4};
   //double prange[nran+1] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4};
   int &nran = trkalg->nran;
   double *prange = &(trkalg->pnode.at(0));

   double &mass = trkalg->mass;
   TTree *data3trkp[nran];
   TTree *data3trkm[nran];
   TTree *data4trkp[nran];
   TTree *data4trkm[nran];
   TTree *dataTtrkp[nran];
   TTree *dataTtrkm[nran];
   char namei[100];
   for (int i=0; i<nran; i++){
     data3trkp[i] = trkalg->data3trkp.at(i);
     data3trkm[i] = trkalg->data3trkm.at(i);
     data4trkp[i] = trkalg->data4trkp.at(i);
     data4trkm[i] = trkalg->data4trkm.at(i);
     dataTtrkp[i] = trkalg->dataTtrkp.at(i);
     dataTtrkm[i] = trkalg->dataTtrkm.at(i);
   }
   TTree* data3trkpAll = trkalg->data3trkpAll; //= new TTree(namei,namei);
   TTree* data3trkmAll = trkalg->data3trkmAll; //= new TTree(namei,namei);
   TTree* data4trkpAll = trkalg->data4trkpAll; //= new TTree(namei,namei);
   TTree* data4trkmAll = trkalg->data4trkmAll; //= new TTree(namei,namei);
   TTree* dataTtrkpAll = trkalg->dataTtrkpAll; //= new TTree(namei,namei);
   TTree* dataTtrkmAll = trkalg->dataTtrkmAll; //= new TTree(namei,namei);

   Long64_t nbytes = 0, nb = 0;
   int counter[30]={0};
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%100000 == 0) cout<<"current entry is "<< jentry << std::endl;
      
      if (Ktag<=3) continue;
      //if (Ktag>0x10) continue;
      trk[0].SetVectMag(TVector3(pippx,pippy,pippz),mpi);
      trk[1].SetVectMag(TVector3(pimpx,pimpy,pimpz),mpi);
      
      // Ktag : Km Kp pip pim, there is a mistake, miss store Km and Kp
    //if (Ktag&0x4) { trk[2].SetVectMag(TVector3(kappx,kappy,kappz),mka);
    //   kapmiss = tot4p - (trk[0]+trk[1]+trk[2]);}
    //if (Ktag&0x8) { trk[3].SetVectMag(TVector3(kampx,kampy,kampz),mka);
    //   kammiss = tot4p - (trk[0]+trk[1]+trk[3]);}
      if (Ktag&0x4) { trk[2].SetVectMag(TVector3(kampx,kampy,kampz),mka);
         kapmiss = tot4p - (trk[0]+trk[1]+trk[2]);}
      if (Ktag&0x8) { trk[3].SetVectMag(TVector3(kappx,kappy,kappz),mka);
         kammiss = tot4p - (trk[0]+trk[1]+trk[3]);}
      
      //else trk[3].SetVectMag((tot4p - (trk[0]+trk[1]+trk[2])).Vect(),mka);
      TLorentzVector trktot;// = trk[0]+trk[1]+trk[2]+trk[3];

      ///////////////////////////////////
      ////    when kap miss       ///////
      ///////////////////////////////////
      if (Ktag&0x4){
        double pkm = kapmiss.Rho();
        // transverse momentum
        double ptkm = kapmiss.Pt();
        double Coskm = kapmiss.CosTheta();
        hpt->Fill(ptkm);
        if (ptkm<prange[parti] || ptkm>prange[parti+1]) goto kapend;
        //if (ptkm<0.0 || ptkm>0.01) continue;
        if (fabs(Coskm)>0.93) goto kapend;
       
        mass = kapmiss.M();
        hMmiss->Fill(mass);
        hMmiss34->Fill(mass);
        //data3trkp[idx]->Fill();
        dataTtrkp[parti]->Fill();
        //data3trkpAll->Fill();
        //dataTtrkpAll->Fill();
        if ((Ktag&0x8) || (Ktag&0x20)) 
        {
              double pprec = trk[3].Rho();
              double pptru = kapmiss.Rho();
              double ppdif = pprec - pptru;
              double pang = trk[3].Angle(kapmiss.Vect())*180/TMath::Pi();
	  if (pang>5) continue;
	  if (ppdif>0.05) continue;

              data4trkp[parti]->Fill();
              //data4trkpAll->Fill();
              hMmiss4->Fill(mass);
              hMrec4->Fill(trk[3].M());
        }
      }
kapend:
   
      ///////////////////////////////////
      ////    when kam miss       ///////
      ///////////////////////////////////
      //cout<<parti<<"\t"<<prange[parti]<<"\t"<<prange[parti+1]<<endl;
      //cout<<"\t"<<jentry<<"\t"<<Ktag<<endl;
      if (Ktag&0x8){
        double pkm = kammiss.Rho();
        // transverse momentum
        double ptkm = kammiss.Pt();
        double Coskm = kammiss.CosTheta();
        //cout<<"\t"<<ptkm<<"\t tagged"<<endl;
        if (ptkm<prange[parti] || ptkm>prange[parti+1]) continue;
        if (fabs(Coskm)>0.93) continue;
             
        //if (jentry>200) break;
        counter[parti]++;
        mass = kammiss.M();
        hMmiss34->Fill(mass);
        //data3trkp[idx]->Fill();
        dataTtrkm[parti]->Fill();
        //data3trkpAll->Fill();
        //dataTtrkmAll->Fill();
        if ((Ktag&0x4)||(Ktag&0x10) )
        {
              double pmrec = trk[2].Rho();
              double pmtru = kammiss.Rho();
              double pmdif = pmrec - pmtru;
              double mang = trk[2].Angle(kammiss.Vect())*180/TMath::Pi();
	  if (mang>5) continue;
	  if (pmdif>0.05) continue;

              data4trkm[parti]->Fill();
              //data4trkmAll->Fill();
              hMmiss4->Fill(mass);
              hMrec4->Fill(trk[2].M());
	  
	  hprec->Fill(trk[2].Rho());
	  if ((mass-0.493677)>0.1) {
	    hpmiss->Fill(kammiss.Rho());
	    hpdiff->Fill(trk[2].Rho()-kammiss.Rho());
	    hptdif->Fill(trk[2].Pt()-kammiss.Pt());
	    hang->Fill(mang);
	  }
	//if (fabs(trk[2].Rho()-kammiss.Rho())<0.1) {
	//  cout<<"\t"<<Ktag<<"\t\t"<<mass<<endl;
	//}
        }
////////////    for (int i=0;i<indexmc;i++)
////////////      cout<<"\t"<<pdgid[i];
////////////      cout<<endl;
      }

   }
   char suffix[1000];
   TCanvas *c1 = new TCanvas();
 //sprintf(suffix,"Trk4_All");
 //FitSpectrum(data4trkAll, 3, suffix);
 //sprintf(suffix,"All_All");
 //FitSpectrum(dataTtrkAll, 3, suffix);
   hMmiss4->SetLineColor(2);
   hMmiss4->Draw();
   hMrec4->SetLineColor(3);
   hMrec4->Draw("same");
   c1->Write("CompreMissAndRecM");
   
   hpt->Draw();
   c1->Print("ptdis.pdf");
   hMmiss->Draw();
   c1->Print("Mmissdis.pdf");
   hprec->SetLineColor(2);
   hprec->Draw();
   hpmiss->SetLineColor(3);
   hpmiss->Draw("same");
   c1->Print("precVspmiss.pdf");
   hpdiff->Draw();
   c1->Print("pdiff.pdf");
   hptdif->Draw();
   c1->Print("ptdif.pdf");
   hang->Draw();
   c1->Print("ang.pdf");
 //hm1->Write();
 //hm2->Write();
 //hm3->Write();
 //hm4->Write();
 //hth->Write();
 //hp->Write();
   
   hMmiss34->Write();
   hMmiss4->Write();
   //hPmiss34->Write();
   //hPmiss4->Write();

   for (int i=0;i<20;i++){
     cout<<"Part "<<i<<": "<< counter[i]<<endl;
   }
   
   delete hm1;
   delete hm2;
   delete hm3;
   delete hm4;
   delete hth;
   delete hp;
   delete hMrec4;
   return;

}


int main(int argc, char** argv)
{
  TFile *ofile = new TFile("output.root","recreate");
  TrackingAlg trking;
  cout<<"KaonTrack"<<endl;
  KaonTrack* t;
  vector<char*> files;
  int partid=-1;
  //if (argc==1) t = new KaonTrack();
  if (argc==1) return -1;
  else {
    for (int i=1;i<argc;i++)
    {
      if (argv[i][0]=='-') partid=atoi(&argv[i][1]);
      else files.push_back(argv[i]);
    }
  }
  
  for (int i=0;i<files.size();i++)
  {
    cout<<"user input file "<< files.at(i)<<endl;
    TFile *f = new TFile(files.at(i));
    TTree* tree = (TTree*)f->Get("KaonTrack");
    t = new KaonTrack(tree);
    cout<<"Loop?"<<endl;
    ofile->cd();
    if (partid==-1) t->Loop(&trking);
    //if (partid==-1) t->LoopH(&trking);
    //if (partid==-1) t->LoopA(&trking);
    else t->Loop(&trking,partid);
    delete t;
  }

  //return 0;

  //trking.ShowInHist();
  if (partid!=-1) trking.SimFitDataSet(partid);
  else {
    //trking.SimFitDataHist();
    trking.SimFitDataSet();
    trking.CreateEffVsPt();
  }
  return 0;
}

#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "TPaveText.h"
#include "TMath.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
using namespace RooFit;

 // 3.08 data
 
double FitSpectrum(TTree *dataraw, double beame, const char* namesfx, double *nsig, double *esig)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double peakvalue = mka;
   double beamlow=0.3;
   double beamup=0.9;
   // try to use roofit
   RooRealVar x("mass","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma("sigma","width of gaussian",0.013+0.01*beame,0.010,0.03);
 //RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   
   RooRealVar co1("co1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2("co2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co3("co3","coefficient #2",   0 ,-100.,100.);
   RooRealVar co4("co4","coefficient #2",   0 ,-100.,100.);
   //RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   RooPolynomial ground("ground","background",x,RooArgList(co1,co2,co3,co4));
   
   RooRealVar signal("signal"," ",1000,0,10000000);//event number
   RooRealVar background("background"," ",600,0,100000000);
     
   //RooRealVar sigma2("sigma2","width of gaussian",0.01,0.008,0.02);
   RooRealVar alpha1("alpha1","#alpha",-1.0,-5.0,5.0);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
   
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_Mmiss_%s",namesfx);
   xframe = x.frame(Title("fit p"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   
   char name[100];
 //char drawopt[100];
 //TH1D *hmKm;
 //sprintf(name,"mass_%f",beame);
 //hmKm = new TH1D(name,name,100,0.1,0.9);
 //sprintf(drawopt,"mass>>%s",name);
 //dataraw->Draw(drawopt);
 //RooDataHist *datahist = new RooDataHist(tmpchr,"data",x,hmKm);
 
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape,ground),RooArgList(signal,background));
   Npar = 9;
   sum->fitTo(*dataset);
   dataset->plotOn(xframe);
   //sum->fitTo(*datahist);
   //datahist->plotOn(xframe);
   sum->plotOn(xframe,Components(cbshape),LineStyle(2),LineColor(2) );
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3)  );
   sum->plotOn(xframe  );
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.65,0.65,0.85,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.png",namesfx);
   c1->Print(tmpchr);

   ofstream fitpar("fitpar",std::ios::app);
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
   fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
   fitpar<<"\t chi: "<<xframe->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   //delete datahist;
   delete sum;
   if (nsig!=0) *nsig = signal.getVal();
   if (esig!=0) *esig = signal.getError();
   return signal.getVal();
}

double FitSpectrum2(TTree *dataraw, double beame, const char* namesfx, double *nsig, double *esig)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double peakvalue = mka;
   double beamlow=0.3;
   double beamup=0.8;
   // try to use roofit
   RooRealVar x("mass","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma("sigma","width of gaussian",0.013+0.01*beame,0.008,0.25);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   
   RooRealVar sigma2("sigma2","width of gaussian",0.03,0.025,0.04);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);

   RooRealVar co1("co1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2("co2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co3("co3","coefficient #2",   0 ,-100.,100.);
   //RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   RooPolynomial ground("ground","background",x,RooArgList(co1,co2,co3));

   RooRealVar signal("signal"," ",2000,0,10000000);//event number
   RooRealVar signal1("signal1"," ",1000,0,10000000);//event number
   RooRealVar signal2("signal2"," ",1000,0,10000000);//event number
   RooRealVar background("background"," ",600,0,100000000);
     
   RooRealVar alpha1("alpha1","#alpha",1.0,-5.0,5.0);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   RooAddPdf sig("sig","signal",RooArgList(gaus,gaus2),RooArgList(signal1,signal2));
   //RooDataHist *data_6pi;
   
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_Mmiss_%s",namesfx);
   xframe = x.frame(Title("fit p"));
  // dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   
   char name[100];
   char drawopt[100];
   TH1D *hmKm;
   sprintf(name,"mass_%f",beame);
   hmKm = new TH1D(name,name,100,0.1,0.9);
   sprintf(drawopt,"mass>>%s",name);
   dataraw->Draw(drawopt);
   RooDataHist *datahist = new RooDataHist(tmpchr,"data",x,hmKm);
 
   
   sum = new RooAddPdf("sum","sum",RooArgList(sig,ground),RooArgList(signal,background));
   Npar = 10;
   sum->fitTo(*datahist);
   //dataset->plotOn(xframe);
   datahist->plotOn(xframe);
   sum->plotOn(xframe,Components(cbshape),LineStyle(2),LineColor(2) );
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3)  );
   sum->plotOn(xframe  );
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.65,0.65,0.85,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.png",namesfx);
   c1->Print(tmpchr);

   ofstream fitpar("fitpar",std::ios::app);
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
   fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
   fitpar<<"\t chi: "<<xframe->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   //delete dataset;
   delete datahist;
   delete sum;
   if (nsig!=0) *nsig = signal.getVal();
   if (esig!=0) *esig = signal.getError();
   return signal.getVal();
}

// cbshape signal
double SimFit(TTree *dataraw1, TTree *dataraw2, double beame, const char* namesfx, double *nsig, double *esig)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double peakvalue = mka;
   double beamlow=0.3;
   double beamup=0.9;
   
   // try to use roofit
   RooRealVar x("mass","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma("sigma","width of gaussian",0.007+0.015*beame,0.005+0.015*beame,0.010+0.015*beame);
 //RooRealVar sigma("sigma","width of gaussian",0.013+0.01*beame,0.010,0.03);
 //RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooRealVar alpha("alpha","#alpha",-1.5,-5.0,5.0);
   RooRealVar nnn("n1","n",1,0.1,10);
   RooCBShape cbshape("cbshape","crystal ball",x,mean,sigma,alpha,nnn);
   
   RooRealVar co1o1("co1o1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co1o2("co1o2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co1o3("co1o3","coefficient #2",   0 ,-100.,100.);
   RooRealVar co1o4("co1o4","coefficient #2",   0 ,-100.,100.);
   RooPolynomial ground1("ground1","background",x,RooArgList(co1o1,co1o2,co1o3,co1o4));
   
   RooRealVar co2o1("co2o1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2o2("co2o2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co2o3("co2o3","coefficient #2",   0 ,-100.,100.);
   RooRealVar co2o4("co2o4","coefficient #2",   0 ,-100.,100.);
   RooPolynomial ground2("ground2","background",x,RooArgList(co2o1,co2o2,co2o3,co2o4));
   
   RooRealVar signal1("signal1"," ",1000,0,100000);//event number
   RooRealVar background1("background1"," ",100,0,100000);
     
   RooRealVar signal2("signal2"," ",1000,0,100000);//event number
   RooRealVar background2("background2"," ",100,0,100000);

   RooAddPdf model1("model1","model1",RooArgList(cbshape,ground1),RooArgList(signal1,background1));
   RooAddPdf model2("model2","model2",RooArgList(cbshape,ground2),RooArgList(signal2,background2));

   //RooAddPdf *sum;
   RooDataSet *dataset1 = new RooDataSet("data1","data1",RooArgSet(x),Import(*dataraw1));
   RooDataSet *dataset2 = new RooDataSet("data2","data2",RooArgSet(x),Import(*dataraw2));
   
   RooCategory cate("cate","cate");
   cate.defineType("ARec");
   cate.defineType("AEvt");
   
   // combined data
   RooDataSet combData("combData","combined data",x,Index(cate),Import("ARec",*dataset1),Import("AEvt",*dataset2));
   
   // construct a simultaneous pdf using category sample as index
   RooSimultaneous simPdf("simPdf","simultaneous pdf",cate);
   simPdf.addPdf(model1,"ARec");
   simPdf.addPdf(model2,"AEvt");

   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   
   simPdf.fitTo(combData);
   
   RooPlot* frame1 = x.frame(Title("All reconstructed"));
   RooPlot* frame2 = x.frame(Title("All events"));
   combData.plotOn(frame1,Cut("cate==cate::ARec"));
   combData.plotOn(frame2,Cut("cate==cate::AEvt"));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),Components("ground1"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kBlue));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),Components("cbshape"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kGreen));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),ProjWData(cate,combData));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),Components("ground2"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kBlue));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),Components("cbshape"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kGreen));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),ProjWData(cate,combData));
   
   TCanvas *c1 = new TCanvas("c1","tracking efficiency",1000,400);
   c1->Divide(2);
   c1->cd(1);gPad->SetLeftMargin(0.15);frame1->GetYaxis()->SetTitleOffset(1.4);frame1->Draw();
   c1->cd(2);gPad->SetLeftMargin(0.15);frame2->GetYaxis()->SetTitleOffset(1.4);frame2->Draw();
   
   char name[100];
   Npar = 10;
   TPaveText *pt1 = new TPaveText(0.6,0.68,0.88,0.88,"BRNDC");
   pt1->SetBorderSize(0); pt1->SetFillStyle(0);//pt1->SetFillStyle(4000);
   pt1->SetTextAlign(12); pt1->SetTextFont(42); pt1->SetTextSize(0.035);
   sprintf(name,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError()); pt1->AddText(name);
   sprintf(name,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError()); pt1->AddText(name);
   sprintf(name,"signal1 = %.1f #pm %.1f",signal1.getVal(),signal1.getError()); pt1->AddText(name);
   sprintf(name,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,frame1->chiSquare(Npar)); pt1->AddText(name);
   c1->cd(1); pt1->Draw();
 
   TPaveText *pt2 = new TPaveText(0.6,0.68,0.88,0.88,"BRNDC");
   pt2->SetBorderSize(0); pt2->SetFillStyle(0);//pt2->SetFillStyle(4000);
   pt2->SetTextAlign(12); pt2->SetTextFont(42); pt2->SetTextSize(0.035);
   sprintf(name,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError()); pt2->AddText(name);
   sprintf(name,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError()); pt2->AddText(name);
   sprintf(name,"signal2 = %.1f #pm %.1f",signal2.getVal(),signal2.getError()); pt2->AddText(name);
   sprintf(name,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,frame2->chiSquare(Npar)); pt2->AddText(name);
   c1->cd(2); pt2->Draw();
   sprintf(name,"p_spectrum_%s",namesfx);
   c1->SetName(name);
   c1->Write();
   sprintf(name,"p_spectrum_%s.png",namesfx);
   c1->Print(name);

  // c1->Print("test.pdf");
   delete pt1;
   delete pt2;
   delete c1;
   
   
   ofstream fitpar("fitpar",std::ios::app);
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<"\t sigNo = " << signal1.getVal() << "\t sigNoE = " << signal1.getError();
   fitpar<<"\t bckNo = " << background1.getVal() << "\t bckNoE = " << background1.getError();
   fitpar<<"\t chi1: "<<frame1->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   fitpar<<"\t sigNo = " << signal2.getVal() << "\t sigNoE = " << signal2.getError();
   fitpar<<"\t bckNo = " << background2.getVal() << "\t bckNoE = " << background2.getError();
   fitpar<<"\t chi2: "<<frame2->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   if (nsig!=0) {nsig[0] = signal1.getVal(); nsig[1] = signal2.getVal(); }
   if (esig!=0) {esig[0] = signal1.getError(); esig[1] = signal2.getError(); }
   delete dataset1;
   delete dataset2;
   return signal1.getVal();
}
double SimFit(TH1 *dataraw1, TH1 *dataraw2, double beame, const char* namesfx, double *nsig, double *esig)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double peakvalue = mka;
   double beamlow=0.3;
   double beamup=0.9;
   
   // try to use roofit
   RooRealVar x("mass","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma("sigma","width of gaussian",0.005+0.017*beame,0.002+0.017*beame,0.008+0.017*beame);
 //RooRealVar sigma("sigma","width of gaussian",0.013+0.01*beame,0.010,0.03);
 //RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooRealVar alpha("alpha","#alpha",-1.0,-5.0,5.0);
   RooRealVar nnn("n1","n",100,1,200);
   RooCBShape cbshape("cbshape","crystal ball",x,mean,sigma,alpha,nnn);
   
   RooRealVar co1o1("co1o1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co1o2("co1o2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co1o3("co1o3","coefficient #2",   0 ,-100.,100.);
   RooRealVar co1o4("co1o4","coefficient #2",   0 ,-100.,100.);
   RooPolynomial ground1("ground1","background",x,RooArgList(co1o1,co1o2,co1o3,co1o4));
   
   RooRealVar co2o1("co2o1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2o2("co2o2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co2o3("co2o3","coefficient #2",   0 ,-100.,100.);
   RooRealVar co2o4("co2o4","coefficient #2",   0 ,-100.,100.);
   RooPolynomial ground2("ground2","background",x,RooArgList(co2o1,co2o2,co2o3,co2o4));
   
   RooRealVar signal1("signal1"," ",1000,0,10000000);//event number
   RooRealVar background1("background1"," ",600,0,100000000);
     
   RooRealVar signal2("signal2"," ",1000,0,10000000);//event number
   RooRealVar background2("background2"," ",600,0,100000000);

   RooAddPdf model1("model1","model1",RooArgList(cbshape,ground1),RooArgList(signal1,background1));
   RooAddPdf model2("model2","model2",RooArgList(cbshape,ground2),RooArgList(signal2,background2));

   //RooAddPdf *sum;
   //RooDataSet *dataset1 = new RooDataSet("data1","data1",RooArgSet(x),Import(*dataraw1));
   //RooDataSet *dataset2 = new RooDataSet("data2","data2",RooArgSet(x),Import(*dataraw2));
   RooDataHist *dataset1 = new RooDataHist("data1","data1",RooArgSet(x),Import(*dataraw1));
   RooDataHist *dataset2 = new RooDataHist("data2","data2",RooArgSet(x),Import(*dataraw2));
   
   RooCategory cate("cate","cate");
   cate.defineType("ARec");
   cate.defineType("AEvt");
   
   // combined data
   //RooDataSet combData("combData","combined data",x,Index(cate),Import("ARec",*dataset1),Import("AEvt",*dataset2));
   RooDataHist combData("combData","combined data",x,Index(cate),Import("ARec",*dataset1),Import("AEvt",*dataset2));
   
   // construct a simultaneous pdf using category sample as index
   RooSimultaneous simPdf("simPdf","simultaneous pdf",cate);
   simPdf.addPdf(model1,"ARec");
   simPdf.addPdf(model2,"AEvt");

   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   
   simPdf.fitTo(combData);
   
   RooPlot* frame1 = x.frame(Title("All reconstructed"));
   RooPlot* frame2 = x.frame(Title("All events"));
   combData.plotOn(frame1,Cut("cate==cate::ARec"));
   combData.plotOn(frame2,Cut("cate==cate::AEvt"));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),Components("ground1"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kBlue));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),Components("cbshape"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kGreen));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),ProjWData(cate,combData));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),Components("ground2"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kBlue));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),Components("cbshape"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kGreen));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),ProjWData(cate,combData));
   
   TCanvas *c1 = new TCanvas("c1","tracking efficiency",1000,400);
   c1->Divide(2);
   c1->cd(1);gPad->SetLeftMargin(0.15);frame1->GetYaxis()->SetTitleOffset(1.4);frame1->Draw();
   c1->cd(2);gPad->SetLeftMargin(0.15);frame2->GetYaxis()->SetTitleOffset(1.4);frame2->Draw();
   
   char name[100];
   Npar = 10;
   TPaveText *pt1 = new TPaveText(0.6,0.68,0.88,0.88,"BRNDC");
   pt1->SetBorderSize(0); pt1->SetFillStyle(0);//pt1->SetFillStyle(4000);
   pt1->SetTextAlign(12); pt1->SetTextFont(42); pt1->SetTextSize(0.035);
   sprintf(name,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError()); pt1->AddText(name);
   sprintf(name,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError()); pt1->AddText(name);
   sprintf(name,"signal1 = %.1f #pm %.1f",signal1.getVal(),signal1.getError()); pt1->AddText(name);
   sprintf(name,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,frame1->chiSquare(Npar)); pt1->AddText(name);
   c1->cd(1); pt1->Draw();
 
   TPaveText *pt2 = new TPaveText(0.6,0.68,0.88,0.88,"BRNDC");
   pt2->SetBorderSize(0); pt2->SetFillStyle(0);//pt2->SetFillStyle(4000);
   pt2->SetTextAlign(12); pt2->SetTextFont(42); pt2->SetTextSize(0.035);
   sprintf(name,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError()); pt2->AddText(name);
   sprintf(name,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError()); pt2->AddText(name);
   sprintf(name,"signal2 = %.1f #pm %.1f",signal2.getVal(),signal2.getError()); pt2->AddText(name);
   sprintf(name,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,frame2->chiSquare(Npar)); pt2->AddText(name);
   c1->cd(2); pt2->Draw();
   sprintf(name,"p_spectrum_%s",namesfx);
   c1->SetName(name);
   c1->Write();
   sprintf(name,"p_spectrum_%s.png",namesfx);
   c1->Print(name);

  // c1->Print("test.pdf");
   delete pt1;
   delete pt2;
   delete c1;
   
   
   ofstream fitpar("fitpar",std::ios::app);
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<"\t sigNo = " << signal1.getVal() << "\t sigNoE = " << signal1.getError();
   fitpar<<"\t bckNo = " << background1.getVal() << "\t bckNoE = " << background1.getError();
   fitpar<<"\t chi1: "<<frame1->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   fitpar<<"\t sigNo = " << signal2.getVal() << "\t sigNoE = " << signal2.getError();
   fitpar<<"\t bckNo = " << background2.getVal() << "\t bckNoE = " << background2.getError();
   fitpar<<"\t chi2: "<<frame2->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   if (nsig!=0) {nsig[0] = signal1.getVal(); nsig[1] = signal2.getVal(); }
   if (esig!=0) {esig[0] = signal1.getError(); esig[1] = signal2.getError(); }
   delete dataset1;
   delete dataset2;
   return signal1.getVal();
}


// 2 gaus signal
double SimFit2(TTree *dataraw1, TTree *dataraw2, double beame, const char* namesfx, double *nsig, double *esig)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double peakvalue = mka;
   double beamlow=0.1;
   double beamup=1.0;
   
   // try to use roofit
   RooRealVar x("mass","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma("sigma","width of gaussian",0.0141+0.01*beame,0.006,0.03);
   RooGaussian gaus1("gaus1","gauss(x,m,s)",x,mean,sigma);
   //RooRealVar mean2("mean2","mean of gaussian",peakvalue,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma2("sigma2","width of gaussian",0.025+0.01*beame,0.020,0.09);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar sigf("sigf","signal1 fraction",0.8,0.5,1);
   RooAddPdf sig("sig","sig",RooArgList(gaus1,gaus2),sigf);
 //RooRealVar alpha("alpha","#alpha",-1.0,-5.0,5.0);
 //RooRealVar nnn("n1","n",100,1,200);
 //RooCBShape cbshape("cbshape","crystal ball",x,mean,sigma,alpha,nnn);
   
   RooRealVar co1o1("co1o1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co1o2("co1o2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co1o3("co1o3","coefficient #2",   0 ,-100.,100.);
   RooRealVar co1o4("co1o4","coefficient #2",   0 ,-100.,100.);
   RooPolynomial ground1("ground1","background",x,RooArgList(co1o1,co1o2,co1o3,co1o4));
   
   RooRealVar co2o1("co2o1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2o2("co2o2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co2o3("co2o3","coefficient #2",   0 ,-100.,100.);
   RooRealVar co2o4("co2o4","coefficient #2",   0 ,-100.,100.);
   RooPolynomial ground2("ground2","background",x,RooArgList(co2o1,co2o2,co2o3,co2o4));
   
   RooRealVar signal1("signal1"," ",1000,0,10000000);//event number
   RooRealVar background1("background1"," ",600,0,100000000);
     
   RooRealVar signal2("signal2"," ",1000,0,10000000);//event number
   RooRealVar background2("background2"," ",600,0,100000000);

   RooAddPdf model1("model1","model1",RooArgList(sig,ground1),RooArgList(signal1,background1));
   RooAddPdf model2("model2","model2",RooArgList(sig,ground2),RooArgList(signal2,background2));

   //RooAddPdf *sum;
   RooDataSet *dataset1 = new RooDataSet("data1","data1",RooArgSet(x),Import(*dataraw1));
   RooDataSet *dataset2 = new RooDataSet("data2","data2",RooArgSet(x),Import(*dataraw2));
   
   RooCategory cate("cate","cate");
   cate.defineType("ARec");
   cate.defineType("AEvt");
   
   // combined data
   RooDataSet combData("combData","combined data",x,Index(cate),Import("ARec",*dataset1),Import("AEvt",*dataset2));
   
   // construct a simultaneous pdf using category sample as index
   RooSimultaneous simPdf("simPdf","simultaneous pdf",cate);
   simPdf.addPdf(model1,"ARec");
   simPdf.addPdf(model2,"AEvt");

   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   
   simPdf.fitTo(combData);
   
   RooPlot* frame1 = x.frame(Title("All reconstructed"));
   RooPlot* frame2 = x.frame(Title("All events"));
   combData.plotOn(frame1,Cut("cate==cate::ARec"));
   combData.plotOn(frame2,Cut("cate==cate::AEvt"));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),Components("ground1"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kBlue));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),Components("sig"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kGreen));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),ProjWData(cate,combData));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),Components("ground2"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kBlue));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),Components("sig"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kGreen));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),ProjWData(cate,combData));
   
   TCanvas *c1 = new TCanvas("c1","tracking efficiency",1000,400);
   c1->Divide(2);
   c1->cd(1);gPad->SetLeftMargin(0.15);frame1->GetYaxis()->SetTitleOffset(1.4);frame1->Draw();
   c1->cd(2);gPad->SetLeftMargin(0.15);frame2->GetYaxis()->SetTitleOffset(1.4);frame2->Draw();
   
   char name[100];
   Npar = 10;
   TPaveText *pt1 = new TPaveText(0.6,0.68,0.88,0.88,"BRNDC");
   pt1->SetBorderSize(0); pt1->SetFillStyle(4000);
   pt1->SetTextAlign(12); pt1->SetTextFont(42); pt1->SetTextSize(0.035);
   sprintf(name,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError()); pt1->AddText(name);
   sprintf(name,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError()); pt1->AddText(name);
   sprintf(name,"signal1 = %.1f #pm %.1f",signal1.getVal(),signal1.getError()); pt1->AddText(name);
   sprintf(name,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,frame1->chiSquare(Npar)); pt1->AddText(name);
   c1->cd(1); pt1->Draw();
 
   TPaveText *pt2 = new TPaveText(0.6,0.68,0.88,0.88,"BRNDC");
   pt2->SetBorderSize(0); pt2->SetFillStyle(4000);
   pt2->SetTextAlign(12); pt2->SetTextFont(42); pt2->SetTextSize(0.035);
   sprintf(name,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError()); pt2->AddText(name);
   sprintf(name,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError()); pt2->AddText(name);
   sprintf(name,"signal2 = %.1f #pm %.1f",signal2.getVal(),signal2.getError()); pt2->AddText(name);
   sprintf(name,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,frame2->chiSquare(Npar)); pt2->AddText(name);
   c1->cd(2); pt2->Draw();
   sprintf(name,"p_spectrum_%s",namesfx);
   c1->SetName(name);
   c1->Write();
   sprintf(name,"p_spectrum_%s.png",namesfx);
   c1->Print(name);

  // c1->Print("test.pdf");
   delete pt1;
   delete pt2;
   delete c1;
   
   
   ofstream fitpar("fitpar",std::ios::app);
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<"\t sigNo = " << signal1.getVal() << "\t sigNoE = " << signal1.getError();
   fitpar<<"\t bckNo = " << background1.getVal() << "\t bckNoE = " << background1.getError();
   fitpar<<"\t chi1: "<<frame1->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   fitpar<<"\t sigNo = " << signal2.getVal() << "\t sigNoE = " << signal2.getError();
   fitpar<<"\t bckNo = " << background2.getVal() << "\t bckNoE = " << background2.getError();
   fitpar<<"\t chi2: "<<frame2->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   if (nsig!=0) {nsig[0] = signal1.getVal(); nsig[1] = signal2.getVal(); }
   if (esig!=0) {esig[0] = signal1.getError(); esig[1] = signal2.getError(); }
   delete dataset1;
   delete dataset2;
   return signal1.getVal();
}

#include "RooFormulaVar.h"
// cbshape signal, efficiency as parameter
double SimFit3(TTree *dataraw1, TTree *dataraw2, double beame, const char* namesfx, double *nsig, double *esig)
{ 
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double peakvalue = mka;
   double beamlow=0.3;
   double beamup=0.9;
   
   // try to use roofit
   RooRealVar x("mass","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma("sigma","width of gaussian",0.013+0.01*beame,0.010,0.03);
 //RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooRealVar alpha("alpha","#alpha",-1.0,-5.0,5.0);
   RooRealVar nnn("n1","n",100,1,200);
   RooCBShape cbshape("cbshape","crystal ball",x,mean,sigma,alpha,nnn);
   
   RooRealVar co1o1("co1o1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co1o2("co1o2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co1o3("co1o3","coefficient #2",   0 ,-100.,100.);
   RooRealVar co1o4("co1o4","coefficient #2",   0 ,-100.,100.);
   RooPolynomial ground1("ground1","background",x,RooArgList(co1o1,co1o2,co1o3,co1o4));
   
   RooRealVar co2o1("co2o1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2o2("co2o2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co2o3("co2o3","coefficient #2",   0 ,-100.,100.);
   RooRealVar co2o4("co2o4","coefficient #2",   0 ,-100.,100.);
   RooPolynomial ground2("ground2","background",x,RooArgList(co2o1,co2o2,co2o3,co2o4));
    
   RooRealVar eff("eff","efficiency",0.9,0,1);
   RooRealVar signal2("signal2"," ",1000,0,10000000);//event number
   RooFormulaVar signal1("signal1","@0*@1",RooArgList(signal2,eff));
   //RooRealVar signal2("signal2"," ",1000,0,10000000);//event number
   RooRealVar background1("background1"," ",600,0,100000000);
   RooRealVar background2("background2"," ",600,0,100000000);

   RooAddPdf model1("model1","model1",RooArgList(cbshape,ground1),RooArgList(signal1,background1));
   RooAddPdf model2("model2","model2",RooArgList(cbshape,ground2),RooArgList(signal2,background2));

   //RooAddPdf *sum;
   RooDataSet *dataset1 = new RooDataSet("data1","data1",RooArgSet(x),Import(*dataraw1));
   RooDataSet *dataset2 = new RooDataSet("data2","data2",RooArgSet(x),Import(*dataraw2));
   
   RooCategory cate("cate","cate");
   cate.defineType("ARec");
   cate.defineType("AEvt");
   
   // combined data
   RooDataSet combData("combData","combined data",x,Index(cate),Import("ARec",*dataset1),Import("AEvt",*dataset2));
   
   // construct a simultaneous pdf using category sample as index
   RooSimultaneous simPdf("simPdf","simultaneous pdf",cate);
   simPdf.addPdf(model1,"ARec");
   simPdf.addPdf(model2,"AEvt");

   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   
   simPdf.fitTo(combData);
   
   RooPlot* frame1 = x.frame(Title("All reconstructed"));
   RooPlot* frame2 = x.frame(Title("All events"));
   combData.plotOn(frame1,Cut("cate==cate::ARec"));
   combData.plotOn(frame2,Cut("cate==cate::AEvt"));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),Components("ground1"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kBlue));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),Components("cbshape"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kGreen));
   simPdf.plotOn(frame1,Slice(cate,"ARec"),ProjWData(cate,combData));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),Components("ground2"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kBlue));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),Components("cbshape"),ProjWData(cate,combData),LineStyle(kDashed),LineColor(kGreen));
   simPdf.plotOn(frame2,Slice(cate,"AEvt"),ProjWData(cate,combData));
   
   TCanvas *c1 = new TCanvas("c1","tracking efficiency",1000,400);
   c1->Divide(2);
   c1->cd(1);gPad->SetLeftMargin(0.15);frame1->GetYaxis()->SetTitleOffset(1.4);frame1->Draw();
   c1->cd(2);gPad->SetLeftMargin(0.15);frame2->GetYaxis()->SetTitleOffset(1.4);frame2->Draw();
   
   double sig1err = sqrt(pow(signal2.getError()*eff.getVal(),2)+pow(signal2.getVal()*eff.getError(),2));
   char name[100];
   Npar = 10;
   TPaveText *pt1 = new TPaveText(0.6,0.68,0.88,0.88,"BRNDC");
   pt1->SetBorderSize(0); pt1->SetFillStyle(4000);
   pt1->SetTextAlign(12); pt1->SetTextFont(42); pt1->SetTextSize(0.035);
   sprintf(name,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError()); pt1->AddText(name);
   sprintf(name,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError()); pt1->AddText(name);
   sprintf(name,"signal1 = %.1f #pm %.1f",signal1.getVal(),sig1err); pt1->AddText(name);
   sprintf(name,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,frame1->chiSquare(Npar)); pt1->AddText(name);
   c1->cd(1); pt1->Draw();
 
   TPaveText *pt2 = new TPaveText(0.6,0.68,0.88,0.88,"BRNDC");
   pt2->SetBorderSize(0); pt2->SetFillStyle(4000);
   pt2->SetTextAlign(12); pt2->SetTextFont(42); pt2->SetTextSize(0.035);
   sprintf(name,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError()); pt2->AddText(name);
   sprintf(name,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError()); pt2->AddText(name);
   sprintf(name,"signal2 = %.1f #pm %.1f",signal2.getVal(),signal2.getError()); pt2->AddText(name);
   sprintf(name,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,frame2->chiSquare(Npar)); pt2->AddText(name);
   c1->cd(2); pt2->Draw();
   sprintf(name,"p_spectrum_%s",namesfx);
   c1->SetName(name);
   c1->Write();
   sprintf(name,"p_spectrum_%s.png",namesfx);
   c1->Print(name);

  // c1->Print("test.pdf");
   delete pt1;
   delete pt2;
   delete c1;
   
   
   ofstream fitpar("fitpar",std::ios::app);
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<"\t sigNo = " << sig1err << "\t sigNoE = " << sig1err;
   fitpar<<"\t bckNo = " << background1.getVal() << "\t bckNoE = " << background1.getError();
   fitpar<<"\t chi1: "<<frame1->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   fitpar<<"\t sigNo = " << signal2.getVal() << "\t sigNoE = " << signal2.getError();
   fitpar<<"\t bckNo = " << background2.getVal() << "\t bckNoE = " << background2.getError();
   fitpar<<"\t chi2: "<<frame2->chiSquare(Npar) <<  std::endl;;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   if (nsig!=0) {nsig[0] = signal1.getVal(); nsig[1] = signal2.getVal(); }
   if (esig!=0) {esig[0] = sig1err; esig[1] = signal2.getError(); }
   delete dataset1;
   delete dataset2;
   return signal1.getVal();
}


