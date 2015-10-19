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

class TrackingAlg
{
public:
   double mass;
   int nran;
   vector<double> pnode; // nran+1
   vector<double> npran; // nran
   vector<double> epran;
   vector<double> neff;
   vector<double> eeff;
   vector<int> ValidRange;

   TH1D *hMmiss34 ;// = new TH1D("hMmiss34" ,"M(K)",100,0.2,0.8);
   TH1D *hMmiss4  ;// = new TH1D("hMmiss4" ,"M(K)",100,0.2,0.8);
   vector<TTree*> data3trk; // nran
   vector<TTree*> data4trk;
   vector<TTree*> dataTtrk;
   TTree* data3trkAll;
   TTree* data4trkAll;
   TTree* dataTtrkAll;

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
       delete data3trk.at(i);
       delete data4trk.at(i);
       delete dataTtrk.at(i);
     }
   }
   
   void init()
   {
     cout<<"Construct tracking alg"<<endl;
     // n pt range
     nran=20;
     cout<<nran<<endl;
     for (int i=0;i<nran+1;i++){
       pnode.push_back(0+i*0.1);
     }

     hMmiss34  = new TH1D("hMmiss34" ,"M(K)",100,0.2,0.8);
     hMmiss4  = new TH1D("hMmiss4" ,"M(K)",100,0.2,0.8);
     TTree* tree;
     char namei[200];
     for (int i=0; i<nran; i++){
       cout<<i<<"th tree constructing"<<endl;
       sprintf(namei,"data3trk_%d",i);
       tree = new TTree(namei,namei);
       tree->Branch("mass", &mass, "mass/D");
       data3trk.push_back(tree);
       sprintf(namei,"data4trk_%d",i);
       tree = new TTree(namei,namei);
       tree->Branch("mass", &mass, "mass/D");
       data4trk.push_back(tree);
       sprintf(namei,"dataTtrk_%d",i); // both 3 and 4 tracks
       tree = new TTree(namei,namei);
       tree->Branch("mass", &mass, "mass/D");
       dataTtrk.push_back(tree);
     }
     sprintf(namei,"data3trk_All");
     data3trkAll = new TTree(namei,namei);
     data3trkAll->Branch("mass", &mass, "mass/D");
     sprintf(namei,"data4trk_All");
     data4trkAll = new TTree(namei,namei);
     data4trkAll->Branch("mass", &mass, "mass/D");
     sprintf(namei,"dataTtrk"); // both 3 and 4 tracks
     dataTtrkAll = new TTree(namei,namei);
     dataTtrkAll->Branch("mass", &mass, "mass/D");

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
       cout<<"\n\n\n\n#################\n"<< i << endl;
       cout<<"size of 4 trk is "<< data4trk.at(i)->GetEntries()<<endl;
       cout<<"size of 3/4 trk is "<< dataTtrk.at(i)->GetEntries()<<endl;
       sprintf(suffix,"Trk4_%.1f", (pnode[i+1]+pnode[i])/2);
       nsig = FitSpectrum(data4trk.at(i), i, suffix, &n4trk, &e4trk);
       sprintf(suffix,"All_%.1f", (pnode[i+1]+pnode[i])/2);
       nsig = FitSpectrum(dataTtrk.at(i), i, suffix, &nTtrk, &eTtrk);
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
 
       npran.push_back( (pnode[i+1]+pnode[i])/2);
       epran.push_back( (pnode[i+1]-pnode[i])/2);
       neff.push_back( TrkEff);
       eeff.push_back( EffErr);
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
       data4trk.at(i)->Draw(drawopt);
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
     TGraphErrors *graph = new TGraphErrors(nran,&npran.at(0),&neff.at(0),&epran.at(0),&eeff.at(0));
     graph->GetXaxis()->SetTitle("p_{t} (GeV/c)");
     graph->GetYaxis()->SetTitle("eff");
     graph->Draw();
     c1->Write("Eff_p");
   }

};

double GetEnergy(int run)
{
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
  return -1;
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

   TLorentzVector trk[4];
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
   TTree *data3trk[nran];
   TTree *data4trk[nran];
   TTree *dataTtrk[nran];
   char namei[100];
   for (int i=0; i<nran; i++){
     data3trk[i] = trkalg->data3trk.at(i);
     data4trk[i] = trkalg->data4trk.at(i);
     dataTtrk[i] = trkalg->dataTtrk.at(i);
   }
   TTree* data3trkAll = trkalg->data3trkAll; //= new TTree(namei,namei);
   TTree* data4trkAll = trkalg->data4trkAll; //= new TTree(namei,namei);
   TTree* dataTtrkAll = trkalg->dataTtrkAll; //= new TTree(namei,namei);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%100000 == 0) cout<<"current entry is "<< jentry << std::endl;
      int pichrg = pi1chrg+pi2chrg+ka1chrg;
      if (pichrg!=1 && pichrg!=-1) { /*cout<<"id: "<< jentry<<endl;*/ continue;}

      trk[0].SetVectMag(TVector3(pi1px,pi1py,pi1pz),mpi);
      trk[1].SetVectMag(TVector3(pi2px,pi2py,pi2pz),mpi);
      trk[2].SetVectMag(TVector3(ka1px,ka1py,ka1pz),mka);
      if (ka_chrg != 0) trk[3].SetVectMag(TVector3(ka_px,ka_py,ka_pz),mka);
      else trk[3].SetVectMag((tot4p - (trk[0]+trk[1]+trk[2])).Vect(),mka);
      TLorentzVector trktot = trk[0]+trk[1]+trk[2]+trk[3];

      double pkm = trk[3].Rho();
      // transverse momentum
      double ptkm = trk[3].Pt();
      double Coskm = trk[3].CosTheta();
      int idx = (int)(ptkm/0.1);
      if (idx>=nran) continue;
      if (fabs(Coskm)>0.93) continue;
      mass = trktot.M();

      if (ka_chrg == 0) {
	hm1->Fill(mass); hm3->Fill(mass);
      }
      else if (pichrg+ka_chrg == 0) {
	hm2->Fill(mass); 
	hm3->Fill(mass);
      
        if (fabs(trktot.M()-Ecm)<0.02) {
          hth->Fill(trk[3].CosTheta());
          hp->Fill(trk[3].Rho());  
        }
      }

      TLorentzVector pmiss = tot4p - (trk[0]+trk[1]+trk[2]);
      mass = pmiss.M();
      hMmiss34->Fill(mass);
      data3trk[idx]->Fill();
      dataTtrk[idx]->Fill();
      data3trkAll->Fill();
      dataTtrkAll->Fill();
      if (ka_chrg!=0) 
      {
	data4trk[idx]->Fill();
	data4trkAll->Fill();
        hMmiss4->Fill(pmiss.M());
	hMrec4->Fill(trk[3].M());
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
   
   hm1->Write();
   hm2->Write();
   hm3->Write();
   hm4->Write();
   hth->Write();
   hp->Write();
   
   hMmiss34->Write();
   hMmiss4->Write();
   //hPmiss34->Write();
   //hPmiss4->Write();

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
  if (argc==1) t = new KaonTrack();
  else {
    for (int i=1;i<argc;i++){
      cout<<"user input file "<< argv[i]<<endl;
      TFile *f = new TFile(argv[i]);
      TTree* tree = (TTree*)f->Get("KaonTrack");
      t = new KaonTrack(tree);
      cout<<"Loop?"<<endl;
      ofile->cd();
      t->Loop(&trking);
      delete t;
    }
  }
  trking.ShowInHist();
  trking.FitDataSet();
  trking.CreateEffVsPt();
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
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
using namespace RooFit;

 // 3.08 data
/* 
double FitSpectrum(TTree *dataraw, double beame, const char* namesfx, double *nsig, double *esig)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double peakvalue = mka;
   double beamlow=0.3;
   double beamup=0.8;
   // try to use roofit
   RooRealVar x("mass","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue-0.001,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma("sigma","width of gaussian",0.0135+0.01/5*beame,0.013,0.03);
 //RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooRealVar mean2("mean2","mean of gaussian",beame,peakvalue,beame+0.1);
 //RooRealVar mean2("mean2","mean of gaussian",peakvalue+0.2,peakvalue+0.05,3.0);
 //RooRealVar sigma2("sigma2","width of gaussian",0.09,0.03,0.15);
 //RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   
   RooRealVar co1("co1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2("co2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co3("co3","coefficient #2",   0 ,-100.,100.);
   RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   
   RooRealVar signal("signal"," ",1000,0,10000000);//event number
   RooRealVar background("background"," ",600,0,100000000);
     
   //RooRealVar sigma2("sigma2","width of gaussian",0.01,0.008,0.02);
   RooRealVar alpha1("alpha1","#alpha",-1.2,-5.0,5.0);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
   
   // fit bhabha sigma and mean
 //RooRealVar meanb("meanb","mean of gaussian",beame,beame-0.015,beame+0.005);
 //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.004,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   
 //RooRealVar alphab("alphab","#alpha",1.0,5);
 //RooRealVar nnnb("nnnb","n",100,1,200);
 //RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
 
 //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
 //RooRealVar sigmae("sigmae","width of gaussian",0.005,0.004,0.02);
 //RooGaussian gause("gause","gauss(x,m,s)",x,meane,sigmae);
   
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_Mmiss_%s",namesfx);
   xframe = x.frame(Title("fit p"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape,ground),RooArgList(signal,background));
   Npar = 9;
   //sum->fitTo(*dataset,Range(peakvalue-0.07,peakvalue+0.07));
   //x.setRange("sigragi",peakvalue-0.18,peakvalue+0.14);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.015);
   sum->fitTo(*dataset);
   //sum->fitTo(*dataset,Range("sigragi"));
   //sum->fitTo(*dataset,Range(peakvalue-0.1,beame+0.01));
   //sum->fitTo(*dataset,"e",Range(peakvalue-0.1,beame));
   //gause.fitTo(*dataset,Range(beame-0.02, beame+0.01));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(cbshape),LineStyle(2),LineColor(2) );
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3)  );
   sum->plotOn(xframe  );
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.15,0.65,0.45,0.90,"BRNDC");
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
 //sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
 //pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.png",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
   fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
   fitpar<< std::endl;;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
 //double intel3 = mean.getVal()-3*sigma.getVal();
 //double inteu3 = mean.getVal()+3*sigma.getVal();
 //double intel5 = mean.getVal()-5*sigma.getVal();
 //double inteu5 = mean.getVal()+5*sigma.getVal();
 ////RooRealVar xint("xint","xint",peakvalue-0.1,beame);
 ////x.setRange("sigragi",peakvalue-0.1,beame);
 ////x.setRange("sigragi",peakvalue-0.1,beame+0.01);
 //x.setRange("sigrag3",intel3, inteu3);
 //x.setRange("sigrag5",intel5, inteu5);
 //RooAbsReal* intsigi = cbshape.createIntegral( x,  NormSet(x),  Range("sigragi"));
 //RooAbsReal* intbcki = ground.createIntegral(x,  NormSet(x),  Range("sigragi"));
 //RooAbsReal* intsig3 = cbshape.createIntegral( x,  NormSet(x),  Range("sigrag3"));
 //RooAbsReal* intbck3 = ground.createIntegral(x,  NormSet(x),  Range("sigrag3"));
 //RooAbsReal* intsig5 = cbshape.createIntegral( x,  NormSet(x),  Range("sigrag5"));
 //RooAbsReal* intbck5 = ground.createIntegral(x,  NormSet(x),  Range("sigrag5"));
 //double signalN3 =   signal.getVal()*intsig3->getVal();
 //double backN3 = background.getVal()*intbck3->getVal();
 //double signalN5 =   signal.getVal()*intsig5->getVal();
 //double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
 //ofstream angSNR("fitpar.dat",std::ios::app);
 //angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
 //std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
 //std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
 //std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   if (nsig!=0) *nsig = signal.getVal();
   if (esig!=0) *esig = signal.getError();
   return signal.getVal();
}
*/

// for 4260 data
/*
double FitSpectrum(TTree *dataraw, double beame, const char* namesfx, double *nsig, double *esig)
{
   cout<<dataraw->GetName()<<"\t"<< dataraw->GetBranch("mass")->GetEntries()<< endl;
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double peakvalue = mka;
   double beamlow=0.2;
   double beamup=0.9;
   // try to use roofit
   RooRealVar x("mass","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue-0.001,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma("sigma","width of gaussian",0.015+0.015/5*beame,0.013,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooRealVar mean2("mean2","mean of gaussian",beame,peakvalue,beame+0.1);
 //RooRealVar mean2("mean2","mean of gaussian",peakvalue+0.2,peakvalue+0.05,3.0);
 //RooRealVar sigma2("sigma2","width of gaussian",0.09,0.03,0.15);
 //RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   
   RooRealVar co1("co1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2("co2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co3("co3","coefficient #2",   0 ,-100.,100.);
   RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   
   RooRealVar signal("signal"," ",1000,0,10000000);//event number
   RooRealVar background("background"," ",600,0,100000000);
     
   //RooRealVar sigma2("sigma2","width of gaussian",0.01,0.008,0.02);
 //RooRealVar alpha1("alpha1","#alpha",0.5,-5.0,5.0);
 //RooRealVar nnn1("n1","n",100,1,200);
 //RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   //RooDataHist *datahist;
   //RooPlot *xframe;
   //RooDataHist *data_6pi;
   
   // fit bhabha sigma and mean
 //RooRealVar meanb("meanb","mean of gaussian",beame,beame-0.015,beame+0.005);
 //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.004,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   
 //RooRealVar alphab("alphab","#alpha",1.0,5);
 //RooRealVar nnnb("nnnb","n",100,1,200);
 //RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
 
 //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
 //RooRealVar sigmae("sigmae","width of gaussian",0.005,0.004,0.02);
 //RooGaussian gause("gause","gauss(x,m,s)",x,meane,sigmae);
   
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_Mmiss_%s",namesfx);
   RooPlot* xframe = x.frame(Title("fit p"));
   //dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   RooDataSet *dataset = new RooDataSet(tmpchr,"data",dataraw,RooArgSet(x));


 //char name[100];
 //char drawopt[100];
 //TH1D *hmKm;
 //sprintf(name,"mass_%f",beame);
 //hmKm = new TH1D(name,name,100,0.1,0.9);
 //sprintf(drawopt,"mass>>%s",name);
 //dataraw->Draw(drawopt);
 //datahist = new RooDataHist(tmpchr,"data",x,hmKm);
    
    
    
    //sum = new RooAddPdf("sum","sum",RooArgList(cbshape,ground),RooArgList(signal,background));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background));
   Npar = 7;
   //sum->fitTo(*datahist);
 //sum->fitTo(*dataset);
   //sum->fitTo(*dataset,Range("sigragi"));
   //datahist->plotOn(xframe);
   dataset->plotOn(xframe);
   //sum->plotOn(xframe,Components(cbshape),LineStyle(2),LineColor(2)  );
 //sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2)  );
 //sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3)  );
 //sum->plotOn(xframe  );
 ////sum->plotOn(xframe);
   xframe->Draw();
 //TPaveText *pt = new TPaveText(0.15,0.65,0.45,0.90,"BRNDC");
 //pt->SetBorderSize(0);
 //pt->SetFillStyle(4000);
 //pt->SetTextAlign(12);
 //pt->SetTextFont(42);
 //pt->SetTextSize(0.035);
 //sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
 //pt->AddText(tmpchr);
 //pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.png",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
   fitpar<<" ene = "<< beame <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
   fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
   fitpar<< std::endl;;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   //delete data_6pi;
   delete xframe;
   delete dataset;
   //delete datahist;
   delete sum;
   if (nsig!=0) *nsig = signal.getVal();
   if (esig!=0) *esig = signal.getError();
   return signal.getVal();
}
*/


 // 3.08 data
 
double FitSpectrum(TTree *dataraw, double beame, const char* namesfx, double *nsig, double *esig)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double peakvalue = mka;
   double beamlow=0.3;
   double beamup=0.8;
   // try to use roofit
   RooRealVar x("mass","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue-0.001,peakvalue-0.01,peakvalue+0.01);
   RooRealVar sigma("sigma","width of gaussian",0.0135+0.01/10*beame,0.013,0.03);
 //RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   
   RooRealVar co1("co1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2("co2","coefficient #2",   0 ,-100.,100.);
   RooRealVar co3("co3","coefficient #2",   0.1 ,-100.,100.);
   RooRealVar co4("co4","coefficient #2",   0 ,-100.,100.);
   //RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   RooPolynomial ground("ground","background",x,RooArgList(co1,co2,co3,co4));
   
   RooRealVar signal("signal"," ",1000,0,10000000);//event number
   RooRealVar background("background"," ",600,0,100000000);
     
   //RooRealVar sigma2("sigma2","width of gaussian",0.01,0.008,0.02);
   RooRealVar alpha1("alpha1","#alpha",-1.2,-5.0,5.0);
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
  // dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   
   char name[100];
   char drawopt[100];
   TH1D *hmKm;
   sprintf(name,"mass_%f",beame);
   hmKm = new TH1D(name,name,100,0.1,0.9);
   sprintf(drawopt,"mass>>%s",name);
   dataraw->Draw(drawopt);
   RooDataHist *datahist = new RooDataHist(tmpchr,"data",x,hmKm);
 
   
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape,ground),RooArgList(signal,background));
   Npar = 9;
   sum->fitTo(*datahist);
   //dataset->plotOn(xframe);
   datahist->plotOn(xframe);
   sum->plotOn(xframe,Components(cbshape),LineStyle(2),LineColor(2) );
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3)  );
   sum->plotOn(xframe  );
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.15,0.65,0.45,0.90,"BRNDC");
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
   fitpar<< std::endl;;
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



