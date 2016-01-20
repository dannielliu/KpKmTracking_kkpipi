#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TF1.h"
#include "TFile.h"
#include "CommonFunc.h"

int EffVsP(const char* filename, TGraphErrors *&graph_mckm)
{
  //ifstream txtfile("EffVsP.txt_mcKm");
  ifstream txtfile(filename);
  if (!txtfile.is_open()) return -1;
  char line[100];
  txtfile.getline(line,100);
  cout<<"first line of "<<filename << "\n\t"<<line<<endl;

  int npran=0;
  double np[100];
  double npe[100];
  double neff[100];
  double neffe[100];
//double ndeff[100];
//double ndeffe[100];
  istringstream iss;
  while (!txtfile.eof()){
    double pave, N4trk, N4err, pave2,NAtrk, NAerr;
    
    txtfile.getline(line,100);
    if (line[0]=='\0' || line[0]=='#') continue;
    iss.clear();
    iss.str(line);
    iss>>pave >> N4trk >> N4err;
    
    txtfile.getline(line,100);
    if (line[0]=='\0' || line[0]=='#') continue;
    iss.clear();
    iss.str(line);
    iss>>pave2 >> NAtrk >> NAerr;

    if (pave!=pave2) {
      cout << "Error: p not equal!!!"<<endl;
      continue;
    }
    cout << pave << "\t" << N4trk << "\t" << N4err<< "\t" << NAtrk << "\t"<< NAerr<< endl;

    np[npran] = pave;
    neff[npran] = N4trk/NAtrk;
    neffe[npran] = 1./NAtrk*sqrt((1-2*neff[npran])*pow(N4err,2)+pow(neff[npran],2)*pow(NAerr,2));

    npran++;
  }
  //extract uncertainty of pt
  for (int i=0;i<npran;i++){
    if (i==npran-1 || i==1) {
      npe[i] = (np[i]-np[i-1])/2.0;
      //cout << "i "<< i << " npe " << npe[i] << " aaaa " <<endl; 
      continue;
    }
    if (i==0 || i==npran-2) {
      npe[i] = (np[i+1]-np[i])/2.0;
      //cout << "i "<< i << " np " << npe[i] << " bbbb "<<endl; 
      continue;
    }

    if (((np[i+1]-np[i])-(np[i]-np[i-1])>1e-6)) 
    {
        if (i!=0 && i!=npran-1) {
	  if (((np[i]-np[i-1])-(np[i-1]-np[i-2]))<1e-6)  npe[i] = (np[i]-np[i-1])/2.0;
	  if (((np[i+2]-np[i+1])-(np[i+1]-np[i]))<1e-6)  npe[i] = (np[i+2]-np[i+1])/2.0;
          //cout << "i "<< i << " np " << npe[i] << " cccc"<<endl; 
	  continue;
	}
    }
    npe[i] = (np[i+1]-np[i])/2.0;
  }
  for (int i=0; i<npran;i++){
    cout << "p error "<< npe[i]<<endl;
  } 
  
  
  TCanvas *c1 = new TCanvas();
  TGraphErrors* graph = new TGraphErrors(npran,np,neff,npe,neffe);
  graph->SetTitle("eff Vs p_{t}");
  graph->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  graph->GetYaxis()->SetTitle("#epsilon_{eff}");
  graph->Draw("AP");

  char output_pdf_name[100];
  sprintf(output_pdf_name,"EffVsP_%s.pdf",getPureName(filename));
  c1->Print(output_pdf_name);
  
  graph_mckm = graph;
  delete c1;
  return 0;
}

void help()
{
  cout<<"Usage:"<<endl;
  cout<<"./EffVsP <args>"<<endl;
}

int main(int argc, char** argv)
{

  if (argc ==2 && strncmp(argv[1],"--help",6)==0) {help(); return 0;}
  
  TGraphErrors* g_EVP;
  if (argc==1) return EffVsP("EffVsP.txt",g_EVP);
  if (argc==2 && argv[1][0]!='-') return EffVsP(argv[1],g_EVP);
  
  char txt_mcKm[] = "3080/EffVsP.txt_mcKm";
  char txt_dataKm[] = "3080/EffVsP.txt_dataKm";
  char* data_mcKm = txt_mcKm;
  char* data_dataKm = txt_dataKm;
  
  if (argc>2) {
   // for (int iarg=1;iarg<argc;iarg++) {
      cout<<"test"<<endl;
      if (argv[1][0]!='-') data_mcKm = argv[1];
      if (argv[2][0]!='-') data_dataKm = argv[2];
    //}
  }

  TGraphErrors* g_mcKm = 0;
  TGraphErrors* g_dataKm = 0;
  int status;
  status = EffVsP(data_mcKm,g_mcKm);
  status = EffVsP(data_dataKm,g_dataKm);
  
  cout<<g_dataKm <<"\t"<<g_mcKm << endl;
  if (g_mcKm==0 || g_dataKm==0) return 1;

  TCanvas *c1 = new TCanvas();
  //c1->SetMargin(0.15,0.1,0.15,0.1);
  //c1->Divide();
  TPad *pad1 = new TPad("pad1","pad1",0.05,0.3,1,0.98);
  TPad *pad2 = new TPad("pad2","pad2",0.05,0.02,1,0.3);
  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0.02);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.3);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  g_mcKm->SetTitle("");
  g_mcKm->GetXaxis()->SetTitleSize(0.05);
  g_mcKm->GetXaxis()->SetLabelSize(0.05);
  g_mcKm->GetXaxis()->SetLabelOffset(1.1);
  g_mcKm->GetXaxis()->SetTitleOffset(1.1);
  g_mcKm->GetXaxis()->SetNdivisions(505);
  g_mcKm->GetYaxis()->SetTitleSize(0.05);
  g_mcKm->GetYaxis()->SetLabelSize(0.05);
  g_mcKm->GetYaxis()->SetTitleOffset(1.1);
  g_mcKm->GetYaxis()->SetNdivisions(505);
  g_mcKm->GetYaxis()->SetRangeUser(0,1.2);
  g_mcKm->SetLineColor(2);
  g_mcKm->SetMarkerColor(2);
  g_mcKm->SetFillColor(0);
  g_mcKm->Draw("AP");
  
  g_dataKm->SetLineColor(1);
  g_dataKm->SetMarkerColor(1);
  g_dataKm->SetFillColor(0);
  g_dataKm->Draw("P");
  
  TLegend *legend = new TLegend(0.5,0.4,0.7,0.55);
  legend->AddEntry(g_mcKm,"MC K^{-}");
  legend->AddEntry(g_dataKm,"DATA K^{-}");
  legend->Draw();

  //c1->Print("EffVsP_MC_Data.pdf");
  
  pad2->cd();
  // compare two graphs
  double* xx;
  double* xe;
  double yy[100];
  double ye[100];
  double yyres[100];
  double yeres[100];
  int np_mc = g_mcKm->GetN();
  int np_data = g_dataKm->GetN();
  if (np_mc != np_data) {
    cout<< "size of two graph are not equal, pleare check it."<<endl;
    cout<< "mc size is "<< np_mc<<", data size is "<< np_data<<endl;
    return -1;
  }
  xx = g_mcKm->GetX();
  xe = g_mcKm->GetEX();
  for (int i=0 ; i<np_mc; i++){
    double y1,y2,ye1,ye2;
    y1 = g_dataKm->GetY()[i]; y2 = g_mcKm->GetY()[i];
    yy[i] = y1 - y2;
    ye[i] = sqrt(pow(g_dataKm->GetEY()[i],2)+pow(g_mcKm->GetEY()[i],2));
    if (fabs(y1)<1e-9) {yyres[i]=0;yeres[i]=0;continue;}
    yyres[i] = (y2 - y1)/y1*100;
    yeres[i] = ye[i]/y1*100;
  }
  TGraphErrors* effcmp = new TGraphErrors(np_mc,xx,yy,xe,ye);
  effcmp->SetTitle("");
  effcmp->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  effcmp->GetYaxis()->SetTitle("#epsilon_{data} - #epsilon_{MC}");
  
  effcmp->GetXaxis()->SetTitleSize(0.12);
  effcmp->GetXaxis()->SetLabelSize(0.12);
  effcmp->GetXaxis()->SetTitleOffset(1.1);
  effcmp->GetXaxis()->SetNdivisions(505);
  effcmp->GetYaxis()->SetRangeUser(-0.1,0.1);
  effcmp->GetYaxis()->SetTitleSize(0.12);
  effcmp->GetYaxis()->SetLabelSize(0.12);
  effcmp->GetYaxis()->SetTitleOffset(0.4);
  effcmp->GetYaxis()->SetNdivisions(502);
 
  //c1->Clear();
  effcmp->Draw("AP");
  TF1 f1("f1","0",0,1.5);
  f1.SetLineStyle(kDashed);
  f1.Draw("same");
  c1->Update();
  c1->Print("effcmp.pdf");
  
  // relatively
  TGraphErrors* effcmp2 = new TGraphErrors(np_mc,xx,yyres,xe,yeres);
  effcmp2->SetTitle("compare data and MC");
  effcmp2->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  effcmp2->GetYaxis()->SetTitle("#Delta#epsilon (%)");
  c1->Clear();
  effcmp2->Draw("AP");
//TF1 f1("f1","0",0,1.5);
//f1.SetLineStyle(kDashed);
  f1.Draw("same");
  c1->Print("effcmp_relative.pdf");
  
  TFile *f = new TFile("output.root","recreate");
  f->WriteTObject(effcmp,"diff");
  f->WriteTObject(effcmp2,"rela_diff");

  return 0;

}

