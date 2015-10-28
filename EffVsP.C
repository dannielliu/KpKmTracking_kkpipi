#include <fstream>
#include <iostream>
#include <sstream>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "CommonFunc.h"
using namespace std;

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
  g_mcKm->SetLineColor(2);
  g_mcKm->SetMarkerColor(2);
  g_mcKm->SetFillColor(0);
  g_mcKm->Draw("AP");
  
  g_dataKm->SetLineColor(1);
  g_dataKm->SetMarkerColor(1);
  g_dataKm->SetFillColor(0);
  g_dataKm->Draw("P");
  
  TLegend *legend = new TLegend(0.5,0.15,0.8,0.35);
  legend->AddEntry(g_mcKm,"MC K^{+}");
  legend->AddEntry(g_dataKm,"DATA K^{+}");
  legend->Draw();

  c1->Print("EffVsP_MC_Data.pdf");
  return 0;

}

