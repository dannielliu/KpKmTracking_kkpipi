#include <fstream>
#include <iostream>
#include <sstream>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
using namespace std;

int main()
{
  ifstream txtfile("EffVsP.txt");
  if (!txtfile.is_open()) return -1;
  char line[100];
  txtfile.getline(line,100);
  cout<<"first line of EffVsP.txt:\n\t"<<line<<endl;

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
      cout << "i "<< i << " npe " << npe[i] << " aaaa " <<endl; 
      continue;
    }
    if (i==0 || i==npran-2) {
      npe[i] = (np[i+1]-np[i])/2.0;
      cout << "i "<< i << " np " << npe[i] << " bbbb "<<endl; 
      continue;
    }

    if (((np[i+1]-np[i])-(np[i]-np[i-1])>1e-6)) 
    {
        if (i!=0 && i!=npran-1) {
	  if (((np[i]-np[i-1])-(np[i-1]-np[i-2]))<1e-6)  npe[i] = (np[i]-np[i-1])/2.0;
	  if (((np[i+2]-np[i+1])-(np[i+1]-np[i]))<1e-6)  npe[i] = (np[i+2]-np[i+1])/2.0;
          cout << "i "<< i << " np " << npe[i] << " cccc"<<endl; 
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
  graph->SetTitle("eff Vs. p_{t}");
  graph->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  graph->GetYaxis()->SetTitle("#epsilon_{eff}");
  graph->Draw("AP");
  c1->Print("EffVsP.pdf");
  return 0;
}
