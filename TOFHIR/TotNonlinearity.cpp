// per compilare: g++ -Wall -o TotNonlinearity `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit TotNonlinearity.cpp

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/MinimizerOptions.h"

#include <sys/stat.h>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>


#include "setTDRStyle.C"

using namespace std;

float PDE(float ov){
  //float pde = 0.4198 * pow(ov,3) - 4.8388 * pow(ov, 2) + 21.277 * ov + 3.6737; //  Yuri
  float pde = 40.1017 * (1.0 - exp(-ov/1.39665)); // mia parametrizzazione basata sui punti di Yuri
  pde/=100;
  return (pde);
}


float Gain(float ov){
  float gain = 97.869 * ov + 36.619;
  gain*=1000;
  return(gain); 
}

float ENF(float ov){
  float enf = 0.00258 * pow(ov,2) -0.00219 * ov + 1.00644;
  return(enf); 
}

//=========== MAIN =================================

int main(int argc, char** argv)
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000); 
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);
 
  setTDRStyle();

  string outdir = "NonLinearityToT/";
  
  int channel = 0;
  
  TFile *file = TFile::Open(argv[1]);
  cout << "Analyzing " << argv[1] <<endl;
  
  TTree* tree = (TTree*)file->Get("data");

  float ov; 
  float step2; 
  float tot;
  unsigned int channelID;

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("step1",1);
  tree->SetBranchStatus("step2",1);
  tree->SetBranchStatus("tot",1);
  tree->SetBranchStatus("channelID",1);
   
  tree->SetBranchAddress("step1",&ov);
  tree->SetBranchAddress("step2",&step2);
  tree->SetBranchAddress("tot",&tot);
  tree->SetBranchAddress("channelID",&channelID);

  cout << "Number of entries : " << tree->GetEntries() <<endl;

  tree -> GetEntry(0);

  int vth1 = int(step2/100) - 1 ;
  int vth2 = int(step2)%100 - 1 ;

  cout << "step2 = " << step2 << "   vth1 = " << vth1 << "  vth2 = "<< vth2 << endl;

  // -- book histograms
  vector<float> Vovs;
  Vovs.clear();
  Vovs.push_back(1.0);
  Vovs.push_back(2.0);
  Vovs.push_back(3.0);
  Vovs.push_back(4.0);
  Vovs.push_back(5.0);
  Vovs.push_back(6.0);

  map<float, TH1F*> h_tot;
  map<float, TH1F*> h_tot_corr;
  for ( int i = 0; i < Vovs.size(); i++){
    h_tot[Vovs[i]] = new TH1F(Form("h_tot_ov%.02f",Vovs[i]), Form("h_tot_ov%.02f",Vovs[i]), 2000, 0, 1000);
    h_tot_corr[Vovs[i]] = new TH1F(Form("h_tot_corr_ov%.02f",Vovs[i]), Form("h_tot_corr_ov%.02f",Vovs[i]), 4000, 0, 1000);
  }


 
  // -- loop over events
  for (int i = 0; i < tree->GetEntries(); i++){
    
    tree -> GetEntry(i);
    
    if (channelID != channel) continue;
    
    for (int iov = 0; iov < Vovs.size(); iov++){
      if (ov == Vovs[iov])
	h_tot[ov] -> Fill(tot/1000.);
    }
    
  }// -- end loop over events




  TGraph *g_pde_vs_ov = new TGraph();
  g_pde_vs_ov->SetName("g_pde_vs_ov");
  g_pde_vs_ov->SetMarkerStyle(20);

  TGraph *g_gain_vs_ov = new TGraph();
  g_gain_vs_ov->SetName("g_gain_vs_ov");
  g_gain_vs_ov->SetMarkerStyle(20);
  
   
  TGraph *g_tot_vs_amp = new TGraph();
  g_tot_vs_amp->SetName("g_tot_vs_amp");
  g_tot_vs_amp->SetMarkerStyle(20);
    
  //TGraph *g_amp_vs_tot = new TGraph();
  TGraphErrors *g_amp_vs_tot = new TGraphErrors();
  g_amp_vs_tot->SetName("g_amp_vs_tot");
  g_amp_vs_tot->SetMarkerStyle(20);


 
  // -- fit 511keV photopeak
  TF1 *fGaus[10];
  for (int iov = 0; iov < Vovs.size(); iov++){
    float ov = Vovs[iov];
    fGaus[iov] = new TF1(Form("fGaus_%.02f",ov),"gaus", 0, 500);
    //fGaus[iov]->SetLineColor(800+20*iov);
    fGaus[iov]->SetLineColor(1);
    h_tot[ov]->GetXaxis()->SetRangeUser(30., 10000.);
    int maxbin = h_tot[ov] -> GetMaximumBin();
    float xmax = h_tot[ov] -> GetBinCenter(maxbin);
    fGaus[iov]->SetParameter(1,xmax);
    fGaus[iov]->SetParameter(2,xmax/10.);
    fGaus[iov]->SetRange(xmax-5, xmax+5);
    h_tot[ov]->GetXaxis()->SetRangeUser(0., 10000.);
    h_tot[ov] ->Fit( fGaus[iov], "QRS");
    fGaus[iov]->SetRange(xmax-0.5*fGaus[iov]->GetParameter(2), xmax+0.5*fGaus[iov]->GetParameter(2));
    h_tot[ov] ->Fit( fGaus[iov], "QRS");
    double peakTot    = fGaus[iov]->GetParameter(1);
    double peakTotErr = fGaus[iov]->GetParError(1);
    double amp = PDE(ov) * Gain(ov) / 1000. * ENF(ov);
    g_pde_vs_ov->SetPoint(iov, ov, PDE(ov));
    g_gain_vs_ov->SetPoint(iov, ov, Gain(ov));
    g_tot_vs_amp->SetPoint(iov, amp, peakTot);
    g_amp_vs_tot->SetPoint(iov, peakTot, amp);
    g_amp_vs_tot->SetPointError(iov, peakTotErr, 0.01*amp); // error on amp is arbitrary (1%)
  }
  
  // -- get correction function for tot
  //TF1 *fCorr = new TF1("fCorr","expo", 0, 1000); // simple exponential
  //TF1 *fCorr = new TF1("fCorr","[0]*x + [1]*x*x + [2]*x*x*x", 0, 1000 ); // pol, passing by 0 for tot = 0
  TF1 *fCorr = new TF1("fCorr","[0] * ( exp([1]*x) - 1 )", 0, 1000); // shifted exponential, passing by 0 for tot = 0
  fCorr->SetParameter(0, 20.);
  fCorr->SetParameter(1, 0.015);
  //TF1 *fCorr = new TF1("fCorr","[0] * ( exp([1]*x) -  exp([2]*x))", 0, 1000); // sum of exponential, passing by 0 for tot = 0
  //fCorr->SetParameter(0, 300);
  //fCorr->SetParameter(1, 0.015);
  //fCorr->SetParameter(2, 0.001);
  //fCorr->SetParLimits(2, 0., 10000.);

  if ( vth2 < 30){
    g_amp_vs_tot-> RemovePoint(0);
  }

  g_amp_vs_tot-> Fit(fCorr,"QSR");

  TGraphErrors *g_totcorr_vs_ov = new TGraphErrors();
  g_totcorr_vs_ov->SetName("g_totcorr_vs_ov");
  g_totcorr_vs_ov->SetMarkerStyle(20);

  // -- apply Tot correction
  for (int i = 0; i < tree->GetEntries(); i++){

    tree -> GetEntry(i);

    if (channelID != channel) continue;
    
    for (int iov = 0; iov < Vovs.size(); iov++){
      if (ov == Vovs[iov]){
	float tot_corr = fCorr->Eval(tot/1000.);
	h_tot_corr[ov] -> Fill(tot_corr);
      }
    }
    
  }


  // -- fit 1275keV photopeak
  string comptonScattering2 = " 1./(1 + exp( [1] * (x -[2]) ) )";
  string emissionPeak2 = "[3] * exp(-pow(x-[4],2)/(2*[5]*[5]))";
  string backScattering2 = "[15] * exp(-pow(x-[16],2)/(2*[17]*[17]))";
  
  string comptonScattering1 = "[6]/(1 + exp( [7] * (x -[8]) ) )";
  string emissionPeak1 = "[9] * exp(-pow(x-[10],2)/(2*[11]*[11]))";
  string backScattering1 = "[12] * exp(-pow(x-[13],2)/(2*[14]*[14]))";
  
  string megafit = "[0]*(" + comptonScattering2 + " + "
    + emissionPeak2 + " + "
    + backScattering2 + " + "
    + comptonScattering1 + " + "
    + emissionPeak1 + " + "
    + backScattering1 + ")" ;

  cout << megafit.c_str() <<endl;

    
  TF1 *fGaus2[10];
  TF1 *f2[10];

  for (int iov = 0; iov < Vovs.size(); iov++){

    float ov = Vovs[iov];

    // -- fit 511 KeV photopeak
    fGaus2[iov] = new TF1(Form("fGaus2_%.02f",ov),"gaus", 0, 500);
    fGaus2[iov]->SetLineColor(800+20*iov);
    int lastbin = h_tot_corr[ov] -> FindLastBinAbove(10);
    float lastx = h_tot_corr[ov] -> GetBinCenter(lastbin);
    h_tot_corr[ov] -> GetXaxis()->SetRangeUser(lastx*0.3, lastx*10);
    int maxbin = h_tot_corr[ov] -> GetMaximumBin();
    double xmax = h_tot_corr[ov] -> GetBinCenter(maxbin);
    fGaus2[iov]->SetParameter(1,xmax);
    fGaus2[iov]->SetParameter(2,0.1*xmax);
    fGaus2[iov]->SetRange(xmax-0.1*xmax, xmax+0.1*xmax);
    h_tot_corr[ov] -> GetXaxis()->SetRangeUser(0, 10000); 
    h_tot_corr[ov] ->Fit( fGaus2[iov], "QR");
    float peak1      =  fGaus2[iov]-> GetParameter(1);
    float peak1sigma =  fGaus2[iov]-> GetParameter(2);
    //cout << "OV = "<< ov << "   xmax = " << xmax << "   gaus mean = " << peak1 << "   gaus sigma = "<< peak1sigma <<endl;

    f2[iov] = new TF1(Form("f2_%.02f",ov), megafit.c_str(), 0, 1000);
    f2[iov]->SetLineColor(800+20*iov);

    //f2[iov]->SetParameter(1, 0.1);
    //f2[iov]->SetParameter(7, 0.1);
    
    f2[iov]->SetParameter(2, peak1*2.5*0.9);
    f2[iov]->SetParameter(8, peak1*0.9);

    f2[iov]->SetParameter(10, peak1);
    f2[iov]->SetParameter(11, peak1sigma);

    f2[iov]->SetParameter(4, peak1*2.5);
    f2[iov]->SetParameter(5, peak1sigma/sqrt(2.5));

    // backscatter peak
    // f2[iov]->FixParameter(12, 0);
    f2[iov]->SetParameter(13, peak1*0.35);
    f2[iov]->SetParameter(14, peak1sigma/sqrt(0.35));
    //f2[iov]->FixParameter(15, 0);
    f2[iov]->SetParameter(16, peak1*0.35);
    f2[iov]->SetParameter(17, peak1sigma/sqrt(0.35));     
        

    if (ov==1.0){
      f2[iov]->SetRange(0.6*peak1, 10000);
      h_tot_corr[ov] ->Fit( f2[iov], "QRS");
    }
    else {
      f2[iov]->SetRange(0.2*peak1, 10000);
      h_tot_corr[ov] ->Fit( f2[iov], "QRS");
    }
    
    f2[iov]->SetRange(peak1-5.0*peak1sigma, 10000);
    TFitResultPtr r = h_tot_corr[ov]->Fit( f2[iov], "QRS");
    
    if (r!=0){
      r = h_tot_corr[ov] ->Fit( f2[iov], "QRSM");
    }
    
    cout << "OV = " << ov 
	 << "  511 keV peak = " << f2[iov]->GetParameter(10) << " +/- " << f2[iov]->GetParError(10) 
	 << " 1275 keV peak = " << f2[iov]->GetParameter(4) << " +/- " << f2[iov]->GetParError(4) 
	 << "      chi2/ndf = " << f2[iov]->GetChisquare()/f2[iov]-> GetNDF() 
	 << "      fitStatus = " << r  <<endl;
        
    float peakCorr1  = f2[iov]->GetParameter(10);
    float peakCorr2  = f2[iov]->GetParameter(4);
    float ratio = peakCorr2/peakCorr1;
    float err = ratio * sqrt(pow(f2[iov]->GetParError(10)/peakCorr1,2) + pow(f2[iov]->GetParError(4)/peakCorr2,2) );

    if ( r==0 && f2[iov]->GetChisquare()/f2[iov]-> GetNDF() < 2 ){
      g_totcorr_vs_ov->SetPoint(iov, ov, ratio);
      g_totcorr_vs_ov->SetPointError(iov, 0, err);
    }
     
  }
  
  g_totcorr_vs_ov->Fit("pol0");

  // == plots
  outdir = outdir + Form("vth1_%d", vth1) + "_" +  Form("vth2_%d", vth2) + "/";
  cout << "Saving plots in " << outdir.c_str() << endl;
  if (mkdir(outdir.c_str(), ACCESSPERMS) == 0){
    cout << "Not able to create directory " << outdir.c_str()  <<endl;  
  }
  
  TLatex *latex = new TLatex(0.15, 0.96, "LYSO:Ce 3x3x3 mm^{3} - HDR2");
  latex->SetNDC();
  latex->SetTextSize(0.030);

  
  TCanvas *c_pde = new TCanvas("c_pde","c_pde",700, 600);
  c_pde->SetGridx();
  c_pde->SetGridy();
  TH2F *hdummy0 = new TH2F("hdummy0","",100, 0, 6.5, 100, 0.0, 0.50);
  hdummy0->GetXaxis()->SetTitle("V_{OV}(V)");
  hdummy0->GetYaxis()->SetTitle("PDE");
  hdummy0->Draw();
  g_pde_vs_ov->Draw("plsame");
  c_pde->SaveAs((outdir+"c_pde_vs_ov.pdf").c_str());
  c_pde->SaveAs((outdir+"c_pde_vs_ov.png").c_str());

  TCanvas *c_gain = new TCanvas("c_gain","c_gain",700, 600);
  c_gain->SetGridx();
  c_gain->SetGridy();
  TH2F *hdummy1 = new TH2F("hdummy1","",100, 0, 6.5, 100, 0.0, 1000000);
  hdummy1->GetXaxis()->SetTitle("V_{OV}(V)");
  hdummy1->GetYaxis()->SetTitle("gain");
  hdummy1->Draw();
  g_gain_vs_ov->Draw("plsame");
  c_gain->SaveAs((outdir+"c_gain_vs_ov.pdf").c_str());
  c_gain->SaveAs((outdir+"c_gain_vs_ov.png").c_str());


  TCanvas *c_corr= new TCanvas();
  c_corr->SetGridx();
  c_corr->SetGridy();
  gStyle->SetOptFit(111);
  TH2F *hdummy2 = new TH2F("hdummy2","",100, 0, 350, 100, 0.0, 300.);
  hdummy2->GetXaxis()->SetTitle("measured ToT (ns)");
  hdummy2->GetYaxis()->SetTitle("amplitude (a.u.)");
  hdummy2->Draw();
  g_amp_vs_tot->Draw("psame");
  latex->Draw();
  gPad->Update();
  TPaveStats *st = (TPaveStats*)g_amp_vs_tot->FindObject("stats");
  st->SetX1NDC(0.65);
  st->SetX2NDC(0.85);
  st->SetY1NDC(0.77);
  st->SetY2NDC(0.92);  
  c_corr->SaveAs((outdir+"c_correction.pdf").c_str());
  c_corr->SaveAs((outdir+"c_correction.png").c_str());
  
  
  TCanvas *c_tot = new TCanvas("c_tot","c_tot",700, 600);
  c_tot->SetTickx();
  c_tot->SetTicky();
  TLegend *leg = new TLegend(0.6, 0.7, 0.80, 0.89);
  leg->SetBorderSize(0);
  for (int iov = 0; iov < Vovs.size(); iov++){
    h_tot[Vovs[iov]]->SetLineColor(800+20*iov);
    h_tot[Vovs[iov]]->GetYaxis()->SetRangeUser(0., h_tot[Vovs[iov]]->GetMaximum()*1.3);
    h_tot[Vovs[iov]]->GetXaxis()->SetRangeUser(0, 350);
    h_tot[Vovs[iov]]->GetXaxis()->SetTitle("ToT (ns)");
    leg->AddEntry(h_tot[Vovs[iov]], Form("V_{OV} = %0.1f",Vovs[iov]),"L");
    if (iov == 0) 
      h_tot[Vovs[iov]]->Draw("");
    else
      h_tot[Vovs[iov]]->Draw("same");
  }
  latex->Draw();
  leg->Draw();
  c_tot->SaveAs((outdir+"c_tot.png").c_str());
  c_tot->SaveAs((outdir+"c_tot.pdf").c_str());



  TCanvas *c_tot_corr_all = new TCanvas("c_tot_corr_all","c_tot_corr_all",700, 600);
  c_tot_corr_all->SetTickx();
  c_tot_corr_all->SetTicky();
  TLegend *leg2 = new TLegend(0.6, 0.7, 0.80, 0.89);
  leg2->SetBorderSize(0);
  for (int iov = 0; iov < Vovs.size(); iov++){
    h_tot_corr[Vovs[iov]]->SetLineColor(800+20*iov);
    h_tot_corr[Vovs[iov]]->GetYaxis()->SetRangeUser(0., h_tot_corr[Vovs[iov]]->GetMaximum()*1.1);
    h_tot_corr[Vovs[iov]]->GetXaxis()->SetRangeUser(0, 700);
    h_tot_corr[Vovs[iov]]->GetXaxis()->SetTitle("amplitude (a.u.)");
    leg2->AddEntry(h_tot_corr[Vovs[iov]], Form("V_{OV} = %0.1f",Vovs[iov]),"L");
    if (iov == 0) 
      h_tot_corr[Vovs[iov]]->Draw("");
    else
      h_tot_corr[Vovs[iov]]->Draw("same");
  }
  latex->Draw();
  leg2->Draw();
  c_tot_corr_all->SaveAs((outdir+"c_tot_corr.png").c_str());
  c_tot_corr_all->SaveAs((outdir+"c_tot_corr.pdf").c_str());


  
  TCanvas *c_tot_corr[10];
  for (int iov = 0; iov < Vovs.size(); iov++){
    c_tot_corr[iov] = new TCanvas(Form("c_tot_corr_Vov%.01f",Vovs[iov]),Form("c_tot_corr_Vov%.01f",Vovs[iov]),700, 600);
    c_tot_corr[iov]->SetLogy();
    h_tot_corr[Vovs[iov]]->GetYaxis()->SetRangeUser(3, h_tot_corr[Vovs[iov]]->GetMaximum()*5);
    h_tot_corr[Vovs[iov]]->GetFunction( Form("f2_%.02f",Vovs[iov]) )-> SetLineColor(1);
    //float xmin  = h_tot_corr[Vovs[iov]]->GetFunction(Form("f2_%.02f",Vovs[iov]) )->GetParameter(10);
    //float xmax  = h_tot_corr[Vovs[iov]]->GetFunction(Form("f2_%.02f",Vovs[iov]) )->GetParameter(4);
    //h_tot_corr[Vovs[iov]]->GetXaxis()->SetRangeUser(xmin*0.6, xmax*1.2);
    int binmax = h_tot_corr[Vovs[iov]]-> FindLastBinAbove(2);
    float xmax = h_tot_corr[Vovs[iov]]-> GetBinCenter(binmax) + 10;
    h_tot_corr[Vovs[iov]]->GetXaxis()->SetRangeUser(0, xmax);
    h_tot_corr[Vovs[iov]]->GetXaxis()->SetTitle("amplitude (a.u.)");
    h_tot_corr[Vovs[iov]]->Draw("");
    latex->Draw();
    c_tot_corr[iov]->SaveAs((outdir+Form("c_tot_corr_Vov%.01f.pdf",Vovs[iov])).c_str());
    c_tot_corr[iov]->SaveAs((outdir+Form("c_tot_corr_Vov%.01f.png",Vovs[iov])).c_str());
  }
  
  
  TCanvas *c_ratio = new TCanvas("c_ratio","c_ratio", 700, 600);
  c_ratio->SetGridx();
  c_ratio->SetGridy();
  TH2F *hdummy = new TH2F("hdummy","",100, 0, 6.5, 100, 1.0, 4.0);
  hdummy->GetXaxis()->SetTitle("V_{OV}(V)");
  hdummy->GetYaxis()->SetTitle("peak_{1275keV}/peak_{511keV}");
  hdummy->Draw();
  latex->Draw();
  g_totcorr_vs_ov->Draw("psame");
  c_ratio->SaveAs((outdir+"c_ratio.pdf").c_str());
  c_ratio->SaveAs((outdir+"c_ratio.png").c_str());
  



  // == save histograms
  TFile *fout = new TFile(Form("output_TotNonlinearity_vth1_%d_vth2_%d.root", vth1, vth2), "recreate");
  for (int iov = 0; iov < Vovs.size(); iov++){
    h_tot[Vovs[iov]]->Write();
    h_tot_corr[Vovs[iov]]->Write();
  }
  g_pde_vs_ov->Write();
  g_tot_vs_amp->Write();
  g_amp_vs_tot->Write();
  g_totcorr_vs_ov->Write();
 
  fout->Close();
  
}


    
