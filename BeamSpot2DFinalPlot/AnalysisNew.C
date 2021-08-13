#define AnalysisNew_cxx
#include "AnalysisNew.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalysisNew::Loop()
{
//   In a ROOT session, you can do:
//      root> .L AnalysisNew.C
//      root> AnalysisNew t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   TString root;
   float mean_Wsum=0;
   float sigma_Wsum=0;
   float sigma_Invsum=0;
   Float_t gmean[300],gmean_Min[300],gmean_Max[300],gmeanY[300],gmean_MinY[300],gmean_MaxY[300] ,grms[300],gFill[300],gdxy[300],gmeanErr[300],grmsErr[300];
   Float_t gmeanX1[300],gmeanXErr1[300],gmeanY1[300],gmeanYErr1[300];
   Float_t gmeanX2[300],gmeanXErr2[300],gmeanY2[300],gmeanYErr2[300];
   Float_t grad1[300], grErr[300];
   int i=0;
   int hi=0;
   int him=0;
   gROOT->cd();
   float x,y;
   TNtuple ntuple("ntuple","data from ascii file","x:y");
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (nFill < 4879 ) continue;//low lumi begining
      if (nFill > 5456 ) continue;//HI run

      //if (nFill < 5150 ) continue;
      //if (nFill > 5250 ) continue;
      //if (nFill < 5024 ) continue;
      //if (nFill > 5183) continue;
      //if (nFill==5196) continue;//no lui fill
      
      //if (nFill == 5424) continue;
      //if (nFill == 5251) continue;
      //if (nFill == 5101) continue;
      if (meanErrX>0.0008) continue;
      if (meanErrY>0.0008) continue;
      
      //	z==0 cuts
      if (meanErrX>0.0004) continue;
      if (meanErrY>0.0004) continue;
      if (meanY>0.02) continue;
      if (meanY<-0.035) continue;
      

      gdxy[i]=0.5;
      gFill[i]=nFill;
      if (abs(meanX)> 0.5){
	//cout<<" The Numbers are"<<gFill[i]<<endl;
      }
      if (sigmaY < 0.25){
	//cout<<" The Numbers less are"<<gFill[i]<<endl;
      }
      
      gmean[i]=meanX;
      gmean_Min[i]=meanX-sigmaX;
      gmean_Max[i]=meanX+sigmaX;
      
      gmeanY[i]=meanY;
      gmean_MinY[i]=meanY-sigmaY;
      gmean_MaxY[i]=meanY+sigmaY;

     
      gmeanErr[i]=meanErrX;
      grmsErr[i]=meanErrY;
     
      if (nFill >= 5024 && nFill <= 5183){
	  //cout<<"The Fill is "<<nFill<<endl;
	  gmeanX1[hi]=meanX;
	  gmeanXErr1[hi]=meanErrX;
	  gmeanY1[hi]=meanY;
	  gmeanYErr1[hi]=meanErrY;
	  //cout<<nFill<<endl;
	  hi++;
      }
      /*
      if (nFill <= 5005 || nFill >= 5521){
      //if (nFill <= 5005 || nFill >= 5424){
	  //cout<<"The Fill is "<<nFill<<endl;
	  gmeanX2[him]=meanX;
	  gmeanXErr2[him]=meanErrX;
	  gmeanY2[him]=meanY;
	  gmeanYErr2[him]=meanErrY;
	  cout<<nFill<<endl;
	  him++;
	}
      */
      if (meanX >= -0.02 && meanX <= 0.00){
	if (meanY >= 0.0 && meanY <= 0.02){
	  cout<<"The Fill for green is "<<nFill<<endl;
	  gmeanX2[him]=meanX;
	  gmeanXErr2[him]=meanErrX;
	  gmeanY2[him]=meanY;
	  gmeanYErr2[him]=meanErrY;
	  him++;
	}
      }
      
      grms[i]=sigmaX;
      mean_Wsum = mean_Wsum + (meanX/(sigmaX*sigmaX));
      sigma_Wsum = sigma_Wsum + (sigmaX/(sigmaX*sigmaX));
      sigma_Invsum = sigma_Invsum + (1 / (sigmaX*sigmaX));
      //r calc
      grad1[i]=pow((meanX*meanX+meanY*meanY),0.5);
      grErr[i]=pow((meanErrX*meanErrX + meanErrY*meanErrY),0.5);
      
      ntuple.Fill(meanX,meanY);
      i++;
      
   }
   //cout<<"The mean is "<<meantot/i<<endl;

   cout<<"The number of fill going to plot are  "<<i<<endl;
   
   //TCanvas *c1 = new TCanvas ("c1","Graph Draw Options",1200,800);
   TCanvas *c1 = new TCanvas ("c1","Graph Draw Options",1000,1000);//for 2D
   //For 2D Plot
   TGraphErrors *gr = new TGraphErrors(i,gmean, gmeanY,gmeanErr,grmsErr);
   TGraphErrors *gr1 = new TGraphErrors(hi, gmeanX1, gmeanY1,gmeanXErr1,gmeanYErr1);
   TGraphErrors *gr2 = new TGraphErrors(him, gmeanX2, gmeanY2,gmeanXErr2,gmeanYErr2);
   
   c1->cd(4);
   //plot styles
   gROOT->SetStyle("Plain");   
   //gStyle->SetOptStat("kFALSE");
   gStyle->SetOptStat(0);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetGrid();
   c1->SetTitle("");
   c1->Modified();

   gStyle->SetOptTitle(kFALSE);
   gPad->SetRightMargin(15/100.);
   gPad->SetLeftMargin(15/100.);
   gPad->SetTopMargin(15/100.);
   gPad->SetBottomMargin(15/100.);
   gStyle->SetPalette(1);
   
   //cout<<"The mean is "<<gr->GetMean(2)<<endl;//2 for y axis
   //cout<<"The Std Dev is "<<gr->GetStdDev(2)<<endl;
   //cout<<"The RMS is "<<gr->GetRMS(2)<<endl;

   gr->SetTitle("");
   gr->SetMarkerColor(1);
   gr->SetMarkerStyle(8);
   gr->SetLineWidth(3);
   gr1->SetLineWidth(3);
   gr1->SetMarkerColor(2);
   gr1->SetMarkerStyle(8);
   gr2->SetLineWidth(3);
   gr2->SetMarkerColor(3);
   gr2->SetMarkerStyle(8);
   gr2->SetLineWidth(3);
   
   gr->SetLineColor(1);
   gr1->SetLineColor(4);
   gr2->SetLineColor(3);
   gr2->SetLineColor(3);
   
   gr->GetYaxis()->SetRangeUser(-0.05,0.05);//for z=0
   gr->GetYaxis()->SetTitleOffset(1.4);
   gr->GetYaxis()->SetTitleSize(0.05);
   gr->GetYaxis()->SetNdivisions(506);
   gr->GetYaxis()->SetTitle("Beamspot y|_{z=0} [cm]");//2D
   

   gr->GetXaxis()->SetLimits(-0.05,0.05);
   gr->GetXaxis()->SetTitleOffset(1.0); 
   gr->GetXaxis()->SetTitleSize(0.05);
   gr->GetXaxis()->SetNdivisions(506);//2D
   gr->GetXaxis()->SetTitle("Beamspot x|_{z=0} [cm]");//2d
   
   
   gr->Draw("AP");
   gr1->Draw("P");
   gr2->Draw("P");
   
      
   TLatex latex;
   //latex.SetTextAlign(32);
   latex.SetTextFont(61);
   latex.DrawTextNDC(.15,.87,"CMS");
   latex.SetTextFont(52);
   latex.DrawTextNDC(.27,.87,"Preliminary");
   latex.SetTextFont(42);
   latex.DrawTextNDC(.58,.87,"2016 (13 TeV)");

   
   TEllipse *el7 = new TEllipse(0.006,-0.004,.03,.03);
   el7->SetFillStyle(0);
   el7->SetLineWidth(2);
   el7->SetLineStyle(7);
   el7->Draw();
      
}
