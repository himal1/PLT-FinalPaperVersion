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
   Float_t gFill1[300], gFill2[300], gdxy1[300],gdxy2[300], gmeanX3[300], gmeanY3[300], gdxy3[300], gFill3[300]; 
   Float_t gmean[300],gmean_Min[300],gmean_Max[300],gmeanY[300],gmean_MinY[300],gmean_MaxY[300] ,grms[300],gFill[300],gdxy[300],gmeanErr[300],grmsErr[300];
   Float_t gmeanX1[300],gmeanXErr1[300],gmeanY1[300],gmeanYErr1[300], gmeanYErr3[300], gmeanXErr3[300];
   Float_t gmeanX2[300],gmeanXErr2[300],gmeanY2[300],gmeanYErr2[300];
   Float_t grad1[300], grErr[300];
   int i=0;
   int hi=0;
   int him=0;
   int hima =0;
   gROOT->cd();
   float x,y;
   TNtuple ntuple("ntuple","data from ascii file","x:y");
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (nFill < 4879 ) continue;//low lumi begining
      //if (nFill > 5183 ) continue;//HI run

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
      
      //cut for track
      //if (meanErrX>0.004) continue;
      //if (meanErrY>0.004) continue;

      //if (meanY>0.08) continue;
      //if (meanX>0.08) continue;
      //if (meanY< -0.1) continue;
      //if (meanX< -0.1) continue;
      
      gdxy[i]=0.5;
      gFill[i]=nFill;
      if (abs(meanX)> 0.5){
	cout<<" The Numbers are"<<gFill[i]<<endl;
      }
      if (sigmaY < 0.25){
	cout<<" The Numbers less are"<<gFill[i]<<endl;
      }
      //gmean[n]=meanY;
      //gmean[n]=r;
      gmean[i]=meanX;
      gmean_Min[i]=meanX-sigmaX;
      gmean_Max[i]=meanX+sigmaX;
      
      gmeanY[i]=meanY;
      gmean_MinY[i]=meanY-sigmaY;
      gmean_MaxY[i]=meanY+sigmaY;

     
      //meantot=meantot+meanX;
      //gmeanErr[n]=meanErrY;
      //gmeanErr[i]=meanErrX;
      gmeanErr[i]=meanErrX;
      grmsErr[i]=meanErrY;
      /* 
      if (meanX >= 0.0 && meanX <= 0.02){
	if (meanY >= -0.03 && meanY <= -0.01){
	  //cout<<"The Fill is "<<nFill<<endl;
	  gmeanX1[hi]=meanX;
	  gmeanXErr1[hi]=meanErrX;
	  gmeanY1[hi]=meanY;
	  gmeanYErr1[hi]=meanErrY;
	  cout<<nFill<<endl;
	  hi++;
	}
      }
      */
      if (nFill >= 5024 && nFill <= 5183){
	  //cout<<"The Fill is "<<nFill<<endl;
	  gmeanX1[hi]=meanX;
	  gmeanXErr1[hi]=meanErrX;
	  gmeanY1[hi]=meanY;
	  gmeanYErr1[hi]=meanErrY;
	  gFill1[hi]=nFill;
	  gdxy1[hi]=0.5;
	  cout<<nFill<<endl;
	  hi++;
	}
      
      if (meanX >= -0.02 && meanX <= 0.00){
	if (meanY >= 0.0 && meanY <= 0.02){
	  //cout<<"The Fill is "<<nFill<<endl;
	  gmeanX2[him]=meanX;
	  gmeanXErr2[him]=meanErrX;
	  gmeanY2[him]=meanY;
	  gmeanYErr2[him]=meanErrY;
	  gFill2[him]=nFill;
	  gdxy2[him]=0.5;
	  him++;
	}
      }

      if (nFill >= 5456){
	//cout<<"The Fill is "<<nFill<<endl;
	gmeanX3[hima]=meanX;
	gmeanXErr3[hima]=meanErrX;
	gmeanY3[hima]=meanY;
	gmeanYErr3[hima]=meanErrY;
	gFill3[hima]=nFill;
	gdxy3[hima]=0.5;
	cout<<nFill<<"He he included"<<endl;
	hima++;
      }
      grms[i]=sigmaX;
      mean_Wsum = mean_Wsum + (meanX/(sigmaX*sigmaX));
      sigma_Wsum = sigma_Wsum + (sigmaX/(sigmaX*sigmaX));
      sigma_Invsum = sigma_Invsum + (1 / (sigmaX*sigmaX));
      //r calc
      grad1[i]=pow((meanX*meanX+meanY*meanY),0.5);
      grErr[i]=pow((meanErrX*meanErrX + meanErrY*meanErrY),0.5);
      //gr[i]=r;
      //grErr[i]=rErr;
      ntuple.Fill(meanX,meanY);
      i++;
      
   }
   //cout<<"The mean is "<<meantot/i<<endl;

   cout<<"The number of fill going to plot are  "<<i<<endl;
   cout<<"The weighted mean is   "<<mean_Wsum/sigma_Invsum<<endl;
   cout<<"The weighted sigma is   "<<sigma_Wsum/sigma_Invsum<<endl;
   
   TCanvas *c1 = new TCanvas ("c1","Graph Draw Options",1200,800);
   //TCanvas *c1 = new TCanvas ("c1","Graph Draw Options",1000,1000);//for 2D
   //ntuple.Draw("y>>K(10,-0.03,0.01)","","");
   //for time evolution plot
   //mean plot for Y
   
   TGraphErrors *gr = new TGraphErrors(i, gFill, gmeanY,gdxy, gmeanErr);
   TGraphErrors *gr1 = new TGraphErrors(hi, gFill1, gmeanY1,gdxy1, gmeanYErr1);
   TGraphErrors *gr2 = new TGraphErrors(him, gFill2, gmeanY2,gdxy2, gmeanYErr2);
   TGraphErrors *gr3 = new TGraphErrors(hima, gFill3, gmeanY3,gdxy3, gmeanYErr3);
   //mean plot for X
   TGraphErrors *gr4 = new TGraphErrors(i, gFill, gmean,gdxy, gmeanErr);
   TGraphErrors *gr5 = new TGraphErrors(hi, gFill1, gmeanX1,gdxy1, gmeanXErr1);
   TGraphErrors *gr6 = new TGraphErrors(him, gFill2, gmeanX2,gdxy2, gmeanXErr2);
   TGraphErrors *gr7 = new TGraphErrors(hima, gFill3, gmeanX3,gdxy3, gmeanXErr3);
    
   c1->cd(4);
   //plot styles
   gROOT->SetStyle("Plain");   
   //gStyle->SetOptStat("kFALSE");
   gStyle->SetOptStat(0);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetTitle("");
   c1->Modified();

   gStyle->SetOptTitle(kFALSE);
   gPad->SetRightMargin(15/100.);
   gPad->SetLeftMargin(15/100.);
   gPad->SetTopMargin(15/100.);
   gPad->SetBottomMargin(15/100.);
   gStyle->SetPalette(1);
   
   cout<<"The mean is "<<gr->GetMean(2)<<endl;//2 for y axis
   //cout<<"The Std Dev is "<<gr->GetStdDev(2)<<endl;
   cout<<"The RMS is "<<gr->GetRMS(2)<<endl;

   gr->SetTitle("");
   gr->SetMarkerColor(1);
   gr->SetMarkerStyle(21);
   gr->SetLineWidth(3);
   gr->GetYaxis()->SetRangeUser(-0.05,0.05);//for z=0
   gr->GetYaxis()->SetTitleOffset(1.5);
   gr->GetYaxis()->SetNdivisions(506); 
   gr->SetLineColor(1);
   gr->GetXaxis()->SetTitleOffset(1.0); 
   gr->GetXaxis()->SetNdivisions(506);//2D
   gr->GetXaxis()->SetTitle("Fill");//2D
   gr->GetYaxis()->SetTitle("Beamspot x(y)|_{z=0} [cm]");//2D
   
   gr1->SetLineWidth(3);
   gr1->SetMarkerColor(2);
   gr1->SetMarkerStyle(21);
   gr1->SetLineColor(2);
   
   gr2->SetLineWidth(3);
   gr2->SetMarkerColor(3);
   gr2->SetMarkerStyle(21);
   gr2->SetLineColor(3);
   
   gr3->SetLineWidth(3);
   gr3->SetMarkerColor(4);
   gr3->SetMarkerStyle(21);
   gr3->SetLineColor(4);

   gr4->SetLineWidth(3);
   gr4->SetMarkerColor(1);
   gr4->SetMarkerStyle(8);
   gr4->SetLineColor(1);

   gr5->SetLineWidth(3);
   gr5->SetMarkerColor(2);
   gr5->SetMarkerStyle(8);
   gr5->SetLineColor(2);

   gr6->SetLineWidth(3);
   gr6->SetMarkerColor(3);
   gr6->SetMarkerStyle(8);
   gr6->SetLineColor(3);

   gr7->SetLineWidth(3);
   gr7->SetMarkerColor(4);
   gr7->SetMarkerStyle(8);
   gr7->SetLineColor(4);
   
   gr->Draw("AP");
   gr1->Draw("P");
   gr2->Draw("P");
   gr3->Draw("P");
   gr4->Draw("P");
   gr5->Draw("P");
   gr6->Draw("P");
   gr7->Draw("P");
 
         
   TLatex latex;
   //latex.SetTextAlign(32);
   latex.SetTextFont(61);
   latex.DrawTextNDC(.15,.87,"CMS");
   latex.SetTextFont(52);
   latex.DrawTextNDC(.27,.87,"Preliminary");
   latex.SetTextFont(42);
   latex.DrawTextNDC(.58,.87,"2016 (13 TeV)");

   //-------------TLegend WOrk---------------
   TLegend *leg = new TLegend(0.70,0.72,0.80,0.80);
   //TLegend *leg1 = new TLegend(0.45,0.75,0.65,0.85);
   leg->SetFillColor(0);
   leg->SetTextAlign(23);
   leg->SetTextFont(42);
   //leg1->SetFillColor(0);
   //leg1->SetTextAlign(23);
   //leg1->SetTextFont(42);
   //leg->SetHeader("Data");
   //leg->AddEntry(gr,"mean_{x} (#bar{x})","l");
   //for r
   
   leg->AddEntry(gr4,"mean_{x}","lep");
   leg->AddEntry(gr,"mean_{y}","lep");
   //leg->AddEntry(gr,"Region 3","lep");
   //leg->AddEntry(gr1,"#bar{x} - #sigma_{x}","l");
   //leg->AddEntry(gr2,"#bar{x} + #sigma_{x}","l");
   //leg->AddEntry(gr3,"mean_{y} (#bar{y})","l");
   //leg1->AddEntry(gr4,"#bar{y} - #sigma_{y}","l");
   //leg1->AddEntry(gr5,"#bar{y} + #sigma_{y}","l");
   leg->SetBorderSize(0);
   leg->Draw();
   //leg1->SetBorderSize(0);
   //leg1->Draw();
   TEllipse *el7 = new TEllipse(0.006,-0.004,.03,.03);
   el7->SetFillStyle(0);
   // el7->Draw();
   //z=0 calc
   TLine *l = new TLine (5456,-0.05,5456,0.025);
   TLine *l1 = new TLine (5024,-0.05,5024,0.025);
   TLine *l2 = new TLine (5183,-0.05,5183,0.025);
   
   /*
   TLine *l = new TLine (5456,-0.2,5456,0.1);
   TLine *l1 = new TLine (5024,-0.2,5024,0.1);
   TLine *l2 = new TLine (5183,-0.2,5183,0.1);
   */
   l->SetLineStyle(2);
   l1->SetLineStyle(2);
   l2->SetLineStyle(2);
   l1->SetLineColor(2);
   l2->SetLineColor(2);
   l->Draw();
   //l1->Draw();
   //l2->Draw();
   
}
