#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <stdio.h>

#include "TRandom.h"
#include "TH1D.h"

#include <TROOT.h>
#include <TStyle.h>
#include "TError.h"

#include "TFile.h"
#include "TCanvas.h"

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"


#include "TDirectory.h"
#include "TDirectoryFile.h"

#include "TGraph.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLine.h"
#include "TTree.h"
#include "TF1.h"
#include "TLegend.h"



#endif

//HLT_Jet60 for data only

//add pileup cut 
//change lower pt limit to 80 GeV/c 
int GetCentBin(int bin);
void drawCentTex(const int centbin);
void drawText(const char *text, float xp, float yp, int size);

void UltimateMCjetQAv2(int algo=7, int month=13, int date=1)  {
	

	
	cout<<" JetID & Quality Check Updated G. Pikul 03 24 14 "<<endl;
	cout<<" ================================================="<<endl;
	
	const char *macroname = "UltimateMCjetQAv2";

	bool reweigh = false;
	bool AAFlag=false;

	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|
	// 		ALGORITHM OPTIONS DEFINED
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|
	const char *algorithm[30] = {"akVs1PF", "akVs2PF", "akVs3PF", "akVs4PF", "akVs5PF",

				  "akPu1PF", "akPu2PF", "akPu3PF", "akPu4PF", "akPu5PF",
	 
				  "ak1PF", "ak2PF", "ak3PF", "ak4PF", "ak5PF",
				
				  "akPu1Calo", "akPu2Calo", "akPu3Calo", "akPu4Calo", "akPu5Calo",
				
				  "ak1Calo", "ak2Calo", "ak3Calo", "ak4Calo", "ak5Calo",
				
				  "icPu1PF", "icPu2PF", "icPu3PF", "icPu4PF", "icPu5PF"}; 

	cout << "algorithm used: " << algorithm[algo] << endl;			  
	//if algo == 0 akPu1PF used, if algo == 1 akPu2PF used, etc.

	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|
	// 		DEFINE SOME COMMONLY USED CUTS
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|

	float global_jtpt_cut = 0;
	float leading_jtpt_cut = 0;
	float jteta_cut = 2;

	
	//~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~
	//BE CAREFUL!!!!!!!
	//NEED TO CHANGE HLT TRIGGER CUT BEFORE RUN WITH CHANGED FILE NAME
	//~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~/~~~
	
	
	//// OutPut file
	//  TFile *f = new TFile("pbpb_Unfo.root","RECREATE");
	
	
	const int NBINS=21;
	Double_t BOUNDARIES[NBINS] = {100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250 ,260, 270, 280, 290, 300};
	
	const int NBINS_=26;
	Double_t BOUNDARIES_[NBINS_] = {50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250 ,260, 270, 280, 290, 300};
	
	const int ncen=6;
	
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	
	
	//! Centrality reweighting function
	TF1* fcen= new TF1("fcen","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
	fcen->SetParameters(2.10653e-02,5.61607,-1.41493e-01,1.00586e-03,-1.32625e-04);
	
	//! Centrality reweighting function
 	// TF1* fcen = new TF1("fcen","exp(-1.0*pow(x+1.11957e+01,2)/pow(1.34120e+01,2)/2)",0,40); 
	
	//! Primary Vertex reweighting function
	TF1 *fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
	//fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05);
	//fVz->SetParameters(6.08192e-02,6.32282e-04, -5.48146e-04,-2.74674e-06 ,1.33737e-06);

	//fVz->SetParameters(8.33922e-01,1.62567e-03,2.91534e-03,-1.14466e-04,5.69407e-05);
	//fVz->SetParameters(8.41697e-01,-6.90553e-03,2.48080e-03,-9.94636e-05,2.57340e-05);
	fVz->SetParameters(1.00991,9.39128e-03,-2.98879e-04,2.69900e-06,9.23987e-07);

	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//
	//		CROSS SECTIONS <--for cross section reweighting 
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//

	//values for PbPb
	float xsection1 = (1.079e-02)-(1.021e-03); //dijet30
	float xsection2 = (1.021e-03)-(9.913e-05); //dijet50
	float xsection3 = (9.913e-05)-(3.069e-05); //dijet80
	float xsection4 = (3.069e-05)-(1.128e-05); //dijet100
	float xsection5 = (1.128e-05); //dijet120
	

	/*float xsection1 = (1.079e-02)-(1.021e-03); //dijet30
	float xsection2 = (1.021e-03)-(9.913e-05); //dijet50
	float xsection3 = (9.913e-05)-(3.069e-05); //dijet80
	float xsection4 = (3.069e-05)-(1.128e-05); //dijet100
	float xsection5 = (1.128e-05)-(1.470e-06); //dijet120
	float xsection6 = (1.470e-06)-(5.310e-07); //dijet170
	float xsection7 = (5.310e-07)-(1.192e-07); //dijet200
	float xsection8 = (1.192e-07)-(3.176e-07); //dijet250
	float xsection9 = (3.176e-07); //dijet300*/

	//values for pp
	/*float xsection1 = (2.034e-01)-(1.075e-02); //dijet15
	float xsection2 = (1.075e-02)-(1.025e-03); //dijet30
	float xsection3 = (1.025e-03)-(9.865e-05); //dijet50
	float xsection4 = (9.865e-05)-(1.129e-05); //dijet80
	float xsection5 = (1.129e-05)-(1.465e-06); //dijet120
	float xsection6 = (1.465e-06)-(2.837e-07); //dijet170
	float xsection7 = (2.837e-07)-(5.323e-08); //dijet220
	float xsection8 = (5.323e-08); //dijet280*/
	
	//	fill in appropriate values here

	//values for Pbp

	//	fill in appropriate values here

	
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//
	//		CREATE HISTOGRAMS
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//

	TH1F* hpT[6];
	TH1F* heta[6];
	TH1F* hphi[6];

	TH1F* hphotonMax[6];
	TH1F* hphotonSum[6];
	TH1F* htrackMax[6];
	TH1F* htrackSum[6];
	TH1F* hneutralMax[6];
	TH1F* hneutralSum[6];
	TH1F* hchargedMax[6];

	TH1F* hvx[6];
	TH1F* hvy[6];
	TH1F* hvz[6];
	TH1F* hbin[6];

	
	TH2F* heta_pT[6];
	TH2F* heta_phi[6];
	
	TH1F* hpthat[6];
	
	
	
	for(int icen=0;icen<ncen;icen++)	{

		//hpT[icen]= new TH1F (Form("hpT%d",icen), Form("hpT%d; Jet p_{T}; entry fraction ",icen),   NBINS_-1,BOUNDARIES_);
		hpT[icen]= new TH1F (Form("hpT%d",icen), Form("hpT%d; jet p_{T} (GeV/c); entry fraction ",icen), 45, 50, 500);
	    heta[icen]= new TH1F (Form("heta%d",icen), Form("heta%d; Jet #eta; entry fraction",icen), 20,-2,+2);
	    hphi[icen]= new TH1F (Form("hphi%d",icen), Form("hphi%d; Jet #phi; entry fraction",icen), 21, -3.5,+3.5);
	    
		hphotonMax[icen] = new TH1F (Form("hphotonMax%d", icen), Form("photonMax%d; max photon p_{T} / jet p_{T}; entry fraction",icen), 20, 0.,1. );
	    hphotonSum[icen]= new TH1F (Form("hphotonSum%d",icen), Form("photonSum%d; sum photon p_{T} / jet p_{T}; entry fraction",icen),   20, 0.,1.);
	    htrackMax[icen]= new TH1F (Form("htrackMax%d",icen), Form("trackMax%d; max track p_{T} / jet p_{T}; entry fraction ",icen),   20, 0.,1.); 
	    htrackSum[icen]= new TH1F (Form("htrackSum%d",icen), Form("trackSum%d; sum track p_{T} / jet p_{T}; entry fraction",icen),   20, 0.,1.); 
	    hneutralMax[icen]= new TH1F (Form("hneutralMax%d",icen), Form("neutralMax%d; max neutral p_{T} / jet p_{T};entry fraction ",icen),  20, 0.,1.); 
	    hneutralSum[icen]= new TH1F (Form("hneutralSum%d",icen), Form("neutralSum%d; sum neutral p_{T} / jet p_{T}; entry fraction ",icen),   20, 0.,1.); 
	    hchargedMax[icen]= new TH1F (Form("hchargedMax%d",icen), Form("hchargedMax%d; max charged p_{T} / jet p_{T}; entry fraction ",icen),   20, 0.,1.); 
	    
	    hvx[icen] = new TH1F (Form("hvx%d",icen), Form("hvx%d; vx (cm); entry fraction ",icen), 40, -20, 20); 
	    hvy[icen] = new TH1F (Form("hvy%d",icen), Form("hvy%d; vy (cm); entry fraction ",icen), 40, -20, 20); 
	    hvz[icen] = new TH1F (Form("hvz%d",icen), Form("hvz%d; vz (cm); entry fraction ",icen), 40, -20, 20);

	    hbin[icen] = new TH1F (Form("hbin%d",icen), Form("hbin%d; centrality bin; entry fraction ",icen), 40, 0, 40);  

	    heta_phi[icen]= new TH2F (Form("heta_phi%d",icen), Form("eta_phi%d; Jet #eta; Jet #phi",icen),20,-2,+2, 20, -3.14,+3.14);
	    heta_pT[icen]= new TH2F (Form("heta_pT%d",icen), Form("eta_phi%d; Jet #eta; Jet p_{T} ",icen),20,-2,+2, 20, 50,500);
	   
	    hpthat[icen]= new TH1F (Form("hpthat%d",icen), Form("hpthat%d; pthat; entry fraction ",icen),   NBINS_-1,BOUNDARIES_);
	    
	  }
	
	// *****************  Ratio Plots *****************
	
	
	TH1F* hphotonMax_Ratio[6];
	TH1F* htrackMax_Ratio[6];
	TH1F* htrackSum_Ratio[6];
	TH1F* hneutralMax_Ratio[6];
	TH1F* hneutralSum_Ratio[6];
	TH1F* hphotonSum_Ratio[6];
	TH1F* hpT_Ratio[6];
	TH1F* heta_Ratio[6];
	TH1F* hphi_Ratio[6];
	
	
	for(int icen=0;icen<ncen;icen++)
	  {   hphotonMax_Ratio[icen] = new TH1F (Form("hphotonMax_Ratio%d",icen), Form("photonMax%d; Max photon p_{T} / Jet p_{T}; Data / MC Ratio",icen), 20, 0.,1. );
	    htrackMax_Ratio[icen]= new TH1F (Form("htrackMax_Ratio%d",icen), Form("trackMax%d; Max track p_{T} / Jet p_{T}; Data / MC Ratio ",icen),   20, 0.,1.); 
	    htrackSum_Ratio[icen]= new TH1F (Form("htrackSum_Ratio%d",icen), Form("trackSum%d; Sum track p_{T} / Jet p_{T};  Data / MC Ratio ",icen),   20, 0.,1.); 
	    hneutralMax_Ratio[icen]= new TH1F (Form("hneutralMax_Ratio%d",icen), Form("neutralMax%d; Max neutral p_{T} / Jet p_{T}; Data / MC Ratio ",icen),  20, 0.,1.); 
	    hneutralSum_Ratio[icen]= new TH1F (Form("hneutralSum_Ratio%d",icen), Form("neutralSum%d; Sum neutral p_{T} / Jet p_{T};  Data / MC Ratio ",icen),   20, 0.,1.); 
	    hphotonSum_Ratio[icen]= new TH1F (Form("hphotonSum_Ratio%d",icen), Form("photonSum%d; Sum photon p_{T} / Jet p_{T};  Data / MC Ratio ",icen),   20, 0.,1.);
	    hpT_Ratio[icen]= new TH1F (Form("hpT_Ratio%d",icen), Form("hpT%d; Jet p_{T};  Data / MC Ration ",icen),   NBINS_-1,BOUNDARIES_);
	    heta_Ratio[icen]= new TH1F (Form("heta_Ratio%d",icen), Form("heta%d; Jet #eta;  Data / MC Ratio ",icen), 20,-2,+2);
	    hphi_Ratio[icen]= new TH1F (Form("hphi_Ratio%d",icen), Form("hphi%d; Jet #phi;  Data / MC Ratio",icen), 20, -3.14,+3.14);
	    
	  }
	
	
	
	
	
	
	
	
	
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|
	//
	//		*****	READ MC FILES ONE-BY-ONE	*****
	//
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|


	//|\\|//|\\|//|\\|//|\\|//|\\|
	//		READ MC 30 FILE
	//|\\|//|\\|//|\\|//|\\|//|\\|
	
	cout<<"read MC 30 file"<<endl;
	
	TString inname1="/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet30_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root";
	//TFile *file1 = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt15/HiForest_v81_merged01/pt15_pp2013_P01_prod22_v81_merged_forest_0.root");

	//TFile* inf1 = new TFile(inname1,"read");

	TFile* inf1 = TFile::Open(inname1,"read");

	cout << "point 1" << endl;
	
	float jtpt1[1000];		
	float jteta1[1000];
	float jtphi1[1000];
	float photonMax1[1000];
	float photonSum1[1000];
	float trackMax1[1000];
	float trackSum1[1000];
	float neutralMax1[1000];
	float neutralSum1[1000];
	float chargedMax1[1000];
	float vx1;
	float vy1;
	float vz1;
	float pthat1;
	int njets1;
	int bin1;
	int pcollisionEventSelection1;
	int pHBHENoiseFilter1;
	
	cout << "point 2" << endl;
	
	TTree* tree1 = (TTree*)inf1->Get(Form("%sJetAnalyzer/t", algorithm[algo]));
	TTree* evt1 = (TTree*)inf1->Get("hiEvtAnalyzer/HiTree");
	TTree* skim1= (TTree*)inf1->Get("skimanalysis/HltTree");
	TTree* hlt1 = (TTree*)inf1->Get("hltanalysis/HltTree");

	cout << "point 3" << endl;
	
	tree1->SetBranchAddress("jtpt",jtpt1);
	cout << "point 3.1" << endl;
	tree1->SetBranchAddress("jteta",jteta1);
	tree1->SetBranchAddress("jtphi",jtphi1);

	cout << "point 3.2" << endl;
	
	tree1->SetBranchAddress("photonMax",photonMax1);
	tree1->SetBranchAddress("photonSum",photonSum1);
	cout << "point 3.3" << endl;
	tree1->SetBranchAddress("trackMax",trackMax1);
	tree1->SetBranchAddress("trackSum",trackSum1);
	cout << "point 3.4" << endl;
	tree1->SetBranchAddress("neutralMax",neutralMax1);
	tree1->SetBranchAddress("neutralSum",neutralSum1);
	cout << "point 3.5" << endl;
	tree1->SetBranchAddress("chargedMax",chargedMax1);
	cout << "point 3.6" << endl;
	evt1->SetBranchAddress("vx",&vx1);
	evt1->SetBranchAddress("vy",&vy1);
	evt1->SetBranchAddress("vz",&vz1);
	cout << "point 3.7" << endl;
	tree1->SetBranchAddress("pthat",&pthat1);
	tree1->SetBranchAddress("nref",&njets1);
	evt1->SetBranchAddress("hiBin",&bin1);

	cout << "point 3.8" << endl;

	skim1->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection1);

	cout << "point 3.9" << endl;

	skim1->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter1);
	
	cout << "point 4" << endl;

//**
	Long64_t nentries1_ = tree1->GetEntries("pthat>30. && pthat<50.");
//**	
	Long64_t nentries1 = tree1->GetEntries();
	
	float scale1=xsection1/nentries1_;
	
	//event loop
	for (Long64_t jentry1=0; jentry1<nentries1;jentry1++)	{
		
  		tree1->GetEntry(jentry1);
  		evt1->GetEntry(jentry1);
		skim1->GetEntry(jentry1);
		hlt1->GetEntry(jentry1);
//**
		if(pthat1<30. || pthat1>50.) continue;
//**		
		if(pcollisionEventSelection1!=1 || pHBHENoiseFilter1!=1 || fabs(vz1)>15) continue; 
		
		double wcen1=1;
		double weight_vz1=1;
		int centb1=0;	

		if (reweigh)	{
			weight_vz1 = 1./fVz->Eval(vz1);
			wcen1=1./fcen->Eval(bin1);
		}	
		//! Centrality Bin Weighting Factor
		if(AAFlag){
			if(bin1<0 || bin1>36)continue;
			int centb1=-1;
				
			centb1=GetCentBin(bin1);
			wcen1 = fcen->Eval(bin1);
			weight_vz1 = fVz->Eval(vz1);
		
			if (centb1 == -1) continue;
		}
		
		hvx[centb1]->Fill(vx1);
		hvy[centb1]->Fill(vy1);
		hvz[centb1]->Fill(vz1, scale1*wcen1*weight_vz1);
		hbin[centb1]->Fill(bin1, scale1*wcen1*weight_vz1);
		hpthat[centb1]->Fill(pthat1,scale1*wcen1*weight_vz1);
		

		//jet loop
		
		for (int i= 0; i < njets1; i++)	{

			if (jtpt1[i]<global_jtpt_cut || jtpt1[0]<leading_jtpt_cut || abs(jteta1[i]) > jteta_cut) continue;
		
			//if((trackMax1[i]/jtpt1[i])<0.01) continue;
			//if((neutralMax1[i]/jtpt1[i])<0.01) continue;

			hpT[centb1]->Fill(jtpt1[i],scale1*wcen1*weight_vz1);
			hphi[centb1]->Fill(jtphi1[i],scale1*wcen1*weight_vz1);
			heta[centb1]->Fill(jteta1[i],scale1*wcen1*weight_vz1);
			
			hphotonMax[centb1]->Fill(photonMax1[i]/jtpt1[i],scale1*wcen1*weight_vz1);
			hphotonSum[centb1]->Fill(photonSum1[i]/jtpt1[i],scale1*wcen1*weight_vz1);
			htrackMax[centb1]->Fill(trackMax1[i]/jtpt1[i],scale1*wcen1*weight_vz1);
			htrackSum[centb1]->Fill(trackSum1[i]/jtpt1[i],scale1*wcen1*weight_vz1);
			hneutralMax[centb1]->Fill(neutralMax1[i]/jtpt1[i],scale1*wcen1*weight_vz1);
			hneutralSum[centb1]->Fill(neutralSum1[i]/jtpt1[i],scale1*wcen1*weight_vz1);
			
			hchargedMax[centb1]->Fill(chargedMax1[i]/jtpt1[i],scale1*wcen1*weight_vz1);

			heta_pT[centb1]->Fill(jteta1[i],jtpt1[i],scale1*wcen1*weight_vz1);
			heta_phi[centb1]->Fill(jteta1[i],jtphi1[i],scale1*wcen1*weight_vz1);
			
			
		}	
	}
	
	inf1->Close();
	

	//|\\|//|\\|//|\\|//|\\|//|\\|
	//		READ MC 50 FILE
	//|\\|//|\\|//|\\|//|\\|//|\\|
	
	cout<<"read MC 50 file"<<endl;
	
	TString inname2="/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root";
	//TFile *file2 = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt30/HiForest_v81_merged01/pt30_pp2013_P01_prod22_v81_merged_forest_0.root");

	TFile* inf2 = new TFile(inname2,"read");
	
	float jtpt2[1000];		
	float jteta2[1000];
	float jtphi2[1000];
	float photonMax2[1000];
	float photonSum2[1000];
	float trackMax2[1000];
	float trackSum2[1000];
	float neutralMax2[1000];
	float neutralSum2[1000];
	float chargedMax2[1000];
	float vx2;
	float vy2;
	float vz2;
	float pthat2;
	int njets2;
	int bin2;
	int pcollisionEventSelection2;
	int pHBHENoiseFilter2;
	
	
	TTree* tree2 = (TTree*)inf2->Get(Form("%sJetAnalyzer/t", algorithm[algo]));
	TTree* evt2 = (TTree*)inf2->Get("hiEvtAnalyzer/HiTree");
	TTree* skim2= (TTree*)inf2->Get("skimanalysis/HltTree");
	TTree* hlt2 = (TTree*)inf2->Get("hltanalysis/HltTree");
	
	tree2->SetBranchAddress("jtpt",jtpt2);
	tree2->SetBranchAddress("jteta",jteta2);
	tree2->SetBranchAddress("jtphi",jtphi2);
	
	tree2->SetBranchAddress("photonMax",photonMax2);
	tree2->SetBranchAddress("photonSum",photonSum2);
	tree2->SetBranchAddress("trackMax",trackMax2);
	tree2->SetBranchAddress("trackSum",trackSum2);
	tree2->SetBranchAddress("neutralMax",neutralMax2);
	tree2->SetBranchAddress("neutralSum",neutralSum2);
	tree2->SetBranchAddress("chargedMax",chargedMax2);
	evt2->SetBranchAddress("vx",&vx2);
	evt2->SetBranchAddress("vy",&vy2);
	evt2->SetBranchAddress("vz",&vz2);
	tree2->SetBranchAddress("pthat",&pthat2);
	tree2->SetBranchAddress("nref",&njets2);
	evt2->SetBranchAddress("hiBin",&bin2);
	skim2->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection2);
	skim2->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter2);
	

//**
	Long64_t nentries2_ = tree2->GetEntries("pthat>50. && pthat<80.");
//**	
	Long64_t nentries2 = tree2->GetEntries();
	
	float scale2=xsection2/nentries2_;
	
	//event loop
	for (Long64_t jentry2=0; jentry2<nentries2;jentry2++)	{
		
  		tree2->GetEntry(jentry2);
  		evt2->GetEntry(jentry2);
		skim2->GetEntry(jentry2);
		hlt2->GetEntry(jentry2);
//**
		if(pthat2<50. || pthat2>80.) continue;
//**		
		if(pcollisionEventSelection2!=1 || pHBHENoiseFilter2!=1 || fabs(vz2)>15) continue; 
		
		double wcen2=1;
		double weight_vz2=1;
		int centb2=0;
		if (reweigh)	{
			weight_vz2 = 1./fVz->Eval(vz2);
			wcen2=1./fcen->Eval(bin2);
		}	

		//! Centrality Bin Weighting Factor
		if(AAFlag){
			if(bin2<0 || bin2>36)continue;
			int centb2=-1;
				
			centb2=GetCentBin(bin2);
			wcen2 = fcen->Eval(bin2);
			weight_vz2 = fVz->Eval(vz2);
		
			if (centb2 == -1) continue;
		}
		
		hvx[centb2]->Fill(vx2);
		hvy[centb2]->Fill(vy2);
		hvz[centb2]->Fill(vz2, scale2*wcen2*weight_vz2);
		hbin[centb2]->Fill(bin2, scale2*wcen2*weight_vz2);
		hpthat[centb2]->Fill(pthat2,scale2*wcen2*weight_vz2);
		

		//jet loop
		
		for (int i= 0; i < njets2; i++)	{

			if (jtpt2[i]<global_jtpt_cut || jtpt2[0]<leading_jtpt_cut || abs(jteta2[i]) > jteta_cut) continue;
		
			//if((trackMax2[i]/jtpt2[i])<0.01) continue;
			//if((neutralMax2[i]/jtpt2[i])<0.01) continue;

			hpT[centb2]->Fill(jtpt2[i],scale2*wcen2*weight_vz2);
			hphi[centb2]->Fill(jtphi2[i],scale2*wcen2*weight_vz2);
			heta[centb2]->Fill(jteta2[i],scale2*wcen2*weight_vz2);
			
			hphotonMax[centb2]->Fill(photonMax2[i]/jtpt2[i],scale2*wcen2*weight_vz2);
			hphotonSum[centb2]->Fill(photonSum2[i]/jtpt2[i],scale2*wcen2*weight_vz2);
			htrackMax[centb2]->Fill(trackMax2[i]/jtpt2[i],scale2*wcen2*weight_vz2);
			htrackSum[centb2]->Fill(trackSum2[i]/jtpt2[i],scale2*wcen2*weight_vz2);
			hneutralMax[centb2]->Fill(neutralMax2[i]/jtpt2[i],scale2*wcen2*weight_vz2);
			hneutralSum[centb2]->Fill(neutralSum2[i]/jtpt2[i],scale2*wcen2*weight_vz2);
			
			hchargedMax[centb2]->Fill(chargedMax2[i]/jtpt2[i],scale2*wcen2*weight_vz2);

			heta_pT[centb2]->Fill(jteta2[i],jtpt2[i],scale2*wcen2*weight_vz2);
			heta_phi[centb2]->Fill(jteta2[i],jtphi2[i],scale2*wcen2*weight_vz2);
			
			
		}	
	}
	
	inf2->Close();

	//|\\|//|\\|//|\\|//|\\|//|\\|
	//		READ MC 80 FILE
	//|\\|//|\\|//|\\|//|\\|//|\\|
	
	cout<<"read MC 80 file"<<endl;
	
	TString inname3="/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root";
	//TFile *file3 = TFile::Open(" /mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt50/HiForest_v81_merged01/pt50_pp2013_P01_prod22_v81_merged_forest_0.root");

	TFile* inf3 = new TFile(inname3,"read");
	
	float jtpt3[1000];		
	float jteta3[1000];
	float jtphi3[1000];
	float photonMax3[1000];
	float photonSum3[1000];
	float trackMax3[1000];
	float trackSum3[1000];
	float neutralMax3[1000];
	float neutralSum3[1000];
	float chargedMax3[1000];
	float vx3;
	float vy3;
	float vz3;
	float pthat3;
	int njets3;
	int bin3;
	int pcollisionEventSelection3;
	int pHBHENoiseFilter3;
	
	
	TTree* tree3 = (TTree*)inf3->Get(Form("%sJetAnalyzer/t", algorithm[algo]));
	TTree* evt3 = (TTree*)inf3->Get("hiEvtAnalyzer/HiTree");
	TTree* skim3= (TTree*)inf3->Get("skimanalysis/HltTree");
	TTree* hlt3 = (TTree*)inf3->Get("hltanalysis/HltTree");
	
	tree3->SetBranchAddress("jtpt",jtpt3);
	tree3->SetBranchAddress("jteta",jteta3);
	tree3->SetBranchAddress("jtphi",jtphi3);
	
	tree3->SetBranchAddress("photonMax",photonMax3);
	tree3->SetBranchAddress("photonSum",photonSum3);
	tree3->SetBranchAddress("trackMax",trackMax3);
	tree3->SetBranchAddress("trackSum",trackSum3);
	tree3->SetBranchAddress("neutralMax",neutralMax3);
	tree3->SetBranchAddress("neutralSum",neutralSum3);
	tree3->SetBranchAddress("chargedMax",chargedMax3);
	evt3->SetBranchAddress("vx",&vx3);
	evt3->SetBranchAddress("vy",&vy3);
	evt3->SetBranchAddress("vz",&vz3);
	tree3->SetBranchAddress("pthat",&pthat3);
	tree3->SetBranchAddress("nref",&njets3);
	evt3->SetBranchAddress("hiBin",&bin3);
	skim3->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection3);
	skim3->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter3);
	

//**
	Long64_t nentries3_ = tree3->GetEntries("pthat>80. && pthat<100.");
//**	
	Long64_t nentries3 = tree3->GetEntries();
	
	float scale3=xsection3/nentries3_;
	
	//event loop
	for (Long64_t jentry3=0; jentry3<nentries3;jentry3++)	{
		
  		tree3->GetEntry(jentry3);
  		evt3->GetEntry(jentry3);
		skim3->GetEntry(jentry3);
		hlt3->GetEntry(jentry3);
//**
		if(pthat3<80. || pthat3>100.) continue;
//**		
		if(pcollisionEventSelection3!=1 || pHBHENoiseFilter3!=1 || fabs(vz3)>15) continue; 
		
		double wcen3=1;
		double weight_vz3=1;
		int centb3=0;	
		if (reweigh)	{
			weight_vz3 = 1./fVz->Eval(vz3);
			wcen3=1./fcen->Eval(bin3);
		}

		//! Centrality Bin Weighting Factor
		if(AAFlag){
			if(bin3<0 || bin3>36)continue;
			int centb3=-1;
				
			centb3=GetCentBin(bin3);
			wcen3 = fcen->Eval(bin3);
			weight_vz3 = fVz->Eval(vz3);
		
			if (centb3 == -1) continue;
		}
		
		hvx[centb3]->Fill(vx3);
		hvy[centb3]->Fill(vy3);
		hvz[centb3]->Fill(vz3, scale3*wcen3*weight_vz3);
		hbin[centb3]->Fill(bin3, scale3*wcen3*weight_vz3);
		hpthat[centb3]->Fill(pthat3,scale3*wcen3*weight_vz3);
		

		//jet loop
		
		for (int i= 0; i < njets3; i++)	{

			if (jtpt3[i]<global_jtpt_cut || jtpt3[0]<leading_jtpt_cut || abs(jteta3[i]) > jteta_cut) continue;
		
			//if((trackMax3[i]/jtpt3[i])<0.01) continue;
			//if((neutralMax3[i]/jtpt3[i])<0.01) continue;

			hpT[centb3]->Fill(jtpt3[i],scale3*wcen3*weight_vz3);
			hphi[centb3]->Fill(jtphi3[i],scale3*wcen3*weight_vz3);
			heta[centb3]->Fill(jteta3[i],scale3*wcen3*weight_vz3);
			
			hphotonMax[centb3]->Fill(photonMax3[i]/jtpt3[i],scale3*wcen3*weight_vz3);
			hphotonSum[centb3]->Fill(photonSum3[i]/jtpt3[i],scale3*wcen3*weight_vz3);
			htrackMax[centb3]->Fill(trackMax3[i]/jtpt3[i],scale3*wcen3*weight_vz3);
			htrackSum[centb3]->Fill(trackSum3[i]/jtpt3[i],scale3*wcen3*weight_vz3);
			hneutralMax[centb3]->Fill(neutralMax3[i]/jtpt3[i],scale3*wcen3*weight_vz3);
			hneutralSum[centb3]->Fill(neutralSum3[i]/jtpt3[i],scale3*wcen3*weight_vz3);
			
			hchargedMax[centb3]->Fill(chargedMax3[i]/jtpt3[i],scale3*wcen3*weight_vz3);

			heta_pT[centb3]->Fill(jteta3[i],jtpt3[i],scale3*wcen3*weight_vz3);
			heta_phi[centb3]->Fill(jteta3[i],jtphi3[i],scale3*wcen3*weight_vz3);
			
			
		}	
	}
	
	inf3->Close();

	//|\\|//|\\|//|\\|//|\\|//|\\|
	//		READ MC 100 FILE
	//|\\|//|\\|//|\\|//|\\|//|\\|
	
	cout<<"read MC 100 file"<<endl;
	
	TString inname4="/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHat_v0/0.root";
	//TFile *file4 = TFile::Open(" /mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt80/HiForest_v81_merged01/pt80_pp2013_P01_prod22_v81_merged_forest_0.root");

	TFile* inf4 = new TFile(inname4,"read");
	
	float jtpt4[1000];		
	float jteta4[1000];
	float jtphi4[1000];
	float photonMax4[1000];
	float photonSum4[1000];
	float trackMax4[1000];
	float trackSum4[1000];
	float neutralMax4[1000];
	float neutralSum4[1000];
	float chargedMax4[1000];
	float vx4;
	float vy4;
	float vz4;
	float pthat4;
	int njets4;
	int bin4;
	int pcollisionEventSelection4;
	int pHBHENoiseFilter4;
	
	
	TTree* tree4 = (TTree*)inf4->Get(Form("%sJetAnalyzer/t", algorithm[algo]));
	TTree* evt4 = (TTree*)inf4->Get("hiEvtAnalyzer/HiTree");
	TTree* skim4= (TTree*)inf4->Get("skimanalysis/HltTree");
	TTree* hlt4 = (TTree*)inf4->Get("hltanalysis/HltTree");
	
	tree4->SetBranchAddress("jtpt",jtpt4);
	tree4->SetBranchAddress("jteta",jteta4);
	tree4->SetBranchAddress("jtphi",jtphi4);
	
	tree4->SetBranchAddress("photonMax",photonMax4);
	tree4->SetBranchAddress("photonSum",photonSum4);
	tree4->SetBranchAddress("trackMax",trackMax4);
	tree4->SetBranchAddress("trackSum",trackSum4);
	tree4->SetBranchAddress("neutralMax",neutralMax4);
	tree4->SetBranchAddress("neutralSum",neutralSum4);
	tree4->SetBranchAddress("chargedMax",chargedMax4);
	evt4->SetBranchAddress("vx",&vx4);
	evt4->SetBranchAddress("vy",&vy4);
	evt4->SetBranchAddress("vz",&vz4);
	tree4->SetBranchAddress("pthat",&pthat4);
	tree4->SetBranchAddress("nref",&njets4);
	evt4->SetBranchAddress("hiBin",&bin4);
	skim4->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection4);
	skim4->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter4);
	

//**
	Long64_t nentries4_ = tree4->GetEntries("pthat>100. && pthat<120.");
//**	
	Long64_t nentries4 = tree4->GetEntries();
	
	float scale4=xsection4/nentries4_;
	
	//event loop
	for (Long64_t jentry4=0; jentry4<nentries4;jentry4++)	{
		
  		tree4->GetEntry(jentry4);
  		evt4->GetEntry(jentry4);
		skim4->GetEntry(jentry4);
		hlt4->GetEntry(jentry4);
//**
		if(pthat4<100. || pthat4>120.) continue;
//**		
		if(pcollisionEventSelection4!=1 || pHBHENoiseFilter4!=1 || fabs(vz4)>15) continue; 
		
		double wcen4=1;
		double weight_vz4=1;
		int centb4=0;	

		if (reweigh)	{
			weight_vz4 = 1./fVz->Eval(vz4);
			wcen4=1./fcen->Eval(bin4);
		}

		//! Centrality Bin Weighting Factor
		if(AAFlag){
			if(bin4<0 || bin4>36)continue;
			int centb4=-1;
				
			centb4=GetCentBin(bin4);
			wcen4 = fcen->Eval(bin4);
			weight_vz4 = fVz->Eval(vz4);
		
			if (centb4 == -1) continue;
		}
		
		hvx[centb4]->Fill(vx4);
		hvy[centb4]->Fill(vy4);
		hvz[centb4]->Fill(vz4, scale4*wcen4*weight_vz4);
		hbin[centb4]->Fill(bin4, scale4*wcen4*weight_vz4);
		hpthat[centb4]->Fill(pthat4,scale4*wcen4*weight_vz4);
		

		//jet loop
		
		for (int i= 0; i < njets4; i++)	{

			if (jtpt4[i]<global_jtpt_cut || jtpt4[0]<leading_jtpt_cut || abs(jteta4[i]) > jteta_cut) continue;
		
			//if((trackMax4[i]/jtpt4[i])<0.01) continue;
			//if((neutralMax4[i]/jtpt4[i])<0.01) continue;

			hpT[centb4]->Fill(jtpt4[i],scale4*wcen4*weight_vz4);
			hphi[centb4]->Fill(jtphi4[i],scale4*wcen4*weight_vz4);
			heta[centb4]->Fill(jteta4[i],scale4*wcen4*weight_vz4);
			
			hphotonMax[centb4]->Fill(photonMax4[i]/jtpt4[i],scale4*wcen4*weight_vz4);
			hphotonSum[centb4]->Fill(photonSum4[i]/jtpt4[i],scale4*wcen4*weight_vz4);
			htrackMax[centb4]->Fill(trackMax4[i]/jtpt4[i],scale4*wcen4*weight_vz4);
			htrackSum[centb4]->Fill(trackSum4[i]/jtpt4[i],scale4*wcen4*weight_vz4);
			hneutralMax[centb4]->Fill(neutralMax4[i]/jtpt4[i],scale4*wcen4*weight_vz4);
			hneutralSum[centb4]->Fill(neutralSum4[i]/jtpt4[i],scale4*wcen4*weight_vz4);
			
			hchargedMax[centb4]->Fill(chargedMax4[i]/jtpt4[i],scale4*wcen4*weight_vz4);

			heta_pT[centb4]->Fill(jteta4[i],jtpt4[i],scale4*wcen4*weight_vz4);
			heta_phi[centb4]->Fill(jteta4[i],jtphi4[i],scale4*wcen4*weight_vz4);
			
			
		}	
	}
	
	inf4->Close();

	//|\\|//|\\|//|\\|//|\\|//|\\|
	//		READ MC 120 FILE
	//|\\|//|\\|//|\\|//|\\|//|\\|
	
	cout<<"read MC 120 file"<<endl;
	
	TString inname5=" /mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0/0.root";
	//TFile *file5 = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt120/HiForest_v81_merged01/pt120_pp2013_P01_prod22_v81_merged_forest_0.root");

	
	TFile* inf5 = new TFile(inname5,"read");
	
	float jtpt5[1000];		
	float jteta5[1000];
	float jtphi5[1000];
	float photonMax5[1000];
	float photonSum5[1000];
	float trackMax5[1000];
	float trackSum5[1000];
	float neutralMax5[1000];
	float neutralSum5[1000];
	float chargedMax5[1000];
	float vx5;
	float vy5;
	float vz5;
	float pthat5;
	int njets5;
	int bin5;
	int pcollisionEventSelection5;
	int pHBHENoiseFilter5;
	
	
	TTree* tree5 = (TTree*)inf5->Get(Form("%sJetAnalyzer/t", algorithm[algo]));
	TTree* evt5 = (TTree*)inf5->Get("hiEvtAnalyzer/HiTree");
	TTree* skim5= (TTree*)inf5->Get("skimanalysis/HltTree");
	TTree* hlt5 = (TTree*)inf5->Get("hltanalysis/HltTree");
	
	tree5->SetBranchAddress("jtpt",jtpt5);
	tree5->SetBranchAddress("jteta",jteta5);
	tree5->SetBranchAddress("jtphi",jtphi5);
	
	tree5->SetBranchAddress("photonMax",photonMax5);
	tree5->SetBranchAddress("photonSum",photonSum5);
	tree5->SetBranchAddress("trackMax",trackMax5);
	tree5->SetBranchAddress("trackSum",trackSum5);
	tree5->SetBranchAddress("neutralMax",neutralMax5);
	tree5->SetBranchAddress("neutralSum",neutralSum5);
	tree5->SetBranchAddress("chargedMax",chargedMax5);
	evt5->SetBranchAddress("vx",&vx5);
	evt5->SetBranchAddress("vy",&vy5);
	evt5->SetBranchAddress("vz",&vz5);
	tree5->SetBranchAddress("pthat",&pthat5);
	tree5->SetBranchAddress("nref",&njets5);
	evt5->SetBranchAddress("hiBin",&bin5);
	skim5->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection5);
	skim5->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter5);
	

//**
	Long64_t nentries5_ = tree5->GetEntries("pthat>120. && pthat<1000.");
//**	
	Long64_t nentries5 = tree5->GetEntries();
	
	float scale5=xsection5/nentries5_;
	
	//event loop
	for (Long64_t jentry5=0; jentry5<nentries5;jentry5++)	{
		
  		tree5->GetEntry(jentry5);
  		evt5->GetEntry(jentry5);
		skim5->GetEntry(jentry5);
		hlt5->GetEntry(jentry5);
//**
		if(pthat5<120. || pthat5>1000.) continue;
//**		
		if(pcollisionEventSelection5!=1 || pHBHENoiseFilter5!=1 || fabs(vz5)>15) continue; 
		
		double wcen5=1;
		double weight_vz5=1;
		int centb5=0;	

		if (reweigh)	{
			weight_vz5 = 1./fVz->Eval(vz5);
			wcen5=1./fcen->Eval(bin5);
		}

		//! Centrality Bin Weighting Factor
		if(AAFlag){
			if(bin5<0 || bin5>36)continue;
			int centb5=-1;
				
			centb5=GetCentBin(bin5);
			wcen5 = fcen->Eval(bin5);
			weight_vz5 = fVz->Eval(vz5);
		
			if (centb5 == -1) continue;
		}
		
		hvx[centb5]->Fill(vx5);
		hvy[centb5]->Fill(vy5);
		hvz[centb5]->Fill(vz5, scale5*wcen5*weight_vz5);
		hbin[centb5]->Fill(bin5, scale5*wcen5*weight_vz5);
		hpthat[centb5]->Fill(pthat5,scale5*wcen5*weight_vz5);
		

		//jet loop
		
		for (int i= 0; i < njets5; i++)	{

			if (jtpt5[i]<global_jtpt_cut || jtpt5[0]<leading_jtpt_cut || abs(jteta5[i]) > jteta_cut) continue;
		
			//if((trackMax5[i]/jtpt5[i])<0.01) continue;
			//if((neutralMax5[i]/jtpt5[i])<0.01) continue;

			hpT[centb5]->Fill(jtpt5[i],scale5*wcen5*weight_vz5);
			hphi[centb5]->Fill(jtphi5[i],scale5*wcen5*weight_vz5);
			heta[centb5]->Fill(jteta5[i],scale5*wcen5*weight_vz5);
			
			hphotonMax[centb5]->Fill(photonMax5[i]/jtpt5[i],scale5*wcen5*weight_vz5);
			hphotonSum[centb5]->Fill(photonSum5[i]/jtpt5[i],scale5*wcen5*weight_vz5);
			htrackMax[centb5]->Fill(trackMax5[i]/jtpt5[i],scale5*wcen5*weight_vz5);
			htrackSum[centb5]->Fill(trackSum5[i]/jtpt5[i],scale5*wcen5*weight_vz5);
			hneutralMax[centb5]->Fill(neutralMax5[i]/jtpt5[i],scale5*wcen5*weight_vz5);
			hneutralSum[centb5]->Fill(neutralSum5[i]/jtpt5[i],scale5*wcen5*weight_vz5);
			
			hchargedMax[centb5]->Fill(chargedMax5[i]/jtpt5[i],scale5*wcen5*weight_vz5);

			heta_pT[centb5]->Fill(jteta5[i],jtpt5[i],scale5*wcen5*weight_vz5);
			heta_phi[centb5]->Fill(jteta5[i],jtphi5[i],scale5*wcen5*weight_vz5);
			
			
		}	
	}
	
	inf5->Close();



	//|\\|//|\\|//|\\|//|\\|//|\\|
	//		READ MC 170 FILE
	//|\\|//|\\|//|\\|//|\\|//|\\|
	/*
	cout<<"read MC 170 file"<<endl;
	
	TString inname6="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt170/HiForest_v81_merged01/pt170_pp2013_P01_prod22_v81_merged_forest_0.root";
	//TFile *file6 = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt170/HiForest_v81_merged01/pt170_pp2013_P01_prod22_v81_merged_forest_0.root");

	TFile* inf6 = new TFile(inname6,"read");
	
	float jtpt6[1000];		
	float jteta6[1000];
	float jtphi6[1000];
	float photonMax6[1000];
	float photonSum6[1000];
	float trackMax6[1000];
	float trackSum6[1000];
	float neutralMax6[1000];
	float neutralSum6[1000];
	float chargedMax6[1000];
	float vx6;
	float vy6;
	float vz6;
	float pthat6;
	int njets6;
	int bin6;
	int pPAcollisionEventSelectionPA6;
	int pHBHENoiseFilter6;
	
	
	TTree* tree6 = (TTree*)inf6->Get(Form("%sJetAnalyzer/t", algorithm[algo]));
	TTree* evt6 = (TTree*)inf6->Get("hiEvtAnalyzer/HiTree");
	TTree* skim6= (TTree*)inf6->Get("skimanalysis/HltTree");
	TTree* hlt6 = (TTree*)inf6->Get("hltanalysis/HltTree");
	
	tree6->SetBranchAddress("jtpt",jtpt6);
	tree6->SetBranchAddress("jteta",jteta6);
	tree6->SetBranchAddress("jtphi",jtphi6);
	
	tree6->SetBranchAddress("photonMax",photonMax6);
	tree6->SetBranchAddress("photonSum",photonSum6);
	tree6->SetBranchAddress("trackMax",trackMax6);
	tree6->SetBranchAddress("trackSum",trackSum6);
	tree6->SetBranchAddress("neutralMax",neutralMax6);
	tree6->SetBranchAddress("neutralSum",neutralSum6);
	tree6->SetBranchAddress("chargedMax",chargedMax6);
	evt6->SetBranchAddress("vx",&vx6);
	evt6->SetBranchAddress("vy",&vy6);
	evt6->SetBranchAddress("vz",&vz6);
	tree6->SetBranchAddress("pthat",&pthat6);
	tree6->SetBranchAddress("nref",&njets6);
	evt6->SetBranchAddress("hiBin",&bin6);
	skim6->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA6);
	skim6->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter6);
	

//**
	Long64_t nentries6_ = tree6->GetEntries("pthat>170. && pthat<220.");
//**	
	Long64_t nentries6 = tree6->GetEntries();
	
	float scale6=xsection6/nentries6_;
	
	//event loop
	for (Long64_t jentry6=0; jentry6<nentries6;jentry6++)	{
		
  		tree6->GetEntry(jentry6);
  		evt6->GetEntry(jentry6);
		skim6->GetEntry(jentry6);
		hlt6->GetEntry(jentry6);
//**
		if(pthat6<170. || pthat6>220.) continue;
//**		
		if(pPAcollisionEventSelectionPA6!=1 || pHBHENoiseFilter6!=1 || fabs(vz6)>15) continue; 
		
		double wcen6=1;
		double weight_vz6=1;
		int centb6=0;	

		//weight_vz6 = 1./fVz->Eval(vz6);

		//! Centrality Bin Weighting Factor
		if(AAFlag){
			if(bin6<0 || bin6>36)continue;
			int centb6=-1;
				
			centb6=GetCentBin(bin6);
			wcen6 = fcen->Eval(bin6);
			weight_vz6 = fVz->Eval(vz6);
		
			if (centb6 == -1) continue;
		}
		
		hvx[centb6]->Fill(vx6);
		hvy[centb6]->Fill(vy6);
		hvz[centb6]->Fill(vz6, scale6*wcen6*weight_vz6);
		hpthat[centb6]->Fill(pthat6,scale6*wcen6*weight_vz6);
		

		//jet loop
		
		for (int i= 0; i < njets6; i++)	{

			if (jtpt6[i]<global_jtpt_cut || jtpt6[0]<leading_jtpt_cut || abs(jteta6[i]) > jteta_cut) continue;
		
			//if((trackMax6[i]/jtpt6[i])<0.01) continue;
			//if((neutralMax6[i]/jtpt6[i])<0.01) continue;

			hpT[centb6]->Fill(jtpt6[i],scale6*wcen6*weight_vz6);
			hphi[centb6]->Fill(jtphi6[i],scale6*wcen6*weight_vz6);
			heta[centb6]->Fill(jteta6[i],scale6*wcen6*weight_vz6);
			
			hphotonMax[centb6]->Fill(photonMax6[i]/jtpt6[i],scale6*wcen6*weight_vz6);
			hphotonSum[centb6]->Fill(photonSum6[i]/jtpt6[i],scale6*wcen6*weight_vz6);
			htrackMax[centb6]->Fill(trackMax6[i]/jtpt6[i],scale6*wcen6*weight_vz6);
			htrackSum[centb6]->Fill(trackSum6[i]/jtpt6[i],scale6*wcen6*weight_vz6);
			hneutralMax[centb6]->Fill(neutralMax6[i]/jtpt6[i],scale6*wcen6*weight_vz6);
			hneutralSum[centb6]->Fill(neutralSum6[i]/jtpt6[i],scale6*wcen6*weight_vz6);
			
			hchargedMax[centb6]->Fill(chargedMax6[i]/jtpt6[i],scale6*wcen6*weight_vz6);

			heta_pT[centb6]->Fill(jteta6[i],jtpt6[i],scale6*wcen6*weight_vz6);
			heta_phi[centb6]->Fill(jteta6[i],jtphi6[i],scale6*wcen6*weight_vz6);
			
			
		}	
	}
	
	inf6->Close();

	//|\\|//|\\|//|\\|//|\\|//|\\|
	//		READ MC 220 FILE
	//|\\|//|\\|//|\\|//|\\|//|\\|
	
	cout<<"read MC 220 file"<<endl;
	
	TString inname7="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt220/HiForest_v81_merged01/pt220_pp2013_P01_prod22_v81_merged_forest_0.root";
	//TFile *file7 = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt220/HiForest_v81_merged01/pt220_pp2013_P01_prod22_v81_merged_forest_0.root");

	TFile* inf7 = new TFile(inname7,"read");
	
	float jtpt7[1000];		
	float jteta7[1000];
	float jtphi7[1000];
	float photonMax7[1000];
	float photonSum7[1000];
	float trackMax7[1000];
	float trackSum7[1000];
	float neutralMax7[1000];
	float neutralSum7[1000];
	float chargedMax7[1000];
	float vx7;
	float vy7;
	float vz7;
	float pthat7;
	int njets7;
	int bin7;
	int pPAcollisionEventSelectionPA7;
	int pHBHENoiseFilter7;
	
	
	TTree* tree7 = (TTree*)inf7->Get(Form("%sJetAnalyzer/t", algorithm[algo]));
	TTree* evt7 = (TTree*)inf7->Get("hiEvtAnalyzer/HiTree");
	TTree* skim7= (TTree*)inf7->Get("skimanalysis/HltTree");
	TTree* hlt7 = (TTree*)inf7->Get("hltanalysis/HltTree");
	
	tree7->SetBranchAddress("jtpt",jtpt7);
	tree7->SetBranchAddress("jteta",jteta7);
	tree7->SetBranchAddress("jtphi",jtphi7);
	
	tree7->SetBranchAddress("photonMax",photonMax7);
	tree7->SetBranchAddress("photonSum",photonSum7);
	tree7->SetBranchAddress("trackMax",trackMax7);
	tree7->SetBranchAddress("trackSum",trackSum7);
	tree7->SetBranchAddress("neutralMax",neutralMax7);
	tree7->SetBranchAddress("neutralSum",neutralSum7);
	tree7->SetBranchAddress("chargedMax",chargedMax7);
	evt7->SetBranchAddress("vx",&vx7);
	evt7->SetBranchAddress("vy",&vy7);
	evt7->SetBranchAddress("vz",&vz7);
	tree7->SetBranchAddress("pthat",&pthat7);
	tree7->SetBranchAddress("nref",&njets7);
	evt7->SetBranchAddress("hiBin",&bin7);
	skim7->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA7);
	skim7->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter7);
	

//**
	Long64_t nentries7_ = tree7->GetEntries("pthat>220. && pthat<280.");
//**	
	Long64_t nentries7 = tree7->GetEntries();
	
	float scale7=xsection7/nentries7_;
	
	//event loop
	for (Long64_t jentry7=0; jentry7<nentries7;jentry7++)	{
		
  		tree7->GetEntry(jentry7);
  		evt7->GetEntry(jentry7);
		skim7->GetEntry(jentry7);
		hlt7->GetEntry(jentry7);
//**
		if(pthat7<220. || pthat7>280.) continue;
//**		
		if(pPAcollisionEventSelectionPA7!=1 || pHBHENoiseFilter7!=1 || fabs(vz7)>15) continue; 
		
		double wcen7=1;
		double weight_vz7=1;
		int centb7=0;	

		weight_vz7 = 1./fVz->Eval(vz7);

		//! Centrality Bin Weighting Factor
		if(AAFlag){
		if(bin7<0 || bin7>36)continue;
		int centb7=-1;
				
		centb7=GetCentBin(bin7);
		wcen7 = fcen->Eval(bin7);
		//weight_vz7 = fVz->Eval(vz7);
		
		if (centb7 == -1) continue;
		}
		
		hvx[centb7]->Fill(vx7);
		hvy[centb7]->Fill(vy7);
		hvz[centb7]->Fill(vz7, scale7*wcen7*weight_vz7);
		hpthat[centb7]->Fill(pthat7,scale7*wcen7*weight_vz7);
		

		//jet loop
		
		for (int i= 0; i < njets7; i++)	{

			if (jtpt7[i]<global_jtpt_cut || jtpt7[0]<leading_jtpt_cut || abs(jteta7[i]) > jteta_cut) continue;
		
			//if((trackMax7[i]/jtpt7[i])<0.01) continue;
			//if((neutralMax7[i]/jtpt7[i])<0.01) continue;

			hpT[centb7]->Fill(jtpt7[i],scale7*wcen7*weight_vz7);
			hphi[centb7]->Fill(jtphi7[i],scale7*wcen7*weight_vz7);
			heta[centb7]->Fill(jteta7[i],scale7*wcen7*weight_vz7);
			
			hphotonMax[centb7]->Fill(photonMax7[i]/jtpt7[i],scale7*wcen7*weight_vz7);
			hphotonSum[centb7]->Fill(photonSum7[i]/jtpt7[i],scale7*wcen7*weight_vz7);
			htrackMax[centb7]->Fill(trackMax7[i]/jtpt7[i],scale7*wcen7*weight_vz7);
			htrackSum[centb7]->Fill(trackSum7[i]/jtpt7[i],scale7*wcen7*weight_vz7);
			hneutralMax[centb7]->Fill(neutralMax7[i]/jtpt7[i],scale7*wcen7*weight_vz7);
			hneutralSum[centb7]->Fill(neutralSum7[i]/jtpt7[i],scale7*wcen7*weight_vz7);
			
			hchargedMax[centb7]->Fill(chargedMax7[i]/jtpt7[i],scale7*wcen7*weight_vz7);

			heta_pT[centb7]->Fill(jteta7[i],jtpt7[i],scale7*wcen7*weight_vz7);
			heta_phi[centb7]->Fill(jteta7[i],jtphi7[i],scale7*wcen7*weight_vz7);
			
			
		}	
	}
	
	inf7->Close();

	//|\\|//|\\|//|\\|//|\\|//|\\|
	//		READ MC 280 FILE
	//|\\|//|\\|//|\\|//|\\|//|\\|
	
	cout<<"read MC 280 file"<<endl;
	
	TString inname8="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt280/HiForest_v81_merged01/pt280_pp2013_P01_prod22_v81_merged_forest_0.root";
	//TFile *file8 = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt280/HiForest_v81_merged01/pt280_pp2013_P01_prod22_v81_merged_forest_0.root");

	TFile* inf8 = new TFile(inname8,"read");
	
	float jtpt8[1000];		
	float jteta8[1000];
	float jtphi8[1000];
	float photonMax8[1000];
	float photonSum8[1000];
	float trackMax8[1000];
	float trackSum8[1000];
	float neutralMax8[1000];
	float neutralSum8[1000];
	float chargedMax8[1000];
	float vx8;
	float vy8;
	float vz8;
	float pthat8;
	int njets8;
	int bin8;
	int pPAcollisionEventSelectionPA8;
	int pHBHENoiseFilter8;
	
	
	TTree* tree8 = (TTree*)inf8->Get(Form("%sJetAnalyzer/t", algorithm[algo]));
	TTree* evt8 = (TTree*)inf8->Get("hiEvtAnalyzer/HiTree");
	TTree* skim8= (TTree*)inf8->Get("skimanalysis/HltTree");
	TTree* hlt8 = (TTree*)inf8->Get("hltanalysis/HltTree");
	
	tree8->SetBranchAddress("jtpt",jtpt8);
	tree8->SetBranchAddress("jteta",jteta8);
	tree8->SetBranchAddress("jtphi",jtphi8);
	
	tree8->SetBranchAddress("photonMax",photonMax8);
	tree8->SetBranchAddress("photonSum",photonSum8);
	tree8->SetBranchAddress("trackMax",trackMax8);
	tree8->SetBranchAddress("trackSum",trackSum8);
	tree8->SetBranchAddress("neutralMax",neutralMax8);
	tree8->SetBranchAddress("neutralSum",neutralSum8);
	tree8->SetBranchAddress("chargedMax",chargedMax8);
	evt8->SetBranchAddress("vx",&vx8);
	evt8->SetBranchAddress("vy",&vy8);
	evt8->SetBranchAddress("vz",&vz8);
	tree8->SetBranchAddress("pthat",&pthat8);
	tree8->SetBranchAddress("nref",&njets8);
	evt8->SetBranchAddress("hiBin",&bin8);
	skim8->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA8);
	skim8->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter8);
	

//**
	Long64_t nentries8_ = tree8->GetEntries("pthat>280. && pthat<1000.");
//**	
	Long64_t nentries8 = tree8->GetEntries();
	
	float scale8=xsection8/nentries8_;
	
	//event loop
	for (Long64_t jentry8=0; jentry8<nentries8;jentry8++)	{
		
  		tree8->GetEntry(jentry8);
  		evt8->GetEntry(jentry8);
		skim8->GetEntry(jentry8);
		hlt8->GetEntry(jentry8);
//**
		if(pthat8<280. || pthat8>1000.) continue;
//**		
		if(pPAcollisionEventSelectionPA8!=1 || pHBHENoiseFilter8!=1 || fabs(vz8)>15) continue; 
		
		double wcen8=1;
		double weight_vz8=1;
		int centb8=0;	

		weight_vz8 = 1./fVz->Eval(vz8);

		//! Centrality Bin Weighting Factor
		if(AAFlag){
		if(bin8<0 || bin8>36)continue;
		int centb8=-1;
				
		centb8=GetCentBin(bin8);
		wcen8 = fcen->Eval(bin8);
		//weight_vz8 = fVz->Eval(vz8);
		
		if (centb8 == -1) continue;
		}
		
		hvx[centb8]->Fill(vx8);
		hvy[centb8]->Fill(vy8);
		hvz[centb8]->Fill(vz8,scale8*wcen8*weight_vz8);
		hpthat[centb8]->Fill(pthat8,scale8*wcen8*weight_vz8);
		

		//jet loop
		
		for (int i= 0; i < njets8; i++)	{

			if (jtpt8[i]<global_jtpt_cut || jtpt8[0]<leading_jtpt_cut || abs(jteta8[i]) > jteta_cut) continue;
		
			//if((trackMax8[i]/jtpt8[i])<0.01) continue;
			//if((neutralMax8[i]/jtpt8[i])<0.01) continue;

			hpT[centb8]->Fill(jtpt8[i],scale8*wcen8*weight_vz8);
			hphi[centb8]->Fill(jtphi8[i],scale8*wcen8*weight_vz8);
			heta[centb8]->Fill(jteta8[i],scale8*wcen8*weight_vz8);
			
			hphotonMax[centb8]->Fill(photonMax8[i]/jtpt8[i],scale8*wcen8*weight_vz8);
			hphotonSum[centb8]->Fill(photonSum8[i]/jtpt8[i],scale8*wcen8*weight_vz8);
			htrackMax[centb8]->Fill(trackMax8[i]/jtpt8[i],scale8*wcen8*weight_vz8);
			htrackSum[centb8]->Fill(trackSum8[i]/jtpt8[i],scale8*wcen8*weight_vz8);
			hneutralMax[centb8]->Fill(neutralMax8[i]/jtpt8[i],scale8*wcen8*weight_vz8);
			hneutralSum[centb8]->Fill(neutralSum8[i]/jtpt8[i],scale8*wcen8*weight_vz8);
			
			hchargedMax[centb8]->Fill(chargedMax8[i]/jtpt8[i],scale8*wcen8*weight_vz8);

			heta_pT[centb8]->Fill(jteta8[i],jtpt8[i],scale8*wcen8*weight_vz8);
			heta_phi[centb8]->Fill(jteta8[i],jtphi8[i],scale8*wcen8*weight_vz8);
			
			
		}	
	}
	
	inf8->Close();	*/

	
	
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//
	//		NORMALIZATION
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//
	
	for(int icen=0;icen<ncen;icen++)	{

		hpT[icen]->Scale(1/hpT[icen]->Integral());
		heta[icen]->Scale(1/heta[icen]->Integral());
		hphi[icen]->Scale(1/hphi[icen]->Integral());
		
		hphotonMax[icen]->Scale(1/hphotonMax[icen]->Integral());
		hphotonSum[icen]->Scale(1/hphotonSum[icen]->Integral());
		htrackMax[icen]->Scale(1/htrackMax[icen]->Integral());
		htrackSum[icen]->Scale(1/htrackSum[icen]->Integral());
		hneutralMax[icen]->Scale(1/hneutralMax[icen]->Integral());
		hneutralSum[icen]->Scale(1/hneutralSum[icen]->Integral());

		hvx[icen]->Scale(1/hvx[icen]->Integral());
		hvy[icen]->Scale(1/hvy[icen]->Integral());
		hvz[icen]->Scale(1/hvz[icen]->Integral());
		hbin[icen]->Scale(1/hbin[icen]->Integral());;
		/*htrackMax[icen]->GetYaxis()->SetTitleOffset(1.5);
		htrackSum[icen]->GetYaxis()->SetTitleOffset(1.5);
		hneutralMax[icen]->GetYaxis()->SetTitleOffset(1.5);
		hneutralSum[icen]->GetYaxis()->SetTitleOffset(1.5);
		hphotonMax[icen]->GetYaxis()->SetTitleOffset(1.5);
		hphotonSum[icen]->GetYaxis()->SetTitleOffset(1.5);*/

		//heta_pT[icen]->Scale(1/heta_pT[icen]->Integral());
		//heta_phi[icen]->Scale(1/heta_phi[icen]->Integral());
	}
	
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//
	//		CREATE CANVASSES AND PLOT HISTOGRAMS
	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//

	gStyle->SetHistFillStyle(3001);
	gStyle->SetHistFillColor(38);
	gStyle->SetHistLineColor(38);
	gStyle->SetOptStat(1111);

	///\///\///\///\///\///\///\///\///\///\///\///\///\///\///\///\
	//
	//	produce SINGLE plots --- no canvas splitting
	//
	///\///\///\///\///\///\///\///\///\///\///\///\///\///\///\///\

	//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|//|\\|
	//~~~//~~~//~~~//~~~//~~~//~~~//~~~//~~~//~~~//~~~//~~~//~~~
	//***
	//standard centralities 0-5,5-10,10-30, 30-50, 50-70, 70-90
	//new centralities 0-10, 10-30, 30-50, 50-100

	TCanvas *c1s = new TCanvas("c1s","c1s",500,500);
	c1s->cd();
	gPad->SetLogy();
	hpT[0]->UseCurrentStyle();
	hpT[0]->Draw("hist");
	hpT[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c2s = new TCanvas("c2s","c2s",500,500);
	c2s->cd();
	heta[0]->UseCurrentStyle();
	heta[0]->SetMinimum(0);
	heta[0]->Draw("hist");
	heta[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c3s = new TCanvas("c3s","c3s",500,500);
	c3s->cd();
	hphi[0]->UseCurrentStyle();
	hphi[0]->SetMinimum(0);
	hphi[0]->Draw("hist");
	hphi[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c4s = new TCanvas("c4s","c4s",500,500);
	c4s->cd();
	hvz[0]->UseCurrentStyle();
	hvz[0]->Draw("hist");
	hvz[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c5s = new TCanvas("c5s","c5s",500,500);
	c5s->cd();
	hbin[0]->UseCurrentStyle();
	hbin[0]->Draw("hist");
	hbin[0]->GetYaxis()->SetTitleOffset(1.3);


	/*TCanvas *c5s = new TCanvas("c5s","c5s",500,500);
	c5s->cd();
	hphotonMax[0]->UseCurrentStyle();
	hphotonMax[0]->Draw("hist");
	hphotonMax[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c6s = new TCanvas("c6s","c6s",500,500);
	c6s->cd();
	hphotonSum[0]->UseCurrentStyle();
	hphotonSum[0]->Draw("hist");
	hphotonSum[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c7s = new TCanvas("c7s","c7s",500,500);
	c7s->cd();
	htrackMax[0]->UseCurrentStyle();
	htrackMax[0]->Draw("hist");
	htrackMax[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c8s = new TCanvas("c8s","c8s",500,500);
	c8s->cd();
	htrackSum[0]->UseCurrentStyle();
	htrackSum[0]->Draw("hist");
	htrackSum[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c9s = new TCanvas("c9s","c9s",500,500);
	c9s->cd();
	hneutralMax[0]->UseCurrentStyle();
	hneutralMax[0]->Draw("hist");
	hneutralMax[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c10s = new TCanvas("c10s","c10s",500,500);
	c10s->cd();
	hneutralSum[0]->UseCurrentStyle();
	hneutralSum[0]->Draw("hist");
	hneutralSum[0]->GetYaxis()->SetTitleOffset(1.3);

	TCanvas *c11s = new TCanvas("c11s","c11s",500,500);
	c11s->cd();
	hchargedMax[0]->UseCurrentStyle();
	hchargedMax[0]->Draw("hist");
	hchargedMax[0]->GetYaxis()->SetTitleOffset(1.3);

	/*TCanvas *c12s = new TCanvas("c12s","c12s",500,500);
	c12s->cd();
	hvz[0]->UseCurrentStyle();
	hvz[0]->Draw("hist");
	hvz[0]->GetYaxis()->SetTitleOffset(1.3);*/

	/*c1s->Print("JetIDv6singleFinalakPu3PFv25.pdf[");
	c1s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	//c2->Print("TransverseMomentumLeading.pdf");
	c2s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c3s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c4s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c5s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c6s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c7s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c8s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c9s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c10s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c11s->Print("JetIDv6singleFinalakPu3PFv25.pdf");
	c11s->Print("JetIDv6singleFinalakPu3PFv25.pdf]");

	printf("Output file closed.\n"); */
	
	TFile *outputFile = new TFile(Form("%s_%sREWEIGHTED.root", macroname, algorithm[algo]), "recreate");

	//TFile *outputFile = new TFile(Form("%s_%s.root", macroname, algorithm[algo]), "recreate");
	hpT[0]->Write();
	heta[0]->Write();
	hphi[0]->Write();
	hvz[0]->Write();
	hbin[0]->Write();
	//hphotonMax[0]->Write();
	//hphotonSum[0]->Write();
	//htrackMax[0]->Write();
	//htrackSum[0]->Write();
	//hneutralMax[0]->Write();
	//hneutralSum[0]->Write();
	//hchargedMax[0]->Write();
	outputFile->Close();
}


int GetCentBin(int bin)
{
	int ibin=-1;
	//! centrality is defined as 2.5% bins of cross section
	//! in 0-39 bins
	
    if(bin<2)ibin=0; //! 0-5%
    else if(bin>=2 && bin<4)ibin=1;    //! 5-10%
    else if(bin>=4 && bin<12)ibin=2;   //! 10-30%
    else if(bin>=12 && bin<20)ibin=3;  //! 30-50%
    else if(bin>=20 && bin<28)ibin=4;  //! 50-70%
    else if(bin>=28 && bin<36)ibin=5;  //! 70-90%
    return ibin;

}

void drawCentTex(const int centbin){
 TLatex *tex=0;
 
if(centbin==0)
  tex = new TLatex(0.7,0.8,"0-5%");
if(centbin==1)
  tex = new TLatex(0.7,0.8,"5-10%");
if(centbin==2)
  tex = new TLatex(0.7,0.8,"10-30%");
if(centbin==3)
  tex = new TLatex(0.7,0.8,"30-50%");
if(centbin==4)
  tex = new TLatex(0.7,0.8,"50-70%");
if(centbin==5)
  tex = new TLatex(0.7,0.8,"70-90%");


  tex->SetTextFont(63);
  tex->SetTextSize(22);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}



void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}




