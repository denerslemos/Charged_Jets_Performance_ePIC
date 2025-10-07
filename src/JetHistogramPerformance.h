TH1D *NEvents = new TH1D("NEvents","",2,0.,2.); NEvents->Sumw2();
TH1D *numRecoJetsEventHist = new TH1D("numRecoJetsEvent","",20,0.,20.); numRecoJetsEventHist->Sumw2();
TH1D *numRecoJetsNoElecEventHist = new TH1D("numRecoJetsNoElecEvent","",20,0.,20.); numRecoJetsNoElecEventHist->Sumw2();
TH1D *numGenJetsEventHist = new TH1D("numGenJetsEventHist","",20,0.,20.); numGenJetsEventHist->Sumw2();
TH1D *numGenJetsNoElecEventHist = new TH1D("numGenJetsNoElecEventHist","",20,0.,20.); numGenJetsNoElecEventHist->Sumw2();

TH1D *JetdR = new TH1D("JetdR","",100,0.,10.); JetdR->Sumw2();
TH1D *JetMindR = new TH1D("JetMindR","",100,0.,10.); JetMindR->Sumw2();

// JER/JES -> {(reco-match)/match, variable: mostly pT }
const int NJESJERAxis = 2;
int	JESJERBins[NJESJERAxis]      =   { 1000  , 2000  };
double JESJERXmin[NJESJERAxis]   =   { -1.0  , -5.0  };
double JESJERXmax[NJESJERAxis]   =   { 1.0  , 195.0 };
THnSparseD *mHistJESJERvsPt = new THnSparseD("mHistJESJERvsPt", "mHistJESJERvsPt", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsPt->Sumw2();
THnSparseD *mHistJESJERvsPtNoElec = new THnSparseD("mHistJESJERvsPtNoElec", "mHistJESJERvsPtNoElec", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsPtNoElec->Sumw2();
THnSparseD *mHistJESJERvsE = new THnSparseD("mHistJESJERvsE", "mHistJESJERvsE", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsE->Sumw2();
THnSparseD *mHistJESJERvsENoElec = new THnSparseD("mHistJESJERvsENoElec", "mHistJESJERvsENoElec", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsENoElec->Sumw2();
THnSparseD *mHistJESJERvsM = new THnSparseD("mHistJESJERvsM", "mHistJESJERvsM", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsM->Sumw2();
THnSparseD *mHistJESJERvsMNoElec = new THnSparseD("mHistJESJERvsMNoElec", "mHistJESJERvsMNoElec", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsMNoElec->Sumw2();
THnSparseD *mHistJESJERvsPhi = new THnSparseD("mHistJESJERvsPhi", "mHistJESJERvsPhi", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsPhi->Sumw2();
THnSparseD *mHistJESJERvsPhiNoElec = new THnSparseD("mHistJESJERvsPhiNoElec", "mHistJESJERvsPhiNoElec", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsPhiNoElec->Sumw2();
THnSparseD *mHistJESJERvsEta = new THnSparseD("mHistJESJERvsEta", "mHistJESJERvsEta", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsEta->Sumw2();
THnSparseD *mHistJESJERvsEtaNoElec = new THnSparseD("mHistJESJERvsEtaNoElec", "mHistJESJERvsEtaNoElec", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsEtaNoElec->Sumw2();
THnSparseD *mHistJESJERvsdR = new THnSparseD("mHistJESJERvsdR", "mHistJESJERvsdR", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsEta->Sumw2();
THnSparseD *mHistJESJERvsdRNoElec = new THnSparseD("mHistJESJERvsdRNoElec", "mHistJESJERvsdRNoElec", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsEtaNoElec->Sumw2();


// Define multidimensional histograms
// Jets -> {pT, eta, phi, mass, energy}
const int NJetAxis = 5;
int	JetBins[NJetAxis]      =   { 100  ,  100 ,   64		  , 100  , 100  };
double JetXmin[NJetAxis]   =   { 0.0  , -5.0 ,   -TMath::Pi() , 0.0 , 0.0  };
double JetXmax[NJetAxis]   =   { 50.0 ,  5.0 ,    TMath::Pi() , 10.0, 100.0};
// all
THnSparseD *mHistJetReco = new THnSparseD("mHistJetReco", "mHistJetReco", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetReco->Sumw2();
THnSparseD *mHistJetMatch = new THnSparseD("mHistJetMatch", "mHistJetMatch", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetMatch->Sumw2();
THnSparseD *mHistJetUnMatch = new THnSparseD("mHistJetUnMatch", "mHistJetUnMatch", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetUnMatch->Sumw2();
THnSparseD *mHistJetGen = new THnSparseD("mHistJetGen", "mHistJetGen", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetGen->Sumw2();
// no Electron
THnSparseD *mHistJetRecoNoElec = new THnSparseD("mHistJetRecoNoElec", "mHistJetRecoNoElec", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetReco_noE->Sumw2();
THnSparseD *mHistJetMatchNoElec = new THnSparseD("mHistJetMatchNoElec", "mHistJetMatchNoElec", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetMatch_noE->Sumw2();
THnSparseD *mHistJetUnMatchNoElec = new THnSparseD("mHistJetUnMatchNoElec", "mHistJetUnMatchNoElec", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetUnMatch_noE->Sumw2();
THnSparseD *mHistJetGenNoElec = new THnSparseD("mHistJetGenNoElec", "mHistJetGenNoElec", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetGen_noE->Sumw2();

void WriteHistos(){

	NEvents->Write();

	numRecoJetsEventHist->Write();
	numRecoJetsNoElecEventHist->Write();
	numGenJetsEventHist->Write();
	numGenJetsNoElecEventHist->Write();
	
	JetdR->Write();
	JetMindR->Write();
	
	mHistJESJERvsPt->Write();
	mHistJESJERvsPtNoElec->Write();
	mHistJESJERvsE->Write();
	mHistJESJERvsENoElec->Write();
	mHistJESJERvsM->Write();
	mHistJESJERvsMNoElec->Write();
	mHistJESJERvsPhi ->Write();
	mHistJESJERvsPhiNoElec->Write();
	mHistJESJERvsEta->Write();
	mHistJESJERvsEtaNoElec->Write();
	mHistJESJERvsdR->Write();
	mHistJESJERvsdRNoElec->Write();

	mHistJetReco->Write();
	mHistJetMatch->Write();
	mHistJetUnMatch->Write();
	mHistJetGen->Write();

	mHistJetRecoNoElec->Write();
	mHistJetMatchNoElec->Write();
	mHistJetUnMatchNoElec->Write();
	mHistJetGenNoElec->Write();
	
}