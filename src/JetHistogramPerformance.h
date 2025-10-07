TH1D *NEvents = new TH1D("NEvents","",2,0.,2.); NEvents->Sumw2();
TH1D *numRecoJetsEventHist = new TH1D("numRecoJetsEvent","",20,0.,20.); numRecoJetsEventHist->Sumw2();
TH1D *numRecoJetsNoElecEventHist = new TH1D("numRecoJetsNoElecEvent","",20,0.,20.); numRecoJetsNoElecEventHist->Sumw2();
TH1D *numRecoJetsNoElecMCPEventHist = new TH1D("numRecoJetsNoElecMCPEventHist","",20,0.,20.); numRecoJetsNoElecMCPEventHist->Sumw2();

TH1D *numGenJetsEventHist = new TH1D("numGenJetsEventHist","",20,0.,20.); numGenJetsEventHist->Sumw2();
TH1D *numGenJetsNoElecEventHist = new TH1D("numGenJetsNoElecEventHist","",20,0.,20.); numGenJetsNoElecEventHist->Sumw2();

TH1D *JetdR = new TH1D("JetdR","",100,0.,10.); JetdR->Sumw2();
TH1D *JetMindR = new TH1D("JetMindR","",100,0.,10.); JetMindR->Sumw2();

// JER/JES -> {reco/match, pT or E, eta}
const int NJESJERAxis = 3;
int	JESJERBins[NJESJERAxis]      =   { 400  , 100  , 100};
double JESJERXmin[NJESJERAxis]   =   { 0.0  , 0.0  , -5.0};
double JESJERXmax[NJESJERAxis]   =   { 2.0  , 100.0, 5.0};
THnSparseD *mHistJESJERvsPt = new THnSparseD("mHistJESJERvsPt", "mHistJESJERvsPt", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsPt->Sumw2();
THnSparseD *mHistJESJERvsE = new THnSparseD("mHistJESJERvsE", "mHistJESJERvsE", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsE->Sumw2();
THnSparseD *mHistJESJERvsPt_noE = new THnSparseD("mHistJESJERvsPt_noE", "mHistJESJERvsPt_noE", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsPt_noE->Sumw2();
THnSparseD *mHistJESJERvsE_noE = new THnSparseD("mHistJESJERvsE_noE", "mHistJESJERvsE_noE", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsE_noE->Sumw2();
THnSparseD *mHistJESJERvsPt_noEMCP = new THnSparseD("mHistJESJERvsPt_noEMCP", "mHistJESJERvsPt_noEMCP", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsPt_noEMCP->Sumw2();
THnSparseD *mHistJESJERvsE_noEMCP = new THnSparseD("mHistJESJERvsE_noEMCP", "mHistJESJERvsE_noEMCP", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsE_noEMCP->Sumw2();


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
THnSparseD *mHistJetReco_noE = new THnSparseD("mHistJetReco_noE", "mHistJetReco_noE", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetReco_noE->Sumw2();
THnSparseD *mHistJetMatch_noE = new THnSparseD("mHistJetMatch_noE", "mHistJetMatch_noE", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetMatch_noE->Sumw2();
THnSparseD *mHistJetUnMatch_noE = new THnSparseD("mHistJetUnMatch_noE", "mHistJetUnMatch_noE", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetUnMatch_noE->Sumw2();
THnSparseD *mHistJetGen_noE = new THnSparseD("mHistJetGen_noE", "mHistJetGen_noE", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetGen_noE->Sumw2();
// no Electron from MC particles
THnSparseD *mHistJetReco_noEMCP = new THnSparseD("mHistJetReco_noEMCP", "mHistJetReco_noEMCP", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetReco_noEMCP->Sumw2();
THnSparseD *mHistJetMatch_noEMCP = new THnSparseD("mHistJetMatch_noEMCP", "mHistJetMatch_noEMCP", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetMatch_noEMCP->Sumw2();
THnSparseD *mHistJetUnMatch_noEMCP = new THnSparseD("mHistJetUnMatch_noEMCP", "mHistJetUnMatch_noEMCP", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetUnMatch_noEMCP->Sumw2();

// Constituents -> {pT, eta, phi}
const int NTrkAxis = 3;
int	TrkBins[NTrkAxis]      =   { 100  ,  100 ,    64		  };
double TrkXmin[NTrkAxis]   =   { 0.0  , -5.0 ,   -TMath::Pi() };
double TrkXmax[NTrkAxis]   =   { 50.0 ,  5.0 ,    TMath::Pi() };
THnSparseD *mHistTrkConstReco = new THnSparseD("mHistTrkConstReco", "mHistTrkConstReco", NTrkAxis, TrkBins, TrkXmin, TrkXmax); mHistTrkConstReco->Sumw2();
THnSparseD *mHistTrkConstReco_noE = new THnSparseD("mHistTrkConstReco_noE", "mHistTrkConstReco_noE", NTrkAxis, TrkBins, TrkXmin, TrkXmax); mHistTrkConstReco_noE->Sumw2();
THnSparseD *mHistTrkConstReco_noEMCP = new THnSparseD("mHistTrkConstReco_noEMCP", "mHistTrkConstReco_noEMCP", NTrkAxis, TrkBins, TrkXmin, TrkXmax); mHistTrkConstReco_noEMCP->Sumw2();
THnSparseD *mHistTrkConstGen = new THnSparseD("mHistTrkConstGen", "mHistTrkConstGen", NTrkAxis, TrkBins, TrkXmin, TrkXmax); mHistTrkConstGen->Sumw2();
THnSparseD *mHistTrkConstGen_noE = new THnSparseD("mHistTrkConstGen_noE", "mHistTrkConstGen_noE", NTrkAxis, TrkBins, TrkXmin, TrkXmax); mHistTrkConstGen_noE->Sumw2();

// Jet Fragmentation
const int NFragAxis = 3;
int	FragBins[NTrkAxis]      =   { 100  ,  100 , 100};
double FragXmin[NTrkAxis]   =   { 0.0  , -5.0 , 0.0};
double FragXmax[NTrkAxis]   =   { 1.0  ,  5.0 , 50.0};

THnSparseD *mHistFrgConstReco = new THnSparseD("mHistFrgConstReco", "mHistFrgConstReco", NFragAxis, FragBins, FragXmin, FragXmax); mHistFrgConstReco->Sumw2();
THnSparseD *mHistFrgConstReco_noE = new THnSparseD("mHistFrgConstReco_noE", "mHistFrgConstReco_noE", NFragAxis, FragBins, FragXmin, FragXmax); mHistFrgConstReco_noE->Sumw2();
THnSparseD *mHistFrgConstReco_noEMCP = new THnSparseD("mHistFrgConstReco_noEMCP", "mHistFrgConstReco_noEMCP", NFragAxis, FragBins, FragXmin, FragXmax); mHistFrgConstReco_noEMCP->Sumw2();
THnSparseD *mHistFrgConstGen = new THnSparseD("mHistFrgConstGen", "mHistFrgConstGen", NFragAxis, FragBins, FragXmin, FragXmax); mHistFrgConstGen->Sumw2();
THnSparseD *mHistFrgConstGen_noE = new THnSparseD("mHistFrgConstGen_noE", "mHistFrgConstGen_noE", NFragAxis, FragBins, FragXmin, FragXmax); mHistFrgConstGen_noE->Sumw2();


// Dijets {etaDijet, Xj, jet pT average, leading jet pT, subleading jet pT}
const int NDiJetAxis = 5;
int	DiJetBins[NDiJetAxis]      =   {  50    , 10 , 200 , 200 , 200 };
double DiJetXmin[NDiJetAxis]   =   {  -5.0  , 0.0, 0.0 , 0.0 , 0.0 };
double DiJetXmax[NDiJetAxis]   =   {  5.0   , 1.0, 50.0, 50.0, 50.0};
THnSparseD *mHistDiJetReco = new THnSparseD("mHistDiJetReco", "mHistDiJetReco", NDiJetAxis, DiJetBins, DiJetXmin, DiJetXmax); mHistDiJetReco->Sumw2();
THnSparseD *mHistDiJetReco_noE = new THnSparseD("mHistDiJetReco_noE", "mHistDiJetReco_noE", NDiJetAxis, DiJetBins, DiJetXmin, DiJetXmax); mHistDiJetReco_noE->Sumw2();
THnSparseD *mHistDiJetReco_noEMCP = new THnSparseD("mHistDiJetReco_noEMCP", "mHistDiJetReco_noEMCP", NDiJetAxis, DiJetBins, DiJetXmin, DiJetXmax); mHistDiJetReco_noEMCP->Sumw2();
THnSparseD *mHistDiJetGen = new THnSparseD("mHistDiJetGen", "mHistDiJetGen", NDiJetAxis, DiJetBins, DiJetXmin, DiJetXmax); mHistDiJetGen->Sumw2();
THnSparseD *mHistDiJetGen_noE = new THnSparseD("mHistDiJetGen_noE", "mHistDiJetGen_noE", NDiJetAxis, DiJetBins, DiJetXmin, DiJetXmax); mHistDiJetGen_noE->Sumw2();


void WriteHistos(){

	NEvents->Write();

	numRecoJetsEventHist->Write();
	numRecoJetsNoElecEventHist->Write();
	numRecoJetsNoElecMCPEventHist->Write();
	numGenJetsEventHist->Write();
	numGenJetsNoElecEventHist->Write();
	
	JetdR->Write();
	JetMindR->Write();
	
	mHistJESJERvsPt->Write();
	mHistJESJERvsE->Write();
	mHistJESJERvsPt_noE->Write();
	mHistJESJERvsE_noE->Write();
	mHistJESJERvsPt_noEMCP->Write();
	mHistJESJERvsE_noEMCP->Write();

	mHistJetReco->Write();
	mHistJetMatch->Write();
	mHistJetUnMatch->Write();
	mHistJetGen->Write();

	mHistJetReco_noE->Write();
	mHistJetMatch_noE->Write();
	mHistJetUnMatch_noE->Write();
	mHistJetGen_noE->Write();
	
	mHistJetReco_noEMCP->Write();
	mHistJetMatch_noEMCP->Write();
	mHistJetUnMatch_noEMCP->Write();
	
	mHistTrkConstReco->Write();
	mHistTrkConstReco_noE->Write();
	mHistTrkConstReco_noEMCP->Write();
	mHistTrkConstGen->Write();
	mHistTrkConstGen_noE->Write();

	mHistFrgConstReco->Write();
	mHistFrgConstReco_noE->Write();
	mHistFrgConstReco_noEMCP->Write();
	mHistFrgConstGen->Write();
	mHistFrgConstGen_noE->Write();

	mHistDiJetReco->Write();
	mHistDiJetReco_noE->Write();
	mHistDiJetReco_noEMCP->Write();
	mHistDiJetGen->Write();
	mHistDiJetGen_noE->Write();
	
}