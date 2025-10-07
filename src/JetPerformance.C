#include "JetTreeReader.h"  // call libraries from ROOT and C++
#include "JetHistogramPerformance.h" // call histograms

// 
void JetPerformance(TString InputFileList, TString OutputFile){

	typedef ROOT::Math::PxPyPzEVector LorentzVector;

	// Read the list of input file(s)
	fstream FileList;
	FileList.open(Form("%s",InputFileList.Data()), ios::in);
	if(!FileList.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << InputFileList.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> FileListVector;
	string FileChain;
	while(getline(FileList, FileChain)){FileListVector.push_back(FileChain.c_str());}
	FileList.close();	
	TChain *mychain = new TChain("events");
	for (std::vector<TString>::iterator listIterator = FileListVector.begin(); listIterator != FileListVector.end(); listIterator++){
		TFile *testfile = TFile::Open(*listIterator,"READ");
		if(testfile && !testfile->IsZombie() && !testfile->TestBit(TFile::kRecovered)){ // safety against corrupted files
			cout << "Adding file " << *listIterator << " to the chains" << endl; // adding files to the chains for each step
			mychain->Add(*listIterator);
		}else{cout << "File: " << *listIterator << " failed!" << endl;}
	}

	// Reading trees
	std::unique_ptr<TTreeReader> tree_reader;
	JetTreeReader(mychain, tree_reader); 
	
	// Make Output
	TFile *OutFile = TFile::Open(Form("%s",OutputFile.Data()),"RECREATE");	
	
	// Loop over events
	int NEVENTS = 0;
	while(tree_reader->Next()) {	
	
		NEvents->Fill(1);
	    if(NEVENTS%50000 == 0) cout << "Events Processed: " << NEVENTS << endl;

	    // Analyze Reconstructed Jets
		int numRecoJetsNoElec = 0;
		numRecoJetsEventHist->Fill(JetRecoType->GetSize());
		for(unsigned int ijet = 0; ijet < JetRecoType->GetSize(); ijet++) {
			// Make 4-vector
			LorentzVector JetReco((*JetRecoPx)[ijet], (*JetRecoPy)[ijet], (*JetRecoPz)[ijet], (*JetRecoE)[ijet]);
    		double RecoJetEta = JetReco.Eta();
    		double RecoJetPhi = JetReco.Phi();
	    	double RecoJetPt  = JetReco.Pt();
	    	double RecoJetE  = JetReco.E();
	    	double RecoJetM  = (*JetRecoM)[ijet];
	    	
	    	double RecoJet[5] = {RecoJetPt, RecoJetEta, RecoJetPhi, RecoJetM, RecoJetE};
	    	mHistJetReco->Fill(RecoJet); // all jets

			// Check if Jet Contains an Electron - Use Particle Matching to Find True PID
			bool hasElectron = false;
			for(unsigned int icjet = (*JetRecoCBegin)[ijet]; icjet < (*JetRecoCEnd)[ijet]; icjet++) {// Loop over jet constituents (particles within the jet)

				// Find electron in a jet
			    int elecIndex = -1;
			    double elecIndexWeight = -1.0;
		   		int chargePartIndex = (*JetRecoCIdx)[icjet]; // ReconstructedChargedParticle Index for m'th Jet Component
		    	for(unsigned int itrkass = 0; itrkass < TrkPartAssocRec->GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations
					if((*TrkPartAssocRec)[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
					    if((*TrkPartAssocWeight)[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
							elecIndex = (*TrkPartAssocSim)[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
							elecIndexWeight = (*TrkPartAssocWeight)[itrkass];
			      		}
			  		}
		      	}
				if((*TrkMCGenPDG)[elecIndex] == 11) hasElectron = true;

	  		}
	  			
			if(!hasElectron){ mHistJetRecoNoElec->Fill(RecoJet); numRecoJetsNoElec++; }
			
			// Reco-Gen matched
		    double minDR = 999.;
		    int matchedGenIdx = -1;
		    for (unsigned int jjet = 0; jjet < JetGenType->GetSize(); jjet++) {
				LorentzVector JetGen((*JetGenPx)[jjet], (*JetGenPy)[jjet], (*JetGenPz)[jjet], (*JetGenE)[jjet]);
        		double GenJetEta = JetGen.Eta();
        		double GenJetPhi = JetGen.Phi();
       			double dEta = RecoJetEta - GenJetEta;
       			double dPhi = TVector2::Phi_mpi_pi(RecoJetPhi - GenJetPhi);  // handles wraparound
				double dR = std::sqrt(dEta * dEta + dPhi * dPhi);
				JetdR->Fill(dR);
		        if (dR < minDR) { // closest dR
        		    minDR = dR;
            		matchedGenIdx = jjet;
        		}
    		}
			JetMindR->Fill(minDR);
			if (minDR < 0.4 && matchedGenIdx != -1) {  // If matched
				LorentzVector JetGen((*JetGenPx)[matchedGenIdx], (*JetGenPy)[matchedGenIdx], (*JetGenPz)[matchedGenIdx], (*JetGenE)[matchedGenIdx]);
				double GenJetPt = JetGen.Pt();
				double GenJetE = JetGen.E();
				double GenJetMass = JetGen.M();
				double GenJetPhi = JetGen.Phi();
				double GenJetEta = JetGen.Eta();

				double resolutionpT = (RecoJetPt - GenJetPt)/GenJetPt;
				double resolutionE = (RecoJetE- GenJetE)/GenJetE;
				double resolutionM = (RecoJetM - GenJetMass)/GenJetMass;
				double resolutionPhi = TVector2::Phi_mpi_pi(RecoJetPhi - GenJetPhi);
				double resolutionEta = (RecoJetEta - GenJetEta);
				double resolutiondR = sqrt(resolutionEta*resolutionEta + resolutionPhi*resolutionPhi);

				double JESJERvspT[2] = {resolutionpT, GenJetPt};
				double JESJERvsE[2] = {resolutionE, GenJetE};
				double JESJERvsM[2] = {resolutionM, GenJetPt}; // as function of pT
				double JESJERvsPhi[2] = {resolutionPhi, GenJetPt}; // as function of pT
				double JESJERvsEta[2] = {resolutionEta, GenJetPt}; // as function of pT
				double JESJERvsdR[2] = {resolutiondR, GenJetPt}; // as function of pT

				mHistJESJERvsPt->Fill(JESJERvspT);
				mHistJESJERvsE->Fill(JESJERvsE);
				mHistJESJERvsM->Fill(JESJERvsM);
				mHistJESJERvsPhi->Fill(JESJERvsPhi);
				mHistJESJERvsEta->Fill(JESJERvsEta);
				mHistJESJERvsdR->Fill(JESJERvsdR);

				mHistJetMatch->Fill(RecoJet);
				if(!hasElectron){ 
					mHistJetMatchNoElec->Fill(RecoJet); 
					mHistJESJERvsPtNoElec->Fill(JESJERvspT);
					mHistJESJERvsENoElec->Fill(JESJERvsE);
					mHistJESJERvsMNoElec->Fill(JESJERvsM);
					mHistJESJERvsPhiNoElec->Fill(JESJERvsPhi);
					mHistJESJERvsEtaNoElec->Fill(JESJERvsEta);
					mHistJESJERvsdRNoElec->Fill(JESJERvsdR);
				}

			} else {
				mHistJetUnMatch->Fill(RecoJet);		
				if(!hasElectron){ mHistJetUnMatchNoElec->Fill(RecoJet); }		
			}
    		
		}
		
		numRecoJetsNoElecEventHist->Fill(numRecoJetsNoElec);
		numRecoJetsNoElecMCPEventHist->Fill(numRecoJetsNoElecMC);
		
		// Analyze Gen Jets
	    int numGenJetsNoElec = 0;
    	numGenJetsEventHist->Fill(JetGenType->GetSize());		
		bool hasGenElectron = false;
		
		for(unsigned int igjet = 0; igjet < JetGenType->GetSize(); igjet++) {

			LorentzVector JetGen((*JetGenPx)[igjet], (*JetGenPy)[igjet], (*JetGenPz)[igjet], (*JetGenE)[igjet]);
    		double GenJetEta = JetGen.Eta();
    		double GenJetPhi = JetGen.Phi();
	    	double GenJetPt  = JetGen.Pt();
	    	double GenJetE  = JetGen.E();
	    	double GenJetM  = (*JetGenM)[igjet];		
		    double GenJet[5] = {GenJetPt, GenJetEta, GenJetPhi, GenJetM, GenJetE};
	    	mHistJetGen->Fill(GenJet);
	    	numGenJets++;
			
			// -> Check for electrons
			// Loop over jet constituents (particles within the jet)
			for(unsigned int icgjet = (*JetGenCBegin)[igjet]; icgjet < (*JetGenCEnd)[igjet]; icgjet++) { 
				int gTrkPDG = (*TrkGenPDG)[(*JetGenCIdx)[igjet]];	    		
			    if(gTrkPDG == 11) hasGenElectron = true;					
			}

			if(!hasGenElectron){ mHistJetGenNoElec->Fill(GenJet); numGenJetsNoElec++; }
			
		}
				
		numGenJetsNoElecEventHist->Fill(numGenJetsNoElec);

		NEVENTS++;
		
	}
	cout << "Total number of events: " << NEVENTS << endl;
	
	//WriteHistos();
	
	OutFile->Close();
	

}

