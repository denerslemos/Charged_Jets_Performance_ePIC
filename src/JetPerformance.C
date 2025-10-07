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
/*
	    // Analyze Reconstructed Jets
		int numRecoJetsNoElec = 0;
		int numRecoJetsSingleParticle = 0;
		numRecoJetsEventHist->Fill(JetRecoType->GetSize());
		for(unsigned int ijet = 0; ijet < JetRecoType->GetSize(); ijet++) {
			LorentzVector JetReco((*JetRecoPx)[ijet], (*JetRecoPy)[ijet], (*JetRecoPz)[ijet], (*JetRecoE)[ijet]);
    		double RecoJetEta = JetReco.Eta();
    		double RecoJetPhi = JetReco.Phi();
	    	double RecoJetPt  = JetReco.Pt();
	    	double RecoJetE  = JetReco.E();
	    	double RecoJetM  = (*JetRecoM)[ijet];
	    	
	    	double RecoJet[5] = {RecoJetPt, RecoJetEta, RecoJetPhi, RecoJetM, RecoJetE};
	    	mHistJetReco->Fill(RecoJet);
	    	DijetReco.push_back(JetReco); // fill histogram for dijet studies	

			// Check if Jet Contains an Electron - Use Particle Matching to Find True PID
			bool noElectron = true;
			bool noElectronMC = true;
			double cesum = 0.0;

			for(unsigned int icjet = (*JetRecoCBegin)[ijet]; icjet < (*JetRecoCEnd)[ijet]; icjet++) {// Loop over jet constituents (particles within the jet)
				// find electron in a jet
				if((*TrkRecoPDG)[(*JetRecoCIdx)[icjet]] == 11) noElectron = false;
			    int elecIndex = -1;
			    double elecIndexWeight = -1.0;
		   		int chargePartIndex = (*JetRecoCIdx)[icjet]; // ReconstructedChargedParticle Index for m'th Jet Component
		    	for(unsigned int itrkass=0; itrkass < TrkPartAssocRec->GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations
					if((*TrkPartAssocRec)[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
					    if((*TrkPartAssocWeight)[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
							elecIndex = (*TrkPartAssocSim)[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
							elecIndexWeight = (*TrkPartAssocWeight)[itrkass];
			      		}
			  		}
		      	}
				if((*TrkMCGenPDG)[elecIndex] == 11) noElectronMC = false;

				LorentzVector TrkReco((*TrkRecoPx)[(*JetRecoCIdx)[icjet]], (*TrkRecoPy)[(*JetRecoCIdx)[icjet]], (*TrkRecoPz)[(*JetRecoCIdx)[icjet]], (*TrkRecoE)[(*JetRecoCIdx)[icjet]]);				
	  			// not looking at PID yet
				double RecoTrk[3] = {TrkReco.Pt(), TrkReco.Eta(), TrkReco.Phi()};	  	
				cesum += (*TrkRecoE)[(*JetRecoCIdx)[icjet]];		
				double zFrag[3] = {TrkReco.Pt()/RecoJetPt, RecoJetEta, RecoJetPt};
				if(RecoJetPt > 5.0 && fabs(RecoJetEta) < 3.0){
					mHistTrkConstReco->Fill(RecoTrk);		
					if(noElectron) mHistTrkConstReco_noE->Fill(RecoTrk);	
					if(noElectronMC) mHistTrkConstReco_noEMCP->Fill(RecoTrk);	
					mHistFrgConstReco->Fill(zFrag);
					if(noElectron) mHistFrgConstReco_noE->Fill(zFrag);
					if(noElectronMC) mHistFrgConstReco_noEMCP->Fill(zFrag);		
			
				}

	  		}
	  		
//	  		if(TMath::Abs(cesum - RecoJetE) > 0.000001) { cout << "bad histos" << endl;}
	
			if(noElectron){ mHistJetReco_noE->Fill(RecoJet); numRecoJetsNoElec++; DijetReco_noE.push_back(JetReco); }
			if(noElectronMC){ mHistJetReco_noEMCP->Fill(RecoJet); numRecoJetsNoElecMC++; DijetReco_noEMC.push_back(JetReco); }
			
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
		        if (dR < minDR) {
        		    minDR = dR;
            		matchedGenIdx = jjet;
        		}
    		}
			JetMindR->Fill(minDR);
			if (minDR < 0.4 && matchedGenIdx != -1) {  // If matched
				LorentzVector JetGen((*JetGenPx)[matchedGenIdx], (*JetGenPy)[matchedGenIdx], (*JetGenPz)[matchedGenIdx], (*JetGenE)[matchedGenIdx]);
				double GenJetPt = JetGen.Pt();
				double GenJetE = JetGen.E();
				double resolutionpT = RecoJetPt/GenJetPt;
				double resolutionE = RecoJetE/GenJetE;
				double JESJERvspT[3] = {resolutionpT, GenJetPt, RecoJetEta};
				mHistJESJERvsPt->Fill(JESJERvspT);
				double JESJERvsE[3] = {resolutionE, GenJetPt, RecoJetEta};
				mHistJESJERvsE->Fill(JESJERvsE);
				mHistJetMatch->Fill(RecoJet);
				if(noElectron){ mHistJetMatch_noE->Fill(RecoJet); mHistJESJERvsPt_noE->Fill(JESJERvspT); mHistJESJERvsE_noE->Fill(JESJERvsE);}
				if(noElectronMC){ mHistJetMatch_noEMCP->Fill(RecoJet); mHistJESJERvsPt_noEMCP->Fill(JESJERvspT); mHistJESJERvsE_noEMCP->Fill(JESJERvsE); }
			}else {
				mHistJetUnMatch->Fill(RecoJet);		
				if(noElectron){ mHistJetUnMatch_noE->Fill(RecoJet); }
				if(noElectronMC){ mHistJetUnMatch_noEMCP->Fill(RecoJet); }		
			}
    		
		}
		
		numRecoJetsNoElecEventHist->Fill(numRecoJetsNoElec);
		numRecoJetsNoElecMCPEventHist->Fill(numRecoJetsNoElecMC);
		
		// Analyze Generator Jets
	    int numGenJetsNoElec = 0;
    	numGenJetsEventHist->Fill(JetGenType->GetSize());		
		for(unsigned int igjet = 0; igjet < JetGenType->GetSize(); igjet++) {
			LorentzVector JetGen((*JetGenPx)[igjet], (*JetGenPy)[igjet], (*JetGenPz)[igjet], (*JetGenE)[igjet]);
    		double GenJetEta = JetGen.Eta();
    		double GenJetPhi = JetGen.Phi();
	    	double GenJetPt  = JetGen.Pt();
	    	double GenJetE  = JetGen.E();
	    	double GenJetM  = (*JetGenM)[igjet];		
		    double GenJet[5] = {GenJetPt, GenJetEta, GenJetPhi, GenJetM, GenJetE};
	    	mHistJetGen->Fill(GenJet);			

			double cesumG = 0.0;
			bool noGenElectron = true;
			
			for(unsigned int icgjet = (*JetGenCBegin)[igjet]; icgjet < (*JetGenCEnd)[igjet]; icgjet++) { // Loop over jet constituents (particles within the jet)

				double gTrkPx = (*TrkGenPx)[(*JetGenCIdx)[igjet]];
				double gTrkPy = (*TrkGenPy)[(*JetGenCIdx)[igjet]];
				double gTrkPz = (*TrkGenPz)[(*JetGenCIdx)[igjet]];
				double gTrkM = (*TrkGenM)[(*JetGenCIdx)[igjet]];
				int gTrkPDG = (*TrkGenPDG)[(*JetGenCIdx)[igjet]];
	    		double gTrkE = TMath::Sqrt(gTrkPx*gTrkPx + gTrkPy*gTrkPy + gTrkPz*gTrkPz + gTrkM*gTrkM);
	    		
	    		cesumG += gTrkE;
			    if(gTrkPDG == 11) noGenElectron = false;

				LorentzVector TrkGen(gTrkPx, gTrkPy, gTrkPz, gTrkE);				
	  			// not looking at PID yet
				double GenTrk[3] = {TrkGen.Pt(), TrkGen.Eta(), TrkGen.Phi()};	  			
				double zFrag[3] = {TrkGen.Pt()/GenJetPt, GenJetEta, GenJetPt};
				if(GenJetPt > 5.0 && fabs(GenJetEta) < 3.0){ 
					mHistTrkConstGen->Fill(GenTrk);		
					if(noGenElectron) mHistTrkConstGen_noE->Fill(GenTrk);	
					mHistFrgConstGen->Fill(zFrag); 
					if(noGenElectron) mHistFrgConstGen_noE->Fill(zFrag); 
				}						
			}

			if(noGenElectron){ mHistJetGen_noE->Fill(GenJet); numGenJetsNoElec++; }
			
		}
				
		numGenJetsNoElecEventHist->Fill(numGenJetsNoElec);
*/
		NEVENTS++;
		
	}
	cout << "Total number of events: " << NEVENTS << endl;
	
	//WriteHistos();
	
	OutFile->Close();
	

}

