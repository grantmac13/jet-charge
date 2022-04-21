//Isaac Mooney, WSU - June 2019
//This file will run the initial analysis on data for the jet mass project.
//It takes in the Picos, performs selections, clusters particles, performs selections on the resulting jets,
//applies the Soft Drop grooming procedure to a copy of the jet population, fills jet trees, and writes to files.
//It is currently capable of doing this for pp or pA and should be easily extendable to AA when I start working on it.

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TFile.h>

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <chrono>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

#include "fastjet/contrib/SoftDrop.hh"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"

#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"

#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

#include "funcs.hh"
#include "params.hh"

using namespace std;
using namespace fastjet;
using namespace Analysis;
typedef fastjet::contrib::SoftDrop SD;

// -------------------------                                                                                                                                   
// Command line arguments:                                                                                                                                    
// (Defaults defined for debugging)                                                                                                                           
// [0]: output directory                                                                                                                                      
// [1]: name for the output root file containing histograms of observables                                                                                    
// [2]: flag determining if we run over ch+ne (i.e. "full") jets or just ch jets. If it's "ch", runs over ch jets. If anything else, runs over full jets.
// [3]: flag determining trigger: ppJP2, ppHT2, ppVPDMB, pAuJP2, pAuHT2, pAuBBCMB, or AA[TBD]. Also used to determine collision species
// [4]: input data: can be a single .root or a .txt or .list of root files - should always be last argument 

// DEF MAIN()
int main ( int argc, const char** argv) {
    
	//Start a timer
	TStopwatch TimeKeeper;
	TimeKeeper.Start( );
	
	//starting timing for function duration
	typedef std::chrono::high_resolution_clock clock;
	
	// Read in command line arguments
	// ------------------------------
	// Defaults
	std::string executable = "./bin/data"; // placeholder                                                                                    
	std::string outputDir = "out/"; // directory where everything will be saved                                                                                  
	std::string outFileName = "test.root"; // histograms will be saved here                                                                                      
	std::string chainList = "list.txt"; // input file: can be .root, .txt, .list                                                                                 
	std::string chainName = "JetTree"; // tree name in input file                                                                                                
	std::string trigger = "pAuJP2"; // trigger name: ppJP2, ppHT2, ppVPDMB, pAuJP2, pAuBBCMB, AA[TBD] 
	double radius = 0.4; //jet radius parameter; input value can range from 0.1 to 9.9.
	bool full = 1; //If TRUE, run over full (ch+ne) jets. If FALSE, run over ch jets.

	// Now check to see if we were given modifying arguments
	switch ( argc ) {
	case 1: // Default case
		__OUT("Using Default Settings");
		break;
	case 7: { // Custom case
		__OUT("Using Custom Settings");
		std::vector<std::string> arguments( argv+1, argv+argc );

		// Set non-default values
		// ----------------------
	
		// output and file names
		outputDir         = arguments[0];
		outFileName       = arguments[1];
		radius            = radius_str_to_double (arguments[2]);
		trigger           = arguments[3]; //ppJP2, ppHT2, ppVPDMB, pAJP2, pAuHT2, pABBCMB, or AA [TBD]
		if (arguments[4] == "ch") {full = 0;} else {full = 1;} //either ch+ne jets (default) or ch jets (if "ch")
		chainList         = arguments[5];
		
		std::cout << "Running analysis of " << arguments[4] << " jets in the " << trigger << "-triggered data. Results will be output to " << outputDir << "." << std::endl;
		std::cout << "The input file is " << chainList << " and the output file is " << outFileName << "." << std::endl;
		
		break;
	}
	default: { // Error: invalid custom settings
		__ERR("Invalid number of command line arguments");
		return -1;
		break;
	}
	}

	//Setting up specifics of analysis based on the flags that were received above!
	string badtows = "", badruns = "";
	double vzdiff = -1;
	if (trigger.find("pp") != string::npos) {badtows = det_badTowers; badruns = dat_bad_run_list; vzdiff = det_vZDiff;}
	if (trigger.find("pA") != string::npos) {badtows = pAu_badTowers; badruns = pAu_bad_run_list; vzdiff = pAu_vZDiff;}
	if (trigger.find("AA") != string::npos) {badtows = ""; badruns = ""; vzdiff = -1;} //TBD

	//in place for now; will encapsulate in a function if it gets much more involved. Hardcodes the trigger IDs.                                                 
	int tID1 = -9999, tID2 = -9999, tID3 = -9999;
	if (trigger == "ppJP2") {tID1 = tppJP2; tID2 = -8888; tID3 = -8888;} //-8888 just ensures it won't accidentally match a trigger
	if (trigger == "ppHT2") {tID1 = tppHT2a; tID2 = tppHT2b; tID3 = tppHT2c;}
	if (trigger == "ppVPDMB") {tID1 = tppVPDMB_nobsmd; tID2 = -8888; tID3 = -8888;}
	if (trigger == "pAuJP2") {tID1 = tpAuJP2a; tID2 = tpAuJP2b; tID3 = -8888;}
	if (trigger == "pAuHT2") {tID1 = tpAuHT2a; tID2 = tpAuHT2b; tID3 = -8888;}
	if (trigger == "pAuBBCMB") {tID1 = tpAuBBCMBa; tID2 = tpAuBBCMBb; tID3 = -8888;}

	// Build our input now
	// --------------------
	TChain* chain = new TChain( chainName.c_str() );

	// Check to see if the input is a .root file or a .txt
	bool inputIsRoot = Analysis::HasEnding( chainList.c_str(), ".root" );
	bool inputIsTxt  = Analysis::HasEnding( chainList.c_str(), ".txt"  );
	bool inputIsList = Analysis::HasEnding( chainList.c_str(), ".list" );

	// If its a recognized file type, build the chain
	// If its not recognized, exit
	if ( inputIsRoot ) { chain->Add( chainList.c_str() ); }
	else if ( inputIsTxt )  { chain = TStarJetPicoUtils::BuildChainFromFileList( chainList.c_str() ); }
	else if ( inputIsList)  { chain = TStarJetPicoUtils::BuildChainFromFileList( chainList.c_str() ); }
	else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }


	// Build the event structure w/ cuts
	// ---------------------------------
	TStarJetPicoReader * reader = new TStarJetPicoReader();
	InitReader(reader, chain, nEvents, "All"/*det_triggerString*/, det_absMaxVz, vzdiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, 0.9999, false, badtows, badruns); //hc = 0.9999

	// Data classes
	// ------------
	TStarJetVectorContainer<TStarJetVector>* container;
	TStarJetVector* sv; // TLorentzVector* would be sufficient
	TStarJetPicoEventHeader* header;
	TStarJetPicoEvent* event;
	TStarJetPicoEventCuts* EventCuts = reader->GetEventCuts();//will use this for hardcoding checks if events passed selections
	
	// Histograms
	// ----------
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TH3::SetDefaultSumw2();

/*
	// Grant's Histograms
	//Event histograms
	TH1D* LeadPt = new TH1D("leading_jet_pt", "pT of leading jets", 40, 0, 100);

	TH2D* cvp_T = new TH2D("Transverse_Vertex", "x,y position of collision vertex", 50, -.1, .1, 50, -.3, .2);
	TH1D* cvp_z = new TH1D("Z_vertex", "z position of collision vertex", 200, -30, 30);
	TH1D* multip = new TH1D("Multiplicity", "Number of Charged Tracks", 100, 0, 50);
	TH1D* refmultip = new TH1D("RefMult", "Reference Multiplicity", 31, 0, 30);
	TH1D* event_zdc = new TH1D("zdc_coincidence", "ZDC Coincidence", 50, 6000.0, 16000.0);
	TH1D* event_bbc = new TH1D("bbc_coincidence", "BBC Coincidence", 50, 300000.0, 500000.0);
	
	TH1D* tow_mult = new TH1D("tow_mult", "tower multiplicity per event", 250, 0, 500);
	TH1D* tr_mult = new TH1D("tr_mult", "track multiplicity per event", 150, 0, 150);    
	
	// Tower histograms
	TH1D* E_hist = new TH1D("Energy", "Energy", 25, 0, 2);
	TH1D* Et_hist = new TH1D("E_T", "Transverse Energy", 30, 0.2, 30);
	TH1D* eta_hist = new TH1D("Eta", "eta", 40, -1, 1);
	TH1D* eta_corr_hist = new TH1D("corrected_eta", "eta_corr", 40, -1, 1);
	TH1D* phi_hist = new TH1D("Phi", "phi", 120, -TMath::Pi(), TMath::Pi());

	TH1D* bad_tow_check = new TH1D("bad_tow_check", "Check for bad towers", 4800, 1, 4800);
	
	// Track histograms
	TH1D* track_pt = new TH1D("track_pt", "track pT", 24, 0.2, 5);

///////////////////////////////////////////////////////////////////////
	TH1D* ptrpt = new TH1D("ptrpt", "track pT", 50, 0.2, 30);
	TH1D* ntrpt = new TH1D("ntrpt", "track pT", 50, 0.2, 30);
*/

	const int nPtBins = 8;
	double ptLo[nPtBins] = {0.2, 0.4, 0.6, 1.0, 2.0, 4.0,  6.0, 10.0};
	double ptHi[nPtBins] = {0.4, 0.6, 1.0, 2.0, 4.0, 6.0, 10.0, 30.0};



	const int nChBins = 2;
	TString chBin[nChBins] = {"_pos", "_neg"};
	
	TString ptBinName[nPtBins] = {"02_04GeV", "04_06GeV", "06_1GeV", "1_2GeV", "2_4GeV", "4_6GeV", "6_10GeV", "10_30GeV"};


	const int nJetPtBins = 9;
	double jetPtLo[nJetPtBins] = { 0.0, 20.0, 20.0, 30.0, 15.0, 10.0, 30.0, 20.0, 25.0};
	double jetPtHi[nJetPtBins] = {50.0, 30.0, 50.0, 50.0, 20.0, 15.0, 40.0, 25.0, 30.0};

	TString jetPtBinName[nJetPtBins] = {"JetPtAll", "JetPtOver20GeVUnder30GeV", "JetPtOver20GeV", "JetPtOver30GeV",
					"JetPtOver15GeVUnder20GeV", "JetPtOver10GeVUnder15GeV", "JetPtOver30GeVUnder40GeV",
					"JetPtOver20GeVUnder25GeV", "JetPtOver25GeVUnder30GeV"};

	// jet histograms
	TH1D* jet_pt_hist[nJetPtBins];
	TH1D* hnjets[nJetPtBins];
	TH1D* hkappa00[nJetPtBins];
	TH1D* hkappa05[nJetPtBins];
	TH1D* hPtNeutrals[nJetPtBins];
	TH1D* hNEF[nJetPtBins];
	
	TH2D* hPt1Pt2[nJetPtBins][nChBins][nChBins];    // pT1-pT2 (++, +-, -+, - -) for each jetPtBinName
	TH1D* hPtCh[nJetPtBins][nChBins];
	TH1D* hZCh[nJetPtBins][nChBins];
	

	for(int i = 0; i < nJetPtBins; i++){
		TString jetPtName;
		hnjets[i] = new TH1D("hnjet" + jetPtBinName[i], "N jets", 60, 0, 50);
		hkappa00[i] = new TH1D("hk00" + jetPtBinName[i], "Net Charge", 9, -4.5, 4.5);
		hkappa05[i] = new TH1D("hk05" + jetPtBinName[i], "Jet Charge", 20, -2.5, 2.5);
		hPtNeutrals[i] = new TH1D("hPtNeutrals" + jetPtBinName[i], "pT neutral constituents", 41, 0, 50);
		hNEF[i] = new TH1D("hNeuEFrac" + jetPtBinName[i], "NEF", 30, 0, 1);

		if(jetPtHi[i] > 45.0){
			jetPtName  = Form("hJetPt%1.0f", jetPtLo[i]);
		}
		else{
			jetPtName  = Form("hJetPt%1.0f_%1.0f", jetPtLo[i], jetPtHi[i]);
		}
		jet_pt_hist[i] = new TH1D(jetPtName, "jet pT", 60, 0, 50);
		
		hPtCh[i][0] = new TH1D("hPtCh" + jetPtBinName[i] + "Pos", "", 31, 0, 30);
		hPtCh[i][1] = new TH1D("hPtCh" + jetPtBinName[i] + "Neg", "", 31, 0, 30);
		hZCh[i][0] = new TH1D("hZCh" + jetPtBinName[i] + "Pos", "", 31, 0, 1);
		hZCh[i][1] = new TH1D("hZCh" + jetPtBinName[i] + "Neg", "", 31, 0, 1);
	}




	// Trees
	// -----
	int dummy_int;
	double dummy_double;

	//variables to link to branches in the tree
	double n_jets; double bbc_east_rate, bbc_east_sum;
	vector<double> Pt; vector<double> Eta; vector<double> Phi; vector<double> M; vector<double> E;


	//contains all (important) jet observables for both groomed and ungroomed jets
	TTree *eventTree = new TTree("event","event");
	eventTree->Branch("n_jets", &n_jets);
	eventTree->Branch("bbc_east_rate", &bbc_east_rate); //~lumi
	eventTree->Branch("bbc_east_sum", &bbc_east_sum); //~centrality
	eventTree->Branch("Pt", &Pt); eventTree->Branch("Eta",&Eta); eventTree->Branch("Phi",&Phi); eventTree->Branch("M",&M); eventTree->Branch("E",&E);


	// JET CHARGE
	// for jet charge need the following:
	// for both particle and detector level
	// jet pt, constituent charge, constituent pt
	// will match jets before filling trees
	//
	vector<double> jetpt; vector<vector<double> > conspt; vector<vector<double> > chcons;

	TTree *jetChargeTree = new TTree("jetChargeTree", "jetChargeTree");
	jetChargeTree->Branch("jetpt", &jetpt);
	jetChargeTree->Branch("conspt", &conspt);
	jetChargeTree->Branch("chcons", &chcons);


	// BALANCE FUNCTION
	// fill tree with constituents' pt, charge
	// think about structure of root macro to calculate R2, A, B2
//	vector<double> jetpt_bf; vector<vector<double> > bf_conspt; vector<vector<double> > bf_chcons;

//	TTree *bfTree = new TTree("bfTree", "bfTree");
//	bfTree->Branch("event", &event); // event number: used for checking if constituents a and b are in same jet
//	bfTree->Branch("jetpt_bf", &jetpt_bf); // jet pt: used for checking if constituents a and b are in same jet
//	bfTree->Branch("bf_conspt", &bf_conspt);
//	bfTree->Branch("bf_chcons",  &bf_chcons);

	// Helpers
	// -------
	vector<PseudoJet> particles;

	unsigned nJets = 0;

	// Constituent selectors
	// ---------------------
	Selector select_track_rap  = fastjet::SelectorAbsEtaMax(max_track_rap);    // was AbsRapMax (4/29/2020)
	Selector select_pt_min     = fastjet::SelectorPtMin( partMinPt );
	Selector select_pt_max     = fastjet::SelectorPtMax( partMaxPt );
	Selector spart = select_track_rap && select_pt_min && select_pt_max;

	// Jet candidate selectors
	// -----------------------
	Selector select_jet_rap     = fastjet::SelectorAbsEtaMax(max_rap);     // was AbsRapMax (4/29/2020)
	Selector select_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
	Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
//	Selector select_jet_m_min   = fastjet::SelectorMassMin( mass_min );
	Selector sjet = select_jet_rap && select_jet_pt_min && select_jet_pt_max;/* && select_jet_m_min;*/

	// Choose a jet and area definition
	// --------------------------------
	JetDefinition jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, radius/*R*/);




	// create an area definition for the clustering
	//----------------------------------------------------------
	// ghosts should go up to the acceptance of the detector or
	// (with infinite acceptance) at least 2R beyond the region
	// where you plan to investigate jets.
	GhostedAreaSpec area_spec = fastjet::GhostedAreaSpec( ghost_maxrap, ghost_repeat, ghost_area );
	AreaDefinition  area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

	//Creating SoftDrop grooming object
	contrib::SoftDrop sd(Beta,z_cut,R0);
	cout << "SoftDrop groomer is: " << sd.description() << endl;

	//for later use looking up PDG masses using particle PID
	//TDatabasePDG *pdg = new TDatabasePDG();

	cout << "Performing analysis." << endl;
	// Cycle through events
	// --------------------  
	int nEventsUsed = 0;
  
	try{
	while ( reader->NextEvent() ) {

		//clearing vectors
		Pt.clear(); Eta.clear(); Phi.clear(); M.clear(); E.clear();
		//initializing variables to -9999
		n_jets = -9999; bbc_east_rate = -9999; bbc_east_sum = -9999;


		jetpt.clear();	// vectors of jet pt doubles
		conspt.clear();	// vector of constituents' pt for each jet
		chcons.clear();	// vector of constituents' charge for each jet

//		jetpt_bf.clear();
//		bf_conspt.clear();
//		bf_chcons.clear();


		reader->PrintStatus(10);

		//get the event header
		event = reader->GetEvent();
		header = event->GetHeader();

		particles.clear();

		// Get the output container from the reader
		// ----------------------------------------
		container = reader->GetOutputContainer();
	
		//     //~~~~~~~~~~~~~~~~~~~~~~~~~~~Skipping undesired events!~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                    

		//if the event lacks the desired trigger, skip it                                                                                                        
		//see function "SetTriggers()" for assignment of tID1, tID2 (or above until I write it)                                                                  
		if ( ! (header->HasTriggerId(tID1) || header->HasTriggerId(tID2) || header->HasTriggerId(tID3) ) ) {//cout << "DEBUG: skipping this event because it lacks appropriate triggers. Does it have trigger ID " << tID1 << "? " << header->HasTriggerId(tID1) << endl;
			continue;
		}
		if (trigger.find("pA") != string::npos) {//removing some runs by hand in pA until we have bad run/tower lists
			//TEMPORARILY SKIPPING THESE RUNS for pA [should define these runIDs somewhere later so they're not magic numbers]                                     
			if (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001) {cout << "DEBUG: should never see this for pp!" << endl; continue;}
			//something weird happened to the towers in run 16135032 (and it looks like it started at the end of run 16135031), so excluding both                  
			if (header->GetRunId() == 16135031 || header->GetRunId() == 16135032) {cout << "DEBUG: should never see this for pp!" << endl; continue;}
			//the event cuts don't check if the vzdiff is acceptable, so I have to hardcode it here. UPDATE: I believe Nick updated the eventstructuredAu to check this condition, so this line should be redundant now
			//	if (!EventCuts->IsVertexZDiffOK(event)) {cout << "DEBUG: shouldn't see this now!" << endl; continue;}
			//Above 64000 seems like detectors saturate (tower multiplicity explodes).                                                                             
			if (header->GetBbcAdcSumEast() >= pAu_BBCE_ADC_sum_max) {continue;}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~STARTING ANALYSIS ON THE EVENT!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      
		// Transform TStarJetVectors into (FastJet) PseudoJets; assign mass to particles
		// ----------------------------------------------------------
		GatherParticles(container, sv, particles, full, 0);  //, pdg); //"pdg" here finds the rest mass for particles with a given PID; 0 means detector-level.

		// Analysis
		// --------
		// Apply selector, spart (defined above, under "Constituent selectors"), to the full particle set
		vector<PseudoJet> good_parts = spart( particles );


/*
		// track QA
//		TStarJetPicoReader track;
		TList *track_list = reader->GetListOfSelectedTracks();
		TIter iterTrack(track_list);
		int ntracks = 0;
		int ntowers = 0;
*/

		// find corresponding jets with soft constituents
		// ----------------------------------------------
		ClusterSequence/*Area*/ csa ( good_parts, jet_def/*, area_def */); // WITHOUT background subtraction

		vector<PseudoJet> good_jets_Initial = fastjet::sorted_by_pt(sjet(csa.inclusive_jets()));//applies jet selector to clustered jets

		vector<PseudoJet> good_jets;

		// apply NEF... as the function's name suggests
		// cut off at 90%
		ApplyNEFSelection(good_jets_Initial, good_jets);

/*
		if (good_jets.size() != 0) {
			//SoftDrop is a groomer not a tagger, so if we have at least one ungroomed jet, we should also have a SoftDrop'd jet.
			nJets += good_jets.size();
			n_jets = good_jets.size();
			bbc_east_rate = header->GetBbcEastRate();
			bbc_east_sum = header->GetBbcAdcSumEast(); 
			for (int i = 0; i < n_jets; ++ i) {
				Pt.push_back(good_jets[i].pt()); Eta.push_back(good_jets[i].eta()); Phi.push_back(good_jets[i].phi());
				M.push_back(good_jets[i].m()); E.push_back(good_jets[i].e());
			}//for loop over jets
		}//if we had jets
*/
		
		// And we're done!
		// -----------------------------
		if (good_jets.size() != 0) {
			nEventsUsed++; //this event was accepted and had at least one jet passing all criteria
			eventTree->Fill();
		}

		if (good_jets.size() > 0){  // if jets found

			///////////////////// fill jetcharge tree and balance function tree or histograms
			for(int i = 0; i < good_jets.size(); i++){
				for(int jj = 0; jj < nJetPtBins; jj++){
					if(jetPtLo[jj] < good_jets[i].pt() && good_jets[i].pt() < jetPtHi[jj]){
						double k00 = 0.0; double k05 = 0.0;
						double djetpt = good_jets[i].pt();
						for(int j = 0; j < good_jets[i].constituents().size(); j++){
							double dch = good_jets[i].constituents()[j].user_index();
							double dconspt = good_jets[i].constituents()[j].pt();

							k00 += dch;
							k05 += (dconspt/djetpt)*dch;
						}
						hkappa00[jj]->Fill(k00);
						hkappa05[jj]->Fill(k05);
					}
				}

				jetpt.push_back(good_jets[i].pt());
//				jetpt_bf.push_back(good_jets[i].pt());

				vector<double> jc_ch; vector<double> jc_pt;
//				vector<double> bf_ch; vector<double> bf_pt;
				for(int j = 0; j < good_jets[i].constituents().size(); j++){
					jc_ch.push_back(good_jets[i].constituents()[j].user_index());
//					bf_ch.push_back(good_jets[i].constituents()[j].user_index());
					jc_pt.push_back(good_jets[i].constituents()[j].pt());
//					bf_pt.push_back(good_jets[i].constituents()[j].pt());
				}
				chcons.push_back(jc_ch);
				conspt.push_back(jc_pt);
				
//				bf_chcons.push_back(bf_ch);
//				bf_conspt.push_back(bf_pt);

			}  // i loop over jets
		} // if jets found


		if(good_jets.size() != 0){
			jetChargeTree->Fill();
//			bfTree->Fill();
		}


/*
				for(int j = 0; j < nJetPtBins; j++){	// j loop only used for histogram array
					bool jetPtPass = false;
					if(jetPtHi[j] > 40.0){	// if jetPtHi = 50GeV (arbitrarily set), always pass
						jetPtPass = true;
					}
					else if(good_jets[i].pt() <= jetPtHi[j]){
						jetPtPass = true;
		 			}
					if(jetPtLo[j] <= good_jets[i].pt() && jetPtPass){	// jetPtPass checks upper bound
						hnjets[j]->Fill(1);
						jet_pt_hist[j]->Fill(good_jets[i].pt());
						double jetch = 0;
						double weighted_jetch = 0;
						double nef = 0;
//						int n_neu = 0;
						for(int k = 0; k < good_jets[i].constituents().size(); k++){
							int q = -1;
							if(good_jets[i].constituents()[k].user_index() == 1){
								q = 0;
								jetch += good_jets[i].constituents()[k].user_index();
								weighted_jetch += (pow(1.0/good_jets[i].pt(), 0.5))*good_jets[i].constituents()[k].user_index()*(pow(good_jets[i].constituents()[k].pt(), 0.5));
								hPtChCons[0]->Fill(good_jets[i].constituents()[k].pt());
							}
							if(good_jets[i].constituents()[k].user_index() == -1){
								q = 1;
								jetch += good_jets[i].constituents()[k].user_index();
								weighted_jetch += (pow(1.0/good_jets[i].pt(), 0.5))*good_jets[i].constituents()[k].user_index()*(pow(good_jets[i].constituents()[k].pt(), 0.5));
								hPtChCons[1]->Fill(good_jets[i].constituents()[k].pt());
							}
							if(good_jets[i].constituents()[k].user_index() == 0){
								nef += good_jets[i].constituents()[k].pt();
								hPtNeutrals[j]->Fill(good_jets[i].constituents()[k].pt());
//								n_neu += 1;
							}
//							jetch += good_jets[i].constituents()[k].user_index();
//							weighted_jetch += (pow(1.0/good_jets[i].pt(), 0.5))*good_jets[i].constituents()[k].user_index()*(pow(good_jets[i].constituents()[k].pt(), 0.5));

							if(q == 0 || q == 1){
								hPtCh[j][q]->Fill(good_jets[i].constituents()[k].pt());
								hZCh[j][q]->Fill(good_jets[i].constituents()[k].pt()/good_jets[i].pt());
							}
							for(int m = 0; m < nPtBins; m++){
								if(ptLo[m] <= good_jets[i].constituents()[k].pt() && good_jets[i].constituents()[k].pt() <= ptHi[m]){
									if(q == 0 || q == 1){
										hnCons[q][m][j]->Fill(good_jets[i].constituents()[k].pt());
									}
									for(int l = 0; l < good_jets[i].constituents().size(); l++){
										int s = -1; // charge index for constituent 2
										if(good_jets[i].constituents()[l].user_index() == 1){
											s = 0;
										}
										if(good_jets[i].constituents()[l].user_index() == -1){
											s = 1;
										}
										if( (q == 0 || q == 1) && (s == 0 || s == 1)  && (k != l)){  // if there is a charged trigger and associated particle
											hpTConsJetPtwithTrigger[j][s]->Fill(good_jets[i].constituents()[l].pt());
											consByCharge[q][s][m][j]->Fill(good_jets[i].constituents()[l].pt());
					    
											hChConsCorr[j][0][q][s]->Fill(abs(good_jets[i].constituents()[k].pt() - good_jets[i].constituents()[l].pt()),
															  good_jets[i].constituents()[k].eta() - good_jets[i].constituents()[l].eta());
											hChConsCorr[j][1][q][s]->Fill(abs(good_jets[i].constituents()[k].pt() - good_jets[i].constituents()[l].pt()),
															  good_jets[i].constituents()[k].delta_phi_to(good_jets[i].constituents()[l]));
											hChConsCorr[j][2][q][s]->Fill(abs(good_jets[i].constituents()[k].pt() - good_jets[i].constituents()[l].pt()),
															  good_jets[i].constituents()[k].delta_R(good_jets[i].constituents()[l]));

											hChConsCorr[j][3][q][s]->Fill(good_jets[i].constituents()[k].eta() - good_jets[i].constituents()[l].eta(),
														      good_jets[i].constituents()[k].delta_phi_to(good_jets[i].constituents()[l]));
											hChConsCorr[j][4][q][s]->Fill(good_jets[i].constituents()[k].eta() - good_jets[i].constituents()[l].eta(),
														      good_jets[i].constituents()[k].delta_R(good_jets[i].constituents()[l]));

											hChConsCorr[j][5][q][s]->Fill(good_jets[i].constituents()[k].delta_phi_to(good_jets[i].constituents()[l]),
                                                                					              good_jets[i].constituents()[k].delta_R(good_jets[i].constituents()[l]));


											//hPt1Pt2[nJetPtBins][nChBins][nChBins]
											hPt1Pt2[j][q][s]->Fill(good_jets[i].constituents()[k].pt(), good_jets[i].constituents()[l].pt());
										}
									} // constituent1 charge
								}
							} // l loop (constituent2)
						} // k loop (constituent1)
						hjetNetCharge[j]->Fill(jetch);
						hjetJetCharge[j]->Fill(weighted_jetch);
						hNEF[j]->Fill(nef/good_jets[i].pt());
					} // fix jetPtHi bound at 50GeV
				} // j loop
*/

	} // Event loop
	}catch ( std::exception& e) {
	std::cerr << "Caught " << e.what() << std::endl;
	return -1;
	}

	// Output
	// ------
	TFile* fout = new TFile((outputDir + outFileName).c_str(), "RECREATE");

	// Close up shop
	// -------------

	// jet

//	hPtChCons[0]->Write();
//	hPtChCons[1]->Write();

	for(int i = 0; i < nJetPtBins; i++){
//		jet_pt_hist[i]->Write();
//		hnjets[i]->Write();
//		hPtNeutrals[i]->Write();
//		hNEF[i]->Write();


		hkappa00[i]->Write();
		hkappa05[i]->Write();

/*
		for(int j = 0; j < nChBins; j++){
			hPtCh[i][j]->Write();
			hZCh[i][j]->Write();
			for(int k = 0; k < nChBins; k++){
				hPt1Pt2[i][j][k]->Write();
			}
		}
*/
	}




	// event tree
	eventTree->Write();
	
	// tree with jet charge information
	jetChargeTree->Write();

	//fout->Write();
	fout->Close();

	cout << "In " << nEventsUsed << " events, found " << endl
	     << nJets << " jets above 5 GeV, with constituents above 0.2 GeV," << endl
	     << "for an average of " << nJets/(double)nEventsUsed << " jets per event" << endl;

	cout << "Wrote to " << fout->GetName() << endl;
	cout << "Bye :)" << endl;

	return 0;
}
  
