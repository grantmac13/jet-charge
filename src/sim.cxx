// sim.cxx
// Grant McNamara, June 2020
// From Isaac Mooney, July 2019

#include "params.hh"
#include "funcs.hh"
#include "TStarJetPicoDefinitions.h"



using namespace fastjet;
using namespace std;
using namespace Analysis;



int main(int argc, const char **argv){
	
	TH1::SetDefaultSumw2( );
	TH2::SetDefaultSumw2( );
	TH3::SetDefaultSumw2( );
	
	
	// Default command line arguments
	std::string outputDir = "out/";
	std::string outFileName = "test.root";
	std::string chainList = "simList.txt";
	double radius = 0.4;
	bool full = 1;
	bool match = 0;
	
	
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
		if (arguments[3] == "ch") {full = 0;} else {full = 1;}
		if (arguments[4] == "matched") {match = 1;} else if (arguments[4] == "unmatched") {match = 0;} else {cerr << "Not a valid flag!" << endl; exit(1);}
		chainList         = arguments[5];
	
		//Printing settings:
		cout << "Outputting to: " << (outputDir+outFileName).c_str() << "\nSettings:\nR = " << radius << " " <<  arguments[3] << " jets;\n match pythia and geant? " << match << ";\n input file: " << chainList << "\n";
		break;
	}
	default: { // Error: invalid custom settings
		__ERR("Invalid number of command line arguments");
		return -1;
		break;
	}
	}
	


	//  Initialize readers and provide chains
	TStarJetPicoReader* P6Reader = new TStarJetPicoReader();
	TChain* P6Chain = new TChain( "JetTreeMc" );
	TStarJetPicoReader* GEANTReader = new TStarJetPicoReader();
	TChain* GEANTChain = new TChain( "JetTree" );


	// Check to see if the input is a .root file or a .txt
	bool inputIsRoot = Analysis::HasEnding( chainList.c_str(), ".root" );
	bool inputIsTxt  = Analysis::HasEnding( chainList.c_str(), ".txt"  );
	bool inputIsList = Analysis::HasEnding( chainList.c_str(), ".list" );
  
	// If its a recognized file type, build the chain
	// If its not recognized, exit
	if ( inputIsRoot ) { P6Chain->Add( chainList.c_str()); GEANTChain->Add( chainList.c_str());}
	else if ( inputIsTxt )  { P6Chain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str()); GEANTChain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
	else if ( inputIsList)  { P6Chain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str()); GEANTChain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
	else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }


	// file with histograms for detector/generator level pT smearing
	// contains detector resolution histograms (meaning?)
	TFile *p6match = new TFile(("~/ppJCandBF/out/sim/hists/p6_hists_det_resolut_R" + (string) argv[3] + ".root").c_str(),"READ");
	//					out/sim/hists/p6_hists_det_resolut_R04.root 

	TH2D *pt_res_py2D = (TH2D*) p6match->Get("deltaPtvPyPt");
	TH2D *pt_res_ge2D = (TH2D*) p6match->Get("deltaPtvGePt");
	
	pt_res_py2D->SetDirectory(0);
	pt_res_ge2D->SetDirectory(0);
	p6match->Close();
	//


	TProfile *pt_res_py = (TProfile*) pt_res_py2D->ProfileX("pt_res_py", 1, 220);
	TProfile *pt_res_ge = (TProfile*) pt_res_ge2D->ProfileX("pt_res_ge", 1, 220);
	

	//define relevant data structures
	TString pythiaFilename;
	TStarJetPicoEventHeader* p_header;
	TStarJetPicoEvent* p_event;
	TStarJetVectorContainer<TStarJetVector>* p_container;
	TStarJetVector* p_sv;


	TString geantFilename;
	TStarJetPicoEventHeader* g_header;
	TStarJetPicoEvent* g_event;
	TStarJetVectorContainer<TStarJetVector>* g_container;
	TStarJetVector* g_sv;

	

	//defining local containers to be linked to tree branches
	int p_EventID;
	double p_n_jets, p_wt;

	// either a vector of doubles or a vector of vectors of doubles...
	//vector<double> part_Pt; vector<double> part_Eta; vector<double> part_Phi; //vector<double> p_consM; vector<double> p_E;
	vector<vector<double> > part_Pt; vector<vector<double> > part_Eta; vector<vector<double> > part_Phi; //vector<double> p_consM; vector<double> p_E;

	vector<double> p_jetMult;
	vector<double> p_Pt; vector<double> p_Eta; vector<double> p_Phi; vector<double> p_M; vector<double> p_E;
	vector<vector<double> > p_consPt; vector<vector<double> > p_consEta; vector<vector<double> > p_consPhi; //vector<double> p_consM; vector<double> p_E;
	vector<double> p_mcd;
	vector<int> p_n_cons;

	int g_EventID;
	double g_n_jets, g_wt;
	vector<double> g_jetMult;
	vector<double> g_Pt; vector<double> g_Eta; vector<double> g_Phi; vector<double> g_M; vector<double> g_E;
	vector<vector<double> > g_consPt; vector<vector<double> > g_consEta; vector<vector<double> > g_consPhi; //vector<double> g_M; vector<double> g_E;
	vector<double> g_mcd;
	


	//tree to hold jet and constituent quantites
	TTree *eventTree = new TTree("event","event");
	eventTree->Branch("p_EventID", &p_EventID);
	eventTree->Branch("p_n_jets", &p_n_jets);
	eventTree->Branch("part_Pt", &part_Pt); eventTree->Branch("part_Eta",&part_Eta); eventTree->Branch("part_Phi",&part_Phi);
//	eventTree->Branch("p_jetMult",&p_jetMult);
	eventTree->Branch("p_Pt", &p_Pt); eventTree->Branch("p_Eta",&p_Eta); eventTree->Branch("p_Phi",&p_Phi);
	eventTree->Branch("p_n_cons", &p_n_cons);
	eventTree->Branch("p_consPt", &p_consPt); eventTree->Branch("p_consEta",&p_consEta); eventTree->Branch("p_consPhi",&p_consPhi);
//	eventTree->Branch("p_mcd",&p_mcd);
	eventTree->Branch("p_weight", &p_wt);

	eventTree->Branch("g_EventID", &g_EventID);
	eventTree->Branch("g_n_jets", &g_n_jets);
//	eventTree->Branch("g_jetMult",&g_jetMult);
	eventTree->Branch("g_Pt", &g_Pt); eventTree->Branch("g_Eta",&g_Eta); eventTree->Branch("g_Phi",&g_Phi);
	eventTree->Branch("g_consPt", &g_consPt); eventTree->Branch("g_consEta",&g_consEta); eventTree->Branch("g_consPhi",&g_consPhi);
//	eventTree->Branch("g_mcd",&g_mcd);
	eventTree->Branch("g_weight", &g_wt);

	// for jet charge need the following:
	// for both particle and detector level
	// jet pt, constituent charge, constituent pt
	// will match jets before filling trees
	//



	// JET CHARGE
	vector<double> part_jetpt; vector<vector<double> > part_conspt; vector<vector<double> > part_chcons;
	vector<double>  det_jetpt; vector<vector<double> >  det_conspt; vector<vector<double> >  det_chcons;
	vector<double> miss_jetpt; vector<vector<double> > miss_conspt; vector<vector<double> > miss_chcons;
	vector<double> fake_jetpt; vector<vector<double> > fake_conspt; vector<vector<double> > fake_chcons;
	
	vector<double> part_jeteta; vector<double> part_jetphi; vector<double> part_jetmass;
	vector<double> det_jeteta; vector<double> det_jetphi; vector<double> det_jetmass;

	vector<double> miss_jeteta; vector<double> miss_jetphi; vector<double> miss_jetmass;
	vector<double> fake_jeteta; vector<double> fake_jetphi; vector<double> fake_jetmass;

	vector<double> part_jetNEF; vector<double> det_jetNEF; vector<double> miss_jetNEF; vector<double> fake_jetNEF;


	TTree *jetChargeTree = new TTree("jetChargeTree", "jetChargeTree");
	jetChargeTree->Branch("p_EventID", &p_EventID);
	jetChargeTree->Branch("p_n_jets", &p_n_jets);
	jetChargeTree->Branch("g_EventID", &g_EventID);
	jetChargeTree->Branch("g_n_jets", &g_n_jets);

	jetChargeTree->Branch("p_n_cons", &p_n_cons);

	jetChargeTree->Branch("part_jetpt", &part_jetpt);
	jetChargeTree->Branch("det_jetpt", &det_jetpt);

	jetChargeTree->Branch("part_jeteta", &part_jeteta);
	jetChargeTree->Branch("part_jetphi", &part_jetphi);
	jetChargeTree->Branch("det_jeteta", &det_jeteta);
	jetChargeTree->Branch("det_jetphi", &det_jetphi);

	jetChargeTree->Branch("part_jetmass", &part_jetmass);
	jetChargeTree->Branch("det_jetmass", &det_jetmass);


	jetChargeTree->Branch("part_conspt", &part_conspt);
	jetChargeTree->Branch("det_conspt", &det_conspt);
	jetChargeTree->Branch("part_chcons", &part_chcons);
	jetChargeTree->Branch("det_chcons", &det_chcons);
	jetChargeTree->Branch("weight", &p_wt);

	jetChargeTree->Branch("miss_jetpt", &miss_jetpt);
	jetChargeTree->Branch("miss_conspt", &miss_conspt);
	jetChargeTree->Branch("miss_chcons", &miss_chcons);

	jetChargeTree->Branch("miss_jeteta", &miss_jeteta);
	jetChargeTree->Branch("miss_jetphi", &miss_jetphi);

	jetChargeTree->Branch("miss_jetmass", &miss_jetmass);


	jetChargeTree->Branch("fake_jetpt", &fake_jetpt);
	jetChargeTree->Branch("fake_conspt", &fake_conspt);
	jetChargeTree->Branch("fake_chcons", &fake_chcons);

	jetChargeTree->Branch("fake_jeteta", &fake_jeteta);
	jetChargeTree->Branch("fake_jetphi", &fake_jetphi);

	jetChargeTree->Branch("fake_jetmass", &fake_jetmass);


	jetChargeTree->Branch("part_jetNEF", &part_jetNEF);
	jetChargeTree->Branch("det_jetNEF", &det_jetNEF);
	jetChargeTree->Branch("miss_jetNEF", &miss_jetNEF);
	jetChargeTree->Branch("fake_jetNEF", &fake_jetNEF);


    
        TTree *towerscale = new TTree("towerscale", "towerscale");

        towerscale->Branch("part_jetpt", &part_jetpt);
        towerscale->Branch("det_jetpt", &det_jetpt);
        towerscale->Branch("part_conspt", &part_conspt);
        towerscale->Branch("det_conspt", &det_conspt);
        towerscale->Branch("part_chcons", &part_chcons);
        towerscale->Branch("det_chcons", &det_chcons);
        towerscale->Branch("weight", &p_wt);

        towerscale->Branch("miss_jetpt", &miss_jetpt);
        towerscale->Branch("miss_conspt", &miss_conspt);
        towerscale->Branch("miss_chcons", &miss_chcons);

        towerscale->Branch("fake_jetpt", &fake_jetpt);
        towerscale->Branch("fake_conspt", &fake_conspt);
        towerscale->Branch("fake_chcons", &fake_chcons);


        TTree *trackingeff = new TTree("trackingeff", "trackingeff");

        trackingeff->Branch("part_jetpt", &part_jetpt);
        trackingeff->Branch("det_jetpt", &det_jetpt);
        trackingeff->Branch("part_conspt", &part_conspt);
        trackingeff->Branch("det_conspt", &det_conspt);
        trackingeff->Branch("part_chcons", &part_chcons);
        trackingeff->Branch("det_chcons", &det_chcons);
        trackingeff->Branch("weight", &p_wt);

        trackingeff->Branch("miss_jetpt", &miss_jetpt);
        trackingeff->Branch("miss_conspt", &miss_conspt);
        trackingeff->Branch("miss_chcons", &miss_chcons);

        trackingeff->Branch("fake_jetpt", &fake_jetpt);
        trackingeff->Branch("fake_conspt", &fake_conspt);
        trackingeff->Branch("fake_chcons", &fake_chcons);
        
        
        
        TTree *hadroncorr50 = new TTree("hadroncorr50", "hadroncorr50");
        
        hadroncorr50->Branch("part_jetpt", &part_jetpt);
        hadroncorr50->Branch("det_jetpt", &det_jetpt);
        hadroncorr50->Branch("part_conspt", &part_conspt);
        hadroncorr50->Branch("det_conspt", &det_conspt);
        hadroncorr50->Branch("part_chcons", &part_chcons);
        hadroncorr50->Branch("det_chcons", &det_chcons);
        hadroncorr50->Branch("weight", &p_wt);

        hadroncorr50->Branch("miss_jetpt", &miss_jetpt);
        hadroncorr50->Branch("miss_conspt", &miss_conspt);
        hadroncorr50->Branch("miss_chcons", &miss_chcons);

        hadroncorr50->Branch("fake_jetpt", &fake_jetpt);
        hadroncorr50->Branch("fake_conspt", &fake_conspt);
        hadroncorr50->Branch("fake_chcons", &fake_chcons);
        
        
        TTree *detsmear = new TTree("detsmear", "detsmear");

        detsmear->Branch("part_jetpt", &part_jetpt);
        detsmear->Branch("det_jetpt", &det_jetpt);
        detsmear->Branch("part_conspt", &part_conspt);
        detsmear->Branch("det_conspt", &det_conspt);
        detsmear->Branch("part_chcons", &part_chcons);
        detsmear->Branch("det_chcons", &det_chcons);
        detsmear->Branch("weight", &p_wt);

        detsmear->Branch("miss_jetpt", &miss_jetpt);
        detsmear->Branch("miss_conspt", &miss_conspt);
        detsmear->Branch("miss_chcons", &miss_chcons);

        detsmear->Branch("fake_jetpt", &fake_jetpt);
        detsmear->Branch("fake_conspt", &fake_conspt);
        detsmear->Branch("fake_chcons", &fake_chcons);
        
        
        TTree *gensmear = new TTree("gensmear", "gensmear");

        gensmear->Branch("part_jetpt", &part_jetpt);
        gensmear->Branch("det_jetpt", &det_jetpt);
        gensmear->Branch("part_conspt", &part_conspt);
        gensmear->Branch("det_conspt", &det_conspt);
        gensmear->Branch("part_chcons", &part_chcons);
        gensmear->Branch("det_chcons", &det_chcons);
        gensmear->Branch("weight", &p_wt);

        gensmear->Branch("miss_jetpt", &miss_jetpt);
        gensmear->Branch("miss_conspt", &miss_conspt);
        gensmear->Branch("miss_chcons", &miss_chcons);

        gensmear->Branch("fake_jetpt", &fake_jetpt);
        gensmear->Branch("fake_conspt", &fake_conspt);
        gensmear->Branch("fake_chcons", &fake_chcons);


	const int nJetBins = 4;
	const double jetEdges[nJetBins+1] = {10.0, 15.0, 20.0, 30.0, 40.0};
	// for the scaling of errors (?) fill histogram with all particle level jets
	TH1D* jetpt_match_plus_miss = new TH1D("jetpt_match_plus_miss", "", nJetBins, jetEdges);
	// figure out binning of this histogram... probably according to the 10-15, 15-20, 20-30, 30-40 bins i would have to assume


	TH1D* jetpt_match = new TH1D("jetpt_match", "", 7, 5.0, 40.0);
	TH1D* jetpt_miss = new TH1D("jetpt_miss", "", 7, 5.0, 40.0);
	TH1D* jetpt_fake = new TH1D("jetpt_fake", "", 7, 5.0, 40.0);

/*
	// BALANCE FUNCTION
	// fill tree with constituents' pt, charge
	// think about structure of root macro to calculate R2, A, B2
	vector<double> part_jetpt_bf; vector<vector<double> > part_bf_conspt; vector<vector<double> > part_bf_chcons;
	vector<double>  det_jetpt_bf; vector<vector<double> >  det_bf_conspt; vector<vector<double> >  det_bf_chcons;

	TTree *bfTree = new TTree("bfTree", "bf");
//	bfTree->Branch("event", &event, "event"); // event number: used for checking if constituents a and b are in same jet
	bfTree->Branch("part_jetpt_bf", &part_jetpt_bf); // jet pt: used for checking if constituents a and b are in same jet
	bfTree->Branch("det_jetpt_bf", &det_jetpt_bf); // jet pt: used for checking if constituents a and b are in same jet
	bfTree->Branch("part_bf_conspt", &part_bf_conspt);
	bfTree->Branch("det_bf_conspt", &det_bf_conspt);
	bfTree->Branch("part_bf_chcons", &part_bf_chcons);
	bfTree->Branch("det_bf_chcons", &det_bf_chcons);
	bfTree->Branch("weight", &p_wt);
*/


	//defining the algorithm and radius parameter for clustering jets
	JetDefinition jet_def(antikt_algorithm, radius/*R*/);


	//SELECTORS
	// Constituent selectors
	// ---------------------
	Selector select_track_rap = fastjet::SelectorAbsEtaMax(max_track_rap);
	Selector select_lopt      = fastjet::SelectorPtMin( partMinPt );
	Selector select_loptmax   = fastjet::SelectorPtMax( partMaxPt );
	Selector spart = select_track_rap * select_lopt * select_loptmax;
	
	// Jet candidate selectors
	// -----------------------
	Selector select_jet_rap     = fastjet::SelectorAbsEtaMax(max_rap); // Cut on eta is correct, cut on rapidity is wrong and legacy--> change variable names to eta here and in src/params.hh
	Selector select_det_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
	Selector select_gen_jet_pt_min = fastjet::SelectorPtMin( jet_ptmin );
	Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
	Selector select_det_jet_m_min = fastjet::SelectorMassMin( 1.0/*0.0*/ /*mass_min*/ ); // for consistency with Isaac: seems that he was still applying this detector level cut, after applying this my particle level matched (to detector level) jets are consistent but when not requiring jets to be matched there remains inconsistency, specifically in particle level jet pt spectrum
	Selector select_gen_jet_m_min = fastjet::SelectorMassMin( 0.0 );
	
	Selector sjet_gen = select_jet_rap && select_gen_jet_pt_min && select_jet_pt_max /*&& select_gen_jet_m_min*/;
	Selector sjet_det = select_jet_rap && select_det_jet_pt_min && select_jet_pt_max && select_det_jet_m_min;
	
	vector<PseudoJet> p_Particles, g_Particles, p_Jets, g_Jets, g_Jets_Initial;
	int p_n_accepted = 0; int g_n_accepted = 0; int p_NJets = 0; int g_NJets = 0;
	int counter_debug = 0;
	double mc_weight = -1;

        double hc = 0.9999;
    
        // loop over systematics
        // start with nominal and tower scaling
        int nSources = 6;

        int njet_events[nSources] = {0, 0, 0, 0, 0, 0};
        for(int iSyst = 0; iSyst < nSources; iSyst++){
            if(iSyst == 0){cout << "RUNNING WITH NOMINAL SETTINGS\n";}
            if(iSyst == 1){cout << "RUNNING WITH TOWER SCALE SETTINGS\n";}
            if(iSyst == 2){cout << "RUNNING WITH TRACKING EFFICIENCY SETTINGS\n";}
            if(iSyst == 3){cout << "RUNNING WITH HADRONIC CORRECTION 50%\n";}
            if(iSyst == 4){cout << "RUNNING WITH DETECTOR SMEARING\n";}
            if(iSyst == 5){cout << "RUNNING WITH GENERATOR SMEARING\n";}
	    
	    
            hc = 0.9999;
            if(iSyst == 3){double hc = 0.5;}
        
    
            //initialize both readers
	    InitReader(P6Reader, P6Chain, nEvents, "All", truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, hc, false, sim_badTowers, sim_bad_run_list);
	    InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, hc, false, det_badTowers, dat_bad_run_list);

            // loop over events to find low statistics jet pt vs constituent pt event numbers
            for(int event = 0; event < P6Chain->GetEntries(); event++){
		P6Reader->ReadEvent(event);
		GEANTReader->ReadEvent(event);

		p_EventID = P6Reader->GetNOfCurrentEvent();
		g_EventID = GEANTReader->GetNOfCurrentEvent();


		//clearing vectors; initializing variables to -9999
		mc_weight = -9999;
		p_n_jets = -9999; p_wt = -9999;
		part_Pt.clear(); part_Eta.clear(); part_Phi.clear();
		p_jetMult.clear();
		p_Pt.clear(); p_Eta.clear(); p_Phi.clear(); p_M.clear(); p_E.clear();
		p_consPt.clear(); p_consEta.clear(); p_consPhi.clear();

		g_n_jets = -9999; g_wt = -9999;
		g_jetMult.clear();
		g_Pt.clear(); g_Eta.clear(); g_Phi.clear(); g_M.clear(); g_E.clear();
		g_consPt.clear(); g_consEta.clear(); g_consPhi.clear();


		// clear vectors that go into trees every event
		part_jetpt.clear(); part_conspt.clear(); part_chcons.clear();
		det_jetpt.clear(); det_conspt.clear(); det_chcons.clear();
		miss_jetpt.clear(); miss_conspt.clear(); miss_chcons.clear();
		fake_jetpt.clear(); fake_conspt.clear(); fake_chcons.clear();
		// vectors that need to be cleared... had not been clearing misses and fakes here...
		// fixed 4/17/22

		part_jeteta.clear(); part_jetphi.clear(); part_jetmass.clear();
		det_jeteta.clear(); det_jetphi.clear(); det_jetmass.clear();

		miss_jeteta.clear(); miss_jetphi.clear(); miss_jetmass.clear();
		fake_jeteta.clear(); fake_jetphi.clear(); fake_jetmass.clear();

		part_jetNEF.clear();
		det_jetNEF.clear();
		miss_jetNEF.clear();
		fake_jetNEF.clear();

		p_n_cons.clear();


//		part_jetpt_bf.clear(); part_bf_conspt.clear(); part_bf_chcons.clear();
//		det_jetpt_bf.clear(); det_bf_conspt.clear(); det_bf_chcons.clear();


		p_Particles.clear(); g_Particles.clear();
		p_Jets.clear(); g_Jets.clear(), g_Jets_Initial.clear();

		if ( match && GEANTReader->ReadEvent(p_EventID) != 1 ) {
			/*cout << "no corresponding geant event...skipping event " << p_EventID << endl;*/
			continue;//goes to the next event
		}
		if ( match && (p_EventID != g_EventID) ) { cerr << "ERROR: READING DIFFERENT EVENTS. EXITING." << endl; exit(1);}


		//filling the data structures that were defined before the event loop:
		p_event = P6Reader->GetEvent();
		p_header = p_event->GetHeader();
		g_event = GEANTReader->GetEvent();
		g_header = g_event->GetHeader();
		
		p_container = P6Reader->GetOutputContainer();
		g_container = GEANTReader->GetOutputContainer();
		
		pythiaFilename =  P6Reader->GetInputChain()->GetCurrentFile()->GetName();
		geantFilename =  GEANTReader->GetInputChain()->GetCurrentFile()->GetName();

		if (match && (pythiaFilename != geantFilename)) {std::cerr << "FILES DON'T MATCH! EXITING." << std::endl; exit(1);}
		
		p_wt = LookupRun12Xsec( pythiaFilename );
		g_wt = LookupRun12Xsec( geantFilename );
		
		// added if statement on 4/20/22
		if (match && (p_wt != g_wt)) {std::cerr << "WEIGHTS DON'T MATCH! EXITING." << std::endl; exit(1);}
		mc_weight = p_wt;

		// from Isaac
	        if (iSyst == 1) {//varying the gain of the towers
        	    for (int i = 0; i < g_container->GetEntries(); ++ i) {
         	        g_sv = g_container->Get(i);
               		if (!(g_sv->IsCharged())) {
                            //cout << "DEBUG: Pre-change: " << g_sv->E() << " " << g_sv->Eta() << " " << g_sv->Phi() << " " << g_sv->M() << endl;
                            double Enew = 1.038*g_sv->E();
                            //g_sv is a shallow copy of g_container->Get(i) so editing it also edits the container, yay
                            g_sv->SetE(Enew);
                            //g_sv->SetPtEtaPhiM(sqrt(Etnew*Etnew - g_sv->M()*g_sv->M()), g_sv->Eta(), g_sv->Phi(), g_sv->M());
                            //cout << "DEBUG: Post-change: " << g_sv->E() << " " << g_sv->Eta() << " " << g_sv->Phi() << " " << g_sv->M() << endl;
                            //cout << "DEBUG: Percent diff: " << ((g_sv->E()/(double)(Enew/(double)1.038)) - 1)*100 << "%" << endl;
                        }
                    }
                }

/////////////////////////// APPLY SELECTORS
		GatherParticles(p_container, p_sv, p_Particles, full, 1); //Pythia; full = 0 => charged-only, 1 => ch+ne
		GatherParticles(g_container, g_sv, g_Particles, full, 0); //Pythia; full = 0 => charged-only, 1 => ch+ne

		// from Isaac
                if (iSyst == 2) {//varying the tracking efficiency randomly by 4%
                    double effic_num;
                    for (int i = 0; i < g_Particles.size(); ++ i) {
                        if (g_Particles[i].user_index() != 0) {
                            effic_num = gRandom->Uniform(0.0, 1.0);
                            if (effic_num > 0.96) {
                                g_Particles.erase(g_Particles.begin() + i);
                                i --; //need to account for the shrinking of the list.
                            }
                        }
                    }
                }

		vector<PseudoJet> p_cut_Particles = spart(p_Particles);
		vector<PseudoJet> g_cut_Particles = spart(g_Particles);

		vector<double> particle_cut_pt; vector<double> particle_cut_eta; vector<double> particle_cut_phi;
		for(int ik = 0; ik < p_cut_Particles.size(); ik++){

			particle_cut_pt.push_back(p_cut_Particles[ik].pt());
			particle_cut_eta.push_back(p_cut_Particles[ik].eta());
			particle_cut_phi.push_back(p_cut_Particles[ik].phi());
		}

		part_Pt.push_back(particle_cut_pt);
		part_Eta.push_back(particle_cut_eta);
		part_Phi.push_back(particle_cut_phi);


/////////////////////////// FIND JETS
		ClusterSequence p_Cluster(p_cut_Particles, jet_def);
		ClusterSequence g_Cluster(g_cut_Particles, jet_def);
		
		p_Jets = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets()));
		g_Jets_Initial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));

//		vector<PseudoJet> g_Jets;
		ApplyNEFSelection(g_Jets_Initial, g_Jets);

		p_n_jets = p_Jets.size();
		g_n_jets = g_Jets.size();

		if (DiscardEvent(pythiaFilename, p_Jets, g_Jets)) { counter_debug ++; cout << "event ID = " << p_EventID << "\n"; continue; }

/////////////////////////// MATCH JETS
		if(match){
			if(iSyst == 0){
				// as of 4/27/22 I only have one match_plus_miss histogram, only fill for nominal settings for now
				for(int kj = 0; kj < p_Jets.size(); kj++){
					jetpt_match_plus_miss->Fill(p_Jets[kj].pt(), mc_weight);

					vector<double> cons_pt; vector<double> cons_eta; vector<double> cons_phi;

					for(int ki = 0; ki < p_Jets[kj].constituents().size(); ki++){
						cons_pt.push_back(p_Jets[kj].constituents()[ki].pt());
						cons_eta.push_back(p_Jets[kj].constituents()[ki].eta());
						cons_phi.push_back(p_Jets[kj].constituents()[ki].phi());
					}
					p_Pt.push_back(p_Jets[kj].pt());
					p_Eta.push_back(p_Jets[kj].eta());
					p_Phi.push_back(p_Jets[kj].phi());

					p_consPt.push_back(cons_pt);
					p_consEta.push_back(cons_eta);
					p_consPhi.push_back(cons_phi);
				}
				for(int kj = 0; kj < g_Jets.size(); kj++){
					vector<double> cons_pt; vector<double> cons_eta; vector<double> cons_phi;

					for(int ki = 0; ki < g_Jets[kj].constituents().size(); ki++){
						cons_pt.push_back(g_Jets[kj].constituents()[ki].pt());
						cons_eta.push_back(g_Jets[kj].constituents()[ki].eta());
						cons_phi.push_back(g_Jets[kj].constituents()[ki].phi());
					}
					g_Pt.push_back(g_Jets[kj].pt());
					g_Eta.push_back(g_Jets[kj].eta());
					g_Phi.push_back(g_Jets[kj].phi());

					g_consPt.push_back(cons_pt);
					g_consEta.push_back(cons_eta);
					g_consPhi.push_back(cons_phi);
				}
			}

			std::vector<fastjet::PseudoJet> g_matches; std::vector<fastjet::PseudoJet> p_matches;
			std::vector<fastjet::PseudoJet> g_matches_for_fakes; std::vector<fastjet::PseudoJet> p_matches_for_fakes;
			std::vector<fastjet::PseudoJet> misses; // only appear in particle level (detector does not find the jet)
			std::vector<fastjet::PseudoJet> fakes; // only appear in detector level (detector jet that does not originate from 'real' jet)
			std::vector<int> match_indices;
			std::vector<int> miss_indices;
			std::vector<int> fake_indices;


			if (p_Jets.size() != 0) {
				g_matches.clear(); p_matches.clear(); misses.clear();
				match_indices.clear(); miss_indices.clear();
				
				match_indices = MatchJets(g_Jets, p_Jets, g_matches, p_matches); //find matches
				
				if (g_matches.size() != p_matches.size()) {std::cerr << "Somehow we have different-sized match vectors. Exiting!" <<std::endl; exit(1);}
				
				if (g_matches.size() < p_Jets.size()) { //then we have misses
					miss_indices = FakesandMisses(p_matches, p_Jets, misses); //find misses
				}
			}
			if(g_Jets.size() != 0){
				g_matches_for_fakes.clear(); p_matches_for_fakes.clear();
				fakes.clear(); fake_indices.clear();
				
				MatchJets(p_Jets, g_Jets, p_matches_for_fakes, g_matches_for_fakes);
				
				if(g_matches_for_fakes.size() != p_matches_for_fakes.size()){std::cerr << "Somehow we have different-sized match vectors. Exiting!" <<std::endl; exit(1);}
				
				if(p_matches_for_fakes.size() < g_Jets.size()){
					fake_indices = FakesandMisses(g_matches_for_fakes, g_Jets, fakes);
				}
			}

			// for constituent quantities to add to trees
			// vector<vector<double> > in trees, clear the vector<double> each jet and add to the outer vector each jet
			vector<double> part_jc_ch; vector<double> part_jc_pt; vector<double> det_jc_ch; vector<double> det_jc_pt;
			
			double prior_adjust = 0;
			
	                if(p_matches.size() != 0){njet_events[iSyst]++;}
	    		for(int i = 0; i < p_matches.size(); i++){
				double NEF_p = 0; double NEF_g = 0;
				double part_pt = p_matches[i].pt();
				double det_pt = g_matches[i].pt();

// /*
				if(iSyst == 4){ // detector level pT smearing
					double res_for_this_jet = pt_res_ge->GetBinContent( pt_res_ge->GetXaxis()->FindBin( det_pt ) );
					prior_adjust = fabs( gRandom->Gaus( 0, fabs( res_for_this_jet * det_pt ) ) );
					det_pt -= prior_adjust;
				}
				if(iSyst == 5){ // generator level pT smearing
					double res_for_this_jet = pt_res_py->GetBinContent( pt_res_py->GetXaxis()->FindBin( part_pt ) );
					prior_adjust = fabs( gRandom->Gaus( 0, fabs( res_for_this_jet * part_pt ) ) );
					part_pt -= prior_adjust;
				}
// */

				// clear vector of constituent quantities for each jet
				part_jc_ch.clear(); part_jc_pt.clear(); det_jc_ch.clear(); det_jc_pt.clear();

				p_n_cons.push_back(p_matches[i].constituents().size());

				for(int j = 0; j < max( g_matches[i].constituents().size(), p_matches[i].constituents().size() ); j++){
					if(j < p_matches[i].constituents().size()){
						double ch_cons = p_matches[i].constituents()[j].user_index();
						double pt_cons = p_matches[i].constituents()[j].pt();

						part_jc_ch.push_back(ch_cons);
						part_jc_pt.push_back(pt_cons);

						if(ch_cons == 0){
							NEF_p += pt_cons/(double) part_pt;
						}
					}
					if(j < g_matches[i].constituents().size()){
						double ch_cons = g_matches[i].constituents()[j].user_index();
						double pt_cons = g_matches[i].constituents()[j].pt();

						det_jc_ch.push_back(ch_cons);
						det_jc_pt.push_back(pt_cons);

						if(ch_cons == 0){
							NEF_g += pt_cons/(double) det_pt;
						}
					}
				}
				// add vector of constituent pt/charge to vector that goes into tree

				// add jet pt to vector that will go into tree
				part_jetpt.push_back(part_pt);
				det_jetpt.push_back(det_pt);
				
				part_jeteta.push_back(p_matches[i].eta());
				part_jetphi.push_back(p_matches[i].phi());
				det_jeteta.push_back(g_matches[i].eta());
				det_jetphi.push_back(g_matches[i].phi());

				if(iSyst == 0){
					part_jetmass.push_back(p_matches[i].m());
					det_jetmass.push_back(g_matches[i].m());

					part_jetNEF.push_back(NEF_p);
					det_jetNEF.push_back(NEF_g);
				}

				part_conspt.push_back(part_jc_pt);
				det_conspt.push_back(det_jc_pt);
				part_chcons.push_back(part_jc_ch);
				det_chcons.push_back(det_jc_ch);


//				part_jetpt_bf.push_back(part_pt);
//				det_jetpt_bf.push_back(det_pt);

//				part_bf_conspt.push_back(part_jc_pt);
//				det_bf_conspt.push_back(det_jc_pt);
//				part_bf_chcons.push_back(part_jc_ch);
//				det_bf_chcons.push_back(det_jc_ch);

			} // i loop over matched jets

			vector<double> miss_jc_ch; vector<double> miss_jc_pt;
			vector<double> fake_jc_ch; vector<double> fake_jc_pt;

			for(int i = 0; i < misses.size(); i++){
				double NEF_m = 0;
				miss_jc_ch.clear(); miss_jc_pt.clear();

				double m_pt = misses[i].pt();
// /*
				if(iSyst == 5){
					double res_for_this_jet = pt_res_py->GetBinContent( pt_res_py->GetXaxis()->FindBin( m_pt ) );
					prior_adjust = fabs( gRandom->Gaus( 0, fabs( res_for_this_jet * m_pt ) ) );
					m_pt -= prior_adjust;
				}
// */

				for(int j = 0; j < misses[i].constituents().size(); j++){
					double pt_cons = misses[i].constituents()[j].pt();
					double ch_cons = misses[i].constituents()[j].user_index();

					miss_jc_ch.push_back(ch_cons);
					miss_jc_pt.push_back(pt_cons);

					if(ch_cons == 0){
						NEF_m += pt_cons/(double) m_pt;
					}
				}

				miss_jetpt.push_back(m_pt);
				miss_conspt.push_back(miss_jc_pt);
				miss_chcons.push_back(miss_jc_ch);

				miss_jeteta.push_back(misses[i].eta());
				miss_jetphi.push_back(misses[i].phi());
				miss_jetmass.push_back(misses[i].m());

				miss_jetNEF.push_back(NEF_m);

			} // loop over missed jets
			for(int i = 0; i < fakes.size(); i++){
				double NEF_f = 0;
				
				fake_jc_ch.clear(); fake_jc_pt.clear();

				double f_pt = fakes[i].pt();
				
// /*
				if(iSyst == 4){
					double res_for_this_jet = pt_res_ge->GetBinContent( pt_res_ge->GetXaxis()->FindBin( f_pt ) );
					prior_adjust = fabs( gRandom->Gaus( 0, fabs( res_for_this_jet * f_pt ) ) );
					f_pt -= prior_adjust;
				}
// */

				for(int j = 0; j < fakes[i].constituents().size(); j++){
					double pt_cons = fakes[i].constituents()[j].pt();
					double ch_cons = fakes[i].constituents()[j].user_index();

					fake_jc_ch.push_back(ch_cons);
					fake_jc_pt.push_back(pt_cons);

					if(ch_cons == 0){
						NEF_f += pt_cons/(double) f_pt;
					}
				}

				fake_jetpt.push_back(f_pt);
				fake_conspt.push_back(fake_jc_pt);
				fake_chcons.push_back(fake_jc_ch);

				fake_jeteta.push_back(fakes[i].eta());
				fake_jetphi.push_back(fakes[i].phi());
				fake_jetmass.push_back(fakes[i].m());

				miss_jetNEF.push_back(NEF_f);

			} // loop over fake jets
			if(iSyst == 0){
				for(int i_match = 0; i_match < part_jetpt.size(); i_match++){
					jetpt_match->Fill(part_jetpt[i_match], mc_weight);
				}
				for(int i_miss = 0; i_miss < miss_jetpt.size(); i_miss++){
					jetpt_miss->Fill(miss_jetpt[i_miss], mc_weight);
				}
				for(int i_fake = 0; i_fake < fake_jetpt.size(); i_fake++){
					jetpt_fake->Fill(fake_jetpt[i_fake], mc_weight);
				}
			}

            		if(iSyst == 0 && p_Jets.size() != 0){
				eventTree->Fill();
                		jetChargeTree->Fill();
            		}
	                else if(iSyst == 1 && p_Jets.size() != 0){
                		towerscale->Fill();
            		}
	                else if(iSyst == 2 && p_Jets.size() != 0){
                		trackingeff->Fill();
	                }
            		else if(iSyst == 3 && p_Jets.size() != 0){
                		hadroncorr50->Fill();
            		}
            		else if(iSyst == 4 && p_Jets.size() != 0){
                		detsmear->Fill();
            		}
            		else if(iSyst == 5 && p_Jets.size() != 0){
                		gensmear->Fill();
            		}

		///////////////////// fill jetcharge tree and balance function tree or histograms

//		if(p_matches.size() != 0){
//			bfTree->Fill();
//		}


		} // if i require jets be matched
		else{
			vector<double> part_jc_ch; vector<double> part_jc_pt; vector<double> det_jc_ch; vector<double> det_jc_pt;
			for(int i = 0; i < max( p_Jets.size(), g_Jets.size() ); i++){
				if(i < p_Jets.size() ){
					double NEF_p = 0;
					double part_pt = p_Jets[i].pt();
				
					part_jc_ch.clear(); part_jc_pt.clear();

					p_n_cons.push_back(p_Jets[i].constituents().size());

					for(int j = 0; j < p_Jets[i].constituents().size(); j++){
						if(p_Jets[i].constituents()[j].user_index() == -9999){continue;}
						double ch_cons = p_Jets[i].constituents()[j].user_index();
						double pt_cons = p_Jets[i].constituents()[j].pt();

						part_jc_ch.push_back(ch_cons);
						part_jc_pt.push_back(pt_cons);

						if(ch_cons == 0){
							NEF_p += pt_cons/(double) part_pt;
						}
					}

					part_jetpt.push_back(part_pt);
					part_conspt.push_back(part_jc_pt);
					part_chcons.push_back(part_jc_ch);

					part_jeteta.push_back(p_Jets[i].eta());
					part_jetphi.push_back(p_Jets[i].phi());
					part_jetmass.push_back(p_Jets[i].m());
					part_jetNEF.push_back(NEF_p);
				} // ith particle jet

				if(i < g_Jets.size() ) {
					double NEF_g = 0;
					double det_pt = g_Jets[i].pt();
					
					det_jc_ch.clear(); det_jc_pt.clear();
					
					for(int j = 0; j < g_Jets[i].constituents().size(); j++){

						double ch_cons = g_Jets[i].constituents()[j].user_index();
						double pt_cons = g_Jets[i].constituents()[j].pt();
						
						det_jc_ch.push_back(ch_cons);
						det_jc_pt.push_back(pt_cons);
						
						if(ch_cons == 0){
							NEF_g += pt_cons/(double) det_pt;
						}
					}

					det_jetpt.push_back(det_pt);
					det_conspt.push_back(det_jc_pt);
					det_chcons.push_back(det_jc_ch);

					det_jeteta.push_back(g_Jets[i].eta());
					det_jetphi.push_back(g_Jets[i].phi());
					det_jetmass.push_back(g_Jets[i].m());
					det_jetNEF.push_back(NEF_g);
				} // ith detector jet
			} // i loop over jets
            		if(iSyst == 0 && p_Jets.size() != 0){
				eventTree->Fill();
                		jetChargeTree->Fill();
            		}
	                else if(iSyst == 1 && p_Jets.size() != 0){
                		towerscale->Fill();
            		}
	                else if(iSyst == 2 && p_Jets.size() != 0){
                		trackingeff->Fill();
	                }
            		else if(iSyst == 3 && p_Jets.size() != 0){
                		hadroncorr50->Fill();
            		}
            		else if(iSyst == 4 && p_Jets.size() != 0){
                		detsmear->Fill();
            		}
            		else if(iSyst == 5 && p_Jets.size() != 0){
                		gensmear->Fill();
            		}

		} // if no matching is required

	    } // event loop
        } // systematic loop

/*
        cout << "\nnumber of matched jets for nominal settings: " << njet_events[0] << "\n\n";
        cout << "\nnumber of matched jets for tower scale variation: " << njet_events[1] << "\n\n";
        cout << "\nnumber of matched jets for tracking efficiency variation: " << njet_events[2] << "\n\n";
        cout << "\nnumber of matched jets for hadronic correction 50%: " << njet_events[3] << "\n\n";
        cout << "\nnumber of matched jets for detector smearing: " << njet_events[4] << "\n\n";
        cout << "\nnumber of matched jets for generator smearing: " << njet_events[5] << "\n\n";
*/

        TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
        fout->cd();
	
    
	if(match){
	    eventTree->Write();
            jetChargeTree->Write();
            towerscale->Write();
	    trackingeff->Write();
            hadroncorr50->Write();
            detsmear->Write();
            gensmear->Write();

            jetpt_match->Write();
            jetpt_miss->Write();
            jetpt_fake->Write();
       	}
	else{
	    eventTree->Write();
            jetChargeTree->Write();
            towerscale->Write();
	    trackingeff->Write();
            hadroncorr50->Write();
            detsmear->Write();
            gensmear->Write();
	}

/////////////////////////////////////////
	cout << endl << "Writing to:  " << fout->GetName() << endl;


	fout->Write();
	fout->Close();
	
	return 0;
}


