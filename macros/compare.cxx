// compare.cxx
//
// Grant McNamara 5/10/22
// take a root file provided by Isaac with a Tree named Youqi and fill histograms
// to compare to histograms filled by jetChargeTree from my own file to find inconsistencies in kinematics, jet selection, etc
// to obtain consistency with his jet mass analysis so that I can perform unfolding
//
// "miss" tree contains events, jets that isaac finds but i "miss"; "fake" tree contains events, jets that i find but isaac does not ("fake" on my side)


#include <stdio.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>


#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TRandom.h>

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "Headers/funcs.hh"

using namespace std;



double str_to_double (std::string str) {
    std::string Num = str.substr(0,1)+"."+str.substr(1,1); //e.g. 0.4
    double str_double = (double) std::stod(Num); //converting the string 0.4 to a double
    std::cout << "DEBUG: string to number = " << str_double << std::endl;

    return str_double;
}



int main(int argc, const char** argv){
    TH1::SetDefaultSumw2();
    
/*
    if (argc != 6) {
    // "$outLocation $outName $inFiles $inputType $kappa"
      cerr << "Should be five arguments: output location, output name, input name, sim, kappa. Received "
       << argc-1 << ". Exiting." << endl;
      exit(1);
    }
    
    // output files from running sim.cxx
    string fin_name = (string) argv[3];
    TFile *fin = new TFile(fin_name.c_str(),"READ");
    cout << "DEBUG: input file name is " << fin->GetName() << endl;
*/

    // to start, only looking at the one file from pthat 11-15 g4
    // file from both isaac and grant
    TFile *fgrant = new TFile("~/ppJCandBF/out/sim/TEST_Cleanpp12Pico_pt11_15_g4_full_unmatched_R04.root", "READ");
    TFile *fisaac = new TFile("~/ppJCandBF/fromisaac.root", "READ");


    
    TFile *fout = new TFile("~/ppJCandBF/out/compare/isaac_vs_grant.root", "RECREATE");
    cout << "DEBUG: output file name is " << fout->GetName() << endl;
    fout->cd();    

    
    int g_EventID = -1;
//    double g_n_jets = -1;
    vector<double> *part_jetpt = 0;
    vector<double> *part_jeteta = 0;
    vector<double> *part_jetphi = 0;
    vector<double> *part_jetm = 0;
    vector<vector<double> > *part_conspt = 0;
    vector<vector<double> > *part_conseta = 0;
    vector<vector<double> > *part_consphi = 0;
//    vector<vector<double> > *part_chcons = 0;


    TTree *grant = (TTree*) fgrant->Get( "jetChargeTree" );

    grant->SetBranchAddress("p_EventID", &g_EventID);
//    grant->SetBranchAddress("p_n_jets", &g_n_jets);
    grant->SetBranchAddress("part_jetpt", &part_jetpt);
    grant->SetBranchAddress("part_jeteta", &part_jeteta);
    grant->SetBranchAddress("part_jetphi", &part_jetphi);
    grant->SetBranchAddress("part_jetmass", &part_jetm);

    grant->SetBranchAddress("part_conspt", &part_conspt);
    grant->SetBranchAddress("part_conseta", &part_conseta);
    grant->SetBranchAddress("part_consphi", &part_consphi);

//    grant->SetBranchAddress("part_chcons", &part_chcons);


    int i_EventID = -1;
//    double i_n_jets = -1;
    vector<double> *p_pt = 0;
    vector<vector<double> > *p_conspt = 0;
//    vector<vector<double> > *part_chcons = 0;


    TTree *isaac = (TTree*) fisaac->Get( "Youqi" );

    isaac->SetBranchAddress("p_EventID", &i_EventID);
    isaac->SetBranchAddress("p_Pt", &p_pt);
    isaac->SetBranchAddress("p_conspT", &p_conspt);
//    isaac->SetBranchAddress("part_chcons", &part_chcons);



    int m_EventID = -1;
    vector<double> miss_jetpt; vector<double> miss_jeteta; vector<double> miss_jetphi; vector<double> miss_jetm;
    vector<vector<double> > miss_conspt; vector<vector<double> > miss_conseta; vector<vector<double> > miss_consphi;

    TTree *miss = new TTree("miss", "miss");
    miss->Branch("m_EventID", &m_EventID);
    miss->Branch("jet_pt", &miss_jetpt);
    miss->Branch("jet_eta", &miss_jeteta);
    miss->Branch("jet_phi", &miss_jetphi);
    miss->Branch("jet_m", &miss_jetm);

    miss->Branch("cons_pt", &miss_conspt);
    miss->Branch("cons_eta", &miss_conseta);
    miss->Branch("cons_phi", &miss_consphi);



    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hists~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    // without gen-level jet mass cut, only discrepancies on my tree, everything in isaac's tree is in grant's tree
    TH2D* diff_pt_vs_eta = new TH2D("diff_pt_vs_eta", "", 21, 0.0, 20.0, 10, -1.0, 1.0);
    TH2D* diff_pt_vs_phi = new TH2D("diff_pt_vs_phi", "", 21, 0.0, 20.0, 10, 0.0, 6.4);
    TH2D* diff_eta_vs_phi = new TH2D("diff_eta_vs_phi", "", 10, -1.0, 1.0, 10, 0.0, 6.4);

    TH2D* diff_pt_vs_m = new TH2D("diff_pt_vs_m", "", 21, 0.0, 20.0, 21, 0.0, 1.0);
    TH2D* diff_m_vs_eta = new TH2D("diff_m_vs_eta", "", 21, 0.0, 1.0, 10, -1.0, 1.0);
    TH2D* diff_m_vs_phi = new TH2D("diff_m_vs_phi", "", 21, 0.0, 1.0, 10, 0.0, 6.4);


    int one_cons_jets = 0;

    // loop over isaac's tree first, since his tree has more events, jets than mine
    for(int i = 0; i < (int) isaac->GetEntries(); i++){
	miss_jetpt.clear(); miss_jeteta.clear(); miss_jetphi.clear(); miss_jetm.clear();
	miss_conspt.clear(); miss_conseta.clear(); miss_consphi.clear();
        // take each event
        isaac->GetEntry(i); // this will set branches appropriately, so e.g. g_EventID will take the proper integer for the i-th entry

//	cerr << "event id: " << i_EventID <</* "\t" << p_EventID <<*/ "\n";
	if( grant->Scan("p_EventID", Form("p_EventID == %i", i_EventID) ) != 0 ){ // if the second tree has the same event as the first tree
	    for(int j = 0; j < (int) p_pt->size(); j++){

//                double i_pt = p_pt->at(j);
		double i_pt = (*p_pt)[j];

		if( (*p_conspt)[j].size() == 1 ){one_cons_jets++; cerr << "ISAAC HAS ONE CONSTITUENT JET\n";}

		// scan returns the number of instances that pass the condition, so anything > 0 means that event matches, and that the jet pt is equal (as float, use of |pt_1 - pt_2| < epsilon)
		if( grant->Scan("p_EventID:part_jetpt", Form("abs(%f - part_jetpt) <= 0.0001 && p_EventID == %i", i_pt, i_EventID) ) != 0 ){
		// if the second tree with the same event id has a jet with the same pt (eta, phi later maybe)
		    
		    // if jet pt (eta, phi), appears in the given event, ok well let's move on
		    continue;
		}
		else{ // if event is in both trees but jet is only in one tree
		    cerr << "filling miss tree with event " << i_EventID << "\n";
//		    fill vectors with event id, jet pt, eta, phi, number of constituents and constituent pt
		    m_EventID = i_EventID;
		    miss_jetpt.push_back( i_pt );
		    miss_jeteta.push_back( (*p_eta)[j] );
		    miss_jetphi.push_back( (*p_phi)[j] );

		    vector<double> pt_cons; vector<double> eta_cons; vector<double> phi_cons;
		    for(int k = 0; k < (int) (*part_conspt)[j].size(); k++){
			pt_cons.push_back( (*part_conspt)[j][k] );
			eta_cons.push_back( (*part_conseta)[j][k] );
			phi_cons.push_back( (*part_consphi)[j][k] );
		    }
		    miss_conspt.push_back( pt_cons ); // do i even need the constituent information? maybe not
		}
	    }
	}
	else{ // if no matching event in second tree, fill "miss" tree vectors
	    cerr << "discrepancy in isaac's event " << i_EventID << "\n";

//	    fill vectors with event id, jet pt (eta, phi), number of constituents and constituent pt
	    m_EventID = i_EventID;
	    for(int j = 0; j < (int) p_pt->size(); j++){
		miss_jetpt.push_back( (*p_pt)[j] );
		miss_jeteta.push_back( (*p_eta)[j] );
		miss_jetphi.push_back( (*p_phi)[j] );

		vector<double> pt_cons; vector<double> eta_cons; vector<double> phi_cons;
		for(int k = 0; k < (int) (*part_conspt)[j].size(); k++){
		    pt_cons.push_back( (*part_conspt)[j][k] );
		    eta_cons.push_back( (*part_conseta)[j][k] );
		    phi_cons.push_back( (*part_consphi)[j][k] );
		}
		miss_conspt.push_back( pt_cons ); // do i even need the constituent information? maybe not
	    }
	}
	if(miss_jetpt.size() != 0){miss->Fill();}
    }





    // then loop over my tree to check to see if i have events, jets that isaac does not have
    int f_EventID = -1;
    vector<double> fake_jetpt; vector<double> fake_jeteta; vector<double> fake_jetphi; vector<double> fake_jetm;
    vector<vector<double> > fake_conspt; vector<vector<double> > fake_conseta; vector<vector<double> > fake_consphi;
    vector<int> fake_ncons;

    TTree *fake = new TTree("fake", "fake");
    fake->Branch("f_EventID", &f_EventID);
    fake->Branch("jet_pt", &fake_jetpt);
    fake->Branch("jet_eta", &fake_jeteta);
    fake->Branch("jet_phi", &fake_jetphi);
    fake->Branch("jet_m", &fake_jetm);

    fake->Branch("ncons", &fake_ncons);
    fake->Branch("cons_pt", &fake_conspt);
    fake->Branch("cons_eta", &fake_conseta);
    fake->Branch("cons_phi", &fake_consphi);


    for(int i = 0; i < (int) grant->GetEntries(); i++){
	fake_jetpt.clear(); fake_jeteta.clear(); fake_jetphi.clear();
	fake_conspt.clear(); fake_conseta.clear(); fake_consphi.clear();
        // take each event
        grant->GetEntry(i); // this will set branches appropriately, so e.g. g_EventID will take the proper integer for the i-th entry
	
	if( isaac->Scan("p_EventID", Form("p_EventID == %i", g_EventID)) != 0 ){ // if the second tree has the same event as the first tree
            for(int j = 0; j < (int) part_jetpt->size(); j++){

//                double g_pt = part_jetpt->at(j);
		double g_pt = (*part_jetpt)[j];

//		loop over jet pt vectors and compare
		if( isaac->Scan("p_EventID:p_Pt", Form("abs(p_Pt - %f) <= 0.0001 && p_EventID == %i", g_pt, g_EventID) ) != 0 ){
		// if the second tree with the same event id has a jet with the same pt (eta, phi later maybe)
//		    compare n_constituents, if i_ncons == g_ncons, loop over and compare pt of constituents

		    continue;
		}
		else{ // if event is in both trees but jet is only in one tree
		    cerr << "filling fake tree with event " << g_EventID << "\n";
//		    fill vectors with event id, jet pt, eta, phi, number of constituents and constituent pt
		    f_EventID = g_EventID;
		    fake_jetpt.push_back( g_pt );
		    fake_jeteta.push_back( (*part_jeteta)[j] );
		    fake_jetphi.push_back( (*part_jetphi)[j] );
		    fake_jetm.push_back( (*part_jetm)[j] );

		    fake_ncons.push_back( (*part_conspt)[j].size() );

		    cerr << "n constituents? " << (*part_conspt)[j].size() << "\n";
		    cerr << "vector of n constituents? " << (*part_conspt).size() << "\n";


		    vector<double> pt_cons; vector<double> eta_cons; vector<double> phi_cons;
		    for(int k = 0; k < (int) (*part_conspt)[j].size(); k++){
			pt_cons.push_back( (*part_conspt)[j][k] );
			eta_cons.push_back( (*part_conseta)[j][k] );
			phi_cons.push_back( (*part_consphi)[j][k] );
		    }

		    fake_conspt.push_back( pt_cons );
		    fake_conseta.push_back( eta_cons );
		    fake_consphi.push_back( phi_cons );


		    diff_pt_vs_eta->Fill( g_pt, (*part_jeteta)[j] );
		    diff_pt_vs_phi->Fill( g_pt, (*part_jetphi)[j] );
		    diff_eta_vs_phi->Fill( (*part_jeteta)[j], (*part_jetphi)[j] );

		    diff_pt_vs_m->Fill( g_pt, (*part_jetm)[j] );
		    diff_m_vs_eta->Fill( (*part_jetm)[j], (*part_jeteta)[j] );
		    diff_m_vs_phi->Fill( (*part_jetm)[j], (*part_jetphi)[j] );
		}
	    }
	}
	else{ // if no matching event in second tree, fill "miss" tree vectors
	    cerr << "discrepancy in grant's event " << g_EventID << "\n";

//	    fill vectors with event id, jet pt (eta, phi), number of constituents and constituent pt
	    f_EventID = g_EventID;
	    for(int j = 0; j < (int) part_jetpt->size(); j++){
		fake_jetpt.push_back( (*part_jetpt)[j] );
		fake_jeteta.push_back( (*part_jeteta)[j] );
		fake_jetphi.push_back( (*part_jetphi)[j] );
//		fake_conspt.push_back( (*part_conspt)[j] ); // do i even need the constituent information? maybe not

		fake_ncons.push_back( (*part_conspt)[j].size() );

		cerr << "n constituents? " << (*part_conspt)[j].size() << "\n";
		cerr << "vector of n constituents? " << (*part_conspt).size() << "\n";

		vector<double> pt_cons; vector<double> eta_cons; vector<double> phi_cons;
		for(int k = 0; k < (int) (*part_conspt)[j].size(); k++){
		    pt_cons.push_back( (*part_conspt)[j][k] );
		    eta_cons.push_back( (*part_conseta)[j][k] );
		    phi_cons.push_back( (*part_consphi)[j][k] );
		}

		fake_conspt.push_back( pt_cons );
		fake_conseta.push_back( eta_cons );
		fake_consphi.push_back( phi_cons );

		diff_pt_vs_eta->Fill( (*part_jetpt)[j], (*part_jeteta)[j] );
		diff_pt_vs_phi->Fill( (*part_jetpt)[j], (*part_jetphi)[j] );
		diff_eta_vs_phi->Fill( (*part_jeteta)[j], (*part_jetphi)[j] );

		diff_pt_vs_m->Fill( (*part_jetpt)[j], (*part_jetm)[j] );
		diff_m_vs_eta->Fill( (*part_jetm)[j], (*part_jeteta)[j] );
		diff_m_vs_phi->Fill( (*part_jetm)[j], (*part_jetphi)[j] );
	    }
	}
	if(fake_jetpt.size() != 0){fake->Fill();}
    }


    grant->ResetBranchAddresses();
    isaac->ResetBranchAddresses();
    

    // tree with events, jets in grant's tree but not isaac's
    miss->Write();
    fake->Write();

    diff_pt_vs_eta->Write();
    diff_pt_vs_phi->Write();
    diff_eta_vs_phi->Write();

    diff_pt_vs_m->Write();
    diff_m_vs_eta->Write();
    diff_m_vs_phi->Write();


    cout << "Wrote to " << fout->GetName() << endl;
    
    //closing file
    fout->Close();
    cout << "Closed " << fout->GetName() << endl;
    
    return 0;
}



