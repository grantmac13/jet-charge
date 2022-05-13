//
//  response.cxx
//  take trees from running sim.cxx on pythia 6 which fill trees with matched, missed, fake jet information and construct response
//
//
//  Created by Grant on 2/19/22.
//

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
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hists~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    const int nJetPtBins = 5; // use 4 jet pt ranges, need bin below 10-15 GeV
    double jetPtLo[nJetPtBins] = { 5.0, 10.0, 15.0, 20.0, 30.0};
    double jetPtHi[nJetPtBins] = {10.0, 15.0, 20.0, 30.0, 40.0};
    double jetptEdges[nJetPtBins+1] = {5.0, 10.0, 15.0, 20.0, 30.0, 40.0};
    
    cout << "DEBUG: this gets passed as kappa--" << argv[5] << "\n";
    // 'sim' is currently passed in the position (argv[4]) for kappa.. need to add if(response) in macro_submit_new.csh
    // set arg = "$outLocation $outName $inFiles $inputType $kappa"
    
    string sim_str = (string) argv[4];
    
    if(sim_str != "sim"){
        cerr << "Response.cxx only to be called for sim. Called for: " << argv[4] << ". Exiting.\n";
        exit(1);
    }
    
    
    const string kappa = (string) argv[5];
    double k = str_to_double(kappa);
     
    
//    int nbins = 25; double binlo = -2.5; double binhi = 2.5; // should work, have extra bin on low and hi end for unfolding
    int nbins = 35; double binlo = -3.5; double binhi = 3.5; // should work, have extra bin on low and hi end for unfolding
//    if(k == 0.0){nbins = 9; binlo = -4.5; binhi = 4.5;}
    if(k == 0.0){nbins = 13; binlo = -6.5; binhi = 6.5;}
    // widen histogram to provide bins in response that are not shown in plot (similar to jet pt range change to include 5-10 GeV jets even though i only show 10-15 GeV and above
    double qBinEdges[nbins];
    for(int i = 0; i < nbins; i++){
        qBinEdges[i] = binlo + (binhi - binlo)*i/nbins;
    }
    
    cout << "what should be the bin edge: " << binlo << "\n";
    cout << "bin edge calculated: " << qBinEdges[1] << "\n";
//    TH2D* partQvPt = new TH2D("partQvPt", ";Q^{jet};p_{T} (GeV/c)", nbins, binlo, binhi, 15,  5.0, 80.0); // previously used jet pt binning
    TH2D* partQvPt = new TH2D("partQvPt", ";Q^{jet};p_{T} (GeV/c)", nbins, binlo, binhi, nJetPtBins, jetptEdges); // want particle level to be binned the same as detector level to compare Q in same jet pt range
    TH2D* detQvPt = new TH2D("detQvPt", ";Q^{jet};p_{T} (GeV/c)", nbins, binlo, binhi, nJetPtBins, jetptEdges);
    
    TH1D* det_pt = new TH1D("det_pt", "", nJetPtBins, jetptEdges);
    TH1D* part_pt = new TH1D("part_pt", "", nJetPtBins, jetptEdges);

    TH1D* partQ = new TH1D("partQ", "", nbins, binlo, binhi);
    TH1D* detQ = new TH1D("detQ", "", nbins, binlo, binhi);
    
    cout << "what the bin edge really is: " << partQ->GetBinLowEdge(1) << "\n";
    
    //responses for systematic uncertainty variation
    RooUnfoldResponse *q_pt_res_nom = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_nom"); //nominal
    RooUnfoldResponse *q_pt_res_TS = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_TS"); //tower scale
    RooUnfoldResponse *q_pt_res_TU = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_TU"); //tracking uncertainty
    RooUnfoldResponse *q_pt_res_HC50 = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_HC50"); //hadronic correction
    RooUnfoldResponse *q_pt_res_DS = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_DS"); //smear detector spectrum
    RooUnfoldResponse *q_pt_res_GS = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_GS"); //shift generator spectrum
    
    
    // jetptEdges[nJetPtBins+1]
    TH1D* match_plus_miss = new TH1D("jetpt_match_plus_miss", "", nJetPtBins, jetptEdges);
    
    RooUnfoldResponse *q_pt_response = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_response");
    RooUnfoldResponse *q_pt_response_counts = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_response_counts");
    RooUnfoldResponse *pt_response = new RooUnfoldResponse(det_pt, part_pt, "pt_response");
    RooUnfoldResponse *q_response = new RooUnfoldResponse(detQ, partQ, "q_response");

    std::vector<RooUnfoldResponse*> syst_res = {q_pt_res_nom, q_pt_res_TS, q_pt_res_TU, q_pt_res_HC50, q_pt_res_DS, q_pt_res_GS};
    

    RooUnfoldResponse *res_1D_1015 = new RooUnfoldResponse(detQ, partQ, "q_response_1015");
    RooUnfoldResponse *res_1D_1520 = new RooUnfoldResponse(detQ, partQ, "q_response_1520");
    RooUnfoldResponse *res_1D_2030 = new RooUnfoldResponse(detQ, partQ, "q_response_2030");
    RooUnfoldResponse *res_1D_3040 = new RooUnfoldResponse(detQ, partQ, "q_response_3040");



    ////////////////////////////////////////////
    // histograms for closure test
    TH1D* sampleA_pt_gen = new TH1D("sampleA_pt_gen", "", nJetPtBins, jetptEdges);
    TH1D* sampleA_pt_det = new TH1D("sampleA_pt_det", "", nJetPtBins, jetptEdges);
    TH1D* sampleA_q_gen = new TH1D("sampleA_q_gen", "", nbins, binlo, binhi);
    TH1D* sampleA_q_det = new TH1D("sampleA_q_det", "", nbins, binlo, binhi);

    TH1D* sampleB_pt_gen = new TH1D("sampleB_pt_gen", "", nJetPtBins, jetptEdges);
    TH1D* sampleB_pt_det = new TH1D("sampleB_pt_det", "", nJetPtBins, jetptEdges);
    TH1D* sampleB_q_gen = new TH1D("sampleB_q_gen", "", nbins, binlo, binhi);
    TH1D* sampleB_q_det = new TH1D("sampleB_q_det", "", nbins, binlo, binhi);


    TH2D* sampleA_q_pt_gen = new TH2D("sampleA_q_pt_gen", "", nbins, binlo, binhi, nJetPtBins, jetptEdges);
    TH2D* sampleA_q_pt_det = new TH2D("sampleA_q_pt_det", "", nbins, binlo, binhi, nJetPtBins, jetptEdges);

    TH2D* sampleB_q_pt_gen = new TH2D("sampleB_q_pt_gen", "", nbins, binlo, binhi, nJetPtBins, jetptEdges);
    TH2D* sampleB_q_pt_det = new TH2D("sampleB_q_pt_det", "", nbins, binlo, binhi, nJetPtBins, jetptEdges);


    // counts histograms for bin dropping
    TH1D* sampleA_pt_gen_counts = new TH1D("sampleA_pt_gen_counts", "", nJetPtBins, jetptEdges);
    TH1D* sampleA_pt_det_counts = new TH1D("sampleA_pt_det_counts", "", nJetPtBins, jetptEdges);
    TH1D* sampleA_q_gen_counts = new TH1D("sampleA_q_gen_counts", "", nbins, binlo, binhi);
    TH1D* sampleA_q_det_counts = new TH1D("sampleA_q_det_counts", "", nbins, binlo, binhi);

    TH1D* sampleB_pt_gen_counts = new TH1D("sampleB_pt_gen_counts", "", nJetPtBins, jetptEdges);
    TH1D* sampleB_pt_det_counts = new TH1D("sampleB_pt_det_counts", "", nJetPtBins, jetptEdges);
    TH1D* sampleB_q_gen_counts = new TH1D("sampleB_q_gen_counts", "", nbins, binlo, binhi);
    TH1D* sampleB_q_det_counts = new TH1D("sampleB_q_det_counts", "", nbins, binlo, binhi);


    TH2D* sampleA_q_pt_gen_counts = new TH2D("sampleA_q_pt_gen_counts", "", nbins, binlo, binhi, nJetPtBins, jetptEdges);
    TH2D* sampleA_q_pt_det_counts = new TH2D("sampleA_q_pt_det_counts", "", nbins, binlo, binhi, nJetPtBins, jetptEdges);

    TH2D* sampleB_q_pt_gen_counts = new TH2D("sampleB_q_pt_gen_counts", "", nbins, binlo, binhi, nJetPtBins, jetptEdges);
    TH2D* sampleB_q_pt_det_counts = new TH2D("sampleB_q_pt_det_counts", "", nbins, binlo, binhi, nJetPtBins, jetptEdges);
    ////////////////////////////////////////////


    ////////////////////////////////////////////
    // responses for closure test
    RooUnfoldResponse *sampleA_pt_response = new RooUnfoldResponse(det_pt, part_pt, "sampleA_pt_response");
    RooUnfoldResponse *sampleA_q_response = new RooUnfoldResponse(nbins, binlo, binhi, "sampleA_q_response");

    RooUnfoldResponse *sampleA_q_pt_response = new RooUnfoldResponse(detQvPt, partQvPt, "sampleA_q_pt_response");
    RooUnfoldResponse *sampleA_q_pt_response_counts = new RooUnfoldResponse(detQvPt, partQvPt, "sampleA_q_pt_response_counts");


    RooUnfoldResponse *sampleB_pt_response = new RooUnfoldResponse(det_pt, part_pt, "sampleB_pt_response");
    RooUnfoldResponse *sampleB_q_response = new RooUnfoldResponse(nbins, binlo, binhi, "sampleB_q_response");

    RooUnfoldResponse *sampleB_q_pt_response = new RooUnfoldResponse(detQvPt, partQvPt, "sampleB_q_pt_response");
    RooUnfoldResponse *sampleB_q_pt_response_counts = new RooUnfoldResponse(detQvPt, partQvPt, "sampleB_q_pt_response_counts");
    ////////////////////////////////////////////

    vector<RooUnfoldResponse*> base_res = {q_response, pt_response, q_pt_response, q_pt_response_counts};

    vector<RooUnfoldResponse*> sampleA_res = {sampleA_q_pt_response, sampleA_q_pt_response_counts, sampleA_q_response, sampleA_pt_response};
    vector<TH1D*> sampleA_h1D = {sampleA_pt_gen, sampleA_pt_det, sampleA_pt_gen_counts, sampleA_pt_det_counts, sampleA_q_gen, sampleA_q_det, sampleA_q_gen_counts, sampleA_q_det_counts};
    vector<TH2D*> sampleA_h2D = {sampleA_q_pt_gen, sampleA_q_pt_det, sampleA_q_pt_gen_counts, sampleA_q_pt_det_counts};

    vector<RooUnfoldResponse*> sampleB_res = {sampleB_q_pt_response, sampleB_q_pt_response_counts, sampleB_q_response, sampleB_pt_response};
    vector<TH1D*> sampleB_h1D = {sampleB_pt_gen, sampleB_pt_det, sampleB_pt_gen_counts, sampleB_pt_det_counts, sampleB_q_gen, sampleB_q_det, sampleB_q_gen_counts, sampleB_q_det_counts};
    vector<TH2D*> sampleB_h2D = {sampleB_q_pt_gen, sampleB_q_pt_det, sampleB_q_pt_gen_counts, sampleB_q_pt_det_counts};


    int nSources = 6;
    
    string tree_name[nSources] = {"jetChargeTree", "towerscale", "trackingeff", "hadroncorr50", "detsmear", "gensmear"};
    
    for(int iSyst = 0; iSyst < nSources; iSyst++){
        
        vector<double> *det_jetpt = 0;
        vector<double> *part_jetpt = 0;
        vector<vector<double> > *det_conspt = 0; vector<vector<double> > *det_chcons = 0;
        vector<vector<double> > *part_conspt = 0; vector<vector<double> > *part_chcons = 0;
        
        vector<double> *miss_jetpt = 0;
        vector<vector<double> > *miss_conspt = 0; vector<vector<double> > *miss_chcons = 0;
        vector<double> *fake_jetpt = 0;
        vector<vector<double> > *fake_conspt = 0; vector<vector<double> > *fake_chcons = 0;
    
        double weight = 1;    

        TTree *t = (TTree*) fin->Get( (tree_name[iSyst]).c_str() ); // should be the same name of the tree for all data, pythia6, pythia8
        
        t->SetBranchAddress("det_jetpt", &det_jetpt);
        t->SetBranchAddress("part_jetpt", &part_jetpt);
        t->SetBranchAddress("det_conspt", &det_conspt);
        t->SetBranchAddress("part_conspt", &part_conspt);
        t->SetBranchAddress("det_chcons", &det_chcons);
        t->SetBranchAddress("part_chcons", &part_chcons);
        
        t->SetBranchAddress("miss_jetpt", &miss_jetpt);
        t->SetBranchAddress("miss_conspt", &miss_conspt);
        t->SetBranchAddress("miss_chcons", &miss_chcons);
        
        t->SetBranchAddress("fake_jetpt", &fake_jetpt);
        t->SetBranchAddress("fake_conspt", &fake_conspt);
        t->SetBranchAddress("fake_chcons", &fake_chcons);
        
        t->SetBranchAddress("weight", &weight);
    
        for(int i = 0; i < (int) t->GetEntries(); i++){
            // take each event
            t->GetEntry(i);
            double clos_rand = gRandom->Uniform(0.0, 1.0); //randomly split events into sampleA, sampleB



            for(int j = 0; j < (int) det_jetpt->size(); j++){
                // loop over matched jets and fill response
                double jc_p = 0.0; double jc_d = 0.0;
                double p_pt = part_jetpt->at(j);
                double d_pt =  det_jetpt->at(j);
                
		if(iSyst == 0){
	                match_plus_miss->Fill(p_pt, weight);
                }
                for(int ii = 0; ii < (int) max( part_conspt->at(j).size(), det_conspt->at(j).size() ); ii++){
                    if(ii < (int) part_conspt->at(j).size()){
                        double p_cons_pt = part_conspt->at(j).at(ii);
                        jc_p += pow( (p_cons_pt/p_pt) , k ) * part_chcons->at(j).at(ii);
                    }
                    if(ii < (int) det_conspt->at(j).size()){
                        double d_cons_pt = det_conspt->at(j).at(ii);
                        jc_d += pow( (d_cons_pt/d_pt) , k ) * det_chcons->at(j).at(ii);
                    }
                }
                syst_res[iSyst]->Fill(jc_d, d_pt, jc_p, p_pt, weight);
                if(iSyst == 0){
                    q_pt_response->Fill(jc_d, d_pt, jc_p, p_pt, weight);
                    q_pt_response_counts->Fill(jc_d, d_pt, jc_p, p_pt);
                    pt_response->Fill(d_pt, p_pt, weight); // fill jet pt response with matched jets

                    if(clos_rand > 0.5){
                        sampleA_pt_gen->Fill(p_pt, weight);
                        sampleA_pt_det->Fill(d_pt, weight);

                        sampleA_pt_gen_counts->Fill(p_pt);
                        sampleA_pt_det_counts->Fill(d_pt);


                        sampleA_q_gen->Fill(jc_p, weight);
                        sampleA_q_det->Fill(jc_d, weight);

                        sampleA_q_gen_counts->Fill(jc_p);
                        sampleA_q_det_counts->Fill(jc_d);


                        sampleA_q_pt_gen->Fill(jc_p, p_pt, weight);
                        sampleA_q_pt_det->Fill(jc_d, d_pt, weight);

                        sampleA_q_pt_gen_counts->Fill(jc_p, p_pt);
                        sampleA_q_pt_det_counts->Fill(jc_d, d_pt);


                        sampleA_pt_response->Fill(d_pt, p_pt, weight);
                        sampleA_q_response->Fill(jc_d, jc_p, weight);

                        sampleA_q_pt_response->Fill(jc_d, d_pt, jc_p, p_pt, weight);
                    }
                    if(clos_rand < 0.5){
                        sampleB_pt_gen->Fill(p_pt, weight);
                        sampleB_pt_det->Fill(d_pt, weight);

                        sampleB_pt_gen_counts->Fill(p_pt);
                        sampleB_pt_det_counts->Fill(d_pt);


                        sampleB_q_gen->Fill(jc_p, weight);
                        sampleB_q_det->Fill(jc_d, weight);

                        sampleB_q_gen_counts->Fill(jc_p);
                        sampleB_q_det_counts->Fill(jc_d);


                        sampleB_q_pt_gen->Fill(jc_p, p_pt, weight);
                        sampleB_q_pt_det->Fill(jc_d, d_pt, weight);

                        sampleB_q_pt_gen_counts->Fill(jc_p, p_pt);
                        sampleB_q_pt_det_counts->Fill(jc_d, d_pt);


                        sampleB_pt_response->Fill(d_pt, p_pt, weight);
                        sampleB_q_response->Fill(jc_d, jc_p, weight);

                        sampleB_q_pt_response->Fill(jc_d, d_pt, jc_p, p_pt, weight);
                    }
                


                    // for 1D unfolding: matched jets that have both particle and detector level in pT range, use ->Fill()
                    // "detector level jets that match with particle level jet with jet pt in pt range but detector jet pt falls outside pT range go into fakes"
                    // "particle level jets that match with detector level jet with jet pt in pt range but particle jet pt falls outside pT range go into misses"
                    if(10.0 < d_pt && d_pt < 15.0 && 10.0 < p_pt && p_pt < 15.0){ // am I requiring both detector and particle level jet pt or only one?
                        res_1D_1015->Fill(jc_d, jc_p, weight);
                    }
                    else if(10.0 < p_pt && p_pt < 15.0){
                        res_1D_1015->Miss(jc_p, weight);
                    }
                    else if(10.0 < d_pt && d_pt < 15.0){
                        res_1D_1015->Fake(jc_d, weight);
                    }

                    if(15.0 < d_pt && d_pt < 20.0 && 15.0 < p_pt && p_pt < 20.0){ // am I requiring both detector and particle level jet pt or only one?
                        res_1D_1520->Fill(jc_d, jc_p, weight);
                    }
                    else if(15.0 < p_pt && p_pt < 20.0){
                        res_1D_1520->Miss(jc_p, weight);
                    }
                    else if(15.0 < d_pt && d_pt < 20.0){
                        res_1D_1520->Fake(jc_d, weight);
                    }

                    if(20.0 < d_pt && d_pt < 30.0 && 20.0 < p_pt && p_pt < 30.0){ // am I requiring both detector and particle level jet pt or only one?
                        res_1D_2030->Fill(jc_d, jc_p, weight);
                    }
                    else if(20.0 < p_pt && p_pt < 30.0){
                        res_1D_2030->Miss(jc_p, weight);
                    }
                    else if(20.0 < d_pt && d_pt < 30.0){
                        res_1D_2030->Fake(jc_d, weight);
                    }

                    if(30.0 < d_pt && d_pt < 40.0 && 30.0 < p_pt && p_pt < 40.0){ // am I requiring both detector and particle level jet pt or only one?
                        res_1D_3040->Fill(jc_d, jc_p, weight);
                    }
                    else if(30.0 < p_pt && p_pt < 40.0){
                        res_1D_3040->Miss(jc_p, weight);
                    }
                    else if(30.0 < d_pt && d_pt < 40.0){
                        res_1D_3040->Fake(jc_d, weight);
                    }
		}
            }
        
            for(int jj = 0; jj < (int) miss_jetpt->size(); jj++){
                // fill response misses with misses from particle level
                double jc_m = 0.0;
                double m_pt = miss_jetpt->at(jj);

		if(iSyst == 0){
                    match_plus_miss->Fill(m_pt, weight);
                }
                for(int ii = 0; ii < (int) miss_conspt->at(jj).size(); ii++){
                    double m_cons_pt = miss_conspt->at(jj).at(ii);
                    jc_m += pow( (m_cons_pt/m_pt), k ) * miss_chcons->at(jj).at(ii);
                }
                syst_res[iSyst]->Miss(jc_m, m_pt, weight);
                if(iSyst == 0){
                    q_pt_response->Miss(jc_m, m_pt, weight);
                    q_pt_response_counts->Miss(jc_m, m_pt);
                    pt_response->Miss(m_pt, weight);
                
                    if(clos_rand > 0.5){
                        sampleA_pt_gen->Fill(m_pt, weight);
                        sampleA_q_gen->Fill(jc_m, weight);

                        sampleA_pt_response->Miss(m_pt, weight);
                        sampleA_q_response->Miss(jc_m, weight);

                        sampleA_q_pt_gen->Fill(jc_m, m_pt, weight);
                        sampleA_q_pt_gen_counts->Fill(jc_m, m_pt);

                        sampleA_q_pt_response->Miss(jc_m, m_pt, weight);
                    }
                    if(clos_rand < 0.5){
                        sampleB_pt_gen->Fill(m_pt, weight);
                        sampleB_q_gen->Fill(jc_m, weight);

                        sampleB_pt_response->Miss(m_pt, weight);
                        sampleB_q_response->Miss(jc_m, weight);

                        sampleB_q_pt_gen->Fill(jc_m, m_pt, weight);
                        sampleB_q_pt_gen_counts->Fill(jc_m, m_pt);

                        sampleB_q_pt_response->Miss(jc_m, m_pt, weight);
                    }
		}
            }
            for(int jj = 0; jj < (int) fake_jetpt->size(); jj++){
                double jc_f = 0.0;
                double f_pt =  fake_jetpt->at(jj);
            
                // fill response fakes with fakes from detector level
                for(int ii = 0; ii < (int) fake_conspt->at(jj).size(); ii++){
                    double f_cons_pt = fake_conspt->at(jj).at(ii);
                    jc_f += pow( (f_cons_pt/f_pt), k ) * fake_chcons->at(jj).at(ii);
                }
                syst_res[iSyst]->Fake(jc_f, f_pt, weight);
                if(iSyst == 0){
                    q_pt_response->Fake(jc_f, f_pt, weight);
                    q_pt_response_counts->Fake(jc_f, f_pt);
                    pt_response->Fake(f_pt, weight);
                
                    if(clos_rand > 0.5){
                        sampleA_pt_response->Fake(f_pt, weight);
                        sampleA_q_response->Fake(jc_f, weight);

                        sampleA_q_pt_det->Fill(jc_f, f_pt, weight);
                        sampleA_q_pt_det_counts->Fill(jc_f, f_pt);

                        sampleA_q_pt_response->Fake(jc_f, f_pt, weight);
                    }
                    if(clos_rand < 0.5){
                        sampleB_pt_response->Fake(f_pt, weight);
                        sampleB_q_response->Fake(jc_f, weight);

                        sampleB_q_pt_det->Fill(jc_f, f_pt, weight);
                        sampleB_q_pt_det_counts->Fill(jc_f, f_pt);

                        sampleB_q_pt_response->Fake(jc_f, f_pt, weight);
		    }
                }
            }
        }

        t->ResetBranchAddresses();
    
    } // end systematic loop
    
    cout << "is the bin edge still where it should be?: " << partQ->GetBinLowEdge(1) << "\n";
    
    TFile *fout = new TFile(((string) argv[1]+(string) argv[2]).c_str(),"RECREATE");
    cout << "DEBUG: output file name is " << fout->GetName() << endl;
    
    //q_pt_response->Write();
    //q_pt_response_counts->Write();
    //q_response->Write();
    //pt_response->Write();
    match_plus_miss->Write();

    res_1D_1015->Write();
    res_1D_1520->Write();
    res_1D_2030->Write();
    res_1D_3040->Write();

    for(int i = 0; i < (int) base_res.size(); i++){
        base_res[i]->Write();
    }
    for(int i = 0; i < (int) sampleA_res.size(); i++){
        sampleA_res[i]->Write();
    }
    for(int i = 0; i < (int) sampleA_h1D.size(); i++){
        sampleA_h1D[i]->Write();
    }
    for(int i = 0; i < (int) sampleA_h2D.size(); i++){
        sampleA_h2D[i]->Write();
    }
    for(int i = 0; i < (int) sampleB_res.size(); i++){
        sampleB_res[i]->Write();
    }
    for(int i = 0; i < (int) sampleB_h1D.size(); i++){
        sampleB_h1D[i]->Write();
    }
    for(int i = 0; i < (int) sampleB_h2D.size(); i++){
        sampleB_h2D[i]->Write();
    }
    for(int i = 0; i < (int) syst_res.size(); i++){
        syst_res[i]->Write();
    }
    

    cout << "Wrote to " << fout->GetName() << endl;
    
    //closing file
    fout->Close();
    cout << "Closed " << fout->GetName() << endl;
    
    return 0;
}



