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
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    
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
    const int nJetPtBins = 4;
//    double jetPtLo[nJetPtBins] = {15.0, 20.0, 25.0, 30.0};
//    double jetPtHi[nJetPtBins] = {20.0, 25.0, 30.0, 40.0};

    // may want to uniform-ize bins to 5 GeV consistently
    double jetptEdges[nJetPtBins+1] = {15.0, 20.0, 25.0, 30.0, 40.0};
//    int jetptHiEdge = 50;

    // temporary to make 5 GeV bins
//    const int nJetPtBins = 9; // use 4 jet pt ranges, need bin below 10-15 GeV
//    double jetptEdges[nJetPtBins+1] = {15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0};

//    const int npartJetPtBins = 15; // use 4 jet pt ranges, need bin below 10-15 GeV
//    double partjetptEdges[npartJetPtBins+1] =
//        {5.0, 10.0, 15.0, 20.0, 25.0,
//        30.0, 35.0, 40.0, 45.0, 50.0,
//        55.0, 60.0, 65.0, 70.0, 75.0, 80.0};
    
//    cout << "DEBUG: this gets passed as kappa--" << argv[5] << "\n";
    // 'sim' is currently passed in the position (argv[4]) for kappa.. need to add if(response) in macro_submit_new.csh
    // set arg = "$outLocation $outName $inFiles $inputType $kappa"
    
    string sim_str = (string) argv[4];
    
    if(sim_str != "sim"){
        cerr << "Response.cxx only to be called for sim. Called for: " << argv[4] << ". Exiting.\n";
        exit(1);
    }
    
    
    const string kappa = (string) argv[5];
    double k = str_to_double(kappa);
    
    
    // the + 1/2 * binwidth is to center a bin at 0

    double qBinWidth = 0.1; //default to ~ CMS binning
    // default to non-k=0.0 bins
    double jetQedge = 4.5 + (1.0/2.0)*qBinWidth; // pre-11/14/22: using 18 bins in [-3.5, 3.5] gives bin width of ~ 7/18 ~= 0.39 [e]
    // what if i use bin width comparable to CMS' choice? ~ 21 bins in [-1, 1] which is roughly 2/21 ~= 0.1 [e]
    if(k == 0.0){qBinWidth = 1.0; jetQedge = 6.0 + (1.0/2.0)*qBinWidth;} // no real choice in binning for k = 0.0 since k = 0.0 corresponds to adding integers together
    else if(k < 0.0){jetQedge = 10.0 + (1.0/2.0)*qBinWidth;} // bin width = 20/100 = 0.2, can even get finer bins; not even as fine as CMS' yet (jet to study effect of bin width in fitting)

    int njetQbins = (2 * jetQedge)/qBinWidth; //


//    cout << "what should be the bin edge: " << binlo << "\n";
//    cout << "bin edge calculated: " << qBinEdges[1] << "\n";
//    TH2D* partQvPt = new TH2D("partQvPt", ";Q^{jet};p_{T} (GeV/c)", nbins, binlo, binhi, 15,  5.0, 80.0); // previously used jet pt binning

    TH2D* partQvPt = new TH2D("partQvPt", ";Q^{jet};p_{T} (GeV/c)", njetQbins, -jetQedge, jetQedge, 15, 5, 80); // want particle level to be binned the same as detector level to compare Q in same jet pt range
    TH2D* detQvPt = new TH2D("detQvPt", ";Q^{jet};p_{T} (GeV/c)", njetQbins, -jetQedge, jetQedge, 9, 15, 60);
    
    TH1D* det_pt = new TH1D("det_pt", "", 9, 15, 60);
    TH1D* part_pt = new TH1D("part_pt", "", 15, 5, 80);


    TH1D* partQ = new TH1D("partQ", "", njetQbins, -jetQedge, jetQedge);
    TH1D* detQ = new TH1D("detQ", "", njetQbins, -jetQedge, jetQedge);


    TH2D* partMvPt = new TH2D("partMvPt", ";M (GeV/c);p_{T} (GeV/c)", 14, 0.0, 14.0, 15, 5, 80); // want particle level to be binned the same as detector level to compare Q in same jet pt range
    TH2D* detMvPt = new TH2D("detMvPt", ";M (GeV/c);p_{T} (GeV/c)", 14, 0.0, 14.0, 9, 15, 60);


    TH1D* partM = new TH1D("partM", "", 14, 0.0, 14.0);
    TH1D* detM = new TH1D("detM", "", 14, 0.0, 14.0);

    
//    cout << "what the bin edge really is: " << partQ->GetBinLowEdge(1) << "\n";


    // TEMPORARY: HARDCODING RADIUS HERE

    // out/sim/hists/unmatchedp6_R04_k00.root
//    TFile *p6unmatch = new TFile("~/ppJCandBF/out/sim/hists/unmatchedp6_R"+rad+"_k"+kapstr+".root", "READ");
    TFile *p6unmatch = new TFile( ("~/ppJCandBF/out/sim/hists/unmatchedp6_R04_k"+kappa+".root").c_str(), "READ");
    TH2D* qp6_2D = (TH2D*) p6unmatch->Get("QvPt_p");
    qp6_2D->SetDirectory(0);
    p6unmatch->Close();
//    TFile *h7unmatch = new TFile("~/for_grant/out/unmatchedh7_R"+rad+"_k"+kapstr+".root", "READ");
    TFile *h7unmatch = new TFile( ("~/for_grant/out/unmatchedh7_R04_k"+kappa+".root").c_str(), "READ");
    TH2D* qh7_2D = (TH2D*) h7unmatch->Get("qvpt");
    qh7_2D->SetDirectory(0);
    h7unmatch->Close();
//    TFile *p8unmatch = new TFile("~/for_grant/out/unmatchedp8_R"+rad+"_k"+kapstr+".root", "READ");
    TFile *p8unmatch = new TFile( ("~/for_grant/out/unmatchedp8_R04_k"+kappa+".root").c_str(), "READ");
    TH2D* qp8_2D = (TH2D*) p8unmatch->Get("qvpt");
    qp8_2D->SetDirectory(0);
    p8unmatch->Close();


    vector<TH1D*> qp6 = {qp6_2D->ProjectionX("qp6_0", qp6_2D->GetYaxis()->FindBin(20), qp6_2D->GetYaxis()->FindBin(25)-1),
			 qp6_2D->ProjectionX("qp6_1", qp6_2D->GetYaxis()->FindBin(25), qp6_2D->GetYaxis()->FindBin(30)-1),
			 qp6_2D->ProjectionX("qp6_2", qp6_2D->GetYaxis()->FindBin(30), qp6_2D->GetYaxis()->FindBin(40)-1),
			};
    vector<TH1D*> qh7 = {qh7_2D->ProjectionX("qh7_0", qh7_2D->GetYaxis()->FindBin(20), qh7_2D->GetYaxis()->FindBin(25)-1),
			 qh7_2D->ProjectionX("qh7_1", qh7_2D->GetYaxis()->FindBin(25), qh7_2D->GetYaxis()->FindBin(30)-1),
			 qh7_2D->ProjectionX("qh7_2", qh7_2D->GetYaxis()->FindBin(30), qh7_2D->GetYaxis()->FindBin(40)-1),
			};
    vector<TH1D*> qp8 = {qp8_2D->ProjectionX("qp8_0", qp8_2D->GetYaxis()->FindBin(20), qp8_2D->GetYaxis()->FindBin(25)-1),
			 qp8_2D->ProjectionX("qp8_1", qp8_2D->GetYaxis()->FindBin(25), qp8_2D->GetYaxis()->FindBin(30)-1),
			 qp8_2D->ProjectionX("qp8_2", qp8_2D->GetYaxis()->FindBin(30), qp8_2D->GetYaxis()->FindBin(40)-1),
			};


    for(int i = 0; i < (int) qp6.size(); i++){
        qh7[i]->Divide(qp6[i]);
        qp8[i]->Divide(qp6[i]);
    }





    //responses for systematic uncertainty variation
    RooUnfoldResponse *q_pt_res_nom = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_nom"); //nominal
    RooUnfoldResponse *q_pt_res_TS = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_TS"); //tower scale
    RooUnfoldResponse *q_pt_res_TU = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_TU"); //tracking uncertainty
    RooUnfoldResponse *q_pt_res_HC50 = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_HC50"); //hadronic correction
    RooUnfoldResponse *q_pt_res_DS = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_DS"); //smear detector spectrum
    RooUnfoldResponse *q_pt_res_GS = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_res_GS"); //shift generator spectrum
    
    
    // jetptEdges[nJetPtBins+1]
    TH1D* match_plus_miss = new TH1D("jetpt_match_plus_miss", "", 15, 5, 80);
    
    RooUnfoldResponse *q_pt_response = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_response");
    RooUnfoldResponse *q_pt_response_counts = new RooUnfoldResponse(detQvPt, partQvPt, "q_pt_response_counts");
    RooUnfoldResponse *pt_response = new RooUnfoldResponse(det_pt, part_pt, "pt_response");
    RooUnfoldResponse *q_response = new RooUnfoldResponse(detQ, partQ, "q_response");


    RooUnfoldResponse *m_pt_response = new RooUnfoldResponse(detMvPt, partMvPt, "m_pt_response");
    RooUnfoldResponse *m_response = new RooUnfoldResponse(detM, partM, "m_response");


    std::vector<RooUnfoldResponse*> syst_res = {q_pt_res_nom, q_pt_res_TS, q_pt_res_TU, q_pt_res_HC50, q_pt_res_DS, q_pt_res_GS};
    

    RooUnfoldResponse *q_res2025_nom = new RooUnfoldResponse(detQ, partQ, "q_res2025_nom");
    RooUnfoldResponse *q_res2530_nom = new RooUnfoldResponse(detQ, partQ, "q_res2530_nom");
    RooUnfoldResponse *q_res3040_nom = new RooUnfoldResponse(detQ, partQ, "q_res3040_nom");

    RooUnfoldResponse *q_res2025_h7smear = new RooUnfoldResponse(detQ, partQ, "q_res2025_h7smear");
    RooUnfoldResponse *q_res2530_h7smear = new RooUnfoldResponse(detQ, partQ, "q_res2530_h7smear");
    RooUnfoldResponse *q_res3040_h7smear = new RooUnfoldResponse(detQ, partQ, "q_res3040_h7smear");

    RooUnfoldResponse *q_res2025_p8smear = new RooUnfoldResponse(detQ, partQ, "q_res2025_p8smear");
    RooUnfoldResponse *q_res2530_p8smear = new RooUnfoldResponse(detQ, partQ, "q_res2530_p8smear");
    RooUnfoldResponse *q_res3040_p8smear = new RooUnfoldResponse(detQ, partQ, "q_res3040_p8smear");



    ////////////////////////////////////////////
    // histograms for closure test
    TH1D* sampleA_pt_gen = new TH1D("sampleA_pt_gen", "", 15, 5, 80);
    TH1D* sampleA_pt_det = new TH1D("sampleA_pt_det", "", 9, 15, 60);
    TH1D* sampleA_q_gen = new TH1D("sampleA_q_gen", "", njetQbins, -jetQedge, jetQedge);
    TH1D* sampleA_q_det = new TH1D("sampleA_q_det", "", njetQbins, -jetQedge, jetQedge);

    TH1D* sampleA_m_gen = new TH1D("sampleA_m_gen", "", 14, 0.0, 14.0);
    TH1D* sampleA_m_det = new TH1D("sampleA_m_det", "", 14, 0.0, 14.0);


    TH1D* sampleB_pt_gen = new TH1D("sampleB_pt_gen", "", 15, 5, 80);
    TH1D* sampleB_pt_det = new TH1D("sampleB_pt_det", "", 9, 15, 60);
    TH1D* sampleB_q_gen = new TH1D("sampleB_q_gen", "", njetQbins, -jetQedge, jetQedge);
    TH1D* sampleB_q_det = new TH1D("sampleB_q_det", "", njetQbins, -jetQedge, jetQedge);

    TH1D* sampleB_m_gen = new TH1D("sampleB_m_gen", "", 14, 0.0, 14.0);
    TH1D* sampleB_m_det = new TH1D("sampleB_m_det", "", 14, 0.0, 14.0);


    TH2D* sampleA_q_pt_gen = new TH2D("sampleA_q_pt_gen", "", njetQbins, -jetQedge, jetQedge, 15, 5, 80);
    TH2D* sampleA_q_pt_det = new TH2D("sampleA_q_pt_det", "", njetQbins, -jetQedge, jetQedge, 9, 15, 60);

    TH2D* sampleB_q_pt_gen = new TH2D("sampleB_q_pt_gen", "", njetQbins, -jetQedge, jetQedge, 15, 5, 80);
    TH2D* sampleB_q_pt_det = new TH2D("sampleB_q_pt_det", "", njetQbins, -jetQedge, jetQedge, 9, 15, 60);


    TH2D* sampleA_m_pt_gen = new TH2D("sampleA_m_pt_gen", "", 14, 0.0, 14.0, 15, 5, 80);
    TH2D* sampleA_m_pt_det = new TH2D("sampleA_m_pt_det", "", 14, 0.0, 14.0, 9, 15, 60);

    TH2D* sampleB_m_pt_gen = new TH2D("sampleB_m_pt_gen", "", 14, 0.0, 14.0, 15, 5, 80);
    TH2D* sampleB_m_pt_det = new TH2D("sampleB_m_pt_det", "", 14, 0.0, 14.0, 9, 15, 60);



    // counts histograms for bin dropping
    TH1D* sampleA_pt_gen_counts = new TH1D("sampleA_pt_gen_counts", "", 15, 5, 80);
    TH1D* sampleA_pt_det_counts = new TH1D("sampleA_pt_det_counts", "", 9, 15, 60);
    TH1D* sampleA_q_gen_counts = new TH1D("sampleA_q_gen_counts", "", njetQbins, -jetQedge, jetQedge);
    TH1D* sampleA_q_det_counts = new TH1D("sampleA_q_det_counts", "", njetQbins, -jetQedge, jetQedge);

    TH1D* sampleB_pt_gen_counts = new TH1D("sampleB_pt_gen_counts", "", 15, 5, 80);
    TH1D* sampleB_pt_det_counts = new TH1D("sampleB_pt_det_counts", "", 9, 15, 60);
    TH1D* sampleB_q_gen_counts = new TH1D("sampleB_q_gen_counts", "", njetQbins, -jetQedge, jetQedge);
    TH1D* sampleB_q_det_counts = new TH1D("sampleB_q_det_counts", "", njetQbins, -jetQedge, jetQedge);


    TH2D* sampleA_q_pt_gen_counts = new TH2D("sampleA_q_pt_gen_counts", "", njetQbins, -jetQedge, jetQedge, 15, 5, 80);
    TH2D* sampleA_q_pt_det_counts = new TH2D("sampleA_q_pt_det_counts", "", njetQbins, -jetQedge, jetQedge, 9, 15, 60);

    TH2D* sampleB_q_pt_gen_counts = new TH2D("sampleB_q_pt_gen_counts", "", njetQbins, -jetQedge, jetQedge, 15, 5, 80);
    TH2D* sampleB_q_pt_det_counts = new TH2D("sampleB_q_pt_det_counts", "", njetQbins, -jetQedge, jetQedge, 9, 15, 60);
    ////////////////////////////////////////////


    ////////////////////////////////////////////
    // responses for closure test
    RooUnfoldResponse *sampleA_pt_response = new RooUnfoldResponse(det_pt, part_pt, "sampleA_pt_response");
    RooUnfoldResponse *sampleA_q_response = new RooUnfoldResponse(njetQbins, -jetQedge, jetQedge, "sampleA_q_response");

    RooUnfoldResponse *sampleA_m_response = new RooUnfoldResponse(detM, partM, "sampleA_m_response");


    RooUnfoldResponse *sampleA_q_pt_response = new RooUnfoldResponse(detQvPt, partQvPt, "sampleA_q_pt_response");
    RooUnfoldResponse *sampleA_q_pt_response_counts = new RooUnfoldResponse(detQvPt, partQvPt, "sampleA_q_pt_response_counts");

    RooUnfoldResponse *sampleA_m_pt_response = new RooUnfoldResponse(detMvPt, partMvPt, "sampleA_m_pt_response");


    RooUnfoldResponse *sampleB_pt_response = new RooUnfoldResponse(det_pt, part_pt, "sampleB_pt_response");
    RooUnfoldResponse *sampleB_q_response = new RooUnfoldResponse(njetQbins, -jetQedge, jetQedge, "sampleB_q_response");

    RooUnfoldResponse *sampleB_m_response = new RooUnfoldResponse(detM, partM, "sampleB_m_response");


    RooUnfoldResponse *sampleB_q_pt_response = new RooUnfoldResponse(detQvPt, partQvPt, "sampleB_q_pt_response");
    RooUnfoldResponse *sampleB_q_pt_response_counts = new RooUnfoldResponse(detQvPt, partQvPt, "sampleB_q_pt_response_counts");

    RooUnfoldResponse *sampleB_m_pt_response = new RooUnfoldResponse(detMvPt, partMvPt, "sampleB_m_pt_response");
    ////////////////////////////////////////////

    vector<RooUnfoldResponse*> base_res = {q_response, pt_response, q_pt_response, q_pt_response_counts};

    vector<RooUnfoldResponse*> sampleA_res = {sampleA_q_pt_response, sampleA_q_pt_response_counts, sampleA_q_response, sampleA_pt_response};
    vector<TH1D*> sampleA_h1D = {sampleA_pt_gen, sampleA_pt_det, sampleA_pt_gen_counts, sampleA_pt_det_counts, sampleA_q_gen, sampleA_q_det, sampleA_q_gen_counts, sampleA_q_det_counts};
    vector<TH2D*> sampleA_h2D = {sampleA_q_pt_gen, sampleA_q_pt_det, sampleA_q_pt_gen_counts, sampleA_q_pt_det_counts};

    vector<RooUnfoldResponse*> sampleB_res = {sampleB_q_pt_response, sampleB_q_pt_response_counts, sampleB_q_response, sampleB_pt_response};
    vector<TH1D*> sampleB_h1D = {sampleB_pt_gen, sampleB_pt_det, sampleB_pt_gen_counts, sampleB_pt_det_counts, sampleB_q_gen, sampleB_q_det, sampleB_q_gen_counts, sampleB_q_det_counts};
    vector<TH2D*> sampleB_h2D = {sampleB_q_pt_gen, sampleB_q_pt_det, sampleB_q_pt_gen_counts, sampleB_q_pt_det_counts};


//    cout << "STARTING SYSTEMATICS LOOP\n";

    int nSources = 7;
    
    string tree_name[nSources] = {"jetChargeTree", "towerscale", "trackingeff", "hadroncorr50", "detsmear", "gensmear", "jetChargeTree"};
    string syst_name[nSources] = {"nom", "TS", "TU", "HC50", "DS", "GS", "QS"}; // Q smear should match nominal, but is split out into 1D responses so i will come back to this 
    string lev[2] = {"part", "det"};

    TH1D* syst_dists[nSources][nJetPtBins][2];
    for(int i = 0; i < nSources; i++){
        for(int j = 0; j < nJetPtBins; j++){
            syst_dists[i][j][0] = new TH1D( Form( (lev[0] + "_q_" + syst_name[i] + "_%1.0f%1.0f").c_str(), jetptEdges[j], jetptEdges[j+1] ), "", njetQbins, -jetQedge, jetQedge );
            syst_dists[i][j][1] = new TH1D( Form( (lev[1] + "_q_" + syst_name[i] + "_%1.0f%1.0f").c_str(), jetptEdges[j], jetptEdges[j+1] ), "", njetQbins, -jetQedge, jetQedge );
        }
    }

    
    for(int iSyst = 0; iSyst < nSources; iSyst++){
        
        vector<double> *det_jetpt = 0;
        vector<double> *part_jetpt = 0;

        vector<double> *det_jetmass = 0;
        vector<double> *part_jetmass = 0;
        vector<vector<double> > *det_conspt = 0; vector<vector<double> > *det_chcons = 0;
        vector<vector<double> > *part_conspt = 0; vector<vector<double> > *part_chcons = 0;
        
        vector<double> *miss_jetpt = 0;
        vector<vector<double> > *miss_conspt = 0; vector<vector<double> > *miss_chcons = 0;
        vector<double> *fake_jetpt = 0;
        vector<vector<double> > *fake_conspt = 0; vector<vector<double> > *fake_chcons = 0;

        vector<double> *miss_jetmass = 0;
        vector<double> *fake_jetmass = 0;

    
        double weight = 1;    

        TTree *t = (TTree*) fin->Get( (tree_name[iSyst]).c_str() ); // should be the same name of the tree for all data, pythia6, pythia8
        
        t->SetBranchAddress("det_jetpt", &det_jetpt);
        t->SetBranchAddress("part_jetpt", &part_jetpt);
        t->SetBranchAddress("det_conspt", &det_conspt);
        t->SetBranchAddress("part_conspt", &part_conspt);
        t->SetBranchAddress("det_chcons", &det_chcons);
        t->SetBranchAddress("part_chcons", &part_chcons);

        if(iSyst == 0){
            t->SetBranchAddress("part_jetmass", &part_jetmass);
            t->SetBranchAddress("det_jetmass", &det_jetmass);
            t->SetBranchAddress("miss_jetmass", &miss_jetmass);
            t->SetBranchAddress("fake_jetmass", &fake_jetmass);
        }
        
        t->SetBranchAddress("miss_jetpt", &miss_jetpt);
        t->SetBranchAddress("miss_conspt", &miss_conspt);
        t->SetBranchAddress("miss_chcons", &miss_chcons);

        
        t->SetBranchAddress("fake_jetpt", &fake_jetpt);
        t->SetBranchAddress("fake_conspt", &fake_conspt);
        t->SetBranchAddress("fake_chcons", &fake_chcons);
        
        t->SetBranchAddress("weight", &weight);


        cout << "starting loop over entries\n";
    
        for(int i = 0; i < (int) t->GetEntries(); i++){
            // take each event
            t->GetEntry(i);
            double clos_rand = gRandom->Uniform(0.0, 1.0); //randomly split events into sampleA, sampleB

            for(int j = 0; j < (int) det_jetpt->size(); j++){
                // loop over matched jets and fill response
                double jc_p = 0.0; double jc_d = 0.0;
                double p_pt = part_jetpt->at(j);
                double d_pt =  det_jetpt->at(j);

                
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
                if(iSyst == 6){
                    double prior_adjust_h7 = 0.0;
                    double prior_adjust_p8 = 0.0;

                    int ji = -1;
                    if(20 < p_pt && p_pt < 25){ji = 0;}
                    else if(25 < p_pt && p_pt < 30){ji = 1;}
                    else if(30 < p_pt && p_pt < 40){ji = 2;}

                    if(ji != -1){
                        prior_adjust_h7 = qh7[ji]->GetBinContent( qh7[ji]->GetXaxis()->FindBin(jc_p) );
                        prior_adjust_p8 = qp8[ji]->GetBinContent( qp8[ji]->GetXaxis()->FindBin(jc_p) );
                    }

                    for(int ij = 0; ij < nJetPtBins; ij++){
                        if( jetptEdges[ij] < p_pt && p_pt < jetptEdges[ij+1] ){
                            syst_dists[iSyst][ij][0]->Fill(jc_p, weight * prior_adjust_h7);
                        }
                    }
                    for(int ij = 0; ij < nJetPtBins; ij++){ // loop over jet bins
                        if( jetptEdges[ij] < d_pt && d_pt < jetptEdges[ij+1] ){ // if jet pt within jet pt bin, fill hist
                            syst_dists[iSyst][ij][1]->Fill(jc_d, weight * prior_adjust_h7);
                        }
                    }


                    if( 20 < d_pt && d_pt <25 && ji == 0 ){
                        q_res2025_nom->Fill(jc_d, jc_p, weight);
                        q_res2025_h7smear->Fill(jc_d, jc_p, weight * prior_adjust_h7);
                        q_res2025_p8smear->Fill(jc_d, jc_p, weight * prior_adjust_p8);
                    }
                    else if( 25 < d_pt && d_pt < 30 && ji == 1 ){
                        q_res2530_nom->Fill(jc_d, jc_p, weight);
                        q_res2530_h7smear->Fill(jc_d, jc_p, weight *   prior_adjust_h7);
                        q_res2530_p8smear->Fill(jc_d, jc_p, weight * prior_adjust_p8);
                    }
                    else if( 30 < d_pt && d_pt < 40 && ji == 2 ){
                        q_res3040_nom->Fill(jc_d, jc_p, weight);
                        q_res3040_h7smear->Fill(jc_d, jc_p, weight * prior_adjust_h7);
                        q_res3040_p8smear->Fill(jc_d, jc_p, weight * prior_adjust_p8);
                    }
                }


                // for systematic loop
                else{
                    for(int ij = 0; ij < nJetPtBins; ij++){
                        if( jetptEdges[ij] < p_pt && p_pt < jetptEdges[ij+1] ){
                            syst_dists[iSyst][ij][0]->Fill(jc_p, weight);
                        }
                    }
                    for(int ij = 0; ij < nJetPtBins; ij++){ // loop over jet bins
                        if( jetptEdges[ij] < d_pt && d_pt < jetptEdges[ij+1] ){ // if jet pt within jet pt bin, fill hist
                            syst_dists[iSyst][ij][1]->Fill(jc_d, weight);
                        }
                    }
                    syst_res[iSyst]->Fill(jc_d, d_pt, jc_p, p_pt, weight);
                }

                if(iSyst == 0){
                    match_plus_miss->Fill(p_pt, weight);
                    q_pt_response->Fill(jc_d, d_pt, jc_p, p_pt, weight);
                    q_pt_response_counts->Fill(jc_d, d_pt, jc_p, p_pt);
                    pt_response->Fill(d_pt, p_pt, weight); // fill jet pt response with matched jets

                    // add mass responses here
                    m_pt_response->Fill( det_jetmass->at(j), d_pt, part_jetmass->at(j), p_pt, weight );
                    m_response->Fill( det_jetmass->at(j), part_jetmass->at(j), weight );

            
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


                        // add mass responses for sample A here
                        sampleA_m_gen->Fill( part_jetmass->at(j), weight );
                        sampleA_m_det->Fill( det_jetmass->at(j), weight );

                        sampleA_m_pt_gen->Fill( part_jetmass->at(j), p_pt, weight );
                        sampleA_m_pt_det->Fill( det_jetmass->at(j), d_pt, weight );

                        sampleA_m_pt_response->Fill( det_jetmass->at(j), d_pt, part_jetmass->at(j), p_pt, weight );
                        sampleA_m_response->Fill( det_jetmass->at(j), part_jetmass->at(j), weight );
                    }
                    else if(clos_rand < 0.5){
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


                        // add mass responses for sample B here

                        sampleB_m_gen->Fill( part_jetmass->at(j), weight );
                        sampleB_m_det->Fill( det_jetmass->at(j), weight );

                        sampleB_m_pt_gen->Fill( part_jetmass->at(j), p_pt, weight );
                        sampleB_m_pt_det->Fill( det_jetmass->at(j), d_pt, weight );

                        sampleB_m_pt_response->Fill( det_jetmass->at(j), d_pt, part_jetmass->at(j), p_pt, weight );
                        sampleB_m_response->Fill( det_jetmass->at(j), part_jetmass->at(j), weight );
                    }
                

/*
                    // for 1D unfolding: matched jets that have both particle and detector level in pT range, use ->Fill()
                    // "detector level jets that match with particle level jet with jet pt in pt range but detector jet pt falls outside pT range go into fakes"
                    // "particle level jets that match with detector level jet with jet pt in pt range but particle jet pt falls outside pT range go into misses"
                    if(15.0 < d_pt && d_pt < 20.0 && 15.0 < p_pt && p_pt < 20.0){ // am I requiring both detector and particle level jet pt or only one?
                        q_res1520_nom->Fill(jc_d, jc_p, weight);
                    }
                    else if(15.0 < p_pt && p_pt < 20.0){
                        q_res1520_nom->Miss(jc_p, weight);
                    }
                    else if(15.0 < d_pt && d_pt < 20.0){
                        q_res1520_nom->Fake(jc_d, weight);
                    }

                    if(20.0 < d_pt && d_pt < 25.0 && 20.0 < p_pt && p_pt < 25.0){ // am I requiring both detector and particle level jet pt or only one?
                        q_res2025_nom->Fill(jc_d, jc_p, weight);
                    }
                    else if(20.0 < p_pt && p_pt < 25.0){
                        q_res2025_nom->Miss(jc_p, weight);
                    }
                    else if(20.0 < d_pt && d_pt < 25.0){
                        q_res2025_nom->Fake(jc_d, weight);
                    }

                    if(25.0 < d_pt && d_pt < 30.0 && 25.0 < p_pt && p_pt < 30.0){ // am I requiring both detector and particle level jet pt or only one?
                        q_res2530_nom->Fill(jc_d, jc_p, weight);
                    }
                    else if(25.0 < p_pt && p_pt < 30.0){
                        q_res2530_nom->Miss(jc_p, weight);
                    }
                    else if(25.0 < d_pt && d_pt < 30.0){
                        q_res2530_nom->Fake(jc_d, weight);
                    }

                    if(30.0 < d_pt && d_pt < 40.0 && 30.0 < p_pt && p_pt < 40.0){ // am I requiring both detector and particle level jet pt or only one?
                        q_res3040_nom->Fill(jc_d, jc_p, weight);
                    }
                    else if(30.0 < p_pt && p_pt < 40.0){
                        q_res3040_nom->Miss(jc_p, weight);
                    }
                    else if(30.0 < d_pt && d_pt < 40.0){
                        q_res3040_nom->Fake(jc_d, weight);
                    }
*/
                }
            }

            for(int jj = 0; jj < (int) miss_jetpt->size(); jj++){
                // fill response misses with misses from particle level
                double jc_m = 0.0;
                double m_pt = miss_jetpt->at(jj);

                for(int ii = 0; ii < (int) miss_conspt->at(jj).size(); ii++){
                    double m_cons_pt = miss_conspt->at(jj).at(ii);
                    jc_m += pow( (m_cons_pt/m_pt), k ) * miss_chcons->at(jj).at(ii);
                }

                if(iSyst == 6){
                    double prior_adjust_h7 = 0.0;
                    double prior_adjust_p8 = 0.0;

                    int ji = -1;
                    if(20 < m_pt && m_pt < 25){ji = 0;}
                    else if(25 < m_pt && m_pt < 30){ji = 1;}
                    else if(30 < m_pt && m_pt < 40){ji = 2;}

                    if(ji != -1){
                        prior_adjust_h7 = qh7[ji]->GetBinContent( qh7[ji]->GetXaxis()->FindBin(jc_m) );
                        prior_adjust_p8 = qp8[ji]->GetBinContent( qp8[ji]->GetXaxis()->FindBin(jc_m) );
                    }

                    if( ji == 0 ){
                        q_res2025_nom->Miss(jc_m, weight);
                        q_res2025_h7smear->Miss(jc_m, weight * prior_adjust_h7);
                        q_res2025_p8smear->Miss(jc_m, weight * prior_adjust_p8);
                    }
                    else if( ji == 1 ){
                        q_res2530_nom->Miss(jc_m, weight);
                        q_res2530_h7smear->Miss(jc_m, weight * prior_adjust_h7);
                        q_res2530_p8smear->Miss(jc_m, weight * prior_adjust_p8);
                    }
                    else if( ji == 2 ){
                        q_res3040_nom->Miss(jc_m, weight);
                        q_res3040_h7smear->Miss(jc_m, weight * prior_adjust_h7);
                        q_res3040_p8smear->Miss(jc_m, weight * prior_adjust_p8);
                    }
                }

                else{
                    syst_res[iSyst]->Miss(jc_m, m_pt, weight);
                }
                if(iSyst == 0){
                    match_plus_miss->Fill(m_pt, weight);
                    
                    q_pt_response->Miss(jc_m, m_pt, weight);
                    q_pt_response_counts->Miss(jc_m, m_pt);
                    pt_response->Miss(m_pt, weight);
                
                    m_pt_response->Miss( miss_jetmass->at(jj), m_pt, weight );
                    m_response->Miss( miss_jetmass->at(jj), weight );
		    
                    if(clos_rand > 0.5){
                        sampleA_pt_gen->Fill(m_pt, weight);
                        sampleA_q_gen->Fill(jc_m, weight);
                    
                        sampleA_q_pt_gen->Fill(jc_m, m_pt, weight);
                        sampleA_q_pt_gen_counts->Fill(jc_m, m_pt);

                        sampleA_pt_response->Miss(m_pt, weight);
                        sampleA_q_response->Miss(jc_m, weight);
                        sampleA_q_pt_response->Miss(jc_m, m_pt, weight);


                        // add mass responses for sample A here
                        sampleA_m_pt_response->Miss( miss_jetmass->at(jj), m_pt, weight );
                        sampleA_m_response->Miss( miss_jetmass->at(jj), weight );

                        sampleA_m_gen->Fill( miss_jetmass->at(jj), weight );
                        sampleA_m_pt_gen->Fill( miss_jetmass->at(jj), m_pt, weight );
                    }
                    else if(clos_rand < 0.5){
                        sampleB_pt_gen->Fill(m_pt, weight);
                        sampleB_q_gen->Fill(jc_m, weight);

                        sampleB_q_pt_gen->Fill(jc_m, m_pt, weight);
                        sampleB_q_pt_gen_counts->Fill(jc_m, m_pt);

                        sampleB_pt_response->Miss(m_pt, weight);
                        sampleB_q_response->Miss(jc_m, weight);
                        sampleB_q_pt_response->Miss(jc_m, m_pt, weight);


                        // add mass responses for sample B here
                        sampleB_m_pt_response->Miss( miss_jetmass->at(jj), m_pt, weight );
                        sampleB_m_response->Miss( miss_jetmass->at(jj), weight );
                    
                        sampleB_m_gen->Fill( miss_jetmass->at(jj), weight );
                        sampleB_m_pt_gen->Fill( miss_jetmass->at(jj), m_pt, weight );
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

                if(iSyst == 6){
                    double prior_adjust_h7 = 0.0;
                    double prior_adjust_p8 = 0.0;

                    int ji = -1;
                    if(20 < f_pt && f_pt < 25){ji = 0;}
                    else if(25 < f_pt && f_pt < 30){ji = 1;}
                    else if(30 < f_pt && f_pt < 40){ji = 2;}
                
                    if(ji != -1){
                        prior_adjust_h7 = qh7[ji]->GetBinContent( qh7[ji]->GetXaxis()->FindBin(jc_f) );
                        prior_adjust_p8 = qp8[ji]->GetBinContent( qp8[ji]->GetXaxis()->FindBin(jc_f) );
                    }

                    if( ji == 0 ){
                        q_res2025_nom->Fake(jc_f, weight);
                        q_res2025_h7smear->Fake(jc_f, weight * prior_adjust_h7);
                        q_res2025_p8smear->Fake(jc_f, weight * prior_adjust_p8);
                    }
                    else if( ji == 1 ){
                        q_res2530_nom->Fake(jc_f, weight);
                        q_res2530_h7smear->Fake(jc_f, weight * prior_adjust_h7);
                        q_res2530_p8smear->Fake(jc_f, weight * prior_adjust_p8);
                    }
                    else if( ji == 2 ){
                        q_res3040_nom->Fake(jc_f, weight);
                        q_res3040_h7smear->Fake(jc_f, weight * prior_adjust_h7);
                        q_res3040_p8smear->Fake(jc_f, weight * prior_adjust_p8);
                    }
                }

                else{
                    syst_res[iSyst]->Fake(jc_f, f_pt, weight);
                }
                if(iSyst == 0){
                    q_pt_response->Fake(jc_f, f_pt, weight);
                    q_pt_response_counts->Fake(jc_f, f_pt);
                    pt_response->Fake(f_pt, weight);
                    
                    m_pt_response->Fake( fake_jetmass->at(jj), f_pt, weight );
                    m_response->Fake( fake_jetmass->at(jj), weight );
                    
                    if(clos_rand > 0.5){
                        sampleA_q_pt_det->Fill(jc_f, f_pt, weight);
                        sampleA_pt_det->Fill(f_pt, weight);

                        sampleA_q_pt_det_counts->Fill(jc_f, f_pt);

                        sampleA_pt_response->Fake(f_pt, weight);
                        sampleA_q_response->Fake(jc_f, weight);
                        sampleA_q_pt_response->Fake(jc_f, f_pt, weight);


                        // add mass responses for sample A here
                        sampleA_m_pt_response->Fake( fake_jetmass->at(jj), f_pt, weight );
                        sampleA_m_response->Fake( fake_jetmass->at(jj), weight );

                        sampleA_m_det->Fill( fake_jetmass->at(jj), weight );
                        sampleA_m_pt_det->Fill( fake_jetmass->at(jj), f_pt, weight );
                    }
                    else if(clos_rand < 0.5){
                        sampleB_q_pt_det->Fill(jc_f, f_pt, weight);
                        sampleB_pt_det->Fill(f_pt, weight);

                        sampleB_q_pt_det_counts->Fill(jc_f, f_pt);
                    
                        sampleB_pt_response->Fake(f_pt, weight);
                        sampleB_q_response->Fake(jc_f, weight);
                        sampleB_q_pt_response->Fake(jc_f, f_pt, weight);


                        // add mass responses for sample B here
                        sampleB_m_pt_response->Fake( fake_jetmass->at(jj), f_pt, weight );
                        sampleB_m_response->Fake( fake_jetmass->at(jj), weight );

                        sampleB_m_det->Fill( fake_jetmass->at(jj), weight );
                        sampleB_m_pt_det->Fill( fake_jetmass->at(jj), f_pt, weight );
                    }
                }
            }
        }

        t->ResetBranchAddresses();
    
    } // end systematic loop
    
    //    cout << "is the bin edge still where it should be?: " << partQ->GetBinLowEdge(1) << "\n";
    
    TFile *fout = new TFile(((string) argv[1]+(string) argv[2]).c_str(),"RECREATE");
    cout << "DEBUG: output file name is " << fout->GetName() << endl;

/*
    for(int i = 0; i < nSources; i++){
        for(int j = 0; j < nJetPtBins; j++){
            syst_dists[i][j][0]->Write();
            syst_dists[i][j][1]->Write();
        }
    }
*/

    
    m_pt_response->Write();
    //m_pt_response_counts->Write();
    m_response->Write();
    //pt_response->Write();
    match_plus_miss->Write();


    q_res2025_nom->Write();
    q_res2530_nom->Write();
    q_res3040_nom->Write();

    q_res2025_h7smear->Write();
    q_res2530_h7smear->Write();
    q_res3040_h7smear->Write();

    q_res2025_p8smear->Write();
    q_res2530_p8smear->Write();
    q_res3040_p8smear->Write();


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
    
//    m_pt_response->Write();
    sampleA_m_response->Write();
    sampleA_m_pt_response->Write();
    sampleB_m_response->Write();
    sampleB_m_pt_response->Write();


    sampleA_m_gen->Write();
    sampleA_m_det->Write();
    sampleB_m_gen->Write();
    sampleB_m_det->Write();

    sampleA_m_pt_gen->Write();
    sampleA_m_pt_det->Write();
    sampleB_m_pt_gen->Write();
    sampleB_m_pt_det->Write();



    cout << "Wrote to " << fout->GetName() << endl;
    
    //closing file
    fout->Close();
    cout << "Closed " << fout->GetName() << endl;
    
    return 0;
}



