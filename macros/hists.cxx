//
//  test_hists.cxx
//  
//
//  Created by Grant McNamara on 12/7/22.
//
//  hists.cxx but without the filling of histograms separated out into functions, just done in the main function all together
//

#include <stdio.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <math.h>

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
#include "Headers/funcs.hh"
#include "Headers/JCparameters.hh"


using namespace std;

using namespace jcAnalysis;



double str_to_double (std::string str) {
    std::string Num = str.substr(0,1)+"."+str.substr(1,1); //e.g. 0.4
    double str_double = (double) std::stod(Num); //converting the string 0.4 to a double
    std::cout << "DEBUG: string to number = " << str_double << std::endl;

    return str_double;
}


int main(int argc, const char ** argv){
    if(argc != 6){
        // $outlocation $outname $infiles $inputtype $kappa
        cout << "Wrong number of arguments!\n Exiting.\n";
        exit(1);
    }
    // arg1: output file location
    // arg2: output file name
    // arg3: infile name
    // arg4: data vs sim
    // arg5: kappa (as string)

    

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    string fin_name = (string) argv[3];
    TFile *fin = new TFile(fin_name.c_str(),"READ");


    string data_string = (string) argv[4];
    bool data_bool = 1;
    if(data_string == "data"){data_bool = 1;}
    else{data_bool = 0;}


    bool match = 0;
    if(fin_name.find("_matched") != string::npos){match = 1;}
    
    
    const int nJetPtBins = 4; // need to match the binning of the detector level in the response
    vector<double> jetPtLo = {15.0, 20.0, 25.0, 30.0};
    vector<double> jetPtHi = {20.0, 25.0, 30.0, 40.0};

    double jetptEdges[nJetPtBins+1] = {15.0, 20.0, 25.0, 30.0, 40.0};
    
    int jetptHiEdge = 50;
    // to use narrower bins in jet pt on det-level, data
    
    
    const string k = (string) argv[5];
    double kappa = str_to_double(k);

    TString k_str(k);

    // to adapt to trying k < 0, kappa argument turns into absolute value, extra argument to designate if not positive
    int k_s = 1;
    
    
    double qBinWidth = 0.1; //default to ~ CMS binning
    // default to non-k=0.0 bins
    double jetQedge = 4.5 + (1.0/2.0)*qBinWidth; // pre-11/14/22: using 18 bins in [-3.5, 3.5] gives bin width of ~ 7/18 ~= 0.39 [e]
    // what if i use bin width comparable to CMS' choice? ~ 21 bins in [-1, 1] which is roughly 2/21 ~= 0.1 [e]
    if(kappa == 0.0){qBinWidth = 1.0; jetQedge = 6.0 + (1.0/2.0)*qBinWidth;} // no real choice in binning for k = 0.0 since k = 0.0 corresponds to adding integers together
    if(kappa < 0.0){jetQedge = 10.0 + (1.0/2.0)*qBinWidth;} // bin width = 20/100 = 0.2, can even get finer bins; not even as fine as CMS' yet (jet to study effect of bin width in fitting)

    int njetQbins = (2 * jetQedge)/qBinWidth; //
    
    
    int nresbins = 40; int resedge = 5;
    if(kappa == 0.0){nresbins = 11;}
    
    
    int nchbins = 11;

    
    /// START MAKING HISTOGRAMS
    
    // pT spectra for data or det- and part-level pythia
    TH1D* njets_pt = new TH1D("jet_pt_hist", "", 60, 0.0, 60.0);
    TH1D* njets_pt_pl = new TH1D("jet_pt_hist_part", "", 60, 0.0, 60.0);
    
    // number of charged constituents, 2D with jet pt
    TH2D* nch_cons =  new TH2D("nch_cons_d", "", nchbins, 0, nchbins - 0.5, 60, 0.0, 60.0);
    TH2D* nch_cons_pl =  new TH2D("nch_cons_p", "", nchbins, -0.5, nchbins - 0.5, 60, 0.0, 60.0);

    // something like the resolution for pythia-6
    TH2D* deltaPtvGePt = new TH2D("deltaPtvGePt", "", 11, 5, 60, 220, -6, 6);
    TH2D* deltaPtvPyPt = new TH2D("deltaPtvPyPt", "", 11, 5, 60, 220, -1, 1);
    
    
    // ( detector level jet Q - particle level jet Q ) / particle level jet Q
    // feels more logical to me to put pythia pt on x-axis
    TH2D* deltaQvPyPt = new TH2D("deltaQvPyPt", ";Gen. p^{jet}_{T};#DeltaQ_{#kappa} (Det - Gen) / Gen. #DeltaQ_{#kappa}", 11, 5, 60, 100, -1, 1);
    TH2D* deltaQvGePt = new TH2D("deltaQvGePt", ";Det. p^{jet}_{T};#DeltaQ_{#kappa} (Det - Gen) / Gen. #DeltaQ_{#kappa}", 11, 5, 60, 100, -6, 6);
    // axis titles following Isaac's convention in jetmass2/macros/hists.cxx
    
    
    TH2D* QvPt = new TH2D("QvPt_d", "", njetQbins, -jetQedge, jetQedge, 9, 5, 60); // data/detector level
    TH2D* QvPt_p = new TH2D("QvPt_p", "", njetQbins, -jetQedge, jetQedge, 15, 5, 80); // particle level (should not be filled for data)
    
    
    
    TH1D* part_posQresolution[nJetPtBins];
    TH1D* part_negQresolution[nJetPtBins];
    TH1D* det_posQresolution[nJetPtBins];
    TH1D* det_negQresolution[nJetPtBins];
    
    
    TH1D* jetQ[nJetPtBins];// = new TH1D(Form("kappa_k%s", k), "", 9, -4.5, 4.5);
    TH1D* jetQ_part[nJetPtBins];// = new TH1D(Form("kappa_part_k%s", k), "", 9, -4.5, 4.5);

    TH2D* q_v_pt_counts[2];
    TH1D* q_counts[2];
    TH1D* pt_counts[2];

    

    for(int ij = 0; ij < nJetPtBins; ij++){
        
        part_posQresolution[ij] = new TH1D( Form( "part_posqResolution_jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
        part_negQresolution[ij] = new TH1D( Form( "part_negqResolution_jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
        det_posQresolution[ij] = new TH1D( Form( "det_posqResolution_jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
        det_negQresolution[ij] = new TH1D( Form( "det_negqResolution_jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
        jetQ[ij] = new TH1D( Form( "jetQ_jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
        jetQ_part[ij] = new TH1D( Form( "jetQ_part_jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
    }
    
    
    
    // [0] is det-level, [1] is part-level; not used in data
    q_v_pt_counts[0] = new TH2D("q_v_pt_counts", "", njetQbins, -jetQedge, jetQedge, 9, 15, 60);
    q_counts[0] = new TH1D("q_counts", "", njetQbins, -jetQedge, jetQedge);
    pt_counts[0] = new TH1D("pt_counts", "", 9, 15, 60);
    if(!data_bool){ // if pythia, need part-level histograms in addition to det-level/data
        q_v_pt_counts[1] = new TH2D("PL_q_v_pt_counts", "", njetQbins, -jetQedge, jetQedge, 15, 5, 80);
        q_counts[1] = new TH1D("PL_q_counts", "", njetQbins, -jetQedge, jetQedge);
        pt_counts[1] = new TH1D("PL_pt_counts", "", 15, 5, 80);
    }
    // HISTOGRAMS ARE MADE
    
    
    //START FILLING HISTOGRAMS
    
    TTree *t = (TTree*) fin->Get("jetChargeTree");
    
    if(data_bool){
        vector<double> *jetpt = 0;
        vector<vector<double> > *conspt = 0; vector<vector<double> > *chcons = 0;
        //double weight = 1.0;
        
        t->SetBranchAddress("jetpt", &jetpt);
        
        t->SetBranchAddress("conspt", &conspt);
        t->SetBranchAddress("chcons", &chcons);
        
        for(int i = 0; i < (int) t->GetEntries(); i++){
            t->GetEntry(i);
            for(int ij = 0; ij < (int) jetpt->size(); ij++){
                
                double jc = 0.0;
                double jpt = jetpt->at(ij);
                int nch = 0;
                
                njets_pt->Fill(jpt);
                
                for(int j = 0; j < (int) conspt->at(ij).size(); j++){
                    double pt_i = conspt->at(ij).at(j);
                    double charge = chcons->at(ij).at(j);
                    
                    jc += pow( (pt_i/jpt), kappa ) * charge;

                }
                // q vs pt 2D
                QvPt->Fill(jc, jpt); // passed a vector of hists (should be max size 2) and for data only fill the first in the vector
                
                // number of charged constituents in jet
                nch_cons->Fill( nch, jpt );
                
                
                for(int kj = 0; kj < nJetPtBins; kj++){ // loop through jet pt bins
                    if( jetPtLo[kj] < jpt && jpt < jetPtHi[kj] ){
                        jetQ[kj]->Fill(jc);
                    }
                }
            }
        }
        t->ResetBranchAddresses();
    }
    else if(match){
        vector<double> *part_jetpt = 0;
        vector<double> *det_jetpt = 0;
        vector<vector<double> > *part_conspt = 0; vector<vector<double> > *part_chcons = 0;
        vector<vector<double> > *det_conspt = 0;  vector<vector<double> > *det_chcons = 0;
        
        double weight = 1.0;
        
        t->SetBranchAddress("part_jetpt", &part_jetpt);
        t->SetBranchAddress("det_jetpt", &det_jetpt);
        
        t->SetBranchAddress("part_conspt", &part_conspt);
        t->SetBranchAddress("det_conspt", &det_conspt);
        t->SetBranchAddress("part_chcons", &part_chcons);
        t->SetBranchAddress("det_chcons", &det_chcons);
        
        t->SetBranchAddress("weight", &weight);
        
        for(int i = 0; i < (int) t->GetEntries(); i++){
            t->GetEntry(i);
            for(int ij = 0; ij < (int) part_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = part_jetpt->at(ij);
                int nch = 0;
                
                njets_pt_pl->Fill(jpt, weight);
                
                for(int j = 0; j < (int) part_conspt->at(ij).size(); j++){
                    double pt_i = part_conspt->at(ij).at(j);
                    double charge = part_chcons->at(ij).at(j);
                    
                    
                    jc += pow( (pt_i/jpt), kappa ) * charge;
                } // part-level jet constituents

                // q vs pt with weights
                QvPt_p->Fill(jc, jpt, weight);

                // q vs pt without weights
                q_v_pt_counts[1]->Fill(jc, jpt);

                // number of charged constituents in jet
                nch_cons_pl->Fill( nch, jpt, weight );
                
                // 1D jet Q hist: find bin corresponding to jet's pt, fill histogram with jet charge
                for(int kj = 0; kj < nJetPtBins; kj++){ // loop through jet pt bins
                    if( jetPtLo[kj] < jpt && jpt < jetPtHi[kj] ){
                        jetQ_part[kj]->Fill(jc, weight);
                        
                        if( jc > 0.0 ){
                            part_posQresolution[kj]->Fill(jc, weight);
                        }
                        else if( jc < 0.0 ){
                            part_negQresolution[kj]->Fill(jc, weight);
                        }
                    }
                }
                
                // 1D counts
                q_counts[1]->Fill(jc);
                pt_counts[1]->Fill(jpt);
                
                // jet pt resolution: delta pt/part pt vs det-/part- pt
                deltaPtvGePt->Fill( det_jetpt->at(ij), ( det_jetpt->at(ij) - part_jetpt->at(ij) ) / (double) part_jetpt->at(ij), weight );
                deltaPtvPyPt->Fill( part_jetpt->at(ij), ( det_jetpt->at(ij) - part_jetpt->at(ij) ) / (double) part_jetpt->at(ij), weight );
                
                
            } // part-level jets
            for(int ij = 0; ij < (int) det_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = det_jetpt->at(ij);
                int nch = 0;
                
                njets_pt->Fill(jpt, weight);
                
                for(int j = 0; j < (int) det_conspt->at(ij).size(); j++){
                    double pt_i = det_conspt->at(ij).at(j);
                    double charge = det_chcons->at(ij).at(j);
                    
                    
                    jc += pow( (pt_i/jpt), kappa ) * charge;
                } // det-level jet constituents

                QvPt->Fill(jc, jpt, weight);
                q_v_pt_counts[0]->Fill(jc, jpt);

                // number of charged constituents in jet
                nch_cons->Fill( nch, jpt, weight );
                
                // 1D jet Q hist: find bin corresponding to jet's pt, fill histogram with jet charge
                for(int kj = 0; kj < nJetPtBins; kj++){ // loop through jet pt bins
                    if( jetPtLo[kj] < jpt && jpt < jetPtHi[kj] ){
                        jetQ[kj]->Fill(jc, weight);
                        
                        if( jc > 0.0 ){
                            det_posQresolution[kj]->Fill(jc, weight);
                        }
                        else if( jc < 0.0 ){
                            det_negQresolution[kj]->Fill(jc, weight);
                        }
                    }
                }
                
                // 1D counts
                q_counts[0]->Fill(jc);
                pt_counts[0]->Fill(jpt);
                
            } // det-level jets
        }// entries
        t->ResetBranchAddresses();
    }// if match
    else if(!match){
        vector<double> *part_jetpt = 0;
        vector<vector<double> > *part_conspt = 0; vector<vector<double> > *part_chcons = 0;
        vector<double> *det_jetpt = 0;
        vector<vector<double> > *det_conspt = 0; vector<vector<double> > *det_chcons = 0;
        
        double weight = 1.0;
        
        
        t->SetBranchAddress("part_jetpt", &part_jetpt);
        t->SetBranchAddress("det_jetpt", &det_jetpt);
        
        t->SetBranchAddress("part_conspt", &part_conspt);
        t->SetBranchAddress("part_chcons", &part_chcons);
        t->SetBranchAddress("det_conspt", &det_conspt);
        t->SetBranchAddress("det_chcons", &det_chcons);
        
        t->SetBranchAddress("weight", &weight);
        
        for(int i = 0; i < (int) t->GetEntries(); i++){
            t->GetEntry(i);
            for(int ij = 0; ij < (int) part_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = part_jetpt->at(ij);
                int nch = 0;
                
                njets_pt_pl->Fill(jpt, weight);
                
                for(int j = 0; j < (int) part_conspt->at(ij).size(); j++){
                    double pt_i = part_conspt->at(ij).at(j);
                    double charge = part_chcons->at(ij).at(j);
                    
                    
                    jc += pow( (pt_i/jpt), kappa ) * charge;
                } // part-level jet constituents
                QvPt_p->Fill(jc, jpt, weight);
                q_v_pt_counts[1]->Fill(jc, jpt);
                
                // number of charged constituents in jet
                nch_cons_pl->Fill( nch, jpt, weight );
                
                // 1D jet Q hist: find bin corresponding to jet's pt, fill histogram with jet charge
                for(int kj = 0; kj < nJetPtBins; kj++){ // loop through jet pt bins
                    if( jetPtLo[kj] < jpt && jpt < jetPtHi[kj] ){
                        jetQ_part[kj]->Fill(jc, weight);
                    }
                }
                
                // 1D counts
                q_counts[1]->Fill(jc);
                pt_counts[1]->Fill(jpt);
                
            } // part-level jets
            for(int ij = 0; ij < (int) det_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = det_jetpt->at(ij);
                int nch = 0;
                
                njets_pt->Fill(jpt, weight);
                
                for(int j = 0; j < (int) det_conspt->at(ij).size(); j++){
                    double pt_i = det_conspt->at(ij).at(j);
                    double charge = det_chcons->at(ij).at(j);
                    
                    
                    jc += pow( (pt_i/jpt), kappa ) * charge;
                } // det-level jet constituents
                QvPt->Fill(jc, jpt, weight);
                q_v_pt_counts[0]->Fill(jc, jpt);
                
                // number of charged constituents in jet
                nch_cons->Fill( nch, jpt );
                
                // 1D jet Q hist: find bin corresponding to jet's pt, fill histogram with jet charge
                for(int kj = 0; kj < nJetPtBins; kj++){ // loop through jet pt bins
                    if( jetPtLo[kj] < jpt && jpt < jetPtHi[kj] ){
                        jetQ[kj]->Fill(jc, weight);
                    }
                }
                
                // 1D counts
                q_counts[0]->Fill(jc);
                pt_counts[0]->Fill(jpt);

            } // det-level jets
        } //entries
        t->ResetBranchAddresses();
    } // if not matched but pythia
    
    
    
    TFile *fout = new TFile(((string) argv[1]+(string) argv[2]).c_str(),"RECREATE");
    cout << "DEBUG: output file name is " << fout->GetName() << endl;
    
    
    // write 1D jet Q histograms, if pythia write both det- and part-level histograms
    for(int i = 0; i < nJetPtBins; i++){
        jetQ[i]->Write();
        if(!data_bool){
            jetQ_part[i]->Write();
        }
    }
    
    // write 2D q vs pt histogram, if pythia write both det- and part-level histogram
    QvPt->Write();
    if(!data_bool){
        // if matched or unmatched sim, write part-level Q vs pT 2D hist
        QvPt_p->Write();
    }
    
    
    // write jet pt spectra, if pythia write both det- and part-level spectra
    njets_pt->Write();
    if(!data_bool){
        // if matched or unmatched sim, write part-level Q vs pT 2D hist
        njets_pt_pl->Write();
    }
    
    
    nch_cons->Write();
    if(!data_bool){
        nch_cons_pl->Write();
    }
    
    
    if(match){ // for matched pythia, write resolution, delta pt, etc
        deltaQvPyPt->Write();
        deltaQvGePt->Write();

        deltaPtvGePt->Write();
        deltaPtvPyPt->Write();

        for(int j = 0; j < (int) nJetPtBins; j++){
            part_posQresolution[j]->Write();
            part_negQresolution[j]->Write();
            det_posQresolution[j]->Write();
            det_negQresolution[j]->Write();
        }
    }
    
    cout << "Wrote to " << fout->GetName() << endl;
    
    fout->Close();
    cout << "Closed " << fout->GetName() << endl;
    
    return 0;
} // main

