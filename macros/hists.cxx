//
//  hists.cxx
//  
//
//  Created by Grant on 2/2/22.
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

using namespace std;



double str_to_double (std::string str) {
    std::string Num = str.substr(0,1)+"."+str.substr(1,1); //e.g. 0.4
    double str_double = (double) std::stod(Num); //converting the string 0.4 to a double
    std::cout << "DEBUG: string to number = " << str_double << std::endl;

    return str_double;
}



void Fill2DHist(TFile* file, vector<TH2D*> hists, double kappa, bool data_bool = 1, bool match = 1){
    // fill 2D histogram needed for unfolding for pythia 6 and data
    
    
    TTree *t = (TTree*) file->Get("jetChargeTree");
    
    
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
                for(int j = 0; j < (int) conspt->at(ij).size(); j++){
                    double pt_i = conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * chcons->at(ij).at(j);
                }
                hists[0]->Fill(jc, jpt); // passed a vector of hists (should be max size 2) and for data only fill the first in the vector
            }
        }
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
                for(int j = 0; j < (int) part_conspt->at(ij).size(); j++){
                    double pt_i = part_conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * part_chcons->at(ij).at(j);
                }
                hists[0]->Fill(jc, jpt, weight); // passed a vector of hists (should be max size 2) and for particle level, fill the first in the vector
            }
            for(int ij = 0; ij < (int) det_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = det_jetpt->at(ij);
                for(int j = 0; j < (int) det_conspt->at(ij).size(); j++){
                    double pt_i = det_conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * det_chcons->at(ij).at(j);
                }
                hists[1]->Fill(jc, jpt, weight); // passed a vector of hists (should be max size 2) and for detector level, fill the second in the vector
            }
        }
    }
    else if(!match){
        vector<double> *part_jetpt = 0;
        vector<vector<double> > *part_conspt = 0; vector<vector<double> > *part_chcons = 0;
        
        double weight = 1.0;
        
        
        t->SetBranchAddress("part_jetpt", &part_jetpt);
        
        t->SetBranchAddress("part_conspt", &part_conspt);
        t->SetBranchAddress("part_chcons", &part_chcons);
        
	t->SetBranchAddress("weight", &weight);

        for(int i = 0; i < (int) t->GetEntries(); i++){
            t->GetEntry(i);
            for(int ij = 0; ij < (int) part_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = part_jetpt->at(ij);
                for(int j = 0; j < (int) part_conspt->at(ij).size(); j++){
                    double pt_i = part_conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * part_chcons->at(ij).at(j);
                }
                hists[1]->Fill(jc, jpt, weight); // passed a vector of hists (should be max size 2) and for particle level, fill the first in the vector
            }
        }
    }
    t->ResetBranchAddresses();
    return;
}




void Fill2DCounts(TFile* file, vector<TH2D*> hists, double kappa, bool data_bool = 1, bool match = 1){
    // fill 2D histogram needed for unfolding for pythia 6 and data
    
    
    TTree *t = (TTree*) file->Get("jetChargeTree");
    
    
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
                for(int j = 0; j < (int) conspt->at(ij).size(); j++){
                    double pt_i = conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * chcons->at(ij).at(j);
                }
                hists[0]->Fill(jc, jpt); // passed a vector of hists (should be max size 2) and for data only fill the first in the vector
            }
        }
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
                for(int j = 0; j < (int) part_conspt->at(ij).size(); j++){
                    double pt_i = part_conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * part_chcons->at(ij).at(j);
                }
                hists[0]->Fill(jc, jpt); // passed a vector of hists (should be max size 2) and for particle level, fill the first in the vector
            }
            for(int ij = 0; ij < (int) det_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = det_jetpt->at(ij);
                for(int j = 0; j < (int) det_conspt->at(ij).size(); j++){
                    double pt_i = det_conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * det_chcons->at(ij).at(j);
                }
                hists[1]->Fill(jc, jpt); // passed a vector of hists (should be max size 2) and for detector level, fill the second in the vector
            }
        }
    }
    else if(!match){
        vector<double> *part_jetpt = 0;
        vector<vector<double> > *part_conspt = 0; vector<vector<double> > *part_chcons = 0;
        
        double weight = 1.0;
        
        
        t->SetBranchAddress("part_jetpt", &part_jetpt);
        
        t->SetBranchAddress("part_conspt", &part_conspt);
        t->SetBranchAddress("part_chcons", &part_chcons);
        
	t->SetBranchAddress("weight", &weight);

        for(int i = 0; i < (int) t->GetEntries(); i++){
            t->GetEntry(i);
            for(int ij = 0; ij < (int) part_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = part_jetpt->at(ij);
                for(int j = 0; j < (int) part_conspt->at(ij).size(); j++){
                    double pt_i = part_conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * part_chcons->at(ij).at(j);
                }
                hists[1]->Fill(jc, jpt); // passed a vector of hists (should be max size 2) and for particle level, fill the first in the vector
            }
        }
    }
    t->ResetBranchAddresses();
    return;
}



void Fill1DCounts(TFile* file, vector<TH1D*> qhists, vector<TH1D*> pthists, double kappa, bool data_bool = 1, bool match = 1){
    // fill 1D histogram needed for unfolding for pythia 6 and data
    
    
    TTree *t = (TTree*) file->Get("jetChargeTree");
    
    
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
                for(int j = 0; j < (int) conspt->at(ij).size(); j++){
                    double pt_i = conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * chcons->at(ij).at(j);
                }
                qhists[0]->Fill(jc); // passed a vector of hists (should be max size 2) and for data only fill the first in the vector
                pthists[0]->Fill(jpt); // passed a vector of hists (should be max size 2) and for data only fill the first in the vector
            }
        }
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
                for(int j = 0; j < (int) part_conspt->at(ij).size(); j++){
                    double pt_i = part_conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * part_chcons->at(ij).at(j);
                }
                qhists[0]->Fill(jc); // passed a vector of hists (should be max size 2) and for particle level, fill the first in the vector
                pthists[0]->Fill(jpt); // passed a vector of hists (should be max size 2) and for particle level, fill the first in the vector
            }
            for(int ij = 0; ij < (int) det_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = det_jetpt->at(ij);
                for(int j = 0; j < (int) det_conspt->at(ij).size(); j++){
                    double pt_i = det_conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * det_chcons->at(ij).at(j);
                }
                qhists[1]->Fill(jc); // passed a vector of hists (should be max size 2) and for detector level, fill the second in the vector
                pthists[1]->Fill(jpt); // passed a vector of hists (should be max size 2) and for detector level, fill the second in the vector
            }
        }
    }
    else if(!match){
        vector<double> *part_jetpt = 0;
        vector<vector<double> > *part_conspt = 0; vector<vector<double> > *part_chcons = 0;
        
        double weight = 1.0;
        
        
        t->SetBranchAddress("part_jetpt", &part_jetpt);
        
        t->SetBranchAddress("part_conspt", &part_conspt);
        t->SetBranchAddress("part_chcons", &part_chcons);
        
	t->SetBranchAddress("weight", &weight);

        for(int i = 0; i < (int) t->GetEntries(); i++){
            t->GetEntry(i);
            for(int ij = 0; ij < (int) part_jetpt->size(); ij++){
                double jc = 0.0;
                double jpt = part_jetpt->at(ij);
                for(int j = 0; j < (int) part_conspt->at(ij).size(); j++){
                    double pt_i = part_conspt->at(ij).at(j);
                    jc += pow( (pt_i/jpt), kappa ) * part_chcons->at(ij).at(j);
                }
                qhists[1]->Fill(jc); // passed a vector of hists (should be max size 2) and for particle level, fill the first in the vector
                pthists[1]->Fill(jpt); // passed a vector of hists (should be max size 2) and for particle level, fill the first in the vector
            }
        }
    }
    t->ResetBranchAddresses();
    return;
}




void JetChargeHist(TFile* file, vector<TH1D*> hjetcharge, TH1D* hjetpt, double kappa, vector<double> jetptLo, vector<double> jetptHi, const bool dataflag = 1){
    //fill histograms for data or pythia particle OR detector level (not matched) with jet charge using provided kappa for jet pt ranges
    
    int N = jetptLo.size(); // should be the same as jetptHi.size() and hjetcharge.size() and length of njets
    
    vector<double> *jetpt = 0;
    vector<vector<double> > *conspt = 0; vector<vector<double> > *chcons = 0;
    double weight = 1.0;// double n_jets[N];
    
    
    TTree *t = (TTree*) file->Get("jetChargeTree");
    
    if(!dataflag){
	t->SetBranchAddress("weight", &weight);

        t->SetBranchAddress("part_jetpt", &jetpt);
        t->SetBranchAddress("part_conspt", &conspt);
        t->SetBranchAddress("part_chcons", &chcons);
    }
    else{
        t->SetBranchAddress("jetpt", &jetpt);
        t->SetBranchAddress("conspt", &conspt);
        t->SetBranchAddress("chcons", &chcons);
    }
    
    for(int i = 0; i < (int) t->GetEntries(); i++){
        t->GetEntry(i);
        for(int j = 0; j < (int) jetpt->size(); j++){
            double jpt = jetpt->at(j);
            hjetpt->Fill(jpt, weight);

            for(int n = 0; n < N; n++){
                if(jetptLo[n] < jpt && jpt < jetptHi[n]){// if jet pt is within n-th bin, fill n-th histogram
                    double jc = 0.0;
                    for(int k = 0; k < (int) conspt->at(j).size(); k++){
                        double z = conspt->at(j).at(k)/jpt;
                        jc += pow(z, kappa)*(chcons->at(j).at(k));
                    }
                    hjetcharge[n]->Fill(jc, weight);
                }
            }
        }
    }
    
    t->ResetBranchAddresses();
    return;
}



void MatchedJetChargeHist(TFile* file, vector<TH1D*> hjetcharge, vector<TH1D*> hjetcharge1, TH1D* hjetpt, TH1D* hjetpt_pl, double kappa, vector<double> jetptLo, vector<double> jetptHi){
    //fill histograms for pythia particle and detector level matched with jet charge using provided kappa for jet pt ranges
    
    int N = jetptLo.size(); // should be the same as jetptHi.size() and hjetcharge.size() and length of njets

    
    vector<double> *det_jetpt = 0; vector<double> *part_jetpt = 0;
    vector<vector<double> > *det_conspt = 0; vector<vector<double> > *part_conspt = 0;
    vector<vector<double> > *det_chcons = 0; vector<vector<double> > *part_chcons = 0;

    double weight = 1.0;// double n_jets = 0;
    
    TTree *t = (TTree*) file->Get("jetChargeTree");
    
//    if(!dataflag){t->SetBranchAddress("weight", &weight);}
    t->SetBranchAddress("det_jetpt", &det_jetpt);
    t->SetBranchAddress("part_jetpt", &part_jetpt);

    t->SetBranchAddress("det_conspt", &det_conspt);
    t->SetBranchAddress("part_conspt", &part_conspt);
    t->SetBranchAddress("det_chcons", &det_chcons);
    t->SetBranchAddress("part_chcons", &part_chcons);

    t->SetBranchAddress("weight", &weight);

    
    for(int i = 0; i < (int) t->GetEntries(); i++){
        t->GetEntry(i);
        for(int j = 0; j < (int) det_jetpt->size(); j++){
            double jpt = det_jetpt->at(j);
            double jpt_part = part_jetpt->at(j);
            hjetpt->Fill(jpt, weight);
            hjetpt_pl->Fill(jpt_part, weight);
            for(int n = 0; n < N; n++){ // loop through relevant jet pt ranges
                if(jetptLo[n] < jpt && jpt < jetptHi[n]){ // if i-th jet's pt is in n-th range, calculate jet Q and fill n-th histogram in vector

                    double jc_det = 0.0;

                    for(int k = 0; k < (int) det_conspt->at(j).size(); k++){
                        double z_cons = det_conspt->at(j).at(k)/jpt;
                        jc_det += pow(z_cons, kappa)*(det_chcons->at(j).at(k));
                    }
                    hjetcharge[n]->Fill(jc_det, weight);
		}
		
		if(jetptLo[n] < jpt_part && jpt_part < jetptHi[n]){
                    double jc_part = 0.0;

                    for(int k = 0; k < (int) part_conspt->at(j).size(); k++){
                        double z_cons = part_conspt->at(j).at(k)/jpt_part;
                        jc_part += pow(z_cons, kappa)*(part_chcons->at(j).at(k));
                    }
                    hjetcharge1[n]->Fill(jc_part, weight);
                }
            }
        }
    }
        
    t->ResetBranchAddresses();
    return;
}



void fillDeltaPtHists(TFile* file, TH2D* partHist, TH2D* detHist){

    vector<double> *det_jetpt = 0; vector<double> *part_jetpt = 0;
    vector<vector<double> > *det_conspt = 0; vector<vector<double> > *part_conspt = 0;
    vector<vector<double> > *det_chcons = 0; vector<vector<double> > *part_chcons = 0;

    double weight = 1.0;// double n_jets = 0;

    TTree *t = (TTree*) file->Get("jetChargeTree");

    t->SetBranchAddress("det_jetpt", &det_jetpt);
    t->SetBranchAddress("part_jetpt", &part_jetpt);

    t->SetBranchAddress("det_conspt", &det_conspt);
    t->SetBranchAddress("part_conspt", &part_conspt);
    t->SetBranchAddress("det_chcons", &det_chcons);
    t->SetBranchAddress("part_chcons", &part_chcons);

    t->SetBranchAddress("weight", &weight);

    for(int i = 0; i < (int) t->GetEntries(); i++){
        t->GetEntry(i);
        for(int j = 0; j < (int) part_jetpt->size(); j++){
            detHist->Fill(det_jetpt->at(j), ( det_jetpt->at(j) - part_jetpt->at(j) ) / (double) part_jetpt->at(j), weight);
            partHist->Fill(part_jetpt->at(j), ( det_jetpt->at(j) - part_jetpt->at(j) ) / (double) part_jetpt->at(j), weight);
        }
    }
    t->ResetBranchAddresses();
    return;
}


void fillJetQResolution(TFile* file, vector<TH1D*> part_pos, vector<TH1D*> part_neg, vector<TH1D*> det_pos, vector<TH1D*> det_neg, double kappa, vector<double> jetptLo, vector<double> jetptHi){
    // vectors of histograms should be [Q < 0, Q > 0] in order
    vector<double> *det_jetpt = 0; vector<double> *part_jetpt = 0;
    vector<vector<double> > *det_conspt = 0; vector<vector<double> > *part_conspt = 0;
    vector<vector<double> > *det_chcons = 0; vector<vector<double> > *part_chcons = 0;

    double weight = 1.0;// double n_jets = 0;

    TTree *t = (TTree*) file->Get("jetChargeTree");

    t->SetBranchAddress("det_jetpt", &det_jetpt);
    t->SetBranchAddress("part_jetpt", &part_jetpt);

    t->SetBranchAddress("det_conspt", &det_conspt);
    t->SetBranchAddress("part_conspt", &part_conspt);
    t->SetBranchAddress("det_chcons", &det_chcons);
    t->SetBranchAddress("part_chcons", &part_chcons);

    t->SetBranchAddress("weight", &weight);

    for(int i = 0; i < (int) t->GetEntries(); i++){
        t->GetEntry(i);
        for(int j = 0; j < (int) part_jetpt->size(); j++){
	    double p_pt = part_jetpt->at(j);
	    double d_pt = det_jetpt->at(j);

	    for(int ji = 0; ji < (int) jetptHi.size(); ji++){
		// only bother to loop over constituents if the jet is in the right pt bin
		if(jetptLo[ji] < p_pt && p_pt < jetptHi[ji]){
		    double jc_p = 0.0;
		    for(int k = 0; k < (int) part_conspt->at(j).size(); k++){
			jc_p += pow( (part_conspt->at(j).at(k)/p_pt ) , kappa) * part_chcons->at(j).at(k);
		    }
		    if(jc_p < 0.0){
			part_neg[ji]->Fill(jc_p, weight);
		    }
		    else if(jc_p > 0.0){
			part_pos[ji]->Fill(jc_p, weight);
		    }
		}
	    
		// only bother to loop over constituents if the jet is in the right pt bin
		if(jetptLo[ji] < d_pt && d_pt < jetptHi[ji]){
		    double jc_d = 0.0;
		    for(int k = 0; k < (int) det_conspt->at(j).size(); k++){
			jc_d += pow( (det_conspt->at(j).at(k)/d_pt ) , kappa) * det_chcons->at(j).at(k);
		    }
		    if(jc_d < 0.0){
			det_neg[ji]->Fill(jc_d, weight);
		    }
		    else if(jc_d > 0.0){
			det_pos[ji]->Fill(jc_d, weight);
		    }
		}
	    }
        }
    }
    t->ResetBranchAddresses();
    return;
}



void fillDeltaQHists(TFile* file, TH2D* parthist, TH2D* dethist, double kappa){
	// this function is the actual resolution (i misunderstood initially)
	// resolution id the difference of the jet charge (det - part level) / part level
	// the mean of this distribution is the "jet charge scale"
	// the width of this distribution is the "jet charge resolution"
	
	// still to consider: do i need to fill a histogram with this quantity split by positive/negative Q?
	// if yes, turn parthist, dethist into an array of size 2 of TH2Ds

    vector<double> *det_jetpt = 0; vector<double> *part_jetpt = 0;
    vector<vector<double> > *det_conspt = 0; vector<vector<double> > *part_conspt = 0;
    vector<vector<double> > *det_chcons = 0; vector<vector<double> > *part_chcons = 0;

    double weight = 1.0;// double n_jets = 0;

    TTree *t = (TTree*) file->Get("jetChargeTree");

    t->SetBranchAddress("det_jetpt", &det_jetpt);
    t->SetBranchAddress("part_jetpt", &part_jetpt);

    t->SetBranchAddress("det_conspt", &det_conspt);
    t->SetBranchAddress("part_conspt", &part_conspt);
    t->SetBranchAddress("det_chcons", &det_chcons);
    t->SetBranchAddress("part_chcons", &part_chcons);

    t->SetBranchAddress("weight", &weight);


    for(int i = 0; i < (int) t->GetEntries(); i++){
        t->GetEntry(i);
        for(int j = 0; j < (int) part_jetpt->size(); j++){
            double p_pt = part_jetpt->at(j);
            double d_pt = det_jetpt->at(j);

            // only bother to loop over constituents if the jet is in the right pt bin
            double jc_p = 0.0;
            double jc_d = 0.0;
            for(int k = 0; k < (int) part_conspt->at(j).size(); k++){
                jc_p += pow( (part_conspt->at(j).at(k)/p_pt ) , kappa) * part_chcons->at(j).at(k);
            }
            for(int k = 0; k < (int) det_conspt->at(j).size(); k++){
                jc_d += pow( (det_conspt->at(j).at(k)/d_pt ) , kappa) * det_chcons->at(j).at(k);
            }
	    parthist->Fill( p_pt, (jc_d - jc_p)/ (double) jc_p , weight );
	    dethist->Fill( d_pt, (jc_d - jc_p)/ (double) jc_p , weight );
	}
    }
    t->ResetBranchAddresses();
    return;
}



void fillDeltaQHists(TFile* file, vector<TH1D*> parthist, vector<TH1D*> dethist, double* jetptedges, double kappa){
	// this function is the actual resolution (i misunderstood initially)
	// resolution id the difference of the jet charge (det - part level) / part level
	// the mean of this distribution is the "jet charge scale"
	// the width of this distribution is the "jet charge resolution"
	
	// still to consider: do i need to fill a histogram with this quantity split by positive/negative Q?
	// if yes, turn parthist, dethist into an array of size 2 of TH2Ds

    vector<double> *det_jetpt = 0; vector<double> *part_jetpt = 0;
    vector<vector<double> > *det_conspt = 0; vector<vector<double> > *part_conspt = 0;
    vector<vector<double> > *det_chcons = 0; vector<vector<double> > *part_chcons = 0;

    double weight = 1.0;// double n_jets = 0;

    TTree *t = (TTree*) file->Get("jetChargeTree");

    t->SetBranchAddress("det_jetpt", &det_jetpt);
    t->SetBranchAddress("part_jetpt", &part_jetpt);

    t->SetBranchAddress("det_conspt", &det_conspt);
    t->SetBranchAddress("part_conspt", &part_conspt);
    t->SetBranchAddress("det_chcons", &det_chcons);
    t->SetBranchAddress("part_chcons", &part_chcons);

    t->SetBranchAddress("weight", &weight);


    for(int i = 0; i < (int) t->GetEntries(); i++){
        t->GetEntry(i);
        for(int j = 0; j < (int) part_jetpt->size(); j++){
            double p_pt = part_jetpt->at(j);
            double d_pt = det_jetpt->at(j);

            // only bother to loop over constituents if the jet is in the right pt bin
            double jc_p = 0.0;
            double jc_d = 0.0;
            for(int k = 0; k < (int) part_conspt->at(j).size(); k++){
                jc_p += pow( (part_conspt->at(j).at(k)/p_pt ) , kappa) * part_chcons->at(j).at(k);
            }
            for(int k = 0; k < (int) det_conspt->at(j).size(); k++){
                jc_d += pow( (det_conspt->at(j).at(k)/d_pt ) , kappa) * det_chcons->at(j).at(k);
            }
	    for(int ji = 0; ji < (int) dethist.size(); ji++){
		if(jetptedges[ji] < p_pt && p_pt < jetptedges[ji+1]){
		    parthist[ji]->Fill( /*p_pt,*/ (jc_d - jc_p)/ (double) jc_p, weight );
		}
		if(jetptedges[ji] < d_pt && d_pt < jetptedges[ji+1]){
		    dethist[ji]->Fill( /*d_pt,*/ (jc_d - jc_p)/ (double) jc_p , weight );
		}
	    }
	}
    }
    t->ResetBranchAddresses();
    return;
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

    TH1D* parthist[nJetPtBins];
    TH1D* dethist[nJetPtBins];

    vector<TH1D*> deltaQbyPartQ_bypartpt;
    vector<TH1D*> deltaQbyPartQ_bydetpt;

/*	// old binning-->now using Isaac's jet mass bins: 20-25, 25-30, 30-40
    vector<double> jetPtLo = {5.0, 10.0, 15.0, 20.0, 30.0};
    vector<double> jetPtHi = {10.0, 15.0, 20.0, 30.0, 40.0};

    double jetEdges[nJetPtBins+1] = {5.0, 10.0, 15.0, 20.0, 30.0, 40.0};
*/

    const string k = (string) argv[5];
    double kappa = str_to_double(k);

    TString k_str(k);

    // jet pt histograms want bin edges at integer values, not bin centers
    // [0-1, 1-2, ... , 49-50]
    // so the proper number of jets is integral [jetptLo + 1, jetptHi]
    TH1D* njets_pt = new TH1D("jet_pt_hist_k" + k_str, "", 60, 0.0, 60.0);
    TH1D* njets_pt_pl = new TH1D("jet_pt_hist_part_k" +  k_str, "", 60, 0.0, 60.0);

    TH2D* deltaPtvGePt = new TH2D("deltaPtvGePt", "", 11, 5, 60, 220, -6, 6);
    TH2D* deltaPtvPyPt = new TH2D("deltaPtvPyPt", "", 11, 5, 60, 220, -1, 1);


    // ( detector level jet Q - particle level jet Q ) / particle level jet Q
    // feels more logical to me to put pythia pt on x-axis
    TH2D* deltaQvPyPt = new TH2D("deltaQvPyPt", ";Gen. p^{jet}_{T};#DeltaQ_{#kappa} (Det - Gen) / Gen. #DeltaQ_{#kappa}", 11, 5, 60, 100, -1, 1);
    TH2D* deltaQvGePt = new TH2D("deltaQvGePt", ";Det. p^{jet}_{T};#DeltaQ_{#kappa} (Det - Gen) / Gen. #DeltaQ_{#kappa}", 11, 5, 60, 100, -6, 6);
    // axis titles following Isaac's convention in jetmass2/macros/hists.cxx


    TH1D* part_posQresolution[nJetPtBins];
    TH1D* part_negQresolution[nJetPtBins];
    TH1D* det_posQresolution[nJetPtBins];
    TH1D* det_negQresolution[nJetPtBins];



    TH1D* jetQ[nJetPtBins];// = new TH1D(Form("kappa_k%s", k), "", 9, -4.5, 4.5);
    TH1D* jetQ_part[nJetPtBins];// = new TH1D(Form("kappa_part_k%s", k), "", 9, -4.5, 4.5);
//    TH1D* jetQ05[nJetPtBins];// = new TH1D("kappa05", "", 25, -2.5, 2.5);
//    TH1D* jetQ05_part[nJetPtBins];// = new TH1D("kappa05_part", "", 25, -2.5, 2.5);
	// range: 5.0, bin width ~0.2? 5.0/0.2 = 25

    TH2D* q_v_pt_counts[2];
    TH1D* q_counts[2];
    TH1D* pt_counts[2];

//    int njetQbins = 25; double jetQedge = 2.5;
    int njetQbins = 35; double jetQedge = 3.5;
//    if( kappa == 0.0 ){njetQbins = 9; jetQedge = 4.5;}
    if( kappa == 0.0 ){njetQbins = 13; jetQedge = 6.5;}
    // widen histogram to provide bins in response that are not shown in plot (similar to jet pt range change to include 5-10 GeV jets even though i only show 10-15 GeV and above
    // if kappa = 0.0, all that happens is addition of integers, don't need finer bins than width = 1.0


    // 2D histogram for response, unfolding--- jet pt in 4 bins [10-15, 15-20, 20-30, 30-40]
    TH2D* QvPt = new TH2D("QvPt_d", "", njetQbins, -jetQedge, jetQedge, nJetPtBins, jetptEdges); // data/detector level
    TH2D* QvPt_p = new TH2D("QvPt_p", "", njetQbins, -jetQedge, jetQedge, nJetPtBins, jetptEdges); // particle level (should not be filled for data)

    cout << "BREAK AT a\n";


    int nresbins = 40; int resedge = 5;
    if(kappa == 0.0){nresbins = 11;}

    for(int ij = 0; ij < nJetPtBins; ij++){
	parthist[ij] = new TH1D( Form("deltaQbypartQ_bypartpt_jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", nresbins, -resedge, resedge);
	dethist[ij] = new TH1D( Form("deltaQbypartQ_bydetpt_jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", nresbins, -resedge, resedge);

	deltaQbyPartQ_bypartpt.push_back(parthist[ij]);
	deltaQbyPartQ_bydetpt.push_back(dethist[ij]);

	part_posQresolution[ij] = new TH1D( Form("part_posqResolution_k" + k_str + "jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
	part_negQresolution[ij] = new TH1D( Form("part_negqResolution_k" + k_str + "jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
	det_posQresolution[ij] = new TH1D( Form("det_posqResolution_k" + k_str + "jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
	det_negQresolution[ij] = new TH1D( Form("det_negqResolution_k" + k_str + "jetpt%1.0f%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
        jetQ[ij] = new TH1D( Form( "jetQ_k" + k_str + "_jetpt%1.0f_%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
        if(match){ // if pythia 6 is being used, jets are matched particle to detector level
            jetQ_part[ij] = new TH1D( Form( "jetQ_k" + k_str + "_part_jetpt%1.0f_%1.0f", jetPtLo[ij], jetPtHi[ij] ), "", njetQbins, -jetQedge, jetQedge);
        }
    }
    q_v_pt_counts[0] = new TH2D("q_v_pt_counts", "", njetQbins, -jetQedge, jetQedge, nJetPtBins, jetptEdges);
    q_counts[0] = new TH1D("q_counts", "", njetQbins, -jetQedge, jetQedge);
    pt_counts[0] = new TH1D("pt_counts", "", nJetPtBins, jetptEdges);
    if(!data_bool){
	q_v_pt_counts[1] = new TH2D("PL_q_v_pt_counts", "", njetQbins, -jetQedge, jetQedge, nJetPtBins, jetptEdges);
	q_counts[1] = new TH1D("PL_q_counts", "", njetQbins, -jetQedge, jetQedge);
	pt_counts[1] = new TH1D("PL_pt_counts", "", nJetPtBins, jetptEdges);
    }
    
    vector<TH1D*> k_vec = {jetQ[0], jetQ[1], jetQ[2], jetQ[3]};//, jetQ[4]};
    vector<TH1D*> k_p = {jetQ_part[0], jetQ_part[1], jetQ_part[2], jetQ_part[3]};//, jetQ_part[4]};
    vector<TH2D*> qvpt = {QvPt, QvPt_p};
    vector<TH2D*> qvpt_counts = {q_v_pt_counts[0], q_v_pt_counts[1]};
    vector<TH1D*> qhists = {q_counts[0], q_counts[1]};
    vector<TH1D*> pthists = {pt_counts[0], pt_counts[1]};

    vector<TH1D*> pos_part_jetQresolution = {part_posQresolution[0], part_posQresolution[1], part_posQresolution[2], part_posQresolution[3]};
    vector<TH1D*> neg_part_jetQresolution = {part_negQresolution[0], part_negQresolution[1], part_negQresolution[2], part_negQresolution[3]};
    vector<TH1D*> pos_det_jetQresolution = {det_posQresolution[0], det_posQresolution[1], det_posQresolution[2], det_posQresolution[3]};
    vector<TH1D*> neg_det_jetQresolution = {det_negQresolution[0], det_negQresolution[1], det_negQresolution[2], det_negQresolution[3]};


    cout << "BREAK AT b\n";
    
    //void JetChargeHist(TFile* file, TH1D* hjetcharge, double *njets, double kappa, double jetptLo = 20.0, double jetptHi = 30.0, const bool dataflag)
    //void MatchedJetChargeHist(TFile* file, TH1D* hjetcharge, TH1D* hjetcharge1, double *njets, double kappa, double jetptLo = 20.0, double jetptHi = 30.0, const bool dataflag)
    //void Fill2DHist(TFile* file, TH2D* hist1, TH2D* hist2, double kappa, bool data_bool = 1)
    //void Fill1DCounts(TFile* file, vector<TH1D*> qhists, vector<TH1D*> pthists, double kappa, bool data_bool = 1)
    // call functions to take trees and fill histograms
    if(match){
        // Only call for non-data
        MatchedJetChargeHist(fin, k_vec, k_p, njets_pt, njets_pt_pl, kappa, jetPtLo, jetPtHi);
	Fill2DHist(fin, qvpt, kappa, data_bool);
	Fill2DCounts(fin, qvpt_counts, kappa, data_bool);
	Fill1DCounts(fin, qhists, pthists, kappa, data_bool);

        fillDeltaPtHists(fin, deltaPtvPyPt, deltaPtvGePt);
	fillJetQResolution(fin, pos_part_jetQresolution, neg_part_jetQresolution, pos_det_jetQresolution, neg_det_jetQresolution, kappa, jetPtLo, jetPtHi);

	fillDeltaQHists(fin, deltaQvPyPt, deltaQvGePt, kappa);
	fillDeltaQHists(fin, deltaQbyPartQ_bypartpt, deltaQbyPartQ_bydetpt, jetptEdges, kappa);
    }
    else{
        JetChargeHist(fin, k_vec, njets_pt, kappa, jetPtLo, jetPtHi, data_bool);
	Fill2DHist(fin, qvpt, kappa, data_bool, match);
	Fill2DCounts(fin, qvpt_counts, kappa, data_bool, match);
	Fill1DCounts(fin, qhists, pthists, kappa, data_bool, match);
    }

    cout << "BREAK AT c\n";
    
    // argv[1] is the output file location, argv[2] is the output file name--from Isaac
    TFile *fout = new TFile(((string) argv[1]+(string) argv[2]).c_str(),"RECREATE");
    cout << "DEBUG: output file name is " << fout->GetName() << endl;
    
    for(int i = 0; i < (int) k_vec.size(); i++){
        k_vec[i]->Write();
        if(match){
            k_p[i]->Write();
        }
    }
    cout << "1D hists get written...\n";
    if(match || data_bool){
        QvPt->Write();
    }
    if(!match && !data_bool){
        QvPt_p->Write();
    }

    cout << "why don't 2D hists though...\n";

    if(match){
    //if(!data_bool){
        QvPt_p->Write();

	deltaQvPyPt->Write();
	deltaQvGePt->Write();

	for(int ji = 0; ji < nJetPtBins; ji++){
	    parthist[ji]->Write();
	    dethist[ji]->Write();
	}

	if(!data_bool){
            deltaPtvGePt->Write();
            deltaPtvPyPt->Write();

cout << "if sim, delta pT and q resolution hists written\n";
	    for(int j = 0; j < (int) pos_part_jetQresolution.size(); j++){
		pos_part_jetQresolution[j]->Write();
		neg_part_jetQresolution[j]->Write();
		pos_det_jetQresolution[j]->Write();
		neg_det_jetQresolution[j]->Write();
	    }
	}
    }
	
cout << "oh no...\n";

    njets_pt->Write();
    njets_pt_pl->Write();

    q_v_pt_counts[0]->Write();
    q_counts[0]->Write();
    pt_counts[0]->Write();
    if(match){
        q_v_pt_counts[1]->Write();
        q_counts[1]->Write();
        pt_counts[1]->Write();
    }

    cout << "Wrote to " << fout->GetName() << endl;
    
    fout->Close();
    cout << "Closed " << fout->GetName() << endl;
    
    return 0;
}

