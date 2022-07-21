//
//  stat_err_scaling.cxx
//  
//
//  Created by Grant on 2/19/22.
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

#include "RooUnfold.h"

using namespace std;



double str_to_double (std::string str) {
    std::string Num = str.substr(0,1)+"."+str.substr(1,1); //e.g. 0.4
    double str_double = (double) std::stod(Num); //converting the string 0.4 to a double
    std::cout << "DEBUG: string to number = " << str_double << std::endl;

    return str_double;
}



int main (int argc, const char** argv) {
    //intro
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //basic argument checking.
    if (argc != 3) {
        cerr << "Should be two argument: jet radius, kappa. Received "
            << argc-1 << ". Exiting." << endl;
	cout << "Received " << argv[0] << ", " << argv[1] << ", " << argv[2] << "\n";
        exit(1);
    }
    
    //argv[1] should be the jet radius e.g. "04".
    string radius = (string) argv[1];
    radius = "_R" + radius;//appending the "_R" used in file naming.
    
    // should pass kappa as argument, which index is TBD
    string kappa = (string) argv[2];
    
//    double k = str_to_double(kappa);
//    str_to_double(kappa);
    
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    const string path = "~/ppJCandBF/out/sim/response/";
    const string file_in = "p6";
    
    TFile *fin = new TFile((path + file_in + radius + "_k" + kappa + "_syst_July20.root").c_str(),"READ");
    // file with responses will be in out/sim/response/p6_R + radius + "_k" + kappa + ".root" or about there
    cout << "DEBUG: input file name is " << fin->GetName() << endl;
    
    RooUnfoldResponse* res = (RooUnfoldResponse*) fin->Get("pt_response"); // pt_response comes from hadded output files from running response.cxx
    TH1D* match_plus_miss = (TH1D*) fin->Get("jetpt_match_plus_miss");
    TH1D* match = (TH1D*) res->Hresponse()->ProjectionY("match");
    
    TH1D* hratio = (TH1D*) match_plus_miss->Clone("hratio");
    TH1D* efficiency = (TH1D*) match->Clone("efficiency");
    
//cout << "ratio binning: " << hratio->GetNbinsX() << " bins from " << hratio->GetBinLowEdge(0) << " to " << hratio->GetBinLowEdge( hratio->GetNbinsX() ) << "\n";
//cout << "efficiency binning: " << efficiency->GetNbinsX() << " bins from " << efficiency->GetBinLowEdge(0) << " to " << efficiency->GetBinLowEdge( efficiency->GetNbinsX() ) << "\n";

    
    hratio->Divide(match);
    efficiency->Divide(match_plus_miss);
    
    const int nBins = hratio->GetNbinsX();
    
    for (int i = 1; i <= nBins; ++ i) {
        hratio->SetBinContent(i, sqrt(hratio->GetBinContent(i)));
    }
    
    TFile *fout = new TFile((path + "stat_err_scaling" + radius + "_k" + kappa + ".root").c_str(),"RECREATE");
    
    fout->cd();
    hratio->Write();
    efficiency->Write();
    
    cout << endl
         << "Wrote statistical error scaling to " << fout->GetName() << endl;
    
    //closing files
    fout->Close();
    cout << endl << "Closed " << fout->GetName() << endl;
    
    return 0;
}
