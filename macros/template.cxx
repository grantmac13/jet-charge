//
//  template.cxx
//  
//
//  Created by Grant McNamara on 8/8/22.
//  This file projects flavor + "_qvpt" histograms onto the jet Q axis, parameterize the distribution (starting with as gaussians, function subject to change) and fit a linear combination of these templates to the unfolded distributions for a given kappa
//
//  treat the same as macros/hists.cxx, macros/reponse.cxx, etc  with input from command line


#include <stdio.h>
#include <string>
#include "math.h"

#include <ctime>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <algorithm>

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


#include "Headers/plot.h"
#include "Headers/JCparameters.hh"


using namespace std;
using namespace jcAnalysis;


int main(int argc, const char** argv){
    
    if(argc != 3){
        cout << "Should receive jet radius, jet charge kappa. Received " << argc-1 << "parameters. Exiting.\n";
        exit(1);
    }
    
    string rad = (string) argv[1];
    string kappa = (string) argv[2];
    
    
//    int njetQbins = 18; double jetQedge = 2.0;
//    if( kappa == 0.0 ){njetQbins = 9; jetQedge = 4.5;}
//    if( kappa == "00" ){njetQbins = 13; jetQedge = 4.5;}
    
    // until pythia 6 root files remade by isaac with parton information included, p8 filepath hardcoded
    TFile* fin = new TFile( ( "~/for_grant/out/p8_decayed_R" + rad + "_k" + kappa + ".root" ).c_str() , "READ");
    
    
    /*
    // 6 quarks * 2 + gluon = 13
    const int nflavors = 13;
    // flavor_names following pdg numbering scheme order for quarks
    //                                 0       1       2
    TString flavor_names[nflavors] = {"d",    "u",    "s",
    //                                 3       4       5       6
                                      "c",    "b",    "t",    "g",
    //                                 0       1       2
                                      "dbar", "ubar", "sbar",
    //                                 3       4       5
                                      "cbar", "bbar", "tbar"};
    */
    
    
    //const int njetbins = 3;

    //double jetPtLo[njetbins] = {/*15.0,*/ 20.0, 25.0, 30.0};
    //double jetPtHi[njetbins] = {/*20.0,*/ 25.0, 30.0, 40.0};
    
    TF1* flavor_fits[n_template_flavors][njetbins]; // will store the fit functions after jetQ_proj[i][j]->Fit() is called
    
    TH2D* jetQ[n_template_flavors];
    TH1D* jetQ_proj[n_template_flavors][njetbins];
    
    TF1* other_fit[njetbins];
    TH1D* other_jetQ[njetbins];
    
    
    for(int i = 0; i < n_template_flavors; i++){
        // get the 2D histogram for the i-th flavor
        jetQ[i] = (TH2D*) fin->Get( template_flavor_names[i] + "_qvpt" );
        
        for(int j = 0; j < njetbins; j++){
            cout << "jet pt bin loop started\n";
            TString fit_name = "fit_" + template_flavor_names[i] + "_" + outFile_pt[j];
            
            //TString fit_name = Form( "fit_" + template_flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] );
            
            // define function to use to fit to template distributions
            //flavor_fits[i][j] = new TF1( fit_name, "[0] * exp( pow( (x-[1])/[2], 2 ) )", -jetQedge, jetQedge);
            flavor_fits[i][j] = new TF1( fit_name, "gaus" );
            
            cout << "tf1 made\n";
            
            // project the 2D by jetpt bin (TH1D::ProjectionX() has weird (non-intuitive) bin convention so first and last bin should always be the same
            // j+2 bin because 0 is underflow, 1 is 15-20 GeV jet pt, only want 20-25 and above
/*
            jetQ_proj[i][j] = (TH1D*) jetQ[i]->ProjectionX( Form( template_flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ), j+2, j+2 );

            TString other_fit_name = Form( "fit_other_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] );
*/
            
            jetQ_proj[i][j] = (TH1D*) jetQ[i]->ProjectionX( template_flavor_names[i] + "_" + outFile_pt[j], j+2, j+2 );

            TString other_fit_name = "fit_other_" + outFile_pt[j];
            
            if(i == 2){ // 0=d, 1=u, 2 != g, so 2 is first to put into "other" template with the rest of the flavors, c, s, t, b combine into "other" for simplicity
                //other_jetQ[j] = (TH1D*) jetQ_proj[i][j]->Clone( Form( "other_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] )  );
                
                other_jetQ[j] = (TH1D*) jetQ_proj[i][j]->Clone( "other_" + outFile_pt[j] );
                
                other_fit[j] = new TF1( other_fit_name, "gaus" );
                cout << "fit name for other: " << other_fit_name << "\n";
                
            }
            else if( ! ( i%7 == 0 || i%7 == 1 || i%7 == 6 ) ){
                other_jetQ[j]->Add( jetQ_proj[i][j] );
            }
            
            if( i + 1 == n_template_flavors ){
                
                other_jetQ[j]->Scale( 1.0 / ( other_jetQ[j]->Integral( 0, other_jetQ[j]->GetNbinsX() + 1 ) ) );
                other_jetQ[j]->Fit( other_fit_name );
            }
            
            jetQ_proj[i][j]->Scale( 1.0 / ( jetQ_proj[i][j]->Integral( 0, jetQ_proj[i][j]->GetNbinsX() + 1 ) ) );

            cout << "about to perform the fit\n";
            jetQ_proj[i][j]->Fit( fit_name );
            // how do i get the parameters?
            // maybe using jetQ_proj[i][j]->GetFunction("func_name");
            
            // store fit function parameterization
            flavor_fits[i][j] = jetQ_proj[i][j]->GetFunction( fit_name );
            
            double* params = flavor_fits[i][j]->GetParameters();
            // gaus defined as f(x) = [0]*exp(-1/2 ( (x-[1]) / [2] )^2 )
            
            // print the parameters of the fit
            cout << "\nfor flavor " << i << " and jet pt " << jetPtLo[j] << "-" << jetPtHi[j] << "GeV\n";
            for(int k = 0; k < flavor_fits[i][j]->GetNpar(); k++){
                cout << "parameter " << k << " = " << params[k] << "\n";
            }
        }
    }
    
    
    TFile *fout = new TFile( ( "~/ppJCandBF/out/template/R" + rad + "_k" + kappa + ".root" ).c_str() , "RECREATE");
    fout->cd();
    
    cout << "output file name: " << ( "~/ppJCandBF/out/template/R" + rad + "_k" + kappa + ".root" ).c_str() << "\n";
    
    
    for(int i = 0; i < n_template_flavors; i++){
        for(int j = 0; j < njetbins; j++){
            //flavor_fits[i][j]->Write(); // template functions (TF1*)
            // no need to write TF1 to output file, it is stored within the TH1D it appears
            jetQ_proj[i][j]->Write(); // projections of qvpt onto jet q by flavor of parton (TH1D*)
        }
    }
    for(int j = 0; j < njetbins; j++){
        other_jetQ[j]->Write();
    }
    
    return 0;
}

