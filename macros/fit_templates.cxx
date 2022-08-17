//
//  fit_templates.cxx
//
//  Created by Grant McNamara on 8/10/22.
//
//  first attempt to define a TF1 as a linear combination of the d, u, g, other templates
//  once working, will then add in getting the unfolded distribution and fit the linear combination TF1 to the unfolded distribution. then extract what should be the free parameters in the sum of gaussians (the coefficients in the linear combination), scale the individual flavor templates by these factors, and plot these later
//
//


#include <stdio.h>
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

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "Headers/funcs.hh"
#include "Headers/JCparameters.hh"


using namespace std;
using namespace jcAnalysis;

int main(int argc, const char** argv){
    if(argc != 3){
        cout << "Should receive jet radius, jet charge kappa. Received " << argc-1 << "parameters. Exiting.\n";
        exit(1);
    }
    
    const int njetbins = 3;
    //double jetEdges[njetbins+1] = {20.0, 25.0, 30.0, 40.0};

    double jetPtLo[njetbins] = {20.0, 25.0, 30.0};
    double jetPtHi[njetbins] = {25.0, 30.0, 40.0};

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();

    string rad = (string) argv[1];
    string kappa = (string) argv[2];
    
    
//    int njetQbins = 18; double jetQedge = 3.5;
//    if( kappa == 0.0 ){njetQbins = 9; jetQedge = 4.5;}
//    if( kappa == "00" ){njetQbins = 13; jetQedge = 6.5;}
    
    
/*
    //const int nflavors = 6;
    const int nflavors = 4;
    TString flavor_names[nflavors] = {"d", "u", "g", "dbar", "ubar", "other"};
                                      //"dbar", "ubar"};
*/
    
    // file with the individual templates
    TFile *ftemp = new TFile( ( "~/ppJCandBF/out/template/R" + rad + "_k" + kappa + ".root" ).c_str() , "READ");
    
    // naming convention for the TH1Ds with associated fitted functions: Form( flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] )
    // naming convention for the associated fitted functions: TString fit_name = Form( "fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] )

    TH1D* qHists[nflavors][njetbins];
    //TF1* fits[nflavors][njetbins];
    double params[nflavors][njetbins][3]; // parameters from the template fits, each template is a gaussian with 3 parameters
    
    TF1* unf_fit[njetbins];
    
    TFile *funf = new TFile( ("~/ppJCandBF/out/unfold/unfolded_R"+ rad + "_k" + kappa + ".root").c_str() , "READ");
    
    TH1D* nom[njetbins];
    
    double params_unf[njetbins][4];
    
    /*
    for(int j = 0; j < njetbins; j++){
        nom[j] = (TH1D*) funf->Get( Form( "unfold_nomx", jetPtLo[j], jetPtHi[j] ) );
        TString fit_name = "";
        for(int i = 0; i < nflavors; i++){
            qHists[i][j] = (TH1D*) ftemp->Get( Form( flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
            
            cout << "histogram breaking? integral = " << qHists[i][j]->Integral() << "\n";
            
            
            cout << "trying to load TF1s from the histograms\n";
            // the below line crashes, segmentation faults
            fits[i][j] = (TF1*) qHists[i][j]->GetFunction( Form( "fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
            cout << "TF1s seemingly loaded successfully...\n";
            
            fit_name += Form( "[%i] * fit_" + flavor_names[i] + "_%1.0f%1.0f(x)", i, jetPtLo[j], jetPtHi[j] );
            //fit_name += Form( "[%i] * fits[%i][%i](x)", i, i, j );
            if( i + 1 != nflavors ){fit_name += " + ";}
        }
        //cout << "DEBUG: formula for sum of fits for jet pt bin " << j << " = " << fit_name << "\n";
        //fit_sum[j] = new TF1( Form( "fit_sum_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ), fit_name);
        
        
    */
    // gaussian formula:
    // [Constant]*exp(-0.5*((x-[Mean])/[Sigma])*((x-[Mean])/[Sigma]))
    // [0]                      [1]      [2]
    
    cout << "DEBUG: getting the unfolded histograms\n";
    for(int j = 0; j < njetbins; j++){
        nom[j] = (TH1D*) funf->Get( "unfold_nomx" + outFile_pt[j] );
        TString fit_name = "";
        cout << "getting the template histogram and fits\n";
        for(int i = 0; i < nflavors; i++){
            qHists[i][j] = (TH1D*) ftemp->Get( flavor_names[i] + "_" + outFile_pt[j] );
            
            TString histfit_name = "fit_" + flavor_names[i] + "_" + outFile_pt[j];
            
            
            /*
            qHists[i][j] = (TH1D*) ftemp->Get( Form( flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
            
            TString histfit_name = Form( "fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] );
            */
            cout << "geting parameters...\n";
            for(int k = 0 ; k < 3; k++){
                params[i][j][k] = qHists[i][j]->GetFunction( histfit_name )->GetParameter( k );
            }
            cout << "parameters successfully(?) gotten\n";
//            fit_name += "[" + i + "] * (" + params[i][j][0] + " * exp( -0.5 * ( (x - " + params[i][j][1] + ") /" + params[i][j][2] + ") * ( (x - " + params[i][j][1] + ") /" + params[i][j][2] + ") ) )", i );
            fit_name += Form( "[%i] * ( %f * exp( -0.5 * ( (x - %f) / %f) * ( (x -  %f) /%f) ) )", i, (float) params[i][j][0], (float) params[i][j][1], (float) params[i][j][2], (float) params[i][j][1], (float) params[i][j][2] );
            if( i + 1 != nflavors ){fit_name += " + ";}
        }
        cout << "flavors done\n";
        
        cout << "DEBUG: new formula for sum of gaussians = " << fit_name << "\n";
        
        unf_fit[j] = new TF1( "fit_unfold_" + outFile_pt[j], fit_name );
        
        for(int k = 0; k < nflavors; k++){
            unf_fit[j]->SetParLimits( k, 0., 1. );
        }

        
        nom[j]->Fit( unf_fit[j] );
        //nom[j]->Fit( Form( "fit_other_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ), fit_name);
        
        //fit_sum[j] = new TF1( Form( "fit_sum_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ), fit_name);
        
        // unfolded distribution->Fit(fit_sum[j])
        //nom[j]->Fit(fit_sum[j], "", "", -jetQedge, jetQedge);
        
        cout << "\n";
        // print out fit_sum[j]->GetParameter(0)
        for(int k = 0 ; k < unf_fit[j]->GetNpar(); k++){
            cout << "parameter " << k << ": " << nom[j]->GetFunction( "fit_unfold_" + outFile_pt[j] )->GetParameter(k) << "\n";
            TString unf_par_name = "fit_unfold_" + outFile_pt[j];
            params_unf[j][k] = nom[j]->GetFunction( unf_par_name )->GetParameter( k );
        }
    }
    
    TFile *fout = new TFile( ( "~/ppJCandBF/out/fit_templates/R" + rad + "_k" + kappa + ".root" ).c_str(), "RECREATE");
    fout->cd();
    
    cout << "output file name: " << ( "~/ppJCandBF/out/fit_templates/R" + rad + "_k" + kappa + ".root" ).c_str() << "\n";
    
    for(int j = 0; j < njetbins; j++){
        nom[j]->Write();
        
    }
    
    return 0;
}


