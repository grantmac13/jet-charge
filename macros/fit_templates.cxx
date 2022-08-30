//
//  fit_templates.cxx
//
//  Created by Grant McNamara on 8/10/22.
//
//  first attempt to define a TF1 as a linear combination of the d, u, g, other templates
//  once working, will then add in getting the unfolded distribution and fit the linear combination TF1 to the unfolded distribution. then extract what should be the free parameters in the sum of gaussians (the coefficients in the linear combination), scale the individual flavor templates by these factors, and plot these later
//
//  8/17/22-- consider fitting all 13 flavors to unfolded result, and THEN constructing "other" as the weighted sum of flavors not u, d, g, ubar, dbar
//      changes: from nflavors, flavor_names -> n_template_flavors, template_flavor_names (6 flavors, u, d, g, ubar, dbar, other to all 13
//      above did not work for some reason, reverting back to nflavors, etc to test the formula expression
//
//  8/20(21)/22 also maybe try to use only u, d, g, ubar, dbar and then a generic gaussian like i had thought about doing to hopefully be able to fit well to unfolded data with relatively strict parameter limits
//
//  8/29/22: two new things to note: now that I am trying to use TFractionFitter, the pythia distributions to use with TFractionFitter should be NOT normalized i believe, but if they are normalized prior to the fitting, the shape does in fact match; second point is that now I will be temporarily looking at pythia-8 as proof of concept that the fitting *can* produce the correct result (as pythia-8 is known to follow the known flavor fractions)
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

#include <TFractionFitter.h>
#include <TObjArray.h>

#include "TStyle.h"


using namespace std;
using namespace jcAnalysis;

int main(int argc, const char** argv){
    if(argc != 3){
        cout << "Should receive jet radius, jet charge kappa. Received " << argc-1 << "parameters. Exiting.\n";
        exit(1);
    }
    
    gStyle->SetOptFit(1111);
    
    string rad = (string) argv[1];
    string kappa = (string) argv[2];
    // add third (and fourth?) argument for fixed/constrained(and to what degree?)/free?
    // instead of having to change it below
    
    
    int n = n_template_flavors/*nflavors*/; // change this variable n between n_template_flavors and nflavors and use n locally in this macro
    const TString *flavs = template_flavor_names;
    
    double tff_params[njetbins][n];
    double tff_errs[njetbins][n];
    
    TH1D* tff_hist[njetbins]; // histograms to fill with weighted sum using parameters from TFracitonFitter fitting to unfolded distribution
    
    double constraint = 0.25; // want to test fitting unfolded distribution to weighted sum of flavor templates, but if the weights are allowed to freely vary during the fitting process the results are unphysical/unrealistic so this is to test if there is a decent fit to the unfolded data where the weights are near (and to see how near this decent fit is to) the pythia8 jet fractions per flavor
    // so we will look at %s of initial value 0.1 (already done), 0.15, 0.2, 0.25
    
    
    // file with the individual templates
    TFile *ftemp = new TFile( ( "~/ppJCandBF/out/template/R" + rad + "_k" + kappa + ".root" ).c_str() , "READ");
    
    // naming convention for the TH1Ds with associated fitted functions: Form( flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] )
    // naming convention for the associated fitted functions: TString fit_name = Form( "fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] )

    TH1D* qHists[n][njetbins];
    //TF1* fits[nflavors][njetbins];
    double params[n][njetbins][3]; // parameters from the template fits, each template is a gaussian with 3 parameters
    
    TF1* unf_fit[njetbins];
    
    //TFile *funf = new TFile( ("~/ppJCandBF/out/unfold/unfolded_R"+ rad + "_k" + kappa + ".root").c_str() , "READ");
    TFile *funf = new TFile( ( "~/ppJCandBF/out/template/R" + rad + "_k" + kappa + ".root" ).c_str() , "READ"); // for fitting to pythia-8
    
    TH1D* nom[njetbins];
    
    double params_unf[njetbins][n];
    
    
    TH1D* rel_abunds[n][njetbins]; // integral is number of jets per flavor per jet pt range in simulation
    TH1D* flav_sum_hist[njetbins]; // integral is total number of jets in simulation
    
    
    TF1* extra_gaus[njetbins];
    
    
    TObjArray *h_fit[njetbins]; // array of histograms used as templates for TFractionFitter
    TH1D* result_tff[njetbins]; // result histograms from TFractionFitter
    //= new TObjArra(nflavors);
    //TH1D* result_mcpredict[njetbins][n];
    
    /*
    for(int j = 0; j < njetbins; j++){
        nom[j] = (TH1D*) funf->Get( Form( "unfold_nomx", jetPtLo[j], jetPtHi[j] ) );
        TString fit_name = "";
        for(int i = 0; i < n_template_flavors; i++){
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
    // reminder of the order of the parameters in the prededined gaussian function
    // gaussian formula:
    // [Constant]*exp(-0.5*((x-[Mean])/[Sigma])*((x-[Mean])/[Sigma]))
    // [0]                      [1]      [2]
    
    cout << "DEBUG: getting the unfolded histograms\n";
    for(int j = 0; j < njetbins; j++){
        h_fit[j] = new TObjArray(n);
        
        //nom[j] = (TH1D*) funf->Get( "unfold_nomx" + outFile_pt[j] );
        nom[j] = (TH1D*) funf->Get( "flav_sum_hist" + outFile_pt[j] ); // for fitting to pythia-8
        
        nom[j]->Scale( 1.0 / nom[j]->Integral() );
        
        TString fit_name = "";
        cout << "getting the template histogram and fits\n";
        
        flav_sum_hist[j] = (TH1D*) ftemp->Get( "flav_sum_hist" + outFile_pt[j] );
        
        extra_gaus[j] = new TF1( "fit_other_" + outFile_pt[j] , "gaus");
        
        for(int i = 0; i < n; i++){
            qHists[i][j] = (TH1D*) ftemp->Get( flavs[i] + "_" + outFile_pt[j] );
            //cout << "BEFORE integral of flav_sum_hist: " << flav_sum_hist[j]->Integral() << "\n";
            
            //qHists[i][j]->Scale( 1.0 / qHists[i][j]->Integral() );

            //cout << "AFTER integral of flav_sum_hist: " << flav_sum_hist[j]->Integral() << "\n";
            
            
            if(i == 0){tff_hist[j] = (TH1D*) qHists[i][j]->Clone( "tff_result_hist" + outFile_pt[j] );} // clone first flavor hist to then scale by param[0], then add qHist[i][j], param[i] for i != 0 to this
            
            h_fit[j]->Add( qHists[i][j] );
            
            rel_abunds[i][j] = (TH1D*) ftemp->Get( "rel_abunds_" + flavs[i] + "_" + outFile_pt[j] );
            
            
            TString histfit_name = "fit_" + flavs[i] + "_" + outFile_pt[j];
            
            /*
            qHists[i][j] = (TH1D*) ftemp->Get( Form( flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
            
            TString histfit_name = Form( "fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] );
            */

            cout << "geting parameters...\n";
            for(int k = 0 ; k < 3; k++){
                params[i][j][k] = qHists[i][j]->GetFunction( histfit_name )->GetParameter( k );
            }
            cout << "parameters successfully(?) gotten\n";

            //if( flavor_names[i] == "t" || flavor_names[i] == "tbar" ){continue;}// t, tbar are not present, so temporary work around to not include in the sum of fits
            //fit_name += Form( "[%i] * ( %f * exp( -0.5 * ( (x - %f) / %f) * ( (x -  %f) /%f) ) )", i, (float) params[i][j][0], (float) params[i][j][1], (float) params[i][j][2], (float) params[i][j][1], (float) params[i][j][2] );

            // testing storing the formula in JCparameters.hh
            fit_name += Form( formula_base, i, (float) params[i][j][0], (float) params[i][j][1], (float) params[i][j][2], (float) params[i][j][1], (float) params[i][j][2] );
            
            if( i + 1 != n ){fit_name += " + ";}
        }
        cout << "flavors done\n";
        
        TFractionFitter* fit = new TFractionFitter( nom[j], h_fit[j] );
        
        // first fit to just u, d, g and maybe ubar, dbar and see how good/bad the fit is
        // then later add an extra gaussian term to act as the contribution from the other 8 flavors (realistically 6, since no t, tbar in pythia8 even)
        //fit_name += " + ";
        //fit_name += Form( gaus_base, nflavors,    nflavors + 1, nflavors + 2, nflavors + 1, nflavors + 2 );
        
        cout << "DEBUG: new formula for sum of gaussians = " << fit_name << "\n";
        
        unf_fit[j] = new TF1( "fit_unfold_" + outFile_pt[j], fit_name );
        
        for(int k = 0; k < n; k++){
            
            // set initial values for parameters to the relative abundances of each flavor of jet
            unf_fit[j]->SetParameter( k, rel_abunds[k][j]->Integral() / flav_sum_hist[j]->Integral() );
            //unf_fit[j]->FixParameter( k, rel_abunds[k][j]->Integral() / flav_sum_hist[j]->Integral() );
            
            // different levels of restriction: parameters definitely need to be between 0 and 1, then try varying levels of strictness ~ 10%, 15%, 20%, 25% to see how far the restriction needs to be lifted to find a decent fit to unfolded data
            //unf_fit[j]->SetParLimits( k, 0., 1. );
            unf_fit[j]->SetParLimits( k,
                                     (1.0 - constraint) * ( rel_abunds[k][j]->Integral() / flav_sum_hist[j]->Integral() ) ,
                                     (1.0 + constraint) * ( rel_abunds[k][j]->Integral() / flav_sum_hist[j]->Integral() ) );
            
            fit->Constrain( k, //0.0, 1.0 );
                          /*(1.0 - constraint) * */ ( rel_abunds[k][j]->Integral() / flav_sum_hist[j]->Integral() ) ,
                          /*(1.0 + constraint) * */ ( rel_abunds[k][j]->Integral() / flav_sum_hist[j]->Integral() )  );
            
            cout << "DEBUG: relative percentage of flavor " + flavs[k] + " in simulation = " << rel_abunds[k][j]->Integral() / flav_sum_hist[j]->Integral() << "\n";
        }
        
        //fit->SetRangeX(4,10);
        int status = fit->Fit();
        cout << "TFractionFitter fit status: " << status << "\n";
        
        result_tff[j] = (TH1D*) fit->GetPlot();
        result_tff[j]->SetName( "frac_fit_" + outFile_pt[j] );
        
        
        // when adding a generic gaussian to u,d,g the first parameter will be > 0, probably < ~.3 (individual templates have constants on the order of .3 to be normalized to integral unity), no constraint on mean or width?
        /*
        if(unf_fit[j]->GetNpar() > nflavors){
            
        }
        */
        
        // no option (chi-square fit) is how i showed things on/before 8/18
        // "WL" to try weighted likelihood to show later
        // ""B" when fixing parameters to show the relative fraction of jets by flavor
        nom[j]->Fit( unf_fit[j], "B" /*, "WL" */);
        
        // either fit to function with dofs or associate histogram with function fixed by relative abundances
        //nom[j]->GetListOfFunctions()->Add( unf_fit[j] );
        
        //nom[j]->Fit( Form( "fit_other_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ), fit_name);
        
        //fit_sum[j] = new TF1( Form( "fit_sum_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ), fit_name);
        
        // unfolded distribution->Fit(fit_sum[j])
        //nom[j]->Fit(fit_sum[j], "", "", -jetQedge, jetQedge);
        
        cout << "\n";
        // print out fit_sum[j]->GetParameter(0)
        double tot_prob = 0;

        TString unf_par_name = "fit_unfold_" + outFile_pt[j];
        
        cout << "X^2 for fit to jet pt bin " << j << " = " << nom[j]->GetFunction( unf_par_name )->GetChisquare() << "\n";
        
        for(int k = 0 ; k < unf_fit[j]->GetNpar(); k++){
            cout << "parameter " << k << ": " << nom[j]->GetFunction( "fit_unfold_" + outFile_pt[j] )->GetParameter(k) << "\n";
            
            params_unf[j][k] = nom[j]->GetFunction( unf_par_name )->GetParameter( k );
            
            //result_mcpredict[j][k] = (TH1D*) fit->GetMCPrediction(k);
            //result_mcpredict[j][k]->SetName( flavs[k] + outFile_pt[j] );
            
            // find the sum of the parameters from TFractionFitter
            fit->GetResult( k, tff_params[j][k], tff_errs[j][k] );
            
            tot_prob += tff_params[j][k];
            
            cout << "fraction " << k << " = " << tff_params[j][k] << "\n";
            
            if(k == 0){ // first flavor is cloned to tff_hist, just need to scale by the fit parameter for this flavor
                tff_hist[j]->Scale( tff_params[j][k] );
            }
            else{ // all other flavors can be added as scaled by their respective parameters using TH1::Add(hist, weight)
                tff_hist[j]->Add( qHists[k][j], tff_params[j][k] );
            }
        }
        cout << "sum of TFractionFitter parameters = " << tot_prob << "\n";
    }
    
    TFile *fout = new TFile( ( "~/ppJCandBF/out/fit_templates/R" + rad + "_k" + kappa + ".root" ).c_str(), "RECREATE");
    fout->cd();
    
    cout << "output file name: " << ( "~/ppJCandBF/out/fit_templates/R" + rad + "_k" + kappa + ".root" ).c_str() << "\n";
    
    for(int j = 0; j < njetbins; j++){
        nom[j]->Write();
        result_tff[j]->Write();
        tff_hist[j]->Write();
        
        /*
        for(int k = 0; k < n; k++){
            result_mcpredict[j][k]->Write();
        }
        */
    }
    
    return 0;
}


