//
//  draw_unfoldedtemplate.C
//  
//
//  Created by Grant McNamara on 8/14/22.
//

#include <string>
#include <iostream>
#include "math.h"
#include "Headers/plot.h"
//#include "Headers/JCparameters.hh"


// plot pythia-6, pythia-6+geant, data and data unfolded for detector effects
void draw_unfoldedtemplate(double k = 0.0, double R = 0.4){
    // ../out/data/hists/ppDataJP2_R0    _R04.root
    //            ppDataJP2_R06.root

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";

    
    const int nflavors = 4;
    TString flavor_names[nflavors] = { "d",     "u",    "g",
                                       /*"dbar", "ubar",*/ "other" };
    
    
    int mark_colors[nflavors] = { 600, 628, 839,
                                  /*600, 628,*/ 617 };


    const int njetbins = 3;

    double jetPtLo[njetbins] = {20.0, 25.0, 30.0};
    double jetPtHi[njetbins] = {25.0, 30.0, 40.0};
    // loop over all jet pt ranges


    double jetQedge = 3.5;
    if( k == 0.0 ){jetQedge = 6.5;}
    
    
    // file containing templates by flavor
    TFile* ftemp = new TFile( "~/ppJCandBF/out/template/R" + rad + "_k" + kappa + ".root" , "READ");
    
    TH1D* temps[njetbins][nflavors];
    TF1* temp_fits[njetbins][nflavors];
    
    
    TFile *funf = new TFile( "~/ppJCandBF/out/fit_templates/R" + rad + "_k" + kappa + ".root" , "READ");
    
    TH1D* unfs[njetbins];
    TF1* unf_fits[njetbins];
    
    
    double temp_params[njetbins][nflavors][3];
    double unf_params[njetbins][nflavors];
    
    TF1* fit_plot[njetbins][nflavors];
    
    TCanvas *c[njetbins];
    TLegend *leg[njetbins];


    for(int j = 0; j < njetbins; j++){
        c[j] = new TCanvas(Form("c_%i", j), "canvas", 600, 600);
        leg[j] = new TLegend(0.18, 0.55, .4, .85);
        leg[j]->SetBorderSize(0);
        //leg[j]->SetNColumns(3);

        c[j]->SetLeftMargin(0.15);
        
        unfs[j] = (TH1D*) funf->Get( Form( "unfold_nomx%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );

        unf_fits[j] = (TF1*) unfs[j]->GetFunction( Form( "fit_unfold_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
        
        //unfs[j]->SetMarkerSize(1);
        unfs[j]->SetMarkerColor(kBlack);
        unfs[j]->SetMarkerStyle(kFullCircle);
        
        
        unfs[j]->SetStats(0);
        unfs[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV", jetPtLo[j], jetPtHi[j] ) );
        unfs[j]->SetXTitle( Form( "Q^{jet}_{#kappa = %1.1f}" , k) );
        unfs[j]->SetYTitle("Prob.");
        
        unfs[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
        unfs[j]->Draw("hist p");
        
        unf_fits[j]->SetLineColor(kBlack);
        unf_fits[j]->SetLineStyle( 9 );
        unf_fits[j]->SetLineWidth( 2 );
        
        
        unf_fits[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
        //unf_fits[j]->Draw("C same");
        
        leg[j]->AddEntry( unfs[j], "unfolded data", "p" );
        leg[j]->AddEntry( unf_fits[j], "unfolded fit", "l" );
        
        for(int i = 0; i < nflavors; i++){
            temps[j][i] = (TH1D*) ftemp->Get( Form( flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
            
            temp_fits[j][i] = (TF1*) temps[j][i]->GetFunction( Form( "fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
            
            for(int k = 0; k < 3; k++){
                temp_params[j][i][k] = temp_fits[j][i]->GetParameter( k );
            }
            //cout << "DEBUG: where do we break?\n";
            
            unf_params[j][i] = unf_fits[j]->GetParameter( i );

            temp_fits[j][i]->SetLineColor( mark_colors[i] );
            
        }
        
        //cout << "DEBUG: EVERYTHING RETRIEVED FROM FILES CORRECTLY\n";
        
        TString fit_formula = "";
        
        for(int i = nflavors-1; i >= 0; i--){
            fit_formula += Form( "(%f) * ( %f * exp( -0.5 * ( (x - %f) / %f) * ( (x -  %f) /%f) ) )", (float) unf_params[j][i], (float) temp_params[j][i][0], (float) temp_params[j][i][1], (float) temp_params[j][i][2], (float) temp_params[j][i][1], (float) temp_params[j][i][2] );
            
            //cout << "DEBUG FIT " + flavor_names[i] +  " FORMULA: " << fit_formula << "\n";
            
            
            fit_plot[j][i] = new TF1( Form( "plot_fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) , fit_formula, -jetQedge, jetQedge );
            
            fit_formula += " + ";
            
            
            //fit_plot[j][i]->SetFillColorAlpha( mark_colors[i], 0.35 );
            fit_plot[j][i]->SetLineColor( mark_colors[i] );
            //fit_plot[j][i]->SetFillStyle( 1001 );

            if( flavor_names[i] == "other" ){
                fit_plot[j][i]->SetLineStyle( 9 );
            }
            else{
                fit_plot[j][i]->SetLineStyle( 1 );
            }

            //temp_fits[j][i]->SetStats(0);
            
        }
        for(int i = 0; i < nflavors; i++){
            leg[j]->AddEntry( fit_plot[j][i], flavor_names[i], "l" );
            
            fit_plot[j][i]->Draw(/*"F*/"C SAME");

        } // flavor loop
        //cout << "\n";
        
        unf_fits[j]->Draw("C same");
        
        leg[j]->Draw();
        // save c[j] here
        drawText( Form( "R = %0.1f", R ), .65, .7, 30 );

        // Save plot here, after points have been plotted
        c[j]->SaveAs( Form( "./plots/unfoldtemplate_jetQ%1.0f%1.0f_R" + rad + "_k" + kappa + ".png", jetPtLo[j], jetPtHi[j] ), "RECREATE" );
        
    } // jet pt loop
}



