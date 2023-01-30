//
//  draw_pythia_and_data.C
//
//  Created by Grant McNamara on 9/19/22.
//
//  plot (total) jet charge distributions by jet pt range for pythia-6, pythia-8, and data to compare to see if pythia-8 describes data as well as pythia-6 as well as other comparisons to data/pythia-6
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
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TRandom.h>

#include "Headers/plot.h"
#include "Headers/JCparameters.hh"




using namespace std;
using namespace jcAnalysis;



void draw_pythia_and_data(double k = 0.0, double R = 0.4){
    
    
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    
    // file with pythia-8 distribution
    TFile *fp8 = new TFile( "~/Desktop/p8_flavorTemplates_R" + rad + "_k" + kappa + ".root" , "READ");

    // file with pythia-6 distribution
    TFile *fp6 = new TFile( "~/Desktop/p6_R" + rad + "_k" + kappa + ".root" , "READ");

    // file with unfolded data distribution
    TFile *fdat = new TFile( "~/Desktop/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    
    

    
    TH1D* unf[njetbins]; // unfolded data jet charge
    TH1D* py8[njetbins]; // pythia-8 jet charge (normalized)
    TH1D* py6[njetbins]; // pythia-6 jet charge (normalized)

        
    
    TCanvas *c[njetbins];
    TLegend *leg[njetbins];
    
    
    //cout << "DEBUG: getting the unfolded histograms\n";
    for(int j = 0; j < njetbins; j++){
        
        // unnormalized pythia-8 truth distribution to fit to
        py8[j] = (TH1D*) fp8->Get( "pythia8_norm_" + outFile_pt[j] );
        
        py6[j] = (TH1D*) fp6->Get( Form( "jetQ_k" + kappa + "_part_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j] ) );
        
        unf[j] = (TH1D*) fdat->Get( "unfold_nomx" + outFile_pt[j] );
        
        // start plotting stuff
        c[j] = new TCanvas(Form("c_%i", j), "canvas", 600, 600);
        leg[j] = new TLegend(0.18, 0.55, .4, .85);
        leg[j]->SetBorderSize(0);
        //leg[j]->SetNColumns(3);
        
        c[j]->SetLeftMargin(0.15);
        
        
        
        unf[j]->SetMarkerStyle(kOpenSquare);
        unf[j]->SetMarkerColor(kBlack);
        unf[j]->SetLineColor(kBlack);

        
        py8[j]->SetMarkerStyle(kOpenCircle);
        py8[j]->SetMarkerColor(kRed);
        py8[j]->SetLineColor(kRed);

        
        py6[j]->SetMarkerStyle(kOpenCircle);
        py6[j]->SetMarkerColor(kBlue);
        py6[j]->SetLineColor(kBlue);
        
                
        
        leg[j]->AddEntry( unf[j], "unfolded data", "p");
        leg[j]->AddEntry( py6[j], "pythia-6", "p");
        leg[j]->AddEntry( py8[j], "pythia-8", "p");
        
        
        unf[j]->SetStats(0);
        
        
        // just in case normalization is wrong (should only be necessary for pythia-6
        unf[j]->Scale( 1.0 / unf[j]->Integral( 0, unf[j]->GetNbinsX() + 1 ) );
        py6[j]->Scale( 1.0 / py6[j]->Integral( 0, py6[j]->GetNbinsX() + 1 ) );
        py8[j]->Scale( 1.0 / py8[j]->Integral( 0, py8[j]->GetNbinsX() + 1 ) );
        
        
        unf[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV", jetPtLo[j], jetPtHi[j] ) );
        unf[j]->SetXTitle( Form( "Q_{jet}  (#kappa = %1.1f)" , k) );
        
        
        unf[j]->SetYTitle("Prob.");
        unf[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
        
        
        
        unf[j]->Draw("p");
        py6[j]->Draw("p same");
        py8[j]->Draw("p same");
        
        leg[j]->Draw();
        drawText( Form( "R = %0.1f", R ), .65, .7, 30 );
        
        c[j]->SaveAs( "/Users/grantmcnamara/Documents/jet-charge/macros/plots/compareDataAndPythia_" + outFile_pt[j] + "_R" + rad + "_k" + kappa + ".png" , "RECREATE" );
        
    }
}


