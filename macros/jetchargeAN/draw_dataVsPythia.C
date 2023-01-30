// 9/29/22 --- Grant McNamara
// macro to recreate all necessary plots for DNP 2022 for jet charge analysis
//
// draw data by itself, data and det-level pythia-6 together, part-level pythia-6 and unfolded data (with systematics) together, then all four together
//


#include <string>
#include <iostream>
#include "math.h"
#include "../Headers/plot.h"
#include "../Headers/JCparameters.hh"

using namespace jcAnalysis;



// draw unfolded data with systematic uncertainties and pythia-6 part-level jet charge for a given kappa, jet R on TCanvas together
void plot_jetcharge_unfPart(double k = 0.0, double R = 0.4){
    cout << "\nPlotting unfolded data and PYTHIA on one canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "../../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // file with pythia-6 distribution
    //p6_unmatched_R04_k00.root
    TFile *fp6 = new TFile( "../../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");

    TFile *fp8 = new TFile( "../../out/p8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" );
    //TFile *fh7 = new TFile( "../../out/h7_undecayed_R" + rad + "_k" + kappa + ".root", "READ" );

    
    TH1D* unfs[njetbins]; // unfolded distributions go here
    TH1D* p6_p[njetbins]; // part-level distriubtions
    TH1D* uerr[njetbins]; // for systematic errors
    
    
    // 2D qvpt hists from pythia-8, herwig-7 that need to be projected by jet pt bin
    TH2D* hp8 = (TH2D*) fp8->Get( "qvpt" );
    //TH2D* hh7 = (TH2D*) fh7->Get( "qvpt" );
    
    TH1D* py8[njetbins];
    TH1D* hw7[njetbins];
    
    
    
    TCanvas *c[njetbins];
    TLegend *leg[njetbins];
    
    
    for(int j = 0; j < njetbins; j++){
        
        c[j] = new TCanvas(Form("unfPart_%i", j), "canvas", 600, 600);
        leg[j] = new TLegend(0.18, 0.55, .4, .85);
        leg[j]->SetBorderSize(0);
        //leg[j]->SetNColumns(3);
        
        c[j]->SetLeftMargin(0.15);
        
        //cout << "GS systematic? " << systs[6][j]->Integral() << "\n";
        
        unfs[j] = (TH1D*) fdat->Get( "unfold_nomx" + outFile_pt[j] );
        
        
        p6_p[j] = (TH1D*) fp6->Get( Form( "jetQ_k" + kappa + "_part_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j] ) );
        
        py8[j] = (TH1D*) hp8->ProjectionX( "p8_" + outFile_pt[j], hp8->GetYaxis()->FindBin( jetPtLo[j] ), hp8->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        //hw7[j] = (TH1D*) hh7->ProjectionX( "p8_" + outFile_pt[j], hh7->GetYaxis()->FindBin( jetPtLo[j] ), hh7->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        //cout << "unfold_errs_" + outFile_pt[j] << "\nEXISTS? " << uerr[j]->Integral() << "\n";
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        unfs[j]->Scale( 1.0 / unfs[j]->Integral( 0, unfs[j]->GetNbinsX() + 1 ) );
        p6_p[j]->Scale( 1.0 / p6_p[j]->Integral( 0, p6_p[j]->GetNbinsX() + 1 ) );
        
        py8[j]->Scale( 1.0 / py8[j]->Integral( 0, py8[j]->GetNbinsX() + 1 ) );
        //hw7[j]->Scale( 1.0 / hw7[j]->Integral( 0, hw7[j]->GetNbinsX() + 1 ) );

        
        
        unfs[j]->SetLineColor(kViolet-3);
        p6_p[j]->SetLineColor(kGreen+2);
        py8[j]->SetLineColor(kGreen+2);
        //hw7[j]->SetLineColor(kRed);
        
        
        unfs[j]->SetLineStyle(kSolid);
        p6_p[j]->SetLineStyle(kSolid);
        py8[j]->SetLineStyle(kDashed);
        //hw7[j]->SetLineStyle(kDashed);
        
        
        leg[j]->AddEntry( unfs[j], "Unfolded Data", "l" );
        leg[j]->AddEntry( p6_p[j], "Pythia-6", "l" );
        leg[j]->AddEntry( py8[j], "Pythia-8 Perugia", "l" );
        //leg[j]->AddEntry( hw7[j], "Herwig-7", "l" );

        
        unfs[j]->SetStats(0);
        unfs[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
        
        unfs[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV", jetPtLo[j], jetPtHi[j] ) );
        unfs[j]->SetXTitle( Form( "Q_{#kappa = %1.1f}^{jet}" , k) );
        unfs[j]->SetYTitle( "1/N_{jets} dN/dQ [1/e]" );
        
        unfs[j]->Draw("C");
        p6_p[j]->Draw("CSAME");
        py8[j]->Draw("CSAME");
        //hw7[j]->Draw("CSAME");
        
        
        leg[j]->Draw();
        
        drawText(Form("R = %0.1f", R), .65, .7, 30);
        
        
        c[j]->SaveAs( "./figs/jetcharge_unfoldedVsPythiaAndHerwig_jetpt" + outFile_pt[j] + "_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
        
    }
}



// MAIN MACRO
// Draw plots by calling macros defined above
void draw_dataVsPythia(double k = 0.0, double R = 0.4){
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    plot_jetcharge_unfPart( k, R );
    
}
