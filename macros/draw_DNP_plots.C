// 9/29/22 --- Grant McNamara
// macro to recreate all necessary plots for DNP 2022 for jet charge analysis
//
// draw data by itself, data and det-level pythia-6 together, part-level pythia-6 and unfolded data (with systematics) together, then all four together
//


#include <string>
#include <iostream>
#include "math.h"
#include "Headers/plot.h"
#include "Headers/JCparameters.hh"

using namespace jcAnalysis;


void getSystematics(TH1D *syst[nsyst][njetbins], double k = 0.0, double R = 0.4){
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    
    
    TFile *fin = new TFile( "../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    
    for(int j = 0; j < njetbins; j++){
        for(int ki = 0; ki < nsyst; ki++){
            syst[ki][j] = (TH1D*) fin->Get( "unfold_" + syst_names[ki] + outFile_pt[j] );
            
            /*if(ki > 4){
                cout << "unfold_" + syst_names[ki] + outFile_pt[j] << "\n";
                cout << syst[ki][j]->Integral() << "\n";
            }*/
        }
    }
    cout << "\n";
}


// calculate quadrature sum of systematic errors from all sources and set bin error of q_err of j-th jet pt bin histogram to the sqrt(quad_sum)
void systErrs( TH1D *syst[nsyst][njetbins], TH1D *q_err[njetbins], int j ){
    // calculate errors within jet pt loop (calculate at each j as the jet pt bin loops
    // set errors on q_err histogram to the quadrature sum of the systematic errors multiplied by the bin content of q_err (or the unf histogram in the main code, should be equivalent)
    
    for(int i = 1; i < (int) syst[0][j]->GetNbinsX(); i++){
        double err_quadsum = 0.0;
        
        for(int ki = 0; ki < nsyst; ki++){
            double syst_err = syst[ki][j]->GetBinContent(i);
            err_quadsum += syst_err*syst_err; // total error is sqrt(sum of individual errors ^2)
        }
        q_err[j]->SetBinError( i, sqrt( err_quadsum ) * q_err[j]->GetBinContent(i) );
        // err_quadsum is a percentage of the bin content, so to calculate absolute error have to multiply sqrt(err_quadsum) * bincontent
    }
    
}



// draw data jet charge for a given kappa, jet R on a TCanvas alone
void plot_jetcharge_dataOnly(double k = 0.0, double R = 0.4){
    cout << "\nPlotting data by itself\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    
    
    
    TH2D* data_2D = (TH2D*) fdat->Get( "QvPt_d" ); // 2D hist used as (raw) data spectrum when unfolding, unfolded dists will be 1Ds from same file (unfold.cxx writes out its input and outputs to same file)
    
    TH1D* data[njetbins];
    
    
    TCanvas *c[njetbins];
    TLegend *leg[njetbins];
    
    for(int j = 0; j < njetbins; j++){
        c[j] = new TCanvas(Form("dataONLY_%i", j), "canvas", 600, 600);
        leg[j] = new TLegend(0.18, 0.55, .4, .85);
        leg[j]->SetBorderSize(0);
        //leg[j]->SetNColumns(3);
        
        c[j]->SetLeftMargin(0.15);
        
        // name projections according to e.g. "data_jetQ_jetpt2025"
        data[j] = (TH1D*) data_2D->ProjectionX( "data_jetQ_jetpt" + outFile_pt[j], j+1, j+1, "e" );
        // or data_2D->GetYaxis()->FindBin( jetPtLo[j] ), data_2D->GetYaxis()->FindBin( jetPtHi[j+1] - 1 ), "e");
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        data[j]->Scale( 1.0 / data[j]->Integral( 0, data[j]->GetNbinsX() + 1 ) );
        data[j]->SetMarkerColor(kViolet-3);
        data[j]->SetMarkerStyle(kOpenSquare);
        
        
        leg[j]->AddEntry( data[j], "pp JP2 Data", "p");
        
        data[j]->SetStats(0);
        data[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
        
        data[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV", jetPtLo[j], jetPtHi[j] ) );
        data[j]->SetXTitle( Form( "Q_{#kappa = %1.1f}^{jet}" , k) );
        data[j]->SetYTitle( "prob." );
        
        data[j]->Draw("P");
        
        
        leg[j]->Draw();
        
        drawText(Form("R = %0.1f", R), .65, .7, 30);
        
        
        c[j]->SaveAs( "./plots/dnp/dataONLY_jetCharge_jetpt" + outFile_pt[j] + "_R" + rad + "_k" + kappa + ".png", "RECREATE");
    }
    
}



// draw data and pythia-6 det-level jet charge for a given kappa, jet R on TCanvas together
void plot_jetcharge_dataDet(double k = 0.0, double R = 0.4){
    cout << "\nPlotting data and PYTHIA+GEANT on one canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // file with pythia-6 distribution
    TFile *fp6 = new TFile( "../out/p6_R" + rad + "_k" + kappa + ".root" , "READ");
    
    
    TH2D* data_2D = (TH2D*) fdat->Get( "QvPt_d" );
    TH1D* data[njetbins];
    TH1D* p6_d[njetbins]; // det-level distributions
    
    
    TCanvas *c[njetbins];
    TLegend *leg[njetbins];
    
    
    for(int j = 0; j < njetbins; j++){
        c[j] = new TCanvas(Form("dataDet_%i", j), "canvas", 600, 600);
        leg[j] = new TLegend(0.18, 0.55, .4, .85);
        leg[j]->SetBorderSize(0);
        //leg[j]->SetNColumns(3);
        
        c[j]->SetLeftMargin(0.15);
    
        
        // name projections according to e.g. "data_jetQ_jetpt2025"
        data[j] = (TH1D*) data_2D->ProjectionX( "data_jetQ_jetpt" + outFile_pt[j], j+1, j+1, "e" );
        // or data_2D->GetYaxis()->FindBin( jetPtLo[j] ), data_2D->GetYaxis()->FindBin( jetPtHi[j+1] - 1 ), "e");
        
        p6_d[j] = (TH1D*) fp6->Get( Form( "jetQ_k" + kappa + "_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j] ) );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        data[j]->Scale( 1.0 / data[j]->Integral( 0, data[j]->GetNbinsX() + 1 ) );
        p6_d[j]->Scale( 1.0 / p6_d[j]->Integral( 0, p6_d[j]->GetNbinsX() + 1 ) );
        
        data[j]->SetMarkerColor(kViolet-3);
        p6_d[j]->SetMarkerColor(kGreen+2);
        
        
        data[j]->SetMarkerStyle(kOpenSquare);
        p6_d[j]->SetMarkerStyle(kOpenSquare);
        
        
        leg[j]->AddEntry( data[j], "pp JP2 Data", "p");
        leg[j]->AddEntry( p6_d[j], "Pythia-6 + GEANT", "p");

        
        data[j]->SetStats(0);
        data[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
        
        data[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV", jetPtLo[j], jetPtHi[j] ) );
        data[j]->SetXTitle( Form( "Q_{#kappa = %1.1f}^{jet}" , k) );
        data[j]->SetYTitle( "prob." );
        
        
        data[j]->Draw("P");
        p6_d[j]->Draw("P SAME");
        
        
        leg[j]->Draw();
        
        drawText(Form("R = %0.1f", R), .65, .7, 30);
        
        
        c[j]->SaveAs( "./plots/dnp/dataAndDetector_jetCharge_jetpt" + outFile_pt[j] + "_R" + rad + "_k" + kappa + ".png", "RECREATE");
    }
    
}



// draw unfolded data with systematic uncertainties and pythia-6 part-level jet charge for a given kappa, jet R on TCanvas together
void plot_jetcharge_unfPart(double k = 0.0, double R = 0.4){
    cout << "\nPlotting unfolded data and PYTHIA on one canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // file with pythia-6 distribution
    TFile *fp6 = new TFile( "../out/p6_R" + rad + "_k" + kappa + ".root" , "READ");

    
    
    TH1D* unfs[njetbins]; // unfolded distributions go here
    TH1D* p6_p[njetbins]; // part-level distriubtions
    TH1D* uerr[njetbins]; // for systematic errors
    
    TH1D* systs[nsyst][njetbins]; // individual systematic histograms as percentages
    
    // get systematics histograms-->separated out into its own function
    getSystematics( systs ); // pass address of 2D array of TH1Ds so that no copy is made, the copy here is actually filled
    
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
        
        //cout << "unfold_errs_" + outFile_pt[j] << "\nEXISTS? " << uerr[j]->Integral() << "\n";
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        unfs[j]->Scale( 1.0 / unfs[j]->Integral( 0, unfs[j]->GetNbinsX() + 1 ) );
        p6_p[j]->Scale( 1.0 / p6_p[j]->Integral( 0, p6_p[j]->GetNbinsX() + 1 ) );

        
        // clone normalized unfolded distributions (to be safe)
        uerr[j] = (TH1D*) unfs[j]->Clone( "unfold_errs_" + outFile_pt[j] );
        
        // calculate quadrature sum of errors from all systematic sources and set bin error of uerr[j] to that sum * bin content (quad sum is % error)
        systErrs(systs, uerr, j);
        
        
        unfs[j]->SetMarkerColor(kViolet-3);
        p6_p[j]->SetMarkerColor(kGreen+2);
        uerr[j]->SetFillColorAlpha(kViolet-3, 0.25);
        
        
        unfs[j]->SetMarkerStyle(kFullCircle);
        p6_p[j]->SetMarkerStyle(kFullCircle);
        uerr[j]->SetFillStyle(3444);
        
        
        leg[j]->AddEntry( unfs[j], "Unfolded Data", "p");
        leg[j]->AddEntry( p6_p[j], "Pythia-6", "p");

        
        unfs[j]->SetStats(0);
        unfs[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
        
        unfs[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV", jetPtLo[j], jetPtHi[j] ) );
        unfs[j]->SetXTitle( Form( "Q_{#kappa = %1.1f}^{jet}" , k) );
        unfs[j]->SetYTitle( "prob." );
        
        unfs[j]->Draw("P");
        p6_p[j]->Draw("P SAME");
        uerr[j]->Draw("E2 SAME");
        
        
        leg[j]->Draw();
        
        drawText(Form("R = %0.1f", R), .65, .7, 30);
        
        
        c[j]->SaveAs( "./plots/dnp/unfoldedAndParticle_jetCharge_jetpt" + outFile_pt[j] + "_R" + rad + "_k" + kappa + ".png", "RECREATE");
        
    }
}


// draw all 4: data, unfolded data (with systematic uncertainties included), det-, and part-level jet charge together on same canvas
void plot_jetcharge_allFour(double k = 0.0, double R = 0.4){
    cout << "\nPlotting data, unfolded data, PYTHIA, and PYTHIA+GEANT on one canvas\n";
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // file with pythia-6 distribution
    TFile *fp6 = new TFile( "../out/p6_R" + rad + "_k" + kappa + ".root" , "READ");
    
    
    TH2D* data_2D = (TH2D*) fdat->Get( "QvPt_d" );
    TH1D* data[njetbins];
    TH1D* unfs[njetbins]; // unfolded distributions go here
    TH1D* p6_p[njetbins]; // part-level distriubtions
    TH1D* p6_d[njetbins]; // det-level distributions
    
    // error-related histograms
    TH1D* uerr[njetbins]; // for systematic errors
    TH1D* systs[nsyst][njetbins]; // individual systematic histograms as percentages
    // get systematics histograms-->separated out into its own function
    getSystematics( systs ); // pass address of 2D array of TH1Ds so that no copy is made, the copy here is actually filled
    
    
    TCanvas *c[njetbins];
    TLegend *leg[njetbins];
    
    
    for(int j = 0; j < njetbins; j++){
        c[j] = new TCanvas(Form("all4_%i", j), "canvas", 600, 600);
        leg[j] = new TLegend(0.18, 0.55, .4, .85);
        leg[j]->SetBorderSize(0);
        
        
        c[j]->SetLeftMargin(0.15);
        
        // name projections according to e.g. "data_jetQ_jetpt2025"
        data[j] = (TH1D*) data_2D->ProjectionX( "data_jetQ_jetpt" + outFile_pt[j], j+1, j+1, "e" );
        // or data_2D->GetYaxis()->FindBin( jetPtLo[j] ), data_2D->GetYaxis()->FindBin( jetPtHi[j+1] - 1 ), "e");
        
        unfs[j] = (TH1D*) fdat->Get( "unfold_nomx" + outFile_pt[j] );
        
        p6_p[j] = (TH1D*) fp6->Get( Form( "jetQ_k" + kappa + "_part_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j] ) );
        p6_d[j] = (TH1D*) fp6->Get( Form( "jetQ_k" + kappa + "_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j] ) );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        data[j]->Scale( 1.0 / data[j]->Integral( 0, data[j]->GetNbinsX() + 1 ) );
        unfs[j]->Scale( 1.0 / unfs[j]->Integral( 0, unfs[j]->GetNbinsX() + 1 ) );
        p6_p[j]->Scale( 1.0 / p6_p[j]->Integral( 0, p6_p[j]->GetNbinsX() + 1 ) );
        p6_d[j]->Scale( 1.0 / p6_d[j]->Integral( 0, p6_d[j]->GetNbinsX() + 1 ) );
     
        
        // clone normalized unfolded distributions (to be safe)
        uerr[j] = (TH1D*) unfs[j]->Clone( "unfold_errs_" + outFile_pt[j] );
        // calculate quadrature sum of errors from all systematic sources and set bin error of uerr[j] to that sum * bin content (quad sum is % error)
        systErrs(systs, uerr, j);
        
        
        
        data[j]->SetMarkerColor(kViolet-3);
        p6_d[j]->SetMarkerColor(kGreen+2);
        unfs[j]->SetMarkerColor(kViolet-3);
        p6_p[j]->SetMarkerColor(kGreen+2);
        uerr[j]->SetFillColorAlpha(kViolet-3, 0.25);
        
        
        data[j]->SetMarkerStyle(kOpenSquare);
        p6_d[j]->SetMarkerStyle(kOpenSquare);
        unfs[j]->SetMarkerStyle(kFullCircle);
        p6_p[j]->SetMarkerStyle(kFullCircle);
        uerr[j]->SetFillStyle(3444);
        
        
        leg[j]->AddEntry( data[j], "pp JP2 Data", "p");
        leg[j]->AddEntry( unfs[j], "Unfolded Data", "p");
        leg[j]->AddEntry( p6_p[j], "Pythia-6", "p");
        leg[j]->AddEntry( p6_d[j], "Pythia-6 + GEANT", "p");

        
        data[j]->SetStats(0);
        data[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
        
        data[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV", jetPtLo[j], jetPtHi[j] ) );
        data[j]->SetXTitle( Form( "Q_{#kappa = %1.1f}^{jet}" , k) );
        data[j]->SetYTitle( "prob." );
        
        
        data[j]->Draw("P");
        unfs[j]->Draw("P SAME");
        p6_p[j]->Draw("P SAME");
        p6_d[j]->Draw("P SAME");
        uerr[j]->Draw("E2 SAME");
        
        
        leg[j]->Draw();
        
        drawText(Form("R = %0.1f", R), .65, .7, 30);
        
        
        c[j]->SaveAs( "./plots/dnp/allfour_jetCharge_jetpt" + outFile_pt[j] + "_R" + rad + "_k" + kappa + ".png", "RECREATE");
    }
    
}



// MAIN MACRO
// Draw plots by calling macros defined above
void draw_DNP_plots(double k = 0.0, double R = 0.4){
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    plot_jetcharge_dataOnly( k, R );
    plot_jetcharge_dataDet( k, R );
    plot_jetcharge_unfPart( k, R );
    plot_jetcharge_allFour( k, R );
    
    
}
