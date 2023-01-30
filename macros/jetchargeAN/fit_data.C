// Grant McNamara 1/13/2023
//
// fit unfolded data with pythia-8 (for now) templates but while keeping s(sbar), t(tbar), c(cbar), b(bbar) fractions constant and u+ubar as a single template, d+dbar as a single template and then g as a template (i think this is how cms is done, testing if this is an improvement due to a dbar, ubar production mechanism messing with varying e.g. u and not ubar)
//
//
// perform fitting to unfolded data using pythia-8 templates
// and draw thstack containing fit results by flavor similar to how it is done in CMS jet charge paper (https://arxiv.org/pdf/2004.00602.pdf)
// and include ratio of fit result/histogram to be fit (pythia-8 here)




#include <string>
#include <iostream>
#include "math.h"
#include "../Headers/plot.h"
#include "../Headers/JCparameters.hh"


using namespace jcAnalysis;



// fit data with pythia-8 templates
void fit_data( double k = 0.0, double R = 0.4 ){
    // perform tfractionfitter fit to pythia-8 total distribution using pythia-8 templates
    cout << "\nfitting unfolded data with pythia-8 templates\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) ( ( k - (int) k ) * 10) );

    
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");//, "Minimize");
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Minimize");
    //ROOT::Math::MinimizerOptions::SetDefaultErrorDef( 0.5 );
    //ROOT::Math::MinimizerOptions::SetDefaultStrategy( 0 );
    //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );
    ROOT::Math::MinimizerOptions::PrintDefault();

    
    // 206, 691 jets in unfolded JP2 data with 20 < pT < 25 GeV
    
    TFile *fdat = new TFile( "../../out/unfolded_R" + rad + "_k" + kappa + ".root", "READ" );
    //unfolded_R04_k07.root
    
    //p8_undecayed_R04_k03.root
    TFile *fp8 = new TFile( "../../out/p8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" );
    
    
    //these files include pythia8_norm_ + outFile_pt[j] and the individual flavor templates
    // templates of the form "u_2025"


    
    const int nflavs = 3;
    TString flavs_names[nflavs+1] = {"d", "u", "g", "other"};

    // 13 = 6*quarks + 6*antiquarks + gluon
    // - nflavs = n flavor templates being used
    // - 1 = omit s as the start of "other" template
    TString other_flav[7] = {"c", "b", "t", "sbar", "cbar", "bbar", "tbar"};


    TH2D* p8_2d = (TH2D*) fp8->Get( "qvpt" ); // 2D histogram with full pythia-8 distribution of q vs pt
    TH1D* hpy8[njetbins]; // p8_2d projected onto 1D histograms


    double temp_fracs[nflavs+1]; // array of pythia-8 flavor fractions to use to constrain fitting if needed/used
    // templates don't get normalized, so this may be redundant?
    TH2D* tmp_2D[nflavs+1]; // 2D histograms by flavor
    TH1D* htemps[njetbins][nflavs+1]; // 2D histograms projected onto 1Ds

    
    //data histograms
    TH1D* hdat[njetbins];
    for(int j = 0; j < njetbins; j++){
        // non-normalized unfolded results (so that statistical errors don't fuck up the fitting)
        hdat[j] = (TH1D*) fdat->Get("unf_scale_" + outFile_pt[j] );
        
    }
    
    for(int i = 0; i < nflavs; i++){ // u,d,g
        
        tmp_2D[i] = (TH2D*) fp8->Get( flavs_names[i] + "_qvpt" );
        
        if( flavs_names[i] != "g" ){
            //cout << "\n" << flavs_names[i] << "bar added to " << flavs_names[i] << "\n\n";
            tmp_2D[i]->Add( (TH2D*) fp8->Get(flavs_names[i] + "bar_qvpt") ); // add e.g. u + ubar templates
        }

        for(int j = 0; j < njetbins; j++){
            if(i == 0){
                hpy8[j] = (TH1D*) p8_2d->ProjectionX( "p8_" + outFile_pt[j], p8_2d->GetYaxis()->FindBin( jetPtLo[j] ), p8_2d->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
                
            }

            htemps[j][i] = (TH1D*) tmp_2D[i]->ProjectionX( flavs_names[i] + "_" + outFile_pt[j], tmp_2D[i]->GetYaxis()->FindBin( jetPtLo[j] ), tmp_2D[i]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        }
    }
    
    
    // "other" 2D histogram, to be projected onto 1D histogram
    tmp_2D[nflavs] = (TH2D*) fp8->Get( "s_qvpt" ); // start "other" template with s, then loop through c,b,t, sbar,cbar,bbar,tbar
    for(int i = 0; i < 7; i++){
        tmp_2D[nflavs]->Add( (TH2D*) fp8->Get( other_flav[i] + "_qvpt" ) );
    }

    // project 2D "other" (aka s template) into 1D "other" list of histogram
    for(int j = 0; j < njetbins; j++){
        htemps[j][nflavs] = (TH1D*) tmp_2D[nflavs]->ProjectionX( "other_" + outFile_pt[j], tmp_2D[nflavs]->GetYaxis()->FindBin( jetPtLo[j] ), tmp_2D[nflavs]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
    }


    // for tfractionfitter
    TObjArray *h_fit[njetbins];
    TFractionFitter *fit[njetbins];

    TH1D* result_tff[njetbins];
    
    THStack *h_result[njetbins]; // stacking flavor histograms from fit result
    TH1D* fit_preds[njetbins][nflavs+1]; // nflavs + other; store fit[j]->GetMCPrediction(i) here
    
    
    TCanvas *c = new TCanvas("c_dat", "canvas", 800, 500);
    c->SetLeftMargin(0.2);

    TLegend *leg1 = new TLegend(0.05, 0.7, 0.9, .85);
    leg1->SetNColumns(2);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.08);

    TLegend *leg3 = new TLegend(0.2, 0.65, 0.9, .75);
    leg3->SetNColumns(3);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.08);

    TLegend *leg2 = new TLegend(0.35, 0.75, 0.65, .85);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.08);

    c->Divide(njetbins, 1, 0, 0);


    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString padup_name = Form( "pdatstack_%i", j );
        TPad* padup = new TPad(padup_name, padup_name, 0, 0.4, 1, 1);

        padup->SetRightMargin(0);
        padup->SetTopMargin(0);
        padup->SetBottomMargin(0.0);
        padup->SetLeftMargin(0);
        if( j == 0 ){padup->SetLeftMargin(0.2);}

        padup->Draw();
        padup->cd();

        h_fit[j] = new TObjArray(nflavs+1); // nflavs: u,d,g; template for all other flavors combined
        
        
        h_result[j] = new THStack("h_result_" + outFile_pt[j], "" );
        
        
        for(int i = 0; i < nflavs+1; i++){
            h_fit[j]->Add( htemps[j][i] );
        }


        double py8_tot = hpy8[j]->Integral( 0, hpy8[j]->GetNbinsX() + 1 );
        double dat_jets = hdat[j]->Integral( 0, hdat[j]->GetNbinsX() + 1 );

        
        hdat[j]->SetMarkerColor(kBlack);
        hdat[j]->SetMarkerStyle(kFullSquare);

        if( j == 0 ){
            leg2->AddEntry( hdat[j], "Corrected STAR", "p" );
        }
        

        // start to perform the fit on unfolded data
        fit[j] = new TFractionFitter( hdat[j], h_fit[j] , "Q");

        
        for(int i = 0; i < nflavs; i++){ // individual flavor fractions only limited to be between 0 and 1
            fit[j]->Constrain( i, 0.0, 1.0 );
        }

        // constrain "other" contribution to very close to pythia fraction
        fit[j]->Constrain( nflavs,
                         ( 1.0 - other_constraint ) * ( htemps[j][nflavs]->Integral() / py8_tot ),
                         ( 1.0 + other_constraint ) * ( htemps[j][nflavs]->Integral() / py8_tot ) ); // fix "other" fraction to known pythia-8 value

        
        int status = fit[j]->Fit();

        cout << "\n\nFor jets with " << jetPtLo[j] << " < pT < " << jetPtHi[j] << " GeV\n";
        cout << "kappa = " << k << "\n\n";

        cout << "TFractionFitter fit status: " << status << "\n";


        result_tff[j] = (TH1D*) fit[j]->GetPlot();
        result_tff[j]->SetName( "frac_fit_" + outFile_pt[j] );
        
        result_tff[j]->SetFillColor(kBlack);
        result_tff[j]->SetMarkerColor(kBlack);
        result_tff[j]->SetMarkerStyle(kOpenSquare);
        result_tff[j]->SetStats(0);
        
        
        
        hdat[j]->Scale( 1.0 / dat_jets );
        
        
        double norm = result_tff[j]->Integral( 0, result_tff[j]->GetNbinsX() + 1 );
        
        result_tff[j]->Scale( 1.0 / norm );

        result_tff[j]->Divide( hdat[j] );

        
        hdat[j]->SetStats(0);
//        hdat[j]->Sumw2(0);
        if(k == 0.0){
            hdat[j]->GetXaxis()->SetRangeUser( -4.5, 4.5 );
            hdat[j]->GetYaxis()->SetRangeUser( 0.001, 0.699 );
        }
        else{
            hdat[j]->GetXaxis()->SetRangeUser( -1, 1 );
            if(k == 0.3){
                hdat[j]->GetYaxis()->SetRangeUser( 0.001, 0.209 );
            }
            if(k == 0.5){
                hdat[j]->GetYaxis()->SetRangeUser( 0.001, 0.329 );
            }
            if(k == 0.7){
                hdat[j]->GetYaxis()->SetRangeUser( 0.001, 0.349 );
            }
        }
        

        hdat[j]->GetYaxis()->SetTitleFont(43);
        hdat[j]->GetYaxis()->SetTitleSize(22);
        hdat[j]->GetYaxis()->SetTitleOffset(1.3);
        hdat[j]->GetYaxis()->SetLabelFont(43);
        hdat[j]->GetYaxis()->SetLabelSize(18);
        hdat[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );


        hdat[j]->GetXaxis()->SetTitleFont(43);
        hdat[j]->GetXaxis()->SetTitleSize(22);
        hdat[j]->GetXaxis()->SetTitleOffset(1.1);
        hdat[j]->GetXaxis()->SetLabelFont(43);
        hdat[j]->GetXaxis()->SetLabelSize(18);
        hdat[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );

        
        hdat[j]->Draw("PSAME");

        
        for(int i = 0; i < nflavs+1; i++){
            fit_preds[j][i] = (TH1D*) fit[j]->GetMCPrediction(i);

            fit_preds[j][i]->SetFillColor( mark_colors[i] );
            fit_preds[j][i]->SetLineColor( mark_colors[i] );

            fit_preds[j][i]->Scale( 1.0 / py8_tot ); // normalize by total number of jets?
            
        }
        
        // add in same order as cms just to be consistent for now -- 11/16/22
        // other, down, up, gluon
        h_result[j]->Add(fit_preds[j][nflavs]);
        h_result[j]->Add(fit_preds[j][0]);
        h_result[j]->Add(fit_preds[j][1]);
        h_result[j]->Add(fit_preds[j][2]);

        h_result[j]->Draw("HISTSAME"); // draw histogram stack
        hdat[j]->DrawCopy("PSAME");
        
        if( j == 0 ){
            drawText( "Corrected Data", .35, .93, 20 );

            leg1->AddEntry( fit_preds[j][0], "d/#bar d", "f" );
            //leg3->AddEntry( (TObject*)0, Form( "%2.0f %%", frac ), "" );
            
            leg1->AddEntry( fit_preds[j][1], "u/#bar u", "f" );
            //leg3->AddEntry( (TObject*)0, Form( "%2.0f %%", frac ), "" );
            
            leg1->AddEntry( fit_preds[j][2], "gluon", "f" );
            //leg3->AddEntry( (TObject*)0, Form( "%2.0f %%", frac ), "" );
            
            leg1->AddEntry( fit_preds[j][3], "other", "f" );
            
        }
        
        if( j == 0 ){
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .36, .87, 15 );
            //drawText( "p_{T}^{track} > " + tr_pt + " MeV", .38, .78, 15 );
        }
        if( j == 1 ){
            leg1->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .22, .87, 15 );
            //drawText( Form( "anti-k_{T}, R = %0.1f", R ), .35, .93, 20 );
        }
        if( j == 2 ){
            leg2->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .22, .87, 15 );
            drawText( "p+p #sqrt{s} = 200 GeV", .25, .93, 15 );
        }
        
        c->cd(j+1);
        
        TString pdown_name = Form( "pbottom_%i", j );
        TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0.05, 1, 0.4);

        p2->SetRightMargin(0);
        p2->SetTopMargin(0);
        p2->SetBottomMargin(0.4);
        p2->SetLeftMargin(0);
        if( j == 0 ){p2->SetLeftMargin(0.2);}

        p2->Draw();
        p2->cd();


        
        result_tff[j]->GetYaxis()->SetRangeUser(0.701, 1.299);
        if(k == 0.0){
            result_tff[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            result_tff[j]->GetXaxis()->SetRangeUser(-1, 1);
        }
        
        result_tff[j]->GetYaxis()->SetTitleFont(43);
        result_tff[j]->GetYaxis()->SetTitleSize(22);
        result_tff[j]->GetYaxis()->SetTitleOffset(1.2);
        result_tff[j]->GetYaxis()->SetNdivisions(505);
        result_tff[j]->GetYaxis()->SetLabelFont(43);
        result_tff[j]->GetYaxis()->SetLabelSize(18);
        result_tff[j]->GetYaxis()->SetTitle( "Fit / Data" );

        result_tff[j]->GetXaxis()->SetTitleFont(43);
        result_tff[j]->GetXaxis()->SetTitleSize(22);
        result_tff[j]->GetXaxis()->SetTitleOffset(1.1);
        result_tff[j]->GetXaxis()->SetLabelFont(43);
        result_tff[j]->GetXaxis()->SetLabelSize(16);
        result_tff[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        result_tff[j]->Draw("E1SAME");
        // ratio of fit result / "data"
        
        
        drawText( Form("fit status: %1.0f", (double) status), 0.45, 0.8, 12 );
        
        for(int i = 0; i < nflavs+1; i++){ // loop over actual flavors (u,d,g)

            
            double frac = 1.0;// set to parameter from the fit result
            double err = 1.0;

            fit[j]->GetResult( i, frac, err );


            cout << "\npythia-8 flavor " << flavs_names[i] << " fraction: " << htemps[j][i]->Integral() / py8_tot << "\n";
            cout << "fit result for " << flavs_names[i] << " fraction: " << frac << " +/- " << err << "\n";
            
        }
        cout << "\n";
    }// j loop over jet pt bins
    c->SaveAs( "./figs/fitUnfoldedData_uplusubar_R" + rad + "_k" + kappa + ".pdf", "RECREATE" );
}

