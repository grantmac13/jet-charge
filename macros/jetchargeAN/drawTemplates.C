// 1/25/2023 Grant McNamara
// draw templates as curves with lines through bins
// solid lines corresponding to Monash (old tune) and dashed lines corresponding to Detroit (new tune) Tune templates to see if there is improvement


#include <string>
#include <iostream>
#include "math.h"
#include "../Headers/plot.h"
#include "../Headers/JCparameters.hh"


using namespace jcAnalysis;


// using the template drawing funciton from DNP macro
void draw_templates_MonashVsDetroit( double k = 0.0, double R = 0.4, int width = 4 ){
    // draw u,d,g templates and pythia-8 total distributions
    cout << "\ncomparing templates from Pythia-8 Monash and Detroit Tunes\n\n";
    
    cout << "what is kappa? k = " << k << "\n\n";
    
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n\n";
    
    double ymax = jc_padmax[0]; // default to k = 0.0 window
    for(int ik = 0; ik < nkaps; ik++){
        if( k == k_vals[ik] ){ ymax = jc_padmax[ik]; }
    }
    
    // adding pythia-8 unmatched to this figure
    // filename = p8_undecayed_R04_k00.root
    TFile *fp8m = new TFile( "../../out/p8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" ); // monash tune
    //this file includes pythia8_norm_ + outFile_pt[j] and the individual flavor templates
    // templates of the form "u_2025"

    TFile *fp8d = new TFile( "../../out/detroit_tune_pythia8_decayed_R" + rad + "_k" + kappa + ".root", "READ" ); // detroit tune, with decays

    TFile *fp8n = new TFile( "../../out/detroit_tune_pythia8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" ); // detroit tune, without decays
    
    /*
    cout << "Monash file name = " << fp8m->GetName() << "\n";
    cout << "Detroit file name = " << fp8d->GetName() << "\n\n";
    */
    
    const int nflavs = 3;
    TString flavs_names[nflavs] = {"d", "u", "g"};
    
    
    // non-normalized
    TH2D* hflavs_m[nflavs];// = (TH2D*) fp8m->Get("");; // flavor_qvpt
    TH1D* htemps_m[njetbins][nflavs]; // projection of hflavs_m onto Q axis by jet pt bin
    TH2D* hpy8_m = (TH2D*) fp8m->Get("qvpt"); // qvpt
    TH1D* p8_monash[njetbins]; // projection of hpy8_m onto Q axis by jet pt bin

    // non-normalized--detroit
    TH2D* hflavs_d[nflavs]; // flavor_qvpt
    TH1D* htemps_d[njetbins][nflavs]; // projection of hflavs_d onto Q axis by jet pt bin
    TH2D* hpy8_d = (TH2D*) fp8d->Get("qvpt");; // qvpt
    TH1D* p8_detroit[njetbins]; // projection of hpy8_d onto Q axis by jet pt bin
    
    // non-normalized--detroit, without decays
    TH2D* hflavs_n[nflavs]; // flavor_qvpt
    TH1D* htemps_n[njetbins][nflavs]; // projection of hflavs_d onto Q axis by jet pt bin
    TH2D* hpy8_n = (TH2D*) fp8n->Get("qvpt");; // qvpt
    TH1D* p8_detroit_n[njetbins]; // projection of hpy8_d onto Q axis by jet pt bin
    
    for(int i = 0; i < nflavs; i++){
        hflavs_m[i] = (TH2D*) fp8m->Get( flavs_names[i] + "_qvpt" );
        hflavs_d[i] = (TH2D*) fp8d->Get( flavs_names[i] + "_qvpt" );
        hflavs_n[i] = (TH2D*) fp8n->Get( flavs_names[i] + "_qvpt" );
    }
    
    
    TCanvas *c = new TCanvas("c_temps", "canvas", 800, 500);
    c->SetLeftMargin(0.25);
    
    TLegend *leg1 = new TLegend(0.35, 0.75, 0.65, .85);
    //leg1->SetNColumns(3);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.08);
    
    TLegend *leg2 = new TLegend(0.35, 0.75, 0.65, .85);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.08);
    
    TLegend *leg3 = new TLegend(0.35, 0.75, 0.65, .85);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.08);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        //cout << "DEBUG: jet pt bin: " << j << " for [ " << jetPtLo[j] << ", " << jetPtHi[j] << "]\n";
        c->cd( j+1 );
        TString padup_name = Form( "ptemps_%i", j );
        TPad* padup = new TPad(padup_name, padup_name, 0, 0, 1, 1);
        
        padup->SetRightMargin(0);
        padup->SetTopMargin(0);
        padup->SetBottomMargin(0.2);
        padup->SetLeftMargin(0);
        if( j == 0 ){padup->SetLeftMargin(0.2);}
        
        padup->Draw();
        padup->cd();
        
        
        // project the total pythia 8 distributions from the qvpt 2D histogram onto slices of jet pt bins
        p8_monash[j] = (TH1D*) hpy8_m->ProjectionX( "monash_pythia8_" + outFile_pt[j], hpy8_m->GetYaxis()->FindBin( jetPtLo[j] ), hpy8_m->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );

        p8_detroit[j] = (TH1D*) hpy8_d->ProjectionX( "detroit_pythia8_" + outFile_pt[j], hpy8_d->GetYaxis()->FindBin( jetPtLo[j] ), hpy8_d->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        p8_detroit_n[j] = (TH1D*) hpy8_n->ProjectionX( "detroit_pythia8_nodecays_" + outFile_pt[j], hpy8_n->GetYaxis()->FindBin( jetPtLo[j] ), hpy8_n->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        
        p8_monash[j]->Scale( 1.0 / p8_monash[j]->Integral( 0, p8_monash[j]->GetNbinsX() + 1 ) );
        p8_detroit[j]->Scale( 1.0 / p8_detroit[j]->Integral( 0, p8_detroit[j]->GetNbinsX() + 1 ) );
        p8_detroit_n[j]->Scale( 1.0 / p8_detroit_n[j]->Integral( 0, p8_detroit_n[j]->GetNbinsX() + 1 ) );
        
        
        // project each flavor for current jet pt bin
        for(int i = 0; i < nflavs; i++){
            //cout << "DEBUG: Flavor = " << flavs_names[i] << "\n";
            
            htemps_m[j][i] = (TH1D*) hflavs_m[i]->ProjectionX( "monash_" + flavs_names[i] + "_" + outFile_pt[j], hflavs_m[i]->GetYaxis()->FindBin( jetPtLo[j] ), hflavs_m[i]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
            htemps_d[j][i] = (TH1D*) hflavs_d[i]->ProjectionX( "detroit_" + flavs_names[i] + "_" + outFile_pt[j], hflavs_d[i]->GetYaxis()->FindBin( jetPtLo[j] ), hflavs_d[i]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
            htemps_n[j][i] = (TH1D*) hflavs_n[i]->ProjectionX( "detroit_nodecays_" + flavs_names[i] + "_" + outFile_pt[j], hflavs_n[i]->GetYaxis()->FindBin( jetPtLo[j] ), hflavs_n[i]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
            
            
            htemps_m[j][i]->Scale( 1.0 / htemps_m[j][i]->Integral( 0, htemps_m[j][i]->GetNbinsX() + 1 ) );
            htemps_d[j][i]->Scale( 1.0 / htemps_d[j][i]->Integral( 0, htemps_d[j][i]->GetNbinsX() + 1 ) );
            htemps_n[j][i]->Scale( 1.0 / htemps_n[j][i]->Integral( 0, htemps_n[j][i]->GetNbinsX() + 1 ) );
        }
        
        
        
        
        p8_monash[j]->SetStats(0);
        p8_monash[j]->GetXaxis()->SetRangeUser( -4.5, 4.5 );
        p8_monash[j]->GetYaxis()->SetRangeUser( 0.001, ymax );
        p8_monash[j]->Sumw2(0);

        p8_detroit[j]->Sumw2(0);
        p8_detroit_n[j]->Sumw2(0);
        
//        p8_monash[j]->GetYaxis()->SetTitleSize(0.1);
        p8_monash[j]->GetYaxis()->SetTitleFont(43);
        p8_monash[j]->GetYaxis()->SetTitleSize(22);
        p8_monash[j]->GetYaxis()->SetTitleOffset(0.8);
//        p8_monash[j]->GetYaxis()->SetLabelSize(0.05);
        p8_monash[j]->GetYaxis()->SetLabelFont(43);
        p8_monash[j]->GetYaxis()->SetLabelSize(22);
        p8_monash[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );

        
        p8_monash[j]->GetXaxis()->SetTitleFont(43);
        p8_monash[j]->GetXaxis()->SetTitleSize(22);
        p8_monash[j]->GetXaxis()->SetTitleOffset(1.2);
        //p8_monash[j]->GetXaxis()->SetNdivisions(505);
        //p8_monash[j]->GetXaxis()->SetLabelSize(0.05);
        p8_monash[j]->GetXaxis()->SetLabelFont(43);
        p8_monash[j]->GetXaxis()->SetLabelSize(22);
        p8_monash[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        p8_monash[j]->SetLineColor(kBlack);
        p8_monash[j]->SetLineStyle(kSolid);
        p8_monash[j]->SetLineWidth(width);
        
        p8_detroit[j]->SetLineColor(kBlack);
        p8_detroit[j]->SetLineStyle(kDashed);
        p8_detroit[j]->SetLineWidth(width);
        
        p8_detroit_n[j]->SetLineColor(kBlack);
        p8_detroit_n[j]->SetLineStyle(kDotted);
        p8_detroit_n[j]->SetLineWidth(width);
        
        if( j == 0 ){
            leg1->AddEntry( p8_monash[j], "Monash", "l" );
            leg2->AddEntry( p8_detroit[j], "Detroit with decays", "l" );
            leg3->AddEntry( p8_detroit_n[j], "Detroit no decays", "l" );
        }
        
        p8_monash[j]->Draw("CSAME");
        p8_detroit[j]->Draw("CSAME");
        p8_detroit_n[j]->Draw("CSAME");
        
        for(int i = 0; i < nflavs; i++){
            htemps_m[j][i]->SetLineColor( mark_colors[i] );
            htemps_m[j][i]->SetLineStyle( kSolid );
            htemps_m[j][i]->SetLineWidth(width);
            htemps_m[j][i]->Sumw2(0);
            
            htemps_d[j][i]->SetLineColor( mark_colors[i] );
            htemps_d[j][i]->SetLineStyle( kDashed );
            htemps_d[j][i]->SetLineWidth(width);
            htemps_d[j][i]->Sumw2(0);

            htemps_n[j][i]->SetLineColor( mark_colors[i] );
            htemps_n[j][i]->SetLineStyle( kDotted );
            htemps_n[j][i]->SetLineWidth(width);
            htemps_n[j][i]->Sumw2(0);
            
            if( j == 0 ){
                leg1->AddEntry( htemps_m[j][i], flavs_names[i], "l" );
                leg2->AddEntry( htemps_d[j][i], flavs_names[i], "l" );
                leg3->AddEntry( htemps_n[j][i], flavs_names[i], "l" );
            }

//            if(i!=2){continue;}
            htemps_m[j][i]->Draw("CSAME");
            htemps_d[j][i]->Draw("CSAME");
            htemps_n[j][i]->Draw("CSAME");

        }
        
        /*
        cout << "DEBUG: DISTRIBUTION CHARACTERISTICS\n";
        for(int i = 0; i < nflavs; i++){
            cout << "MONASH\n";
            cout << flavs_names[i] << " mean = " << htemps_m[j][i]->GetMean() << " and std dev. = " << htemps_m[j][i]->GetStdDev() << "\n";
            
            cout << "DETROIT\n";
            cout << flavs_names[i] << " mean = " << htemps_d[j][i]->GetMean() << " and std dev. = " << htemps_d[j][i]->GetStdDev() << "\n";
        }
        */
        
        
        if( j == 0 ){
            leg3->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .36, .87, 15 );
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
        
        
        TLine *line = new TLine( 0, 0.001, 0, ymax );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        
        //line->Draw();

    }
    c->SaveAs( "./figs/templatesMonashVsDetroit_R" + rad + "_k" + kappa + ".pdf", "RECREATE" );
}



// from /Users/grantmcnamara/Documents/jet-charge/macros/jetchargeAN/draw_templateFitting.C
/// ----------------------------------------------------------------------- ///
// perform various fitting processes
void fit_templates( double k = 0.0, TString tune = "Detroit", TString opt1 = "Minuit", TString opt2 = "Minimize", double R = 0.4 ){
    // perform tfractionfitter fit to pythia-8 total distribution using pythia-8 templates
    cout << "\nfitting pythia-8 with pythia-8 templates\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n\n";
    
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");//, "Minimize");
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer( opt1, opt2 );
    //ROOT::Math::MinimizerOptions::SetDefaultErrorDef( 0.5 );
    //ROOT::Math::MinimizerOptions::SetDefaultStrategy( 0 );
    //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );
    ROOT::Math::MinimizerOptions::PrintDefault();

    double ymax = jc_padmax[0]; // default to k = 0.0 window
    for(int ik = 0; ik < nkaps; ik++){
        if( k == k_vals[ik] ){ ymax = jc_padmax[ik]; }
    }
    
    // 40 000 000 000 jets in pythia-8 with 20 < pT < 25 GeV
    
    //p8_undecayed_R04_k03.root

    // original binning for all kappas, only binning available for k = 0.0 due to integer addition making finer bins worthless
    TFile *fp8;//
    if( tune == "Monash" ){ fp8 = new TFile( "../../out/p8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" ); }
    else if( tune == "Detroit" ){ fp8 = new TFile( "../../out/detroit_tune_pythia8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" ); }
    // templates of the form "u_2025"

    
    const int nflavs = 3;
    TString flavs_names[nflavs+1] = {"d", "u", "g", "other"};

    // 13 = 6*quarks + 6*antiquarks + gluon
    // - nflavs = n flavor templates being used
    // - 1 = omit s as the start of "other" template
    TString other_flav[13-nflavs-1] = {"c", "b", "t", "dbar", "ubar", "sbar", "cbar", "bbar", "tbar"};


    TH2D* p8_2d = (TH2D*) fp8->Get( "qvpt" ); // 2D histogram with full pythia-8 distribution of q vs pt
    TH1D* hpy8[njetbins]; // p8_2d projected onto 1D histograms


    double temp_fracs[nflavs+1]; // array of pythia-8 flavor fractions to use to constrain fitting if needed/used
    // templates don't get normalized, so this may be redundant?
    TH2D* tmp_2D[nflavs+1]; // 2D histograms by flavor
    TH1D* htemps[njetbins][nflavs+1]; // 2D histograms projected onto 1Ds

    
    for(int i = 0; i < nflavs; i++){
        tmp_2D[i] = (TH2D*) fp8->Get( flavs_names[i] + "_qvpt" );

        for(int j = 0; j < njetbins; j++){
            if(i == 0){
                hpy8[j] = (TH1D*) p8_2d->ProjectionX( "p8_" + outFile_pt[j], p8_2d->GetYaxis()->FindBin( jetPtLo[j] ), p8_2d->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );

                cout << "pythia jets = " << hpy8[j]->Integral() << "\n\n";
                //hpy8[j]->Scale( 1.0 / 1000.0 );
            }

            htemps[j][i] = (TH1D*) tmp_2D[i]->ProjectionX( flavs_names[i] + "_" + outFile_pt[j], tmp_2D[i]->GetYaxis()->FindBin( jetPtLo[j] ), tmp_2D[i]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
            
        }
    }
    // "other" 2D histogram, to be projected onto 1D histogram
    tmp_2D[nflavs] = (TH2D*) fp8->Get( "s_qvpt" ); // start "other" template with s, then loop through c,b,t, sbar,cbar,bbar,tbar
    for(int i = 0; i < 13-nflavs-1; i++){
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
    
    
    TCanvas *c = new TCanvas("c_fits", "canvas", 800, 500);
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
        TString padup_name = Form( "pstack_%i", j );
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
        
        
        hpy8[j]->SetFillColor(kBlack);
        hpy8[j]->SetMarkerColor(kBlack);
        hpy8[j]->SetMarkerStyle(kFullSquare);
        //        hpy8[j]->SetLineWidth(2);
        
        if( j == 0 ){
            leg2->AddEntry( hpy8[j], "PYTHIA-8", "p" );
        }
        
        
        fit[j] = new TFractionFitter( hpy8[j], h_fit[j] /*, "V"*/ );
        
        for(int i = 0; i < nflavs; i++){
            fit[j]->Constrain( i, 0.0, 1.0 );
            //( 1.0 - 0.5 ) * ( frac_temps[j][i]->Integral() / frac_py8[j]->Integral() ),
            //( 1.0 + 0.5 ) * ( frac_temps[j][i]->Integral() / frac_py8[j]->Integral() ) );
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
        result_tff[j]->SetMarkerStyle(kFullSquare);
        result_tff[j]->SetStats(0);
        
        //hpy8[j]->Scale( 1000.0 );
        hpy8[j]->Scale( 1.0 / py8_tot );
        
        cout << "pythia integral = " << hpy8[j]->Integral() << "\n\n";
        
        result_tff[j]->Scale( 1.0 / result_tff[j]->Integral( 0, result_tff[j]->GetNbinsX() + 1 ) );
        
        result_tff[j]->Divide( hpy8[j] );
        
        
        hpy8[j]->SetStats(0);
        
        hpy8[j]->GetYaxis()->SetRangeUser( 0.001, ymax );
        if(k == 0.0){
            hpy8[j]->GetXaxis()->SetRangeUser( -4.5, 4.5 );
        }
        else{
            hpy8[j]->GetXaxis()->SetRangeUser( -1, 1 );
        }


        hpy8[j]->GetYaxis()->SetTitleFont(43);
        hpy8[j]->GetYaxis()->SetTitleSize(22);
        hpy8[j]->GetYaxis()->SetTitleOffset(1.3);
        hpy8[j]->GetYaxis()->SetLabelFont(43);
        hpy8[j]->GetYaxis()->SetLabelSize(18);
        hpy8[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );


        hpy8[j]->GetXaxis()->SetTitleFont(43);
        hpy8[j]->GetXaxis()->SetTitleSize(22);
        hpy8[j]->GetXaxis()->SetTitleOffset(1.1);
        hpy8[j]->GetXaxis()->SetLabelFont(43);
        hpy8[j]->GetXaxis()->SetLabelSize(18);
        hpy8[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );

        
        hpy8[j]->Draw("PSAME");

        
        for(int i = 0; i < nflavs+1; i++){
            fit_preds[j][i] = (TH1D*) fit[j]->GetMCPrediction(i);

            fit_preds[j][i]->SetFillColor( mark_colors[i] );
            fit_preds[j][i]->SetLineColor( mark_colors[i] );

            //fit_preds[j][i]->Scale( 1.0 / 1000.0 );
            fit_preds[j][i]->Scale( 1.0 / py8_tot ); // normalize by total number of jets
        }
        
        // add in same order as cms just to be consistent for now -- 11/16/22
        // other, down, up, gluon
        h_result[j]->Add(fit_preds[j][nflavs]);
        h_result[j]->Add(fit_preds[j][0]);
        h_result[j]->Add(fit_preds[j][1]);
        h_result[j]->Add(fit_preds[j][2]);

        h_result[j]->Draw("HISTSAME"); // draw histogram stack
        hpy8[j]->DrawCopy("PSAME");
        
        if( j == 0 ){
            drawText( "PYTHIA-8 " + tune, .35, .93, 20 );
            //drawText( "p_{T}^{track} > " + tr_pt + " MeV", .38, .78, 15 );
            
            leg1->AddEntry( fit_preds[j][0], "down quark", "f" );
            //leg3->AddEntry( (TObject*)0, Form( "%2.0f %%", frac ), "" );
            
            leg1->AddEntry( fit_preds[j][1], "up quark", "f" );
            //leg3->AddEntry( (TObject*)0, Form( "%2.0f %%", frac ), "" );
            
            leg1->AddEntry( fit_preds[j][2], "gluon", "f" );
            //leg3->AddEntry( (TObject*)0, Form( "%2.0f %%", frac ), "" );
            
            leg1->AddEntry( fit_preds[j][3], "other", "f" );
        }
        
        if( j == 0 ){
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .36, .87, 15 );
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
        result_tff[j]->GetYaxis()->SetTitle( "Fit / PYTHIA-8" );

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

            //double frac = frac_temps[j][i]->Integral( 0, frac_temps[j][i]->GetNbinsX() + 1 ) / py8_tot;
            double frac = 1.0;// set to parameter from the fit result
            double err = 1.0;

            fit[j]->GetResult( i, frac, err );


            cout << "\npythia-8 flavor " << flavs_names[i] << " fraction: " << htemps[j][i]->Integral() / py8_tot << "\n";
            cout << "fit result for u fraction: " << frac << "\n";
        }
        cout << "\n";
    }// j loop over jet pt bins
    c->SaveAs( "./figs/" + tune + "/fits/templateFitP8_R" + rad + "_k" + kappa + ".pdf", "RECREATE" );
}




// fit data with pythia-8 templates
void fit_data( double k = 0.0, TString tune = "Detroit", TString opt1 = "Minuit", TString opt2 = "Minimize", double R = 0.4 ){
    // perform tfractionfitter fit to pythia-8 total distribution using pythia-8 templates
    cout << "\nfitting unfolded data with pythia-8 templates\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n\n";
    
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");//, "Minimize");
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer( opt1, opt2 );
    //ROOT::Math::MinimizerOptions::SetDefaultErrorDef( 0.5 );
    //ROOT::Math::MinimizerOptions::SetDefaultStrategy( 0 );
    //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );
    ROOT::Math::MinimizerOptions::PrintDefault();

    
    double ymax = jc_padmax[0]; // default to k = 0.0 window
    for(int ik = 0; ik < nkaps; ik++){
        if( k == k_vals[ik] ){ ymax = jc_padmax[ik]; }
    }
    
    // 206, 691 jets in unfolded JP2 data with 20 < pT < 25 GeV
    
    TFile *fdat = new TFile( "../../out/unfolded_R" + rad + "_k" + kappa + ".root", "READ" );
    //unfolded_R04_k07.root
    
    //p8_undecayed_R04_k03.root
    TFile *fp8;//
    if( tune == "Monash" ){ fp8 = new TFile( "../../out/p8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" ); }
    else if( tune == "Detroit" ){ fp8 = new TFile( "../../out/detroit_tune_pythia8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" ); }
    
    //these files include pythia8_norm_ + outFile_pt[j] and the individual flavor templates
    // templates of the form "u_2025"


    
    const int nflavs = 3;
    TString flavs_names[nflavs+1] = {"d", "u", "g", "other"};

    // 13 = 6*quarks + 6*antiquarks + gluon
    // - nflavs = n flavor templates being used
    // - 1 = omit s as the start of "other" template
    TString other_flav[13-nflavs-1] = {"c", "b", "t", "dbar", "ubar", "sbar", "cbar", "bbar", "tbar"};


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
        
        cout << "\nunfolded data statistics: " << hdat[j]->Integral() << " jets\n\n";
        
        hdat[j]->Scale( 1000.0 );
    }
    
    for(int i = 0; i < nflavs; i++){
        
        tmp_2D[i] = (TH2D*) fp8->Get( flavs_names[i] + "_qvpt" );

        for(int j = 0; j < njetbins; j++){
            if(i == 0){
                hpy8[j] = (TH1D*) p8_2d->ProjectionX( "p8_" + outFile_pt[j], p8_2d->GetYaxis()->FindBin( jetPtLo[j] ), p8_2d->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
                
                //cout << "binning in pythia-8: " << hpy8[j]->GetNbinsX() << " bins in [" << hpy8[j]->GetBinLowEdge(1) << ", " << hpy8[j]->GetBinLowEdge(hpy8[j]->GetNbinsX()) << "]\n";
            }

            htemps[j][i] = (TH1D*) tmp_2D[i]->ProjectionX( flavs_names[i] + "_" + outFile_pt[j], tmp_2D[i]->GetYaxis()->FindBin( jetPtLo[j] ), tmp_2D[i]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );

            //htemps[j][i]->Scale( 1.0 / 1000.0 );
        }
    }
    
    
    // "other" 2D histogram, to be projected onto 1D histogram
    tmp_2D[nflavs] = (TH2D*) fp8->Get( "s_qvpt" ); // start "other" template with s, then loop through c,b,t, sbar,cbar,bbar,tbar
    for(int i = 0; i < 13-nflavs-1; i++){
        tmp_2D[nflavs]->Add( (TH2D*) fp8->Get( other_flav[i] + "_qvpt" ) );
    }

    // project 2D "other" (aka s template) into 1D "other" list of histogram
    for(int j = 0; j < njetbins; j++){
        htemps[j][nflavs] = (TH1D*) tmp_2D[nflavs]->ProjectionX( "other_" + outFile_pt[j], tmp_2D[nflavs]->GetYaxis()->FindBin( jetPtLo[j] ), tmp_2D[nflavs]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        //htemps[j][nflavs]->Scale( 1.0 / 1000.0 );
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
        fit[j] = new TFractionFitter( hdat[j], h_fit[j] /*, "V"*/);

        
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

        hdat[j]->GetYaxis()->SetRangeUser( 0.001, ymax );
        if(k == 0.0){
            hdat[j]->GetXaxis()->SetRangeUser( -4.5, 4.5 );
        }
        else{
            hdat[j]->GetXaxis()->SetRangeUser( -1, 1 );
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

            //fit_preds[j][i]->Scale( 1000.0 );
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

            leg1->AddEntry( fit_preds[j][0], "down quark", "f" );
            //leg3->AddEntry( (TObject*)0, Form( "%2.0f %%", frac ), "" );
            
            leg1->AddEntry( fit_preds[j][1], "up quark", "f" );
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
    c->SaveAs( "./figs/" + tune + "/fits/unfoldedData_R" + rad + "_k" + kappa + ".pdf", "RECREATE" );
}



// include fitting to each systematic source's unfolded distributions
void fit_systematics( double k = 0.0, TString tune = "Detroit", TString opt1 = "Minuit", TString opt2 = "Minimize", double R = 0.4 ){
    cout << "\ndrawing fit result for systematics, with ratio of fit result and systematic distributions\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n\n";

    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"/*, "Minimize"*/);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer( opt1, opt2 );
    //ROOT::Math::MinimizerOptions::SetDefaultErrorDef( 0.5 );
    //ROOT::Math::MinimizerOptions::SetDefaultStrategy( 0 );
    //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );
    ROOT::Math::MinimizerOptions::PrintDefault();

    double ymax = jc_padmax[0]; // default to k = 0.0 window
    for(int ik = 0; ik < nkaps; ik++){
        if( k == k_vals[ik] ){ ymax = jc_padmax[ik]; }
    }
    
    
    TFile *fdat = new TFile( "../../out/unfolded_R" + rad + "_k" + kappa + ".root", "READ" );
    //unfolded_R04_k07.root
    
    //p8_undecayed_R04_k03.root
    TFile *fp8;//
    if( tune == "Monash" ){ fp8 = new TFile( "../../out/p8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" ); }
    else if( tune == "Detroit" ){ fp8 = new TFile( "../../out/detroit_tune_pythia8_undecayed_R" + rad + "_k" + kappa + ".root", "READ" ); }
    // templates of the form "u_2025"
    

    // systematic sources naming
    const int nsources = 7; // 9 with 1D Q smearing (P8, H7)
    TString syst_name[nsources] = {"IP2", "IP6", "TS", "TU", "HC50", "DS", "GS"/*, "P8", "H7"*/};
    
    const int nflavs = 3;
    TString flavs_names[nflavs+1] = {"d", "u", "g", "other"};

    // 13 = 6*quarks + 6*antiquarks + gluon
    // - nflavs = n flavor templates being used
    // - 1 = omit s as the start of "other" template
    TString other_flav[13-nflavs-1] = {"c", "b", "t", "dbar", "ubar", "sbar", "cbar", "bbar", "tbar"};
    
    // cleaner way to fill the above is to have a list saved of the 13 flavors and fill the above with flavors that do not match any in flavs_names, but this will do for now, since settled on u,d,g, other
    
    
    TH2D* p8_2d = (TH2D*) fp8->Get( "qvpt" ); // 2D histogram with full pythia-8 distribution of q vs pt
    TH1D* hpy8[njetbins]; // p8_2d projected onto 1D histograms
    
    
    double temp_fracs[nflavs+1]; // array of pythia-8 flavor fractions to use to constrain fitting if needed/used
    // templates don't get normalized, so this may be redundant?
    TH2D* tmp_2D[nflavs+1]; // 2D histograms by flavor
    TH1D* htemps[njetbins][nflavs+1]; // 2D histograms projected onto 1Ds
    
    
    for(int i = 0; i < nflavs; i++){
    
        tmp_2D[i] = (TH2D*) fp8->Get( flavs_names[i] + "_qvpt" );

        for(int j = 0; j < njetbins; j++){
            if(i == 0){
                hpy8[j] = (TH1D*) p8_2d->ProjectionX( "p8_" + outFile_pt[j], p8_2d->GetYaxis()->FindBin( jetPtLo[j] ), p8_2d->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
            }
            
            htemps[j][i] = (TH1D*) tmp_2D[i]->ProjectionX( flavs_names[i] + "_" + outFile_pt[j], tmp_2D[i]->GetYaxis()->FindBin( jetPtLo[j] ), tmp_2D[i]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        }
    }
    
    // "other" 2D histogram, to be projected onto 1D histogram
    tmp_2D[nflavs] = (TH2D*) fp8->Get( "s_qvpt" ); // start "other" template with s, then loop through c,b,t, sbar,cbar,bbar,tbar
    for(int i = 0; i < 13-nflavs-1; i++){
        tmp_2D[nflavs]->Add( (TH2D*) fp8->Get( other_flav[i] + "_qvpt" ) );
    }
    
    // project 2D "other" (aka s template) into 1D "other" list of histogram
    for(int j = 0; j < njetbins; j++){
        htemps[j][nflavs] = (TH1D*) tmp_2D[nflavs]->ProjectionX( "other_" + outFile_pt[j], tmp_2D[nflavs]->GetYaxis()->FindBin( jetPtLo[j] ), tmp_2D[nflavs]->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
    }
    

    // loop over nsources to fit pythia to each
    for(int isys = 0; isys < nsources; isys++){
        
        cout << "\n\nsystematic: " << syst_name[isys] << "\n\n";
    
        
        
        //data histograms
        TH1D* hdat[njetbins];
        for(int j = 0; j < njetbins; j++){
            // non-normalized unfolded results (so that statistical errors don't fuck up the fitting)
            // need to figure out the naming convention used for systematics
            // syst_name[i] + "_copy_" + outFile_pt[j]
            // where syst_name[] = {"IP2", "IP6", "TS", "TU", "HC50", "DS", "GS"};//, and then something for Q smeared systematics
            hdat[j] = (TH1D*) fdat->Get( syst_name[isys] + "_copy_" + outFile_pt[j] );
            
            cout << "\nsystematic integral = " << hdat[j]->Integral() << "\n\n";
            // still normalized, so fitting will have to wait until i remake unfolded files
        }
    
        
        
        
        // for tfractionfitter
        TObjArray *h_fit[njetbins];
        TFractionFitter *fit[njetbins];
        
        TH1D* result_tff[njetbins];
        
        THStack *h_result[njetbins]; // stacking flavor histograms from fit result
        TH1D* fit_preds[njetbins][nflavs+1]; // nflavs + other; store fit[j]->GetMCPrediction(i) here
        
        
        TCanvas *c = new TCanvas("c_" + syst_name[isys], "canvas", 800, 500);
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
                leg2->AddEntry( hdat[j], syst_name[isys], "p" );
            }
            
            
            // start to perform the fit on unfolded data
            fit[j] = new TFractionFitter( hdat[j], h_fit[j] /*, "V"*/);
            
            
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
//          hdat[j]->Sumw2(0);
            
            hdat[j]->GetYaxis()->SetRangeUser( 0.001, ymax );
            if(k == 0.0){
                hdat[j]->GetXaxis()->SetRangeUser( -4.5, 4.5 );
            }
            else{
                hdat[j]->GetXaxis()->SetRangeUser( -1, 1 );
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
                // should be ~ (non-normalized) Template * fraction
                
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
                drawText( "Corrected Systematic", .35, .93, 20 );
                
                leg1->AddEntry( fit_preds[j][0], "down quark", "f" );
                //leg3->AddEntry( (TObject*)0, Form( "%2.0f %%", frac ), "" );
                
                leg1->AddEntry( fit_preds[j][1], "up quark", "f" );
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
            result_tff[j]->GetYaxis()->SetTitle( "Fit / Systematic" );
            
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
        c->SaveAs( "./figs/" + tune + "fits/systematics_" + syst_name[isys] + "_R" + rad + "_k" + kappa + ".pdf", "RECREATE" );
    } // loop over systematic sources
}



// run all above macros
void drawTemplates( double k = 0.0, TString tune = "Detroit", TString opt1 = "Minuit", TString opt2 = "Minimize", double R = 0.4 ){
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n\n";
    
    
    
    draw_templates_MonashVsDetroit( k , R , 3 );
    //fit_templates( k , tune , opt1 , opt2 , R  );
    //fit_data( k , tune, opt1 , opt2 , R  );
    //fit_systematics( k , tune, opt1 , op2 , R  );
    
    
}



