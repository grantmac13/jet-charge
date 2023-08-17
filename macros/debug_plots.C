// 8/8/2023 Grant McNamara
// copied from ~/macros/jetchargeAN/thesis/ to ~/macros/jetchargeAN/thesis/debug
// to test that jet pt, constituent pt, nef spectra as well as unfolding are consistent.
// first look at raw data, det-level comparison to check that this is still good agreement between data and pythia-6+geant
// then plot pt spectra; then the constituent pt, nef as function of jet pt that i report
// then the response, matched jets; and finally the unfolded distribution to pythia-6 part-level


#include <string>
#include <iostream>
#include "math.h"
#include "../../../Headers/plot.h"
#include "../../../Headers/JCparameters.hh"


using namespace jcAnalysis;



// draw data and pythia-6 det-level jet charge for a given kappa, jet R on TCanvas together
void plot_unfPart_jetcharge( double k = 0.0, double R = 0.4 ){
    cout << "\nPlotting unfolded data and PYTHIA on one canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    
    //TFile *fdat = new TFile( "./out/data_unfold_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    TFile *fdat = new TFile( "./out/unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    
    TFile *fp6 = new TFile( "./out/py6_unmatched_hists_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../../../out/p6_matched_R" + rad + "_k" + kappa + ".root" , "READ" );

    
    TH2D* p62 = (TH2D*) fp6->Get( "QvPt_p" );
    
    
    TH1D* unfs[njetbins]; // unfolded distributions go here
    TH1D* p6_p[njetbins]; // part-level distriubtions from pythia-6 (soon to be unmatched)
    
    
    
    TH1D* p6_r[njetbins]; // ratio of pythia-6/data
    TH1D* er_r[njetbins]; // ratio of uncertainty histogram/unfolded data
    
    TH1D* syst_errs_envs[njetbins];
    
    
    TCanvas *c = new TCanvas("c_unf", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);

    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    TLegend *leg3 = new TLegend(0.3, 0.6, .7, .75);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString padup_name = Form( "pup_%i", j );
        TPad* padup = new TPad(padup_name, padup_name, 0, 0.4, 1, 1.0);
        
        padup->SetRightMargin(0);
        padup->SetTopMargin(0);
        padup->SetBottomMargin(0);
        padup->SetLeftMargin(0);
        if( j == 0 ){padup->SetLeftMargin(0.2);}
        
        padup->Draw();
        padup->cd();

        
        
        unfs[j] = (TH1D*) fdat->Get( "unfold_nomx" + outFile_pt[j] );
        syst_errs_envs[j] = (TH1D*) fdat->Get( Form( "w_systs_%i", j+1 ) );
        
        p6_p[j] = (TH1D*) p62->ProjectionX( "p6_part_" + outFile_pt[j], p62->GetYaxis()->FindBin( jetPtLo[j] ), p62->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );

        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        unfs[j]->Scale( 1.0 / unfs[j]->Integral( 0, unfs[j]->GetNbinsX() + 1 ) );
        p6_p[j]->Scale( 1.0 / p6_p[j]->Integral( 0, p6_p[j]->GetNbinsX() + 1 ) );
        unfs[j]->Scale( 1.0 / unfs[j]->Integral( 0, unfs[j]->GetNbinsX() + 1 ) );
        
        
        unfs[j]->SetMarkerColor(kRed);
        p6_p[j]->SetLineColor(kBlue);
        
        syst_errs_envs[j]->SetFillColor(kRed - 10);
        
        
        
        unfs[j]->SetMarkerStyle(kFullStar);
//        unfs[j]->SetMarkerSize(2);
        
        p6_p[j]->SetLineStyle(kSolid);
        p6_p[j]->SetLineWidth(3);
        
        
        syst_errs_envs[j]->SetFillStyle(1001);
        
        
        // to take ratios of pythia-6(8)/unfolded data
        p6_r[j] = (TH1D*) p6_p[j]->Clone( "p6ratio_" + outFile_pt[j] );
        
        er_r[j] = (TH1D*) syst_errs_envs[j]->Clone( "syst_errs_envs_" + outFile_pt[j] );
        
        
        p6_r[j]->Divide( unfs[j] );
        er_r[j]->Divide( unfs[j] );
        
        p6_r[j]->GetYaxis()->SetRangeUser(0.001, 1.999);
        er_r[j]->GetYaxis()->SetRangeUser(0.001, 1.999);
        
        if(k == 0.0){
            p6_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
            er_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            p6_r[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            er_r[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
        if( j == 0 ){
            leg1->AddEntry( unfs[j], "Corrected STAR", "p" );
        }
        else if( j == 1 ){
            leg2->AddEntry( p6_p[j], "PYTHIA-6", "l" ); // only data with markers, all else with lines only
        }
        
        
        syst_errs_envs[j]->SetStats(0);
        syst_errs_envs[j]->GetYaxis()->SetRangeUser(0.001, 0.699);
        syst_errs_envs[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        
                
//      syst_errs_envs[j]->GetYaxis()->SetTitleSize(0.1);
        syst_errs_envs[j]->GetYaxis()->SetTitleFont(43);
        syst_errs_envs[j]->GetYaxis()->SetTitleSize(22);
        syst_errs_envs[j]->GetYaxis()->SetTitleOffset(1.2);
        syst_errs_envs[j]->GetYaxis()->SetNdivisions(505);
//        syst_errs_envs[j]->GetYaxis()->SetLabelSize(0.05);
        syst_errs_envs[j]->GetYaxis()->SetLabelFont(43);
        syst_errs_envs[j]->GetYaxis()->SetLabelSize(22);
        syst_errs_envs[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        
        
        p6_p[j]->SetStats(0);
        p6_p[j]->GetYaxis()->SetRangeUser(0.001, 0.699);
        
        if(k == 0.0){
            p6_p[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            p6_p[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
        
        p6_p[j]->GetYaxis()->SetTitleFont(43);
        p6_p[j]->GetYaxis()->SetTitleSize(22);
        p6_p[j]->GetYaxis()->SetTitleOffset(1.2);
        p6_p[j]->GetYaxis()->SetNdivisions(505);
        p6_p[j]->GetYaxis()->SetLabelFont(43);
        p6_p[j]->GetYaxis()->SetLabelSize(22);
        p6_p[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        

        
        p6_p[j]->Sumw2(0);
        
        
        syst_errs_envs[j]->Draw("E3SAME");
        p6_p[j]->Draw("CSAME");
        unfs[j]->Draw("PSAME");
        //line1->Draw();
        
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            //leg2->Draw();
            drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
            
        }
        else if( j == 1){
            leg2->Draw();
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            leg3->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( 0, 0.001, 0, 0.599 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        //line2->Draw();
        
        
        c->cd(j+1);
        TString paddown_name = Form( "paddown_%i", j );
        TPad *paddown = new TPad( paddown_name, paddown_name, 0, 0.05, 1, 0.4);
        
        paddown->SetTopMargin(0);
        paddown->SetBottomMargin(0.4);
        paddown->SetRightMargin(0);
        paddown->SetLeftMargin(0);
        if( j == 0 ){
            paddown->SetLeftMargin(0.2);
        }
        paddown->Draw();
        paddown->cd();
        
        er_r[j]->SetStats(0);
        
        er_r[j]->GetYaxis()->SetTitleFont(43);
        er_r[j]->GetYaxis()->SetTitleSize(22);
        er_r[j]->GetYaxis()->SetTitleOffset(1.2);
        er_r[j]->GetYaxis()->SetNdivisions(505);
        er_r[j]->GetYaxis()->SetLabelFont(43);
        er_r[j]->GetYaxis()->SetLabelSize(22);
        er_r[j]->GetYaxis()->SetTitle( "MC / data" );
        
        er_r[j]->GetXaxis()->SetTitleFont(43);
        er_r[j]->GetXaxis()->SetTitleSize(22);
        er_r[j]->GetXaxis()->SetTitleOffset(1.1);
        er_r[j]->GetXaxis()->SetLabelSize(0.1);
        er_r[j]->GetXaxis()->SetLabelFont(43);
        er_r[j]->GetXaxis()->SetLabelSize(22);
        
        er_r[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        er_r[j]->Draw("E3SAME");
        p6_r[j]->Draw("E1SAME");
        
        
        cout << "\n";
        cout << jetPtLo[j] << " < p_{T}^{jet} < " << jetPtHi[j] << " jet charge means\n";
        cout << "data: " << unfs[j]->GetMean() << "\n";
        cout << "pythia-6: " << p6_p[j]->GetMean() << "\n";
        cout << "\n\n";
        
    }

    c->SaveAs( "./plots/unfAndPart_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}


void plot_dataDet_jetcharge( double k = 0.0, double R = 0.4 ){
    cout << "\nPlotting raw data and PYTHIA+GEANT on one canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    // file with both raw and unfolded data distributions for all jet pT bins
    //TFile *fdat = new TFile( "./out/data_hists_ppjp2_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    TFile *fdat = new TFile( "./out/unfolded_R" + rad + "_k" + kappa + ".root" , "READ");

    // file with pythia-6 distribution
    // unmatched pythia-6 --- only has part-level distribution...
    TFile *fp6 = new TFile( "./out/py6_unmatched_hists_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");

    
    
    TH1D* data[njetbins]; // unfolded distributions go here
    TH1D* p6_d[njetbins]; // part-level distriubtions
    
    
    TH2D* dat = (TH2D*) fdat->Get( "QvPt_d" );
    TH2D* p62 = (TH2D*) fp6->Get( "QvPt_d" );


    // for ratios
    TH1D* det_r[njetbins]; // ratio pythia-6/data
    
    
    //fdat->Close();
    //fp6->Close();
    
    
    TCanvas *c = new TCanvas("c_dat", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString pad_name = Form( "p_%i", j );
        TPad* p1 = new TPad(pad_name, pad_name, 0, 0.4, 1, 1.0);
        
        p1->SetRightMargin(0);
        p1->SetTopMargin(0);
        p1->SetBottomMargin(0);
        p1->SetLeftMargin(0);
        if( j == 0 ){p1->SetLeftMargin(0.2);}
        
        p1->Draw();
        p1->cd();

        
        
        data[j] = (TH1D*) dat->ProjectionX( "data_raw_" + outFile_pt[j], dat->GetYaxis()->FindBin( jetPtLo[j] ), dat->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        p6_d[j] = (TH1D*) p62->ProjectionX( "p6det_" + outFile_pt[j], p62->GetYaxis()->FindBin( jetPtLo[j] ), p62->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        data[j]->Scale( 1.0 / data[j]->Integral( 0, data[j]->GetNbinsX() + 1 ) );
        p6_d[j]->Scale( 1.0 / p6_d[j]->Integral( 0, p6_d[j]->GetNbinsX() + 1 ) );

        
        data[j]->SetMarkerColor(kBlack);
        p6_d[j]->SetMarkerColor(kBlue);
        
        
        data[j]->SetMarkerStyle(kFullStar);
        data[j]->SetMarkerSize(2);
        
        p6_d[j]->SetMarkerStyle(kOpenCircle);
        //p6_d[j]->SetLineStyle(kSolid);
        p6_d[j]->SetMarkerSize(2);
        //p6_d[j]->SetLineWidth(3);
        
        
        // to take ratios of pythia-6(8)/unfolded data
        det_r[j] = (TH1D*) p6_d[j]->Clone( "p6_detlevel_ratio_" + outFile_pt[j] );
        //p8_r[j] = (TH1D*) p8_p[j]->Clone( "p8ratio_" + outFile_pt[j] );
        
        det_r[j]->Divide( data[j] );
        //p8_r[j]->Divide( unfs[j] );
        
        det_r[j]->GetYaxis()->SetRangeUser(0.001, 1.999);
        
        if(k == 0.0){
            det_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            det_r[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }//p8_r[j]->GetYaxis()->SetRangeUser(0,2);
        
        
        if( j == 0 ){
            leg1->AddEntry( data[j], "STAR (uncorrected)", "p");

        }
        else if( j == 1 ){
            leg2->AddEntry( p6_d[j], "PYTHIA-6+GEANT", "p"); // only data with markers, all else with lines only
        }
        
        data[j]->SetStats(0);
        data[j]->GetYaxis()->SetRangeUser(0.001, 0.699);
        
        if(k == 0.0){
            data[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            data[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        data[j]->GetYaxis()->SetTitleSize(0.1);
        data[j]->GetYaxis()->SetTitleFont(43);
        data[j]->GetYaxis()->SetTitleSize(22);
        data[j]->GetYaxis()->SetTitleOffset(1.2);
        data[j]->GetYaxis()->SetNdivisions(505);
//        data[j]->GetYaxis()->SetLabelSize(0.05);
        data[j]->GetYaxis()->SetLabelFont(43);
        data[j]->GetYaxis()->SetLabelSize(22);
        data[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        
        p6_d[j]->SetStats(0);
        p6_d[j]->GetYaxis()->SetRangeUser(0.001, 0.699);
        if(k == 0.0){
            p6_d[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            p6_d[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        p6_d[j]->GetYaxis()->SetTitleSize(0.1);
        p6_d[j]->GetYaxis()->SetTitleFont(43);
        p6_d[j]->GetYaxis()->SetTitleSize(22);
        p6_d[j]->GetYaxis()->SetTitleOffset(1.2);
        p6_d[j]->GetYaxis()->SetNdivisions(505);
//        p6_d[j]->GetYaxis()->SetLabelSize(0.05);
        p6_d[j]->GetYaxis()->SetLabelFont(43);
        p6_d[j]->GetYaxis()->SetLabelSize(22);
        p6_d[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        
        
        p6_d[j]->Sumw2(0);
        
        p6_d[j]->Draw("PSAME");
        data[j]->Draw("PSAME");
        
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            //leg2->Draw();
            drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
        }
        else if( j == 1){
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            leg2->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( 0, 0.001, 0, 0.599 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        //line2->Draw();
        
        
        c->cd( j+1 );
        
        TString pdown_name = Form( "pdown_%i", j );
        TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0.05, 1, 0.4);
        
        p2->SetRightMargin(0);
        p2->SetTopMargin(0);
        p2->SetBottomMargin(0.4);
        p2->SetLeftMargin(0);
        if( j == 0 ){p2->SetLeftMargin(0.2);}
        
        p2->Draw();
        p2->cd();
        
        
        det_r[j]->SetStats(0);
        
        det_r[j]->SetMarkerSize(1);
        
//        det_r[j]->GetYaxis()->SetTitleSize(0.2);
        det_r[j]->GetYaxis()->SetTitleFont(43);
        det_r[j]->GetYaxis()->SetTitleSize(22);
        det_r[j]->GetYaxis()->SetTitleOffset(1.2);
        det_r[j]->GetYaxis()->SetNdivisions(505);
//        det_r[j]->GetYaxis()->SetLabelSize(0.1);
        det_r[j]->GetYaxis()->SetLabelFont(43);
        det_r[j]->GetYaxis()->SetLabelSize(22);
        det_r[j]->GetYaxis()->SetTitle( "MC / data" );
        
        //det_r[j]->GetXaxis()->SetTitleSize(0.2);
        det_r[j]->GetXaxis()->SetTitleFont(43);
        det_r[j]->GetXaxis()->SetTitleSize(22);
        det_r[j]->GetXaxis()->SetTitleOffset(1.1);
        //det_r[j]->GetXaxis()->SetNdivisions(505);
        //det_r[j]->GetXaxis()->SetLabelSize(0.1);
        det_r[j]->GetXaxis()->SetLabelFont(43);
        det_r[j]->GetXaxis()->SetLabelSize(22);
        det_r[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        //er_r[j]->Draw("E3SAME");
        det_r[j]->Draw("E1SAME");
        
        
        TLine *line = new TLine( -4.5, 1, 4.5, 1 );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        
        line->Draw();
        
        
        cout << "\n";
        cout << jetPtLo[j] << " < p_{T}^{jet} < " << jetPtHi[j] << " jet charge means\n";
        cout << "data: " << data[j]->GetMean() << "\n";
        cout << "pythia-6: " << p6_d[j]->GetMean() << "\n";
        cout << "\n\n";
        
    }

    c->SaveAs( "./plots/dataAndDet_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}



// attempt at recreating det-level prelim figure
// testing whether the file from 10/10, 9/8 2022 is the correct one
void plot_prelim_det( double k = 0.0, double R = 0.4 ){
    cout << "\nRECREATING PRELIMINARY FIGURE?\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "../../../../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // file from ~09/08/2022 (pre-DNP)
    

    // file with pythia-6 distribution
    // unmatched pythia-6 --- only has part-level distribution...
    TFile *fp6 = new TFile( "./out/unmatchedp6_101022_R" + rad + "_k" + kappa + ".root" , "READ");
    // from 10/10/2022 (pre-DNP)
    
    // can i reproduce DNP det-level? -- maybe can figure out what i have changed since then?
    
    
    
    TH1D* data[njetbins]; // unfolded distributions go here
    TH1D* p6_d[njetbins]; // part-level distriubtions
    
    
    TH2D* dat = (TH2D*) fdat->Get( "QvPt_d" );
    TH2D* p62 = (TH2D*) fp6->Get( "QvPt_d" );


    // for ratios
    TH1D* det_r[njetbins]; // ratio pythia-6/data
    
    
    //fdat->Close();
    //fp6->Close();
    
    
    TCanvas *c = new TCanvas("c_prelim_test", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString pad_name = Form( "p_%i", j );
        TPad* p1 = new TPad(pad_name, pad_name, 0, 0.4, 1, 1.0);
        
        p1->SetRightMargin(0);
        p1->SetTopMargin(0);
        p1->SetBottomMargin(0);
        p1->SetLeftMargin(0);
        if( j == 0 ){p1->SetLeftMargin(0.2);}
        
        p1->Draw();
        p1->cd();

        
        
        data[j] = (TH1D*) dat->ProjectionX( "data_raw_" + outFile_pt[j], dat->GetYaxis()->FindBin( jetPtLo[j] ), dat->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        p6_d[j] = (TH1D*) p62->ProjectionX( "p6det_" + outFile_pt[j], p62->GetYaxis()->FindBin( jetPtLo[j] ), p62->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        data[j]->Scale( 1.0 / data[j]->Integral( 0, data[j]->GetNbinsX() + 1 ) );
        p6_d[j]->Scale( 1.0 / p6_d[j]->Integral( 0, p6_d[j]->GetNbinsX() + 1 ) );

        
        data[j]->SetMarkerColor(kBlack);
        p6_d[j]->SetMarkerColor(kBlue);
        
        
        data[j]->SetMarkerStyle(kFullStar);
        data[j]->SetMarkerSize(2);
        
        p6_d[j]->SetMarkerStyle(kOpenCircle);
        //p6_d[j]->SetLineStyle(kSolid);
        p6_d[j]->SetMarkerSize(2);
        //p6_d[j]->SetLineWidth(3);
        
        
        // to take ratios of pythia-6(8)/unfolded data
        det_r[j] = (TH1D*) p6_d[j]->Clone( "p6_detlevel_ratio_" + outFile_pt[j] );
        //p8_r[j] = (TH1D*) p8_p[j]->Clone( "p8ratio_" + outFile_pt[j] );
        
        det_r[j]->Divide( data[j] );
        //p8_r[j]->Divide( unfs[j] );
        
        det_r[j]->GetYaxis()->SetRangeUser(0.001, 1.999);
        
        if(k == 0.0){
            det_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            det_r[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }//p8_r[j]->GetYaxis()->SetRangeUser(0,2);
        
        
        if( j == 0 ){
            leg1->AddEntry( data[j], "STAR (uncorrected)", "p");

        }
        else if( j == 1 ){
            leg2->AddEntry( p6_d[j], "PYTHIA-6+GEANT", "p"); // only data with markers, all else with lines only
        }
        
        data[j]->SetStats(0);
        data[j]->GetYaxis()->SetRangeUser(-0.099, 0.599);
        
        if(k == 0.0){
            data[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            data[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        data[j]->GetYaxis()->SetTitleSize(0.1);
        data[j]->GetYaxis()->SetTitleFont(43);
        data[j]->GetYaxis()->SetTitleSize(22);
        data[j]->GetYaxis()->SetTitleOffset(1.2);
        data[j]->GetYaxis()->SetNdivisions(505);
//        data[j]->GetYaxis()->SetLabelSize(0.05);
        data[j]->GetYaxis()->SetLabelFont(43);
        data[j]->GetYaxis()->SetLabelSize(22);
        data[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        
        p6_d[j]->SetStats(0);
        p6_d[j]->GetYaxis()->SetRangeUser(-0.099, 0.599);
        if(k == 0.0){
            p6_d[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            p6_d[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        p6_d[j]->GetYaxis()->SetTitleSize(0.1);
        p6_d[j]->GetYaxis()->SetTitleFont(43);
        p6_d[j]->GetYaxis()->SetTitleSize(22);
        p6_d[j]->GetYaxis()->SetTitleOffset(1.2);
        p6_d[j]->GetYaxis()->SetNdivisions(505);
//        p6_d[j]->GetYaxis()->SetLabelSize(0.05);
        p6_d[j]->GetYaxis()->SetLabelFont(43);
        p6_d[j]->GetYaxis()->SetLabelSize(22);
        p6_d[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        
        
        p6_d[j]->Sumw2(0);
        
        p6_d[j]->Draw("PSAME");
        data[j]->Draw("PSAME");
        
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            //leg2->Draw();
            drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
        }
        else if( j == 1){
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            leg2->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( 0, 0.001, 0, 0.599 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        //line2->Draw();
        
        
        c->cd( j+1 );
        
        TString pdown_name = Form( "pdown_%i", j );
        TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0.05, 1, 0.4);
        
        p2->SetRightMargin(0);
        p2->SetTopMargin(0);
        p2->SetBottomMargin(0.4);
        p2->SetLeftMargin(0);
        if( j == 0 ){p2->SetLeftMargin(0.2);}
        
        p2->Draw();
        p2->cd();
        
        
        det_r[j]->SetStats(0);
        
        det_r[j]->SetMarkerSize(1);
        
//        det_r[j]->GetYaxis()->SetTitleSize(0.2);
        det_r[j]->GetYaxis()->SetTitleFont(43);
        det_r[j]->GetYaxis()->SetTitleSize(22);
        det_r[j]->GetYaxis()->SetTitleOffset(1.2);
        det_r[j]->GetYaxis()->SetNdivisions(505);
//        det_r[j]->GetYaxis()->SetLabelSize(0.1);
        det_r[j]->GetYaxis()->SetLabelFont(43);
        det_r[j]->GetYaxis()->SetLabelSize(22);
        det_r[j]->GetYaxis()->SetTitle( "MC / data" );
        
        //det_r[j]->GetXaxis()->SetTitleSize(0.2);
        det_r[j]->GetXaxis()->SetTitleFont(43);
        det_r[j]->GetXaxis()->SetTitleSize(22);
        det_r[j]->GetXaxis()->SetTitleOffset(1.1);
        //det_r[j]->GetXaxis()->SetNdivisions(505);
        //det_r[j]->GetXaxis()->SetLabelSize(0.1);
        det_r[j]->GetXaxis()->SetLabelFont(43);
        det_r[j]->GetXaxis()->SetLabelSize(22);
        det_r[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        //er_r[j]->Draw("E3SAME");
        det_r[j]->Draw("E1SAME");
        
        
        TLine *line = new TLine( -4.5, 1, 4.5, 1 );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        
        line->Draw();
        
        
        cout << "\n";
        cout << jetPtLo[j] << " < p_{T}^{jet} < " << jetPtHi[j] << " jet charge means\n";
        cout << "data: " << data[j]->GetMean() << "\n";
        cout << "pythia-6: " << p6_d[j]->GetMean() << "\n";
        cout << "\n\n";
        
    }

    c->SaveAs( "./plots/TESTprelim_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}



void plot_conspt( double k = 0.0, double R = 0.4 ){
    cout << "\nPlotting constituent pt for raw data; part-, det-level on one canvas\n";
    
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "./out/data_hists_ppjp2_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    //TFile *fdat = new TFile( "./out/data_unfold_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    
    // file with pythia-6 distribution
    // unmatched pythia-6 --- only has part-level distribution...
    TFile *fp6 = new TFile( "./out/py6_unmatched_hists_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");
    
    
    //TH2D* unf_qpt = (TH2D*) funf->Get( "unfold_nom" );
    
    TH1D* dat_jets = (TH1D*) fdat->Get( "jet_pt_hist" );
    TH1D* p6_det_jets = (TH1D*) fp6->Get( "jet_pt_hist" );
    TH1D* p6_part_jets = (TH1D*) fp6->Get( "jet_pt_hist_part" );
    
    
    
    TH1D* data[njetbins]; // unfolded distributions go here
    TH1D* p6_d[njetbins]; // part-level distriubtions
    TH1D* p6_p[njetbins]; // part-level distriubtions
    
    
    TH2D* dat = (TH2D*) fdat->Get( "det_cons_vs_jetpt" );
    TH2D* p62_d = (TH2D*) fp6->Get( "det_cons_vs_jetpt" );
    TH2D* p62_p = (TH2D*) fp6->Get( "part_cons_vs_jetpt" );
    
    
    // for ratios
    TH1D* det_r[njetbins]; // ratio pythia-6/data
    TH1D* part_r[njetbins]; // ratio pythia-6/data
    
    
    
    
    TCanvas *c = new TCanvas("c_conspt", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString pad_name = Form( "p_%i", j );
        TPad* p1 = new TPad(pad_name, pad_name, 0, 0.4, 1, 1.0);
        
        p1->SetRightMargin(0);
        p1->SetTopMargin(0);
        p1->SetBottomMargin(0);
        p1->SetLeftMargin(0);
        if( j == 0 ){p1->SetLeftMargin(0.2);}
        
        p1->Draw();
        p1->cd();
        
        p1->SetLogy();
        
        data[j] = (TH1D*) dat->ProjectionX( "data_raw_" + outFile_pt[j], dat->GetYaxis()->FindBin( jetPtLo[j] ), dat->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        p6_d[j] = (TH1D*) p62_d->ProjectionX( "p6_det_" + outFile_pt[j], p62_d->GetYaxis()->FindBin( jetPtLo[j] ), p62_d->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        p6_p[j] = (TH1D*) p62_p->ProjectionX( "p6_part_" + outFile_pt[j], p62_p->GetYaxis()->FindBin( jetPtLo[j] ), p62_p->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        data[j]->Scale( 1.0 / data[j]->Integral() );
        p6_d[j]->Scale( 1.0 / p6_d[j]->Integral() );
        p6_p[j]->Scale( 1.0 / p6_p[j]->Integral() );

        /*
        data[j]->Scale( 1.0 / dat_jets->Integral( dat_jets->GetXaxis()->FindBin( jetPtLo[j] ), dat_jets->GetXaxis()->FindBin( jetPtHi[j] ) - 1 ) );
        p6_d[j]->Scale( 1.0 / p6_det_jets->Integral( p6_det_jets->GetXaxis()->FindBin( jetPtLo[j] ), p6_det_jets->GetXaxis()->FindBin( jetPtHi[j] ) - 1 ) );
        p6_p[j]->Scale( 1.0 / p6_part_jets->Integral( p6_part_jets->GetXaxis()->FindBin( jetPtLo[j] ), p6_part_jets->GetXaxis()->FindBin( jetPtHi[j] ) - 1 ) );
        */
        
        data[j]->SetMarkerColor(kBlack);
        p6_d[j]->SetMarkerColor(kBlue);
        p6_p[j]->SetMarkerColor(kGreen);
        
        
        data[j]->SetMarkerStyle(kFullStar);
        data[j]->SetMarkerSize(2);
        
        p6_d[j]->SetMarkerStyle(kOpenCircle);
        p6_d[j]->SetMarkerSize(2);
        
        p6_p[j]->SetMarkerStyle(kOpenSquare);
        p6_p[j]->SetMarkerSize(2);
        
        
        
        // to take ratios of pythia-6(8)/unfolded data
        det_r[j] = (TH1D*) p6_d[j]->Clone( "p6_detlevel_ratio_" + outFile_pt[j] );
        part_r[j] = (TH1D*) p6_p[j]->Clone( "p6_partlevel_ratio_" + outFile_pt[j] );
        
        det_r[j]->Divide( data[j] );
        part_r[j]->Divide( data[j] );
        
        det_r[j]->GetYaxis()->SetRangeUser(0.501, 1.499);
        
        det_r[j]->GetXaxis()->SetRangeUser(0.0, 1.0);
        
        
        if( j == 0 ){
            leg1->AddEntry( data[j], "STAR (uncorrected)", "p");
            
        }
        else if( j == 1 ){
            //leg2->AddEntry( p6_p[j], "PYTHIA-6", "p");
            leg2->AddEntry( p6_d[j], "PYTHIA-6+GEANT", "p"); // only data with markers, all else with lines only
        }
        
        data[j]->SetStats(0);
//        data[j]->GetYaxis()->SetRangeUser(-0.099, 0.599);
        
        //data[j]->GetXaxis()->SetRangeUser(0.0, 10.0);
        
        //        data[j]->GetYaxis()->SetTitleSize(0.1);
        data[j]->GetYaxis()->SetTitleFont(43);
        data[j]->GetYaxis()->SetTitleSize(22);
        data[j]->GetYaxis()->SetTitleOffset(1.2);
        data[j]->GetYaxis()->SetNdivisions(505);
        //        data[j]->GetYaxis()->SetLabelSize(0.05);
        data[j]->GetYaxis()->SetLabelFont(43);
        data[j]->GetYaxis()->SetLabelSize(22);
        data[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dp_{T}^{jet} [1/(GeV/c)]" );
        
        
        p6_d[j]->SetStats(0); p6_p[j]->SetStats(0);
        //p6_d[j]->GetYaxis()->SetRangeUser(0.0, 0.6);
        p6_d[j]->GetXaxis()->SetRangeUser(0.0, 1.0);
        
        
        //        p6_d[j]->GetYaxis()->SetTitleSize(0.1);
        p6_d[j]->GetYaxis()->SetTitleFont(43);
        p6_d[j]->GetYaxis()->SetTitleSize(18);
        p6_d[j]->GetYaxis()->SetTitleOffset(1.2);
        p6_d[j]->GetYaxis()->SetNdivisions(505);
        //        p6_d[j]->GetYaxis()->SetLabelSize(0.05);
        p6_d[j]->GetYaxis()->SetLabelFont(43);
        p6_d[j]->GetYaxis()->SetLabelSize(22);
        p6_d[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dp_{T}^{jet} [1/(GeV/c))]" );
        
        
        
        //p6_d[j]->Sumw2(0);
        
        p6_d[j]->Draw("PSAME");
        data[j]->Draw("PSAME");
        //p6_p[j]->Draw("PSAME");
        
        
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            //leg2->Draw();
            drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
        }
        else if( j == 1){
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            leg2->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        
        c->SetLogy();
        
        
        TLine *line2 = new TLine( 0.2, 0.0, 0.2, 1.0 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        line2->Draw();
        
        
        c->cd( j+1 );
        
        TString pdown_name = Form( "pdown_%i", j );
        TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0.05, 1, 0.4);
        
        p2->SetRightMargin(0);
        p2->SetTopMargin(0);
        p2->SetBottomMargin(0.4);
        p2->SetLeftMargin(0);
        if( j == 0 ){p2->SetLeftMargin(0.2);}
        
        p2->Draw();
        p2->cd();
        
        
        det_r[j]->SetStats(0);
        
        det_r[j]->SetMarkerSize(1);
        
        //        det_r[j]->GetYaxis()->SetTitleSize(0.2);
        det_r[j]->GetYaxis()->SetTitleFont(43);
        det_r[j]->GetYaxis()->SetTitleSize(22);
        det_r[j]->GetYaxis()->SetTitleOffset(1.2);
        det_r[j]->GetYaxis()->SetNdivisions(505);
        //        det_r[j]->GetYaxis()->SetLabelSize(0.1);
        det_r[j]->GetYaxis()->SetLabelFont(43);
        det_r[j]->GetYaxis()->SetLabelSize(22);
        det_r[j]->GetYaxis()->SetTitle( "MC / data" );
        
        
        //det_r[j]->GetXaxis()->SetTitleSize(0.2);
        det_r[j]->GetXaxis()->SetTitleFont(43);
        det_r[j]->GetXaxis()->SetTitleSize(22);
        det_r[j]->GetXaxis()->SetTitleOffset(1.1);
        //det_r[j]->GetXaxis()->SetNdivisions(505);
        //det_r[j]->GetXaxis()->SetLabelSize(0.1);
        det_r[j]->GetXaxis()->SetLabelFont(43);
        det_r[j]->GetXaxis()->SetLabelSize(22);
        det_r[j]->GetXaxis()->SetTitle( "p_{T}^{cons}" );
        
        
        
        //er_r[j]->Draw("E3SAME");
        det_r[j]->Draw("E1SAME");
        //part_r[j]->Draw("E1SAME");
        
        
        TLine *line = new TLine( 0.0, 1, 2.0, 1 );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        
        line->Draw();
        
    }
    
    c->SaveAs( "./plots/constituentPtSpectrum_R" + rad + "_k" + kappa + "_" + ".pdf", "RECREATE");
}



void plot_nef( double k = 0.0, double R = 0.4 ){
    cout << "\nPlotting NEF comparison on a canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "./out/data_hists_ppjp2_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    //TFile *fdat = new TFile( "./out/data_unfold_correct_R" + rad + "_k" + kappa + ".root" , "READ");

    // file with pythia-6 distribution
    // unmatched pythia-6 --- only has part-level distribution...
    TFile *fp6 = new TFile( "./out/py6_unmatched_hists_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");

    
    
    TH1D* data[njetbins]; // unfolded distributions go here
    TH1D* p6_d[njetbins]; // part-level distriubtions
    
    
    TH2D* dat = (TH2D*) fdat->Get( "det_nef_vs_jetpt" );
    TH2D* p62 = (TH2D*) fp6->Get( "det_nef_vs_jetpt" );


    // for ratios
    TH1D* det_r[njetbins]; // ratio pythia-6/data
    
    
    
    TCanvas *c = new TCanvas("c_nef", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString pad_name = Form( "p_%i", j );
        TPad* p1 = new TPad(pad_name, pad_name, 0, 0.4, 1, 1.0);
        
        p1->SetRightMargin(0);
        p1->SetTopMargin(0);
        p1->SetBottomMargin(0);
        p1->SetLeftMargin(0);
        if( j == 0 ){p1->SetLeftMargin(0.2);}
        
        p1->Draw();
        p1->cd();

        
        
        data[j] = (TH1D*) dat->ProjectionX( "data_raw_" + outFile_pt[j], dat->GetYaxis()->FindBin( jetPtLo[j] ), dat->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        p6_d[j] = (TH1D*) p62->ProjectionX( "p6det_" + outFile_pt[j], p62->GetYaxis()->FindBin( jetPtLo[j] ), p62->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        data[j]->Scale( 1.0 / data[j]->Integral( ) );
        p6_d[j]->Scale( 1.0 / p6_d[j]->Integral( ) );

        
        data[j]->SetMarkerColor(kBlack);
        p6_d[j]->SetMarkerColor(kBlue);
        
        
        data[j]->SetMarkerStyle(kFullStar);
        data[j]->SetMarkerSize(2);
        
        p6_d[j]->SetMarkerStyle(kOpenCircle);
        //p6_d[j]->SetLineStyle(kSolid);
        p6_d[j]->SetMarkerSize(2);
        //p6_d[j]->SetLineWidth(3);
        
        
        // to take ratios of pythia-6(8)/unfolded data
        det_r[j] = (TH1D*) p6_d[j]->Clone( "p6_detlevel_ratio_" + outFile_pt[j] );
        //p8_r[j] = (TH1D*) p8_p[j]->Clone( "p8ratio_" + outFile_pt[j] );
        
        det_r[j]->Divide( data[j] );
        //p8_r[j]->Divide( unfs[j] );
        
        det_r[j]->GetYaxis()->SetRangeUser(0.001, 1.999);
        
        data[j]->GetXaxis()->SetRangeUser(0.0, 1.0);
        data[j]->GetYaxis()->SetRangeUser(0.0, 0.1);

        
        if( j == 0 ){
            leg1->AddEntry( data[j], "STAR (uncorrected)", "p");

        }
        else if( j == 1 ){
            leg2->AddEntry( p6_d[j], "PYTHIA-6+GEANT", "p"); // only data with markers, all else with lines only
        }
        
        data[j]->SetStats(0);
        


//        data[j]->GetYaxis()->SetTitleSize(0.1);
        data[j]->GetYaxis()->SetTitleFont(43);
        data[j]->GetYaxis()->SetTitleSize(16);
        data[j]->GetYaxis()->SetTitleOffset(1.2);
//        data[j]->GetYaxis()->SetNdivisions(505);
        data[j]->GetYaxis()->SetLabelSize(0.05);
        data[j]->GetYaxis()->SetLabelFont(43);
       data[j]->GetYaxis()->SetLabelSize(16);
        data[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );

        
        p6_d[j]->SetStats(0);
        
/*
//        p6_d[j]->GetYaxis()->SetTitleSize(0.1);
        p6_d[j]->GetYaxis()->SetTitleFont(43);
        p6_d[j]->GetYaxis()->SetTitleSize(22);
        p6_d[j]->GetYaxis()->SetTitleOffset(1.2);
        p6_d[j]->GetYaxis()->SetNdivisions(505);
//        p6_d[j]->GetYaxis()->SetLabelSize(0.05);
        p6_d[j]->GetYaxis()->SetLabelFont(43);
        p6_d[j]->GetYaxis()->SetLabelSize(22);
        p6_d[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
*/
        
        
        data[j]->Draw("PSAME");
        p6_d[j]->Draw("PSAME");
        
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            //leg2->Draw();
            drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
        }
        else if( j == 1){
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            leg2->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( 0, 0.001, 0, 0.599 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        //line2->Draw();
        
        
        c->cd( j+1 );
        
        TString pdown_name = Form( "pdown_%i", j );
        TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0.05, 1, 0.4);
        
        p2->SetRightMargin(0);
        p2->SetTopMargin(0);
        p2->SetBottomMargin(0.4);
        p2->SetLeftMargin(0);
        if( j == 0 ){p2->SetLeftMargin(0.2);}
        
        p2->Draw();
        p2->cd();
        
        
        det_r[j]->SetStats(0);
        
        det_r[j]->SetMarkerSize(1);
        
//        det_r[j]->GetYaxis()->SetTitleSize(0.2);
        det_r[j]->GetYaxis()->SetTitleFont(43);
        det_r[j]->GetYaxis()->SetTitleSize(22);
        det_r[j]->GetYaxis()->SetTitleOffset(1.2);
        det_r[j]->GetYaxis()->SetNdivisions(505);
//        det_r[j]->GetYaxis()->SetLabelSize(0.1);
        det_r[j]->GetYaxis()->SetLabelFont(43);
        det_r[j]->GetYaxis()->SetLabelSize(22);
        det_r[j]->GetYaxis()->SetTitle( "MC / data" );
        
        
        //det_r[j]->GetXaxis()->SetTitleSize(0.2);
        det_r[j]->GetXaxis()->SetTitleFont(43);
        det_r[j]->GetXaxis()->SetTitleSize(18);
//        det_r[j]->GetXaxis()->SetTitleOffset(1.1);
        //det_r[j]->GetXaxis()->SetNdivisions(505);
        //det_r[j]->GetXaxis()->SetLabelSize(0.1);
        det_r[j]->GetXaxis()->SetLabelFont(43);
        det_r[j]->GetXaxis()->SetLabelSize(18);
        det_r[j]->GetXaxis()->SetTitle( "Neutral Energy Fraction" );
        
        
        //er_r[j]->Draw("E3SAME");
        det_r[j]->Draw("E1SAME");
        
        
    }

    c->SaveAs( "./plots/nef_debug_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}



void plot_jetspectrum( double k = 0.0, double R = 0.4 ){
    cout << "\nPlotting jet spectra for all four on one canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "./out/data_hists_ppjp2_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    TFile *funf = new TFile( "./out/data_unfold_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    // file with pythia-6 distribution
    // unmatched pythia-6 --- only has part-level distribution...
    TFile *fp6 = new TFile( "./out/py6_unmatched_hists_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../../../out/p6_matched_R" + rad + "_k" + kappa + ".root" , "READ" );

    
    
    TH2D* unf_qpt = (TH2D*) funf->Get( "unfold_nom" );
    
    TH1D* data = (TH1D*) fdat->Get( "jet_pt_hist" );
    TH1D* p6_d = (TH1D*) fp6->Get( "jet_pt_hist" );
    TH1D* p6_p = (TH1D*) fp6->Get( "jet_pt_hist_part" );
    TH1D* unfs = (TH1D*) unf_qpt->ProjectionY( "unfold_pt" );

    
    cout << "\nraw data pt spectrum... integral = " << data->Integral() << "\n\n";
    cout << "\nraw data pt spectrum... mean = " << data->GetMean() << "\n\n";
    
    cout << "\nunfolded pt spectrum... integral = " << unfs->Integral() << "\n\n";
    cout << "\nunfolded pt spectrum... mean = " << unfs->GetMean() << "\n\n";
    
    
    cout << "\npythia+geant pt spectrum... integral = " << p6_d->Integral() << "\n\n";
    cout << "\npythia+geant pt spectrum... mean = " << p6_d->GetMean() << "\n\n";
    
    cout << "\npythia pt spectrum... integral = " << p6_p->Integral() << "\n\n";
    cout << "\npythia pt spectrum... mean = " << p6_p->GetMean() << "\n\n";

    
    TCanvas *c = new TCanvas("c_jetspectra", "canvas", 800, 800);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.55, 0.75, .87, .87);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.03);

    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    TLegend *leg3 = new TLegend(0.3, 0.6, .7, .75);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.05);
    
    
    // normalize all to a probability distribution (integral = 1) including over/underflow bins
    data->Scale( 1.0 / data->Integral(  ) );
    p6_d->Scale( 1.0 / p6_d->Integral(  ) );
    p6_p->Scale( 1.0 / p6_p->Integral(  ) );
    unfs->Scale( 1.0 / unfs->Integral(  ) );
    
    
    unfs->SetMarkerColor(kRed);
    p6_p->SetMarkerColor(kGreen);
    data->SetMarkerColor(kBlack);
    p6_d->SetMarkerColor(kBlue);
    
    
    
    
    data->SetMarkerStyle(kFullStar);
//    data->SetMarkerSize(2);
    
    p6_d->SetMarkerStyle(kOpenCircle);
    
    unfs->SetMarkerStyle(kFullStar);
//    unfs->SetMarkerSize(2);
    
    p6_p->SetLineStyle(kOpenSquare);
//    p6_p->SetLineWidth(3);
    
    
    
    leg1->AddEntry( data, "Raw STAR", "p" );
    //leg1->AddEntry( unfs, "Corrected STAR", "p" );
    //leg2->AddEntry( p6_p, "PYTHIA-6", "p" );
    leg1->AddEntry( p6_d, "PYTHIA-6+GEANT", "p" );
    
    
    
    p6_d->SetStats(0);
    
    p6_d->GetYaxis()->SetTitleFont(43);
    p6_d->GetYaxis()->SetTitleSize(22);
//    p6_d->GetYaxis()->SetTitleOffset(1.2);
    p6_d->GetYaxis()->SetNdivisions(505);
    p6_d->GetYaxis()->SetLabelFont(43);
    p6_d->GetYaxis()->SetLabelSize(22);
    p6_d->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dp_{T}^{jet} [1/(GeV/c)]" );
    
    
    //p6_d->Sumw2(0);
    //p6_p->Sumw2(0);
    
    p6_d->Draw("PSAME");
    data->Draw("PSAME");
    //p6_p->Draw("PSAME");
    //unfs->Draw("PSAME");
    //line1->Draw();
    
    //leg1->Draw();
    //leg2->Draw();
    
    leg1->Draw();
    c->SetLogy();
    
    drawText( "p+p #sqrt{s} = 200 GeV", .25, .36, 18 );
    drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .25, .3, 18 );
    drawText( "p_{T}^{cons} > 200 MeV", .25, .24, 18 );
    
    /*
    leg1->Draw();
    //leg2->Draw();
    drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
    drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
    
    leg2->Draw();
    drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
    drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
    leg3->Draw();
    drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
    drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
    */
    

    c->SaveAs( "./plots/jetptspectra_debug_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}


void compare_to_preDNP_unfold( double k = 0.0, double R = 0.4 )
{
    cout << "\nComparing current unfolded distribution to as old as i can find\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    
    TFile *fdat = new TFile( "./out/unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    
    TFile *fp6 = new TFile( "../../../../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../../../out/p6_matched_R" + rad + "_k" + kappa + ".root" , "READ" );

    
    //TH2D* p62 = (TH2D*) fp6->Get( "QvPt_p" );
    
    
    TH1D* unfs[njetbins]; // unfolded distributions go here
    TH1D* p6_p[njetbins]; // part-level distriubtions from pythia-6 (soon to be unmatched)
    
    
    TH2D* new_tmp = (TH2D*) fdat->Get( "unfold_nom" );
    TH2D* old_tmp = (TH2D*) fp6->Get( "unfold_nom" );

    
    TH1D* new_q[njetbins];
    TH1D* old_q[njetbins];
    
    
    TH1D* p6_r[njetbins]; // ratio of pythia-6/data
    TH1D* er_r[njetbins]; // ratio of uncertainty histogram/unfolded data
    
    TH1D* syst_errs_envs[njetbins];
    
    
    TCanvas *c = new TCanvas("c_unf_dnp", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);

    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    TLegend *leg3 = new TLegend(0.3, 0.6, .7, .75);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString padup_name = Form( "pup_%i", j );
        TPad* padup = new TPad(padup_name, padup_name, 0, 0.4, 1, 1.0);
        
        padup->SetRightMargin(0);
        padup->SetTopMargin(0);
        padup->SetBottomMargin(0);
        padup->SetLeftMargin(0);
        if( j == 0 ){padup->SetLeftMargin(0.2);}
        
        padup->Draw();
        padup->cd();

        
        
        unfs[j] = (TH1D*) fdat->Get( "unfold_nomx" + outFile_pt[j] );
        syst_errs_envs[j] = (TH1D*) fdat->Get( Form( "w_systs_%i", j+1 ) );
        
        p6_p[j] = (TH1D*) fp6->Get( "unfold_nomx" + outFile_pt[j] );

        
        new_q[j] = (TH1D*) new_tmp->ProjectionX( "new_" + outFile_pt[j] , new_tmp->GetYaxis()->FindBin( jetPtLo[j] ) , new_tmp->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        old_q[j] = (TH1D*) old_tmp->ProjectionX( "old_" + outFile_pt[j] , old_tmp->GetYaxis()->FindBin( jetPtLo[j] ) , old_tmp->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        
        cout << "jet pt: [" << jetPtLo[j] << ", " << jetPtHi[j] << "]\n";
        cout << "\tcurrent integral = " << new_q[j]->Integral( 0, new_q[j]->GetNbinsX() + 1 ) << "\n";
        cout << "\tDNP integral = " << old_q[j]->Integral( 0, old_q[j]->GetNbinsX() + 1 ) << "\n";

        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        unfs[j]->Scale( 1.0 / unfs[j]->Integral( 0, unfs[j]->GetNbinsX() + 1 ) );
        p6_p[j]->Scale( 1.0 / p6_p[j]->Integral( 0, p6_p[j]->GetNbinsX() + 1 ) );
        
        
        
        unfs[j]->SetMarkerColor(kRed);
        p6_p[j]->SetLineColor(kGreen);
        
        syst_errs_envs[j]->SetFillColor(kRed - 10);
        
        
        
        unfs[j]->SetMarkerStyle(kFullStar);
//        unfs[j]->SetMarkerSize(2);
        
        p6_p[j]->SetMarkerStyle(kOpenCircle);
        p6_p[j]->SetLineWidth(3);
        
        
        syst_errs_envs[j]->SetFillStyle(1001);
        
        
        // to take ratios of pythia-6(8)/unfolded data
        p6_r[j] = (TH1D*) p6_p[j]->Clone( "p6ratio_" + outFile_pt[j] );
        
        er_r[j] = (TH1D*) syst_errs_envs[j]->Clone( "syst_errs_envs_" + outFile_pt[j] );
        
        
        p6_r[j]->Divide( unfs[j] );
        er_r[j]->Divide( unfs[j] );
        
        p6_r[j]->GetYaxis()->SetRangeUser(0.801, 1.199);
        er_r[j]->GetYaxis()->SetRangeUser(0.801, 1.199);
        
        if(k == 0.0){
            p6_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
            er_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            p6_r[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            er_r[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
        if( j == 0 ){
            leg1->AddEntry( unfs[j], "Current (8/23)", "p" );
        }
        else if( j == 1 ){
            leg2->AddEntry( p6_p[j], "From DNP (10/22)", "l" ); // only data with markers, all else with lines only
        }
        
        
        syst_errs_envs[j]->SetStats(0);
        syst_errs_envs[j]->GetYaxis()->SetRangeUser(0.001, 0.699);
        syst_errs_envs[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        
        unfs[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        
                
//      syst_errs_envs[j]->GetYaxis()->SetTitleSize(0.1);
        syst_errs_envs[j]->GetYaxis()->SetTitleFont(43);
        syst_errs_envs[j]->GetYaxis()->SetTitleSize(22);
        syst_errs_envs[j]->GetYaxis()->SetTitleOffset(1.2);
        syst_errs_envs[j]->GetYaxis()->SetNdivisions(505);
//        syst_errs_envs[j]->GetYaxis()->SetLabelSize(0.05);
        syst_errs_envs[j]->GetYaxis()->SetLabelFont(43);
        syst_errs_envs[j]->GetYaxis()->SetLabelSize(22);
        syst_errs_envs[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        
        
        p6_p[j]->SetStats(0);
        p6_p[j]->GetYaxis()->SetRangeUser(0.001, 0.699);
        
        if(k == 0.0){
            p6_p[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            p6_p[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
        
        p6_p[j]->GetYaxis()->SetTitleFont(43);
        p6_p[j]->GetYaxis()->SetTitleSize(22);
        p6_p[j]->GetYaxis()->SetTitleOffset(1.2);
        p6_p[j]->GetYaxis()->SetNdivisions(505);
        p6_p[j]->GetYaxis()->SetLabelFont(43);
        p6_p[j]->GetYaxis()->SetLabelSize(22);
        p6_p[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        

        
        //p6_p[j]->Sumw2(0);
        
        
        syst_errs_envs[j]->Draw("E3SAME");
        p6_p[j]->Draw("PSAME");
        unfs[j]->Draw("PSAME");
        //line1->Draw();
        
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            //leg2->Draw();
            drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
            
        }
        else if( j == 1){
            leg2->Draw();
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            leg3->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( 0, 0.001, 0, 0.599 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        //line2->Draw();
        
        
        c->cd(j+1);
        TString paddown_name = Form( "paddown_%i", j );
        TPad *paddown = new TPad( paddown_name, paddown_name, 0, 0.05, 1, 0.4);
        
        paddown->SetTopMargin(0);
        paddown->SetBottomMargin(0.4);
        paddown->SetRightMargin(0);
        paddown->SetLeftMargin(0);
        if( j == 0 ){
            paddown->SetLeftMargin(0.2);
        }
        paddown->Draw();
        paddown->cd();
        
        er_r[j]->SetStats(0);
        
        er_r[j]->GetYaxis()->SetTitleFont(43);
        er_r[j]->GetYaxis()->SetTitleSize(22);
        er_r[j]->GetYaxis()->SetTitleOffset(1.2);
        er_r[j]->GetYaxis()->SetNdivisions(505);
        er_r[j]->GetYaxis()->SetLabelFont(43);
        er_r[j]->GetYaxis()->SetLabelSize(22);
        er_r[j]->GetYaxis()->SetTitle( "MC / data" );
        
        er_r[j]->GetXaxis()->SetTitleFont(43);
        er_r[j]->GetXaxis()->SetTitleSize(22);
        er_r[j]->GetXaxis()->SetTitleOffset(1.1);
        er_r[j]->GetXaxis()->SetLabelSize(0.1);
        er_r[j]->GetXaxis()->SetLabelFont(43);
        er_r[j]->GetXaxis()->SetLabelSize(22);
        
        er_r[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        er_r[j]->Draw("E3SAME");
        p6_r[j]->Draw("PSAME");
        
        
        cout << "\n";
        cout << jetPtLo[j] << " < p_{T}^{jet} < " << jetPtHi[j] << " jet charge means\n";
        cout << "Current: " << unfs[j]->GetMean() << "\n";
        cout << "DNP version: " << p6_p[j]->GetMean() << "\n";
        cout << "\n\n";
        
    }

    c->SaveAs( "./plots/unfoldCompToDNP_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
    
}



void compare_to_preDNP_raw( double k = 0.0, double R = 0.4 ){
    cout << "\nComparing raw data from current (8/23) to DNP (10/22)\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "./out/unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    //TFile *fdat = new TFile( "./out/data_unfold_correct_R" + rad + "_k" + kappa + ".root" , "READ");

    // file with pythia-6 distribution
    // unmatched pythia-6 --- only has part-level distribution...
    TFile *fp6 = new TFile( "../../../../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");

    
    
    TH1D* data[njetbins]; // unfolded distributions go here
    TH1D* p6_d[njetbins]; // part-level distriubtions
    
    
    TH2D* dat = (TH2D*) fdat->Get( "QvPt_d" );
    TH2D* p62 = (TH2D*) fp6->Get( "QvPt_d" );


    // for ratios
    TH1D* det_r[njetbins]; // ratio pythia-6/data
    
    
    //fdat->Close();
    //fp6->Close();
    
    
    TCanvas *c = new TCanvas("c_raw_dnp", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString pad_name = Form( "p_%i", j );
        TPad* p1 = new TPad(pad_name, pad_name, 0, 0.4, 1, 1.0);
        
        p1->SetRightMargin(0);
        p1->SetTopMargin(0);
        p1->SetBottomMargin(0);
        p1->SetLeftMargin(0);
        if( j == 0 ){p1->SetLeftMargin(0.2);}
        
        p1->Draw();
        p1->cd();

        
        
        data[j] = (TH1D*) dat->ProjectionX( "data_raw_" + outFile_pt[j], dat->GetYaxis()->FindBin( jetPtLo[j] ), dat->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        p6_d[j] = (TH1D*) p62->ProjectionX( "p6det_" + outFile_pt[j], p62->GetYaxis()->FindBin( jetPtLo[j] ), p62->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        //data[j]->Scale( 1.0 / data[j]->Integral( 0, data[j]->GetNbinsX() + 1 ) );
        //p6_d[j]->Scale( 1.0 / p6_d[j]->Integral( 0, p6_d[j]->GetNbinsX() + 1 ) );

        
        data[j]->SetMarkerColor(kBlack);
        p6_d[j]->SetMarkerColor(kBlue);
        
        
        data[j]->SetMarkerStyle(kFullStar);
        data[j]->SetMarkerSize(2);
        
        p6_d[j]->SetMarkerStyle(kOpenCircle);
        //p6_d[j]->SetLineStyle(kSolid);
        p6_d[j]->SetMarkerSize(2);
        //p6_d[j]->SetLineWidth(3);
        
        
        // to take ratios of pythia-6(8)/unfolded data
        det_r[j] = (TH1D*) p6_d[j]->Clone( "p6_detlevel_ratio_" + outFile_pt[j] );
        //p8_r[j] = (TH1D*) p8_p[j]->Clone( "p8ratio_" + outFile_pt[j] );
        
        det_r[j]->Divide( data[j] );
        //p8_r[j]->Divide( unfs[j] );
        
        det_r[j]->GetYaxis()->SetRangeUser(0.801, 1.199);
        
        if(k == 0.0){
            det_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            det_r[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }//p8_r[j]->GetYaxis()->SetRangeUser(0,2);
        
        
        if( j == 0 ){
            leg1->AddEntry( data[j], "Current", "p");

        }
        else if( j == 1 ){
            leg2->AddEntry( p6_d[j], "DNP", "p"); // only data with markers, all else with lines only
        }
        
        double ymax = 0.699*data[j]->Integral( 0, data[j]->GetNbinsX() + 1 );
        
        data[j]->SetStats(0);
        data[j]->GetYaxis()->SetRangeUser(0.001, ymax);
        
        if(k == 0.0){
            data[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            data[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        data[j]->GetYaxis()->SetTitleSize(0.1);
        data[j]->GetYaxis()->SetTitleFont(43);
        data[j]->GetYaxis()->SetTitleSize(22);
        data[j]->GetYaxis()->SetTitleOffset(1.2);
        data[j]->GetYaxis()->SetNdivisions(505);
//        data[j]->GetYaxis()->SetLabelSize(0.05);
        data[j]->GetYaxis()->SetLabelFont(43);
        data[j]->GetYaxis()->SetLabelSize(22);
        data[j]->GetYaxis()->SetTitle( "dN_{jet}/dQ_{jet} [1/e]" );
        
        
        p6_d[j]->SetStats(0);
        p6_d[j]->GetYaxis()->SetRangeUser(0.001, ymax);
        if(k == 0.0){
            p6_d[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            p6_d[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        p6_d[j]->GetYaxis()->SetTitleSize(0.1);
        p6_d[j]->GetYaxis()->SetTitleFont(43);
        p6_d[j]->GetYaxis()->SetTitleSize(22);
        p6_d[j]->GetYaxis()->SetTitleOffset(1.2);
        p6_d[j]->GetYaxis()->SetNdivisions(505);
//        p6_d[j]->GetYaxis()->SetLabelSize(0.05);
        p6_d[j]->GetYaxis()->SetLabelFont(43);
        p6_d[j]->GetYaxis()->SetLabelSize(22);
        p6_d[j]->GetYaxis()->SetTitle( "dN_{jet}/dQ_{jet} [1/e]" );
        
        
        
        //p6_d[j]->Sumw2(0);
        
        p6_d[j]->Draw("PSAME");
        data[j]->Draw("PSAME");
        
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            //leg2->Draw();
            drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
        }
        else if( j == 1){
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            leg2->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( 0, 0.001, 0, 0.599 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        //line2->Draw();
        
        
        c->cd( j+1 );
        
        TString pdown_name = Form( "pdown_%i", j );
        TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0.05, 1, 0.4);
        
        p2->SetRightMargin(0);
        p2->SetTopMargin(0);
        p2->SetBottomMargin(0.4);
        p2->SetLeftMargin(0);
        if( j == 0 ){p2->SetLeftMargin(0.2);}
        
        p2->Draw();
        p2->cd();
        
        
        det_r[j]->SetStats(0);
        
        det_r[j]->SetMarkerSize(1);
        
//        det_r[j]->GetYaxis()->SetTitleSize(0.2);
        det_r[j]->GetYaxis()->SetTitleFont(43);
        det_r[j]->GetYaxis()->SetTitleSize(22);
        det_r[j]->GetYaxis()->SetTitleOffset(1.2);
        det_r[j]->GetYaxis()->SetNdivisions(505);
//        det_r[j]->GetYaxis()->SetLabelSize(0.1);
        det_r[j]->GetYaxis()->SetLabelFont(43);
        det_r[j]->GetYaxis()->SetLabelSize(22);
        det_r[j]->GetYaxis()->SetTitle( "DNP / new" );
        
        //det_r[j]->GetXaxis()->SetTitleSize(0.2);
        det_r[j]->GetXaxis()->SetTitleFont(43);
        det_r[j]->GetXaxis()->SetTitleSize(22);
        det_r[j]->GetXaxis()->SetTitleOffset(1.1);
        //det_r[j]->GetXaxis()->SetNdivisions(505);
        //det_r[j]->GetXaxis()->SetLabelSize(0.1);
        det_r[j]->GetXaxis()->SetLabelFont(43);
        det_r[j]->GetXaxis()->SetLabelSize(22);
        det_r[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        //er_r[j]->Draw("E3SAME");
        det_r[j]->Draw("E1SAME");
        
        
        TLine *line = new TLine( -4.5, 1, 4.5, 1 );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        
        line->Draw();
        
        
        cout << "\n";
        cout << jetPtLo[j] << " < p_{T}^{jet} < " << jetPtHi[j] << " jet charge means\n";
        cout << "Current: " << data[j]->GetMean() << "\n";
        cout << "DNP: " << p6_d[j]->GetMean() << "\n";
        cout << "\n\n";
        
    }

    c->SaveAs( "./plots/compareDNP_raw_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}


// plot DNP pythia/pythia+geant and current to see how much has changed--potentially leave non-normalized to see statistics change too
void compare_to_preDNP_pythia( double k = 0.0, double R = 0.4 ){
    cout << "\nComparing pythia particle level from current (8/23) to DNP (10/22)\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    // file current pythia-6 distributions
    TFile *fdat = new TFile( "./out/py6_unmatched_hists_R" + rad + "_k" + kappa + ".root" , "READ");

    // old ("DNP") pythia
    //TFile *fp6 = new TFile( "../../../../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");
    TFile *fp6 = new TFile( "./out/unmatchedp6_101022_R" + rad + "_k" + kappa + ".root" , "READ");
    
    // matched pythia-6
    //TFile *fp6 = new TFile( "../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");

    
    
    TH1D* data[njetbins]; // unfolded distributions go here
    TH1D* p6_d[njetbins]; // part-level distriubtions
    
    
    TH2D* dat = (TH2D*) fdat->Get( "QvPt_p" );
    TH2D* p62 = (TH2D*) fp6->Get( "QvPt_p" );


    // for ratios
    TH1D* det_r[njetbins]; // ratio pythia-6/data
    
    
    //fdat->Close();
    //fp6->Close();
    
    
    TCanvas *c = new TCanvas("c_py6_dnp", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString pad_name = Form( "p_%i", j );
        TPad* p1 = new TPad(pad_name, pad_name, 0, 0.4, 1, 1.0);
        
        p1->SetRightMargin(0);
        p1->SetTopMargin(0);
        p1->SetBottomMargin(0);
        p1->SetLeftMargin(0);
        if( j == 0 ){p1->SetLeftMargin(0.2);}
        
        p1->Draw();
        p1->cd();

        
        
        data[j] = (TH1D*) dat->ProjectionX( "data_raw_" + outFile_pt[j], dat->GetYaxis()->FindBin( jetPtLo[j] ), dat->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        p6_d[j] = (TH1D*) p62->ProjectionX( "p6det_" + outFile_pt[j], p62->GetYaxis()->FindBin( jetPtLo[j] ), p62->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        //data[j]->Scale( 1.0 / data[j]->Integral( 0, data[j]->GetNbinsX() + 1 ) );
        //p6_d[j]->Scale( 1.0 / p6_d[j]->Integral( 0, p6_d[j]->GetNbinsX() + 1 ) );

        
        data[j]->SetMarkerColor(kBlack);
        p6_d[j]->SetMarkerColor(kBlue);
        
        
        data[j]->SetMarkerStyle(kFullStar);
        data[j]->SetMarkerSize(2);
        
        p6_d[j]->SetMarkerStyle(kOpenCircle);
        //p6_d[j]->SetLineStyle(kSolid);
        p6_d[j]->SetMarkerSize(2);
        //p6_d[j]->SetLineWidth(3);
        
        
        // to take ratios of pythia-6(8)/unfolded data
        det_r[j] = (TH1D*) p6_d[j]->Clone( "p6_detlevel_ratio_" + outFile_pt[j] );
        //p8_r[j] = (TH1D*) p8_p[j]->Clone( "p8ratio_" + outFile_pt[j] );
        
        det_r[j]->Divide( data[j] );
        //p8_r[j]->Divide( unfs[j] );
        
        double ymax = 0.699 * data[j]->Integral( 0, data[j]->GetNbinsX() + 1 );

        det_r[j]->GetYaxis()->SetRangeUser(0.801, 1.199);
        
        if(k == 0.0){
            det_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            det_r[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }//p8_r[j]->GetYaxis()->SetRangeUser(0,2);
        
        
        if( j == 0 ){
            leg1->AddEntry( data[j], "Current", "p");

        }
        else if( j == 1 ){
            leg2->AddEntry( p6_d[j], "DNP", "p"); // only data with markers, all else with lines only
        }
        
        data[j]->SetStats(0);
        data[j]->GetYaxis()->SetRangeUser(0.001, ymax);
        
        if(k == 0.0){
            data[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            data[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        data[j]->GetYaxis()->SetTitleSize(0.1);
        data[j]->GetYaxis()->SetTitleFont(43);
        data[j]->GetYaxis()->SetTitleSize(22);
        data[j]->GetYaxis()->SetTitleOffset(1.2);
        data[j]->GetYaxis()->SetNdivisions(505);
//        data[j]->GetYaxis()->SetLabelSize(0.05);
        data[j]->GetYaxis()->SetLabelFont(43);
        data[j]->GetYaxis()->SetLabelSize(22);
        data[j]->GetYaxis()->SetTitle( "dN_{jet}/dQ_{jet} [1/e]" );
        
        
        p6_d[j]->SetStats(0);
        p6_d[j]->GetYaxis()->SetRangeUser(0.001, ymax);
        if(k == 0.0){
            p6_d[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            p6_d[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        p6_d[j]->GetYaxis()->SetTitleSize(0.1);
        p6_d[j]->GetYaxis()->SetTitleFont(43);
        p6_d[j]->GetYaxis()->SetTitleSize(22);
        p6_d[j]->GetYaxis()->SetTitleOffset(1.2);
        p6_d[j]->GetYaxis()->SetNdivisions(505);
//        p6_d[j]->GetYaxis()->SetLabelSize(0.05);
        p6_d[j]->GetYaxis()->SetLabelFont(43);
        p6_d[j]->GetYaxis()->SetLabelSize(22);
        p6_d[j]->GetYaxis()->SetTitle( "dN_{jet}/dQ_{jet} [1/e]" );
        
        
        
        //p6_d[j]->Sumw2(0);
        
        p6_d[j]->Draw("PSAME");
        data[j]->Draw("PSAME");
        
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            //leg2->Draw();
            drawText( " PYTHIA-6 p+p"/* #sqrt{s} = 200 GeV"*/, .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
        }
        else if( j == 1){
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            leg2->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( 0, 0.001, 0, ymax );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        //line2->Draw();
        
        
        c->cd( j+1 );
        
        TString pdown_name = Form( "pdown_%i", j );
        TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0.05, 1, 0.4);
        
        p2->SetRightMargin(0);
        p2->SetTopMargin(0);
        p2->SetBottomMargin(0.4);
        p2->SetLeftMargin(0);
        if( j == 0 ){p2->SetLeftMargin(0.2);}
        
        p2->Draw();
        p2->cd();
        
        
        det_r[j]->SetStats(0);
        
        det_r[j]->SetMarkerSize(1);
        
//        det_r[j]->GetYaxis()->SetTitleSize(0.2);
        det_r[j]->GetYaxis()->SetTitleFont(43);
        det_r[j]->GetYaxis()->SetTitleSize(22);
        det_r[j]->GetYaxis()->SetTitleOffset(1.2);
        det_r[j]->GetYaxis()->SetNdivisions(505);
//        det_r[j]->GetYaxis()->SetLabelSize(0.1);
        det_r[j]->GetYaxis()->SetLabelFont(43);
        det_r[j]->GetYaxis()->SetLabelSize(22);
        det_r[j]->GetYaxis()->SetTitle( "DNP / new" );
        
        //det_r[j]->GetXaxis()->SetTitleSize(0.2);
        det_r[j]->GetXaxis()->SetTitleFont(43);
        det_r[j]->GetXaxis()->SetTitleSize(22);
        det_r[j]->GetXaxis()->SetTitleOffset(1.1);
        //det_r[j]->GetXaxis()->SetNdivisions(505);
        //det_r[j]->GetXaxis()->SetLabelSize(0.1);
        det_r[j]->GetXaxis()->SetLabelFont(43);
        det_r[j]->GetXaxis()->SetLabelSize(22);
        det_r[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        //er_r[j]->Draw("E3SAME");
        det_r[j]->Draw("E1SAME");
        
        
        TLine *line = new TLine( -4.5, 1, 4.5, 1 );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        
        line->Draw();
        
        
        cout << "\n";
        cout << jetPtLo[j] << " < p_{T}^{jet} < " << jetPtHi[j] << " jet charge means\n";
        cout << "Current: " << data[j]->GetMean() << "\n";
        cout << "DNP: " << p6_d[j]->GetMean() << "\n";
        cout << "\n\n";
        
    }

    c->SaveAs( "./plots/compareDNP_py6_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}


void compare_to_preDNP_geant( double k = 0.0, double R = 0.4 ){
    cout << "\nComparing pythia+geant detector level from current (8/23) to DNP (10/22)\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    // file current pythia-6 distributions
    TFile *fdat = new TFile( "./out/py6_unmatched_hists_R" + rad + "_k" + kappa + ".root" , "READ");

    // old ("DNP") pythia
    //TFile *fp6 = new TFile( "../../../../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");
    TFile *fp6 = new TFile( "./out/unmatchedp6_101022_R" + rad + "_k" + kappa + ".root" , "READ");
    
    // matched pythia-6
    //TFile *fp6 = new TFile( "../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");

    
    
    TH1D* data[njetbins]; // unfolded distributions go here
    TH1D* p6_d[njetbins]; // part-level distriubtions
    
    
    TH2D* dat = (TH2D*) fdat->Get( "QvPt_d" );
    TH2D* p62 = (TH2D*) fp6->Get( "QvPt_d" );


    // for ratios
    TH1D* det_r[njetbins]; // ratio pythia-6/data
    
    
    //fdat->Close();
    //fp6->Close();
    
    
    TCanvas *c = new TCanvas("c_gea_dnp", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString pad_name = Form( "p_%i", j );
        TPad* p1 = new TPad(pad_name, pad_name, 0, 0.4, 1, 1.0);
        
        p1->SetRightMargin(0);
        p1->SetTopMargin(0);
        p1->SetBottomMargin(0);
        p1->SetLeftMargin(0);
        if( j == 0 ){p1->SetLeftMargin(0.2);}
        
        p1->Draw();
        p1->cd();

        
        
        data[j] = (TH1D*) dat->ProjectionX( "data_raw_" + outFile_pt[j], dat->GetYaxis()->FindBin( jetPtLo[j] ), dat->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        p6_d[j] = (TH1D*) p62->ProjectionX( "p6det_" + outFile_pt[j], p62->GetYaxis()->FindBin( jetPtLo[j] ), p62->GetYaxis()->FindBin( jetPtHi[j] ) - 1 );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        //data[j]->Scale( 1.0 / data[j]->Integral( 0, data[j]->GetNbinsX() + 1 ) );
        //p6_d[j]->Scale( 1.0 / p6_d[j]->Integral( 0, p6_d[j]->GetNbinsX() + 1 ) );

        
        data[j]->SetMarkerColor(kBlack);
        p6_d[j]->SetMarkerColor(kBlue);
        
        
        data[j]->SetMarkerStyle(kFullStar);
        data[j]->SetMarkerSize(2);
        
        p6_d[j]->SetMarkerStyle(kOpenCircle);
        //p6_d[j]->SetLineStyle(kSolid);
        p6_d[j]->SetMarkerSize(2);
        //p6_d[j]->SetLineWidth(3);
        
        
        // to take ratios of pythia-6(8)/unfolded data
        det_r[j] = (TH1D*) p6_d[j]->Clone( "p6_detlevel_ratio_" + outFile_pt[j] );
        //p8_r[j] = (TH1D*) p8_p[j]->Clone( "p8ratio_" + outFile_pt[j] );
        
        det_r[j]->Divide( data[j] );
        //p8_r[j]->Divide( unfs[j] );
        
        double ymax = 0.699 * data[j]->Integral( 0, data[j]->GetNbinsX() + 1 );

        
        det_r[j]->GetYaxis()->SetRangeUser(0.801, 1.199);
        
        if(k == 0.0){
            det_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            det_r[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }//p8_r[j]->GetYaxis()->SetRangeUser(0,2);
        
        
        if( j == 0 ){
            leg1->AddEntry( data[j], "Current", "p");

        }
        else if( j == 1 ){
            leg2->AddEntry( p6_d[j], "DNP", "p"); // only data with markers, all else with lines only
        }
        
        data[j]->SetStats(0);
        data[j]->GetYaxis()->SetRangeUser(0.001, ymax);
        
        if(k == 0.0){
            data[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            data[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        data[j]->GetYaxis()->SetTitleSize(0.1);
        data[j]->GetYaxis()->SetTitleFont(43);
        data[j]->GetYaxis()->SetTitleSize(22);
        data[j]->GetYaxis()->SetTitleOffset(1.2);
        data[j]->GetYaxis()->SetNdivisions(505);
//        data[j]->GetYaxis()->SetLabelSize(0.05);
        data[j]->GetYaxis()->SetLabelFont(43);
        data[j]->GetYaxis()->SetLabelSize(22);
        data[j]->GetYaxis()->SetTitle( "dN_{jet}/dQ_{jet} [1/e]" );
        
        
        p6_d[j]->SetStats(0);
        p6_d[j]->GetYaxis()->SetRangeUser(0.001, ymax);
        if(k == 0.0){
            p6_d[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        }
        else{
            p6_d[j]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        }
        
//        p6_d[j]->GetYaxis()->SetTitleSize(0.1);
        p6_d[j]->GetYaxis()->SetTitleFont(43);
        p6_d[j]->GetYaxis()->SetTitleSize(22);
        p6_d[j]->GetYaxis()->SetTitleOffset(1.2);
        p6_d[j]->GetYaxis()->SetNdivisions(505);
//        p6_d[j]->GetYaxis()->SetLabelSize(0.05);
        p6_d[j]->GetYaxis()->SetLabelFont(43);
        p6_d[j]->GetYaxis()->SetLabelSize(22);
        p6_d[j]->GetYaxis()->SetTitle( "dN_{jet}/dQ_{jet} [1/e]" );
        
        
        
        //p6_d[j]->Sumw2(0);
        
        p6_d[j]->Draw("PSAME");
        data[j]->Draw("PSAME");
        
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            //leg2->Draw();
            drawText( " PYTHIA-6+GEANT"/* p+p #sqrt{s} = 200 GeV"*/, .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
        }
        else if( j == 1){
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            leg2->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( 0, 0.001, 0, 0.599 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        //line2->Draw();
        
        
        c->cd( j+1 );
        
        TString pdown_name = Form( "pdown_%i", j );
        TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0.05, 1, 0.4);
        
        p2->SetRightMargin(0);
        p2->SetTopMargin(0);
        p2->SetBottomMargin(0.4);
        p2->SetLeftMargin(0);
        if( j == 0 ){p2->SetLeftMargin(0.2);}
        
        p2->Draw();
        p2->cd();
        
        
        det_r[j]->SetStats(0);
        
        det_r[j]->SetMarkerSize(1);
        
//        det_r[j]->GetYaxis()->SetTitleSize(0.2);
        det_r[j]->GetYaxis()->SetTitleFont(43);
        det_r[j]->GetYaxis()->SetTitleSize(22);
        det_r[j]->GetYaxis()->SetTitleOffset(1.2);
        det_r[j]->GetYaxis()->SetNdivisions(505);
//        det_r[j]->GetYaxis()->SetLabelSize(0.1);
        det_r[j]->GetYaxis()->SetLabelFont(43);
        det_r[j]->GetYaxis()->SetLabelSize(22);
        det_r[j]->GetYaxis()->SetTitle( "DNP / new" );
        
        //det_r[j]->GetXaxis()->SetTitleSize(0.2);
        det_r[j]->GetXaxis()->SetTitleFont(43);
        det_r[j]->GetXaxis()->SetTitleSize(22);
        det_r[j]->GetXaxis()->SetTitleOffset(1.1);
        //det_r[j]->GetXaxis()->SetNdivisions(505);
        //det_r[j]->GetXaxis()->SetLabelSize(0.1);
        det_r[j]->GetXaxis()->SetLabelFont(43);
        det_r[j]->GetXaxis()->SetLabelSize(22);
        det_r[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        //er_r[j]->Draw("E3SAME");
        det_r[j]->Draw("E1SAME");
        
        
        TLine *line = new TLine( -4.5, 1, 4.5, 1 );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        
        line->Draw();
        
        
        cout << "\n";
        cout << jetPtLo[j] << " < p_{T}^{jet} < " << jetPtHi[j] << " jet charge means\n";
        cout << "Current: " << data[j]->GetMean() << "\n";
        cout << "DNP: " << p6_d[j]->GetMean() << "\n";
        cout << "\n\n";
        
    }

    c->SaveAs( "./plots/compareDNP_py6_geant_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}



void compare_jet_pt_counts( double k = 0.0, double R = 0.4 )
{
    cout << "\nComparing raw data jet pt counts from current (8/23) to DNP (10/22)\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    /*
    TFile *fdat = new TFile( "./out/data_hists_ppjp2_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    //TFile *fdat = new TFile( "./out/data_unfold_correct_R" + rad + "_k" + kappa + ".root" , "READ");

    // file with pythia-6 distribution
    // unmatched pythia-6 --- only has part-level distribution...
    TFile *fp6 = new TFile( "../../../../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");
    */
    
    
    
    // file current pythia-6 distributions
    TFile *fdat = new TFile( "./out/py6_unmatched_hists_R" + rad + "_k" + kappa + ".root" , "READ");

    // old ("DNP") pythia
    //TFile *fp6 = new TFile( "../../../../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");
    TFile *fp6 = new TFile( "./out/unmatchedp6_101022_R" + rad + "_k" + kappa + ".root" , "READ");
    
    
    TH2D* dat = (TH2D*) fdat->Get( "QvPt_d" );
    TH2D* p62 = (TH2D*) fp6->Get( "QvPt_d" );
    
    
    TH1D* data = (TH1D*) dat->ProjectionY( "data_pt" ); // unfolded distributions go here
    TH1D* p6_d = (TH1D*) p62->ProjectionY( "py6_pt" ); // part-level distriubtions
    
    
    
    
    TCanvas *c = new TCanvas("c_jetcounts", "canvas", 800, 800);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.55, 0.75, .87, .87);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.03);

    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    TLegend *leg3 = new TLegend(0.3, 0.6, .7, .75);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.05);
    

    
    
    data->SetMarkerColor(kBlack);
    p6_d->SetMarkerColor(kBlue);
    
    
    
    
    data->SetMarkerStyle(kFullStar);
//    data->SetMarkerSize(2);
    
    p6_d->SetMarkerStyle(kOpenCircle);
    

    
    
    leg1->AddEntry( data, "Current", "p" );
    leg1->AddEntry( p6_d, "DNP", "p" );
    
    
    
    data->SetStats(0);
    
    data->GetYaxis()->SetTitleFont(43);
    data->GetYaxis()->SetTitleSize(22);
//    data->GetYaxis()->SetTitleOffset(1.2);
    data->GetYaxis()->SetNdivisions(505);
    data->GetYaxis()->SetLabelFont(43);
    data->GetYaxis()->SetLabelSize(22);
    data->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dp_{T}^{jet} [1/(GeV/c)]" );
    
    
    //p6_d->Sumw2(0);
    //p6_p->Sumw2(0);
    
    data->Draw("PSAME");
    p6_d->Draw("PSAME");
    //line1->Draw();
    
    //leg1->Draw();
    //leg2->Draw();
    
    leg1->Draw();
    c->SetLogy();
    
    drawText( "p+p #sqrt{s} = 200 GeV", .25, .36, 18 );
    drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .25, .3, 18 );
    drawText( "p_{T}^{cons} > 200 MeV", .25, .24, 18 );
    drawText( "PYTHIA-6+GEANT", .25, .18, 18 );
    
    
    cout << "integral of jet spectrum:\n";
    cout << "20-25 GeV\n";
    cout << "\tCurrent = " << data->Integral( data->GetXaxis()->FindBin( 20 ), data->GetXaxis()->FindBin( 25 ) - 1 ) << "\n";
    cout << "\tDNP = " << p6_d->Integral( p6_d->GetXaxis()->FindBin( 20 ), p6_d->GetXaxis()->FindBin( 25 ) - 1 ) << "\n";

    cout << "25-30 GeV\n";
    cout << "\tCurrent = " << data->Integral( data->GetXaxis()->FindBin( 25 ), data->GetXaxis()->FindBin( 30 ) - 1 ) << "\n";
    cout << "\tDNP = " << p6_d->Integral( p6_d->GetXaxis()->FindBin( 25 ), p6_d->GetXaxis()->FindBin( 30 ) - 1 ) << "\n";

    cout << "30-40 GeV\n";
    cout << "\tCurrent = " << data->Integral( data->GetXaxis()->FindBin( 30 ), data->GetXaxis()->FindBin( 40 ) - 1 ) << "\n";
    cout << "\tDNP = " << p6_d->Integral( p6_d->GetXaxis()->FindBin( 30 ), p6_d->GetXaxis()->FindBin( 40 ) - 1 ) << "\n";

    c->SaveAs( "./plots/compare_pt_spectrum_raw_R" + rad + "_k" + kappa + ".pdf", "RECREATE" );
}



void compare_ptspectrum_unfold( double k = 0.0, double R = 0.4 )
{
    cout << "\nComparing current unfolded pt spectrum to DNP\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    
    TFile *fdat = new TFile( "./out/unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    
    TFile *fp6 = new TFile( "../../../../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../../../out/p6_matched_R" + rad + "_k" + kappa + ".root" , "READ" );

    
    //TH2D* p62 = (TH2D*) fp6->Get( "QvPt_p" );
    
    
    TH1D* unfs; // unfolded distributions go here
    TH1D* p6_p; // part-level distriubtions from pythia-6 (soon to be unmatched)
    
    
    TH2D* new_tmp = (TH2D*) fdat->Get( "unfold_nom" );
    TH2D* old_tmp = (TH2D*) fp6->Get( "unfold_nom" );

    
    
    
    TCanvas *c = new TCanvas("c_unf_pt_dnp", "canvas", 800, 800);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.45, 0.7, .7, .85);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);

    
    TString pad_name = "pup";
    TPad* p1 = new TPad(pad_name, pad_name, 0, 0.4, 1, 1.0);
    
    
    p1->SetRightMargin(0);
    p1->SetTopMargin(0);
    p1->SetBottomMargin(0);
    p1->SetLeftMargin(0.15);
    
    p1->Draw();
    p1->cd();
    
    
    
    unfs = (TH1D*) new_tmp->ProjectionY( "new_unfolded_ptspectrum" );
    p6_p = (TH1D*) old_tmp->ProjectionY( "old_unfolded_ptspectrum" );
    
    
    cout << "DEBUG: Number of bins\n";
    cout << "\tnew = " << unfs->GetXaxis()->GetNbins() << "\n";
    cout << "\told = " << p6_p->GetXaxis()->GetNbins() << "\n";
    
    
    for(int ik = 0; ik < max( unfs->GetXaxis()->GetNbins(), p6_p->GetXaxis()->GetNbins() ) + 1; ik++ ){
        if( ik < unfs->GetXaxis()->GetNbins() + 1 ){
            cout << "\tnew bin " << ik << " low edge-- " << unfs->GetXaxis()->GetBinLowEdge(ik) << "\n";
        }
        else{
            cout << "\tno more bins in NEW version.\n";
        }
        if( ik < p6_p->GetXaxis()->GetNbins() + 1 ){
            cout << "\told bin " << ik << " low edge-- " << p6_p->GetXaxis()->GetBinLowEdge(ik) << "\n";
        }
        else{
            cout << "\tno more bins in OLD version.\n";
        }
    }
    
    
    unfs->SetMarkerColor(kRed);
    p6_p->SetLineColor(kGreen);
    
    
    TH1D* rat = (TH1D*) unfs->Clone( "new_to_old_ratio"  );
    
    
    unfs->SetMarkerStyle(kFullStar);
//    unfs->SetMarkerSize(2);
    
    p6_p->SetMarkerStyle(kOpenCircle);
//    p6_p->SetLineWidth(3);
    
    
    
    leg1->AddEntry( unfs, "Current (8/23)", "p" );
    leg1->AddEntry( p6_p, "From DNP (10/22)", "p" ); // only data with markers, all else with lines only
    
    
    unfs->SetStats(0);
    //unfs->GetYaxis()->SetRangeUser(0.001, 0.699);
    //unfs->GetXaxis()->SetRangeUser(-4.5, 4.5);
    
            
//     unfs->GetYaxis()->SetTitleSize(0.1);
    unfs->GetYaxis()->SetTitleFont(43);
    unfs->GetYaxis()->SetTitleSize(22);
    unfs->GetYaxis()->SetTitleOffset(1.2);
    unfs->GetYaxis()->SetNdivisions(505);
//       unfs->GetYaxis()->SetLabelSize(0.05);
    unfs->GetYaxis()->SetLabelFont(43);
    unfs->GetYaxis()->SetLabelSize(22);
    unfs->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dp_{T}^{jet} [1/(GeV/c)]" );
    
    
    
    p6_p->SetStats(0);
    
    
    
    p6_p->GetYaxis()->SetTitleFont(43);
    p6_p->GetYaxis()->SetTitleSize(22);
    p6_p->GetYaxis()->SetTitleOffset(1.2);
    p6_p->GetYaxis()->SetNdivisions(505);
    p6_p->GetYaxis()->SetLabelFont(43);
    p6_p->GetYaxis()->SetLabelSize(22);
    p6_p->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dp_{T}^{jet} [1/(GeV/c)]" );
    
    
    
    //p6_p->Sumw2(0);
    p1->SetLogy();
    
    p6_p->Draw("PSAME");
    unfs->Draw("PSAME");
    //line1->Draw();
    
    
    leg1->Draw();
    //leg2->Draw();
    /*
    drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
    drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
        
    leg2->Draw();
    drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
    drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
    
    leg3->Draw();
    drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
    drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
    
    TLine *line2 = new TLine( 0, 0.001, 0, 0.599 );
    line2->SetLineColor(kBlack);
    line2->SetLineStyle(kDashed);
    line2->SetLineWidth(1);
    */
    
    
    c->cd();
    
    TString pdown_name = "pdown";
    TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0.05, 1, 0.4);
    
    p2->SetRightMargin(0);
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.2);
    p2->SetLeftMargin(0.15);
    
    p2->Draw();
    p2->cd();
    
    
    rat->Divide( p6_p );
    rat->GetYaxis()->SetRangeUser(0.501, 1.499);
    rat->SetStats(0);
    rat->SetMarkerStyle( kFullStar );
    rat->SetMarkerColor( kRed );
    
    
    rat->SetTitle( "" );
    
    rat->GetYaxis()->SetTitleFont(43);
    rat->GetYaxis()->SetTitleSize(22);
    rat->GetYaxis()->SetTitleOffset(1.2);
    rat->GetYaxis()->SetNdivisions(505);
    rat->GetYaxis()->SetLabelFont(43);
    rat->GetYaxis()->SetLabelSize(22);
    rat->GetYaxis()->SetTitle( "new / DNP" );
    
    
    
    rat->GetXaxis()->SetTitleFont(43);
    rat->GetXaxis()->SetTitleSize(22);
    rat->GetXaxis()->SetTitleOffset(0.9);
    rat->GetXaxis()->SetNdivisions(505);
    rat->GetXaxis()->SetLabelFont(43);
    rat->GetXaxis()->SetLabelSize(22);
    
    rat->GetXaxis()->SetTitle( "p_{T}^{jet} [GeV/c]" );
   
    
    rat->Draw("PSAME");
    

    c->SaveAs( "./plots/unfoldjetPt_compToDNP_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
    
}



// take ratio of pythia to pythia+geant to see rough estimate of the correction we would expect to see from the unfolding
void estimate_correction( double k = 0.0, double R = 0.4 )
{
    
    
}




// produce jet charge resolution, scale "delta Q vs pT" (similar to Isaac's jet mass)
void plot_resolution_jetcharge( double k = 0.0, double R = 0.4 )
{
    
}




// make plot with different iterations
void plot_closure_jetcharge( double k = 0.0, double R = 0.4 )
{
    cout << "\nPlotting jet charge closure\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    
    //TFile *fdat = new TFile( "./out/data_unfold_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    TFile *fdat = new TFile( "./out/unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    
    
    
    TH1D* clos[njetbins];
    TH1D* same[njetbins];
    
    
    TCanvas *c = new TCanvas("c_closure", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.4, 0.65, .8, .8);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    
    
    TLegend *leg2 = new TLegend(0.4, 0.65, .8, .8);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    
    TLegend *leg3 = new TLegend(0.4, 0.65, .8, .8);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.05);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString padup_name = Form( "pup_%i", j );
        TPad* padup = new TPad(padup_name, padup_name, 0, 0, 1, 1);
        
        padup->SetRightMargin(0);
        padup->SetTopMargin(0);
        padup->SetBottomMargin(0.15);
        padup->SetLeftMargin(0);
        if( j == 0 ){padup->SetLeftMargin(0.2);}
        
        padup->Draw();
        padup->cd();

        
        
        clos[j] = (TH1D*) fdat->Get( "opp_4iterx" + outFile_pt[j] );
        same[j] = (TH1D*) fdat->Get( "same_4iterx" + outFile_pt[j] );
        
        
        clos[j]->SetMarkerColor(kRed);
        clos[j]->SetMarkerStyle(kOpenCircle);
        
        
        same[j]->SetMarkerColor(kBlack);
        same[j]->SetMarkerStyle(kOpenCircle);
        
        
        if( j == 0 ){
            leg1->AddEntry(clos[j], "opp-side", "p");
            leg1->AddEntry(same[j], "same-side", "p");
        }
        
        
        same[j]->SetStats(0);
        clos[j]->SetStats(0);
        clos[j]->GetYaxis()->SetRangeUser(0.501, 1.499);
            
        
        clos[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        
        
        clos[j]->GetYaxis()->SetTitleFont(43);
        clos[j]->GetYaxis()->SetTitleSize(22);
        clos[j]->GetYaxis()->SetTitleOffset(1.2);
        clos[j]->GetYaxis()->SetNdivisions(505);
        clos[j]->GetYaxis()->SetLabelFont(43);
        clos[j]->GetYaxis()->SetLabelSize(22);
        clos[j]->GetYaxis()->SetTitle( "Closure" );
        
        
        
        clos[j]->Draw("PSAME");
        same[j]->Draw("PSAME");
            //line1->Draw();
        
        //TLatex *tprelim = new TLatex();
        //tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
            
        }
        else if( j == 1){
            //leg2->Draw();
            drawText( "nominal 4 iterations" , .05, 0.85, 16 );
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            //leg3->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( -4.5, 1.0, 4.5, 1.0 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        line2->Draw();
        
    }

    c->SaveAs( "./plots/closure_jetcharge_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}

// make closure plot again to diagnose unfolding procedure
void plot_iterCheck_jetcharge( TString side = "opp", double k = 0.0, double R = 0.4 )
{
    cout << "\nPlotting jet charge closure\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    
    //TFile *fdat = new TFile( "./out/data_unfold_correct_R" + rad + "_k" + kappa + ".root" , "READ");
    TFile *fdat = new TFile( "./out/unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    
    int nIters = 10; // closure is done for [1, 10] iterations (4 = nominal unfolding)
    int nleg_cols = 3;
    
    TH1D* clos[nIters][njetbins];
    
    
    TCanvas *c = new TCanvas("c_itercheck", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.4, 0.65, .8, .8);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    leg1->SetNColumns(nleg_cols);

    
    TLegend *leg2 = new TLegend(0.4, 0.65, .8, .8);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    leg2->SetNColumns(nleg_cols);
    
    
    TLegend *leg3 = new TLegend(0.4, 0.65, .8, .8);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.05);
    leg3->SetNColumns(nleg_cols);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString padup_name = Form( "pup_%i", j );
        TPad* padup = new TPad(padup_name, padup_name, 0, 0, 1, 1);
        
        padup->SetRightMargin(0);
        padup->SetTopMargin(0);
        padup->SetBottomMargin(0.15);
        padup->SetLeftMargin(0);
        if( j == 0 ){padup->SetLeftMargin(0.2);}
        
        padup->Draw();
        padup->cd();

        
        for(int ii = 0; ii < nIters; ii++){
            clos[ii][j] = (TH1D*) fdat->Get( Form( side + "_%iiterx" + outFile_pt[j], ii+1 ) );
            
            clos[ii][j]->SetMarkerColor(ii+1);
            clos[ii][j]->SetMarkerStyle(kOpenCircle);
            
            if( j == 0 ){
                if(ii%nleg_cols == 0){
                    leg1->AddEntry( clos[ii][j], Form( "%i", ii+1) , "p" );
                }
                else if(ii%nleg_cols == 1){
                    leg2->AddEntry( clos[ii][j], Form( "%i", ii+1) , "p" );
                }
                else if(ii%nleg_cols == 2){
                    leg3->AddEntry( clos[ii][j], Form( "%i", ii+1) , "p" );
                }
            }
            
            
            
            clos[ii][j]->SetStats(0);
            clos[ii][j]->GetYaxis()->SetRangeUser(0.001, 1.999);
            
            
            clos[ii][j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
            
            
            clos[ii][j]->GetYaxis()->SetTitleFont(43);
            clos[ii][j]->GetYaxis()->SetTitleSize(22);
            clos[ii][j]->GetYaxis()->SetTitleOffset(1.2);
            clos[ii][j]->GetYaxis()->SetNdivisions(505);
            clos[ii][j]->GetYaxis()->SetLabelFont(43);
            clos[ii][j]->GetYaxis()->SetLabelSize(22);
            clos[ii][j]->GetYaxis()->SetTitle( "Closure" );
            
            
            
            clos[ii][j]->Draw("PSAME");
            //line1->Draw();
        }
        
        //TLatex *tprelim = new TLatex();
        //tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
        
        if( j == 0 ){
            leg1->Draw();
            drawText( "p+p #sqrt{s} = 200 GeV", .36, .9, 20 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .34, .81, 15 );
            
        }
        else if( j == 1){
            leg2->Draw();
            drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .05, .9, 18 );
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
        }
        else if( j == 2 ){
            leg3->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
            drawText( "p_{T}^{cons} > 200 MeV/c", .05, .9, 18 );
        }
        
        TLine *line2 = new TLine( -4.5, 1.0, 4.5, 1.0 );
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(1);
        
        
        line2->Draw();
        
    }

    c->SaveAs( "./plots/" + side + "_itercheck_jetcharge_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
}



// look at the response, as well as the fakes and misses distributions





void debug_plots( double k = 0.0, double R = 0.4 )
{
    //plot_unfPart_jetcharge( k , R );
    //plot_dataDet_jetcharge( k , R );
    //plot_conspt( k , R );
    //plot_nef( k , R );
    //plot_jetspectrum( k , R );
    
    
    compare_to_preDNP_unfold( k , R );
    compare_to_preDNP_raw( k , R );
    compare_to_preDNP_pythia( k , R );
    compare_to_preDNP_geant( k , R );
    
    //plot_closure_jetcharge( k , R );
    
    
    //plot_prelim_det( k , R );
    

    //compare_jet_pt_counts( k , R );
    
    compare_ptspectrum_unfold( k , R );
    
    
}

