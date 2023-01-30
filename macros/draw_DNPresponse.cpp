//
//  draw_DNPresponse.C
//  
//
//  Created by Grant McNamara on 10/9/22.
//
//
// draw Response matrix used for jet charge unfolding for DNP 2022 slides

#include <string>
#include <iostream>
#include "math.h"
#include "Headers/plot.h"
#include "Headers/JCparameters.hh"


using namespace jcAnalysis;


void draw_DNPresponse( double k = 0.0, double R = 0.4 ){
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    gROOT->ForceStyle();
    gStyle->SetPalette(kPastel);
    
    
    // need to change file path, name to where the output file of unfold.cxx is on the grid
    TFile* fres = new TFile( "../out/unfold/unfolded_R" + rad + "_k" + kappa + ".root" , "READ" );
    
    
    RooUnfoldResponse *res = (RooUnfoldResponse*) fres->Get( "q_pt_response" );
    
    TLatex *p = new TLatex();
    
    TCanvas *cres = new TCanvas("cres", "", 600, 600);
    TH1D* res_mat = (TH1D*) res->Hresponse();
    
    res_mat->SetStats(0);
    
    cres->SetLeftMargin(0.2);
    cres->SetBottomMargin(0.2);
    cres->SetRightMargin(0.2);
    cres->SetTopMargin(0.15);
    cres->SetLogz();
    
    
    res_mat->GetZaxis()->SetRangeUser( pow(10, -10), pow(10, -5) );
    
//    res_mat->GetYaxis()->SetTitleSize(0.05);
//    res_mat->GetYaxis()->SetTitleOffset(0.);
    // /*
    res_mat->GetYaxis()->SetLabelSize(0.0);
    res_mat->GetYaxis()->SetTicks("U");
    res_mat->GetYaxis()->SetTickLength(0.0);
    res_mat->GetXaxis()->SetLabelSize(0.0);
    res_mat->GetXaxis()->SetTicks("U");
    res_mat->GetXaxis()->SetTickLength(0.0);
    // */
//    res_mat->GetXaxis()->SetTitleSize(0.05);
//    res_mat->GetXaxis()->SetTitleOffset(0.1);
//    res_mat->GetXaxis()->SetLabelSize(0.05);
    
    
    double ymax = res_mat->GetYaxis()->GetBinUpEdge( res_mat->GetYaxis()->GetNbins() );
    double ymin = res_mat->GetYaxis()->GetBinLowEdge( 1 );
    double xmax = res_mat->GetXaxis()->GetBinUpEdge( res_mat->GetXaxis()->GetNbins() );
    double xmin = res_mat->GetXaxis()->GetBinLowEdge( 1 );
    
    // x-axis
    TF1 *f1 = new TF1("f1", "x", 15, 25 );
    TGaxis *A1 = new TGaxis( 0 , 0, (13.0*4.0)/2.0, 0, "f1", 4 );
    //A1->SetTitle( "p^{det}_{T}" );
    //A1->SetLabelSize(0.03);
    //A1->SetTitleSize(0.03);
    //A1->SetTitleOffset(1.2);

    TF1 *f2 = new TF1("f2", "x", 30, 40 );
    TGaxis *A2 = new TGaxis( 3.0/4.0 * (13.0*4.0) , 0, (13.0*4.0), 0, "f2", 2 );
    A2->SetTitle( "p^{jet, det}_{T} [GeV/c]" );
    //A2->SetLabelSize(0.03);
    A2->SetTitleSize(0.08);
    A2->SetTitleOffset(0.8);
    //A2->Set
    
    //                                 max - min / 5 GeV bins * 13 Q bins per jet bin
    res_mat->GetYaxis()->SetRangeUser(0, (50.0-5.0) * (13.0/5.0) ); // want to show 5-10, 10-15, 15-20, 20-25, 25-30, 30-35, 35-40, 40-45, 45-50?
    //res_mat->GetXaxis()->SetRangeUser(0, 5.0 * 13.0 ); // want to show 15-20, 20-25, 25-30, 30-40, 40-50?
    
    // y-axis
    TF1 *f3 = new TF1("f3", "x", 5, 50 );
    TGaxis *A3 = new TGaxis( 0, 0 , 0, (50.0-5.0) * (13.0/5.0), "f3", 10 );
    A3->SetTitle( "p^{jet,part}_{T} [GeV/c]" );
    //A3->SetLabelSize(0.03);
    A3->SetTitleSize(0.08);
    A3->SetTitleOffset(0.8);
    
    

    res_mat->GetZaxis()->SetTitle( "1/#sigma d#sigma/dQdp_{T} [1/e c/GeV]" );
    res_mat->GetZaxis()->SetTitleSize(0.05);
    res_mat->GetZaxis()->SetTitleOffset(1.3);
    res_mat->SetTitle( "4D jet charge response matrix" );
    res_mat->SetTitleOffset(0.05);

    res_mat->Draw("colz");

    A1->Draw();
    A2->Draw();
    A3->Draw();
    
    //p->DrawLatexNDC(0.2,0.8, "pp #sqrt{s} = 200 GeV");
    //p->DrawLatexNDC(0.2,0.75, "anti-k_{T}, R = 0.4, |#eta| < 0.4");
    
    TLatex *tprelim = new TLatex();
    tprelim->SetTextSize(0.05); tprelim->SetTextColor(kRed);
    tprelim->DrawLatexNDC(0.31, 0.87, "STAR Simulation");
    
    drawText( "pp #sqrt{s} = 200 GeV", 0.23, 0.8, 20);
    drawText("anti-k_{T}, R = 0.4, |#eta| < 0.4", 0.23, 0.75, 20);
    

    
    cres->SaveAs( "./plots/dnp_jetcharge_responseMatrix.png" , "RECREATE" );
}

