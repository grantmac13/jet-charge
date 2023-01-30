//
//  temp_DNPsystematicPlotting.C
//  
//
//  Created by Grant McNamara on 10/8/22.
//


#include <string>
#include <iostream>
#include "math.h"
#include "Headers/plot.h"
#include "Headers/JCparameters.hh"


using namespace jcAnalysis;

/*
// from Isaac to keep consistency
//projects a 2D histogram in desired ranges and returns an array of the (1D) projections on desired axis
std::vector<TH1D*> Projection2D (TH2D * hist2D, const int nBins, double * ranges, const std::string axis) {
  std::vector<TH1D*> proj1Ds;
  for (int i = 0; i < nBins; ++ i) {
    std::string low = std::to_string(ranges[i]);
    std::string high = std::to_string(ranges[i+1]);
    std::string low_rough = low.substr(0,2);
    std::string high_rough = high.substr(0,2);
    if (low_rough.substr(1,2) == ".") {low_rough = low_rough.substr(0,1);}
    if (high_rough.substr(1,2) == ".") {high_rough = high_rough.substr(0,1);}
    if (axis == "x" || axis == "X" || axis == "1") {
      cout << "now including e option" << endl;
      proj1Ds.push_back(hist2D->ProjectionX((hist2D->GetName() + axis + low_rough + high_rough).c_str(),ranges[i],ranges[i+1]- 1,"e"));
    }
    else if (axis == "y" || axis == "Y" || axis == "2") {
      cout << "now including e option" << endl;
      proj1Ds.push_back(hist2D->ProjectionY((hist2D->GetName() + axis + low_rough + high_rough).c_str(),ranges[i],ranges[i+1] - 1,"e"));
    }
    else {
      std::cerr << "Improper axis given for projections. Exiting." << std::endl; exit(1);
    }
    proj1Ds[i]->SetTitle("");
  }
  return proj1Ds;
}
*/


// from Grant's unfold.cxx
//projects a 2D histogram in desired ranges and returns an array of the (1D) projections on desired axis
std::vector<TH1D*> Projection2D (TH2D * hist2D, const int nBins, double * ranges, const std::string axis) {
    std::vector<TH1D*> proj1Ds;
    for (int i = 0; i < nBins; ++ i) {
        std::string low = std::to_string(ranges[i]);
        std::string high = std::to_string(ranges[i+1]);
        std::string low_rough = low.substr(0,2);
        std::string high_rough = high.substr(0,2);
        if (low_rough.substr(1,2) == ".") {low_rough = low_rough.substr(0,1);}
        if (high_rough.substr(1,2) == ".") {high_rough = high_rough.substr(0,1);}
        //cerr << "low edge: " << low_rough << "\n";
        //cerr << "high edge: " << high_rough << "\n";
        if (axis == "x" || axis == "X" || axis == "1") {
            //cout << "now including e option" << endl;
            proj1Ds.push_back(hist2D->ProjectionX((hist2D->GetName() + axis + low_rough + high_rough).c_str(), hist2D->GetYaxis()->FindBin(ranges[i]), hist2D->GetYaxis()->FindBin(ranges[i+1]) - 1, "e"));
            //cout << "for bin: [" << ranges[i] << ", " << ranges[i+1] << "]\n";
            //cout << "project onto y bins [" << hist2D->GetYaxis()->FindBin(ranges[i]) << ", " << hist2D->GetYaxis()->FindBin(ranges[i+1]) - 1 << "]\n";
            // testing above:
//            proj1Ds.push_back(hist2D->ProjectionX((hist2D->GetName() + axis + low_rough + high_rough).c_str(), ranges[i], ranges[i+1] - 1, "e"));
        }
        else if (axis == "y" || axis == "Y" || axis == "2") {
            cout << "now including e option" << endl;
            proj1Ds.push_back(hist2D->ProjectionY((hist2D->GetName() + axis + low_rough + high_rough).c_str(), ranges[i], ranges[i+1] - 1, "e"));
        }
        else {
            std::cerr << "Improper axis given for projections. Exiting." << std::endl; exit(1);
        }
        proj1Ds[i]->SetTitle("");
        // debug

    }
    return proj1Ds;
}



void systematics( double k = 0.0, double R = 0.4 ){
    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    gROOT->ForceStyle();
    gStyle->SetPalette(kPastel);
    
    
    // file will be loaded within function based on R, kappa
    TFile* fres = new TFile( "../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ" );
    
    // don't load roounfold objects on local machine due to roounfold not being setup, just pull already unfolded histograms by systematic -- unfold_Xx2025 histograms already in form of percentage difference compared to nominal
  
    TH1D* reco_noms[njetbins];
    TH1D* reco_IP2s[njetbins];
    TH1D* reco_IP6s[njetbins];
    TH1D* reco_TSs[njetbins];
    TH1D* reco_TUs[njetbins];
    TH1D* reco_HC50s[njetbins];
    TH1D* reco_DSs[njetbins];
    TH1D* reco_GSs[njetbins];
    
    TH1D* reco_h7smear1D[njetbins];
    TH1D* reco_p8smear1D[njetbins];
    
    
    vector<TH1D*> reco_noms_copy;
    
    
    for(int j = 0; j < njetbins; j++){
        reco_noms[j] = (TH1D*) fres->Get( "unfold_nomx" + outFile_pt[j] );
        reco_IP2s[j] = (TH1D*) fres->Get( "unfold_IP2x" + outFile_pt[j] );
        reco_IP6s[j] = (TH1D*) fres->Get( "unfold_IP6x" + outFile_pt[j] );
        reco_TSs[j] = (TH1D*) fres->Get( "unfold_TSx" + outFile_pt[j] );
        reco_TUs[j] = (TH1D*) fres->Get( "unfold_TUx" + outFile_pt[j] );
        reco_HC50s[j] = (TH1D*) fres->Get( "unfold_HC50x" + outFile_pt[j] );
        reco_DSs[j] = (TH1D*) fres->Get( "unfold_DSx" + outFile_pt[j] );
        reco_GSs[j] = (TH1D*) fres->Get( "unfold_GSx" + outFile_pt[j] );
        
        reco_h7smear1D[j] = (TH1D*) fres->Get( "unfold_h7smear1D_" + outFile_pt[j] );
        reco_p8smear1D[j] = (TH1D*) fres->Get( "unfold_p8smear1D_" + outFile_pt[j] );
    
        
        reco_noms_copy.push_back((TH1D*) reco_noms[j]->Clone( "nom_copy" + outFile_pt[j] ) );
        reco_noms_copy[j]->Scale( 1 / (double)reco_noms_copy[j]->Integral() );
        
        
        reco_h7smear1D[j]->SetLineWidth(2);
        reco_h7smear1D[j]->SetLineColor(11); reco_h7smear1D[j]->SetLineStyle(kSolid);
        reco_h7smear1D[j]->SetFillColor(11); reco_h7smear1D[j]->SetFillStyle(3690);
        
        reco_p8smear1D[j]->SetLineWidth(2);
        reco_p8smear1D[j]->SetLineColor(11); reco_p8smear1D[j]->SetLineStyle(kSolid);
        reco_p8smear1D[j]->SetFillColor(11); reco_p8smear1D[j]->SetFillStyle(3444);
        

        
        reco_IP2s[j]->SetFillColor(2); reco_IP2s[j]->SetFillStyle(3305);
        reco_IP6s[j]->SetFillColor(3); reco_IP6s[j]->SetFillStyle(3395);
        reco_TSs[j]->SetFillColor(kGreen+2); reco_TSs[j]->SetFillStyle(3490);
        reco_TUs[j]->SetFillColor(kMagenta+1); reco_TUs[j]->SetFillStyle(3436);
        reco_HC50s[j]->SetFillColor(6); reco_HC50s[j]->SetFillStyle(3335);
        reco_DSs[j]->SetFillColor(8); reco_DSs[j]->SetFillStyle(3944);
        reco_GSs[j]->SetFillColor(9); reco_GSs[j]->SetFillStyle(3544);
        
        
    }
    
    TLegend *l1 = new TLegend(0.5,0.3,0.7,0.6); l1->SetBorderSize(0);
    l1->AddEntry(reco_IP2s[0],"IP2","f");
    l1->AddEntry(reco_IP6s[0],"IP6","f");
    l1->AddEntry(reco_TSs[0],"TS","f");
    l1->AddEntry(reco_TUs[0],"TU","f");

    TLegend *l2 = new TLegend(0.5,0.3,0.7,0.6); l2->SetBorderSize(0);
    l2->AddEntry(reco_HC50s[0],"HC50","f");
    l2->AddEntry(reco_DSs[0],"DS","f");
    l2->AddEntry(reco_GSs[0],"GS","f");
    //l2->AddEntry(reco_MSs[0],"MS","f");
    l2->AddEntry(reco_p8smear1D[0],"P8","f");
    l2->AddEntry(reco_h7smear1D[0],"H7","f");
    

    TLatex *slice = new TLatex();

    TCanvas *csys = new TCanvas("csys", "csys", 800, 500);
//    DivideCanvas(csys,"0", nBins,1);
    csys->SetLeftMargin(0);
    
    csys->Divide(njetbins, 1, 0, 0);
    

    TH1D* ldummy = new TH1D("ldummy", ";;relative uncertainty", 13, -6.5, 6.5);
    TH1D* rdummy = new TH1D("rdummy", ";;relative uncertainty", 13, -6.5, 6.5);
    
    
    ldummy->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
    rdummy->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
    
    double pad_max = 0.599;
    
    ldummy->GetXaxis()->SetRangeUser(-4.5,4.5);
    ldummy->GetYaxis()->SetRangeUser(0.001,pad_max);
    rdummy->GetXaxis()->SetRangeUser(-4.5,4.5);
    rdummy->GetYaxis()->SetRangeUser(0.001,pad_max);
    
    
//    ldummy->GetYaxis()->SetTitleSize(0.08);
    ldummy->GetYaxis()->SetTitleFont(43);
    ldummy->GetYaxis()->SetTitleSize(22);
    ldummy->GetYaxis()->SetTitleOffset(1.2);
    ldummy->GetYaxis()->SetNdivisions(505);
    rdummy->GetYaxis()->SetNdivisions(505);
//    ldummy->GetYaxis()->SetLabelSize(0.05);
    ldummy->GetYaxis()->SetLabelFont(43);
    ldummy->GetYaxis()->SetLabelSize(22);

    
    ldummy->GetXaxis()->SetTitleFont(43);
    ldummy->GetXaxis()->SetTitleSize(22);
    ldummy->GetXaxis()->SetTitleOffset(1.1);
    ldummy->GetXaxis()->SetLabelFont(43);
    ldummy->GetXaxis()->SetLabelSize(22);
    //ldummy->GetXaxis()->SetLabelOffset(0.06);
    
    rdummy->GetXaxis()->SetTitleFont(43);
    rdummy->GetXaxis()->SetTitleSize(22);
    rdummy->GetXaxis()->SetTitleOffset(1.1);
    rdummy->GetXaxis()->SetLabelFont(43);
    rdummy->GetXaxis()->SetLabelSize(22);
    //rdummy->GetXaxis()->SetLabelOffset(0.06);
    
    
    ldummy->SetStats(0);
    rdummy->SetStats(0);
    
    
    TLatex *p = new TLatex();
    //p->SetTextAlign(11);
    //p->SetTextSize(0.07);

    TLatex *t = new TLatex();
    for (int i = 0; i < njetbins; ++ i) {
        csys->cd(i+1);
        if (i == 0) {ldummy->Draw();} if (i != 0) {rdummy->Draw();}
        reco_HC50s[i]->Draw("LF2same");
        reco_TSs[i]->Draw("LF2same");
        reco_TUs[i]->Draw("LF2same");
        reco_DSs[i]->Draw("LF2same");
        reco_GSs[i]->Draw("LF2same");
        reco_IP2s[i]->Draw("LF2same");
        reco_IP6s[i]->Draw("LF2same");
        reco_p8smear1D[i]->Draw("LF2same");
        reco_h7smear1D[i]->Draw("LF2same");
        //reco_MSs[i]->Draw("LF2same");
        if (i == 0) {
            p->DrawLatexNDC(0.3,0.9, "pp 200 GeV run12 JP2");
            p->DrawLatexNDC(0.3,0.85, "anti-k_{T}, R = 0.4");
            p->DrawLatexNDC(0.3,0.8, "Ch+Ne jets, |#eta| < 0.4");
        }
        if (i==0) {l1->Draw("same");}
        if (i==1) {
            //reco_p8smear1D_3045->Draw("LF2same");
            //reco_h7smear1D_3045->Draw("LF2same");
            l2->Draw("same");
        }
        slice->DrawLatexNDC( 0.3, 0.7, Form( "%1.0f < p_{T} < %1.0f GeV/c", jetPtLo[i], jetPtHi[i] ) );
    }

    //csys->SaveAs(("~/jetmass2_11-10-2020_11_10-2020/plots/DNP_talk/"+fstart+"systematics_new.pdf").c_str());
    
    
    double edges[njetbins+1];
    for(int j = 0; j < njetbins+1; j++){
        edges[j] = jetEdges[j];
    }
    
    
    //taking maximum envelopes!
    TH2D* env_HC = new TH2D("env_HC","",13,-6.5,6.5,15,5,80);
    vector<TH1D*> env_HCs = Projection2D(env_HC,njetbins,edges,"x");
    TH2D* env_un = new TH2D("env_un","",13,-6.5,6.5,15,5,80);
    vector<TH1D*> env_uns = Projection2D(env_un,njetbins,edges,"x");
    TH2D* net = new TH2D("net","",13,-6.5,6.5,15,5,80);
    vector<TH1D*> nets = Projection2D(net,njetbins,edges,"x");
    TH2D* stat2D = new TH2D("stat","",13,-6.5,6.5,15,5,80);
    vector<TH1D*> stats = Projection2D(stat2D,njetbins,edges,"x");
    
    //  const int nBins_m = 14;
    
    vector<vector< double> > syst_errs2D;
    
    for (int i = 0; i < njetbins; ++ i) {
        vector<double> syst_errs1D;

        for (int j = 1; j < 13 + 1; ++ j) {
            //hadronic correction envelope - using an ordered set here to automatically get the largest value.
            double hcs [1] = {reco_HC50s[i]->GetBinContent(j)};
            set<double> hc_sort (hcs, hcs+1);
            set<double>::iterator hc = hc_sort.end(); hc --;
            double hc_envelope = *hc;
            env_HCs[i]->SetBinContent(j, hc_envelope);
            
            
            //unfolding envelope - using an ordered set here to automatically get the largest value
            double uns [6] = {reco_DSs[i]->GetBinContent(j), reco_GSs[i]->GetBinContent(j), /*reco_MSs[i]->GetBinContent(j),*/ reco_h7smear1D[i]->GetBinContent(j), reco_p8smear1D[i]->GetBinContent(j), reco_IP2s[i]->GetBinContent(j), reco_IP6s[i]->GetBinContent(j)};
            
            set<double> un_sort (uns, uns+6);
            set<double>::iterator un = un_sort.end(); un --;
            double un_envelope = *un;
            if(un_envelope > pad_max){
                //cout << "jet bin: " << i << " and histogram bin: " << j << "\n";
                //cout << "uncertainty > 100%: " << un_envelope << "\n";
                un_envelope = pad_max;
            }
            env_uns[i]->SetBinContent(j, un_envelope);
            //cout << "unfolding uncertainty value for bin = " << un_envelope << "\n";
            
            //total uncertainty = TU + TS + un envelope + hc envelope
            double square = pow(hc_envelope,2) + pow(un_envelope,2) + pow(reco_TUs[i]->GetBinContent(j),2) + pow(reco_TSs[i]->GetBinContent(j),2);
            nets[i]->SetBinContent(j,sqrt(square));
            stats[i]->SetBinContent(j,reco_noms_copy[i]->GetBinError(j));
            syst_errs1D.push_back(nets[i]->GetBinContent(j));
        }
        
        syst_errs2D.push_back(syst_errs1D);
    }
    
    for (int i = 0; i < njetbins; ++ i) {
        
        env_HCs[i]->SetLineColor(kRed); env_HCs[i]->SetLineStyle(kSolid); env_HCs[i]->SetLineWidth(2);
        env_HCs[i]->SetFillColor(kRed); env_HCs[i]->SetFillStyle(3353);
        
        env_uns[i]->SetLineColor(kRed); env_uns[i]->SetLineStyle(kSolid); env_uns[i]->SetLineWidth(2);
        env_uns[i]->SetFillColor(kBlue); env_uns[i]->SetFillStyle(3305);
        
        
        stats[i]->SetLineColor(kBlack); stats[i]->SetLineStyle(kDashed); stats[i]->SetLineWidth(2);
        nets[i]->SetLineColor(kBlack); nets[i]->SetLineStyle(kSolid); nets[i]->SetLineWidth(2);
        
    }
    
    TLegend *tenvs = new TLegend(0.2,0.35,0.8,0.85); tenvs->SetBorderSize(0);
    tenvs->AddEntry(env_HCs[0],"Hadronic correction","f");
    tenvs->AddEntry(reco_TSs[0],"Tower energy scale","f");
    tenvs->AddEntry(reco_TUs[0],"Tracking","f");
    tenvs->AddEntry(env_uns[0],"Unfolding","f");
    tenvs->AddEntry(nets[0],"Total systematic","l");
    
    tenvs->SetTextSize(0.07);
    
    
    TCanvas *cenvs = new TCanvas("cenvs","cenvs",800,500);
    cenvs->SetLeftMargin(0.2);
    cenvs->SetBottomMargin(0.15);
    cenvs->Divide(njetbins, 1, 0, 0);
    
    
    for (int i = 0; i < njetbins; ++ i) {
        cenvs->cd(i+1);
        TString pdown_name = Form( "psyst_%i", i );
        TPad* p2 = new TPad(pdown_name, pdown_name, 0, 0, 1, 1);
        
        p2->SetRightMargin(0);
        p2->SetTopMargin(0);
        p2->SetBottomMargin(0.15);
        p2->SetLeftMargin(0);
        if( i == 0 ){p2->SetLeftMargin(0.2);}
        
        p2->Draw();
        p2->cd();
        
        
        if (i == 0) {ldummy->Draw();} if (i != 0) {rdummy->DrawCopy();}
        env_uns[i]->Draw("lf2same PLC PFC PMC");
        reco_TUs[i]->Draw("lf2same PLC PFC PMC");
        env_HCs[i]->Draw("lf2same PLC PFC PMC");
        reco_TSs[i]->Draw("lf2same PLC PFC PMC");
        nets[i]->Draw("lf2same");
        
        if (i == 0) {
            //slice->DrawLatexNDC( 0.34, 0.7, Form( "%1.0f GeV/c < p_{T} < %1.0f GeV/c", jetPtLo[i], jetPtHi[i] ) );
            //slice->DrawLatexNDC(0.3,0.7,(jetPtLo[i]+" < p_{T} < "+jetPtHi[i]+" GeV/c").c_str());
            drawText( Form( "%1.0f GeV/c < p_{T} < %1.0f GeV/c", jetPtLo[i], jetPtHi[i] ), 0.34, 0.7, 15 );

            drawText( "pp 200 GeV 2012", 0.4, 0.9, 20 );
            drawText( "anti-k_{T}, R = 0.4", 0.42, 0.85, 20 );
            drawText( "|#eta| < 0.6", 0.48, 0.8, 20 );

            //p->DrawLatexNDC(0.43,0.9, "pp 200 GeV 2012");
            //p->DrawLatexNDC(0.45,0.85, "anti-k_{T}, R = 0.4");
            //p->DrawLatexNDC(0.52,0.8, "|#eta| < 0.6");
        }
        if (i == 1) {
            tenvs->Draw("same");
            //slice->DrawLatexNDC( 0.22, 0.9, Form( "%1.0f GeV/c < p_{T} < %1.0f GeV/c", jetPtLo[i], jetPtHi[i] ) );
            drawText( Form( "%1.0f GeV/c < p_{T} < %1.0f GeV/c", jetPtLo[i], jetPtHi[i] ), 0.22, 0.9, 15 );
        }
        else if (i == 2){
            TLatex *tprelim = new TLatex();
            tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
            tprelim->DrawLatex(-2.0, 0.55, "STAR Preliminary");
            //slice->DrawLatexNDC( 0.22, 0.7, Form( "%1.0f GeV/c < p_{T} < %1.0f GeV/c", jetPtLo[i], jetPtHi[i] ) );
            drawText( Form( "%1.0f GeV/c < p_{T} < %1.0f GeV/c", jetPtLo[i], jetPtHi[i] ), 0.21, 0.8, 15 );
        }
        TLine *line = new TLine( 0, 0.001, 0, 0.599 );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        
        //line->Draw();
        
    }
  
    cenvs->SaveAs( "./plots/jetQ_systematic_envelopes_R" + rad + "_k" + kappa + ".pdf" );
  
}



void temp_DNPsystematicPlotting(){
    systematics(0.0, 0.4);
}


