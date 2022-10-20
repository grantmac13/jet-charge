
#include <string>
#include <iostream>
#include "math.h"
#include "Headers/plot.h"
#include "Headers/JCparameters.hh"

using namespace jcAnalysis;


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


//isaac's method of systematic calculation (envelopes) return histograms to use in plot with bin content = unfold_nom bin content, bin error = systematic error
vector<TH1D*> systematics( double k = 0.0, double R = 0.4 ){
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
        reco_noms_copy[j]->Scale( 1 / (double)reco_noms_copy[j]->Integral( 0, reco_noms_copy[j]->GetNbinsX() + 1 ) );
        
        
    }
    
    
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
            if(un_envelope > 0.999){
                //cout << "jet bin: " << i << " and histogram bin: " << j << "\n";
                //cout << "uncertainty > 100%: " << un_envelope << "\n"; un_envelope = 0.999;
            }
            env_uns[i]->SetBinContent(j, un_envelope);
            //cout << "unfolding uncertainty value for bin = " << un_envelope << "\n";
            
            //total uncertainty = TU + TS + un envelope + hc envelope
            double square = pow(hc_envelope,2) + pow(un_envelope,2) + pow(reco_TUs[i]->GetBinContent(j),2) + pow(reco_TSs[i]->GetBinContent(j),2);
            nets[i]->SetBinContent(j, reco_noms_copy[i]->GetBinContent(j));
            nets[i]->SetBinError(j, sqrt(square) * reco_noms_copy[i]->GetBinContent(j) );
            stats[i]->SetBinContent(j,reco_noms_copy[i]->GetBinError(j));
            syst_errs1D.push_back(nets[i]->GetBinContent(j));
        }
    }
    return nets;
}



// calculate quadrature sum of systematic errors from all sources and set bin error of q_err of j-th jet pt bin histogram to the sqrt(quad_sum)
void systErrs( TH1D *syst[nsyst][njetbins], TH1D *q_err[njetbins], int j ){
    // calculate errors within jet pt loop (calculate at each j as the jet pt bin loops
    // set errors on q_err histogram to the quadrature sum of the systematic errors multiplied by the bin content of q_err (or the unf histogram in the main code, should be equivalent)
    
    // systematics order: IP2, IP6, TS, TU, HC50, DS, GS, h7, p8
    
    
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


// calculate maximum envelope for unfolding uncertainty, quadrature sum TS, TU, unfolding, HC50 systematic errors
void systEnvs( TH1D *syst[nsyst][njetbins], TH1D *q_err[njetbins], int j ){
    // calculate errors within jet pt loop (calculate at each j as the jet pt bin loops
    // set errors on q_err histogram to the quadrature sum of the systematic errors multiplied by the bin content of q_err (or the unf histogram in the main code, should be equivalent)
    
    // systematics order: IP2, IP6, TS, TU, HC50, DS, GS, h7, p8

    
    for(int i = 1; i < (int) syst[0][j]->GetNbinsX(); i++){
        double un_envelope = 0.0; // initialize unfolding maximum envelope to 0 for each histogram bin
        
        double err_quadsum = 0.0; // total systematic error to apply to bin
        
        
        for(int ki = 0; ki < nsyst; ki++){
            if(ki == 0 || ki == 1 || ki >= 5){
                //cout << "SANITY CHECK: systematic part of unfolding uncertainty envelope? " << syst_names[ki] << "\n";
                if(syst[ki][j]->GetBinContent(i) > err_quadsum){ // if ki-th syst is larger than current largest systematic for bin, replace un_envelope with largest syst
                    un_envelope = syst[ki][j]->GetBinContent(i);
                }
            }
            else{ // if ki-th systematic is not IP2(6), D(G)S, h7/p8
                double syst_err = syst[ki][j]->GetBinContent(i);
                err_quadsum += syst_err*syst_err; // total error is sqrt(sum of individual errors ^2)
            }
        }
        err_quadsum += un_envelope*un_envelope;
        q_err[j]->SetBinError( i, sqrt( err_quadsum ) * q_err[j]->GetBinContent(i) );
        // err_quadsum is a percentage of the bin content, so to calculate absolute error have to multiply sqrt(err_quadsum) * bincontent
    }
    
}



void plot_unfPart_jetcharge(double k = 0.0, double R = 0.4){
    cout << "\nPlotting unfolded data and PYTHIA on one canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    // file with pythia-6 distribution
    // unmatched pythia-6 --- only has part-level distribution...
    TFile *fp6 = new TFile( "../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    //TFile *fp6 = new TFile( "../out/p6_R" + rad + "_k" + kappa + ".root" , "READ" );

    // adding pythia-8 unmatched to this figure
    TFile *fp8 = new TFile( "../out/p8_flavorTemplates_nodecays_R" + rad + "_k" + kappa + ".root", "READ" );
    //this file includes pythia8_norm_ + outFile_pt[j] and the individual flavor templates
    
    
    vector<TH1D*> syst_errs_envs = systematics( k, R ); // get systematics, quad sum the envelopes according to Isaac's procedure
    
    
    TH1D* unfs[njetbins]; // unfolded distributions go here
    TH1D* p6_p[njetbins]; // part-level distriubtions from pythia-6 (soon to be unmatched)
    TH1D* p8_p[njetbins]; // part-level distriubtions from pythia-8
    
    TH1D* uerr[njetbins]; // for systematic errors
    
    
    TH1D* p6_r[njetbins]; // ratio of pythia-6/data
    TH1D* p8_r[njetbins]; // ratio of pythia-8/data
    
    TH1D* er_r[njetbins]; // ratio of uncertainty histogram/unfolded data
    
    TH1D* systs[nsyst][njetbins]; // individual systematic histograms as percentages
    
    // get systematics histograms-->separated out into its own function
    getSystematics( systs ); // pass address of 2D array of TH1Ds so that no copy is made, the copy here is actually filled
    
    
    TCanvas *c = new TCanvas("c_unf", "canvas", 800, 500);
    c->SetLeftMargin(0.2);
    
    TLegend *leg1 = new TLegend(0.35, 0.7, .65, .75);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);

    
    TLegend *leg2 = new TLegend(0.3, 0.7, .7, .75);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    
    TLegend *leg3 = new TLegend(0.3, 0.7, .7, .75);
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
        
        p6_p[j] = (TH1D*) fp6->Get( Form( "jetQ_k" + kappa + "_part_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j] ) );

        p8_p[j] = (TH1D*) fp8->Get( "pythia8_norm_" + outFile_pt[j] );
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        unfs[j]->Scale( 1.0 / unfs[j]->Integral( 0, unfs[j]->GetNbinsX() + 1 ) );
        p6_p[j]->Scale( 1.0 / p6_p[j]->Integral( 0, p6_p[j]->GetNbinsX() + 1 ) );
        p8_p[j]->Scale( 1.0 / p8_p[j]->Integral( 0, p8_p[j]->GetNbinsX() + 1 ) );
        
        
        // clone normalized unfolded distributions (to be safe)
        uerr[j] = (TH1D*) unfs[j]->Clone( "unfold_errs_" + outFile_pt[j] );
        
        // calculate quadrature sum of errors from all systematic sources and set bin error of uerr[j] to that sum * bin content (quad sum is % error)
        //systErrs(systs, uerr, j);
        systEnvs(systs, uerr, j);
        
        
        // draw line at mean of unfolded distribution to show shift of mean
        TLine *line1 = new TLine( unfs[j]->GetMean(), 0, unfs[j]->GetMean(), 0.3 );
        line1->SetLineColor(kRed);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        
        
        
        
        unfs[j]->SetMarkerColor(kRed);
        p6_p[j]->SetLineColor(kBlue);
        p8_p[j]->SetLineColor(kBlack);
        
        syst_errs_envs[j]->SetFillColor(kRed - 10);
        uerr[j]->SetFillColor(kRed - 10);
        
        
        
        unfs[j]->SetMarkerStyle(kFullStar);
        unfs[j]->SetMarkerSize(2);
        
        //p6_p[j]->SetMarkerStyle(kFullStar);
        p6_p[j]->SetLineStyle(kSolid);
        //p6_p[j]->SetMarkerSize(2);
        p6_p[j]->SetLineWidth(3);
        
        
        p8_p[j]->SetLineStyle(kDashed);
        p8_p[j]->SetLineWidth(3);
        
        syst_errs_envs[j]->SetFillStyle(1001);
        uerr[j]->SetFillStyle(1001);
        
        
        // to take ratios of pythia-6(8)/unfolded data
        p6_r[j] = (TH1D*) p6_p[j]->Clone( "p6ratio_" + outFile_pt[j] );
        p8_r[j] = (TH1D*) p8_p[j]->Clone( "p8ratio_" + outFile_pt[j] );
        
        er_r[j] = (TH1D*) syst_errs_envs[j]->Clone( "syst_errs_envs_" + outFile_pt[j] );
        //er_r[j] = (TH1D*) uerr[j]->Clone( "uncertainty_ratio_" + outFile_pt[j] );
        
        p8_r[j]->SetLineStyle(kSolid);
        
        p6_r[j]->Divide( unfs[j] );
        p8_r[j]->Divide( unfs[j] );
        er_r[j]->Divide( unfs[j] );
        
        p6_r[j]->GetYaxis()->SetRangeUser(0.001, 1.999);
        p8_r[j]->GetYaxis()->SetRangeUser(0.001, 1.999);
        er_r[j]->GetYaxis()->SetRangeUser(0.001, 1.999);
        
        p6_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        er_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);

        
        if( j == 0 ){
            leg1->AddEntry( unfs[j], "STAR", "p");
        }
        else if( j == 1 ){
            leg2->AddEntry( p6_p[j], "PYTHIA-6 Perugia"/* 2012"*/, "l"); // only data with markers, all else with lines only
        }
        else if( j == 2 ){
            leg3->AddEntry( p8_p[j], "PYTHIA-8 Monash", "l");
        }
        
        
        syst_errs_envs[j]->SetStats(0);
        syst_errs_envs[j]->GetYaxis()->SetRangeUser(-0.099, 0.599);
        syst_errs_envs[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        
        
        syst_errs_envs[j]->GetYaxis()->SetTitleSize(0.1);
        syst_errs_envs[j]->GetYaxis()->SetTitleOffset(0.8);
        syst_errs_envs[j]->GetYaxis()->SetNdivisions(505);
        syst_errs_envs[j]->GetYaxis()->SetLabelSize(0.05);
        syst_errs_envs[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        //uerr[j]->SetXTitle( "Q_{jet}" );
        
        
        uerr[j]->SetStats(0);
        uerr[j]->GetYaxis()->SetRangeUser(-0.099, 0.599);
        uerr[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        
        
//        uerr[j]->GetYaxis()->SetTitleSize(0.1);
        uerr[j]->GetYaxis()->SetTitleFont(43);
        uerr[j]->GetYaxis()->SetTitleSize(22);
        uerr[j]->GetYaxis()->SetTitleOffset(1.2);
        uerr[j]->GetYaxis()->SetNdivisions(505);
//        uerr[j]->GetYaxis()->SetLabelSize(0.05);
        uerr[j]->GetYaxis()->SetLabelFont(43);
        uerr[j]->GetYaxis()->SetLabelSize(22);
        uerr[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        //uerr[j]->SetXTitle( "Q_{jet}" );
        
        
        cout << "jet bin: [" << jetPtLo[j] << ", " << jetPtHi[j] << "] mean = " << unfs[j]->GetMean() << "\n";
        
        
        p6_p[j]->Sumw2(0);
        p8_p[j]->Sumw2(0);
        
        //syst_errs_envs[j]->Draw("E3SAME");
        uerr[j]->Draw("E3SAME");
        unfs[j]->Draw("PSAME");
        p6_p[j]->Draw("CSAME");
        p8_p[j]->Draw("CSAME");
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
            //drawText( "", .28, .9, 20 );
            tprelim->DrawLatex(-2.3, 0.55, "STAR Preliminary");
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
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
        
//        er_r[j]->GetYaxis()->SetTitleSize(0.2);
        er_r[j]->GetYaxis()->SetTitleFont(43);
        er_r[j]->GetYaxis()->SetTitleSize(22);
        er_r[j]->GetYaxis()->SetTitleOffset(1.2);
        er_r[j]->GetYaxis()->SetNdivisions(505);
//        er_r[j]->GetYaxis()->SetLabelSize(0.1);
        er_r[j]->GetYaxis()->SetLabelFont(43);
        er_r[j]->GetYaxis()->SetLabelSize(22);
        er_r[j]->GetYaxis()->SetTitle( "MC / data" );
        
        //er_r[j]->GetXaxis()->SetTitleSize(0.2);
        er_r[j]->GetXaxis()->SetTitleFont(43);
        er_r[j]->GetXaxis()->SetTitleSize(22);
        er_r[j]->GetXaxis()->SetTitleOffset(1.1);
        //er_r[j]->GetXaxis()->SetNdivisions(505);
        er_r[j]->GetXaxis()->SetLabelSize(0.1);
        er_r[j]->GetXaxis()->SetLabelFont(43);
        er_r[j]->GetXaxis()->SetLabelSize(22);
        
        er_r[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        er_r[j]->Draw("E3SAME");
        p6_r[j]->Draw("E1SAME");
        p8_r[j]->Draw("E1SAME");
        
    }

    c->SaveAs( "./plots/dnp/TMPunfoldedAndParticle_jetCharge_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
    
}


void plot_dataDet_jetcharge(double k = 0.0, double R = 0.4){
    cout << "\nPlotting raw data and PYTHIA+GEANT on one canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    // file with both raw and unfolded data distributions for all jet pT bins
    TFile *fdat = new TFile( "../out/ppJP2data_unfolded_R" + rad + "_k" + kappa + ".root" , "READ");
    
    // file with pythia-6 distribution
    // unmatched pythia-6 --- only has part-level distribution...
    //TFile *fp6 = new TFile( "../out/p6_unmatched_R" + rad + "_k" + kappa + ".root" , "READ");
    // matched pythia-6
    TFile *fp6 = new TFile( "../out/p6_R" + rad + "_k" + kappa + ".root" , "READ");

    
    
    TH1D* data[njetbins]; // unfolded distributions go here
    TH1D* p6_d[njetbins]; // part-level distriubtions

    // for ratios
    TH1D* det_r[njetbins]; // ratio pythia-6/data
    
    
    
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

        
        
        data[j] = (TH1D*) fdat->Get( Form( "jetQ_k" + kappa + "_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j] ) );
        
        p6_d[j] = (TH1D*) fp6->Get( Form( "jetQ_k" + kappa + "_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j] ) );
        
        
        
        // normalize all to a probability distribution (integral = 1) including over/underflow bins
        data[j]->Scale( 1.0 / data[j]->Integral( 0, data[j]->GetNbinsX() + 1 ) );
        p6_d[j]->Scale( 1.0 / p6_d[j]->Integral( 0, p6_d[j]->GetNbinsX() + 1 ) );

        
        data[j]->SetMarkerColor(kBlack);
        p6_d[j]->SetMarkerColor(kBlue);
        
//        uerr[j]->SetFillColor(kViolet-3 - 1);
        
        
        
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
        det_r[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        //p8_r[j]->GetYaxis()->SetRangeUser(0,2);
        
        
        if( j == 0 ){
            leg1->AddEntry( data[j], "STAR (uncorrected)", "p");

        }
        else if( j == 1 ){
            leg2->AddEntry( p6_d[j], "PYTHIA-6+GEANT", "p"); // only data with markers, all else with lines only
        }
        
        data[j]->SetStats(0);
        data[j]->GetYaxis()->SetRangeUser(-0.099, 0.599);
        data[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        
        
//        data[j]->GetYaxis()->SetTitleSize(0.1);
        data[j]->GetYaxis()->SetTitleFont(43);
        data[j]->GetYaxis()->SetTitleSize(22);
        data[j]->GetYaxis()->SetTitleOffset(1.2);
        data[j]->GetYaxis()->SetNdivisions(505);
//        data[j]->GetYaxis()->SetLabelSize(0.05);
        data[j]->GetYaxis()->SetLabelFont(43);
        data[j]->GetYaxis()->SetLabelSize(22);
        data[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        //data[j]->GetXaxis()->SetTitleSize(0.1);
        //data[j]->GetXaxis()->SetTitleOffset(0.8);
        //data[j]->GetXaxis()->SetNdivisions(505);
        //data[j]->GetXaxis()->SetLabelSize(0.05);
        //data[j]->GetXaxis()->SetTitle( "Q_{jet}" );
        
        
        p6_d[j]->Sumw2(0);
        
        
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
            //drawText( "STAR Preliminary", .3, .9, 20 );
            tprelim->DrawLatex(-2.3, 0.55, "STAR Preliminary");
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .17, .81, 15 );
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
        
        
    }

    c->SaveAs( "./plots/dnp/TMPdataAndDetector_jetCharge_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
    
}


void plot_jetChargeResolution( double k = 0.0, double R = 0.4 ){
    cout << "\nPlotting jet charge resolution in pythia-6 / pythia-6 + geant on one canvas\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";
    
    
    // file with pythia-6 distribution
    TFile *fp6 = new TFile( "../out/p6_matched_R" + rad + "_k" + kappa + ".root" , "READ");
    
    TH1D* jc_res[njetbins];
    
    TCanvas *c;
    TLegend *leg;
    
    c = new TCanvas("c_res", "canvas", 800, 600);
    //c->SetLeftMargin(0.2);
    
    leg = new TLegend(0.55, 0.6, 0.85, .85);
    leg->SetBorderSize(0);
    
    for(int j = 0; j < njetbins; j++){
        // from hists.cxx, this histogram is ( Q_det - Q_part ) / Q_part for e.g. det jet pT = [20,25]
        jc_res[j] = (TH1D*) fp6->Get( "deltaQbypartQ_bydetpt_jetpt" + outFile_pt[j] );
        
        jc_res[j]->Scale( 1.0 / jc_res[j]->Integral() );
        jc_res[j]->SetStats(0);
        jc_res[j]->GetYaxis()->SetRangeUser(0.0, 0.7);
        jc_res[j]->GetYaxis()->SetTitleOffset(1.3);
        jc_res[j]->GetYaxis()->SetNdivisions(505);
        jc_res[j]->GetYaxis()->SetTitle( "1/N_{jet} dN/dr" );


        jc_res[j]->SetMarkerSize(2);
        
        jc_res[j]->GetXaxis()->SetRangeUser(-4.5, 4.5);
        jc_res[j]->SetXTitle( "r = ( Q^{det} - Q^{part} ) / Q^{part}" );
        
        if(j == 0){
            jc_res[j]->SetMarkerColor(kBlack);
            jc_res[j]->SetMarkerStyle(kOpenCircle);
            jc_res[j]->SetMarkerSize(2);
            
            jc_res[j]->Draw("P");
        }
        else if(j == 1){
            jc_res[j]->SetMarkerColor(kRed);
            jc_res[j]->SetMarkerStyle(kOpenSquare);

            jc_res[j]->Draw("PSAME");
        }
        else if(j == 2){
            jc_res[j]->SetMarkerColor(kBlue);
            jc_res[j]->SetMarkerStyle(kOpenCross);

            jc_res[j]->Draw("PSAME");
        }
        leg->AddEntry( jc_res[j], Form( "%1.0f GeV/c < p_{T}^{det.} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), "p" );
        
        
        drawText( "p+p #sqrt{s} = 200 GeV", .15, .81, 20 );
        drawText( "PYTHIA-6+GEANT", 0.15, 0.75, 20);
        drawText( Form( "anti-k_{T}, R = %0.1f, |#eta_{jet}| < 1-R", R ), .15, .69, 20 );
        
    }
    TLine *line = new TLine(jc_res[0]->GetMean() , 0, jc_res[0]->GetMean(), 0.6);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(2);
    
    
    leg->Draw();
    line->Draw();
    
    c->SaveAs( "./plots/dnp/TMPjetChargeRes_R" + rad + "_k" + kappa + ".pdf", "RECREATE");
    
}


void drawTemplates( double k = 0.0, double R = 0.4 ){
    // draw u,d,g templates and pythia-8 total distributions
    cout << "\ndrawing pythia-8 with pythia-8 templates\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n\n";
    
    
    // adding pythia-8 unmatched to this figure
    TFile *fp8 = new TFile( "../out/p8_flavorTemplates_nodecays_R" + rad + "_k" + kappa + ".root", "READ" );
    //this file includes pythia8_norm_ + outFile_pt[j] and the individual flavor templates
    // templates of the form "u_2025"

    const int nflavs = 3;
    TString flavs_names[nflavs] = {"d", "u", "g"};
    
    
    // normalized
    TH1D* htemps[njetbins][nflavs];
    TH1D* hpy8[njetbins];
    
    // non-normalized -- to calculate fractions
    TH1D* frac_temps[njetbins][nflavs]; // rel_abunds_u_2025
    TH1D* frac_py8[njetbins]; // pythia8_2025
    
    
    TCanvas *c = new TCanvas("c_temps", "canvas", 800, 500);
    c->SetLeftMargin(0.25);
    
    TLegend *leg1 = new TLegend(0.2, 0.75, 0.9, .85);
    leg1->SetNColumns(3);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.08);
    
    TLegend *leg2 = new TLegend(0.35, 0.75, 0.65, .85);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.08);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
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
        
        
        
        hpy8[j] = (TH1D*) fp8->Get( "pythia8_norm_" +outFile_pt[j] );
        frac_py8[j] = (TH1D*) fp8->Get( "pythia8_" +outFile_pt[j] );
        
        double py8_tot = frac_py8[j]->Integral( 0, frac_py8[j]->GetNbinsX() + 1 );
        
        
        hpy8[j]->SetStats(0);
        hpy8[j]->GetXaxis()->SetRangeUser( -4.5, 4.5 );
        hpy8[j]->GetYaxis()->SetRangeUser( 0.001, 0.599 );
        hpy8[j]->Sumw2(0);
        
//        hpy8[j]->GetYaxis()->SetTitleSize(0.1);
        hpy8[j]->GetYaxis()->SetTitleFont(43);
        hpy8[j]->GetYaxis()->SetTitleSize(22);
        hpy8[j]->GetYaxis()->SetTitleOffset(0.8);
//        hpy8[j]->GetYaxis()->SetLabelSize(0.05);
        hpy8[j]->GetYaxis()->SetLabelFont(43);
        hpy8[j]->GetYaxis()->SetLabelSize(22);
        hpy8[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );

        
        hpy8[j]->GetXaxis()->SetTitleFont(43);
        hpy8[j]->GetXaxis()->SetTitleSize(22);
        hpy8[j]->GetXaxis()->SetTitleOffset(1.2);
        //hpy8[j]->GetXaxis()->SetNdivisions(505);
        //hpy8[j]->GetXaxis()->SetLabelSize(0.05);
        hpy8[j]->GetXaxis()->SetLabelFont(43);
        hpy8[j]->GetXaxis()->SetLabelSize(22);
        hpy8[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        hpy8[j]->SetLineColor(kBlack);
        hpy8[j]->SetLineStyle(kSolid);
        hpy8[j]->SetLineWidth(2);
        
        if( j == 0 ){
            leg2->AddEntry( hpy8[j], "PYTHIA-8", "l" );
        }
        
        //hpy8[j]->Draw("CSAME");
        
        /*
        for(int i = 0; i < n_template_flavors; i++){
            if( template_flavor_names[i] == "u" || template_flavor_names[i] == "d" || template_flavor_names[i] == "g" ){continue;}
            else if( template_flavor_names[i] == "s" ){ // is is first flavor after d,u
                h_other[j] = (TH1D*) fp8->Get( "rel_abunds_" + template_flavor_names[i] + "_" + outFile_pt[j] );
            }
            else{ // if flavor is not u,d,g and not s (first non-u,d,g flavor) add non-normalized template to h_other
                h_other[j]->Add( (TH1D*) fp8->Get( "rel_abunds_" + template_flavor_names[i] + "_" + outFile_pt[j] ) );
            }
        }
        */
        
        
        for(int k = 0; k < nflavs; k++){
            htemps[j][k] = (TH1D*) fp8->Get( flavs_names[k] + "_" + outFile_pt[j] );
            
            frac_temps[j][k] = (TH1D*) fp8->Get( "rel_abunds_" + flavs_names[k] + "_" + outFile_pt[j] );
            
        
            double frac = frac_temps[j][k]->Integral( 0, frac_temps[j][k]->GetNbinsX() + 1 ) / py8_tot;
            
            htemps[j][k]->SetStats(0);
            htemps[j][k]->GetXaxis()->SetRangeUser( -4.5, 4.5 );
            htemps[j][k]->GetYaxis()->SetRangeUser( 0.001, 0.599 );
            htemps[j][k]->Sumw2(0);
            
            
            
//            htemps[j][k]->GetYaxis()->SetTitleSize(0.1);
            htemps[j][k]->GetYaxis()->SetTitleFont(43);
            htemps[j][k]->GetYaxis()->SetTitleSize(22);
            htemps[j][k]->GetYaxis()->SetTitleOffset(1.3);
//            htemps[j][k]->GetYaxis()->SetLabelSize(0.05);
            htemps[j][k]->GetYaxis()->SetLabelFont(43);
            htemps[j][k]->GetYaxis()->SetLabelSize(22);
            htemps[j][k]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
            
            //htemps[j][k]->GetXaxis()->SetTitleSize(0.1);
            htemps[j][k]->GetXaxis()->SetTitleFont(43);
            htemps[j][k]->GetXaxis()->SetTitleSize(22);
            htemps[j][k]->GetXaxis()->SetTitleOffset(1.1);
            //frac_temps[j][k]->GetXaxis()->SetNdivisions(505);
            //htemps[j][k]->GetXaxis()->SetLabelSize(0.05);
            htemps[j][k]->GetXaxis()->SetLabelFont(43);
            htemps[j][k]->GetXaxis()->SetLabelSize(22);
            htemps[j][k]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
            
            //htemps[j][k]->Scale( frac );
            
            htemps[j][k]->SetLineColor( mark_colors[k] );
            htemps[j][k]->SetLineStyle( line_styles[k] );
            htemps[j][k]->SetLineWidth(2);
            
            htemps[j][k]->Sumw2(0);
            
            if( j == 0 ){
                drawText( "PYTHIA-8 Monash", .35, .93, 20 );
                
                if( k == 0 ){
                    leg1->AddEntry( htemps[j][k], "d", "l" );
                }
                else if( k == 1 ){
                    leg1->AddEntry( htemps[j][k], "u", "l" );
                }
                else if( k == 2 ){
                    leg1->AddEntry( htemps[j][k], "g", "l" );
                }
            }
            if( k == 0 ){
                htemps[j][k]->Draw("C");
            }
            else{
                htemps[j][k]->Draw("CSAME");
            }
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
            //leg2->Draw();
            drawText( Form( "%1.0f GeV/c < p_{T}^{jet} < %1.0f GeV/c", jetPtLo[j], jetPtHi[j] ), .22, .87, 15 );
            drawText( "p+p #sqrt{s} = 200 GeV", .25, .93, 15 );
        }
        
        
        TLine *line = new TLine( 0, 0.001, 0, 0.599 );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        /*
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.05); tprelim->SetTextColor(kRed);
        tprelim->DrawLatexNDC(0.31, 0.87, "STAR Simulation");
        tprelim->Draw();
        */
        
        //line->Draw();

    }
    c->SaveAs( "./plots/dnp/p8templates_udg_R04_k00.pdf", "RECREATE" );
}




void templates_weighted( double k = 0.0, double R = 0.4 ){
    // perform tfractionfitter fit to pythia-8 total distribution using pythia-8 templates
    cout << "\nfitting pythia-8 with pythia-8 templates\n";

    TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n\n";
    
    
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Minimize");
    //ROOT::Math::MinimizerOptions::SetDefaultErrorDef( 1. );
    
    //ROOT::Math::MinimizerOptions::SetDefaultStrategy( 0 );
    
    //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );
    
    ROOT::Math::MinimizerOptions::PrintDefault();
    
    
    
    // adding pythia-8 unmatched to this figure
    TFile *fp8 = new TFile( "../out/p8_flavorTemplates_nodecays_R" + rad + "_k" + kappa + ".root", "READ" );
    //this file includes pythia8_norm_ + outFile_pt[j] and the individual flavor templates
    // templates of the form "u_2025"

    const int nflavs = 3;
    TString flavs_names[nflavs] = {"d", "u", "g"};
    
    
    // normalized
    TH1D* htemps[njetbins][nflavs];
    TH1D* hpy8[njetbins];
    
    // non-normalized -- to calculate fractions
    TH1D* frac_temps[njetbins][nflavs]; // rel_abunds_u_2025
    TH1D* frac_py8[njetbins]; // pythia8_2025
    
    
    // for tfractionfitter
    TObjArray *h_fit[njetbins];
    TFractionFitter *fit[njetbins];
    TH1D* h_other[njetbins];
    
    TH1D* result_tff[njetbins];
    
    
    TCanvas *c = new TCanvas("c_fits", "canvas", 800, 500);
    c->SetLeftMargin(0.25);
    
    TLegend *leg1 = new TLegend(0.2, 0.75, 0.9, .85);
    leg1->SetNColumns(3);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.08);
    
    TLegend *leg2 = new TLegend(0.35, 0.75, 0.65, .85);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.08);
    
    c->Divide(njetbins, 1, 0, 0);
    
    
    for(int j = 0; j < njetbins; j++){
        c->cd( j+1 );
        TString padup_name = Form( "pfits_%i", j );
        TPad* padup = new TPad(padup_name, padup_name, 0, 0, 1, 1);
        
        padup->SetRightMargin(0);
        padup->SetTopMargin(0);
        padup->SetBottomMargin(0.2);
        padup->SetLeftMargin(0);
        if( j == 0 ){padup->SetLeftMargin(0.2);}
        
        padup->Draw();
        padup->cd();
        
        h_fit[j] = new TObjArray(nflavs+1); // nflavs: u,d,g; template for all other flavors combined
        
        
        hpy8[j] = (TH1D*) fp8->Get( "pythia8_norm_" +outFile_pt[j] );
        frac_py8[j] = (TH1D*) fp8->Get( "pythia8_" +outFile_pt[j] );
        
        double py8_tot = frac_py8[j]->Integral( 0, frac_py8[j]->GetNbinsX() + 1 );
        
        hpy8[j]->SetStats(0);
        hpy8[j]->Sumw2(0);
        hpy8[j]->GetXaxis()->SetRangeUser( -4.5, 4.5 );
        hpy8[j]->GetYaxis()->SetRangeUser( 0.001, 0.599 );

        
//        hpy8[j]->GetYaxis()->SetTitleSize(0.1);
        hpy8[j]->GetYaxis()->SetTitleFont(43);
        hpy8[j]->GetYaxis()->SetTitleSize(22);
        hpy8[j]->GetYaxis()->SetTitleOffset(1.3);
//        hpy8[j]->GetYaxis()->SetLabelSize(0.05);
        hpy8[j]->GetYaxis()->SetLabelFont(43);
        hpy8[j]->GetYaxis()->SetLabelSize(22);
        hpy8[j]->GetYaxis()->SetTitle( "1/N_{jet} dN_{jet}/dQ_{jet} [1/e]" );
        
        /*
        hpy8[j]->GetXaxis()->SetTitleFont(63);
        hpy8[j]->GetXaxis()->SetTitleSize(20); // labels will be 14 pixels
        hpy8[j]->GetXaxis()->SetLabelFont(63);
        hpy8[j]->GetXaxis()->SetLabelSize(16); // labels will be 14 pixels
        */
        
        //hpy8[j]->GetXaxis()->SetTitleSize(0.1);
        hpy8[j]->GetXaxis()->SetTitleFont(43);
        hpy8[j]->GetXaxis()->SetTitleSize(22);
        hpy8[j]->GetXaxis()->SetTitleOffset(1.1);
        //hpy8[j]->GetXaxis()->SetNdivisions(505);
        //hpy8[j]->GetXaxis()->SetLabelSize(0.05);
        hpy8[j]->GetXaxis()->SetLabelFont(43);
        hpy8[j]->GetXaxis()->SetLabelSize(22);
        hpy8[j]->GetXaxis()->SetTitle( Form( "Q^{jet}_{#kappa = %1.1f} [e]", k ) );
        
        hpy8[j]->SetLineColor(kBlack);
        hpy8[j]->SetLineStyle(kSolid);
        hpy8[j]->SetLineWidth(2);
        
        if( j == 0 ){
            leg2->AddEntry( hpy8[j], "PYTHIA-8", "l" );
        }
        
        hpy8[j]->Draw("CSAME");
        
        // loop over all 11 flavors with templates (no t, tbar)
        for(int i = 0; i < n_template_flavors; i++){
            if( template_flavor_names[i] == "u" || template_flavor_names[i] == "d" || template_flavor_names[i] == "g" ){continue;} // skip u,d,g for "other" template, these will have their own templates
            else if( template_flavor_names[i] == "s" ){ // is is first flavor after d,u
                cout << "s template used as OTHER\n";
                h_other[j] = (TH1D*) fp8->Get( "rel_abunds_" + template_flavor_names[i] + "_" + outFile_pt[j] );
            }
            else{ // if flavor is not u,d,g and not s (first non-u,d,g flavor) add non-normalized template to h_other
                cout << template_flavor_names[i] << " template added to OTHER\n";
                h_other[j]->Add( (TH1D*) fp8->Get( "rel_abunds_" + template_flavor_names[i] + "_" + outFile_pt[j] ) );
            }
        }
        
        
        // loop over just u,d,g to get their individual templates
        for(int k = 0; k < nflavs; k++){
            cout << "BEFORE FITTING...\n";
            cout << flavs_names[k] << "\n";
            htemps[j][k] = (TH1D*) fp8->Get( flavs_names[k] + "_" + outFile_pt[j] );
            
            frac_temps[j][k] = (TH1D*) fp8->Get( "rel_abunds_" + flavs_names[k] + "_" + outFile_pt[j] );

            cout << frac_temps[j][k]->Integral() / frac_py8[j]->Integral() << "\n";
            
            h_fit[j]->Add( frac_temps[j][k] );
        }
        h_fit[j]->Add( h_other[j] ); // add "other" after u,d,g
        
        fit[j] = new TFractionFitter( hpy8[j], h_fit[j] /*, "V"*/ );
        
        for(int k = 0; k < nflavs; k++){
            fit[j]->Constrain( k, 0.0, 1.0 );
                              //( 1.0 - 0.5 ) * ( frac_temps[j][k]->Integral() / frac_py8[j]->Integral() ),
                              //( 1.0 + 0.5 ) * ( frac_temps[j][k]->Integral() / frac_py8[j]->Integral() ) );
        }
        
        fit[j]->Constrain( nflavs, ( 1.0 - 0.001 ) * ( h_other[j]->Integral() / frac_py8[j]->Integral() ), ( 1.0 + 0.001 ) * ( h_other[j]->Integral() / frac_py8[j]->Integral() ) ); // fix "other" fraction to known pythia-8 value
        
        int status = fit[j]->Fit();
        cout << "TFractionFitter fit status: " << status << "\n";
        
        
        result_tff[j] = (TH1D*) fit[j]->GetPlot();
        result_tff[j]->SetName( "frac_fit_" + outFile_pt[j] );
        
        
        // end fitting procedure, now can plot
        
        for(int k = 0; k < nflavs; k++){ // loop over actual flavors (u,d,g)
        
            //double frac = frac_temps[j][k]->Integral( 0, frac_temps[j][k]->GetNbinsX() + 1 ) / py8_tot;
            double frac = 1.0;// set to parameter from the fit result
            double err = 1.0;
            
            fit[j]->GetResult( k, frac, err );
            
            if(k == 1){
                cout << "\npythia-8 flavor u fraction: " << frac_temps[j][k]->Integral() / frac_py8[j]->Integral() << "\n";
                cout << "fit result for u fraction: " << frac << "\n";
            }
            else if(k == 0){
                cout << "\npythia-8 flavor d fraction: " << frac_temps[j][k]->Integral() / frac_py8[j]->Integral() << "\n";
                cout << "fit result for d fraction: " << frac << "\n";
            }
            else if( k == 2){
                cout << "\npythia-8 flavor g fraction: " << frac_temps[j][k]->Integral() / frac_py8[j]->Integral() << "\n";
                cout << "fit result for g fraction: " << frac << "\n";
                
                
                double frac_other = 0.; double err_other = 0.;
                fit[j]->GetResult(nflavs, frac_other, err_other);
                
                cout << "\npythia-8 flavor other fraction: " << h_other[j]->Integral() / frac_py8[j]->Integral() << "\n";
                cout << "fit result for other fraction: " << frac_other << "\n";
            }
            
            htemps[j][k]->Scale( frac );
            
            htemps[j][k]->SetLineColor( mark_colors[k] );
            htemps[j][k]->SetLineStyle( line_styles[k] );
            htemps[j][k]->SetLineWidth(2);
            
            htemps[j][k]->Sumw2(0);
            
            if( j == 0 ){
                drawText( "PYTHIA-8 Monash", .35, .93, 20 );
                
                if( k == 0 ){
                    leg1->AddEntry( htemps[j][k], "d", "l" );
                }
                else if( k == 1 ){
                    leg1->AddEntry( htemps[j][k], "u", "l" );
                }
                else if( k == 2 ){
                    leg1->AddEntry( htemps[j][k], "g", "l" );
                }
            }
            
            htemps[j][k]->Draw("CSAME");
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
        TLine *line = new TLine( 0, 0.001, 0, 0.599 );
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(1);
        
        /*
        TLatex *tprelim = new TLatex();
        tprelim->SetTextSize(0.05); tprelim->SetTextColor(kRed);
        tprelim->DrawLatexNDC(0.31, 0.87, "STAR Simulation");
        tprelim->Draw();
        */
        
        //line->Draw();

    }
    c->SaveAs( "./plots/dnp/FitTemplatesToP8_udg_R04_k00.pdf", "RECREATE" );
}


void draw_plotsForDNP( double k = 0.0, double R = 0.4 ){
    plot_unfPart_jetcharge( k, R );
    plot_dataDet_jetcharge( k, R );
    plot_jetChargeResolution( k, R );
    drawTemplates( k, R );
    templates_weighted( k, R );
}


