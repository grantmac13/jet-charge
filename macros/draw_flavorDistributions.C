// 8/13/22 Grant McNamara
// draw distributions from pythia-8, not actual templates
// draw_flavorDistributions.C

#include <string>
#include <iostream>
#include "math.h"
#include "Headers/plot.h"
//#include "Headers/JCparameters.hh"


// plot pythia-6, pythia-6+geant, data and data unfolded for detector effects
void draw_flavorDistributions(double k = 0.0, TString decay = "decayed", double R = 0.4){
	// ../out/data/hists/ppDataJP2_R0	_R04.root
	//			ppDataJP2_R06.root

	TString rad = Form("%i%i", (int) R, (int) (R*10) );
        TString kappa = Form("%i%i", (int) k, (int) (k*10) );
        cout << "converting double to string: radius = " << rad << "\n";
        cout << "converting double to string: kappa = " << kappa << "\n";

    const int nflavors = 13;
    TString flavor_names[nflavors] = {"d",    "u",    "s",
                                      "c",    "b",    "t",  "g",
                                      "dbar", "ubar", "sbar",
                                      "cbar", "bbar", "tbar"};
    
    
    int mark_colors[(int) (nflavors+1)/2] = { 600, 628, 423,
                                              891, 874, 1, 839};


	const int njetbins = 3;

	double jetPtLo[njetbins] = {/*15.0,*/ 20.0, 25.0, 30.0};
	double jetPtHi[njetbins] = {/*20.0,*/ 25.0, 30.0, 40.0};
	// loop over all jet pt ranges



	TFile* fin = new TFile( "~/for_grant/out/p8_" + decay + "_R" + rad + "_k" + kappa + ".root" , "READ");
	

	TH2D* jetQ[nflavors];
    TH1D* jetQ_proj[nflavors][njetbins];

    
    for(int i = 0; i < nflavors; i++){
        jetQ[i] = (TH2D*) fin->Get(  flavor_names[i] + "_qvpt" );
    }
    

    TCanvas *c[njetbins];
	TLegend *leg[njetbins];


	for(int j = 0; j < njetbins; j++){
		c[j] = new TCanvas(Form("c_%i", j), "canvas", 600, 600);
		leg[j] = new TLegend(0.18, 0.65, .45, .85);
		leg[j]->SetBorderSize(0);
//        leg[j]->SetNColumns(3);

        c[j]->SetLeftMargin(0.15);
        
        // jetpt axis includes 15-20 bin, so j+2 = [2, 4] corresponding to 20-25, 25-30, 30-40
        for(int i = 0; i < nflavors; i++){
            jetQ_proj[i][j] = (TH1D*) jetQ[i]->ProjectionX( Form( flavor_names[i] + "jetQ_pt%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ), j+2, j+2 );



            // normalize histograms i guess?
            jetQ_proj[i][j]->Scale( 1.0 / ( jetQ_proj[i][j]->Integral() ) );


            jetQ_proj[i][j]->SetStats(0);
            jetQ_proj[i][j]->GetYaxis()->SetRangeUser(0.0, 0.5);
/*		if(kappa != "00"){
			jetQ_unfold[j]->GetYaxis()->SetRangeUser(0.0, 0.3);
		}
*/


            c[j]->SetLeftMargin(0.15);

            jetQ_proj[i][j]->SetLineColor( mark_colors[i%7] ); // remainder when i / 7 should go from 0 to 6 , i%7 + 1 ->[1, 7]
            jetQ_proj[i][j]->SetLineWidth(2);
            if(i < 7){ // 0-6 should be quarks and gluon
                jetQ_proj[i][j]->SetLineStyle(1); // solid
            }
            else{ // 7-12 should be antiquarks
                jetQ_proj[i][j]->SetLineStyle(9); // long dashes
            }
            
            if( i%7 == 0 || i%7 == 1 || i%7 == 6 ){ // if u, d, g, ubar, dbar
                TString leg_name = "";// = flavor_names[i];
                if(6 < i && i < 12){leg_name = "#bar";}
                leg_name += flavor_names[i%7];
            
                leg[j]->AddEntry(jetQ_proj[i][j], leg_name, "l");
                //leg[j]->AddEntry((TObject*)0, Form("mean: %1.3f", jetQ_proj[i][j]->GetMean()), "");

                jetQ_proj[i][j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV" /*	R = %0.1f*/, jetPtLo[j], jetPtHi[j] /*, R */) );
                jetQ_proj[i][j]->SetXTitle( Form( "Q^{jet}_{#kappa = %1.1f}" , k ) );

                jetQ_proj[i][j]->Draw("HIST C SAME"); // option to use: either "C" or "L" ?
            }
        }
		leg[j]->Draw();

//		drawText("p+p #sqrt{s_{NN}} = 200 GeV}", .25, .85, 15);
		drawText( Form( "R = %0.1f", R ), .65, .7, 30 );

		// Save plot here, after points have been plotted
		c[j]->SaveAs( Form( "./plots/flavorDistributions_jetQ%1.0f%1.0f_R" + rad + "_k" + kappa + ".png", jetPtLo[j], jetPtHi[j] ), "RECREATE" );
	}
}
