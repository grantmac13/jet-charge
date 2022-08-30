// 2/4/2022 Grant McNamara
// take the hadd-ed output files from running data, pythia6, pythia8 over hists.cxx (which makes kappa = 0.0, 0.5 jet charge hists)
// and plot the jet charge histograms for a given kappa (default k = 0.5)

#include <string>
#include <iostream>
#include "math.h"
#include "Headers/plot.h"
//#include "Headers/JCparameters.hh"


// plot pythia-6, pythia-6+geant, data and data unfolded for detector effects
void draw_unfoldedFits(double k = 0.0, TString decay = "decayed", double R = 0.4){
	// ../out/data/hists/ppDataJP2_R0	_R04.root
	//			ppDataJP2_R06.root

	TString rad = Form("%i%i", (int) R, (int) (R*10) );
    TString kappa = Form("%i%i", (int) k, (int) (k*10) );
    cout << "converting double to string: radius = " << rad << "\n";
    cout << "converting double to string: kappa = " << kappa << "\n";

    const int nflavors = 6;
    TString flavor_names[nflavors] = {"d", "u", "g", "dbar", "ubar", "other"};
    
    
    //                             d,   u,   g,  other
    int mark_colors[nflavors/2 + 1] = { 600, 628, 839,   1};
    


	const int njetbins = 3;

	double jetPtLo[njetbins] = {20.0, 25.0, 30.0};
	double jetPtHi[njetbins] = {25.0, 30.0, 40.0};
	// loop over all jet pt ranges


    // file containing template fits by flavor
	TFile* ftemp = new TFile( "~/ppJCandBF/out/template/R" + rad + "_k" + kappa + ".root" , "READ");

    // file containing unfolded result with fit
    TFile* funf = new TFile( "~/ppJCandBF/out/fit_templates/R" + rad + "_k" + kappa + ".root" , "READ");


    TH1D* temps[nflavors][njetbins]; // naming convention for template histograms: Form( flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] )
    TF1* temp_fits[nflavors][njetbins]; // naming convention for the template fit functionsTString histfit_name = Form( "fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] );
    
    
    TH1D* unf[njetbins]; // unfolded distribution, with fit named according to:  Form( "fit_unfold_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] )
    
    TF1* unf_fit[njetbins]; // fits associated to unfolded distribution from fit_templates.cxx with naming convention Form( "fit_unfold_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] )
    
    double params[njetbins][nflavors]; // one remaining parameter not contained within the temp_fits functions, just the relative weight between flavors in the sum of gaussians fit
    
    
    
    TF1* plot_fits[nflavors][njetbins];
    

    TCanvas *c[njetbins];
	TLegend *leg[njetbins];


	for(int j = 0; j < njetbins; j++){
		c[j] = new TCanvas(Form("c_%i", j), "canvas", 600, 600);
		leg[j] = new TLegend(0.18, 0.65, .45, .85);
		leg[j]->SetBorderSize(0);
        leg[j]->SetNColumns(3);

        // jetpt axis includes 15-20 bin, so j+2 = [2, 4] corresponding to 20-25, 25-30, 30-40
        
        TString plot_fit_name = ""; // for largest params[j][i1] value, start with corresponding i1, plot by itself, next largest, add largest i*temp_fits[i1][j] + i2*temp_fits[i2][j], etc
        
        
        for(int i = 0; i < nflavors; i++){
            temps[i][j] = (TH1D*) ftemp->Get( Form( flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
            
            temp_fits[i][j] = (TF1*) temps[i][j]->GetFunction( Form( "fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );


            unf[j] = (TH1D*) funf->Get( Form( "unfold_nomx%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
            
            TString unf_fit_name = Form( "fit_unfold_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] );
            
            unf_fit[j] = (TF1*) unf[j]->GetFunction( unf_fit_name );
            
            params[j][i] = unf_fit[j]->GetParameter( i );


            /*
            // what i want to do here plotting:
            first, plot d weight in sum fit of unfolded distribution (params[j][0]) times the d template from template.cxx stored in temps[0][j]
            then, for u: plot params[j][0]*temps[0][j] + params[j][1]*temps[1][j]
            then, for g: plot params[j][0]*temps[0][j] + params[j][1]*temps[1][j] + params[j][2]*temps[2][j]
            then, for other: plot params[j][0]*temps[0][j] + params[j][1]*temps[1][j] + params[j][2]*temps[2][j] + params[j][3]*temps[3][j]
        
             but plot in reverse order so that the smallest (since all functions are positive definite), first distribution for d is on top and all are visible simulatneously
            */
            
            
            // Form( "fit_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] )
            plot_fit_name += Form( "%f * fit_" + flavor_names[i] + "_%1.0f%1.0f(x)" , (float) params[j][i], jetPtLo[j], jetPtHi[j] );
            
            
            plot_fits[i][j] = new TF1( Form( "plot_" + flavor_names[i] + "_%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ), plot_fit_name );
            
            if(i != nflavors){plot_fit_name += " + ";}

            
            if(i%3 != 2){plot_fits[i][j]->SetLineColor( mark_colors[i%3] );}
            else{plot_fits[i][j]->SetLineColor( mark_colors[i/3 + 2] );}
            
            if( !(i/3) ){ // i is 0-5 so i/3 if i is int is 0 or 1
                plot_fits[i][j]->SetLineStyle(1); // solid
            }
            else{
                plot_fits[i][j]->SetLineStyle(9); // long dashes
            }
            
            
            plot_fits[i][j]->Draw("same");
            
            
        } // flavor loop
        
    } // jet pt loop
            
/*

            jetQ_proj[i][j]->SetStats(0);
            jetQ_proj[i][j]->GetYaxis()->SetRangeUser(0.0, 0.5);
//
//            if(kappa != "00"){
//                jetQ_unfold[j]->GetYaxis()->SetRangeUser(0.0, 0.3);
//            }
//


            c[j]->SetLeftMargin(0.15);

            jetQ_proj[i][j]->SetLineColor( mark_colors[i%7] ); // remainder when i / 7 should go from 0 to 6 , i%7 + 1 ->[1, 7]
            jetQ_proj[i][j]->SetLineWidth(2);
            if(i < 7){ // 0-6 should be quarks and gluon
                jetQ_proj[i][j]->SetLineStyle(1); // solid
            }
            else{ // 7-12 should be antiquarks
                jetQ_proj[i][j]->SetLineStyle(9); // long dashes
            }
            
            if( i%7 == 0 || i%7 == 1 || i%7 == 6 ){ //
                TString leg_name = "";// = flavor_names[i];
                if(6 < i && i < 12){leg_name = "#bar";}
                leg_name += flavor_names[i%7];
            
                leg[j]->AddEntry(jetQ_proj[i][j], leg_name, "l");
                //leg[j]->AddEntry((TObject*)0, Form("mean: %1.3f", jetQ_proj[i][j]->GetMean()), "");

                jetQ_proj[i][j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV", jetPtLo[j], jetPtHi[j] ) );
                jetQ_proj[i][j]->SetXTitle( Form( "Q^{jet}_{#kappa = %1.1f}" , k ) );

                jetQ_proj[i][j]->Draw("HIST C SAME"); // option to use: either "C" or "L" ?
            }
        }
		leg[j]->Draw();

//		drawText("p+p #sqrt{s_{NN}} = 200 GeV}", .25, .85, 15);
		drawText( Form( "R = %0.1f", R ), .65, .7, 30 );

		// Save plot here, after points have been plotted
		c[j]->SaveAs( Form( "./plots/Template_jetQ%1.0f%1.0f_R" + rad + "_k" + kappa + ".png", jetPtLo[j], jetPtHi[j] ), "RECREATE" );
	}
*/
    

}



