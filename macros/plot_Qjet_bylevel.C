// 2/4/2022 Grant McNamara
// take the hadd-ed output files from running data, pythia6, pythia8 over hists.cxx (which makes kappa = 0.0, 0.5 jet charge hists)
// and plot the jet charge histograms for a given kappa (default k = 0.5)

#include <string>
#include <iostream>
#include "math.h"
#include "Headers/plot.h"


// plot 2 of 4 jet charge distributions at a time: either particle level and corrected(unfolded) or detector level and raw data
void plot_Qjet_bylevel(double R = 0.4, double k = 0.0, bool level = 0){
	// level: if 0->data, det-level; 1->part-level, unfolded with systematics

	// ../out/data/hists/ppDataJP2_R0	_R04.root
	//			ppDataJP2_R06.root

	TString rad = Form("%i%i", (int) R, (int) (R*10) );
        TString kappa = Form("%i%i", (int) k, (int) (k*10) );
        cout << "converting double to string: radius = " << rad << "\n";
        cout << "converting double to string: kappa = " << kappa << "\n";


	const int nsyst = 9;
	string syst_names[nsyst] = {"IP2x", "IP6x", "TSx", "TUx", "HC50x",
				     "DSx", "GSx", "h7smear1D_", "p8smear1D_"};


	const int njetbins = 3;

	double jetPtLo[njetbins] = {/*15.0,*/ 20.0, 25.0, 30.0};
	double jetPtHi[njetbins] = {/*20.0,*/ 25.0, 30.0, 40.0};
	// loop over all jet pt ranges

	TString lev[2] = {"data", "unf"};


	TString hist_name = "jetQ_k" + kappa; //hists for p8, data
	TString hist_part = "jetQ_k" + kappa + "_part"; // for p6
	// then add jetpt%1.0f_%1.0f, jetPtLo[], jetPtHi[] for the corresponding jet pt range
	// 

	TFile* data; TFile* p6; TFile* unfold;

	TH1D* d_jets; TH1D* p6_jets; TH1D* p6_jets_p;
	
	
	data = new TFile("../out/data/hists/ppDataJP2_R" + rad + "_k" + kappa + ".root", "READ");
	unfold = new TFile("../out/unfold/unfolded_R" + rad + "_k" + kappa + ".root", "READ");
	p6 = new TFile("../out/sim/hists/p6_R" + rad + "_k" + kappa + "_June7.root", "READ");

	d_jets = (TH1D*) data->Get("jet_pt_hist_k" + kappa);
	p6_jets = (TH1D*) p6->Get("jet_pt_hist_k" + kappa);
	p6_jets_p = (TH1D*) p6->Get("jet_pt_hist_part_k" + kappa);

	TH1D* jetQ_data[njetbins];
	TH1D* jetQ_p6_p[njetbins];
	TH1D* jetQ_p6_d[njetbins];
	TH1D* jetQ_unfold[njetbins];
	TH1D* jetQ_err[njetbins];

	TH1D* data_2Ds[njetbins];

	TH2D* jetQ_data_2D = (TH2D*) data->Get( "QvPt_d" );

	TH1D* systs[njetbins][nsyst];


	TCanvas *c[njetbins];
	TLegend *leg[njetbins];


	for(int j = 0; j < njetbins; j++){
		c[j] = new TCanvas(Form("c_%i", j), "canvas", 600, 600);
		leg[j] = new TLegend(0.18, 0.65, .45, .85);
		leg[j]->SetBorderSize(0);

		// only take histogram corresponding to kappa that is passed to function
		// jetQ_k05_jetpt10_15

		data_2Ds[j] = (TH1D*) jetQ_data_2D->ProjectionX( Form("q_2D_jetpt%1.0f%1.0f", jetPtLo[j], jetPtHi[j]), j+1, j+1, "e" );
		cout << "PROJECTION COMPLETE\n";

		jetQ_data[j] = (TH1D*) data->Get(hist_name + Form("_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j]));
		jetQ_p6_d[j] = (TH1D*) p6->Get(hist_name + Form("_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j]));
		jetQ_p6_p[j] = (TH1D*) p6->Get(hist_part + Form("_jetpt%1.0f_%1.0f", jetPtLo[j], jetPtHi[j]));
//		jetQ_unfold[j] = (TH1D*) unfold->Get( Form( "unfold_objx%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
		jetQ_unfold[j] = (TH1D*) unfold->Get( Form( "unfold_nomx%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );
		jetQ_err[j] = (TH1D*) jetQ_unfold[j]->Clone( Form( "unfold_errsx%1.0f%1.0f", jetPtLo[j], jetPtHi[j] ) );

		for(int k = 0; k < nsyst; k++){
			systs[j][k] = (TH1D*) unfold->Get( Form( ("unfold_" + syst_names[k] + "%1.0f%1.0f").c_str(), jetPtLo[j], jetPtHi[j]  ) );
		}
		// loaded all systematics hists, so now I can loop over n-th bin for all systs and quad sum the error

		for(int i = 1; i < (int) systs[j][0]->GetNbinsX(); i++){
			double err_quadsum = 0.0;
			for(int k = 0; k < nsyst; k++){
				double syst_err = systs[j][k]->GetBinContent(i);
				err_quadsum += syst_err*syst_err; // total error is sqrt(sum of individual errors ^2)
			}
			jetQ_err[j]->SetBinError( i, sqrt( err_quadsum ) * jetQ_unfold[j]->GetBinContent(i) );
			// err_quadsum is a percentage of the bin content, so to calculate absolute error have to multiply sqrt(err_quadsum) * bincontent
		}


		// scale histograms by number of jets in relevant jet pt range (maybe unnecessary, mean doesn't depend on scaling?)
		jetQ_data[j]->Scale(1.0/(d_jets->Integral( (int) jetPtLo[j] + 1, (int) jetPtHi[j] )  ) );
		jetQ_p6_d[j]->Scale(1.0/(p6_jets->Integral( (int) jetPtLo[j] + 1, (int) jetPtHi[j] )  ) );
		
		// need to change the normalization here
		jetQ_p6_p[j]->Scale(1.0/(p6_jets_p->Integral( (int) jetPtLo[j] + 1, (int) jetPtHi[j] )  ) );

		cout << "integral of projection: " << data_2Ds[j]->Integral() << "\n";
		data_2Ds[j]->Scale( 1.0 / data_2Ds[j]->Integral() );
		cout << "NORMALIZATION COMPLETE\n";


		c[j]->SetLeftMargin(0.15);


		if(level){
			jetQ_unfold[j]->SetStats(0);
			jetQ_unfold[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
			if(kappa != "00"){
				jetQ_unfold[j]->GetYaxis()->SetRangeUser(0.0, 0.3);
			}
			jetQ_err[j]->SetStats(0);
			jetQ_err[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
			if(kappa != "00"){
				jetQ_err[j]->GetYaxis()->SetRangeUser(0.0, 0.3);
			}

	
			jetQ_unfold[j]->SetMarkerColor(kViolet-3);
			jetQ_err[j]->SetFillColorAlpha(kViolet-3, 0.1);
			jetQ_p6_p[j]->SetMarkerColor(kGreen+2);


			jetQ_unfold[j]->SetMarkerStyle(kFullCircle);
			jetQ_p6_p[j]->SetMarkerStyle(kFullCircle);
			jetQ_err[j]->SetFillStyle(3444);
//			jetQ_err[j]->SetFillStyle(4020);

			leg[j]->AddEntry(jetQ_p6_p[j], "Pythia 6 particle level", "p");
			leg[j]->AddEntry((TObject*)0, Form("mean: %1.3f", jetQ_p6_p[j]->GetMean()), "");

			leg[j]->AddEntry(jetQ_unfold[j], "unfolded", "p");
			leg[j]->AddEntry((TObject*)0, Form("mean: %1.3f", jetQ_unfold[j]->GetMean()), "");


			jetQ_unfold[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV" /*	R = %0.1f*/, jetPtLo[j], jetPtHi[j] /*, R */) );
			jetQ_unfold[j]->SetXTitle( Form( "#LTQ_{jet}#GT  (#kappa = %1.1f)" , k) );

			jetQ_err[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV" /*	R = %0.1f*/, jetPtLo[j], jetPtHi[j] /*, R */) );
			jetQ_err[j]->SetXTitle( Form( "#LTQ_{jet}#GT  (#kappa = %1.1f)" , k) );


			jetQ_err[j]->Draw("E2");
			jetQ_unfold[j]->Draw("p same");

			jetQ_p6_p[j]->Draw("p same");

			leg[j]->Draw();

			drawText(Form("R = %0.1f", R), .65, .7, 30);

		}
		else{
			data_2Ds[j]->SetStats(0);
			data_2Ds[j]->GetYaxis()->SetRangeUser(0.0, 0.5);
			if(kappa != "00"){
				data_2Ds[j]->GetYaxis()->SetRangeUser(0.0, 0.3);
			}
			
			data_2Ds[j]->SetMarkerColor(kViolet-3);
			jetQ_p6_d[j]->SetMarkerColor(kGreen+2);

			data_2Ds[j]->SetMarkerStyle(kOpenSquare);
			jetQ_p6_d[j]->SetMarkerStyle(kOpenSquare);


			leg[j]->AddEntry(data_2Ds[j], "Data", "p");
			leg[j]->AddEntry((TObject*)0, Form("mean: %1.3f", data_2Ds[j]->GetMean()), "");

			leg[j]->AddEntry(jetQ_p6_d[j], "Pythia 6 detector level", "p");
			leg[j]->AddEntry((TObject*)0, Form("mean: %1.3f", jetQ_p6_d[j]->GetMean()), "");


			data_2Ds[j]->SetTitle( Form( "%1.0f GeV < p_{T}^{jet} < %1.0f GeV" /*	R = %0.1f*/, jetPtLo[j], jetPtHi[j] /*, R */) );
			data_2Ds[j]->SetXTitle( Form( "#LTQ_{jet}#GT  (#kappa = %1.1f)" , k) );


			cout << "absolute error on Q=0 bin = " << data_2Ds[j]->GetBinError( data_2Ds[j]->FindBin(0.0) ) / data_2Ds[j]->GetBinContent( data_2Ds[j]->FindBin(0.0) ) << "\n";


			data_2Ds[j]->Draw("p same");
			jetQ_p6_d[j]->Draw("p same");


			leg[j]->Draw();

//			drawText("p+p #sqrt{s_{NN}} = 200 GeV}", .25, .85, 15);
			drawText(Form("R = %0.1f", R), .65, .7, 30);

		}
		

		// Save plot here, after points have been plotted
		c[j]->SaveAs(Form("./plots/unfoldedJetQ_jetpt_%1.0f_%1.0f_k" + kappa + "_R" + rad + "_" + lev[level]+ ".png", jetPtLo[j], jetPtHi[j]), "RECREATE");
	}
}
