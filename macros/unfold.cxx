#include <ctime>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <algorithm>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TRandom.h>

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "Headers/funcs.hh"


using namespace std;


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
		cerr << "low edge: " << low_rough << "\n";
		cerr << "high edge: " << high_rough << "\n";
		if (axis == "x" || axis == "X" || axis == "1") {
			cout << "now including e option" << endl;
			proj1Ds.push_back(hist2D->ProjectionX((hist2D->GetName() + axis + low_rough + high_rough).c_str(), hist2D->GetYaxis()->FindBin(ranges[i]), hist2D->GetYaxis()->FindBin(ranges[i+1]) - 1, "e"));
			// testing above: 
//			proj1Ds.push_back(hist2D->ProjectionX((hist2D->GetName() + axis + low_rough + high_rough).c_str(), ranges[i], ranges[i+1] - 1, "e"));
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
        
        cout << "DEBUG projection? err/value of Q=0 bin for jet pt bin above " << ranges[i] << " = " << proj1Ds[i]->GetBinError( proj1Ds[i]->FindBin(0.0) ) / proj1Ds[i]->GetBinContent( proj1Ds[i]->FindBin(0.0) ) << "\n";
	}
	return proj1Ds;
}



int main(int argc, const char** argv){
	if(argc != 3){
		cout << "Should receive jet radius, jet charge kappa. Received " << argc-1 << "parameters. Exiting.\n";
		exit(1);
	}

	const int nBins = 4;
	double jetEdges[nBins+1] = {15.0, 20.0, 25.0, 30.0, 40.0};

	double jetPtLo[nBins] = {15.0, 20.0, 25.0, 30.0};
//	double jetPtHi[nBins] = {20.0, 25.0, 30.0, 40.0};

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TH3::SetDefaultSumw2();

	string radius = (string) argv[1];
	string kappa = (string) argv[2];


	// will have to think about how this will be called in macro_submit.csh or how else to call this if this does not work
	// will hadd output files from response.cxx so that response matrix will be across pthat bins
	// 	and apply the "full" response matrix to the full data 2D histogram (Q_jet vs jet pT)
	// May 26 is latest version, with updated jet pt det-level cut at 15 GeV, no mass assignment, no mass cut
	//	and after achieving consistency with Isaac's jet mass analysis, so it is from a "known good" state
	//out/sim/hists/p6_R04_k05_June7.root
	TFile *file1 = new TFile(("~/ppJCandBF/out/sim/response/p6_R" + radius + "_k" + kappa + "_syst_July20.root").c_str(), "READ"); // file with response in the form of RooUnfoldResponse object
	TFile *file2 = new TFile(("~/ppJCandBF/out/data/hists/ppDataJP2_R" + radius + "_k" + kappa + "_July20.root").c_str(), "READ"); // file to get dat_spectrum
	TFile *file3 = new TFile(("~/ppJCandBF/out/sim/response/stat_err_scaling_R" + radius + "_k" + kappa + ".root").c_str(), "READ"); // scale factors file--- output from stat_err_scaling.cxx
	// ppJCandBF/out/sim/response/stat_err_scaling_R06_k00.root


	//(a) get response object and data spectrum
	RooUnfoldResponse *res = (RooUnfoldResponse*) file1->Get("q_pt_response"); // either named "QvPtResponse" or "q_pt_response"
	TH2D* dat_spectrum = (TH2D*) file2->Get("QvPt_d"); // check what the jet charge histogram name is
	// 2D histogram in data is named QvPt_d
	// 2D histograms in pythia 6 are QvPt_d and QvPt_p for detector/particle level, respectively

/*
	RooUnfoldResponse *res_1D_1520 = (RooUnfoldResponse*) file1->Get("q_res1520_nom"); // either named "QvPtResponse" or "q_pt_response"
	RooUnfoldResponse *res_1D_2025 = (RooUnfoldResponse*) file1->Get("q_res2025_nom"); // either named "QvPtResponse" or "q_pt_response"
	RooUnfoldResponse *res_1D_2530 = (RooUnfoldResponse*) file1->Get("q_res2530_nom"); // either named "QvPtResponse" or "q_pt_response"
	RooUnfoldResponse *res_1D_3040 = (RooUnfoldResponse*) file1->Get("q_res3040_nom"); // either named "QvPtResponse" or "q_pt_response"
*/


	TH1D* dat_1D_1520 = (TH1D*) file2->Get(("jetQ_k" + kappa + "_jetpt15_20").c_str()); // check what the jet charge histogram name is
	TH1D* dat_1D_2025 = (TH1D*) file2->Get(("jetQ_k" + kappa + "_jetpt20_25").c_str()); // check what the jet charge histogram name is
	TH1D* dat_1D_2530 = (TH1D*) file2->Get(("jetQ_k" + kappa + "_jetpt25_30").c_str()); // check what the jet charge histogram name is
	TH1D* dat_1D_3040 = (TH1D*) file2->Get(("jetQ_k" + kappa + "_jetpt30_40").c_str()); // check what the jet charge histogram name is
	// Form( "jetQ_k" + k_str + "_jetpt%1.0f_%1.0f", jetPtLo[ij], jetPtHi[ij] )


	// started adding this 4/8/22
	// systematics responses:
	// TS: tower scaling; TU: tracking uncertainty; HC50: hadronic correction = .50; DS,GS: detector, generator spectra smearing;
	// e.g. h7smear: smearing prior with e.g. herwig 7;
	// nom: nominal (to which to compare)
	RooUnfoldResponse *rnom = (RooUnfoldResponse*) file1->Get("q_pt_res_nom");
	RooUnfoldResponse *rTS = (RooUnfoldResponse*) file1->Get("q_pt_res_TS"); 
	RooUnfoldResponse *rTU = (RooUnfoldResponse*) file1->Get("q_pt_res_TU");
	RooUnfoldResponse *rHC50 = (RooUnfoldResponse*) file1->Get("q_pt_res_HC50");
	RooUnfoldResponse *rDS = (RooUnfoldResponse*) file1->Get("q_pt_res_DS");
	RooUnfoldResponse *rGS = (RooUnfoldResponse*) file1->Get("q_pt_res_GS");
	//  RooUnfoldResponse *rMS = (RooUnfoldResponse*) file1->Get("q_pt_res_QS");

//cout << "1D Q responses gotten from response file\n";


	RooUnfoldResponse *rQS_2025_nom = (RooUnfoldResponse*) file1->Get("q_res2025_nom");
	RooUnfoldResponse *rQS_2025_h7smear = (RooUnfoldResponse*) file1->Get("q_res2025_h7smear");
	RooUnfoldResponse *rQS_2025_p8smear = (RooUnfoldResponse*) file1->Get("q_res2025_p8smear");
	RooUnfoldResponse *rQS_2530_nom = (RooUnfoldResponse*) file1->Get("q_res2530_nom");
	RooUnfoldResponse *rQS_2530_h7smear = (RooUnfoldResponse*) file1->Get("q_res2530_h7smear");
	RooUnfoldResponse *rQS_2530_p8smear = (RooUnfoldResponse*) file1->Get("q_res2530_p8smear");
	RooUnfoldResponse *rQS_3040_nom = (RooUnfoldResponse*) file1->Get("q_res3040_nom");
	RooUnfoldResponse *rQS_3040_h7smear = (RooUnfoldResponse*) file1->Get("q_res3040_h7smear");
	RooUnfoldResponse *rQS_3040_p8smear = (RooUnfoldResponse*) file1->Get("q_res3040_p8smear");

//cout << "systematic q smearing responses gotten from response file\n";


	//(b) create an object which dictates how the unfolding will happen (e.g. what data spectrum will be unfolded, how many iterations of unfolding, etc.)
	RooUnfoldBayes *unfold_nom = new RooUnfoldBayes(rnom, dat_spectrum, 4, false, "unfold_nom", "");
	RooUnfoldBayes *unfold_IP2 = new RooUnfoldBayes(rnom, dat_spectrum, 2, false, "unfold_IP2", "");
	RooUnfoldBayes *unfold_IP6 = new RooUnfoldBayes(rnom, dat_spectrum, 6, false, "unfold_IP6", "");
	RooUnfoldBayes *unfold_TS = new RooUnfoldBayes(rTS, dat_spectrum, 4, false, "unfold_TS", "");
	RooUnfoldBayes *unfold_TU = new RooUnfoldBayes(rTU, dat_spectrum, 4, false, "unfold_TU", "");
	RooUnfoldBayes *unfold_HC50 = new RooUnfoldBayes(rHC50, dat_spectrum, 4, false, "unfold_HC50", "");
	RooUnfoldBayes *unfold_DS = new RooUnfoldBayes(rDS, dat_spectrum, 4, false, "unfold_DS", "");
	RooUnfoldBayes *unfold_GS = new RooUnfoldBayes(rGS, dat_spectrum, 4, false, "unfold_GS", "");
	

	RooUnfoldBayes *unfold_nom1D_2025 = new RooUnfoldBayes(rQS_2025_nom, dat_1D_2025, 4, false, "unfold_nom1D_2025", "");
	RooUnfoldBayes *unfold_h7smear1D_2025 = new RooUnfoldBayes(rQS_2025_h7smear, dat_1D_2025, 4, false, "unfold_h7smear1D_2025", "");
	RooUnfoldBayes *unfold_p8smear1D_2025 = new RooUnfoldBayes(rQS_2025_p8smear, dat_1D_2025, 4, false, "unfold_p8smear1D_2025", "");

	RooUnfoldBayes *unfold_nom1D_2530 = new RooUnfoldBayes(rQS_2530_nom, dat_1D_2530, 4, false, "unfold_nom1D_2530", "");
	RooUnfoldBayes *unfold_h7smear1D_2530 = new RooUnfoldBayes(rQS_2530_h7smear, dat_1D_2530, 4, false, "unfold_h7smear1D_2530", "");
	RooUnfoldBayes *unfold_p8smear1D_2530 = new RooUnfoldBayes(rQS_2530_p8smear, dat_1D_2530, 4, false, "unfold_p8smear1D_2530", "");

	RooUnfoldBayes *unfold_nom1D_3040 = new RooUnfoldBayes(rQS_3040_nom, dat_1D_3040, 4, false, "unfold_nom1D_3040", "");
	RooUnfoldBayes *unfold_h7smear1D_3040 = new RooUnfoldBayes(rQS_3040_h7smear, dat_1D_3040, 4, false, "unfold_h7smear1D_3040", "");
	RooUnfoldBayes *unfold_p8smear1D_3040 = new RooUnfoldBayes(rQS_3040_p8smear, dat_1D_3040, 4, false, "unfold_p8smear1D_3040", "");

//cout << "Q smearing systematic unfolding performed\n";


	// 4 is nominal number of iterations

/*
	// 1D jet Q unfolding
	RooUnfoldBayes *unfolded_1D_1520 = new RooUnfoldBayes(res_1D_1520, dat_1D_1520, 4, false, "unfolded_1D_1520", "");
	RooUnfoldBayes *unfolded_1D_2025 = new RooUnfoldBayes(res_1D_2025, dat_1D_2025, 4, false, "unfolded_1D_2025", "");
	RooUnfoldBayes *unfolded_1D_2530 = new RooUnfoldBayes(res_1D_2530, dat_1D_2530, 4, false, "unfolded_1D_2530", "");
	RooUnfoldBayes *unfolded_1D_3040 = new RooUnfoldBayes(res_1D_3040, dat_1D_3040, 4, false, "unfolded_1D_3040", "");
*/
//cout << "1D Q unfolding perfromed\n";


	//(c) [optional] adjustment of unfolding procedure from nominal (e.g. "IncludeSystematics" which makes the unfolding consider the statistical uncertainties of the response)
//	unfold_nom->IncludeSystematics(1); // Isaac didn't do it for pp but important when statistics are low
//	unfold_nom->SetNToys(1000); // slows it down a lot, Isaac also did not do for pp
	// come back to this...

	//(d) do the unfolding
	TH2D* reco_nom = (TH2D*) unfold_nom->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_IP2 = (TH2D*) unfold_IP2->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_IP6 = (TH2D*) unfold_IP6->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_TS = (TH2D*) unfold_TS->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_TU = (TH2D*) unfold_TU->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_HC50 = (TH2D*) unfold_HC50->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_DS = (TH2D*) unfold_DS->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_GS = (TH2D*) unfold_GS->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options


	TH1D* reco_nom1D_2025 = (TH1D*) unfold_nom1D_2025->Hreco((RooUnfold::ErrorTreatment) 3);
	TH1D* reco_h7smear1D_2025 = (TH1D*) unfold_h7smear1D_2025->Hreco((RooUnfold::ErrorTreatment) 3);
	TH1D* reco_p8smear1D_2025 = (TH1D*) unfold_p8smear1D_2025->Hreco((RooUnfold::ErrorTreatment) 3);
	TH1D* reco_nom1D_2530 = (TH1D*) unfold_nom1D_2530->Hreco((RooUnfold::ErrorTreatment) 3);
	TH1D* reco_h7smear1D_2530 = (TH1D*) unfold_h7smear1D_2530->Hreco((RooUnfold::ErrorTreatment) 3);
	TH1D* reco_p8smear1D_2530 = (TH1D*) unfold_p8smear1D_2530->Hreco((RooUnfold::ErrorTreatment) 3);
	TH1D* reco_nom1D_3040 = (TH1D*) unfold_nom1D_3040->Hreco((RooUnfold::ErrorTreatment) 3);
	TH1D* reco_h7smear1D_3040 = (TH1D*) unfold_h7smear1D_3040->Hreco((RooUnfold::ErrorTreatment) 3);
	TH1D* reco_p8smear1D_3040 = (TH1D*) unfold_p8smear1D_3040->Hreco((RooUnfold::ErrorTreatment) 3);



/*
	TH2D* reco_1D_1520 = (TH2D*) unfolded_1D_1520->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_1D_2025 = (TH2D*) unfolded_1D_2025->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_1D_2530 = (TH2D*) unfolded_1D_2530->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_1D_3040 = (TH2D*) unfolded_1D_3040->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
*/




	//(e) project onto 1D histograms
	vector<TH1D*> reco_noms = Projection2D(reco_nom, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_IP2s = Projection2D(reco_IP2, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_IP6s = Projection2D(reco_IP6, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_TSs = Projection2D(reco_TS, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_TUs = Projection2D(reco_TU, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_HC50s = Projection2D(reco_HC50, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_DSs = Projection2D(reco_DS, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_GSs = Projection2D(reco_GS, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	// project the 2D onto 1Ds for the pT ranges I select. see functions "Projection2D" from unfold.cxx on Isaac's github

	//vector<TH1D*> reco_slices_1D = Projection2D(reco_1D, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms


	// stop here to see if things are working..


	//(f) since RooUnfold doesn't handle the statistical uncertainties correctly due to the misses, we scale them up by-hand afterwards
	TH1D* scalefactors = (TH1D*) file3->Get("hratio"); // each bin of this histogram will be a scale factor for the errors of each pT slice of your unfolded 2D hist
	// histograms in file3: "hratio" or "efficiency"----Isaac's code suggests the proper histogram for scalefactors is "hratio"

	// it will scale all bins of jet charge equally, so e.g. if my jet charge distribution for pt = 20-30 GeV has GetBinError() of 0.3, 1.0, 0.09, 0.07, 0.07, 0.1 for Q = -3:-2, -2:-1, -1:0, 0:1, 1:2, 2:3, and the scaling (bin content of the scaling histogram) is 1.5 for 20-30 then you get 0.45, 0.15, 0.135, 0.105, 0.105, 0.15
	for(int i = 0; i < nBins; i++){
        double scaling = -1;
        if(i == 1){scaling = scalefactors->GetBinContent(scalefactors->GetXaxis()->FindBin(20));}
        if(i == 2){scaling = scalefactors->GetBinContent(scalefactors->GetXaxis()->FindBin(25));} //if 15 is the lowest bin
        if(i == 3){scaling = max( scalefactors->GetBinContent(scalefactors->GetXaxis()->FindBin(30)), scalefactors->GetBinContent(scalefactors->GetXaxis()->FindBin(35)) );}

        cout << "\njet pt bin low end: " << jetEdges[i] << "\n";

		for(int j = 1; j <= reco_noms[i]->GetNbinsX(); j++){
            double binerr = reco_noms[i]->GetBinError(j);
			reco_noms[i]->SetBinError(j, (double) binerr*scaling);
            
            cout << "bin center: " << reco_noms[i]->GetBinCenter(j) << "\n";
            cout << "bin content: " << reco_noms[i]->GetBinContent(j) << "\n";
            cout << "bin error before scaling = " << binerr << "\n";
            cout << "bin content/bin error before scaling: " << binerr/reco_noms[i]->GetBinContent(j) << "\n";
            cout << "scaling = " << scaling << "\n";
            cout << "bin error after scaling = " << binerr*scaling << "\n";
            cout << "bin content/bin error after scaling: " << binerr*scaling/reco_noms[i]->GetBinContent(j) << "\n";
        }
	}
	// reco_slices[i] now have the fully corrected data for each pT range, and we're done!


	for(int i = 0; i < nBins; i++){
		reco_noms[i]->Scale( 1.0 / ( reco_noms[i]->Integral(0, reco_noms[i]->GetNbinsX() + 1) ) );
		reco_IP2s[i]->Scale( 1.0 / ( reco_IP2s[i]->Integral(0, reco_IP2s[i]->GetNbinsX() + 1) ) );
		reco_IP6s[i]->Scale( 1.0 / ( reco_IP6s[i]->Integral(0, reco_IP6s[i]->GetNbinsX() + 1) ) );
		reco_TSs[i]->Scale( 1.0 / ( reco_TSs[i]->Integral(0, reco_TSs[i]->GetNbinsX() + 1) ) );
		reco_TUs[i]->Scale( 1.0 / ( reco_TUs[i]->Integral(0, reco_TUs[i]->GetNbinsX() + 1) ) );
		reco_HC50s[i]->Scale( 1.0 / ( reco_HC50s[i]->Integral(0, reco_HC50s[i]->GetNbinsX() + 1) ) );
		reco_DSs[i]->Scale( 1.0 / ( reco_DSs[i]->Integral(0, reco_DSs[i]->GetNbinsX() + 1) ) );
		reco_GSs[i]->Scale( 1.0 / ( reco_GSs[i]->Integral(0, reco_GSs[i]->GetNbinsX() + 1) ) );


		// divide by nominal to see effect as percent of nominal
		reco_IP2s[i]->Divide(reco_noms[i]);
		reco_IP6s[i]->Divide(reco_noms[i]);
		reco_TSs[i]->Divide(reco_noms[i]);
		reco_TUs[i]->Divide(reco_noms[i]);
		reco_HC50s[i]->Divide(reco_noms[i]);
		reco_DSs[i]->Divide(reco_noms[i]);
		reco_GSs[i]->Divide(reco_noms[i]);
	}

//cout << "scaling 1D q smearing unfolded distributions\n";

	reco_nom1D_2025->Scale( 1.0 / reco_nom1D_2025->Integral() );
	reco_nom1D_2530->Scale( 1.0 / reco_nom1D_2530->Integral() );
	reco_nom1D_3040->Scale( 1.0 / reco_nom1D_3040->Integral() );

	reco_h7smear1D_2025->Scale( 1.0 / reco_h7smear1D_2025->Integral() );
	reco_h7smear1D_2530->Scale( 1.0 / reco_h7smear1D_2530->Integral() );
	reco_h7smear1D_3040->Scale( 1.0 / reco_h7smear1D_3040->Integral() );

	reco_p8smear1D_2025->Scale( 1.0 / reco_p8smear1D_2025->Integral() );
	reco_p8smear1D_2530->Scale( 1.0 / reco_p8smear1D_2530->Integral() );
	reco_p8smear1D_3040->Scale( 1.0 / reco_p8smear1D_3040->Integral() );


//cout << "taking 1D q smearing unfolded distribution ratio to nominal\n";

	reco_h7smear1D_2025->Divide( reco_nom1D_2025 );
//cout << "point a\n";
	reco_h7smear1D_2530->Divide( reco_nom1D_2530 );
//cout << "point b\n";
	reco_h7smear1D_3040->Divide( reco_nom1D_3040 );
//cout << "point c\n";

	reco_p8smear1D_2025->Divide( reco_nom1D_2025 );
//cout << "point d\n";
	reco_p8smear1D_2530->Divide( reco_nom1D_2530 );
//cout << "point e\n";
	reco_p8smear1D_3040->Divide( reco_nom1D_3040 );
//cout << "point f\n";



	double stats[5] = {0, 0, 0, 0, 0};

	reco_h7smear1D_2025->PutStats(stats);
	reco_h7smear1D_2530->PutStats(stats);
	reco_h7smear1D_3040->PutStats(stats);
	reco_p8smear1D_2025->PutStats(stats);
	reco_p8smear1D_2530->PutStats(stats);
	reco_p8smear1D_3040->PutStats(stats);

	reco_h7smear1D_2025->Sumw2(0);
	reco_h7smear1D_2530->Sumw2(0);
	reco_h7smear1D_3040->Sumw2(0);
	reco_p8smear1D_2025->Sumw2(0);
	reco_p8smear1D_2530->Sumw2(0);
	reco_p8smear1D_3040->Sumw2(0);



cout << "setting histogram content to % difference compared to nominal\n";

	for(int j = 1; j <= (int) reco_h7smear1D_2025->GetNbinsX(); j++){
		reco_h7smear1D_2025->SetBinContent( j, fabs( reco_h7smear1D_2025->GetBinContent(j) - 1 ) );
		reco_p8smear1D_2025->SetBinContent( j, fabs( reco_p8smear1D_2025->GetBinContent(j) - 1 ) );

		reco_h7smear1D_2530->SetBinContent( j, fabs( reco_h7smear1D_2530->GetBinContent(j) - 1 ) );
		reco_p8smear1D_2530->SetBinContent( j, fabs( reco_p8smear1D_2530->GetBinContent(j) - 1 ) );

		reco_h7smear1D_3040->SetBinContent( j, fabs( reco_h7smear1D_3040->GetBinContent(j) - 1 ) );
		reco_p8smear1D_3040->SetBinContent( j, fabs( reco_p8smear1D_3040->GetBinContent(j) - 1 ) );
		
	}


	for(int i = 0; i < nBins; i++){

		reco_IP2s[i]->PutStats(stats);
		reco_IP6s[i]->PutStats(stats);
		reco_TSs[i]->PutStats(stats);
		reco_TUs[i]->PutStats(stats);
		reco_HC50s[i]->PutStats(stats);
		reco_DSs[i]->PutStats(stats);
		reco_GSs[i]->PutStats(stats);

		reco_IP2s[i]->Sumw2(0);
		reco_IP6s[i]->Sumw2(0);
		reco_TSs[i]->Sumw2(0);
		reco_TUs[i]->Sumw2(0);
		reco_HC50s[i]->Sumw2(0);
		reco_DSs[i]->Sumw2(0);
		reco_GSs[i]->Sumw2(0);


		for(int j = 1; j <= (int) reco_IP2s[i]->GetXaxis()->GetNbins(); j++){
			reco_IP2s[i]->SetBinContent( j, fabs( reco_IP2s[i]->GetBinContent(j) - 1 ) );
			reco_IP6s[i]->SetBinContent( j, fabs( reco_IP6s[i]->GetBinContent(j) - 1 ) );
			reco_TSs[i]->SetBinContent( j, fabs( reco_TSs[i]->GetBinContent(j) - 1 ) );
			reco_TUs[i]->SetBinContent( j, fabs( reco_TUs[i]->GetBinContent(j) - 1 ) );
			reco_HC50s[i]->SetBinContent( j, fabs( reco_HC50s[i]->GetBinContent(j) - 1 ) );
			reco_DSs[i]->SetBinContent( j, fabs( reco_DSs[i]->GetBinContent(j) - 1 ) );
			reco_GSs[i]->SetBinContent( j, fabs( reco_GSs[i]->GetBinContent(j) - 1 ) );
			
		}
	}


//cout << "done setting histogram content to % difference compared to nominal\n";
//cout << "ready to write to output file\n";

	
	TFile *fout = new TFile( ("~/ppJCandBF/out/unfold/unfolded_R"+ radius + "_k" + kappa + ".root").c_str() , "RECREATE");
	fout->cd();

	cout << "output file name: " << ("~/ppJCandBF/out/unfold/unfolded_R"+ radius + "_k" + kappa + ".root").c_str() << "\n";
	
	for(int i = 0; i < nBins; i++){
		reco_noms[i]->Write();
		reco_IP2s[i]->Write();
		reco_IP6s[i]->Write();
		reco_TSs[i]->Write();
		reco_TUs[i]->Write();
		reco_HC50s[i]->Write();
		reco_DSs[i]->Write();
		reco_GSs[i]->Write();
	}

	reco_h7smear1D_2025->Write();
	reco_h7smear1D_2530->Write();
	reco_h7smear1D_3040->Write();

	reco_p8smear1D_2025->Write();
	reco_p8smear1D_2530->Write();
	reco_p8smear1D_3040->Write();


/*
	for(int i = 0; i < nBins; i++){
		reco_slices_1D[i]->Write();
	}
*/	
	res->Write();
	dat_spectrum->Write();
	reco_nom->Write();
	

/*
	res_1D_1520->Write();
	res_1D_2025->Write();
	res_1D_2530->Write();
	res_1D_3040->Write();
*/

	dat_1D_1520->Write();
	dat_1D_2025->Write();
	dat_1D_2530->Write();
	dat_1D_3040->Write();

/*
	reco_1D_1520->Write();
	reco_1D_2025->Write();
	reco_1D_2530->Write();
	reco_1D_3040->Write();
*/

	fout->Close();
	
	return 0;
}

