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
		cerr << "DEBUG: integral of projection histogram = " << proj1Ds[i]->Integral() << "\n";
	}
	return proj1Ds;
}



int main(int argc, const char** argv){
	if(argc != 3){
		cout << "Should receive jet radius, jet charge kappa. Received " << argc-1 << "parameters. Exiting.\n";
		exit(1);
	}

	const int nBins = 4;
	double jetEdges[nBins+1] = {10.0, 15.0, 20.0, 30.0, 40.0};

	double jetPtLo[nBins] = {10.0, 15.0, 20.0, 30.0};
//	double jetPtHi[nBins] = {15.0, 20.0, 30.0, 40.0};

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TH3::SetDefaultSumw2();

	string radius = (string) argv[1];
	string kappa = (string) argv[2];


	// will have to think about how this will be called in macro_submit.csh or how else to call this if this does not work
	// will hadd output files from response.cxx so that response matrix will be across pthat bins
	// 	and apply the "full" response matrix to the full data 2D histogram (Q_jet vs jet pT)
	TFile *file1 = new TFile(("~/ppJCandBF/out/sim/response/p6_R" + radius + "_k" + kappa + "_syst.root").c_str(), "READ"); // file with response in the form of RooUnfoldResponse object
	TFile *file2 = new TFile(("~/ppJCandBF/out/data/hists/ppDataJP2_R" + radius + "_k" + kappa + ".root").c_str(), "READ"); // file to get dat_spectrum
	TFile *file3 = new TFile(("~/ppJCandBF/out/sim/response/stat_err_scaling_R" + radius + "_k" + kappa + ".root").c_str(), "READ"); // scale factors file--- output from stat_err_scaling.cxx
	// ppJCandBF/out/sim/response/stat_err_scaling_R06_k00.root


	//(a) get response object and data spectrum
	RooUnfoldResponse *res = (RooUnfoldResponse*) file1->Get("q_pt_response"); // either named "QvPtResponse" or "q_pt_response"
	TH2D* dat_spectrum = (TH2D*) file2->Get("QvPt_d"); // check what the jet charge histogram name is
	// 2D histogram in data is named QvPt_d
	// 2D histograms in pythia 6 are QvPt_d and QvPt_p for detector/particle level, respectively

	RooUnfoldResponse *res_1D_1015 = (RooUnfoldResponse*) file1->Get("q_response_1015"); // either named "QvPtResponse" or "q_pt_response"
	RooUnfoldResponse *res_1D_1520 = (RooUnfoldResponse*) file1->Get("q_response_1520"); // either named "QvPtResponse" or "q_pt_response"
	RooUnfoldResponse *res_1D_2030 = (RooUnfoldResponse*) file1->Get("q_response_2030"); // either named "QvPtResponse" or "q_pt_response"
	RooUnfoldResponse *res_1D_3040 = (RooUnfoldResponse*) file1->Get("q_response_3040"); // either named "QvPtResponse" or "q_pt_response"

	TH2D* dat_1D_1015 = (TH2D*) file2->Get(("jetQ_k" + kappa + "_jetpt10_15").c_str()); // check what the jet charge histogram name is
	TH2D* dat_1D_1520 = (TH2D*) file2->Get(("jetQ_k" + kappa + "_jetpt15_20").c_str()); // check what the jet charge histogram name is
	TH2D* dat_1D_2030 = (TH2D*) file2->Get(("jetQ_k" + kappa + "_jetpt20_30").c_str()); // check what the jet charge histogram name is
	TH2D* dat_1D_3040 = (TH2D*) file2->Get(("jetQ_k" + kappa + "_jetpt30_40").c_str()); // check what the jet charge histogram name is
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

/*
	RooUnfoldResponse *rMS_2025_nom = (RooUnfoldResponse*) fres->Get((hname+"_res2025_nom").c_str());
	RooUnfoldResponse *rMS_2025_h7smear = (RooUnfoldResponse*) fres->Get((hname+"_res2025_h7smear").c_str());
	RooUnfoldResponse *rMS_2025_p8smear = (RooUnfoldResponse*) fres->Get((hname+"_res2025_p8smear").c_str());
	RooUnfoldResponse *rMS_2530_nom = (RooUnfoldResponse*) fres->Get((hname+"_res2530_nom").c_str());
	RooUnfoldResponse *rMS_2530_h7smear = (RooUnfoldResponse*) fres->Get((hname+"_res2530_h7smear").c_str());
	RooUnfoldResponse *rMS_2530_p8smear = (RooUnfoldResponse*) fres->Get((hname+"_res2530_p8smear").c_str());
	RooUnfoldResponse *rMS_3040_nom = (RooUnfoldResponse*) fres->Get((hname+"_res3040_nom").c_str());
	RooUnfoldResponse *rMS_3040_h7smear = (RooUnfoldResponse*) fres->Get((hname+"_res3040_h7smear").c_str());
	RooUnfoldResponse *rMS_3040_p8smear = (RooUnfoldResponse*) fres->Get((hname+"_res3040_p8smear").c_str());
*/

	//(b) create an object which dictates how the unfolding will happen (e.g. what data spectrum will be unfolded, how many iterations of unfolding, etc.)
	RooUnfoldBayes *unfold_nom = new RooUnfoldBayes(rnom, dat_spectrum, 4, false, "unfold_nom", "");
	RooUnfoldBayes *unfold_IP2 = new RooUnfoldBayes(rnom, dat_spectrum, 2, false, "unfold_IP2", "");
	RooUnfoldBayes *unfold_IP6 = new RooUnfoldBayes(rnom, dat_spectrum, 6, false, "unfold_IP6", "");
	RooUnfoldBayes *unfold_TS = new RooUnfoldBayes(rTS, dat_spectrum, 4, false, "unfold_TS", "");
	RooUnfoldBayes *unfold_TU = new RooUnfoldBayes(rTU, dat_spectrum, 4, false, "unfold_TU", "");
	RooUnfoldBayes *unfold_HC50 = new RooUnfoldBayes(rHC50, dat_spectrum, 4, false, "unfold_HC50", "");
	RooUnfoldBayes *unfold_DS = new RooUnfoldBayes(rDS, dat_spectrum, 4, false, "unfold_DS", "");
	RooUnfoldBayes *unfold_GS = new RooUnfoldBayes(rGS, dat_spectrum, 4, false, "unfold_GS", "");
	
	// 4 is nominal number of iterations

	// 1D jet Q unfolding
	RooUnfoldBayes *unfolded_1D_1015 = new RooUnfoldBayes(res_1D_1015, dat_1D_1015, 4, false, "unfolded_1D_1015", "");
	RooUnfoldBayes *unfolded_1D_1520 = new RooUnfoldBayes(res_1D_1520, dat_1D_1520, 4, false, "unfolded_1D_1520", "");
	RooUnfoldBayes *unfolded_1D_2030 = new RooUnfoldBayes(res_1D_2030, dat_1D_2030, 4, false, "unfolded_1D_2030", "");
	RooUnfoldBayes *unfolded_1D_3040 = new RooUnfoldBayes(res_1D_3040, dat_1D_3040, 4, false, "unfolded_1D_3040", "");



	//(c) [optional] adjustment of unfolding procedure from nominal (e.g. "IncludeSystematics" which makes the unfolding consider the statistical uncertainties of the response)
//	jetQ_unfold->IncludeSystematics(1); // Isaac didn't do it for pp but important when statistics are low
//	jetQ_unfold->SetNToys(1000); // slows it down a lot, Isaac also did not do for pp
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


	TH2D* reco_1D_1015 = (TH2D*) unfolded_1D_1015->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_1D_1520 = (TH2D*) unfolded_1D_1520->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_1D_2030 = (TH2D*) unfolded_1D_2030->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_1D_3040 = (TH2D*) unfolded_1D_3040->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options


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

//	vector<TH1D*> reco_slices_1D = Projection2D(reco_1D, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms


	// stop here to see if things are working..

/*
	//(f) since RooUnfold doesn't handle the statistical uncertainties correctly due to the misses, we scale them up by-hand afterwards
	TH1D* scalefactors = (TH1D*) file3->Get("hratio"); // each bin of this histogram will be a scale factor for the errors of each pT slice of your unfolded 2D hist
	// histograms in file3: "hratio" or "efficiency"----Isaac's code suggests the proper histogram for scalefactors is "hratio"

	// it will scale all bins of jet charge equally, so e.g. if my jet charge distribution for pt = 20-30 GeV has GetBinError() of 0.3, 1.0, 0.09, 0.07, 0.07, 0.1 for Q = -3:-2, -2:-1, -1:0, 0:1, 1:2, 2:3, and the scaling (bin content of the scaling histogram) is 1.5 for 20-30 then you get 0.45, 0.15, 0.135, 0.105, 0.105, 0.15
	for(int i = 0; i < nBins; i++){
//		double scaling = scalefactors->GetBinContent( scalefactors->GetXaxis()->FindBin( jetPtLo[i] ) );
//		if(i == 0){scaling = scalefactors->GetBinContent(scalefactors->GetXaxis()->FindBin(12));} //if 15 is the lowest bin
		
		for(int j = 1; j <= reco_slices[i]->GetNbinsX(); j++){
			reco_noms[i]->SetBinError(j, (double) reco_noms[i]->GetBinError(j)*scaling);
		}
	}
	// reco_slices[i] now have the fully corrected data for each pT range, and we're done!
*/

	for(int i = 0; i < nBins; i++){
		reco_noms[i]->Scale( 1.0/ ( reco_noms[i]->Integral(0, reco_noms[i]->GetNbinsX() + 1) ) );
		reco_IP2s[i]->Scale( 1.0/ ( reco_IP2s[i]->Integral(0, reco_IP2s[i]->GetNbinsX() + 1) ) );
		reco_IP6s[i]->Scale( 1.0/ ( reco_IP6s[i]->Integral(0, reco_IP6s[i]->GetNbinsX() + 1) ) );
		reco_TSs[i]->Scale( 1.0/ ( reco_TSs[i]->Integral(0, reco_TSs[i]->GetNbinsX() + 1) ) );
		reco_TUs[i]->Scale( 1.0/ ( reco_TUs[i]->Integral(0, reco_TUs[i]->GetNbinsX() + 1) ) );
		reco_HC50s[i]->Scale( 1.0/ ( reco_HC50s[i]->Integral(0, reco_HC50s[i]->GetNbinsX() + 1) ) );
		reco_DSs[i]->Scale( 1.0/ ( reco_DSs[i]->Integral(0, reco_DSs[i]->GetNbinsX() + 1) ) );
		reco_GSs[i]->Scale( 1.0/ ( reco_GSs[i]->Integral(0, reco_GSs[i]->GetNbinsX() + 1) ) );

		cerr << "nominal integral(0, nBins + 1) = " << reco_noms[i]->Integral(0, reco_noms[i]->GetNbinsX() + 1) << "\n";
		cerr << "nominal integral(0, nBins) + 1 = " << reco_noms[i]->Integral(0, reco_noms[i]->GetNbinsX()) + 1 << "\n";

		cerr << "integral of reco_noms before scaling = " << reco_noms[i]->Integral() << "\n";


//		reco_noms_1D[i]->Scale( 1.0/ ( reco_slices_1D[i]->Integral(0, reco_slices_1D[i]->GetNbinsX()) + 1) );

		// if normalizing to exclude under/overflow: integral [1, nbins], aka just call ->Integral()
//		reco_slices[i]->Scale( 1.0/ ( reco_slices[i]->Integral() ) );


		// divide by nominal to see effect as percent of nominal
		reco_IP2s[i]->Divide(reco_noms[i]);
		reco_IP6s[i]->Divide(reco_noms[i]);
		reco_TSs[i]->Divide(reco_noms[i]);
		reco_TUs[i]->Divide(reco_noms[i]);
		reco_HC50s[i]->Divide(reco_noms[i]);
		reco_DSs[i]->Divide(reco_noms[i]);
		reco_GSs[i]->Divide(reco_noms[i]);
	}


	
	TFile *fout = new TFile( ("~/ppJCandBF/out/unfold/unfolded_R"+ radius + "_k" + kappa + ".root").c_str() , "RECREATE");
	fout->cd();
	
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

/*
	for(int i = 0; i < nBins; i++){
		reco_slices_1D[i]->Write();
	}
*/	
	res->Write();
	dat_spectrum->Write();
	reco_nom->Write();
	

	res_1D_1015->Write();
	res_1D_1520->Write();
	res_1D_2030->Write();
	res_1D_3040->Write();

	dat_1D_1015->Write();
	dat_1D_1520->Write();
	dat_1D_2030->Write();
	dat_1D_3040->Write();

	reco_1D_1015->Write();
	reco_1D_1520->Write();
	reco_1D_2030->Write();
	reco_1D_3040->Write();


	fout->Close();
	
	return 0;
}

