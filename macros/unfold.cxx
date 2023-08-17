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
			proj1Ds.push_back(hist2D->ProjectionX((hist2D->GetName() + axis + low_rough + high_rough).c_str(), hist2D->GetYaxis()->FindBin( ranges[i] ), hist2D->GetYaxis()->FindBin( ranges[i+1] ) - 1, "e"));
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

        //cout << "DEBUG projection? err/value of Q=0 bin for jet pt bin above " << ranges[i] << " = " << proj1Ds[i]->GetBinError( proj1Ds[i]->FindBin(0.0) ) / proj1Ds[i]->GetBinContent( proj1Ds[i]->FindBin(0.0) ) << "\n";
	}
	return proj1Ds;
}




int main(int argc, const char** argv){
	if(argc != 3){
		cout << "Should receive jet radius, jet charge kappa. Received " << argc-1 << "parameters. Exiting.\n";
		exit(1);
	}


	const int nBins = 4;
    int nBins_Q = 13; // change for different kappas, but for k = 0.0, 13 bins from -6.5 to 6.5

	double jetEdges[nBins+1] = {15.0, 20.0, 25.0, 30.0, 40.0};

	double jetPtLo[nBins] = {15.0, 20.0, 25.0, 30.0};
	double jetPtHi[nBins] = {20.0, 25.0, 30.0, 40.0};

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TH3::SetDefaultSumw2();

	string radius = (string) argv[1];
	string kappa = (string) argv[2];


        // add _bindropped at end of file name after implementing bin_drop.cxx from Isaac to remove bins with < 20 counts from response, histograms

	TFile *file1 = new TFile(("./out/sim/response/p6_R" + radius + "_k" + kappa + /*"_bindropped" + */ ".root").c_str(), "READ"); // file with response in the form of RooUnfoldResponse object

	TFile *file2 = new TFile(("./out/data/hists/ppJP2data_R" + radius + "_k" + kappa + /*"_bindropped" */ ".root").c_str(), "READ"); // file to get dat_spectrum

	TFile *file3 = new TFile(("./out/sim/response/stat_err_scaling_R" + radius + "_k" + kappa + ".root").c_str(), "READ"); // scale factors file--- output from stat_err_scaling.cxx
	// ~/thesis/out/sim/response/stat_err_scaling_R06_k00.root




	//(a) get response object and data spectrum
	RooUnfoldResponse *res = (RooUnfoldResponse*) file1->Get("q_pt_response"); // either named "QvPtResponse" or "q_pt_response"
	TH2D* dat_spectrum = (TH2D*) file2->Get("QvPt_d"); // check what the jet charge histogram name is
	// 2D histogram in data is named QvPt_d
	// 2D histograms in pythia 6 are QvPt_d and QvPt_p for detector/particle level, respectively


    // for opposite-side closure, unfold A with response from B
    TH2D* A_qvpt_d = (TH2D*) file1->Get( "sampleA_q_pt_det" );
    TH2D* A_qvpt_p = (TH2D*) file1->Get( "sampleA_q_pt_gen" );

    vector<TH1D*> A_det = Projection2D(A_qvpt_d, nBins, jetEdges, "x");
    vector<TH1D*> A_gen = Projection2D(A_qvpt_p, nBins, jetEdges, "x");


    RooUnfoldResponse* clos_res = (RooUnfoldResponse*) file1->Get( "sampleB_q_pt_response" );


    // for same-side closure, unfold B with response from B
    TH2D* B_qvpt_d = (TH2D*) file1->Get( "sampleB_q_pt_det" );
    TH2D* B_qvpt_p = (TH2D*) file1->Get( "sampleB_q_pt_gen" );

    vector<TH1D*> B_det = Projection2D(B_qvpt_d, nBins, jetEdges, "x");
    vector<TH1D*> B_gen = Projection2D(B_qvpt_p, nBins, jetEdges, "x");


    nBins_Q = dat_spectrum->GetNbinsX();
    double Q_max = -1.0 * dat_spectrum->GetXaxis()->GetBinLowEdge(1); // low edge of bin 1 is < 0, q range is from this value to its negative (which will be > 0)

/*
	RooUnfoldResponse *res_1D_1520 = (RooUnfoldResponse*) file1->Get("q_res1520_nom"); // either named "QvPtResponse" or "q_pt_response"
	RooUnfoldResponse *res_1D_2025 = (RooUnfoldResponse*) file1->Get("q_res2025_nom"); // either named "QvPtResponse" or "q_pt_response"
	RooUnfoldResponse *res_1D_2530 = (RooUnfoldResponse*) file1->Get("q_res2530_nom"); // either named "QvPtResponse" or "q_pt_response"
	RooUnfoldResponse *res_1D_3040 = (RooUnfoldResponse*) file1->Get("q_res3040_nom"); // either named "QvPtResponse" or "q_pt_response"
*/


	//TH1D* dat_1D_1520 = (TH1D*) file2->Get( "jetQ_jetpt1520" ); // check what the jet charge histogram name is
	TH1D* dat_1D_2025 = (TH1D*) file2->Get( "jetQ_jetpt2025" ); // check what the jet charge histogram name is
	TH1D* dat_1D_2530 = (TH1D*) file2->Get( "jetQ_jetpt2530" ); // check what the jet charge histogram name is
	TH1D* dat_1D_3040 = (TH1D*) file2->Get( "jetQ_jetpt3040" ); // check what the jet charge histogram name is
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
    RooUnfoldBayes *unfold_test = new RooUnfoldBayes(res, dat_spectrum, 4, false, "unfold_test", "");
    RooUnfoldBayes *unfold_nom = new RooUnfoldBayes(rnom, dat_spectrum, 4, false, "unfold_nom", "");
	RooUnfoldBayes *unfold_IP2 = new RooUnfoldBayes(rnom, dat_spectrum, 2, false, "unfold_IP2", "");
	RooUnfoldBayes *unfold_IP6 = new RooUnfoldBayes(rnom, dat_spectrum, 6, false, "unfold_IP6", "");
	RooUnfoldBayes *unfold_TS = new RooUnfoldBayes(rTS, dat_spectrum, 4, false, "unfold_TS", "");
	RooUnfoldBayes *unfold_TU = new RooUnfoldBayes(rTU, dat_spectrum, 4, false, "unfold_TU", "");
	RooUnfoldBayes *unfold_HC50 = new RooUnfoldBayes(rHC50, dat_spectrum, 4, false, "unfold_HC50", "");
	RooUnfoldBayes *unfold_DS = new RooUnfoldBayes(rDS, dat_spectrum, 4, false, "unfold_DS", "");
	RooUnfoldBayes *unfold_GS = new RooUnfoldBayes(rGS, dat_spectrum, 4, false, "unfold_GS", "");


    int N_iters = 10;

    RooUnfoldBayes* same_side[N_iters]; RooUnfoldBayes* opp_side[N_iters];
    for(int ii = 0; ii < N_iters; ii++){
        same_side[ii] = new RooUnfoldBayes(clos_res, B_qvpt_d, ii+1, false, Form( "same_%iiter", ii+1 ), "");
        opp_side[ii] = new RooUnfoldBayes(clos_res, A_qvpt_d, ii+1, false, Form( "opp_%iiter", ii+1 ), "");
    }



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
    
    cout << "\nNOMINAL UNFOLDING: systematics included? " << unfold_nom->SystematicsIncluded() << "\n";
    
    
	//(d) do the unfolding
    TH2D* reco_test = (TH2D*) unfold_test->Hreco((RooUnfold::ErrorTreatment) 3);
    
	TH2D* reco_nom = (TH2D*) unfold_nom->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_IP2 = (TH2D*) unfold_IP2->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_IP6 = (TH2D*) unfold_IP6->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_TS = (TH2D*) unfold_TS->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_TU = (TH2D*) unfold_TU->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_HC50 = (TH2D*) unfold_HC50->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_DS = (TH2D*) unfold_DS->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options
	TH2D* reco_GS = (TH2D*) unfold_GS->Hreco((RooUnfold::ErrorTreatment) 3); // 3 is nominal for our group, see link in email for reference of other options


    TH2D* same_reco[N_iters]; TH2D* opp_reco[N_iters];
    for(int ii = 0; ii < N_iters; ii++){
        same_reco[ii] = (TH2D*) same_side[ii]->Hreco((RooUnfold::ErrorTreatment) 3);
        opp_reco[ii] = (TH2D*) opp_side[ii]->Hreco((RooUnfold::ErrorTreatment) 3);
    }



//    cout << "DEBUG: Where is the break? 1\n";

	TH1D* reco_nom1D_2025 = (TH1D*) unfold_nom1D_2025->Hreco((RooUnfold::ErrorTreatment) 3);
    // this line is the issue???


//    cout << "DEBUG: Where is the break? 2a\n";

    TH1D* reco_h7smear1D_2025 = (TH1D*) unfold_h7smear1D_2025->Hreco((RooUnfold::ErrorTreatment) 3);
//    cout << "DEBUG: Where is the break? 2b\n";
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
    vector<TH1D*> reco_tests = Projection2D(reco_test, nBins, jetEdges, "x");
    
	vector<TH1D*> reco_noms = Projection2D(reco_nom, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_IP2s = Projection2D(reco_IP2, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_IP6s = Projection2D(reco_IP6, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_TSs = Projection2D(reco_TS, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_TUs = Projection2D(reco_TU, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_HC50s = Projection2D(reco_HC50, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_DSs = Projection2D(reco_DS, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	vector<TH1D*> reco_GSs = Projection2D(reco_GS, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms
	// project the 2D onto 1Ds for the pT ranges I select. see functions "Projection2D" from unfold.cxx on Isaac's github

//    cout << "DEBUG: Where is the break? 3\n";


    vector<TH1D*> recos_same[N_iters]; vector<TH1D*> recos_opp[N_iters];
    for(int ii = 0; ii < N_iters; ii++){
        recos_same[ii] = Projection2D(same_reco[ii], nBins, jetEdges, "x");
        recos_opp[ii] = Projection2D(opp_reco[ii], nBins, jetEdges, "x");
    }


//    cout << "DEBUG: Where is the break? 4\n";


	//vector<TH1D*> reco_slices_1D = Projection2D(reco_1D, nBins, jetEdges, "x"); // i think--- Isaac's custom function will return a vector of histograms

    // vector to store non-normalized unfolded distribution since template fitting requires statistics as counts, so needs non-normalized
    vector<TH1D*> reco_noms_scale;

	// jet charge distributions for the systematics, don't turn these into % errors
    vector<TH1D*> reco_tests_copy;
    
	vector<TH1D*> reco_noms_copy;
	vector<TH1D*> reco_IP2s_copy;
	vector<TH1D*> reco_IP6s_copy;
	vector<TH1D*> reco_TSs_copy;
	vector<TH1D*> reco_TUs_copy;
	vector<TH1D*> reco_HC50s_copy;
	vector<TH1D*> reco_DSs_copy;
	vector<TH1D*> reco_GSs_copy;

    vector<TH1D*> reco_P8s_copy;
    vector<TH1D*> reco_H7s_copy;


//    cout << "DEBUG: Where is the break? 5\n";

    for(int i = 0; i < nBins; i++){
        reco_noms_scale.push_back( (TH1D*) reco_noms[i]->Clone( Form("unf_scale_%1.0f%1.0f", jetEdges[i], jetEdges[i+1] ) ) );

        reco_tests[i]->Scale( 1.0 / ( reco_tests[i]->Integral(0, reco_tests[i]->GetNbinsX() + 1) ) );
        
		reco_noms[i]->Scale( 1.0 / ( reco_noms[i]->Integral(0, reco_noms[i]->GetNbinsX() + 1) ) );
		reco_IP2s[i]->Scale( 1.0 / ( reco_IP2s[i]->Integral(0, reco_IP2s[i]->GetNbinsX() + 1) ) );
		reco_IP6s[i]->Scale( 1.0 / ( reco_IP6s[i]->Integral(0, reco_IP6s[i]->GetNbinsX() + 1) ) );
		reco_TSs[i]->Scale( 1.0 / ( reco_TSs[i]->Integral(0, reco_TSs[i]->GetNbinsX() + 1) ) );
		reco_TUs[i]->Scale( 1.0 / ( reco_TUs[i]->Integral(0, reco_TUs[i]->GetNbinsX() + 1) ) );
		reco_HC50s[i]->Scale( 1.0 / ( reco_HC50s[i]->Integral(0, reco_HC50s[i]->GetNbinsX() + 1) ) );
		reco_DSs[i]->Scale( 1.0 / ( reco_DSs[i]->Integral(0, reco_DSs[i]->GetNbinsX() + 1) ) );
		reco_GSs[i]->Scale( 1.0 / ( reco_GSs[i]->Integral(0, reco_GSs[i]->GetNbinsX() + 1) ) );


		// clone unfolded projections after normalization to simplify things (instead of cloning before and normalizing both e.g. reco_noms and reco_noms_copy)
		reco_noms_copy.push_back( (TH1D*) reco_noms[i]->Clone( Form("nom_copy_%1.0f%1.0f", jetEdges[i], jetEdges[i+1] ) ) );
		reco_IP2s_copy.push_back( (TH1D*) reco_IP2s[i]->Clone( Form("IP2_copy_%1.0f%1.0f", jetEdges[i], jetEdges[i+1] ) ) );
		reco_IP6s_copy.push_back( (TH1D*) reco_IP6s[i]->Clone( Form("IP6_copy_%1.0f%1.0f", jetEdges[i], jetEdges[i+1] ) ) );
		reco_TSs_copy.push_back( (TH1D*) reco_TSs[i]->Clone( Form("TS_copy_%1.0f%1.0f", jetEdges[i], jetEdges[i+1] ) ) );
		reco_TUs_copy.push_back( (TH1D*) reco_TUs[i]->Clone( Form("TU_copy_%1.0f%1.0f", jetEdges[i], jetEdges[i+1] ) ) );
		reco_HC50s_copy.push_back( (TH1D*) reco_HC50s[i]->Clone( Form("HC50_copy_%1.0f%1.0f", jetEdges[i], jetEdges[i+1] ) ) );
		reco_DSs_copy.push_back( (TH1D*) reco_DSs[i]->Clone( Form("DS_copy_%1.0f%1.0f", jetEdges[i], jetEdges[i+1] ) ) );
		reco_GSs_copy.push_back( (TH1D*) reco_GSs[i]->Clone( Form("GS_copy_%1.0f%1.0f", jetEdges[i], jetEdges[i+1] ) ) );



		// divide by nominal to see effect as percent of nominal
		reco_IP2s[i]->Divide(reco_noms[i]);
		reco_IP6s[i]->Divide(reco_noms[i]);
		reco_TSs[i]->Divide(reco_noms[i]);
		reco_TUs[i]->Divide(reco_noms[i]);
		reco_HC50s[i]->Divide(reco_noms[i]);
		reco_DSs[i]->Divide(reco_noms[i]);
		reco_GSs[i]->Divide(reco_noms[i]);




        //cout << "DEBUG: Where is the break? normalization " << same_1iters[i]->Integral( 0 , same_1iters[i]->GetNbinsX() + 1 ) << "\n";

        for(int jj = 0; jj < N_iters; jj++){
            recos_same[jj][i]->Scale( 1.0 / recos_same[jj][i]->Integral( 0 , recos_same[jj][i]->GetNbinsX() + 1 ) );
            recos_opp[jj][i]->Scale( 1.0 / recos_opp[jj][i]->Integral( 0 , recos_opp[jj][i]->GetNbinsX() + 1 ) );
        }

        // to always be able to divide by the 4iteration same-side histogram
        //TH1D* norm = (TH1D*) recos_same[3][i]->Clone( Form( "nominal_4iters_%i", i ) );

        B_gen[i]->Scale( 1.0 / B_gen[i]->Integral( 0 , B_gen[i]->GetNbinsX() + 1 ) );
        A_gen[i]->Scale( 1.0 / A_gen[i]->Integral( 0 , A_gen[i]->GetNbinsX() + 1 ) );

        for(int jj = 0; jj < N_iters; jj++){
            recos_same[jj][i]->Divide( B_gen[i] );
            recos_opp[jj][i]->Divide( A_gen[i] );
        }

	}





	reco_nom1D_2025->Scale( 1.0 / reco_nom1D_2025->Integral() );
	reco_nom1D_2530->Scale( 1.0 / reco_nom1D_2530->Integral() );
	reco_nom1D_3040->Scale( 1.0 / reco_nom1D_3040->Integral() );

	reco_h7smear1D_2025->Scale( 1.0 / reco_h7smear1D_2025->Integral() );
	reco_h7smear1D_2530->Scale( 1.0 / reco_h7smear1D_2530->Integral() );
	reco_h7smear1D_3040->Scale( 1.0 / reco_h7smear1D_3040->Integral() );

	reco_p8smear1D_2025->Scale( 1.0 / reco_p8smear1D_2025->Integral() );
	reco_p8smear1D_2530->Scale( 1.0 / reco_p8smear1D_2530->Integral() );
	reco_p8smear1D_3040->Scale( 1.0 / reco_p8smear1D_3040->Integral() );



    reco_P8s_copy.push_back( (TH1D*) reco_p8smear1D_2025->Clone( "P8_copy_2025" ) );
    reco_P8s_copy.push_back( (TH1D*) reco_p8smear1D_2530->Clone( "P8_copy_2530" ) );
    reco_P8s_copy.push_back( (TH1D*) reco_p8smear1D_3040->Clone( "P8_copy_3040" ) );

    reco_H7s_copy.push_back( (TH1D*) reco_h7smear1D_2025->Clone( "H7_copy_2025" ) );
    reco_H7s_copy.push_back( (TH1D*) reco_h7smear1D_2530->Clone( "H7_copy_2530" ) );
    reco_H7s_copy.push_back( (TH1D*) reco_h7smear1D_3040->Clone( "H7_copy_3040" ) );




	reco_h7smear1D_2025->Divide( reco_nom1D_2025 );
	reco_h7smear1D_2530->Divide( reco_nom1D_2530 );
	reco_h7smear1D_3040->Divide( reco_nom1D_3040 );

	reco_p8smear1D_2025->Divide( reco_nom1D_2025 );
	reco_p8smear1D_2530->Divide( reco_nom1D_2530 );
	reco_p8smear1D_3040->Divide( reco_nom1D_3040 );




	//(f) since RooUnfold doesn't handle the statistical uncertainties correctly due to the misses, we scale them up by-hand afterwards
	TH1D* scalefactors = (TH1D*) file3->Get("hratio"); // each bin of this histogram will be a scale factor for the errors of each pT slice of your unfolded 2D hist
	// histograms in file3: "hratio" or "efficiency"----Isaac's code suggests the proper histogram for scalefactors is "hratio"

	// it will scale all bins of jet charge equally, so e.g. if my jet charge distribution for pt = 20-30 GeV has GetBinError() of 0.3, 1.0, 0.09, 0.07, 0.07, 0.1 for Q = -3:-2, -2:-1, -1:0, 0:1, 1:2, 2:3, and the scaling (bin content of the scaling histogram) is 1.5 for 20-30 then you get 0.45, 0.15, 0.135, 0.105, 0.105, 0.15
	for(int i = 0; i < nBins; i++){
        double scaling = -1;
        if(i == 1){scaling = scalefactors->GetBinContent(scalefactors->GetXaxis()->FindBin(20));}
        if(i == 2){scaling = scalefactors->GetBinContent(scalefactors->GetXaxis()->FindBin(25));} //if 15 is the lowest bin
        if(i == 3){scaling = max( scalefactors->GetBinContent(scalefactors->GetXaxis()->FindBin(30)), scalefactors->GetBinContent(scalefactors->GetXaxis()->FindBin(35)) );}


		for(int j = 1; j <= reco_noms[i]->GetNbinsX(); j++){
			double binerr = reco_noms[i]->GetBinError(j);
			reco_noms[i]->SetBinError(j, (double) binerr*scaling);


            reco_noms_scale[i]->SetBinError( j, (double) reco_noms_scale[i]->GetBinError(j)*scaling );
		}
	}
	// reco_slices[i] now have the fully corrected data for each pT range, and we're done!



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



//cout << "setting histogram content to % difference compared to nominal\n";

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



        for(int jj = 0; jj < N_iters; jj++){
            recos_same[jj][i]->PutStats(stats);
//            recos_same[jj][i]->Sumw2(0);

            recos_opp[jj][i]->PutStats(stats);
//            recos_opp[jj][i]->Sumw2(0);
        }


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


    // 2/7/23 envelope procedure taken from Isaac as well

    //taking maximum envelopes!
    TH2D* env_HC = new TH2D("env_HC", "", nBins_Q, -Q_max, Q_max, 15, 5, 80);
    vector<TH1D*> env_HCs = Projection2D(env_HC, nBins, jetEdges, "x");
    TH2D* env_un = new TH2D("env_un", "", nBins_Q, -Q_max, Q_max, 15, 5, 80);
    vector<TH1D*> env_uns = Projection2D(env_un, nBins, jetEdges, "x");
    TH2D* net = new TH2D("net", "", nBins_Q, -Q_max, Q_max, 15, 5, 80);
    vector<TH1D*> nets = Projection2D(net, nBins, jetEdges, "x");
    TH2D* stat = new TH2D("stat", "", nBins_Q, -Q_max, Q_max, 15, 5, 80);
    vector<TH1D*> stat_s = Projection2D(stat, nBins, jetEdges, "x");

    vector<vector< double> > syst_errs2D;


    for (int i = 0; i < nBins; ++ i) {
        vector<double> syst_errs1D;
        reco_noms_copy[i]->Scale(1/(double)reco_noms_copy[i]->Integral());
        //for each mass bin, determine the largest effect in each category
        for (int j = 1; j < nBins_Q + 1; ++ j) {
            //hadronic correction envelope
            double hcs [1] = {reco_HC50s[i]->GetBinContent(j)};
            set<double> hc_sort (hcs, hcs+1);
            set<double>::iterator hc = hc_sort.end(); hc --;
            double hc_envelope = *hc;
            env_HCs[i]->SetBinContent(j, hc_envelope);
            //a clunky way of doing this for both current ranges of 1D mass responses: 20-30, 30-45:
            TH1D* reco_h7_1D = new TH1D("reco_h7_1D","", nBins_Q, -Q_max, Q_max);
            TH1D* reco_p8_1D = new TH1D("reco_p8_1D","", nBins_Q, -Q_max, Q_max);
            if (i == 0) {reco_h7_1D = reco_h7smear1D_2025; reco_p8_1D = reco_p8smear1D_2025;}
            if (i == 1) {reco_h7_1D = reco_h7smear1D_2530; reco_p8_1D = reco_p8smear1D_2530;}
            if (i == 2) {reco_h7_1D = reco_h7smear1D_3040; reco_p8_1D = reco_p8smear1D_3040;}
            //unfolding envelope - using an ordered set here to automatically get the largest value
            double uns [6] = {reco_DSs[i]->GetBinContent(j), reco_GSs[i]->GetBinContent(j),  reco_h7_1D->GetBinContent(j), reco_p8_1D->GetBinContent(j), reco_IP2s[i]->GetBinContent(j), reco_IP6s[i]->GetBinContent(j)};
            set<double> un_sort (uns, uns+6);
            set<double>::iterator un = un_sort.end(); un --;
            double un_envelope = *un;
            env_uns[i]->SetBinContent(j, un_envelope);
            //total uncertainty = TU + TS + un envelope + hc envelope
            double square = pow(hc_envelope,2) + pow(un_envelope,2) + pow(reco_TUs[i]->GetBinContent(j),2) + pow(reco_TSs[i]->GetBinContent(j),2);
            
            nets[i]->SetBinContent(j,sqrt(square));
            syst_errs1D.push_back(nets[i]->GetBinContent(j));

            /*
            if( reco_noms[i]->GetBinContent(j) < 0.00001 ){
                stat_s[i]->SetBinContent( j, 0 );
            }
            else{*/
                stat_s[i]->SetBinContent( j, reco_noms[i]->GetBinError(j) / reco_noms[i]->GetBinContent(j) );
//            }
        }

        syst_errs2D.push_back(syst_errs1D);
    }

    //unfolded result with systematic errors!
    vector<TH1D*> w_systs;
    for (int i = 0; i < nBins; ++ i) {
        w_systs.push_back( (TH1D*) reco_noms_copy[i]->Clone( ("w_systs_" + to_string(i) ).c_str() ) );
        for (int j = 1; j < nBins_Q + 1; ++ j) {
            w_systs[i]->SetBinError( j, syst_errs2D[i][j-1]*w_systs[i]->GetBinContent(j) );
        }
    }




	TFile *fout = new TFile( ("./out/unfold/unfolded_R"+ radius + "_k" + kappa + ".root").c_str() , "RECREATE");
	fout->cd();

	//cout << "output file name: " << ("./out/unfold/unfolded_R"+ radius + "_k" + kappa + ".root").c_str() << "\n";

    cout << "output file name: " << fout->GetName() << "\n";


	for(int i = 0; i < nBins; i++){
        reco_noms_scale[i]->Write();

        reco_tests[i]->Write();
        
		reco_noms[i]->Write();
		reco_IP2s[i]->Write();
		reco_IP6s[i]->Write();
		reco_TSs[i]->Write();
		reco_TUs[i]->Write();
		reco_HC50s[i]->Write();
		reco_DSs[i]->Write();
		reco_GSs[i]->Write();

		reco_noms_copy[i]->Write();
		reco_IP2s_copy[i]->Write();
		reco_IP6s_copy[i]->Write();
		reco_TSs_copy[i]->Write();
		reco_TUs_copy[i]->Write();
		reco_HC50s_copy[i]->Write();
		reco_DSs_copy[i]->Write();
		reco_GSs_copy[i]->Write();


        if( i < nBins - 1 ){
            reco_P8s_copy[i]->Write();
            reco_H7s_copy[i]->Write();
        }

        env_HCs[i]->Write(); //systematic envelopes
        env_uns[i]->Write();
        nets[i]->Write(); // full systematic envelope

        //reco_noms_copy[i]->Write(); //unfolded datapoints
        //dats[i]->Write(); //raw datapoints
        w_systs[i]->Write(); //unfolded data with errors equal to net systematic uncertainty

        stat_s[i]->Write();

	}

	reco_h7smear1D_2025->Write();
	reco_h7smear1D_2530->Write();
	reco_h7smear1D_3040->Write();

	reco_p8smear1D_2025->Write();
	reco_p8smear1D_2530->Write();
	reco_p8smear1D_3040->Write();


	res->Write();
	dat_spectrum->Write();
	reco_nom->Write();

    //cout << "write out objects\n";

	//dat_1D_1520->Write();
	dat_1D_2025->Write();
	dat_1D_2530->Write();
	dat_1D_3040->Write();

    //cout << "what did i add that breaks things??\n";


    //qres2D_part->Write();
    //qres2D_det->Write();

    //cout << "is it the 2D distributions??\n";

    //ptres2D_part->Write();
    //ptres2D_det->Write();

    cout << "Done.\n";

    for(int i = 0; i < nBins; i++){
        for(int jj = 0; jj < N_iters; jj++){
            recos_same[jj][i]->Write();
            recos_opp[jj][i]->Write();
        }

        A_det[i]->Write();
        A_gen[i]->Write();
        B_det[i]->Write();
        B_gen[i]->Write();
    }


    A_qvpt_p->Write();
    A_qvpt_d->Write();
    B_qvpt_p->Write();
    B_qvpt_d->Write();



	fout->Close();

	return 0;
}
