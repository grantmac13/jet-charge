// From Isaac Mooney-- adapted by Grant McNamara, March 2022
// performs closure tests for 2D and 1D responses


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

using namespace std;

int main (int argc, const char** argv) {
	//intro
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//basic argument checking.
	if (argc < 2) {
		cerr << "Should be at least two arguments: jet radius and kappa. Received "
		     << argc-1 << ". Exiting." << endl;
		exit(1);
	}

	cout << "A\n";
	//argv[1] should be the jet radius e.g. "04".
	string radius = (string) argv[1];
	string kappa = (string) argv[2];


	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TH3::SetDefaultSumw2();


	//input files
	// May 25th was day when we used mass cut, pdg mass process to create plots of mass closure to compare to Isaac to debug my closure process
	// May 26th was day after walking through with Isaac, confirming that mass cut, pdg mass process was not used anymore to create charge closure
	TFile *f = new TFile(("~/ppJCandBF/out/sim/response/p6_R" + radius + "_k" + kappa + "_syst_May26.root").c_str(),"READ");
	//responses
	RooUnfoldResponse *res1D_q_A = (RooUnfoldResponse*) f->Get("sampleA_q_response");
	RooUnfoldResponse *res1D_pt_A = (RooUnfoldResponse*) f->Get("sampleA_pt_response");
	RooUnfoldResponse *res2D_A = (RooUnfoldResponse*) f->Get("sampleA_q_pt_response");

	RooUnfoldResponse *res1D_m_A = (RooUnfoldResponse*) f->Get("sampleA_m_response");
	RooUnfoldResponse *res2D_mpt_A = (RooUnfoldResponse*) f->Get("sampleA_m_pt_response");


	RooUnfoldResponse *res1D_q_B = (RooUnfoldResponse*) f->Get("sampleB_q_response");
	RooUnfoldResponse *res1D_pt_B = (RooUnfoldResponse*) f->Get("sampleB_pt_response");
	RooUnfoldResponse *res2D_B = (RooUnfoldResponse*) f->Get("sampleB_q_pt_response");

	RooUnfoldResponse *res1D_m_B = (RooUnfoldResponse*) f->Get("sampleB_m_response");
	RooUnfoldResponse *res2D_mpt_B = (RooUnfoldResponse*) f->Get("sampleB_m_pt_response");


	//pseudo-data
	//pop A
	TH2D* q_pt_det_A = (TH2D*) f->Get("sampleA_q_pt_det");
	TH2D* q_pt_gen_A = (TH2D*) f->Get("sampleA_q_pt_gen");
	TH1D* q_det_A = (TH1D*) f->Get("sampleA_q_det");
	TH1D* q_gen_A = (TH1D*) f->Get("sampleA_q_gen");
	TH1D* pt_det_A = (TH1D*) f->Get("sampleA_pt_det");
	TH1D* pt_gen_A = (TH1D*) f->Get("sampleA_pt_gen");

	TH2D* m_pt_det_A = (TH2D*) f->Get("sampleA_m_pt_det");
	TH2D* m_pt_gen_A = (TH2D*) f->Get("sampleA_m_pt_gen");
	TH1D* m_det_A = (TH1D*) f->Get("sampleA_m_det");
	TH1D* m_gen_A = (TH1D*) f->Get("sampleA_m_gen");


	//pop B
	TH2D* q_pt_det_B = (TH2D*) f->Get("sampleB_q_pt_det");
	TH2D* q_pt_gen_B = (TH2D*) f->Get("sampleB_q_pt_gen");
	TH1D* q_det_B = (TH1D*) f->Get("sampleB_q_det");
	TH1D* q_gen_B = (TH1D*) f->Get("sampleB_q_gen");
	TH1D* pt_det_B = (TH1D*) f->Get("sampleB_pt_det");
	TH1D* pt_gen_B = (TH1D*) f->Get("sampleB_pt_gen");

	TH2D* m_pt_det_B = (TH2D*) f->Get("sampleB_m_pt_det");
	TH2D* m_pt_gen_B = (TH2D*) f->Get("sampleB_m_pt_gen");
	TH1D* m_det_B = (TH1D*) f->Get("sampleB_m_det");
	TH1D* m_gen_B = (TH1D*) f->Get("sampleB_m_gen");


	cout << "B" << endl;
	/////////////////
	// 3/21/22: test through 2 iterations, 3,4,5,6,7,8,9,10 iterations

	RooUnfoldBayes *unfold_opp_q1D = new RooUnfoldBayes(res1D_q_A, q_det_B, 4, false, "unfold_opp_q1D","");
	RooUnfoldBayes *unfold_opp_pt1D = new RooUnfoldBayes(res1D_pt_A, pt_det_B, 4, false, "unfold_opp_pt1D","");
 	RooUnfoldBayes *unfold_opp_2D = new RooUnfoldBayes(res2D_A, q_pt_det_B, 4, false, "unfold_opp_2D","");

	RooUnfoldBayes *unfold_opp_m1D = new RooUnfoldBayes(res1D_m_A, m_det_B, 6, false, "unfold_opp_m1D","");

	// COME BACK TO THIS GRANT
 	RooUnfoldBayes *unfold_opp_mpt_2D = new RooUnfoldBayes(res2D_mpt_A, m_pt_det_B, 6, false, "unfold_opp_mpt_2D","");
	cout << "j" << endl;


	// for testing number of iterations
 	RooUnfoldBayes *unfold_opp_2D_2 = new RooUnfoldBayes(res2D_A, q_pt_det_B, 2, false, "unfold_opp_2D_2","");
 	RooUnfoldBayes *unfold_opp_2D_3 = new RooUnfoldBayes(res2D_A, q_pt_det_B, 3, false, "unfold_opp_2D_3","");
 	RooUnfoldBayes *unfold_opp_2D_4 = new RooUnfoldBayes(res2D_A, q_pt_det_B, 4, false, "unfold_opp_2D_4","");
 	RooUnfoldBayes *unfold_opp_2D_5 = new RooUnfoldBayes(res2D_A, q_pt_det_B, 5, false, "unfold_opp_2D_5","");
 	RooUnfoldBayes *unfold_opp_2D_6 = new RooUnfoldBayes(res2D_A, q_pt_det_B, 6, false, "unfold_opp_2D_6","");
 	RooUnfoldBayes *unfold_opp_2D_7 = new RooUnfoldBayes(res2D_A, q_pt_det_B, 7, false, "unfold_opp_2D_7","");
 	RooUnfoldBayes *unfold_opp_2D_8 = new RooUnfoldBayes(res2D_A, q_pt_det_B, 8, false, "unfold_opp_2D_8","");
 	RooUnfoldBayes *unfold_opp_2D_9 = new RooUnfoldBayes(res2D_A, q_pt_det_B, 9, false, "unfold_opp_2D_9","");
 	RooUnfoldBayes *unfold_opp_2D_10 = new RooUnfoldBayes(res2D_A, q_pt_det_B, 10, false, "unfold_opp_2D_10","");


	cout << "k" << endl;

	RooUnfoldBayes *unfold_same_q1D = new RooUnfoldBayes(res1D_q_A, q_det_A, 4, false, "unfold_same_q1D","");
	RooUnfoldBayes *unfold_same_pt1D = new RooUnfoldBayes(res1D_pt_A, pt_det_A, 4, false, "unfold_same_pt1D","");
	RooUnfoldBayes *unfold_same_2D = new RooUnfoldBayes(res2D_A, q_pt_det_A, 4, false, "unfold_same_2D","");


	RooUnfoldBayes *unfold_same_m1D = new RooUnfoldBayes(res1D_m_A, m_det_A, 6, false, "unfold_same_m1D","");
	RooUnfoldBayes *unfold_same_mpt_2D = new RooUnfoldBayes(res2D_mpt_A, m_pt_det_A, 6, false, "unfold_same_mpt_2D","");


	cout << "C" << endl;

	TH1D *reco_opp_q1D = (TH1D*) unfold_opp_q1D->Hreco((RooUnfold::ErrorTreatment) 3);

	cout << "l" << endl;

	TH1D *reco_opp_pt1D = (TH1D*) unfold_opp_pt1D->Hreco((RooUnfold::ErrorTreatment) 3);

	cout << "i" << endl;

	TH2D *reco_opp_2D = (TH2D*) unfold_opp_2D->Hreco((RooUnfold::ErrorTreatment) 3);

	cout << "d" << endl;


	TH1D *reco_opp_m1D = (TH1D*) unfold_opp_m1D->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D *reco_opp_mpt_2D = (TH2D*) unfold_opp_mpt_2D->Hreco((RooUnfold::ErrorTreatment) 3);

	cout << "e" << endl;



	TH2D* reco_opp_2D_2 = (TH2D*) unfold_opp_2D_2->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D* reco_opp_2D_3 = (TH2D*) unfold_opp_2D_3->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D* reco_opp_2D_4 = (TH2D*) unfold_opp_2D_4->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D* reco_opp_2D_5 = (TH2D*) unfold_opp_2D_5->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D* reco_opp_2D_6 = (TH2D*) unfold_opp_2D_6->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D* reco_opp_2D_7 = (TH2D*) unfold_opp_2D_7->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D* reco_opp_2D_8 = (TH2D*) unfold_opp_2D_8->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D* reco_opp_2D_9 = (TH2D*) unfold_opp_2D_9->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D* reco_opp_2D_10 = (TH2D*) unfold_opp_2D_10->Hreco((RooUnfold::ErrorTreatment) 3);

	cout << "f" << endl;


	TH1D *reco_same_q1D = (TH1D*) unfold_same_q1D->Hreco((RooUnfold::ErrorTreatment) 3);
	TH1D *reco_same_pt1D = (TH1D*) unfold_same_pt1D->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D *reco_same_2D = (TH2D*) unfold_same_2D->Hreco((RooUnfold::ErrorTreatment) 3);


	TH1D *reco_same_m1D = (TH1D*) unfold_same_m1D->Hreco((RooUnfold::ErrorTreatment) 3);
	TH2D *reco_same_mpt_2D = (TH2D*) unfold_same_mpt_2D->Hreco((RooUnfold::ErrorTreatment) 3);



	reco_opp_q1D->Divide(q_gen_A);
	reco_opp_pt1D->Divide(pt_gen_A);

	reco_opp_m1D->Divide(m_gen_A);



	TH1D* opp_pt2D = (TH1D*) reco_opp_2D->ProjectionY("opp_pt2D");
	vector<TH1D*> opp_q_projs =
		{(TH1D*) reco_opp_2D->ProjectionX("opp_q1015", reco_opp_2D->GetYaxis()->FindBin(10), reco_opp_2D->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D->ProjectionX("opp_q1520", reco_opp_2D->GetYaxis()->FindBin(15), reco_opp_2D->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D->ProjectionX("opp_q2030", reco_opp_2D->GetYaxis()->FindBin(20), reco_opp_2D->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D->ProjectionX("opp_q3040", reco_opp_2D->GetYaxis()->FindBin(30), reco_opp_2D->GetYaxis()->FindBin(40)-1)
		};

	vector<TH1D*> opp_q_projs_2 =
		{(TH1D*) reco_opp_2D_2->ProjectionX("opp_q1015_2iter", reco_opp_2D_2->GetYaxis()->FindBin(10), reco_opp_2D_2->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D_2->ProjectionX("opp_q1520_2iter", reco_opp_2D_2->GetYaxis()->FindBin(15), reco_opp_2D_2->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D_2->ProjectionX("opp_q2030_2iter", reco_opp_2D_2->GetYaxis()->FindBin(20), reco_opp_2D_2->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D_2->ProjectionX("opp_q3040_2iter", reco_opp_2D_2->GetYaxis()->FindBin(30), reco_opp_2D_2->GetYaxis()->FindBin(40)-1)
		};
	vector<TH1D*> opp_q_projs_3 =
		{(TH1D*) reco_opp_2D_3->ProjectionX("opp_q1015_3iter", reco_opp_2D_3->GetYaxis()->FindBin(10), reco_opp_2D_3->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D_3->ProjectionX("opp_q1520_3iter", reco_opp_2D_3->GetYaxis()->FindBin(15), reco_opp_2D_3->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D_3->ProjectionX("opp_q2030_3iter", reco_opp_2D_3->GetYaxis()->FindBin(20), reco_opp_2D_3->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D_3->ProjectionX("opp_q3040_3iter", reco_opp_2D_3->GetYaxis()->FindBin(30), reco_opp_2D_3->GetYaxis()->FindBin(40)-1)
		};
	vector<TH1D*> opp_q_projs_4 =
		{(TH1D*) reco_opp_2D_4->ProjectionX("opp_q1015_4iter", reco_opp_2D_4->GetYaxis()->FindBin(10), reco_opp_2D_4->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D_4->ProjectionX("opp_q1520_4iter", reco_opp_2D_4->GetYaxis()->FindBin(15), reco_opp_2D_4->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D_4->ProjectionX("opp_q2030_4iter", reco_opp_2D_4->GetYaxis()->FindBin(20), reco_opp_2D_4->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D_4->ProjectionX("opp_q3040_4iter", reco_opp_2D_4->GetYaxis()->FindBin(30), reco_opp_2D_4->GetYaxis()->FindBin(40)-1)
		};
	vector<TH1D*> opp_q_projs_5 =
		{(TH1D*) reco_opp_2D_5->ProjectionX("opp_q1015_5iter", reco_opp_2D_5->GetYaxis()->FindBin(10), reco_opp_2D_5->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D_5->ProjectionX("opp_q1520_5iter", reco_opp_2D_5->GetYaxis()->FindBin(15), reco_opp_2D_5->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D_5->ProjectionX("opp_q2030_5iter", reco_opp_2D_5->GetYaxis()->FindBin(20), reco_opp_2D_5->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D_5->ProjectionX("opp_q3040_5iter", reco_opp_2D_5->GetYaxis()->FindBin(30), reco_opp_2D_5->GetYaxis()->FindBin(40)-1)
		};
	vector<TH1D*> opp_q_projs_6 =
		{(TH1D*) reco_opp_2D_6->ProjectionX("opp_q1015_6iter", reco_opp_2D_6->GetYaxis()->FindBin(10), reco_opp_2D_6->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D_6->ProjectionX("opp_q1520_6iter", reco_opp_2D_6->GetYaxis()->FindBin(15), reco_opp_2D_6->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D_6->ProjectionX("opp_q2030_6iter", reco_opp_2D_6->GetYaxis()->FindBin(20), reco_opp_2D_6->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D_6->ProjectionX("opp_q3040_6iter", reco_opp_2D_6->GetYaxis()->FindBin(30), reco_opp_2D_6->GetYaxis()->FindBin(40)-1)
		};
	vector<TH1D*> opp_q_projs_7 =
		{(TH1D*) reco_opp_2D_7->ProjectionX("opp_q1015_7iter", reco_opp_2D_7->GetYaxis()->FindBin(10), reco_opp_2D_7->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D_7->ProjectionX("opp_q1520_7iter", reco_opp_2D_7->GetYaxis()->FindBin(15), reco_opp_2D_7->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D_7->ProjectionX("opp_q2030_7iter", reco_opp_2D_7->GetYaxis()->FindBin(20), reco_opp_2D_7->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D_7->ProjectionX("opp_q3040_7iter", reco_opp_2D_7->GetYaxis()->FindBin(30), reco_opp_2D_7->GetYaxis()->FindBin(40)-1)
		};
	vector<TH1D*> opp_q_projs_8 =
		{(TH1D*) reco_opp_2D_8->ProjectionX("opp_q1015_8iter", reco_opp_2D_8->GetYaxis()->FindBin(10), reco_opp_2D_8->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D_8->ProjectionX("opp_q1520_8iter", reco_opp_2D_8->GetYaxis()->FindBin(15), reco_opp_2D_8->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D_8->ProjectionX("opp_q2030_8iter", reco_opp_2D_8->GetYaxis()->FindBin(20), reco_opp_2D_8->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D_8->ProjectionX("opp_q3040_8iter", reco_opp_2D_8->GetYaxis()->FindBin(30), reco_opp_2D_8->GetYaxis()->FindBin(40)-1)
		};
	vector<TH1D*> opp_q_projs_9 =
		{(TH1D*) reco_opp_2D_9->ProjectionX("opp_q1015_9iter", reco_opp_2D_9->GetYaxis()->FindBin(10), reco_opp_2D_9->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D_9->ProjectionX("opp_q1520_9iter", reco_opp_2D_9->GetYaxis()->FindBin(15), reco_opp_2D_9->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D_9->ProjectionX("opp_q2030_9iter", reco_opp_2D_9->GetYaxis()->FindBin(20), reco_opp_2D_9->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D_9->ProjectionX("opp_q3040_9iter", reco_opp_2D_9->GetYaxis()->FindBin(30), reco_opp_2D_9->GetYaxis()->FindBin(40)-1)
		};
	vector<TH1D*> opp_q_projs_10 =
		{(TH1D*) reco_opp_2D_10->ProjectionX("opp_q1015_10iter", reco_opp_2D_10->GetYaxis()->FindBin(10), reco_opp_2D_10->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_2D_10->ProjectionX("opp_q1520_10iter", reco_opp_2D_10->GetYaxis()->FindBin(15), reco_opp_2D_10->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_opp_2D_10->ProjectionX("opp_q2030_10iter", reco_opp_2D_10->GetYaxis()->FindBin(20), reco_opp_2D_10->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_2D_10->ProjectionX("opp_q3040_10iter", reco_opp_2D_10->GetYaxis()->FindBin(30), reco_opp_2D_10->GetYaxis()->FindBin(40)-1)
		};



	vector<TH1D*> opp_m_projs =
		{(TH1D*) reco_opp_mpt_2D->ProjectionX("opp_m1015", reco_opp_mpt_2D->GetYaxis()->FindBin(10), reco_opp_mpt_2D->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_opp_mpt_2D->ProjectionX("opp_m2025", reco_opp_mpt_2D->GetYaxis()->FindBin(20), reco_opp_mpt_2D->GetYaxis()->FindBin(25)-1),
		 (TH1D*) reco_opp_mpt_2D->ProjectionX("opp_m2530", reco_opp_mpt_2D->GetYaxis()->FindBin(25), reco_opp_mpt_2D->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_opp_mpt_2D->ProjectionX("opp_m3040", reco_opp_mpt_2D->GetYaxis()->FindBin(30), reco_opp_mpt_2D->GetYaxis()->FindBin(40)-1)
		};


	reco_same_q1D->Divide(q_gen_A);
	reco_same_pt1D->Divide(pt_gen_A);
	TH1D* same_pt2D = (TH1D*) reco_same_2D->ProjectionY("same_pt2D");

	reco_same_m1D->Divide(m_gen_A);
	TH1D* same_mpt2D = (TH1D*) reco_same_mpt_2D->ProjectionY("same_mpt2D");


	vector<TH1D*> same_q_projs =
		{(TH1D*) reco_same_2D->ProjectionX("same_q1015", reco_same_2D->GetYaxis()->FindBin(10), reco_same_2D->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_same_2D->ProjectionX("same_q1520", reco_same_2D->GetYaxis()->FindBin(15), reco_same_2D->GetYaxis()->FindBin(20)-1),
		 (TH1D*) reco_same_2D->ProjectionX("same_q2030", reco_same_2D->GetYaxis()->FindBin(20), reco_same_2D->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_same_2D->ProjectionX("same_q3040", reco_same_2D->GetYaxis()->FindBin(30), reco_same_2D->GetYaxis()->FindBin(40)-1)
		};

	vector<TH1D*> gen_q_A_projs =
		{(TH1D*) q_pt_gen_A->ProjectionX("gen_q1015", q_pt_gen_A->GetYaxis()->FindBin(10), q_pt_gen_A->GetYaxis()->FindBin(15)-1),
		 (TH1D*) q_pt_gen_A->ProjectionX("gen_q1520", q_pt_gen_A->GetYaxis()->FindBin(15), q_pt_gen_A->GetYaxis()->FindBin(20)-1),
		 (TH1D*) q_pt_gen_A->ProjectionX("gen_q2030", q_pt_gen_A->GetYaxis()->FindBin(20), q_pt_gen_A->GetYaxis()->FindBin(30)-1),
		 (TH1D*) q_pt_gen_A->ProjectionX("gen_q3040", q_pt_gen_A->GetYaxis()->FindBin(30), q_pt_gen_A->GetYaxis()->FindBin(40)-1)
		};


	vector<TH1D*> same_m_projs =
		{(TH1D*) reco_same_mpt_2D->ProjectionX("same_m1015", reco_same_mpt_2D->GetYaxis()->FindBin(10), reco_same_mpt_2D->GetYaxis()->FindBin(15)-1),
		 (TH1D*) reco_same_mpt_2D->ProjectionX("same_m2025", reco_same_mpt_2D->GetYaxis()->FindBin(20), reco_same_mpt_2D->GetYaxis()->FindBin(25)-1),
		 (TH1D*) reco_same_mpt_2D->ProjectionX("same_m2530", reco_same_mpt_2D->GetYaxis()->FindBin(25), reco_same_mpt_2D->GetYaxis()->FindBin(30)-1),
		 (TH1D*) reco_same_mpt_2D->ProjectionX("same_m3040", reco_same_mpt_2D->GetYaxis()->FindBin(30), reco_same_mpt_2D->GetYaxis()->FindBin(40)-1)
		};

	vector<TH1D*> gen_m_A_projs =
		{(TH1D*) m_pt_gen_A->ProjectionX("gen_m1015", m_pt_gen_A->GetYaxis()->FindBin(10), m_pt_gen_A->GetYaxis()->FindBin(15)-1),
		 (TH1D*) m_pt_gen_A->ProjectionX("gen_m2025", m_pt_gen_A->GetYaxis()->FindBin(20), m_pt_gen_A->GetYaxis()->FindBin(25)-1),
		 (TH1D*) m_pt_gen_A->ProjectionX("gen_m2530", m_pt_gen_A->GetYaxis()->FindBin(25), m_pt_gen_A->GetYaxis()->FindBin(30)-1),
		 (TH1D*) m_pt_gen_A->ProjectionX("gen_m3040", m_pt_gen_A->GetYaxis()->FindBin(30), m_pt_gen_A->GetYaxis()->FindBin(40)-1)
		};



	const int nBins = 4; //10-15, 15-20, 20-30, 30-40
	//		for mass: 10-15, 20-25, 25-30, 30-40
	for (int i = 0; i < nBins; ++ i) {
		opp_q_projs[i]->Scale( 1.0 / ( opp_q_projs[i]->Integral() ) );
		same_q_projs[i]->Scale( 1.0 / ( same_q_projs[i]->Integral() ) );
		gen_q_A_projs[i]->Scale( 1.0 / ( gen_q_A_projs[i]->Integral() ) );


		opp_m_projs[i]->Scale( 1.0 / ( opp_m_projs[i]->Integral() ) );
		same_m_projs[i]->Scale( 1.0 / ( same_m_projs[i]->Integral() ) );
		gen_m_A_projs[i]->Scale( 1.0 / ( gen_m_A_projs[i]->Integral() ) );


		opp_q_projs_2[i]->Scale( 1.0 / ( opp_q_projs_2[i]->Integral() ) );
		opp_q_projs_3[i]->Scale( 1.0 / ( opp_q_projs_3[i]->Integral() ) );
		opp_q_projs_4[i]->Scale( 1.0 / ( opp_q_projs_4[i]->Integral() ) );
		opp_q_projs_5[i]->Scale( 1.0 / ( opp_q_projs_5[i]->Integral() ) );
		opp_q_projs_6[i]->Scale( 1.0 / ( opp_q_projs_6[i]->Integral() ) );
		opp_q_projs_7[i]->Scale( 1.0 / ( opp_q_projs_7[i]->Integral() ) );
		opp_q_projs_8[i]->Scale( 1.0 / ( opp_q_projs_8[i]->Integral() ) );
		opp_q_projs_9[i]->Scale( 1.0 / ( opp_q_projs_9[i]->Integral() ) );
		opp_q_projs_10[i]->Scale( 1.0 / ( opp_q_projs_10[i]->Integral() ) );

		opp_q_projs[i]->Divide(gen_q_A_projs[i]); //det B unfolded w/ A, divided by gen A
		same_q_projs[i]->Divide(gen_q_A_projs[i]);//det A unfolded w/ A, divided by gen A


		opp_m_projs[i]->Divide(gen_m_A_projs[i]); //det B unfolded w/ A, divided by gen A
		same_m_projs[i]->Divide(gen_m_A_projs[i]);//det A unfolded w/ A, divided by gen A


		opp_q_projs_2[i]->Divide(gen_q_A_projs[i]);
		opp_q_projs_3[i]->Divide(gen_q_A_projs[i]);
		opp_q_projs_4[i]->Divide(gen_q_A_projs[i]);
		opp_q_projs_5[i]->Divide(gen_q_A_projs[i]);
		opp_q_projs_6[i]->Divide(gen_q_A_projs[i]);
		opp_q_projs_7[i]->Divide(gen_q_A_projs[i]);
		opp_q_projs_8[i]->Divide(gen_q_A_projs[i]);
		opp_q_projs_9[i]->Divide(gen_q_A_projs[i]);
		opp_q_projs_10[i]->Divide(gen_q_A_projs[i]);

	}


	string ftitle = "closure";

	TFile *fout = new TFile(("~/ppJCandBF/out/closure/" + ftitle + "_R" + radius + "_k" + kappa + "_new.root").c_str(),"RECREATE");
	fout->cd();

	reco_opp_q1D->Write();
	reco_opp_pt1D->Write();
	reco_same_q1D->Write();
	reco_same_pt1D->Write();

	reco_opp_m1D->Write();
	reco_same_m1D->Write();

	same_pt2D->Write();
	opp_pt2D->Write();

	for (int i = 0; i < nBins; ++ i) {
		opp_q_projs[i]->Write();
		same_q_projs[i]->Write();
		gen_q_A_projs[i]->Write();

		opp_m_projs[i]->Write();
		same_m_projs[i]->Write();
		gen_m_A_projs[i]->Write();
	}
	for (int i = 0; i < nBins; ++ i) {
		opp_q_projs_2[i]->Write();
		opp_q_projs_3[i]->Write();
		opp_q_projs_4[i]->Write();
		opp_q_projs_5[i]->Write();
		opp_q_projs_6[i]->Write();
		opp_q_projs_7[i]->Write();
		opp_q_projs_8[i]->Write();
		opp_q_projs_9[i]->Write();
		opp_q_projs_10[i]->Write();
	}

	fout->Close();

	return 0;
}
