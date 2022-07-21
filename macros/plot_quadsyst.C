void plot_quadsyst(double R = 0.4, double kappa = 0.0){

	TString rad = Form("%i%i", (int) R, (int) (R*10) );
	TString k = Form("%i%i", (int) kappa, (int) (kappa*10) );
	cout << "converting double to string: radius = " << rad << "\n";
	cout << "converting double to string: kappa = " << k << "\n";
	TFile* fin = new TFile("~/ppJCandBF/out/unfold/unfolded_R" + rad + "_k" + k + ".root", "READ");

	
	const int nJetBins = 3;
	int jetPt[nJetBins+1] = {20, 25, 30, 40};
	//get systematics unfolded histograms
	TH1D* IP2s[nJetBins];
	TH1D* IP6s[nJetBins];
	TH1D* TSs[nJetBins];
	TH1D* TUs[nJetBins];
	TH1D* HC50s[nJetBins];
	TH1D* DSs[nJetBins];
	TH1D* GSs[nJetBins];
	TH1D* h7s[nJetBins];
	TH1D* p8s[nJetBins];


	// histogram to fill with the quadrature sums of the systematics bin-by-bin
	TH1D* noms[nJetBins];


	TCanvas *c[nJetBins];
	TLegend *leg[nJetBins];

	for(int i = 0; i < nJetBins; i++){
		c[i] = new TCanvas(Form("c%i", i), "", 600, 600);
		c[i]->SetLeftMargin(0.15);
		leg[i] = new TLegend(0.28, 0.7, 0.88, 0.88);
		leg[i]->SetBorderSize(0);
		leg[i]->SetNColumns(4);


		IP2s[i] = (TH1D*) fin->Get(Form("unfold_IP2x%i%i", jetPt[i], jetPt[i+1]));
		IP6s[i] = (TH1D*) fin->Get(Form("unfold_IP6x%i%i", jetPt[i], jetPt[i+1]));
		TSs[i] = (TH1D*) fin->Get(Form("unfold_TSx%i%i", jetPt[i], jetPt[i+1]));
		TUs[i] = (TH1D*) fin->Get(Form("unfold_TUx%i%i", jetPt[i], jetPt[i+1]));
		HC50s[i] = (TH1D*) fin->Get(Form("unfold_HC50x%i%i", jetPt[i], jetPt[i+1]));
		DSs[i] = (TH1D*) fin->Get(Form("unfold_DSx%i%i", jetPt[i], jetPt[i+1]));
		GSs[i] = (TH1D*) fin->Get(Form("unfold_GSx%i%i", jetPt[i], jetPt[i+1]));
		h7s[i] = (TH1D*) fin->Get(Form("unfold_h7smear1D_%i%i", jetPt[i], jetPt[i+1]));
		p8s[i] = (TH1D*) fin->Get(Form("unfold_p8smear1D_%i%i", jetPt[i], jetPt[i+1]));



		// get Q bin procedure from macros/hists.cxx or macros/response.cxx
		// or initialize after getting histograms from fin and use the bin info from one of the systematics histograms
		int nBins = IP2s[i]->GetNbinsX();
		double loEdge = IP2s[i]->GetXaxis()->GetBinLowEdge(1);
		double hiEdge = IP2s[i]->GetXaxis()->GetBinLowEdge(nBins+1);
		// loEdge should == -hiEdge
		if(loEdge != -1*hiEdge){cout << "ERROR: ASYMMETRIC BINS SOMEHOW?\n";}

		noms[i] = new TH1D( Form( "total_systs_jetQ_jetpt%i%i", jetPt[i], jetPt[i+1] ), "", nBins, loEdge, hiEdge);
		

		for(int j = 1; j < nBins; j++){
			double quad_sum2 = pow( IP2s[i]->GetBinContent(j), 2 )
					 + pow( IP6s[i]->GetBinContent(j), 2 )
					 + pow( TSs[i]->GetBinContent(j), 2 )
					 + pow( TUs[i]->GetBinContent(j), 2 )
					 + pow( HC50s[i]->GetBinContent(j), 2 )
					 + pow( DSs[i]->GetBinContent(j), 2 )
					 + pow( GSs[i]->GetBinContent(j), 2 )
					 + pow( h7s[i]->GetBinContent(j), 2 )
					 + pow( p8s[i]->GetBinContent(j), 2 );
			
			noms[i]->SetBinContent( j, sqrt(quad_sum2) );
		}

		noms[i]->SetLineStyle(kSolid);
		noms[i]->SetLineColor(kGreen+2);
		noms[i]->SetLineWidth(2);
//		noms[i]->SetFillColor(kGreen+2);
		noms[i]->Sumw2(0);

		noms[i]->GetYaxis()->SetRangeUser(0.0, .999);

		noms[i]->SetStats(0);

		noms[i]->SetTitle( Form("%i < p_{T}^{jet} < %i", jetPt[i], jetPt[i+1]) );
		noms[i]->SetXTitle("Q^{jet}");
		noms[i]->SetYTitle("total % relative error");

		noms[i]->Draw();

//		leg[i]->Draw();

		c[i]->SaveAs(Form("./plots/totalsystematicsJetQ_jetpt%i%i_R" + rad + "_k" + k + ".png", jetPt[i], jetPt[i+1]), "RECREATE");
	}
}
