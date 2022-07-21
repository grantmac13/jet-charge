void plot_systematics(double R = 0.4, double kappa = 0.0){

	TString rad = Form("%i%i", (int) R, (int) (R*10) );
	TString k = Form("%i%i", (int) kappa, (int) (kappa*10) );
	cout << "converting double to string: radius = " << rad << "\n";
	cout << "converting double to string: kappa = " << k << "\n";
	TFile* fin = new TFile("~/ppJCandBF/out/unfold/unfolded_R" + rad + "_k" + k + ".root", "READ");

	
	const int nJetBins = 3;
	int jetPt[nJetBins+1] = {20, 25, 30, 40};
	//get systematics unfolded histograms
	TH1D* noms[nJetBins];
	TH1D* IP2s[nJetBins];
	TH1D* IP6s[nJetBins];
	TH1D* TSs[nJetBins];
	TH1D* TUs[nJetBins];
	TH1D* HC50s[nJetBins];
	TH1D* DSs[nJetBins];
	TH1D* GSs[nJetBins];
	TH1D* h7s[nJetBins];
	TH1D* p8s[nJetBins];


	TCanvas *c[nJetBins];
	TLegend *leg[nJetBins];

	for(int i = 0; i < nJetBins; i++){
		c[i] = new TCanvas(Form("c%i", i), "", 600, 600);
		c[i]->SetLeftMargin(0.15);
		leg[i] = new TLegend(0.28, 0.7, 0.8, 0.88);
		leg[i]->SetBorderSize(0);
		leg[i]->SetNColumns(4);


		noms[i] = (TH1D*) fin->Get(Form("unfold_nomx%i%i", jetPt[i], jetPt[i+1]));
		IP2s[i] = (TH1D*) fin->Get(Form("unfold_IP2x%i%i", jetPt[i], jetPt[i+1]));
		IP6s[i] = (TH1D*) fin->Get(Form("unfold_IP6x%i%i", jetPt[i], jetPt[i+1]));
		TSs[i] = (TH1D*) fin->Get(Form("unfold_TSx%i%i", jetPt[i], jetPt[i+1]));
		TUs[i] = (TH1D*) fin->Get(Form("unfold_TUx%i%i", jetPt[i], jetPt[i+1]));
		HC50s[i] = (TH1D*) fin->Get(Form("unfold_HC50x%i%i", jetPt[i], jetPt[i+1]));
		DSs[i] = (TH1D*) fin->Get(Form("unfold_DSx%i%i", jetPt[i], jetPt[i+1]));
		GSs[i] = (TH1D*) fin->Get(Form("unfold_GSx%i%i", jetPt[i], jetPt[i+1]));
		h7s[i] = (TH1D*) fin->Get(Form("unfold_h7smear1D_%i%i", jetPt[i], jetPt[i+1]));
		p8s[i] = (TH1D*) fin->Get(Form("unfold_p8smear1D_%i%i", jetPt[i], jetPt[i+1]));


//////////////////////////////// SETTING MARKERS FOR SYSTEMATICS PLOT
		// may be useful to make a function like Isaac's "Prettify1DwLineStyle" to do this more in a more organized way

		IP2s[i]->SetLineStyle(kSolid);
		IP2s[i]->SetLineColor(kBlack);
		IP2s[i]->SetLineWidth(2);
//		IP2s[i]->SetFillColorAlpha(kBlack, 0.1);
//		IP2s[i]->SetFillStyle(3305);
		IP2s[i]->Sumw2(0);


		IP6s[i]->SetLineStyle(kSolid);
		IP6s[i]->SetLineColor(kBlue);
		IP6s[i]->SetLineWidth(2);
//		IP6s[i]->SetFillColorAlpha(kBlue, 0.1);
//		IP6s[i]->SetFillStyle(3395);
		IP6s[i]->Sumw2(0);


		TSs[i]->SetLineStyle(kSolid);
		TSs[i]->SetLineColor(kRed);
		TSs[i]->SetLineWidth(2);
//		TSs[i]->SetFillColorAlpha(kRed, 0.1);
//		TSs[i]->SetFillStyle(3490);
		TSs[i]->Sumw2(0);


		TUs[i]->SetLineStyle(kSolid);
		TUs[i]->SetLineColor(kGreen+2);
		TUs[i]->SetLineWidth(2);
//		TUs[i]->SetFillColorAlpha(kGreen+2, 0.1);
//		TUs[i]->SetFillStyle(3436);
		TUs[i]->Sumw2(0);


		HC50s[i]->SetLineStyle(kSolid);
		HC50s[i]->SetLineColor(kCyan);
		HC50s[i]->SetLineWidth(2);
//		HC50s[i]->SetFillColorAlpha(kCyan, 0.1);
//		HC50s[i]->SetFillStyle(3335);
		HC50s[i]->Sumw2(0);
		

		DSs[i]->SetLineStyle(kSolid);
		DSs[i]->SetLineColor(kOrange+4);
		DSs[i]->SetLineWidth(2);
//		DSs[i]->SetFillColorAlpha(kOrange+4, 0.1);
//		DSs[i]->SetFillStyle(3944);
		DSs[i]->Sumw2(0);


		GSs[i]->SetLineStyle(kSolid);
		GSs[i]->SetLineColor(kViolet-3);
		GSs[i]->SetLineWidth(2);
//		GSs[i]->SetFillColorAlpha(kViolet-3, 0.1);
//		GSs[i]->SetFillStyle(3544);
		GSs[i]->Sumw2(0);


		h7s[i]->SetLineStyle(kSolid);
		h7s[i]->SetLineColor(kGray+1);
		h7s[i]->SetLineWidth(2);
//		h7s[i]->SetFillColorAlpha(kGray+1, 0.1);
//		h7s[i]->SetFillStyle(3690);
		h7s[i]->Sumw2(0);


		p8s[i]->SetLineStyle(kSolid);
		p8s[i]->SetLineColor(kMagenta-9);
		p8s[i]->SetLineWidth(2);
//		p8s[i]->SetFillColorAlpha(kMagenta-9, 0.1);
//		p8s[i]->SetFillStyle(3444);
		p8s[i]->Sumw2(0);




		leg[i]->AddEntry(IP2s[i], "IP2", "f");
		leg[i]->AddEntry(IP6s[i], "IP6", "f");
		leg[i]->AddEntry(TSs[i], "TS", "f");
		leg[i]->AddEntry(TUs[i], "TU", "f");
		leg[i]->AddEntry(HC50s[i], "HC50", "f");
		leg[i]->AddEntry(DSs[i], "DS", "f");
		leg[i]->AddEntry(GSs[i], "GS", "f");
		leg[i]->AddEntry(h7s[i], "H7", "f");
		leg[i]->AddEntry(p8s[i], "P8", "f");


		if(kappa != 0.0){
			IP2s[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
		}
		IP2s[i]->GetYaxis()->SetRangeUser(0.0, .999);

		IP2s[i]->SetStats(0);

		IP2s[i]->SetTitle( Form("%i < p_{T}^{jet} < %i", jetPt[i], jetPt[i+1]) );
		IP2s[i]->SetXTitle("Q^{jet}");
		IP2s[i]->SetYTitle("% relative error");

		IP2s[i]->Draw("L");
		IP6s[i]->Draw("L same");
		TSs[i]->Draw("L same");
		TUs[i]->Draw("L same");
		HC50s[i]->Draw("L same");
		DSs[i]->Draw("L same");
		GSs[i]->Draw("L same");
		h7s[i]->Draw("L same");
		p8s[i]->Draw("L same");

		leg[i]->Draw();

		c[i]->SaveAs(Form("./plots/systematicsJetQ_jetpt%i%i_R" + rad + "_k" + k + ".png", jetPt[i], jetPt[i+1]), "RECREATE");
	}
}
