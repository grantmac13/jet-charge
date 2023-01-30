// Standard library
#include <math.h>
#include <iostream>
#include <fstream>

// ROOT Library
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TCut.h>
#include <TPad.h>
#include "TLatex.h"
#include <TColor.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>




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
//            proj1Ds.push_back(hist2D->ProjectionY((hist2D->GetName() + axis + low_rough + high_rough).c_str(), ranges[i], ranges[i+1] - 1, "e"));
            proj1Ds.push_back(hist2D->ProjectionY((hist2D->GetName() + axis + low_rough + high_rough).c_str(), hist2D->GetXaxis()->FindBin(ranges[i]), hist2D->GetXaxis()->FindBin(ranges[i+1]) - 1, "e"));
        }
        else {
            std::cerr << "Improper axis given for projections. Exiting." << std::endl; exit(1);
        }
        proj1Ds[i]->SetTitle("");
        // debug

    }
    return proj1Ds;
}



void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}
