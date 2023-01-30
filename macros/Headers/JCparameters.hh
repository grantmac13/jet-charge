//
//  JCparameters.hh
//  
//
//  Created by Grant McNamara on 8/6/22.
//  following the advice of Veronica Verkest, this file will contain plot labels, bin information, marker style/color information and name segments stored in const vectors, arrays, etc.
//  this file will be used to make plotting macros much neater as the information will be self-contained within this file. Loops will also become cleaner and, more or less, standardized--may even break certain loops out into standalone functions themselves and store them in another file similarly named to this one.
// use #include "JCparameters.hh" to use in macro, .cxx file

// as of end of day 8/17/22: fitting to all 11 flavors present in pythia-8 simulation (excluding t, tbar) is able to perform the fit, in one bin there is either the quark or the anti-quark present in a meaningful quantity (u, c, g, dbar, sbar, bbar are weighted >= .1% and the rest are < 10^-6)
//      but another pt bin has g ~30%, dbar ~65%, cbar ~5% with everything else < 10^-10
//      and 30-40 bin is u ~22%, g ~24%, dbar ~55%
//
//



#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"


#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>



#ifndef JCparameters_hh
#define JCparameters_hh


namespace jcAnalysis
{

    const int njetbins = 3;
    const double jetPtLo[njetbins] = {20.0, 25.0, 30.0};
    const double jetPtHi[njetbins] = {25.0, 30.0, 40.0};
    const double jetedges[njetbins+1] = {20.0, 25.0, 30.0, 40.0};

    const TString outFile_pt[njetbins] = {"2025", "2530", "3040"}; // for names of histograms, output filenames, etc


    // plotting colors, markers for det-level and part-level data,unfold plots of jet Q
    const int nLevels = 2; // det-, part-level
    const int nTypes = 2; // pythia, pythia+geant and data, unfolded
    const int color_jetQ[nLevels] = {877, 419};
    const Style_t marker[nLevels*nTypes] = {kOpenCircle, kOpenSquare, kFullCircle, kFullSquare};


    // flavors for plotting
    // e.g. draw_unfoldedtemplate.C, draw_template.C, etc.
    // as of 8/22/22 i want to account for only u,d,g,ubar,dbar--and then "other" when i want it will be made separately because it will be a totally free gaussian to account for the minor contribution of the rest of the quarks/anti-quarks
    const int nflavors = 3;
    const TString flavor_names[nflavors] = { "d",    "u",     "g"
                                            // , "dbar", "ubar"
                                        };

//                                        blue, red, teal, purple
    const int mark_colors[nflavors+1] = { 600, 628, 839, 875
                                        //, 600, 628
                                      };

    const int line_styles[nflavors] = { 1, 1, 1
                                       //, 9, 9
                                      };



    // info for when all flavors are required (e.g. template.cxx)
    const int n_template_flavors = 11; // 13;
    const TString template_flavor_names[n_template_flavors] =
                    { "d",     "u",    "s",
                      "c",     "b"/*,  "t"*/,    "g",
                      "dbar", "ubar", "sbar",
                      "cbar", "bbar"/*, "tbar"*/
                    };
    
    // one color per quark/gluon:           d,    u,   s,    c,    b,   g
    const int template_mark_colors[6] = {880+1, 416+3, 600, 632, 800+8, 1};
    // one marker per particle vs anti-particle
    //                          q/g, qbar
    const int mark_styles[2] = { 8,  4};//full circle, open circle

    // systematics naming styles
    const int nsyst = 8;
    const TString systNames[nsyst] = {"nom", "TS", "TU", "HC50", "DS", "GS", "P8", "H7"};

    
    const double other_constraint = 0.000001;

    const int nkaps = 4;
    const double k_vals[nkaps] = {0.0, 0.3, 0.5, 0.7};
    // y-max value for different kappa values for 3-panel plots, etc
    const double jc_padmax[nkaps] = { 0.599, 0.209, 0.329, 0.349};

}

#endif /* JCparameters_hh */
