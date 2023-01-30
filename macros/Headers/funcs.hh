// Isaac Mooney, WSU - June 2019
// functions.hh


#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"


//#include "ktTrackEff.hh"
#include <string>
#include <utility>      // std::pair
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <memory>
#include <set>
#include <fstream>

#ifndef funcs_hh
#define funcs_hh

using std::unordered_map; using std::make_shared; using std::shared_ptr; 

namespace Analysis {


  // IO/OS MANIP Functions from Nick
  
  // Helper to build the TChain, used to decide which input format                                         
  bool HasEnding (std::string const &, std::string const &);

  template <typename T>
  std::set<T> ParseCSV(std::string);

  template<typename T>
  bool CanCast(std::string);

  template<typename T>
  T CastTo(std::string);  

  double radius_str_to_double (std::string radius_str);
  double str_to_double (std::string str);

}

#endif
