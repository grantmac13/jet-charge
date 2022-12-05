# jet-charge
Performs the jet charge analysis: for pp data and PYTHIA 6, PYTHIA 6 + GEANT
data analyzed using src/data.cxx, then histograms filled from trees filled in data.cxx using macros/hists.cxx
pythia analyzed using src/sim.cxx, histograms filled with entries in trees filled in sim.cxx,
        response matrices for unfolding filled in macros/response.cxx from these same trees filled in sim.cxx,
        and then the outputs from response.cxx are hadded to unfold using unfold.cxx, and closure.cxx takes the output from response.cxx and performs a closure test


To perform initial jet-finding for a given radius R:
        submit/analyze.csh data(sim) trigger species R full/charged (matched)

To then create histograms using these jets and a particular kappa value k:
        csh submit/macro_submit.csh hists sim species R matched(unmatched)(anything else for data) k
        
        hadd into single file for data (pre-unfolding), and unmatched pythia to use in macros that plot distributions
        
        need unmatched pythia-6 and then pythia-8 and herwig-7 histograms hadded into a single file each to do Q-smearing systematic uncertainty
        

To create response matrices:
        csh submit/macro_submit.csh response sim species R matched k
        
        hadd into a file to use in macros, unfolding procedure
        
        
To scale statistical uncertainties properly:
        csh submit/macro_submit.csh stat_err_scaling sim species R matched k
        
        
To perform unfolding:
        csh submit/macro_submit.csh unfold data species R whatever k


