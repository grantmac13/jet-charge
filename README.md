# jet-charge
Performs the jet charge analysis: for pp data and PYTHIA 6, PYTHIA 6 + GEANT
data analyzed using src/data.cxx, then histograms filled from trees filled in data.cxx using macros/hists.cxx
pythia analyzed using src/sim.cxx, histograms filled with entries in trees filled in sim.cxx,
        response matrices for unfolding filled in macros/response.cxx from these same trees filled in sim.cxx,
        and then the outputs from response.cxx are hadded to unfold using unfold.cxx, and closure.cxx takes the output from response.cxx and performs a closure test
