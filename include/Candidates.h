// Candidates.h
#ifndef CANDIDATES_H
#define CANDIDATES_H

#include <ROOT/RVec.hxx>  // Include the ROOT header for RVec

struct Candidates {
    ROOT::RVec<float> pt;
    ROOT::RVec<float> eta;
    ROOT::RVec<float> phi;
    ROOT::RVec<int> type;
    ROOT::RVec<int> charge;
    ROOT::RVec<int> dxy;
    ROOT::RVec<float> iso;
    ROOT::RVec<float> assoc_jet_pt;
    ROOT::RVec<float> assoc_rel_jet_pt;
    ROOT::RVec<float> assoc_muon_pt;
    ROOT::RVec<float> assoc_rel_muon_pt;
    ROOT::RVec<float> assoc_eg_pt;
    ROOT::RVec<float> assoc_rel_eg_pt;
    ROOT::RVec<int> type_original;
};

#endif // CANDIDATES_H