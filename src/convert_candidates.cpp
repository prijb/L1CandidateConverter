//This script tries to convert the L1 objects into a generic candidate object
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <string>

#include "Candidates.h" 

auto preprocess(ROOT::RDF::RNode df){

    using Vbool = ROOT::RVec<bool>;
    using Vint = ROOT::RVec<int>;
    using Vfloat = ROOT::RVec<float>;

    /*
    struct Candidates{
        Vfloat pt;
        Vfloat eta;
        Vfloat phi;
        Vint type;
        //Exclusive to muons
        Vint charge;
        Vint dxy;
        //Exclusive to EGammas
        Vfloat iso;
        //The following variables use associations between objects
        Vfloat assoc_jet_pt;
        Vfloat assoc_rel_jet_pt;
        Vfloat assoc_muon_pt;
        Vfloat assoc_rel_muon_pt;
        Vfloat assoc_eg_pt;
        Vfloat assoc_rel_eg_pt;
        //To track original class
        Vint type_original;
    };
    */
     
    //Match obj1 to obj2 (This determines the purity of the object)
    auto getObjMatch = [](ROOT::RVec<float> Obj1_eta, ROOT::RVec<float> Obj1_phi, 
    ROOT::RVec<float> Obj2_eta, ROOT::RVec<float> Obj2_phi){
        ROOT::RVec<int> isMatched;
        std::vector<int> matchedIndices;
        for (int i=0; i<Obj1_eta.size(); i++){
            float mindR = 999;
            int matchedIndex = -1;
            for (int j=0; j<Obj2_eta.size(); j++){
                //Skip already matched objects
                if(std::find(matchedIndices.begin(), matchedIndices.end(), j) != matchedIndices.end()) continue;
                float deta = Obj1_eta[i] - Obj2_eta[j];
                float dphi = TVector2::Phi_mpi_pi(Obj1_phi[i] - Obj2_phi[j]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if((dR < 0.4)&&(dR < mindR)){
                    mindR = dR;
                    matchedIndex = j;
                }
            }
            if(matchedIndex != -1){
                isMatched.push_back(1);
                matchedIndices.push_back(matchedIndex);
            }
            else{
                isMatched.push_back(0);
            }
        }
        return isMatched;
    };

    //Make a list of candidates from the event starting with muons, EG and then jets
    //Type: Muon=0, EGamma=1, LightJet=2, HeavyJet=3, Fake=-1
    auto makeCandidates = [](ROOT::RVec<float> Muon_pt, ROOT::RVec<float> Muon_eta, ROOT::RVec<float> Muon_phi, ROOT::RVec<int> Muon_charge, ROOT::RVec<int> Muon_dxy, ROOT::RVec<int> Muon_isGenmatched,
    ROOT::RVec<float> EGamma_pt, ROOT::RVec<float> EGamma_eta, ROOT::RVec<float> EGamma_phi, ROOT::RVec<int> EGamma_iso, ROOT::RVec<int> EGamma_isGenmatched,
    ROOT::RVec<float> Jet_pt, ROOT::RVec<float> Jet_eta, ROOT::RVec<float> Jet_phi, ROOT::RVec<int> Jet_isGenmatched, ROOT::RVec<int> Jet_partonFlav){
        Candidates candidates;
        //Start with muons
        for(int iMu=0; iMu<Muon_pt.size(); iMu++){
            //Calculate associated variables
            //Muon activity
            TLorentzVector mu_activity;
            mu_activity.SetPtEtaPhiM(0., 0., 0., 0.);
            for(int jMu=0; jMu<Muon_pt.size(); jMu++){
                if(iMu == jMu) continue;
                //Only consider muons within dR < 0.4
                float deta = Muon_eta[iMu] - Muon_eta[jMu];
                float dphi = TVector2::Phi_mpi_pi(Muon_phi[iMu] - Muon_phi[jMu]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if(dR < 0.4){
                    TLorentzVector mu;
                    mu.SetPtEtaPhiM(Muon_pt[jMu], Muon_eta[jMu], Muon_phi[jMu], 0.1055);
                    mu_activity += mu;
                }
            }
            //EGamma activity
            TLorentzVector eg_activity;
            eg_activity.SetPtEtaPhiM(0., 0., 0., 0.);
            for(int jEGamma=0; jEGamma<EGamma_pt.size(); jEGamma++){
                //Only consider EGamma within dR < 0.4
                float deta = Muon_eta[iMu] - EGamma_eta[jEGamma];
                float dphi = TVector2::Phi_mpi_pi(Muon_phi[iMu] - EGamma_phi[jEGamma]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if(dR < 0.4){
                    TLorentzVector eg;
                    eg.SetPtEtaPhiM(EGamma_pt[jEGamma], EGamma_eta[jEGamma], EGamma_phi[jEGamma], 0.);
                    eg_activity += eg;
                }
            }
            //Jet activity
            TLorentzVector jet_activity;
            jet_activity.SetPtEtaPhiM(0., 0., 0., 0.);
            for(int jJet=0; jJet<Jet_pt.size(); jJet++){
                //Only consider jets within dR < 0.4
                float deta = Muon_eta[iMu] - Jet_eta[jJet];
                float dphi = TVector2::Phi_mpi_pi(Muon_phi[iMu] - Jet_phi[jJet]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if(dR < 0.4){
                    TLorentzVector jet;
                    jet.SetPtEtaPhiM(Jet_pt[jJet], Jet_eta[jJet], Jet_phi[jJet], 0.);
                    jet_activity += jet;
                }
            }
            //Get the type
            int type = -1;
            int charge = 0;
            if(Muon_isGenmatched[iMu] == 1) type = 0;
            //Convert hwCharge
            if(Muon_charge[iMu] == 1) charge = 1;
            else if(Muon_charge[iMu] == 0) charge = -1;
            //Fill candidate variables
            candidates.pt.push_back(Muon_pt[iMu]);
            candidates.eta.push_back(Muon_eta[iMu]);
            candidates.phi.push_back(Muon_phi[iMu]);
            candidates.type.push_back(type);
            candidates.charge.push_back(charge);
            candidates.dxy.push_back(Muon_dxy[iMu]);
            candidates.iso.push_back(-999);
            candidates.assoc_muon_pt.push_back(mu_activity.Pt());
            candidates.assoc_rel_muon_pt.push_back(mu_activity.Pt()/Muon_pt[iMu]);
            candidates.assoc_eg_pt.push_back(eg_activity.Pt());
            candidates.assoc_rel_eg_pt.push_back(eg_activity.Pt()/Muon_pt[iMu]);
            candidates.assoc_jet_pt.push_back(jet_activity.Pt());
            candidates.assoc_rel_jet_pt.push_back(jet_activity.Pt()/Muon_pt[iMu]);
            candidates.type_original.push_back(0);
        }
        //Next EGamma
        for(int iEGamma=0; iEGamma<EGamma_pt.size(); iEGamma++){
            //Calculate associated variables
            //Muon activity
            TLorentzVector mu_activity;
            mu_activity.SetPtEtaPhiM(0., 0., 0., 0.);
            for(int jMu=0; jMu<Muon_pt.size(); jMu++){
                //Only consider muons within dR < 0.4
                float deta = EGamma_eta[iEGamma] - Muon_eta[jMu];
                float dphi = TVector2::Phi_mpi_pi(EGamma_phi[iEGamma] - Muon_phi[jMu]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if(dR < 0.4){
                    TLorentzVector mu;
                    mu.SetPtEtaPhiM(Muon_pt[jMu], Muon_eta[jMu], Muon_phi[jMu], 0.1055);
                    mu_activity += mu;
                }
            }
            //EGamma activity
            TLorentzVector eg_activity;
            eg_activity.SetPtEtaPhiM(0., 0., 0., 0.);
            for(int jEGamma=0; jEGamma<EGamma_pt.size(); jEGamma++){
                if(iEGamma == jEGamma) continue;
                //Only consider EGamma within dR < 0.4
                float deta = EGamma_eta[iEGamma] - EGamma_eta[jEGamma];
                float dphi = TVector2::Phi_mpi_pi(EGamma_phi[iEGamma] - EGamma_phi[jEGamma]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if(dR < 0.4){
                    TLorentzVector eg;
                    eg.SetPtEtaPhiM(EGamma_pt[jEGamma], EGamma_eta[jEGamma], EGamma_phi[jEGamma], 0.);
                    eg_activity += eg;
                }
            }
            //Jet activity
            TLorentzVector jet_activity;
            jet_activity.SetPtEtaPhiM(0., 0., 0., 0.);
            for(int jJet=0; jJet<Jet_pt.size(); jJet++){
                //Only consider jets within dR < 0.4
                float deta = EGamma_eta[iEGamma] - Jet_eta[jJet];
                float dphi = TVector2::Phi_mpi_pi(EGamma_phi[iEGamma] - Jet_phi[jJet]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if(dR < 0.4){
                    TLorentzVector jet;
                    jet.SetPtEtaPhiM(Jet_pt[jJet], Jet_eta[jJet], Jet_phi[jJet], 0.);
                    jet_activity += jet;
                }
            }
            //Get the type
            int type = -1;
            if(EGamma_isGenmatched[iEGamma] == 1) type = 1;
            //Fill candidate variables
            candidates.pt.push_back(EGamma_pt[iEGamma]);
            candidates.eta.push_back(EGamma_eta[iEGamma]);
            candidates.phi.push_back(EGamma_phi[iEGamma]);
            candidates.type.push_back(type);
            candidates.charge.push_back(0);
            candidates.dxy.push_back(-999);
            candidates.iso.push_back(EGamma_iso[iEGamma]);
            candidates.assoc_muon_pt.push_back(mu_activity.Pt());
            candidates.assoc_rel_muon_pt.push_back(mu_activity.Pt()/EGamma_pt[iEGamma]);
            candidates.assoc_eg_pt.push_back(eg_activity.Pt());
            candidates.assoc_rel_eg_pt.push_back(eg_activity.Pt()/EGamma_pt[iEGamma]);
            candidates.assoc_jet_pt.push_back(jet_activity.Pt());
            candidates.assoc_rel_jet_pt.push_back(jet_activity.Pt()/EGamma_pt[iEGamma]);
            candidates.type_original.push_back(1);
        }
        //Finally jets
        for(int iJet=0; iJet<Jet_pt.size(); iJet++){
            //Calculate associated variables
            //Muon activity
            TLorentzVector mu_activity;
            mu_activity.SetPtEtaPhiM(0., 0., 0., 0.);
            for(int jMu=0; jMu<Muon_pt.size(); jMu++){
                //Only consider muons within dR < 0.4
                float deta = Jet_eta[iJet] - Muon_eta[jMu];
                float dphi = TVector2::Phi_mpi_pi(Jet_phi[iJet] - Muon_phi[jMu]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if(dR < 0.4){
                    TLorentzVector mu;
                    mu.SetPtEtaPhiM(Muon_pt[jMu], Muon_eta[jMu], Muon_phi[jMu], 0.1055);
                    mu_activity += mu;
                }
            }
            //EGamma activity
            TLorentzVector eg_activity;
            eg_activity.SetPtEtaPhiM(0., 0., 0., 0.);
            for(int jEGamma=0; jEGamma<EGamma_pt.size(); jEGamma++){
                //Only consider EGamma within dR < 0.4
                float deta = Jet_eta[iJet] - EGamma_eta[jEGamma];
                float dphi = TVector2::Phi_mpi_pi(Jet_phi[iJet] - EGamma_phi[jEGamma]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if(dR < 0.4){
                    TLorentzVector eg;
                    eg.SetPtEtaPhiM(EGamma_pt[jEGamma], EGamma_eta[jEGamma], EGamma_phi[jEGamma], 0.);
                    eg_activity += eg;
                }
            }
            //Jet activity
            TLorentzVector jet_activity;
            jet_activity.SetPtEtaPhiM(0., 0., 0., 0.);
            for(int jJet=0; jJet<Jet_pt.size(); jJet++){
                if(iJet == jJet) continue;
                //Only consider jets within dR < 0.4
                float deta = Jet_eta[iJet] - Jet_eta[jJet];
                float dphi = TVector2::Phi_mpi_pi(Jet_phi[iJet] - Jet_phi[jJet]);
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if(dR < 0.4){
                    TLorentzVector jet;
                    jet.SetPtEtaPhiM(Jet_pt[jJet], Jet_eta[jJet], Jet_phi[jJet], 0.);
                    jet_activity += jet;
                }
            }
            //Get the type
            int type = -1;
            int type_original = -1;
            if((Jet_isGenmatched[iJet] == 1)&&(Jet_partonFlav[iJet] != 0)){
                if(Jet_partonFlav[iJet] == 5) type = 3;
                else type = 2;
            }
            //For the original type
            if(Jet_partonFlav[iJet] == 5) type_original = 3;
            else type_original = 2;
            //Fill candidate variables
            candidates.pt.push_back(Jet_pt[iJet]);
            candidates.eta.push_back(Jet_eta[iJet]);
            candidates.phi.push_back(Jet_phi[iJet]);
            candidates.type.push_back(type);
            candidates.charge.push_back(0);
            candidates.dxy.push_back(-999);
            candidates.iso.push_back(-999);
            candidates.assoc_muon_pt.push_back(mu_activity.Pt());
            candidates.assoc_rel_muon_pt.push_back(mu_activity.Pt()/Jet_pt[iJet]);
            candidates.assoc_eg_pt.push_back(eg_activity.Pt());
            candidates.assoc_rel_eg_pt.push_back(eg_activity.Pt()/Jet_pt[iJet]);
            candidates.assoc_jet_pt.push_back(jet_activity.Pt());
            candidates.assoc_rel_jet_pt.push_back(jet_activity.Pt()/Jet_pt[iJet]);
            candidates.type_original.push_back(type_original);
        }

        return candidates;
    };

    //Filter dataframe to meet selection stream requirements
    auto df_trig = df.Filter("nJet > 1")
    .Filter("Jet_pt[0] > 30")
    .Filter("Jet_pt[1] > 30");

    //Get list of gen particles by type
    auto df_gen = df_trig.Define("GenMuon_pt", "GenPart_pt[abs(GenPart_pdgId) == 13]")
    .Define("GenMuon_eta", "GenPart_eta[abs(GenPart_pdgId) == 13]")
    .Define("GenMuon_phi", "GenPart_phi[abs(GenPart_pdgId) == 13]")
    .Define("GenEGamma_pt", "GenPart_pt[abs(GenPart_pdgId) == 11 || abs(GenPart_pdgId) == 22]")
    .Define("GenEGamma_eta", "GenPart_eta[abs(GenPart_pdgId) == 11 || abs(GenPart_pdgId) == 22]")
    .Define("GenEGamma_phi", "GenPart_phi[abs(GenPart_pdgId) == 11 || abs(GenPart_pdgId) == 22]")
    //Perform acceptance cuts on all collections
    .Define("GenMuon_cut", "GenMuon_pt > 4 && abs(GenMuon_eta) < 2.4")
    .Define("Muon_cut", "Muon_pt > 4 && abs(Muon_eta) < 2.4")
    .Define("GenEGamma_cut", "GenEGamma_pt > 10 && abs(GenEGamma_eta) < 3.0")
    .Define("EGamma_cut", "EGamma_pt > 10 && abs(EGamma_eta) < 3.0")
    .Define("GenJet_cut", "GenJet_pt > 10 && abs(GenJet_eta) < 3.0")
    .Define("Jet_cut", "Jet_pt > 10 && abs(Jet_eta) < 3.0")
    .Redefine("GenMuon_pt", "GenMuon_pt[GenMuon_cut]")
    .Redefine("GenMuon_eta", "GenMuon_eta[GenMuon_cut]")
    .Redefine("GenMuon_phi", "GenMuon_phi[GenMuon_cut]")
    .Redefine("GenEGamma_pt", "GenEGamma_pt[GenEGamma_cut]")
    .Redefine("GenEGamma_eta", "GenEGamma_eta[GenEGamma_cut]")
    .Redefine("GenEGamma_phi", "GenEGamma_phi[GenEGamma_cut]")
    .Redefine("GenJet_pt", "GenJet_pt[GenJet_cut]")
    .Redefine("GenJet_eta", "GenJet_eta[GenJet_cut]")
    .Redefine("GenJet_phi", "GenJet_phi[GenJet_cut]")
    .Redefine("Muon_pt", "Muon_pt[Muon_cut]")
    .Redefine("Muon_eta", "Muon_eta[Muon_cut]")
    .Redefine("Muon_phi", "Muon_phi[Muon_cut]")
    .Redefine("Muon_etaAtVtx", "Muon_etaAtVtx[Muon_cut]")
    .Redefine("Muon_phiAtVtx", "Muon_phiAtVtx[Muon_cut]")
    .Redefine("EGamma_pt", "EGamma_pt[EGamma_cut]")
    .Redefine("EGamma_eta", "EGamma_eta[EGamma_cut]")
    .Redefine("EGamma_phi", "EGamma_phi[EGamma_cut]")
    .Redefine("Jet_pt", "Jet_pt[Jet_cut]")
    .Redefine("Jet_eta", "Jet_eta[Jet_cut]")
    .Redefine("Jet_phi", "Jet_phi[Jet_cut]")
    .Redefine("Jet_partonFlav", "Jet_partonFlav[Jet_cut]");

    //Perform genmatching
    auto df_match = df_gen.Define("Muon_GenMuon_Match", getObjMatch, {"Muon_etaAtVtx", "Muon_phiAtVtx", "GenMuon_eta", "GenMuon_phi"})
    .Define("Muon_GenEGamma_Match", getObjMatch, {"Muon_etaAtVtx", "Muon_phiAtVtx", "GenEGamma_eta", "GenEGamma_phi"})
    .Define("Muon_GenJet_Match", getObjMatch, {"Muon_etaAtVtx", "Muon_phiAtVtx", "GenJet_eta", "GenJet_phi"})
    .Define("EGamma_GenEGamma_Match", getObjMatch, {"EGamma_eta", "EGamma_phi", "GenEGamma_eta", "GenEGamma_phi"})
    .Define("EGamma_GenJet_Match", getObjMatch, {"EGamma_eta", "EGamma_phi", "GenJet_eta", "GenJet_phi"})
    .Define("Jet_GenJet_Match", getObjMatch, {"Jet_eta", "Jet_phi", "GenJet_eta", "GenJet_phi"})
    .Define("Jet_Genmatched_partonFlav", "Jet_partonFlav[Jet_GenJet_Match == 1]")
    .Define("Jet_Genmatched_Heavy", "Jet_Genmatched_partonFlav == 5")
    .Define("Jet_Genmatched_Light", "(Jet_Genmatched_partonFlav!=5)&&(Jet_Genmatched_partonFlav!=0)");

    //Make candidates
    auto df_cands = df_match.Define("Candidate", makeCandidates, {"Muon_pt", "Muon_etaAtVtx", "Muon_phiAtVtx", "Muon_hwCharge", "Muon_hwDXY", "Muon_GenMuon_Match", 
    "EGamma_pt", "EGamma_eta", "EGamma_phi", "EGamma_Iso", "EGamma_GenEGamma_Match",
    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_GenJet_Match", "Jet_partonFlav"})
    .Define("nCandidate", "Candidate.pt.size()")
    .Define("Candidate_pt", "Candidate.pt")
    .Define("Candidate_eta", "Candidate.eta")
    .Define("Candidate_phi", "Candidate.phi")
    .Define("Candidate_type", "Candidate.type")
    .Define("Candidate_charge", "Candidate.charge")
    .Define("Candidate_dxy", "Candidate.dxy")
    .Define("Candidate_iso", "Candidate.iso")
    .Define("Candidate_assoc_jet_pt", "Candidate.assoc_jet_pt")
    .Define("Candidate_assoc_rel_jet_pt", "Candidate.assoc_rel_jet_pt")
    .Define("Candidate_assoc_muon_pt", "Candidate.assoc_muon_pt")
    .Define("Candidate_assoc_rel_muon_pt", "Candidate.assoc_rel_muon_pt")
    .Define("Candidate_assoc_eg_pt", "Candidate.assoc_eg_pt")
    .Define("Candidate_assoc_rel_eg_pt", "Candidate.assoc_rel_eg_pt")
    .Define("Candidate_type_original", "Candidate.type_original");

    return df_cands;
}

//Main function
int main(int argc, char** argv){
    char* input_dir = argv[1];
    float weight = std::stof(argv[2]);
    char* output_file = argv[3];

    gSystem->Load("build/CandidatesDict.so");

    TChain chain;
    chain.Add((std::string(input_dir)+std::string("/output_1.root/scNtuplizer/Events")).c_str());

    ROOT::EnableImplicitMT(8);
    ROOT::RDataFrame df(chain);

    auto entries = df.Count();
    int num_entries = *entries;

    auto df_processed = preprocess(df);
    auto df_processed_weighted = df_processed.Define("weight", [&num_entries, &weight]() { return (weight)/float(num_entries);});

    //Save processed file as a snapshot
    std::cout << "Saving snapshot" << std::endl;
    ROOT::EnableImplicitMT(0);

    df_processed_weighted.Snapshot("Events", output_file, {"weight", "nCandidate", "Candidate_pt", "Candidate_eta", "Candidate_phi", "Candidate_type", "Candidate_charge", "Candidate_dxy", "Candidate_iso", "Candidate_assoc_jet_pt", "Candidate_assoc_rel_jet_pt", "Candidate_assoc_muon_pt", "Candidate_assoc_rel_muon_pt", "Candidate_assoc_eg_pt", "Candidate_assoc_rel_eg_pt", "Candidate_type_original"});

    return 0;
}