import os
import sys
import json
import ROOT

# load custom analyzer with collections of tracks
analyzer_path = os.path.join(os.path.dirname(__file__), '../../analyzers', 'analyzer_tracktools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

analyzer_path = os.path.join(os.path.dirname(__file__), '../../analyzers', 'analyzer_recotomctools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# main analyzer class
class RDFanalysis():

    def analysers(df):

        # initialization
        dfout = df

        # getting input data and aliases
        dfout = dfout.Alias("Particle", "MCParticles")
        dfout = dfout.Alias("ReconstructedParticles", "RecoParticles")
        dfout = dfout.Alias("ParticleIDs", "ParticleID")
        dfout = dfout.Alias("EFlowTrack", "Tracks")
        dfout = dfout.Define("EFlowTrack_1", "TrackTools::getModifiedTrackStates(_Tracks_trackStates)")
        dfout = dfout.Define("EFlowTrack_2", "1.0 / ReconstructedParticle::get_p(ReconstructedParticles)")
        dfout = dfout.Alias("Reco2TrackLinks", "_RecoParticles_tracks")

        # store all gen-particle info
        dfout = dfout.Define("GenParticle_pdgId", "MCParticle::get_pdg(Particle)")
        dfout = dfout.Define("GenParticle_genStatus", "MCParticle::get_genStatus(Particle)")
        dfout = dfout.Define("GenParticle_px", "MCParticle::get_px(Particle)")
        dfout = dfout.Define("GenParticle_py", "MCParticle::get_py(Particle)")
        dfout = dfout.Define("GenParticle_pz", "MCParticle::get_pz(Particle)")
        dfout = dfout.Define("GenParticle_x", "ROOT::VecOps::RVec<float> out; for (const auto& t : Particle) out.emplace_back(t.vertex.x); return out;")
        dfout = dfout.Define("GenParticle_y", "ROOT::VecOps::RVec<float> out; for (const auto& t : Particle) out.emplace_back(t.vertex.y); return out;")
        dfout = dfout.Define("GenParticle_z", "ROOT::VecOps::RVec<float> out; for (const auto& t : Particle) out.emplace_back(t.vertex.z); return out;")

        # store all track info
        dfout = dfout.Define("Tracks_d0", "ROOT::VecOps::RVec<float> out; for (const auto& t : EFlowTrack_1) out.emplace_back(t.D0); return out;")
        dfout = dfout.Define("Tracks_z0", "ROOT::VecOps::RVec<float> out; for (const auto& t : EFlowTrack_1) out.emplace_back(t.Z0); return out;")
        dfout = dfout.Define("Tracks_phi0", "ROOT::VecOps::RVec<float> out; for (const auto& t : EFlowTrack_1) out.emplace_back(t.phi); return out;")
        dfout = dfout.Define("Tracks_omega", "ROOT::VecOps::RVec<float> out; for (const auto& t : EFlowTrack_1) out.emplace_back(t.omega); return out;")
        dfout = dfout.Define("Tracks_tanlambda", "ROOT::VecOps::RVec<float> out; for (const auto& t : EFlowTrack_1) out.emplace_back(t.tanLambda); return out;")        

        # store links
        dfout = dfout.Define("TrackToMCMap", "RecoToMCTools::makeTrackToMCMapping(EFlowTrack, _trackMCLink_to, _trackMCLink_from)")

        return dfout

    def output():

        branchList = []

        # gen-level stuff
        branchList += [
            'GenParticle_pdgId',
            'GenParticle_genStatus',
            'GenParticle_px',
            'GenParticle_py',
            'GenParticle_pz',
            'GenParticle_x',
            'GenParticle_y',
            'GenParticle_z',

            'Tracks_d0',
            'Tracks_z0',
            'Tracks_phi0',
            'Tracks_omega',
            'Tracks_tanlambda',

            'TrackToMCMap'
        ]

        return branchList    
