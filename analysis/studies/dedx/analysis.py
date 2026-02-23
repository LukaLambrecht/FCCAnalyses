import os
import sys
import json
import ROOT

analyzers_dir = '../../analyzers'

# load custom analyzer for particle ID retrieval
analyzer_path = os.path.join(os.path.dirname(__file__), analyzers_dir, 'analyzer_particleid.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with dE/dx tools
analyzer_path = os.path.join(os.path.dirname(__file__), analyzers_dir, 'analyzer_dEdx.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with dE/dx tools
analyzer_path = os.path.join(os.path.dirname(__file__), analyzers_dir, 'analyzer_recotomctools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with collections of tracks
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_tracktools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# main analyzer class
class RDFanalysis():

    def analysers(df):

        # initialization
        dfout = df
        det = 'aleph'

        # automatically determine data type based on the presence of some branches
        # (maybe make more robust later)
        dtype = 'data'
        try:
            if det=='fcc': df.Alias("_", "Particle")
            if det=='aleph': df.Alias("_", "MCParticles")
            dtype = 'sim'
        except: pass

        # for Aleph simulation, the collection of MC particles is called "MCParticles",
        # while for FCC simulation it is called "Particle", so define an alias here.
        # and same for "RecoParticles" vs "ReconstructedParticles".
        # and same for "ParticleID" vs "ParticleIDs" (note: not well defined for FCC).
        if det=='aleph':
            if dtype=='sim': dfout = dfout.Alias("Particle", "MCParticles")
            dfout = dfout.Alias("ReconstructedParticles", "RecoParticles")
            dfout = dfout.Alias("ParticleIDs", "ParticleID")

        # for Aleph simulation, the collections EFlowTrack, EFlowTrack_1 and EFlowTrack_2 do not seem to exist,
        # so we need to alias them with other collections.
        # note: the alias for EFlowTrack_2 has not yet been validated, no guarantee that it is correct.
        if det=='aleph':
            dfout = (
                dfout

                .Alias("EFlowTrack", "Tracks")
                # (must be an object of type ROOT::VecOps::RVec<edm4hep::TrackData>)
                .Alias("EFlowTrack_1", "TrackTools::getModifiedTrackStates(_Tracks_trackStates)")
                # (must be an object of type ROOT::VecOps::RVec<edm4hep::TrackState>)
                .Define("EFlowTrack_2", "1.0 / ReconstructedParticle::get_p(ReconstructedParticles)")
                # (must be an object of type rv::RVec<edm4hep::Quantity>)
            )

        # for Aleph simulation, the link between RecoParticles and Tracks is indirect,
        # via an intermediate collection _RecoParticles_tracks;
        # for FCC simulation, we must define it as a transparent mapping for the syntax
        if det=='aleph': dfout = dfout.Alias("Reco2TrackLinks", "_RecoParticles_tracks")

        # do jet clustering
        dfout = (
            dfout

            # define the momentum, energy, mass and charge of all reconstructed particles.
            .Define("RP_px",          "ReconstructedParticle::get_px(ReconstructedParticles)")
            .Define("RP_py",          "ReconstructedParticle::get_py(ReconstructedParticles)")
            .Define("RP_pz",          "ReconstructedParticle::get_pz(ReconstructedParticles)")
            .Define("RP_e",           "ReconstructedParticle::get_e(ReconstructedParticles)")
            .Define("RP_m",           "ReconstructedParticle::get_mass(ReconstructedParticles)")
            .Define("RP_q",           "ReconstructedParticle::get_charge(ReconstructedParticles)")
            
            # build "pseudo-jets", meaning each particle is converted to a jet
            # consisting of only that one particle (which can then be clustered in the next step)
            .Define("pseudo_jets",    "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
            
            # note: the arguments are (in order):
            # - exclusive
            # - cut
            # - sorted
            # - recombination
            .Define("FCCAnalysesJets_ee_genkt", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")

            # get the jets out of the struct
            .Define("jets_ee_genkt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_genkt)")

            # find the reconstructed particles grouped per jet
            # (output struct is a vector of vectors of ReconstructedParticle objects, one vector of ReconstructedParticles for each jet)
            .Define("jetconstituents_ee_genkt", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_genkt)")
            .Define("JetsConstituents", "JetConstituentsUtils::build_constituents_cluster(ReconstructedParticles, jetconstituents_ee_genkt)")
            .Define("Jets_nConstituents", "JetConstituentsUtils::count_consts(JetsConstituents)")

        )

        # store the PDG ID of every jet constituent
        # (or at least every particle for which it is available, see more details in the helper function)
        if dtype=='sim':
            dfout = (
                dfout
                .Define("TrackToMCMap", "RecoToMCTools::makeTrackToMCMapping(EFlowTrack, _trackMCLink_to, _trackMCLink_from)")
                .Define("JetsConstituents_pdgId", "RecoToMCTools::get_pdgid(JetsConstituents, EFlowTrack, Particle, Reco2TrackLinks, TrackToMCMap)")
            )

        else:
            # set dummies (maybe later find out how to avoid the need for this)
            dfout = (
                dfout
                .Define("JetsConstituents_pdgId", "RecoToMCTools::get_pdgidDummy(JetsConstituents)")
            )

        # rest of the analysis
        dfout = (
            dfout

            # Extract ParticleID types for all jet constituents
            # ParticleID.type legend: 0:Track, 1:Electron, 2:Muon, 3:Track from V0, 
            #                         4:EM, 5:Ecal hadron/residual, 6:Hcal element, 7:Lcal element
            .Define("JetsConstituents_Types", "getParticleIDTypes(ParticleIDs, jetconstituents_ee_genkt)")
            
            # Convert the ParticleID types to binary flags for particle classification
            # Muon: type 2
            # Electron: type 1
            # Gamma: type 4 (EM)
            # Charged hadron: types 0 (Track), 3 (Track from V0)
            # Neutral hadron: types 5 (Ecal hadron/residual), 6 (Hcal element), 7 (Lcal element)
            .Define("JetsConstituents_isMu", "get_isMu_from_type(JetsConstituents_Types)")
            .Define("JetsConstituents_isEl", "get_isEl_from_type(JetsConstituents_Types)")
            .Define("JetsConstituents_isChargedHad", "get_isChargedHad_from_type(JetsConstituents_Types)")
            .Define("JetsConstituents_isGamma", "get_isGamma_from_type(JetsConstituents_Types)")
            .Define("JetsConstituents_isNeutralHad", "get_isNeutralHad_from_type(JetsConstituents_Types)")

            # basic kinematics
            .Define("JetsConstituents_e", "JetConstituentsUtils::get_e(JetsConstituents)")
            .Define("JetsConstituents_pt", "JetConstituentsUtils::get_pt(JetsConstituents)")
            .Define("JetsConstituents_px", "JetConstituentsUtils::get_px(JetsConstituents)")
            .Define("JetsConstituents_py", "JetConstituentsUtils::get_py(JetsConstituents)")
            .Define("JetsConstituents_pz", "JetConstituentsUtils::get_pz(JetsConstituents)")
            .Define("JetsConstituents_theta", "JetConstituentsUtils::get_theta(JetsConstituents)")
            .Define("JetsConstituents_phi", "JetConstituentsUtils::get_phi(JetsConstituents)")
            .Define("JetsConstituents_charge", "JetConstituentsUtils::get_charge(JetsConstituents)")
        )

        # dE/dx values
        dfout = (
            dfout
            .Define("dEdxPadsValue" , "dEdxPads.dQdx.value")
            .Define("dEdxPadsError" , "dEdxPads.dQdx.error")
            .Define("dEdxWiresValue" , "dEdxWires.dQdx.value")
            .Define("dEdxWiresError" , "dEdxWires.dQdx.error")

            # basic variables for pads
            .Define("jet_constituents_dEdx_pads_objs", "dEdxTools::build_constituents_dEdx()(RecoParticles, Reco2TrackLinks.index, dEdxPads, _dEdxPads_track.index, jetconstituents_ee_genkt, false)")
            .Define("jet_constituents_dEdx_pads_dedx", "jet_constituents_dEdx_pads_objs.dedx_constituents")
            .Define("JetsConstituents_dEdx_pads_type", "dEdxTools::get_dEdx_type(jet_constituents_dEdx_pads_dedx)")
            .Define("JetsConstituents_dEdx_pads_value", "dEdxTools::get_dEdx_value(jet_constituents_dEdx_pads_dedx)")
            .Define("JetsConstituents_dEdx_pads_error", "dEdxTools::get_dEdx_error(jet_constituents_dEdx_pads_dedx)")

            # PID hypotheses compatibility for pads
            .Define("jet_constituents_dEdx_pads_pvals", "jet_constituents_dEdx_pads_objs.pid_array_constituents")
            .Define("JetsConstituents_PID_pval_pads_ele", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 0)")
            .Define("JetsConstituents_PID_pval_pads_mu", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 1)")
            .Define("JetsConstituents_PID_pval_pads_pi", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 2)")
            .Define("JetsConstituents_PID_pval_pads_kaon", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 3)")
            .Define("JetsConstituents_PID_pval_pads_proton", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 4)")

            # basic variables for wires
            .Define("jet_constituents_dEdx_wires_objs", "dEdxTools::build_constituents_dEdx()(RecoParticles, Reco2TrackLinks.index, dEdxWires, _dEdxWires_track.index, jetconstituents_ee_genkt, true)")
            .Define("jet_constituents_dEdx_wires_dedx", "jet_constituents_dEdx_wires_objs.dedx_constituents")
            .Define("JetsConstituents_dEdx_wires_type", "dEdxTools::get_dEdx_type(jet_constituents_dEdx_wires_dedx)")
            .Define("JetsConstituents_dEdx_wires_value", "dEdxTools::get_dEdx_value(jet_constituents_dEdx_wires_dedx)")
            .Define("JetsConstituents_dEdx_wires_error", "dEdxTools::get_dEdx_error(jet_constituents_dEdx_wires_dedx)")

            # PID hypotheses compatibility for wires
            .Define("jet_constituents_dEdx_wires_pvals", "jet_constituents_dEdx_wires_objs.pid_array_constituents")
            .Define("JetsConstituents_PID_pval_wires_ele", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 0)")
            .Define("JetsConstituents_PID_pval_wires_mu", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 1)")
            .Define("JetsConstituents_PID_pval_wires_pi", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 2)")
            .Define("JetsConstituents_PID_pval_wires_kaon", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 3)")
            .Define("JetsConstituents_PID_pval_wires_proton", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 4)")        

        )

        return dfout

    def output():

        # define what output to store

        # define output branches
        branchList = []

        # jet-constituent-level variables
        branchList += [

            'JetsConstituents_pdgId',

            'JetsConstituents_e', 'JetsConstituents_pt',
            'JetsConstituents_px', 'JetsConstituents_py', 'JetsConstituents_pz',
            'JetsConstituents_theta', 'JetsConstituents_phi',
            'JetsConstituents_charge',

            'JetsConstituents_dEdx_pads_type',
            'JetsConstituents_dEdx_pads_value',
            'JetsConstituents_dEdx_pads_error',
            'JetsConstituents_dEdx_wires_type',
            'JetsConstituents_dEdx_wires_value',
            'JetsConstituents_dEdx_wires_error',

            'JetsConstituents_PID_pval_pads_ele',
            'JetsConstituents_PID_pval_pads_mu',
            'JetsConstituents_PID_pval_pads_pi',
            'JetsConstituents_PID_pval_pads_kaon',
            'JetsConstituents_PID_pval_pads_proton',

            'JetsConstituents_PID_pval_wires_ele',
            'JetsConstituents_PID_pval_wires_mu',
            'JetsConstituents_PID_pval_wires_pi',
            'JetsConstituents_PID_pval_wires_kaon',
            'JetsConstituents_PID_pval_wires_proton',

            'JetsConstituents_isMu', 
            'JetsConstituents_isEl', 
            'JetsConstituents_isChargedHad',
            'JetsConstituents_isGamma', 
            'JetsConstituents_isNeutralHad',
        ]

        return branchList    
