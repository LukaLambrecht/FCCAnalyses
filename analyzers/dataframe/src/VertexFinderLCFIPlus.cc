// reference: https://arxiv.org/pdf/1506.08371.pdf
// contact: kunal.gautam@cern.ch

#include "FCCAnalyses/VertexFinderLCFIPlus.h"
#include <iostream>

namespace FCCAnalyses{

namespace VertexFinderLCFIPlus{

/*
Global settings and constants
*/

bool debug_me = false;
// if particle masses defined in a dedicated file, call those rather than defining here
const double m_pi = 0.13957039; // pi+- mass [GeV]
const double m_p  = 0.93827208; // p+- mass [GeV]
const double m_e  = 0.00051099; // e+- mass [GeV]


/*
Secondary vertex finding per jet.
*/


ROOT::VecOps::RVec<ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>>
get_SV_jets(
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>& tracksPerJet,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& allTracks,
        const VertexingUtils::FCCAnalysesVertex& PV,
        bool V0_veto,
        double chi2_cut, double invariant_mass_cut, double trackChi2_cut){
    // Base function for doing per-jet secondary vertex finding.
    // Input arguments:
    // - tracksPerJet: tracks grouped per jet.
    //   No further filtering is performed on this collection,
    //   the tracks are assumed to be already selected
    //   (e.g. sufficient quality, incompatible with primary vertex, etc).
    //   Exception: a V0 veto can still be performed if V0_veto is set to true.
    // - allTracks: all tracks in event; not sure why this is needed.
    // - PV: primary vertex, not sure why this is needed.
    // Returns:
    // - a collection of secondary vertices per jet.

    // initialize result
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>> result;

    // loop over jets
    for(unsigned int jetIdx = 0; jetIdx < tracksPerJet.size(); jetIdx++){

        // get tracks for this jet
        ROOT::VecOps::RVec<edm4hep::TrackState> tracksInThisJet = tracksPerJet.at(jetIdx);

        // do V0 veto if requested
        ROOT::VecOps::RVec<edm4hep::TrackState> tracksToUse = V0rejection_tight(tracksInThisJet, PV, V0_veto);

        // find secondary vertices
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> verticesInThisJet;
        verticesInThisJet = findSVfromTracks(tracksToUse, allTracks, PV, chi2_cut, invariant_mass_cut, trackChi2_cut);

        // add to result
        result.push_back(verticesInThisJet);
  }
  return result;    
}


ROOT::VecOps::RVec<ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>>
get_SV_jets(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recoparticles,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& allTracks,
	    const VertexingUtils::FCCAnalysesVertex& PV,
	    ROOT::VecOps::RVec<bool> isInPrimary,
	    ROOT::VecOps::RVec<fastjet::PseudoJet> jets,
	    std::vector<std::vector<int>> jet_constituent_indices,
	    bool V0_rej,
	    double chi2_cut, double invM_cut, double chi2Tr_cut){
    // Same as above (base function) but with more built-in input handling.
    
    // initializations
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> selectedTracksPerJet;

    // retrieve tracks from reco particles & get a vector with their indices in the reco collection
    ROOT::VecOps::RVec<edm4hep::TrackState> tracks = ReconstructedParticle2Track::getRP2TRK( recoparticles, allTracks );
    ROOT::VecOps::RVec<int> reco_ind_tracks = ReconstructedParticle2Track::get_recoindTRK( recoparticles, allTracks );

    // sanity checks
    if(tracks.size() != reco_ind_tracks.size()){
        std::cout<<"ERROR: reco index vector not the same size as no of tracks"<<std::endl;
    }
    if(tracks.size() != isInPrimary.size()){
        std::cout<<"ERROR: isInPrimary vector size not the same as no of tracks"<<std::endl;
    }
    if(debug_me) std::cout<<"tracks extracted from the reco particles"<<std::endl;

    // loop over jets
    for(unsigned int jetIdx = 0; jetIdx < jets.size(); jetIdx++){

        // get all secondary tracks for this jet
        ROOT::VecOps::RVec<edm4hep::TrackState> tracksInThisJet;
        std::vector<int> this_jet_constituent_indices = jet_constituent_indices[jetIdx];
        for (int ctr=0; ctr<tracks.size(); ctr++) {
            // skip primary tracks
            if(isInPrimary[ctr]) continue;
            // skip tracks that are not belonging to this jet
            // note: this seems to be a bug? should be != instead of == ?
            if(std::find(this_jet_constituent_indices.begin(), this_jet_constituent_indices.end(), reco_ind_tracks[ctr])
                == this_jet_constituent_indices.end()){
	            tracksInThisJet.push_back(tracks[ctr]);
            }
        }
    
        if(debug_me){
            std::cout << "primary tracks removed; there are " << tracksInThisJet.size() << " non-primary tracks";
            std::cout << " in jet n. " << jetIdx+1 << std::endl;
        }
    
        // V0 rejection (tight) - perform V0 rejection with tight constraints if user chooses
        ROOT::VecOps::RVec<edm4hep::TrackState> tracksToUse = V0rejection_tight(tracksInThisJet, PV, V0_rej);
    
        if(debug_me) {
            std::cout<<tracksInThisJet.size()-tracksToUse.size()<<" V0 tracks removed"<<std::endl;
            std::cout<<"now starting to find secondary vertices..."<<std::endl;
        }

        // add to selected tracks
        selectedTracksPerJet.push_back(tracksToUse);
    }
    
    // find secondary vertices using the selected tracks
    return get_SV_jets(selectedTracksPerJet, allTracks, PV, V0_rej, chi2_cut, invM_cut, chi2Tr_cut); 
}


/*
Secondary vertex finding per event
*/


ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
get_SV_event(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& allTracks,
        const VertexingUtils::FCCAnalysesVertex& PV,
        bool V0_veto,
        double chi2_cut, double invariant_mass_cut, double trackChi2_cut){
    // Base function for doing per-event secondary vertex finding.
    // Input arguments:
    // - tracks: tracks to use.
    //   No further filtering is performed on this collection,
    //   the tracks are assumed to be already selected
    //   (e.g. sufficient quality, incompatible with primary vertex, etc).
    //   Exception: a V0 veto can still be performed if V0_veto is set to true.
    // - allTracks: all tracks in event; not sure why this is needed.
    // - PV: primary vertex, not sure why this is needed.
    // Returns:
    // - a collection of secondary vertices.

    // do V0 veto if requested
    ROOT::VecOps::RVec<edm4hep::TrackState> tracksToUse = V0rejection_tight(tracks, PV, V0_veto);

    if(debug_me) {
        std::cout<<tracks.size()-tracksToUse.size()<<" V0 tracks removed"<<std::endl;
        std::cout<<"now starting to find secondary vertices..."<<std::endl;
    }

    // start finding secondary vertices
    return findSVfromTracks(tracksToUse, allTracks, PV, chi2_cut, invariant_mass_cut, trackChi2_cut);
}


ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
get_SV_event(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recoparticles,
		const ROOT::VecOps::RVec<edm4hep::TrackState>& allTracks,
		VertexingUtils::FCCAnalysesVertex PV,
		ROOT::VecOps::RVec<bool> isInPrimary,
		bool V0_rej,
		double chi2_cut, double invM_cut, double chi2Tr_cut){
    // Same as above (base function) but with more built-in input handling.
    
    // retrieve the tracks associated to the recoparticles
    ROOT::VecOps::RVec<edm4hep::TrackState> tracks = ReconstructedParticle2Track::getRP2TRK( recoparticles, allTracks );
    if(debug_me) std::cout<<"tracks extracted from the reco particles"<<std::endl;

    // sanity checks
    if(tracks.size() != isInPrimary.size()){
        std::cout<<"ISSUE: track vector and primary-nonprimary vector of diff sizes"<<std::endl;
    }

    // remove primary tracks
    ROOT::VecOps::RVec<edm4hep::TrackState> secondaryTracks;
    for(unsigned int i=0; i<tracks.size(); i++) {
        if (!isInPrimary[i]) secondaryTracks.push_back(tracks[i]);
    }

    if(debug_me) std::cout<<"primary tracks removed; there are "<<secondaryTracks.size()<<" non-primary tracks in the event"<<std::endl;

    // V0 rejection (tight) - perform V0 rejection with tight constraints if user chooses
    ROOT::VecOps::RVec<edm4hep::TrackState> tracksToUse = V0rejection_tight(secondaryTracks, PV, V0_rej);

    if(debug_me) {
        std::cout<<secondaryTracks.size()-tracksToUse.size()<<" V0 tracks removed"<<std::endl;
        std::cout<<"now starting to find secondary vertices..."<<std::endl;
    }
  
    // find secondary vertices
    return get_SV_event(tracksToUse, allTracks, PV, V0_rej, chi2_cut, invM_cut, chi2Tr_cut);
}


//** internal functions for SV finder **//

ROOT::VecOps::RVec<int> VertexSeed_best(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
					VertexingUtils::FCCAnalysesVertex PV,
					double chi2_cut, double invM_cut) {

  // gives indices of the best pair of tracks

  ROOT::VecOps::RVec<int> result;
  int isel = 0;
  int jsel = 1;
  
  int nTr = tracks.size();
  // push empty tracks to make a size=2 vector
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  edm4hep::TrackState tr_i, tr_j;
  tr_pair.push_back(tr_i);
  tr_pair.push_back(tr_j);
  VertexingUtils::FCCAnalysesVertex vtx_seed;
  double chi2_min = 99;
 
  // loop over pairs of tracks 
  for(unsigned int i=0; i<nTr-1; i++) {
    tr_pair[0] = tracks[i];
    TVector3 p_i( TMath::Cos(tracks[i].phi), TMath::Sin(tracks[i].phi), tracks[i].tanLambda );
    for(unsigned int j=i+1; j<nTr; j++) {
      tr_pair[1] = tracks[j];
      TVector3 p_j( TMath::Cos(tracks[j].phi), TMath::Sin(tracks[j].phi), tracks[j].tanLambda );

      // pair preselection
      // - opposite charge
      if(tr_pair[0].omega * tr_pair[1].omega > 0) continue;
      // - more or less same direction
      if( p_i.DeltaR(p_j) > 0.8 ) continue;
 
      // V0 rejection (tight)
      ROOT::VecOps::RVec<bool> isInV0 = isV0(tr_pair, PV, true);
      if(isInV0[0] && isInV0[1]) continue;
      
      // fit a common vertex (fast method) and do basic checks
      // (done only for speedup, but disabled since it doesn't gain that much and might have some unwanted effects)
      vtx_seed = VertexFitterSimple::VertexFitter_Tk(2, tr_pair, false, 0, 0, 0, 0, 0, 0, true);
      double chi2 = vtx_seed.vertex.chi2; // normalised, but ndof is 1 anyway
      if(chi2 >= chi2_cut) continue;

      // fit a common vertex (full method)
      vtx_seed = VertexFitterSimple::VertexFitter_Tk(2, tr_pair);

      // constraints check
      bool pass = check_constraints(vtx_seed, tr_pair, PV, true, chi2_cut, invM_cut);
      if(!pass) continue;
      
      // if a pair passes all constraints compare chi2, store lowest chi2
      double chi2_seed = vtx_seed.vertex.chi2; // normalised but nDOF=1 for nTr=2      
      if(chi2_seed < chi2_min) {
	    isel = i; jsel =j;
	    chi2_min = chi2_seed;
      }
    }
  }

  if(chi2_min != 99){
    result.push_back(isel); 
    result.push_back(jsel);
  }
  return result;
}

ROOT::VecOps::RVec<int> addTrack_best(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
				      ROOT::VecOps::RVec<int> vtx_tr,
				      VertexingUtils::FCCAnalysesVertex PV,
				      double chi2_cut, double invM_cut, double chi2Tr_cut) {
  // adds index of the best track to the (seed) vtx
  
  ROOT::VecOps::RVec<int> result = vtx_tr;
  if(tracks.size() == vtx_tr.size()) return result;
  
  // initializations
  int isel = -1;
  int nTr = tracks.size();
  VertexingUtils::FCCAnalysesVertex vtx;
  double chi2_min = 99;

  // add tracks of the previously formed vtx to a vector
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx;
  TVector3 p_tracks_sum;
  for(int tr : vtx_tr){
    tr_vtx.push_back(tracks[tr]);
    p_tracks_sum += TVector3( TMath::Cos(tracks[tr].phi), TMath::Sin(tracks[tr].phi), tracks[tr].tanLambda );
  }
  int iTr = tr_vtx.size();
  // add an empty track to increase vector size by 1
  edm4hep::TrackState tr_i;
  tr_vtx.push_back(tr_i);

  // loop over tracks
  for(unsigned int i=0; i<nTr; i++) {
    TVector3 p_i( TMath::Cos(tracks[i].phi), TMath::Sin(tracks[i].phi), tracks[i].tanLambda );

    // track preselection
    // - do not consider tracks already in the vertex (?)
    if(std::find(vtx_tr.begin(), vtx_tr.end(), i) != vtx_tr.end()) continue;
    // - more or less same direction as other tracks in vertex
    if( p_i.DeltaR(p_tracks_sum) > 0.8 ) continue;

    tr_vtx[iTr] = tracks[i];
    double nDOF = 2*tr_vtx.size() - 3; // nDOF 

    // fit a common vertex (fast method) and do basic checks
    // (done only for speedup, but disabled since it doesn't gain that much and might have some unwanted effects)
    vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_vtx, false, 0, 0, 0, 0, 0, 0, true);
    double chi2 = vtx.vertex.chi2; // normalised
    chi2 = chi2 * nDOF;
    if(chi2 >= chi2_cut) continue;
    ROOT::VecOps::RVec<float> chi2_tr = vtx.reco_chi2;
    if(chi2_tr[tr_vtx.size()-1] >= chi2Tr_cut) continue;

    // fit a common vertex (full method)
    vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_vtx);

    // check constraints
    bool pass = check_constraints(vtx, tr_vtx, PV, false, chi2_cut, invM_cut, chi2Tr_cut);
    if(!pass) continue;
    
    // if a track passes all constraints compare chi2, store lowest chi2
    double chi2_vtx = vtx.vertex.chi2; // normalised
    chi2_vtx = chi2_vtx * nDOF;
    if(chi2_vtx < chi2_min) {
      isel = i;
      chi2_min = chi2_vtx;
    }    
  }
  if(isel>=0) result.push_back(isel);
  return result;
}

ROOT::VecOps::RVec<edm4hep::TrackState> V0rejection_tight(
    ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
	VertexingUtils::FCCAnalysesVertex PV,
	bool V0_rej) {
  // perform V0 rejection with tight constraints if user chooses
  ROOT::VecOps::RVec<edm4hep::TrackState> result;
  if(V0_rej) {
    bool tight = true;
    ROOT::VecOps::RVec<bool> isInV0 = isV0(tracks, PV, tight);
    for(unsigned int i=0; i<isInV0.size(); i++) {
      if (!isInV0[i]) result.push_back(tracks[i]);
    }
  }
  else result = tracks;
  return result;
}

ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> findSVfromTracks(
    const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks_fin,
    const ROOT::VecOps::RVec<edm4hep::TrackState>&  alltracks,
    VertexingUtils::FCCAnalysesVertex PV,
    double chi2_cut, double invM_cut, double chi2Tr_cut) {
  // find SVs (only if there are 2 or more tracks)
  
  // initialize result
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;

  // make a modifiable copy of input tracks
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks_leftover = tracks_fin;

  // keep iterating as long as there are at least 2 leftover tracks
  while(tracks_leftover.size() > 1){
    
    // find vertex seed and break if none can be found
    ROOT::VecOps::RVec<int> vtx_seed_indices = VertexSeed_best(tracks_leftover, PV, chi2_cut, invM_cut);
    if(vtx_seed_indices.size() == 0) break;
    
    // add tracks to the seed
    ROOT::VecOps::RVec<int> vtx_indices = vtx_seed_indices;
    int vtx_indices_size = 0; // dummy to start the loop
    while(vtx_indices_size != vtx_indices.size()){
      vtx_indices_size = vtx_indices.size();
      vtx_indices = addTrack_best(tracks_leftover, vtx_indices, PV, chi2_cut, invM_cut, chi2Tr_cut);
    }
    
    // fit secondary vertex to tracks
    ROOT::VecOps::RVec<edm4hep::TrackState> tracks_in_vertex;
    for(int trackidx : vtx_indices){
      tracks_in_vertex.push_back(tracks_leftover[trackidx]);
    }
    VertexingUtils::FCCAnalysesVertex secondary_vertex;
    secondary_vertex = VertexFitterSimple::VertexFitter_Tk(2, tracks_in_vertex); // flag 2 for SVs
    result.push_back(secondary_vertex);

    // remove used tracks from collection of leftover tracks
    ROOT::VecOps::RVec<edm4hep::TrackState> temp = tracks_leftover;
    tracks_leftover.clear();
    for(unsigned int trackidx = 0; trackidx < temp.size(); trackidx++) {
      if(std::find(vtx_indices.begin(), vtx_indices.end(), trackidx) == vtx_indices.end()){
        tracks_leftover.push_back(temp[trackidx]);
      }
    }
  }
  return result;
}

bool check_constraints(VertexingUtils::FCCAnalysesVertex vtx,
		       ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
		       VertexingUtils::FCCAnalysesVertex PV,
		       bool seed,
		       double chi2_cut, double invM_cut, double chi2Tr_cut) {
  // if all constraints pass -> true
  // if any constraint fails -> false

  bool result = true;

  int nTr = tracks.size();  // no of tracks
  
  // Constraints
  // chi2 < cut (9)
  double chi2 = vtx.vertex.chi2; // normalised
  double nDOF = 2*nTr - 3; // nDOF
  chi2 = chi2 * nDOF;
  if(chi2 >= chi2_cut) result = false;

  // invM < cut (10GeV)
  double invM = VertexingUtils::get_invM(vtx);
  if(invM >= invM_cut) result = false;

  // invM < sum of energy
  double E_tracks = 0.;
  for(edm4hep::TrackState tr_e : tracks) E_tracks += VertexingUtils::get_trackE(tr_e);
  if(invM >= E_tracks) result = false;

  // momenta sum & vtx r on same side
  double angle = VertexingUtils::get_PV2vtx_angle(tracks, vtx, PV);
  if(angle<0) result = false;

  if(!seed) {
    // chi2_contribution(track) < threshold
    ROOT::VecOps::RVec<float> chi2_tr = vtx.reco_chi2;
    if(chi2_tr[nTr-1] >= chi2Tr_cut) result = false;
  }
  return result;
}

ROOT::VecOps::RVec<bool> isV0(
    ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
    VertexingUtils::FCCAnalysesVertex PV,
	bool tight){
  // find which tracks belong to V0 candidates (e.g. for V0 rejection)
  //
  // take all non-primary tracks & assign "true" to pairs that form V0
  // if(tight)  -> tight constraints
  // if(!tight) -> loose constraints

  int nTr = np_tracks.size();

  ROOT::VecOps::RVec<bool> result(nTr, false);
  // true -> forms a V0, false -> doesn't form a V0
  if(nTr<2) return result;

  // set constraints (if(tight==true) tight_set)
  ROOT::VecOps::RVec<double> isKs = constraints_Ks(tight);
  ROOT::VecOps::RVec<double> isLambda0 = constraints_Lambda0(tight);
  ROOT::VecOps::RVec<double> isGamma = constraints_Gamma(tight);
  
  // make a dummy track pair
  ROOT::VecOps::RVec<edm4hep::TrackState> t_pair;
  edm4hep::TrackState tr_i, tr_j;
  t_pair.push_back(tr_i);
  t_pair.push_back(tr_j);
  VertexingUtils::FCCAnalysesVertex V0;

  // loop over pairs of tracks
  for(unsigned int i=0; i<nTr-1; i++) {
    // do not consider tracks that are already paired in a V0 candidate.
    // update: disable this skip to be maximally inclusive in V0 reconstruction.
    //if(result[i] == true) continue;
    t_pair[0] = np_tracks[i];
    TVector3 p_i( TMath::Cos(np_tracks[i].phi), TMath::Sin(np_tracks[i].phi), np_tracks[i].tanLambda );
    for(unsigned int j=i+1; j<nTr; j++) {
      // do not consider tracks that are already paired in a V0 candidate.
      // update: disable this skip to be maximally inclusive in V0 reconstruction.
      //if(result[j] == true) continue;
      t_pair[1] = np_tracks[j];
      TVector3 p_j( TMath::Cos(np_tracks[j].phi), TMath::Sin(np_tracks[j].phi), np_tracks[j].tanLambda );

      // pair preselection
      // - opposite charge
      if(t_pair[0].omega * t_pair[1].omega > 0) continue;
      // - more or less same direction
      //if( p_i.DeltaR(p_j) > 0.4 ) continue;

      // fit common vertex
      double chi2max = 10.;
      ROOT::VecOps::RVec<double> V0_cand = get_V0candidate(V0, t_pair, PV, true, chi2max);
      if( V0_cand[0] < 0 ){ continue; }

      // Ks
      if(V0_cand[0]>isKs[0] && V0_cand[0]<isKs[1] && V0_cand[4]>isKs[2] && V0_cand[5]>isKs[3]) {
	    result[i] = true;
	    result[j] = true;
      }

      // Lambda0
      if(V0_cand[1]>isLambda0[0] && V0_cand[1]<isLambda0[1] && V0_cand[4]>isLambda0[2] && V0_cand[5]>isLambda0[3]) {
	    result[i] = true;
	    result[j] = true;
      }
      if(V0_cand[2]>isLambda0[0] && V0_cand[2]<isLambda0[1] && V0_cand[4]>isLambda0[2] && V0_cand[5]>isLambda0[3]) {
	    result[i] = true;
	    result[j] = true;
      }

      // photon conversion
      if(V0_cand[3]<isGamma[1] && V0_cand[4]>isGamma[2] && V0_cand[5]>isGamma[3]) {
	    result[i] = true;
	    result[j] = true;
      }	  
    }
  }

  return result;
}


///////////////////////////
//** V0 Reconstruction **//
///////////////////////////

VertexingUtils::FCCAnalysesV0 get_V0s(
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
	    VertexingUtils::FCCAnalysesVertex PV,
	    const ROOT::VecOps::RVec<double>& constraints_ks,
        const ROOT::VecOps::RVec<double>& constraints_lambda0,
        const ROOT::VecOps::RVec<double>& constraints_gamma,
		double chi2_cut) {
  // Base function for V0 reconstruction in a set of provided tracks.

  // initializations
  VertexingUtils::FCCAnalysesV0 result;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vtx; // FCCAnalyses vertex object
  ROOT::VecOps::RVec<int> pdgAbs;                            // absolute PDG ID
  ROOT::VecOps::RVec<double> invM;                           // invariant mass
  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;
  VertexingUtils::FCCAnalysesVertex V0_vtx;
  
  // check number of tracks
  int nTr = tracks.size();
  if(nTr<2) return result;
  ROOT::VecOps::RVec<bool> isInV0(nTr, false);

  // initialize (empty) track pair
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  edm4hep::TrackState tr_i, tr_j;
  tr_pair.push_back(tr_i);
  tr_pair.push_back(tr_j);
  
  // loop over pairs of tracks
  for(unsigned int i=0; i<nTr-1; i++) {
    // do not consider tracks that are already paired in a V0 candidate.
    // update: disable this skip to be maximally inclusive in V0 reconstruction.
    //if(isInV0[i] == true) continue;
    tr_pair[0] = tracks[i];
    TVector3 p_i( TMath::Cos(tracks[i].phi), TMath::Sin(tracks[i].phi), tracks[i].tanLambda );
    for(unsigned int j=i+1; j<nTr; j++) {
      // do not consider tracks that are already paired in a V0 candidate.
      // update: disable this skip to be maximally inclusive in V0 reconstruction.
      //if(isInV0[j] == true) continue;
      tr_pair[1] = tracks[j];
      TVector3 p_j( TMath::Cos(tracks[j].phi), TMath::Sin(tracks[j].phi), tracks[j].tanLambda );

      // pair preselection
      // - opposite charge
      if(tr_pair[0].omega * tr_pair[1].omega > 0) continue;
      // - more or less same direction
      //if( p_i.DeltaR(p_j) > 0.4 ) continue;

      // try to make a V0 candidate
      ROOT::VecOps::RVec<double> V0_cand = get_V0candidate(V0_vtx, tr_pair, PV, true, chi2_cut);
      if(V0_cand[0] == -1) continue;

      // decide which one it is
      // update: now allow a candidate to be identified as more than one type,
      // e.g. a kaon under the pi-pi hypothesis but also a lambda under the pi-p hypothesis.
      
      // Ks
      if(V0_cand[0]>constraints_ks[0] && V0_cand[0]<constraints_ks[1] && V0_cand[4]>constraints_ks[2] && V0_cand[5]>constraints_ks[3]) {
	    if(debug_me) std::cout<<"Found a Ks"<<std::endl;
        isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(310);
	    invM.push_back(V0_cand[0]);
      }
      
      // Lambda0
      if(V0_cand[1]>constraints_lambda0[0] && V0_cand[1]<constraints_lambda0[1] && V0_cand[4]>constraints_lambda0[2] && V0_cand[5]>constraints_lambda0[3]) {
	    if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(3122);
	    invM.push_back(V0_cand[1]);
      }
      if(V0_cand[2]>constraints_lambda0[0] && V0_cand[2]<constraints_lambda0[1] && V0_cand[4]>constraints_lambda0[2] && V0_cand[5]>constraints_lambda0[3]) {
	    if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(3122);
	    invM.push_back(V0_cand[2]);
      }
	
      // photon conversion
      if(V0_cand[3]<constraints_gamma[1] && V0_cand[4]>constraints_gamma[2] && V0_cand[5]>constraints_gamma[3]) {
	    if(debug_me) std::cout<<"Found a Photon coversion"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(22);
	    invM.push_back(V0_cand[3]);
      }
    }
  } // end of loop over track pairs

  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;
  
  return result;
}

VertexingUtils::FCCAnalysesV0 get_V0s(
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
        VertexingUtils::FCCAnalysesVertex PV,
        bool tight,
        double chi2_cut) {
    // same as above, but with predefined constraints
    ROOT::VecOps::RVec<double> isKs = constraints_Ks(tight);
    ROOT::VecOps::RVec<double> isLambda0 = constraints_Lambda0(tight);
    ROOT::VecOps::RVec<double> isGamma = constraints_Gamma(tight);
    return get_V0s(tracks, PV,
             isKs,
             isLambda0,
             isGamma,
             chi2_cut);
}

VertexingUtils::FCCAnalysesV0 get_V0s(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
				      VertexingUtils::FCCAnalysesVertex PV,
				      double Ks_invM_low, double Ks_invM_high, double Ks_dis, double Ks_cosAng,
				      double Lambda_invM_low, double Lambda_invM_high, double Lambda_dis, double Lambda_cosAng,
				      double Gamma_invM_low, double Gamma_invM_high, double Gamma_dis, double Gamma_cosAng,
				      double chi2_cut) {
  // same as above, but with constraints given in different format
  ROOT::VecOps::RVec<double> isKs      = constraints_Ks(Ks_invM_low, Ks_invM_high, Ks_dis, Ks_cosAng);
  ROOT::VecOps::RVec<double> isLambda0 = constraints_Lambda0(Lambda_invM_low, Lambda_invM_high, Lambda_dis, Lambda_cosAng);
  ROOT::VecOps::RVec<double> isGamma   = constraints_Gamma(Gamma_invM_low, Gamma_invM_high, Gamma_dis, Gamma_cosAng);
  return get_V0s(tracks, PV,
             isKs,
             isLambda0,
             isGamma,
             chi2_cut);
}

VertexingUtils::FCCAnalysesV0 get_V0s_jet(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
					  ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
					  ROOT::VecOps::RVec<bool> isInPrimary,
					  ROOT::VecOps::RVec<fastjet::PseudoJet> jets,
					  std::vector<std::vector<int>> jet_consti,
					  VertexingUtils::FCCAnalysesVertex PV,
					  bool tight,
					  double chi2_cut) {
  // V0 reconstruction after jet clustering
  // if(tight)  -> tight constraints
  // if(!tight) -> loose constraints

  VertexingUtils::FCCAnalysesV0 result;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vtx; // FCCAnalyses vertex object
  ROOT::VecOps::RVec<int> pdgAbs;                            // absolute PDG ID
  ROOT::VecOps::RVec<double> invM;                           // invariant mass
  ROOT::VecOps::RVec<int> nSV_jet(jets.size(),0);
  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;
  result.nSV_jet = nSV_jet;

  VertexingUtils::FCCAnalysesVertex V0_vtx;

  int n_par = recoparticles.size();
  if(n_par<2) return result;

  // set constraints (if(tight==true) tight_set)
  ROOT::VecOps::RVec<double> isKs      = constraints_Ks(tight);
  ROOT::VecOps::RVec<double> isLambda0 = constraints_Lambda0(tight);
  ROOT::VecOps::RVec<double> isGamma   = constraints_Gamma(tight);
  
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  // push empty tracks to make a size=2 vector
  edm4hep::TrackState tr_i, tr_j;
  tr_pair.push_back(tr_i);
  tr_pair.push_back(tr_j);

  // find V0s inside the jet loop (only from non-primary tracks)
  ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks;
  
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks   = ReconstructedParticle2Track::getRP2TRK( recoparticles, thetracks );
  ROOT::VecOps::RVec<int> reco_ind_tracks     = ReconstructedParticle2Track::get_recoindTRK( recoparticles, thetracks );
  if(tracks.size() != reco_ind_tracks.size()) std::cout<<"ERROR: reco index vector not the same size as no of tracks"<<std::endl;

  if(tracks.size() != isInPrimary.size()) std::cout<<"ERROR: isInPrimary vector size not the same as no of tracks"<<std::endl;
  
  if(debug_me) std::cout<<"tracks extracted from the reco particles"<<std::endl;
  
  //
  for (unsigned int j=0; j<jets.size(); j++) {

    int i_nSV = 0;
    
    // remove primary tracks & separate non-primary tracks by jet
    std::vector<int> i_jetconsti = jet_consti[j];
    for (int ctr=0; ctr<tracks.size(); ctr++) {
      if(isInPrimary[ctr]) continue; // remove primary tracks
      if(std::find(i_jetconsti.begin(), i_jetconsti.end(), reco_ind_tracks[ctr]) == i_jetconsti.end()) {
	np_tracks.push_back(tracks[ctr]); // separate tracks by jet
      }
    }
    
    if(debug_me) std::cout<<"primary tracks removed; there are "<<np_tracks.size()<<" non-primary tracks in jet#"<<j+1<<std::endl;

    int nTr = np_tracks.size();
    if(nTr<2) continue;    
    ROOT::VecOps::RVec<bool> isInV0(nTr, false);
    //
    for(unsigned int i=0; i<nTr-1; i++) {
      if(isInV0[i] == true) continue; // don't pair a track if it already forms a V0
      tr_pair[0] = np_tracks[i];
      
      for(unsigned int j=i+1; j<nTr; j++) {
	if(isInV0[j] == true) continue; // don't pair a track if it already forms a V0
	if(tr_pair[0].omega * np_tracks[j].omega > 0) continue; // don't pair tracks with same charge (same sign curvature = same sign charge)
	tr_pair[1] = np_tracks[j];
	
	ROOT::VecOps::RVec<double> V0_cand = get_V0candidate(V0_vtx, tr_pair, PV, true, chi2_cut);
	if(V0_cand[0] == -1) continue;
	
	// Ks
	if(V0_cand[0]>isKs[0] && V0_cand[0]<isKs[1] && V0_cand[4]>isKs[2] && V0_cand[5]>isKs[3]) {
	  if(debug_me) std::cout<<"Found a Ks"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(310);
	  invM.push_back(V0_cand[0]);
	  i_nSV++;
	  break;
	}
	  
	// Lambda0
	else if(V0_cand[1]>isLambda0[0] && V0_cand[1]<isLambda0[1] && V0_cand[4]>isLambda0[2] && V0_cand[5]>isLambda0[3]) {
	  if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(V0_cand[1]);
	  i_nSV++;
	  break;
	}
	else if(V0_cand[2]>isLambda0[0] && V0_cand[2]<isLambda0[1] && V0_cand[4]>isLambda0[2] && V0_cand[5]>isLambda0[3]) {
	  if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(V0_cand[2]);
	  i_nSV++;
	  break;
	}
	
	// photon conversion
	else if(V0_cand[3]<isGamma[1] && V0_cand[4]>isGamma[2] && V0_cand[5]>isGamma[3]) {
	  if(debug_me) std::cout<<"Found a Photon coversion"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(22);
	  invM.push_back(V0_cand[3]);
	  i_nSV++;
	  break;
	}
	//      
      }
    }

    nSV_jet[j] = i_nSV;
    // clean-up
    np_tracks.clear();
  } // jet loop ends

  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;
  result.nSV_jet = nSV_jet;
  //
  return result;
}

ROOT::VecOps::RVec<double> get_V0candidate(VertexingUtils::FCCAnalysesVertex &V0_vtx,
					   ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair,
					   VertexingUtils::FCCAnalysesVertex PV,
					   bool chi2,
					   double chi2_cut)
{
  // get invariant mass, distance from PV, and colliniarity variables for all V0 candidates
  
  // [0] -> invM_Ks [GeV]
  // [1] -> invM_Lambda1 [GeV]
  // [2] -> invM_Lambda2 [GeV]
  // [3] -> invM_Gamma [GeV]
  // [4] -> r (distance from PV) [mm]
  // [5] -> r.p (colinearity) [r & p - unit vectors]
  // skip the candidate with output entries = -1
  
  ROOT::VecOps::RVec<double> result(6, -1);
  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm
  
  // fit a common vertex (fast method) and do basic checks
  /*if( chi2 ){
    V0_vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_pair, false, 0, 0, 0, 0, 0, 0, true);
    double chi2 = V0_vtx.vertex.chi2; // normalised but ndof = 1
    if(chi2 >= chi2_cut) return result;
  }*/

  // fit a common vertex (full method)
  V0_vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_pair);

  if(chi2) {
    // constraint on chi2: chi2 < cut (9)
    double chi2_V0 = V0_vtx.vertex.chi2; // normalised but nDOF=1
    if(chi2_V0 >= chi2_cut) return result;
  }
  
  // invariant masses for V0 candidates
  result[0] = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_pi);
  result[1] = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_p);
  result[2] = VertexingUtils::get_invM_pairs(V0_vtx, m_p, m_pi);
  result[3] = VertexingUtils::get_invM_pairs(V0_vtx, m_e, m_e);

  // V0 candidate distance from PV
  edm4hep::Vector3f r_V0 = V0_vtx.vertex.position;
  TVector3 r_V0_PV(r_V0[0] - r_PV[0], r_V0[1] - r_PV[1], r_V0[2] - r_PV[2]);
  result[4] = r_V0_PV.Mag();

  // angle b/n V0 candidate momentum & PV-V0 displacement vector
  result[5] = VertexingUtils::get_PV2V0angle(V0_vtx, PV);

  return result;
}

// functions to fill constraint thresholds
// tight  -> tight constraints
// !tight -> loose constraints
//
// [0] -> invariant mass lower limit [GeV]
// [1] -> invariant mass upper limit [GeV]
// [2] -> distance from PV lower limit [mm] (note: [cm] for ALEPH data (?))
// [3] -> colinearity lower limit

ROOT::VecOps::RVec<double> constraints_Ks(bool tight) {

  ROOT::VecOps::RVec<double> result(4, 0);

  if(tight) {
    result[0] = 0.453; // widened based on observed Ks resolution in ALEPH; originally 0.493;
    result[1] = 0.553; // widened based on observed Ks resolution in ALEPH; orginally 0.503;
    result[2] = 0.1; // originally 0.5 [mm]
    result[3] = 0.999;
  }

  else {
    result[0] = 0.1; // widended to include sideband; originally 0.488;
    result[1] = 1.4; // widened to include sideband; originally 0.508;
    result[2] = 0.1; // originally 0.3 [mm]
    result[3] = 0.0; // widened based on observed values in ALEPH; originally 0.999;
  }
  
  return result;
}

ROOT::VecOps::RVec<double> constraints_Lambda0(bool tight) {

  ROOT::VecOps::RVec<double> result(4, 0);

  if(tight) {
    result[0] = 1.06; // widened based on observed Lambda resolution in ALEPH; originally 1.111;
    result[1] = 1.16; // widened based on observed Lambda resolution in ALEPH; originally 1.121;
    result[2] = 0.1; // originally 0.5 [mm]
    result[3] = 0.99995;
  }

  else {
    result[0] = 0.1; // widened to include sideband; originally 1.106;
    result[1] = 1.4; // widened to include sideband; originally 1.126;
    result[2] = 0.1; // originally 0.3 [mm]
    result[3] = 0.0; // widened based on observed values in ALEPH; originally 0.999;
  }

  return result;
}

ROOT::VecOps::RVec<double> constraints_Gamma(bool tight) {

  ROOT::VecOps::RVec<double> result(4, 0);

  if(tight) {
    result[1] = 0.005;
    result[2] = 0.9; // originally 9 [mm]
    result[3] = 0.99995;
  }

  else {
    result[1] = -1; // originally 0.01;
    result[2] = 0.9; // originally 9 [mm]
    result[3] = 0.999;
  }

  return result;
}

// User set constraints
ROOT::VecOps::RVec<double> constraints_Ks(double invM_low, double invM_high, double dis, double cosAng) {

  ROOT::VecOps::RVec<double> result(4, 0);

  result[0] = invM_low;
  result[1] = invM_high;
  result[2] = dis;
  result[3] = cosAng;
  //
  return result;
}

ROOT::VecOps::RVec<double> constraints_Lambda0(double invM_low, double invM_high, double dis, double cosAng) {

  ROOT::VecOps::RVec<double> result(4, 0);

  result[0] = invM_low;
  result[1] = invM_high;
  result[2] = dis;
  result[3] = cosAng;
  //
  return result;
}

ROOT::VecOps::RVec<double> constraints_Gamma(double invM_low, double invM_high, double dis, double cosAng) {

  ROOT::VecOps::RVec<double> result(4, 0);

  result[0] = invM_low;
  result[1] = invM_high;
  result[2] = dis;
  result[3] = cosAng;
  //
  return result;
}

// ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> VertexSeed_all(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
// 					  		      VertexingUtils::FCCAnalysesVertex PV,
// 							      double chi2_cut, double invM_cut) {

  // // gives indices of the all pairs of tracks which pass the constraints

  // ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
  // ROOT::VecOps::RVec<int> ij_sel;
  
  // int nTr = tracks.size();
  // ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  // // push empty tracks to make a size=2 vector
  // edm4hep::TrackState tr_i;
  // edm4hep::TrackState tr_j;
  // tr_pair.push_back(tr_i);
  // tr_pair.push_back(tr_j);
  // VertexingUtils::FCCAnalysesVertex vtx_seed;
  
//   for(unsigned int i=0; i<nTr-1; i++) {
//     if(i!=0) tr_pair[0] = tracks[i];

//     for(unsigned int j=i+1; j<nTr; j++) {
//       if(j!=1) tr_pair[1] = tracks[j];

//       // V0 rejection (loose)
//       ROOT::VecOps::RVec<bool> isInV0 = isV0(tr_pair, PV, false);
//       if(isInV0[0] && isInV0[1]) continue;
      
//       vtx_seed = VertexFitterSimple::VertexFitter_Tk(2, tr_pair);

//       // Constraints
//       bool pass = check_constraints(vtx_seed, tr_pair, PV, true, chi2_cut, invM_cut, chi2Tr_cut);
//       if(!pass) continue;

//       // if a pair passes all constraints, store indices
//       ij_sel.push_back(i); ij_sel.push_back(j);
//       result.push_back(ij_sel);
//       ij_sel.clear();
//     }
//   }

//   return result;
// }


// ROOT::VecOps::RVec<int> addTrack_multi(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
//				          ROOT::VecOps::RVec<int> vtx_tr,
// 					  VertexingUtils::FCCAnalysesVertex PV,
// 					  double chi2_cut, double invM_cut, double chi2Tr_cut) {
//   // adds indices of all tracks passing constraints to the (seed) vtx
  
//   ROOT::VecOps::RVec<int> result = vtx_tr;
//   if(tracks.size() == vtx_tr.size()) return result;

//   int nTr = tracks.size();
//   ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx;
//   VertexingUtils::FCCAnalysesVertex vtx;

//   // tracks from the previously formed vtx
//   for(int tr : vtx_tr) {
//     tr_vtx.push_back(tracks[tr]);
//   }
//   int iTr = tr_vtx.size();
  
//   // find best track to add to the vtx
//   for(unsigned int i=0; i<nTr; i++) {
//     if(std::find(vtx_tr.begin(), vtx_tr.end(), i) != vtx_tr.end()) continue;

//     if(iTr != tr_vtx.size()) tr_vtx[iTr] = tracks[i];
//     else tr_vtx.push_back(tracks[i]);
    
//     vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_vtx);

//     // Constraints
//     bool pass = check_constraints(vtx_seed, tr_pair, PV, false, chi2_cut, invM_cut, chi2Tr_cut);
//     if(!pass) continue;
    
//     // if a pair passes all constraints add to the vtx
//     result.push_back(i);
//     iTr++;
//   }

//   return result;
// }

}//end NS VertexFinderLCFIPlus

}//end NS FCCAnalyses
