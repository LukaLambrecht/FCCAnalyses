#include "FCCAnalyses/ReconstructedParticle2Track.h"
#include "FCCAnalyses/VertexingUtils.h"

namespace FCCAnalyses{

namespace ReconstructedParticle2Track{

  /*
  Define a switch between different detectors (fcc or aleph);
  not the cleanest solution to hard-code it here, but good enough for now...
  */
  const std::string detector = "aleph"; // choose from "fcc" or "aleph"


  /*
  Indexing methods
  */ 

  size_t getTrackIndex(const edm4hep::ReconstructedParticleData& rp,
                       const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    /*
    Get index of track in track collection corrsponding to a reco particle.
    Note: depending on the edm4hep version, this is NOT simply RecoParticle.tracks_begin...
    */
    if(rp.tracks_begin < reco2track_links.size()) {
      const auto &oid = reco2track_links.at(rp.tracks_begin);
      size_t trackIndex = oid.index;
      return trackIndex;
    }
    return 9999;
  }


  /*
  General
  */

  ROOT::VecOps::RVec<edm4hep::TrackState>
  getRP2TRK_trackState(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
    const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    /*
    Get the full TrackState object for each input particle.
    Note: contrary to most other functions of this template,
          no dummy value is returned for neutral particles;
          hence the output vector contains all the tracks of the set of input particles,
          but the length might be different and there is no one-to-one correspondence!
          this function is useful e.g. for getting all tracks of a given jet
          (or other collection of particles) but not for extracting per-particle properties.
    */

    // initializations
    ROOT::VecOps::RVec<edm4hep::TrackState> result;

    // loop over reco particles
    for (auto & p: rps) {
        // find track
        size_t trackIndex = getTrackIndex(p, reco2track_links);
        if(trackIndex < tracks.size()){
            edm4hep::TrackState tr = tracks.at(trackIndex);
            result.push_back(tr);
        }
    }
    return result;
  }

  ROOT::VecOps::RVec<edm4hep::TrackData>
  getRP2TRK_track(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
    const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
    const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    /*
    Get the full Track object for each input particle.
    Note: contrary to most other functions of this template,
          no dummy value is returned for neutral particles;
          hence the output vector contains all the tracks of the set of input particles,
          but the length might be different and there is no one-to-one correspondence!
          this function is useful e.g. for getting all tracks of a given jet
          (or other collection of particles) but not for extracting per-particle properties.
    */

    // initializations
    ROOT::VecOps::RVec<edm4hep::TrackData> result;

    // loop over reco particles
    for (auto & p: rps) {
        // find track
        size_t trackIndex = getTrackIndex(p, reco2track_links);
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            result.push_back(tr);
        }
    }
    return result;
  } 

  ROOT::VecOps::RVec<float> 
  getRP2TRK_mom(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
                ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<float> result;
    for (auto & p: in) {
      if (p.tracks_begin<tracks.size())
        result.push_back(VertexingUtils::get_trackMom(tracks.at(p.tracks_begin)));
      else result.push_back(std::nan(""));
    }
    return result;
  }

  ROOT::VecOps::RVec<float> 
  getRP2TRK_charge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
                   ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<float> result;
    for (auto & p: in) {
      if (p.tracks_begin<tracks.size())
        result.push_back(p.charge);
      else result.push_back(std::nan(""));
    }
    return result;
  }

  ROOT::VecOps::RVec<float>
  getRP2TRK_chi2(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
                   const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
                   const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    // get chi2 of track fit
    
    // initializations
    ROOT::VecOps::RVec<float> result;

    // loop over reco particles
    for (auto & p: in) {
        bool valid = false;
        // find track
        size_t trackIndex = getTrackIndex(p, reco2track_links);
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            float chi2 = tr.chi2;
            result.push_back(chi2);
            valid = true;
        }
        if(!valid){ result.push_back( -1 ); }
    }
    return result;
  }

  ROOT::VecOps::RVec<int>
  getRP2TRK_ndof(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
                   const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
                   const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    // get number of degrees of freedom of track fit

    // initializations
    ROOT::VecOps::RVec<int> result;

    // loop over reco particles
    for (auto & p: in) {
        bool valid = false;
        // find track
        size_t trackIndex = getTrackIndex(p, reco2track_links);
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            int ndof = tr.ndf;
            result.push_back(ndof);
            valid = true;
        }
        if(!valid){ result.push_back( -1 ); }
    }
    return result;
  }

  ROOT::VecOps::RVec<float>
  getRP2TRK_chi2Normalized(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
                   const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
                   const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    // get chi2 of track fit

    // initializations
    ROOT::VecOps::RVec<float> result;

    // loop over reco particles
    for (auto & p: in) {
        bool valid = false;
        // find track
        size_t trackIndex = getTrackIndex(p, reco2track_links);
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            float chi2 = tr.chi2;
            int ndof = tr.ndf;
            if(ndof > 0){
                result.push_back(chi2 / (float)ndof);
                valid = true;
            }
        }
        if(!valid){ result.push_back( -1 ); }
    }
    return result;
  }


  /*
  Get number of hits in subdetectors.
  Note: custom addition for Aleph files, not sure how to use for FCC.
  */
  
  ROOT::VecOps::RVec<int> getRP2TRK_nTrackHits(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<int>& subdetectorHitNumbers,
        const int subdetectorNumber,
        const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    // Get number of track hits for a given set of reco particles
    // for a given subdetector number

    // initializations
    ROOT::VecOps::RVec<float> out;

    // loop over reco particles
    for(auto & p: rps) {
        bool valid = false;
        // find track
        size_t trackIndex = getTrackIndex(p, reco2track_links);
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            // find index in hit collection
            size_t hitIdx = tr.subdetectorHitNumbers_begin + subdetectorNumber;
            if(hitIdx < subdetectorHitNumbers.size()){
                int nHits = subdetectorHitNumbers.at(hitIdx);
                out.push_back(nHits);
                valid = true;
            }
        }
        if(!valid){ out.push_back( -1 ); }
    }
    return out;
  }

  ROOT::VecOps::RVec<int> getRP2TRK_nTrackHits_VDET(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<int>& subdetectorHitNumbers,
        const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    // Get number of VDET (vertex detector) hits for a given set of reco particles.
    // Note: assumes inside-out numbering of the subdetectors,
    // such that VDET is at index 0 (to check).
    return getRP2TRK_nTrackHits(rps, tracks, subdetectorHitNumbers, 0, reco2track_links);
  }

  ROOT::VecOps::RVec<int> getRP2TRK_nTrackHits_ITC(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<int>& subdetectorHitNumbers,
        const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    // Get number of ITC (inner tracking chamber) hits for a given set of reco particles.
    // Note: assumes inside-out numbering of the subdetectors,
    // such that ITC is at index 1 (to check).
    return getRP2TRK_nTrackHits(rps, tracks, subdetectorHitNumbers, 1, reco2track_links);
  }

  ROOT::VecOps::RVec<int> getRP2TRK_nTrackHits_TPC(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<int>& subdetectorHitNumbers,
        const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    // Get number of TPC (time projection chamber) hits for a given set of reco particles.
    // Note: assumes inside-out numbering of the subdetectors,
    // such that TPC is at index 2 (to check).
    return getRP2TRK_nTrackHits(rps, tracks, subdetectorHitNumbers, 2, reco2track_links);
  }


  /*
  Magnetic field from track parameters
  */

  ROOT::VecOps::RVec<float> getRP2TRK_Bz(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
        const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links) {
    // Get the magnetic field strength from the momentum of a particle
    // and the curvature of its associated track.
    // Based on the formula p[GeV/c] = c[= 0.3 in these units] B[T] r[m]
    // with p the momentum, c the speed of light, B the magnetic field strength,
    // and r the radius of curvature.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units
    // (see formula above, and also accounting for r in mm instead of in m).
    ROOT::VecOps::RVec<float> out;

    // loop over reco particles
    for(auto & p: rps) {
        bool valid = false;
        size_t trackIndex = getTrackIndex(p, reco2track_links);
        if(trackIndex < trackStates.size()) {
            edm4hep::TrackState trst = trackStates.at(trackIndex); 
	        double pt = sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y);
            double omega = trst.omega;
            if( detector=="aleph" ){
                omega *= 0.1;
                // extra numerical factor needed for Aleph data
                // (because of the curvature being expressed in 1/cm instead of 1/mm)
            }
	        double Bz = omega / cSpeed * pt * std::copysign(1.0, p.charge);
	        out.push_back(Bz);
            valid = true;
        }
        if(!valid){ out.push_back(-9.); }
    }
    return out;
  }

  float Bz(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
           const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
           const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links) {
    // Get the magnetic field strength for an event
    // from the momenta of the particles and the curvatures of their associated tracks.
    // Based on the formula p[GeV/c] = c[= 0.3 in these units] B[T] r[m]
    // with p the momentum, c the speed of light, B the magnetic field strength,
    // and r the radius of curvature.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units
    // (see formula above, and also accounting for r in mm instead of in m).
    double Bz = -9; // dummy value if no valid result found

    // loop over reco particles
    for(auto & p: rps) {
      size_t trackIndex = getTrackIndex(p, reco2track_links);
      if(trackIndex >= trackStates.size()){ continue; }
      edm4hep::TrackState trst = trackStates.at(trackIndex);
      double pt = sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y);
      double omega = trst.omega;
      if( detector=="aleph" ){
        omega *= 0.1;
        // extra numerical factor needed for Aleph data
        // (because of the curvature being expressed in 1/cm instead of 1/mm)
      }
      Bz = omega / cSpeed * pt * std::copysign(1.0, p.charge);
      // break after the first particle with valid result
      break;
    }
    return Bz;
  }


  /*
  Impact parameters
  */

  ROOT::VecOps::RVec<float> XPtoPar_dxy(
      const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
	  const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
      const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links,
	  const TLorentzVector& PV,
	  const float& Bz) {
    // Get dxy with respect to primary vertex,
    // given the track state parameters calculated with respect to the origin.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units (see Bz calculation)
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {
      size_t trackIndex = getTrackIndex(rp, reco2track_links);
      if(trackIndex < trackStates.size()){
        edm4hep::TrackState trst = trackStates.at(trackIndex);

        float D0_wrt0 = trst.D0;
        float Z0_wrt0 = trst.Z0;
        float phi0_wrt0 = trst.phi;

        // note: phi0 is not the position vector azimuth,
        // but the azimuth of the momentum vector at the point of closest approach!
        TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
        TVector3 x = X - PV.Vect();
        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
        double pt = p.Pt();

        double a = - rp.charge * Bz * cSpeed;
        if( detector=="aleph" ){
            // extra numerical factor for aleph, see Bz calculation
            a *= -10;
        }
        double r2 = x(0) * x(0) + x(1) * x(1);
        double cross = x(0) * p(1) - x(1) * p(0);
        double D=-9;
        if (pt * pt - 2 * a * cross + a * a * r2 > 0) {
          double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
      	  if (pt < 10.0) D = (T - pt) / a;
          else D = (-2 * cross + a * r2) / (T + pt);
        }
	    out.push_back(D);
      } else { out.push_back(-9.); }
    }
    return out;
  }


  ROOT::VecOps::RVec<float> XPtoPar_dz(
      const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
      const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
      const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links,
      const TLorentzVector& PV,
      const float& Bz) {
    // Get dz with respect to primary vertex,
    // given the track state parameters calculated with respect to the origin.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {
      size_t trackIndex = getTrackIndex(rp, reco2track_links);
      if(trackIndex < trackStates.size()){
        edm4hep::TrackState trst = trackStates.at(trackIndex);

        float D0_wrt0 = trst.D0;
        float Z0_wrt0 = trst.Z0;
        float phi0_wrt0 = trst.phi;

        // note: phi0 is not the position vector azimuth,
        // but the azimuth of the momentum vector at the point of closest approach!
        TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
        TVector3 x = X - PV.Vect();
        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
        double pt = p.Pt();

        double a = - rp.charge * Bz * cSpeed;
        if( detector=="aleph" ){
            // extra numerical factor for aleph, see Bz calculation
            a *= -10;
        }
        double C = a/(2 * pt);
        double r2 = x(0) * x(0) + x(1) * x(1);
        double cross = x(0) * p(1) - x(1) * p(0);
        double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
        double D;
        if (pt < 10.0) D = (T - pt) / a;
        else D = (-2 * cross + a * r2) / (T + pt);
        double B = C * TMath::Sqrt(TMath::Max(r2 - D * D, 0.0) / (1 + 2 * C * D));
        if ( TMath::Abs(B) > 1.) B = TMath::Sign(1, B);
        double st = TMath::ASin(B) / C;
        double ct = p(2) / pt;
        double z0;
        double dot = x(0) * p(0) + x(1) * p(1);
        if (dot > 0.0) z0 = x(2) - ct * st;
        else z0 = x(2) + ct * st;

        out.push_back(z0);
      } else { out.push_back(-9.); }
    }
    return out;
  }

  ROOT::VecOps::RVec<float> XPtoPar_phi(
      const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
      const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
      const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links,
      const TLorentzVector& PV,
      const float& Bz) {
    // Get phi with respect to primary vertex,
    // given the track state parameters calculated with respect to the origin.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {
      size_t trackIndex = getTrackIndex(rp, reco2track_links);
      if(trackIndex < trackStates.size()){
        edm4hep::TrackState trst = trackStates.at(trackIndex);

        float D0_wrt0 = trst.D0;
        float Z0_wrt0 = trst.Z0;
        float phi0_wrt0 = trst.phi;

        // note: phi0 is not the position vector azimuth,
        // but the azimuth of the momentum vector at the point of closest approach!
        TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
        TVector3 x = X - PV.Vect();
        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
        double pt = p.Pt();

        double a = - rp.charge * Bz * cSpeed;
        if( detector=="aleph" ){
            // extra numerical factor for aleph, see Bz calculation
            a *= -10;
        }
        double r2 = x(0) * x(0) + x(1) * x(1);
        double cross = x(0) * p(1) - x(1) * p(0);
        double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
        double phi0 = TMath::ATan2((p(1) - a * x(0)) / T, (p(0) + a * x(1)) / T);

	    out.push_back(phi0);

      } else { out.push_back(-9.); }
    }
    return out;
  }

  ROOT::VecOps::RVec<float> XPtoPar_C(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
				       const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
				       const float& Bz) {

    const double cSpeed = 2.99792458e8 * 1.0e3 * 1.0e-15;
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {

      if( rp.tracks_begin < tracks.size()) {

        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);

        double a = std::copysign(1.0, rp.charge) * Bz * cSpeed;
	    double pt = p.Pt();
        double C = a/(2 * pt);

	    out.push_back(C);

      } else { out.push_back(-9.); }
    }
    return out;
  }

  ROOT::VecOps::RVec<float> XPtoPar_ct(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
				       const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
				       const float& Bz) {

    const double cSpeed = 2.99792458e8 * 1.0e-9;
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {

      if( rp.tracks_begin < tracks.size()) {

        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
	    double pt = p.Pt();

        double ct = p(2) / pt;

	    out.push_back(ct);

      } else { out.push_back(-9.); }
    }
    return out;
  }


/*
Simple track parameter getter functions
*/

ROOT::VecOps::RVec<float>
getRP2TRK_D0(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
             ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
             ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex < trackStates.size()){
      float D0_wrt0 = trackStates.at(trackIndex).D0;
      result.push_back(D0_wrt0);
    }
    else result.push_back(-9.);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_D0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
			     ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
                 ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex < trackStates.size()){
      result.push_back(trackStates.at(trackIndex).covMatrix[0]);
    }
    else result.push_back(-9.);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_D0_sig(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
		         ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
                 ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex < trackStates.size()){
      float D0_wrt0 = trackStates.at(trackIndex).D0;
      result.push_back(D0_wrt0/sqrt(trackStates.at(trackIndex).covMatrix[0]));
    }
    else result.push_back(-9.);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_Z0(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
		     ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
             ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).Z0);
    else result.push_back(-9.);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_Z0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[9]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_Z0_sig(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).Z0/sqrt(trackStates.at(trackIndex).covMatrix[9]));
    else result.push_back(std::nan(""));
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).phi);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[2]);
    else result.push_back(-9);
  }
  return result;
}


ROOT::VecOps::RVec<float>
getRP2TRK_omega(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size()){
      float omega = trackStates.at(trackIndex).omega;
      result.push_back(omega);
    }
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_omega_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[5]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_tanLambda(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).tanLambda);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_tanLambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[14]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_d0_phi0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[1]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_d0_omega_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        	    ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[3]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_d0_z0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        	 ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[6]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_d0_tanlambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        		ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[10]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi0_omega_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        	      ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[4]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi0_z0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        	   ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[7]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi0_tanlambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[11]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_omega_z0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        	    ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[8]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_omega_tanlambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        		   ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[12]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_z0_tanlambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
        		ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
        ROOT::VecOps::RVec<podio::ObjectID> reco2track_links) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    size_t trackIndex = getTrackIndex(p, reco2track_links);
    if (trackIndex<trackStates.size())
      result.push_back(trackStates.at(trackIndex).covMatrix[13]);
    else result.push_back(-9);
  }
  return result;
}


/*
Other undocumented garbage
*/

ROOT::VecOps::RVec<edm4hep::TrackState>
getRP2TRK( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					ROOT::VecOps::RVec<edm4hep::TrackState> tracks )
{

  ROOT::VecOps::RVec<edm4hep::TrackState> result ;
  result.reserve( in.size() );

  for (auto & p: in) {
    if (p.tracks_begin >= 0 && p.tracks_begin<tracks.size()) {
	result.push_back(tracks.at(p.tracks_begin) ) ;
    }
  }
  return result;
}

// returns reco indices of tracks
ROOT::VecOps::RVec<int> 
get_recoindTRK( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, 
		ROOT::VecOps::RVec<edm4hep::TrackState> tracks )
{

  ROOT::VecOps::RVec<int> result ;
  
  for (unsigned int ctr=0; ctr<in.size(); ctr++) {
    edm4hep::ReconstructedParticleData p = in.at(ctr);
    if (p.tracks_begin >= 0 && p.tracks_begin<tracks.size()) result.push_back(ctr) ;
  }
 return result ;
}

int getTK_n(ROOT::VecOps::RVec<edm4hep::TrackState> x) {
  int result =  x.size();
  return result;
}

///
ROOT::VecOps::RVec<bool> 
hasTRK( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in ) {

  ROOT::VecOps::RVec<bool> result ;
  result.reserve( in.size() );
  
  for (auto & p: in) {
    if (p.tracks_begin >= 0 && p.tracks_begin != p.tracks_end) result.push_back(true) ;
    else result.push_back(false);
  }
 return result ;
}

}//end NS ReconstructedParticle2Track

}//end NS FCCAnalyses
