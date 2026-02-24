#include "FCCAnalyses/VertexFitterSimple.h"
#include "FCCAnalyses/MCParticle.h"

#include "edm4hep/EDM4hepVersion.h"

#include <iostream>

#include "TFile.h"
#include "TString.h"

#include <fcntl.h>


namespace FCCAnalyses {

namespace VertexFitterSimple {

int supress_stdout() {
  fflush(stdout);

  int ret = dup(1);
  int nullfd = open("/dev/null", O_WRONLY);
  // check nullfd for error omitted
  dup2(nullfd, 1);
  close(nullfd);

  return ret;
}

void resume_stdout(int fd) {
  fflush(stdout);
  dup2(fd, 1);
  close(fd);
  std::cout << std::flush;
}

// -----------------------------------------------------------------------------

VertexingUtils::FCCAnalysesVertex VertexFitter(
    int Primary,
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
    ROOT::VecOps::RVec<edm4hep::TrackState> thetracks, bool BeamSpotConstraint,
    double bsc_sigmax, double bsc_sigmay, double bsc_sigmaz, double bsc_x,
    double bsc_y, double bsc_z) {

  // input = a collection of recoparticles (in case one want to make
  // associations to RecoParticles ?) and thetracks = the collection of all
  // TrackState in the event

  VertexingUtils::FCCAnalysesVertex thevertex;

  // retrieve the tracks associated to the recoparticles
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks =
      ReconstructedParticle2Track::getRP2TRK(recoparticles, thetracks);

  // and run the vertex fitter

  // FCCAnalysesVertex thevertex = VertexFitter_Tk( Primary, tracks, thetracks)
  // ;
  thevertex =
      VertexFitter_Tk(Primary, tracks, thetracks, BeamSpotConstraint,
                      bsc_sigmax, bsc_sigmay, bsc_sigmaz, bsc_x, bsc_y, bsc_z);

  // fill the indices of the tracks
  ROOT::VecOps::RVec<int> reco_ind;
  int Ntr = tracks.size();
  for (auto &p : recoparticles) {
    // std::cout << " in VertexFitter:  a recoparticle with charge = " <<
    // p.charge << std::endl;
    if (p.tracks_begin >= 0 && p.tracks_begin < thetracks.size()) {
      reco_ind.push_back(p.tracks_begin);
    }
  }
  if (reco_ind.size() != Ntr)
    std::cout << " ... problem in Vertex, size of reco_ind != Ntr "
              << std::endl;

  thevertex.reco_ind = reco_ind;

  return thevertex;
}

// ---------------------------------------------------------------------------------------------------------------------------

VertexingUtils::FCCAnalysesVertex
VertexFitter_Tk(int Primary, ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
                bool BeamSpotConstraint, double bsc_sigmax, double bsc_sigmay,
                double bsc_sigmaz, double bsc_x, double bsc_y, double bsc_z,
                bool fast) {

  ROOT::VecOps::RVec<edm4hep::TrackState> dummy;
  return VertexFitter_Tk(Primary, tracks, dummy, BeamSpotConstraint, bsc_sigmax,
                         bsc_sigmay, bsc_sigmaz, bsc_x, bsc_y, bsc_z,
                         fast);
}

// ---------------------------------------------------------------------------------------------------------------------------

VertexingUtils::FCCAnalysesVertex
VertexFitter_Tk(int Primary, ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
                const ROOT::VecOps::RVec<edm4hep::TrackState> &alltracks,
                bool BeamSpotConstraint, double bsc_sigmax, double bsc_sigmay,
                double bsc_sigmaz, double bsc_x, double bsc_y, double bsc_z,
                bool fast) {
  // Base function for fitting a vertex to a collection of tracks

  // Suppressing printf() output from TMatrixBase:
  // https://github.com/root-project/root/blob/722eb4652bfc79149df00c8b0e92d0837caf054c/math/matrix/src/TMatrixTBase.cxx#L662
  // The solution found here:
  // https://stackoverflow.com/questions/46728680/how-to-temporarily-suppress-output-from-printf
  int fd = supress_stdout();

  // Units for the beam-spot : mum
  // See https://github.com/HEP-FCC/FCCeePhysicsPerformance/tree/master/
  // General#generating-events-under-realistic-fcc-ee-environment-conditions

  // initialization of result
  VertexingUtils::FCCAnalysesVertex TheVertex;

  // other initializations
  edm4hep::VertexData result;
  ROOT::VecOps::RVec<float> reco_chi2;
  ROOT::VecOps::RVec<TVectorD> updated_track_parameters;
  ROOT::VecOps::RVec<int> reco_ind;
  ROOT::VecOps::RVec<float> final_track_phases;
  ROOT::VecOps::RVec<TVector3> updated_track_momentum_at_vertex;

  // if the collection of all tracks has been passed, keep trace of the indices
  // of the tracks that are used to fit this vertex
  if (alltracks.size() > 0) {
    for (int i = 0; i < tracks.size(); i++) { // the fitted tracks
      edm4hep::TrackState tr1 = tracks[i];
      for (int j = 0; j < alltracks.size();
           j++) { // the collection of all tracks
        edm4hep::TrackState tr2 = alltracks[j];
        if (VertexingUtils::compare_Tracks(tr1, tr2)) {
          reco_ind.push_back(j);
          break;
        }
      }
    }
  }

  // initialize properties of the resulting FCCAnalysesVertex
  TheVertex.vertex = result;
  TheVertex.reco_chi2 = reco_chi2;
  TheVertex.updated_track_parameters = updated_track_parameters;
  TheVertex.reco_ind = reco_ind;
  TheVertex.final_track_phases = final_track_phases;
  TheVertex.updated_track_momentum_at_vertex = updated_track_momentum_at_vertex;

  // safety for too few provided tracks
  int Ntr = tracks.size();
  TheVertex.ntracks = Ntr;
  if (Ntr <= 1) {
    resume_stdout(fd);
    return TheVertex;
  }

  TVectorD **trkPar = new TVectorD *[Ntr];
  TMatrixDSym **trkCov = new TMatrixDSym *[Ntr];

  bool Units_mm = true;

  for (Int_t i = 0; i < Ntr; i++) {
    edm4hep::TrackState t = tracks[i];
    TVectorD par = VertexingUtils::get_trackParam(t, Units_mm);
    trkPar[i] = new TVectorD(par);
    TMatrixDSym Cov = VertexingUtils::get_trackCov(t, Units_mm);
    trkCov[i] = new TMatrixDSym(Cov);
  }

  // case of fast fit
  if( fast ){

    // initialize the fitter object (does not perform the fit yet)
    VertexFitFast theVertexFit(Ntr, trkPar, trkCov);

    // perform the fit
    TVectorD x = theVertexFit.GetVtx();

    // get position (in mm)
    result.position = edm4hep::Vector3f(x(0), x(1), x(2));

    // store the results in an edm4hep::VertexData object
    float Chi2 = theVertexFit.GetVtxChi2();
    float Ndof = 2.0 * Ntr - 3.0;
    result.chi2 = Chi2 / Ndof;

    // the chi2 of all the tracks
    TVectorD tracks_chi2 = theVertexFit.GetVtxChi2List();
    for (int it = 0; it < Ntr; it++) {
        reco_chi2.push_back(tracks_chi2[it]);
    }
  } // end of fast fit case

  else{

    // initialize the fitter object (does not perform the fit yet)
    VertexFit theVertexFit(Ntr, trkPar, trkCov);

    // add beamspot constraints to fitter object
    if (BeamSpotConstraint) {
        float conv_BSC = 1e-3;
        // (conversion factor from micrometer to millimeter;
        // the beamspot positions and widths are expected to be given in micrometers,
        // while the track- and vertex properties are typically in millimeter.)
        TVectorD xv_BS(3);
        xv_BS[0] = bsc_x * conv_BSC;
        xv_BS[1] = bsc_y * conv_BSC;
        xv_BS[2] = bsc_z * conv_BSC;
        TMatrixDSym cov_BS(3);
        cov_BS[0][0] = pow(bsc_sigmax * conv_BSC, 2);
        cov_BS[1][1] = pow(bsc_sigmay * conv_BSC, 2);
        cov_BS[2][2] = pow(bsc_sigmaz * conv_BSC, 2);
        theVertexFit.AddVtxConstraint(xv_BS, cov_BS);
    }

    // perform the fit
    TVectorD x = theVertexFit.GetVtx();

    // get position (in mm)
    result.position = edm4hep::Vector3f(x(0), x(1), x(2));

    // store the results in an edm4hep::VertexData object
    float Chi2 = theVertexFit.GetVtxChi2();
    float Ndof = 2.0 * Ntr - 3.0;
    result.chi2 = Chi2 / Ndof;

    // the chi2 of all the tracks
    TVectorD tracks_chi2 = theVertexFit.GetVtxChi2List();
    for (int it = 0; it < Ntr; it++) {
        reco_chi2.push_back(tracks_chi2[it]);
    }

    TMatrixDSym covX = theVertexFit.GetVtxCov();
    std::array<float, 6> covMatrix; // covMat in edm4hep is a LOWER-triangle matrix.
    covMatrix[0] = covX(0, 0);
    covMatrix[1] = covX(1, 0);
    covMatrix[2] = covX(1, 1);
    covMatrix[3] = covX(2, 0);
    covMatrix[4] = covX(2, 1);
    covMatrix[5] = covX(2, 2);
    result.covMatrix = covMatrix;

    // get more information from the vertex fit
    // note: up until this point, units didn't matter much as long as they are internally consistent;
    //       but the code below assumes mm, and this is important for correct momentum propagation.
    //       need to get the fit with correct units, but I don't want to break the code above,
    //       so best solution seems to re-do the vertex fit with modified input parameters (although this is slow)...
    TVectorD **trkPar_2 = new TVectorD *[Ntr];
    TMatrixDSym **trkCov_2 = new TMatrixDSym *[Ntr];
    bool Units_mm = true;
    for (Int_t i = 0; i < Ntr; i++) {
        edm4hep::TrackState t = tracks[i];
        TVectorD par = VertexingUtils::get_trackParam(t, Units_mm);
        // modify units
        par[0] *= 10; // d0 in cm instead of mm
        par[2] /= 10; // omega in 1/cm instead of 1/mm
        par[3] *= 10; // z0 in cm instead of mm
        // extra: the VertexMore class seems to have a magnetic field of 2T hard-coded, so modify omega accordingly...
        par[2] *= (2/1.5);
        trkPar_2[i] = new TVectorD(par);
        TMatrixDSym Cov = VertexingUtils::get_trackCov(t, Units_mm);
        trkCov_2[i] = new TMatrixDSym(Cov);
    }
    // note: no beamspot constraints yet as this part is typically only used for secondary vertices
    VertexingUtils::FCCAnalysesVertex TheVertex_2;
    VertexFit theVertexFit_2(Ntr, trkPar_2, trkCov_2);
    TVectorD x_2 = theVertexFit_2.GetVtx();
    // make extra info
    VertexMore theVertexMore(&theVertexFit_2, Units_mm);

    // copy properties
    for (Int_t i = 0; i < Ntr; i++) {
        TVectorD updated_par = theVertexFit.GetNewPar(i); // use the old fit (in original units)
        TVectorD updated_par_edm4hep = VertexingUtils::Delphes2Edm4hep_TrackParam(updated_par, Units_mm);
        updated_track_parameters.push_back(updated_par_edm4hep);
        TVector3 ptrack_at_vertex = theVertexMore.GetMomentum(i); // use the new fit (in correct units)
        updated_track_momentum_at_vertex.push_back(ptrack_at_vertex);
    }
  
    TheVertex.updated_track_parameters = updated_track_parameters;
    TheVertex.updated_track_momentum_at_vertex = updated_track_momentum_at_vertex;
    TheVertex.final_track_phases = final_track_phases;

    // memory cleanup
    for (Int_t i = 0; i < Ntr; i++) {
        delete trkPar_2[i];
        delete trkCov_2[i];
    }
    delete[] trkPar_2;
    delete[] trkCov_2;
  } // end case of full fit

#if EDM4HEP_BUILD_VERSION <= EDM4HEP_VERSION(0, 10, 5)
  result.primary = Primary;
#else
  result.type = Primary; // NOTE: Here we are relying on users passing in the correct value
#endif

  // add final parameters
  result.algorithmType = 1;
  TheVertex.vertex = result;
  TheVertex.reco_chi2 = reco_chi2;

  // memory cleanup
  for (Int_t i = 0; i < Ntr; i++) {
    delete trkPar[i];
    delete trkCov[i];
  }
  delete[] trkPar;
  delete[] trkCov;

  resume_stdout(fd);

  return TheVertex;
}

// ---------------------------------------------------------------------------------------------------------------------------

ROOT::VecOps::RVec<edm4hep::TrackState>
get_PrimaryTracks(ROOT::VecOps::RVec<edm4hep::TrackState> tracks, double chi2max,
                  bool BeamSpotConstraint, double bsc_sigmax, double bsc_sigmay,
                  double bsc_sigmaz, double bsc_x, double bsc_y, double bsc_z) {

  // iterative procedure to determine the primary vertex - and the primary tracks

  // Units for the beam-spot : mum
  // See
  // https://github.com/HEP-FCC/FCCeePhysicsPerformance/tree/master/General#generating-events-under-realistic-fcc-ee-environment-conditions

  // bool debug  = true ;
  bool debug = false;

  if (debug) {
    std::cout << " ... enter in VertexFitterSimple::get_PrimaryTracks   Ntr = "
              << tracks.size() << std::endl;
  }

  ROOT::VecOps::RVec<edm4hep::TrackState> seltracks = tracks;

  if (seltracks.size() <= 1)
    return seltracks;

  int Ntr = tracks.size();

  TVectorD **trkPar = new TVectorD *[Ntr];
  TMatrixDSym **trkCov = new TMatrixDSym *[Ntr];

  for (Int_t i = 0; i < Ntr; i++) {
    edm4hep::TrackState t = tracks[i];
    TVectorD par = VertexingUtils::get_trackParam(t);
    trkPar[i] = new TVectorD(par);
    TMatrixDSym Cov = VertexingUtils::get_trackCov(t);
    trkCov[i] = new TMatrixDSym(Cov);
  }

  VertexFit theVertexFit(Ntr, trkPar, trkCov);

  if (BeamSpotConstraint) {
    float conv_BSC = 1e-3;
    // (conversion factor from micrometer to millimeter;
    // the beamspot positions and widths are expected to be given in micrometers,
    // while the track- and vertex properties are typically in millimeter.)
    TVectorD xv_BS(3);
    xv_BS[0] = bsc_x * conv_BSC;
    xv_BS[1] = bsc_y * conv_BSC;
    xv_BS[2] = bsc_z * conv_BSC;
    TMatrixDSym cov_BS(3);
    cov_BS[0][0] = pow(bsc_sigmax * conv_BSC, 2);
    cov_BS[1][1] = pow(bsc_sigmay * conv_BSC, 2);
    cov_BS[2][2] = pow(bsc_sigmaz * conv_BSC, 2);
    theVertexFit.AddVtxConstraint(xv_BS, cov_BS);
  }

  TVectorD x = theVertexFit.GetVtx(); // this actually runs the fit

  float chi2_max = 1e30;

  while (chi2_max >= chi2max) {

    TVectorD tracks_chi2 = theVertexFit.GetVtxChi2List();
    chi2_max = tracks_chi2.Max();

    int n_removed = 0;
    for (int i = 0; i < theVertexFit.GetNtrk(); i++) {
      float track_chi2 = tracks_chi2[i];
      if (track_chi2 >= chi2_max) {
        theVertexFit.RemoveTrk(i);
        seltracks.erase(seltracks.begin() + i);
        n_removed++;
      }
    }
    if (n_removed > 0) {
      if (theVertexFit.GetNtrk() > 1) {
        // run the fit again:
        x = theVertexFit.GetVtx();
        TVectorD new_tracks_chi2 = theVertexFit.GetVtxChi2List();
        chi2_max = new_tracks_chi2.Max();
      } else {
        chi2_max = 0; // exit from the loop w/o crashing..
      }
    }
  } // end while

  // memory cleanup :
  for (Int_t i = 0; i < Ntr; i++) {
    delete trkPar[i];
    delete trkCov[i];
  }
  delete[] trkPar;
  delete[] trkCov;

  return seltracks;
}

// ---------------------------------------------------------------------------------------------------------------------------

ROOT::VecOps::RVec<edm4hep::TrackState>
get_NonPrimaryTracks(ROOT::VecOps::RVec<edm4hep::TrackState> allTracks,
                     ROOT::VecOps::RVec<edm4hep::TrackState> primaryTracks) {

  ROOT::VecOps::RVec<edm4hep::TrackState> result;
  for (auto &track : allTracks) {
    bool isInPrimary = false;
    for (auto &primary : primaryTracks) {
      if (VertexingUtils::compare_Tracks(track, primary)) {
        isInPrimary = true;
        break;
      }
    }
    if (!isInPrimary)
      result.push_back(track);
  }

  return result;
}

// ---------------------------------------------------------------------------------------------------------------------------

ROOT::VecOps::RVec<bool>
IsPrimary_forTracks(ROOT::VecOps::RVec<edm4hep::TrackState> allTracks,
                    ROOT::VecOps::RVec<edm4hep::TrackState> primaryTracks) {

  ROOT::VecOps::RVec<bool> result;
  for (auto &track : allTracks) {
    bool isInPrimary = false;
    for (auto &primary : primaryTracks) {
      if (VertexingUtils::compare_Tracks(track, primary)) {
        isInPrimary = true;
        break;
      }
    }
    result.push_back(isInPrimary);
  }
  return result;
}

} // namespace VertexFitterSimple

} // namespace FCCAnalyses
