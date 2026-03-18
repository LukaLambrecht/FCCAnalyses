# Measure the beamspot position in data

Notes on synchronization:
The results in `output_legacy` which are currently used in the analysis are slightly outdated in 2 ways:
- The chi2 cutoff value in the primary track selection is 25 (old) instead of 5 (new).
  The value of 25 is still used in the vertex fit applied in this analysis script;
  it can be changed at any point (for consistency with the main analysis),
  but then the results might come out slightly different...
- They were calculated before propagating the sign flip of D0 and omega to the track covariance matrix.
  This propagation is currently enabled by default, hence new runs will give slightly different results,
  unless this propagation is temporarily disabled again (verified that in this way the legacy results can be reproduced).
