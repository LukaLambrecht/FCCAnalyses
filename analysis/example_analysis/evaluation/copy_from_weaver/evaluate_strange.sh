# Keep track of commands run for final plots

# Compare following models:
# - Latest model (28/02) with V0 and dEdx (but masked in order to have good sim to data agreement).
# - Previous iteration (26/02) with dEdx (masked as above) but no V0.
# - Previous iteration (26/02) trained without dEdx or particle type variables.
# All of this in the conventional training and testing selection
# (does not correspond 100% with R_b or A_FB selection, but that's fine
# because only a relative comparison is made between these models, which are consistent)
python evaluate_strange.py -i ../models/output_20260228_withks_withdedx_masked/output.root ../models/output_20260226_withstrange_withdedx_masked/output.root ../models/output_20260226_withstrange_nodedx_noptype/output.root
