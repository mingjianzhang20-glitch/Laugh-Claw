"""Utilities for reporting Boltz2 continuous affinity outputs."""


def predicted_nominal_ic50_nm(affinity_pred_value):
    """Convert log10(IC50 in micromolar) to predicted nominal IC50 in nM."""
    return 10 ** float(affinity_pred_value) * 1000.0
