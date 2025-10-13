import numpy as np

def westfall_young_fwer(stats, null_stats_matrix, alpha=0.05):
    """
    stats: array of observed statistics (larger is more significant)
    null_stats_matrix: (B, m) matrix, each row is a permutation's stats
    """
    stats = np.asarray(stats)
    B, m = null_stats_matrix.shape
    # max-T method
    max_null = null_stats_matrix.max(axis=1)
    pvals = np.array([(np.sum(max_null >= s) + 1) / (B + 1) for s in stats])
    return pvals, pvals <= alpha
