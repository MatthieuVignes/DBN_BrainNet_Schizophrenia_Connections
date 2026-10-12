"""Toy calibration check of group-comparison tests on bootstrap edge counts.
Under a true null (both groups drawn from the same population) we compare:
 (a) Wilcoxon / mean-of-replicates on bootstrap counts  -> anti-conservative
 (b) percentile interval of all pairwise bootstrap differences -> near nominal
Run: python3 26_toy_bootstrap_calibration.py   (numpy, scipy)"""
import numpy as np
from scipy.stats import mannwhitneyu
rng = np.random.default_rng(1)
N, B, SIMS = 30, 100, 400
def boot_counts(x):  # statistic = mean (stand-in for edge count)
    return np.array([rng.choice(x, len(x)).mean() for _ in range(B)])
fp_w = fp_i = 0
for _ in range(SIMS):
    a, b = rng.normal(size=N), rng.normal(size=N)
    ba, bb = boot_counts(a), boot_counts(b)
    if mannwhitneyu(ba, bb).pvalue < .05: fp_w += 1
    dif = (ba[:, None] - bb[None, :]).ravel()
    lo, hi = np.percentile(dif, [2.5, 97.5])
    if lo > 0 or hi < 0: fp_i += 1
print(f"Wilcoxon on replicates   : false-positive rate {fp_w/SIMS:.3f}")
print(f"Pairwise-difference interval: false-positive rate {fp_i/SIMS:.3f}")
