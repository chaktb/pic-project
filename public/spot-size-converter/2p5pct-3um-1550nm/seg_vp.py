"""
seg_vp.py -- widening, variable-(minimum-)pitch segmented SSC  (2.5% Delta, 3x3 um, 1550 nm)
=========================================================================================
Structure (same family as the 1.5% 5x5 um page, adapted to a high-contrast core):

  chip WG 3 um solid  ->  segments whose WIDTH grows gradually  3 um -> W_facet
                          while the DUTY falls  d_start -> d_facet
                          on the NARROWEST pitch the mask rule allows locally

Local pitch rule (period k, local duty d):
      P_k = max( P_min ,  g_min / (1-d) ,  t_min / d )     (t_min = g_min = 0.5 um)
  then nudged upward out of the contra-directional Bragg bands
      m * lambda / (2 n_eff),  m = 2..6,  n_eff in [n_clad, n_eff(chip)]  (+-0.03 um guard)
  -> dense segmentation where the duty is mid-range (largest scattering otherwise),
     every tooth and every gap >= 0.5 um.

Duty / width laws (s = normalized position 0..1 along the segmented section):
      d(s) = d_start + (d_facet - d_start) * s**p_d
      w(s) = W0      + (W_facet - W0)      * s**p_w      (W_facet > W0: widening)

Duck-types SSCDesign (teeth(), lead_um, drc(), to_dict()) so ssc_core.segmented_slices
and the BPM work unchanged.
"""
import math
from dataclasses import dataclass
import numpy as np

LAM = 1.55
N_LO, N_HI = 1.4440, 1.4625


def bragg_bands(m_max=6, guard=0.03):
    return [(m * LAM / (2 * N_HI) - guard, m * LAM / (2 * N_LO) + guard) for m in range(2, m_max + 1)]


BANDS = bragg_bands()


def avoid_bragg(P):
    for lo, hi in BANDS:
        if lo <= P <= hi:
            return hi + 1e-3
    return P


def ceil_grid(v, g=0.001):
    return math.ceil(v / g - 1e-9) * g


@dataclass
class WideningSSC:
    L_um: float = 600.0          # segmented-section length
    d_start: float = 0.80
    d_facet: float = 0.22
    w_facet_um: float = 7.0
    p_d: float = 1.0
    p_w: float = 1.0
    P_min_um: float = 1.20       # global pitch floor (> 2nd-order Bragg band)
    lead_um: float = 20.0
    W0_um: float = 3.0
    t_min: float = 0.5
    g_min: float = 0.5

    def _law(self, s):
        s = min(max(s, 0.0), 1.0)
        d = self.d_start + (self.d_facet - self.d_start) * s ** self.p_d
        w = self.W0_um + (self.w_facet_um - self.W0_um) * s ** self.p_w
        return d, w

    def periods(self):
        """List of (z_start, pitch, tooth, w0, w1) along the segmented section."""
        out, z = [], 0.0
        while z < self.L_um - 1e-9:
            d, _ = self._law(z / self.L_um)
            P = max(self.P_min_um, self.g_min / max(1 - d, 1e-3), self.t_min / max(d, 1e-3))
            P = avoid_bragg(P)
            t = max(self.t_min, d * P)
            g = P - t
            if g < self.g_min:
                g = self.g_min; P = t + g
            t, P = ceil_grid(t), ceil_grid(P)
            if P - t < self.g_min:
                P = ceil_grid(t + self.g_min)
            _, w0 = self._law(z / self.L_um)
            _, w1 = self._law((z + t) / self.L_um)
            out.append((z, P, t, w0, w1))
            z += P
        return out

    def teeth(self, grid_snap=0.001):
        segs = [(0.0, self.lead_um, self.W0_um, self.W0_um)]
        for z, P, t, w0, w1 in self.periods():
            za = self.lead_um + z + (P - t)     # gap first, then tooth
            segs.append((round(za, 3), round(za + t, 3), round(w0, 3), round(w1, 3)))
        return segs                 # facet = end of the last tooth

    def drc(self):
        segs = self.teeth()
        tooth = min(s[1] - s[0] for s in segs[1:])
        gaps = [segs[i + 1][0] - segs[i][1] for i in range(len(segs) - 1)]
        P = [p[1] for p in self.periods()]
        return dict(min_tooth=round(tooth, 4), min_gap=round(min(gaps), 4),
                    min_width=round(min(min(s[2], s[3]) for s in segs), 4),
                    max_width=round(max(max(s[2], s[3]) for s in segs), 4),
                    n_teeth=len(segs) - 1, pitch_min=round(min(P), 4), pitch_max=round(max(P), 4),
                    length=round(segs[-1][1], 3))

    def to_dict(self):
        return dict(type="widening variable-pitch", L_um=self.L_um, d_start=self.d_start,
                    d_facet=self.d_facet, w_facet_um=self.w_facet_um, p_d=self.p_d, p_w=self.p_w,
                    P_min_um=self.P_min_um, lead_um=self.lead_um, W0_um=self.W0_um,
                    t_min=self.t_min, g_min=self.g_min)


NAMES = ["L_um", "d_start", "d_facet", "w_facet_um", "p_d", "p_w", "P_min_um"]


def from_genes(v):
    return WideningSSC(L_um=v[0], d_start=v[1], d_facet=v[2], w_facet_um=v[3], p_d=v[4], p_w=v[5],
                       P_min_um=v[6])
