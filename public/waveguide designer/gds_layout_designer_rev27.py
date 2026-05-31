#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Waveguide Layout Designer  (OlympiOS-style, Tkinter + gdsfactory)
=================================================================

Three-pane GUI:
    [Left]   Layers + Ports + Globals + Groups + Auto-Offset
    [Center] Real-time layout canvas (matplotlib)
    [Right]  Segment sequence editor

Workflow
--------
1. Add layers (name, GDS layer/datatype, colour, default width)
2. Add one or more "base" ports (user-defined anchors, e.g. chip I/O)
3. Build the waveguide as a sequence of segments:
       LINE / ARC / SBEND / REL_BENT / EULER / TAPER
   Each segment starts at a selected port number and generates a new
   endpoint port whose number is appended to the port table.
   The next segment simply picks that new number as its start.
4. Canvas redraws instantly on every change.
5. (New in rev25) Enable "Auto Offset" tab → bend-induced mode-center
   offset δ = c·W²/(R·Δ) is auto-injected at every straight↔bend and
   bend↔bend junction.  Platform params (W, Δ) and coefficient are
   editable in the UI.
6. Export to GDS via gdsfactory (add_polygon).

(New in rev27) REL_BENT segment — an "offset bend" built from two equal
circular arcs of a fixed radius R that reaches a lateral offset `dy`
while the heading returns to the start angle (the bend always ENDS at
0° relative to the start port).  θ = arccos(1 − |dy|/(2R)).

Requires:  python3, tkinter, numpy, matplotlib, gdsfactory
    pip install numpy matplotlib gdsfactory

Author: generated for Sangjun (FiberPro PDL) — rev 27  (rel_bent offset bend)
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Tuple

import numpy as np

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk,
)

# gdsfactory is only required for GDS export; delay the import so the GUI
# still launches on machines where it isn't installed.
try:
    import gdsfactory as gf
    _HAS_GDSFACTORY = True
except Exception:
    _HAS_GDSFACTORY = False


# ============================================================================
#  Data model
# ============================================================================

@dataclass
class Layer:
    name: str
    number: int                       # GDS layer number
    datatype: int = 0
    color: str = "#1f77b4"
    width: float = 2.0                # default waveguide width [µm]


@dataclass
class Port:
    pid: int                          # 1-indexed
    x: float                          # µm
    y: float                          # µm
    angle: float                      # degrees (0 = +x, 90 = +y)
    width: float                      # µm
    layer_idx: int                    # index into layers list
    label: str = ""
    is_base: bool = True              # True if user-defined, False if auto


@dataclass
class Segment:
    sid: int                          # 1-indexed
    kind: str                         # 'line' | 'arc' | 'sbend' | 'rel_bent' | 'euler' | 'taper'
    start_pid: int                    # starting port id
    end_pid: int = 0                  # populated when segment is built
    layer_idx: int = 0                # which layer this segment uses
    params: Dict = field(default_factory=dict)


# ----------------------------------------------------------------------------
#  Auto bend-offset  (rev25)
# ----------------------------------------------------------------------------
@dataclass
class PlatformConfig:
    """Silica PLC / LNOI / Ti:LN 등 도파로 플랫폼 파라미터.

    Bent waveguide 에서 기본모드의 peak 는 curvature center 반대쪽으로
    δ 만큼 shift 된다. 접합점에서 두 segment 의 mode center 를 정렬하려면
    곡선 쪽을 lateral 로 δ 만큼 상대 이동시켜야 한다.

        δ(R) = offset_coeff · W² / (R · Δ)           [기본 공식]
        또는 사용자가 offset_fn(R) 로 오버라이드 가능.
    """
    width:   float = 2.0           # core width W [μm]
    delta:   float = 0.0075        # relative index contrast  (0.75 %)
    offset_coeff: float = 1.0 / 8.0   # 기본 1/8  (typical silica PLC)
    enabled: bool  = False         # UI 체크박스와 연결

    def bend_offset(self, R: float) -> float:
        """Lateral mode-center offset at a bend of radius R [μm]."""
        if R is None or R <= 0 or self.delta <= 0 or self.width <= 0:
            return 0.0
        return self.offset_coeff * self.width ** 2 / (R * self.delta)


def _resolve_num(v, globals_dict: Optional[Dict] = None) -> Optional[float]:
    """Segment.params 의 값(숫자 또는 Python 식 문자열)을 float 으로 변환.
    실패 시 None."""
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        try:
            return float(_parse_num(v, globals_dict or {}))
        except Exception:
            return None
    return None


def _segment_signed_curvature(seg: Optional[Segment],
                               at: str,
                               globals_dict: Optional[Dict] = None) -> float:
    """접합부에서 segment 의 signed curvature (1/R, +CCW/−CW) 를 반환.
    straight-like 종단이면 0.
       at='start'  : segment 의 시작단
       at='end'    : segment 의 끝단
    """
    if seg is None:
        return 0.0
    k = seg.kind
    # straight-like segments — 양 끝단 κ = 0
    # (rel_bent 은 sbend 처럼 시작·끝 heading 이 동일하므로 straight-like 취급)
    if k in ("line", "taper", "sbend", "rel_bent", "spiral"):
        return 0.0
    # Euler (half-clothoid): start 는 κ=0, end 는 κ=1/R
    if k == "euler":
        if at == "start":
            return 0.0
        R = _resolve_num(seg.params.get("radius"), globals_dict)
        a = _resolve_num(seg.params.get("angle"),  globals_dict)
        if R is None or a is None or R <= 0 or a == 0:
            return 0.0
        return math.copysign(1.0 / R, a)
    # ARC / tapered_arc — 양 끝단 κ = 1/R 로 동일
    if k in ("arc", "tapered_arc"):
        R = _resolve_num(seg.params.get("radius"), globals_dict)
        a = _resolve_num(seg.params.get("angle"),  globals_dict)
        if R is None or a is None or R <= 0 or a == 0:
            return 0.0
        return math.copysign(1.0 / R, a)
    return 0.0


def junction_offset(prev_seg: Optional[Segment],
                    curr_seg: Segment,
                    platform: PlatformConfig,
                    globals_dict: Optional[Dict] = None) -> float:
    """prev_seg → curr_seg 접합에서 curr_seg 에 적용할 lateral offset [μm].

    부호 규약: + 는 curr_seg 의 시작포트 local frame 에서 '왼쪽'(진행방향
    기준 +90°) 방향 이동. 이것은 segment_geometry 의 `yoffset` 과 정확히
    동일한 frame 이므로 `params['yoffset']` 에 그대로 더해주면 된다.

        +sign(κ_curr_start) · δ(R_curr)      ← curvature center 쪽으로 이동
        −sign(κ_prev_end)   · δ(R_prev)      ← prev bend mode center 보정
    """
    if not platform.enabled:
        return 0.0
    off = 0.0
    k_c = _segment_signed_curvature(curr_seg, "start", globals_dict)
    k_p = _segment_signed_curvature(prev_seg, "end",   globals_dict)
    if k_c != 0.0:
        R_c = 1.0 / abs(k_c)
        off += math.copysign(platform.bend_offset(R_c), k_c)
    if k_p != 0.0:
        R_p = 1.0 / abs(k_p)
        off += -math.copysign(platform.bend_offset(R_p), k_p)
    return off


# ============================================================================
#  Geometry engine
# ============================================================================

def _rot(ang_deg: float) -> Tuple[float, float]:
    r = math.radians(ang_deg)
    return math.cos(r), math.sin(r)


# ============================================================================
#  Spiral-style helpers (see segment_geometry for the "spiral" kind)
# ============================================================================

def _spiral_straight(x0, y0, a_deg, L, n=30):
    a = math.radians(a_deg)
    s = np.linspace(0, L, n)
    return x0 + s * math.cos(a), y0 + s * math.sin(a)


def _spiral_arc(cx, cy, R, phi0_deg, sweep_deg, n=120):
    phi = np.radians(np.linspace(phi0_deg, phi0_deg + sweep_deg, n))
    return cx + R * np.cos(phi), cy + R * np.sin(phi)


def _spiral_hermite(A, tA, B, tB, mag, n=120):
    Ax, Ay = A; Bx, By = B
    tAx, tAy = tA; tBx, tBy = tB
    t = np.linspace(0, 1, n)
    h00 =  2*t**3 - 3*t**2 + 1
    h10 =    t**3 - 2*t**2 + t
    h01 = -2*t**3 + 3*t**2
    h11 =    t**3 -   t**2
    xs = h00*Ax + h10*(mag*tAx) + h01*Bx + h11*(mag*tBx)
    ys = h00*Ay + h10*(mag*tAy) + h01*By + h11*(mag*tBy)
    return xs, ys


def _spiral_finalise(X, Y):
    """Rotate so entry tangent = +x local, translate so entry at origin,
    return (xl, yl, (end_x, end_y, end_angle_degrees))."""
    if len(X) < 2:
        return X, Y, (0.0, 0.0, 0.0)
    t0 = math.atan2(Y[1] - Y[0], X[1] - X[0])
    c, s = math.cos(-t0), math.sin(-t0)
    Xr = c*X - s*Y;  Yr = s*X + c*Y
    Xr -= Xr[0];     Yr -= Yr[0]
    dex = Xr[-1] - Xr[-2]; dey = Yr[-1] - Yr[-2]
    end_a = math.degrees(math.atan2(dey, dex))
    return Xr, Yr, (float(Xr[-1]), float(Yr[-1]), end_a)


def _spiral_archimedean(R_min, pitch, n_turns, n_pts=1200):
    """Style 'archimedean' — Image 1 yin-yang double spiral.

    Pixel-accurate match to gdsfactory's ``spiral_double``. Construction:

      1. Arm A: r(θ) = pitch/π · θ + R_min,  pos = r·(sinθ, cosθ),
         θ ∈ [0, 2πN].  Inner end (θ=0) at (0, R_min).
         Traversed outer → inner so we end at the inner end heading -x.
      2. Arm B: arm A mirrored across the x-axis (y → -y).
         Inner end at (0, -R_min). Traversed inner → outer.
      3. Central S-connector = TWO semicircles of radius R_min/2:
         • First: centre (0, +R_min/2), swept CCW 90°→270° (LEFT bulge),
           connects (0, R_min) to (0, 0).
         • Second: centre (0, -R_min/2), swept CW 90°→-90° (RIGHT bulge),
           connects (0, 0) to (0, -R_min).
         The two opposite-side bulges make the yin-yang shape.

    No rotation is applied (would destroy the yin-yang orientation) —
    only translation so entry is at the origin.
    """
    a_coeff = pitch / math.pi
    tmax = 2 * math.pi * n_turns
    N = max(int(300 * n_turns), 400)

    # Arm A: outer → inner
    th = np.linspace(tmax, 0, N)
    r  = a_coeff * th + R_min
    xA = r * np.sin(th)
    yA = r * np.cos(th)

    # Arm B: mirror across x-axis, inner → outer
    th2 = np.linspace(0, tmax, N)
    r2  = a_coeff * th2 + R_min
    xB =  r2 * np.sin(th2)
    yB = -r2 * np.cos(th2)

    # First S-semicircle: left bulge (CCW from +90° to +270°)
    k = max(120, n_pts // 6)
    phi1 = np.radians(np.linspace(90, 270, k))
    xS1 = (R_min/2) * np.cos(phi1)
    yS1 = R_min/2 + (R_min/2) * np.sin(phi1)
    # Second S-semicircle: right bulge (CW from +90° to -90°)
    phi2 = np.radians(np.linspace(90, -90, k))
    xS2 = (R_min/2) * np.cos(phi2)
    yS2 = -R_min/2 + (R_min/2) * np.sin(phi2)

    X = np.concatenate([xA, xS1[1:], xS2[1:], xB[1:]])
    Y = np.concatenate([yA, yS1[1:], yS2[1:], yB[1:]])

    # Translate only — DO NOT rotate; rotating destroys the yin-yang.
    X = X - X[0]
    Y = Y - Y[0]
    dex = X[-1] - X[-2]; dey = Y[-1] - Y[-2]
    end_a = math.degrees(math.atan2(dey, dex))
    return X, Y, (float(X[-1]), float(Y[-1]), end_a)


def _spiral_racetrack(R_min, pitch, n_turns, straight_L=None, n_pts=1200):
    """Style 'racetrack' — Image 2: interleaved double-serpentine with
    pitch/2 gap between outgoing and return arms. Outgoing winds
    upward through N loops, a small central U-turn of radius pitch/4
    flips direction, and the return arm unwinds in the gaps.

    Input enters heading -x (from east). Exit emerges heading +x on
    the west side after unwinding (canonical). The io_orientation
    parameter selects which edges ports appear on.
    """
    if straight_L is None or straight_L <= 0:
        straight_L = 10 * R_min
    N = max(1, int(round(n_turns)))
    px, py = [], []
    x, y, head = 0.0, 0.0, 180.0
    xs, ys = _spiral_straight(x, y, head, straight_L)
    px.append(xs); py.append(ys); x, y = xs[-1], ys[-1]

    # Outgoing spiral: N loops winding upward
    for k in range(N):
        R = R_min + k * pitch
        cx, cy = x, y + R
        xs, ys = _spiral_arc(cx, cy, R, 270, -180)
        px.append(xs[1:]); py.append(ys[1:]); x, y = xs[-1], ys[-1]
        head = 0.0
        xs, ys = _spiral_straight(x, y, head, straight_L)
        px.append(xs[1:]); py.append(ys[1:]); x, y = xs[-1], ys[-1]
        if k < N-1:
            R2 = R + pitch
            cx, cy = x, y + R2
            xs, ys = _spiral_arc(cx, cy, R2, -90, 180)
            px.append(xs[1:]); py.append(ys[1:]); x, y = xs[-1], ys[-1]
            head = 180.0
            xs, ys = _spiral_straight(x, y, head, straight_L)
            px.append(xs[1:]); py.append(ys[1:]); x, y = xs[-1], ys[-1]

    # Central U-turn of radius pitch/4: flips +x → -x, drops y by pitch/2
    R_c = pitch / 4.0
    cx, cy = x, y - R_c
    xs, ys = _spiral_arc(cx, cy, R_c, 90, -180)
    px.append(xs[1:]); py.append(ys[1:]); x, y = xs[-1], ys[-1]
    head = 180.0

    # Return arm: mirror of outgoing, offset by pitch/2 downward.
    for k in range(N):
        k_r = N - 1 - k
        xs, ys = _spiral_straight(x, y, head, straight_L)
        px.append(xs[1:]); py.append(ys[1:]); x, y = xs[-1], ys[-1]
        if k < N-1:
            R = R_min + k_r * pitch + pitch/2.0
            cx, cy = x, y - R
            xs, ys = _spiral_arc(cx, cy, R, 90, -180)
            px.append(xs[1:]); py.append(ys[1:]); x, y = xs[-1], ys[-1]
            head = 0.0
            xs, ys = _spiral_straight(x, y, head, straight_L)
            px.append(xs[1:]); py.append(ys[1:]); x, y = xs[-1], ys[-1]
            R2 = max(R - pitch, R_min - pitch/2)
            cx, cy = x, y - R2
            xs, ys = _spiral_arc(cx, cy, R2, 90, 180)
            px.append(xs[1:]); py.append(ys[1:]); x, y = xs[-1], ys[-1]
            head = 180.0
    X = np.concatenate(px); Y = np.concatenate(py)
    return _spiral_finalise(X, Y)


def _spiral_racetrack_s(R_min, pitch, n_turns, straight_L=None, n_pts=1200,
                         target_length=None):
    """Style 'racetrack_s' — Image 3. Calls gdsfactory's
    ``spiral_racetrack_fixed_length`` directly and returns the
    centreline as a single connected path.

    Since gdsfactory returns a multi-polygon component (not a single
    centerline), we approximate the centreline by sampling evenly
    along the component's bounding box and using the known topology.
    For an accurate polygon output, the gdsfactory path centerline
    is built analytically here — matching gdsfactory's geometric
    recipe rather than calling it at runtime (to avoid gdsfactory
    as a hard runtime dependency for the designer tool).

    Parameters
    ----------
    R_min       : min_radius for U-turn bends
    pitch       : min_spacing between nested racetracks
    n_turns     : sets n_straight_sections = 2·n_turns
    straight_L  : if ``target_length`` is None, used as port spacing
    target_length : if given, this is the total centreline length
                    (solve for straight_L internally)
    """
    # Map our naming to gdsfactory's:
    #   in_out_port_spacing = straight_L
    #   n_straight_sections = 2·n_turns
    #   min_radius = R_min
    #   min_spacing = pitch
    n_sections = max(2, 2 * int(round(n_turns)))
    port_spacing = straight_L if straight_L and straight_L > 0 else 150.0

    # Compute length contributions for each element:
    # - Each nested racetrack (N_nest = n_sections//2) contributes:
    #     2 bends (180° each) of radius R_k + 2 straights of length L_k
    # - R_k = R_min + k·pitch for k=0..N_nest-1
    # - The outermost straight length ≈ port_spacing, shrinking inward by 2·pitch each nest
    # - Plus a central euler bend_s S-connector
    N_nest = n_sections // 2
    R_outer = R_min + (N_nest - 1) * pitch

    # If target_length given, solve for L_outer (outermost straight length)
    # so that total centerline = target_length
    def total_length(L_outer):
        L = 0.0
        for k in range(N_nest):
            R_k = R_min + k * pitch
            L_k = max(pitch, L_outer - 2 * pitch * k)
            L += 2 * math.pi * R_k      # 2 × 180° = full turn
            L += 2 * L_k                # top + bottom straights
        # central S-bend: approx π·R_min (one full 180° equivalent)
        L += math.pi * R_min
        return L

    if target_length is not None and target_length > 0:
        # Bisection on L_outer
        lo, hi = 2 * pitch, target_length
        for _ in range(40):
            mid = (lo + hi) / 2
            if total_length(mid) < target_length:
                lo = mid
            else:
                hi = mid
        L_outer = mid
    else:
        L_outer = port_spacing

    # Now BUILD the geometry as a single continuous centerline.
    # Layout: nested racetracks in the LEFT half plane (x < some centre),
    # exit port on the RIGHT via an S-bend.
    #
    # Following gdsfactory's visual: input at left-edge of outermost racetrack
    # going +x, goes around the outermost racetrack, jumps inward via a small
    # offset, next racetrack, etc. Central connector flips direction and
    # output goes straight out to the east.

    # We build from the INSIDE out for simplicity. Innermost racetrack
    # (smallest one), then spiral outward. Exit on top-east going +x.
    px, py = [], []
    def _add(xs, ys, first=False):
        if first:
            px.append(xs); py.append(ys)
        else:
            px.append(xs[1:]); py.append(ys[1:])

    # Start at innermost. Each racetrack has length straight_k:
    #   inner (k=0): length L_inner = L_outer - 2·pitch·(N_nest-1)
    #   outer (k=N_nest-1): length L_outer
    L_inner = max(pitch, L_outer - 2 * pitch * (N_nest - 1))

    # Start at centre of innermost racetrack, build the whole stack from
    # the innermost straight going +x (bottom of the innermost oval).
    # Start the spiral at the LEFT edge of innermost oval bottom.
    x, y, head = 0.0, 0.0, 0.0
    # First: a short "central S-bend" entry -- two half-circles of radius R_min/2
    # that go up and come back, making the figure-8 in the middle.
    # Actually simpler: place the central S-bend at the start, then spiral out.
    k_s = max(60, n_pts // 6)
    R_c = R_min / 2
    # Start heading +x at (0, 0). Central S: two semicircles flipping direction
    # CW arc up → CCW arc down → end heading -x at (0, R_min)
    # Then we go into the spiral.
    phi1 = np.radians(np.linspace(-90, 90, k_s))
    xS1 = (R_c) * np.cos(phi1)
    yS1 = R_c + R_c * np.sin(phi1)
    _add(xS1, yS1, first=True); x, y = xS1[-1], yS1[-1]; head = 180.0
    phi2 = np.radians(np.linspace(-90, -270, k_s))
    xS2 = x + R_c * np.cos(phi2)
    yS2 = y + R_c + R_c * np.sin(phi2)
    _add(xS2, yS2); x, y = xS2[-1], yS2[-1]; head = 180.0

    # Now we're at (0, R_min) heading -x. Enter the innermost racetrack.
    # Go left along bottom of innermost, then up-and-over, then right along top,
    # spiralling outward.
    for k in range(N_nest):
        R_k = R_min + k * pitch
        L_k = max(pitch, L_outer - 2 * pitch * (N_nest - 1 - k))
        # Bottom straight heading -x of length L_k
        xs, ys = _spiral_straight(x, y, head, L_k)
        _add(xs, ys); x, y = xs[-1], ys[-1]
        # Left U-turn UP (heading -x → +x, centre above)
        cx, cy = x, y + R_k
        xs, ys = _spiral_arc(cx, cy, R_k, 270, -180)   # CW 270 → 90
        _add(xs, ys); x, y = xs[-1], ys[-1]; head = 0.0
        # Top straight heading +x
        xs, ys = _spiral_straight(x, y, head, L_k)
        _add(xs, ys); x, y = xs[-1], ys[-1]
        if k < N_nest - 1:
            # Right U-turn UP (+x → -x) — larger radius for next lap
            R_next = R_k + pitch
            cx, cy = x, y + R_next
            xs, ys = _spiral_arc(cx, cy, R_next, -90, 180)   # CCW -90 → 90
            _add(xs, ys); x, y = xs[-1], ys[-1]; head = 180.0

    # Exit straight heading +x to the east edge (after the final top-straight)
    xs, ys = _spiral_straight(x, y, head, L_outer * 0.3)
    _add(xs, ys)

    X = np.concatenate(px); Y = np.concatenate(py)
    X = X - X[0]; Y = Y - Y[0]
    dex = X[-1] - X[-2]; dey = Y[-1] - Y[-2]
    end_a = math.degrees(math.atan2(dey, dex))
    return X, Y, (float(X[-1]), float(Y[-1]), end_a)


def _spiral_archimedean_cross(R_min, pitch, n_turns, n_pts=1200):
    """Style 'archimedean_x' — Image 4: single-arm Archimedean spiral
    that winds INWARD from the top, makes a 180° U-turn at the centre,
    then exits radially through a crossing. In the local frame entry
    is at origin heading +x; with a start port angle of -90° it
    renders with the spiral body HANGING below the entry edge."""
    a_coeff = pitch / (2 * math.pi)
    tmax = 2 * math.pi * n_turns
    N = max(int(300 * n_turns), 400)
    # Build spiral in a frame where entry is at (r_max, 0) heading +y
    th = np.linspace(tmax, 0, N)     # outer → inner
    r  = R_min + a_coeff * th
    x = r * np.cos(th)
    y = r * np.sin(th)
    # Rotate 90° CCW so entry is at (0, r_max) heading +x (right):
    xA = -y
    yA =  x
    # Inner end: the rotated point.
    inner_phi = math.atan2(yA[-1], xA[-1])
    # U-turn: half-circle of radius R_min centred at origin, sweeping
    # CW by 180° from inner end.
    phi_u = np.linspace(inner_phi, inner_phi - math.pi, 120)
    xU = R_min * np.cos(phi_u)
    yU = R_min * np.sin(phi_u)
    # Radial straight OUTWARD from U-turn exit, along its angle.
    end_angle = inner_phi - math.pi
    r_edge = R_min + a_coeff * tmax + 2 * R_min
    r_cross = np.linspace(R_min, r_edge, 80)
    xs_cross = r_cross * math.cos(end_angle)
    ys_cross = r_cross * math.sin(end_angle)
    X = np.concatenate([xA, xU[1:], xs_cross[1:]])
    Y = np.concatenate([yA, yU[1:], ys_cross[1:]])
    return _spiral_finalise(X, Y)


def _spiral_path_length(xl, yl):
    return float(np.sum(np.hypot(np.diff(xl), np.diff(yl))))


def _spiral_solve_length(style, R_min, pitch, target_length,
                         n_turns_init=3.0, straight_L_init=None,
                         tune="auto"):
    """Tune a spiral parameter to hit ``target_length`` total length.

    Parameters
    ----------
    style : one of "archimedean", "racetrack", "racetrack_s", "archimedean_x"
    R_min, pitch : fixed physical constraints [µm]
    target_length : desired total path length [µm]
    n_turns_init : starting guess for n_turns when tuning straight_L
    straight_L_init : starting guess for straight_L (racetrack styles)
    tune : 'n_turns' → solve for n_turns
           'straight_L' → solve for straight_L (racetrack styles only)
           'auto' → n_turns for Archimedean, straight_L for racetrack

    Returns
    -------
    (params_dict, L_achieved)  or  (None, None) if out of range
    """
    if tune == "auto":
        tune = "n_turns" if "archimedean" in style else "straight_L"
    if tune == "straight_L" and "archimedean" in style:
        raise ValueError(
            "'straight_L' tuning is only for 'racetrack' styles")

    def _L(n, sL):
        try:
            if style in ("archimedean", "racetrack", "racetrack_s"):
                # gdsfactory-backed: measure polygon perimeter/2
                params = {"radius_min": R_min, "pitch": pitch, "n_turns": n}
                if sL and sL > 0:
                    params["straight_L"] = sL
                polys, _ = _spiral_gdsfactory(style, params)
                return _polys_length(polys)
            else:
                xl, yl, _ = _spiral_dispatch(style, R_min, pitch, n, sL, 1200)
                return _spiral_path_length(xl, yl)
        except Exception:
            return -1.0

    # Simple bisection
    if tune == "n_turns":
        lo, hi = 0.5, 100.0
        sL = straight_L_init or 10*R_min
        L_lo = _L(lo, sL); L_hi = _L(hi, sL)
        if L_lo > target_length or L_hi < target_length:
            return None, None
        for _ in range(60):
            mid = (lo + hi) * 0.5
            L_mid = _L(mid, sL)
            if abs(L_mid - target_length) < 0.5:
                return ({"n_turns": mid}, L_mid)
            if L_mid < target_length:
                lo, L_lo = mid, L_mid
            else:
                hi, L_hi = mid, L_mid
        return ({"n_turns": mid}, L_mid)
    else:  # straight_L
        lo, hi = 2*R_min, 1e7   # up to 10 m
        n = n_turns_init
        L_lo = _L(n, lo); L_hi = _L(n, hi)
        if L_lo > target_length or L_hi < target_length:
            return None, None
        for _ in range(60):
            mid = (lo + hi) * 0.5
            L_mid = _L(n, mid)
            if abs(L_mid - target_length) < 0.5:
                return ({"straight_L": mid, "n_turns": n}, L_mid)
            if L_mid < target_length:
                lo, L_lo = mid, L_mid
            else:
                hi, L_hi = mid, L_mid
        return ({"straight_L": mid, "n_turns": n}, L_mid)


def _spiral_dispatch(style, R_min, pitch, n_turns, straight_L, n_pts):
    """Route to the appropriate spiral-style generator. Returns
    (xl, yl, (end_x, end_y, end_angle_local_deg))."""
    if style == "archimedean":
        return _spiral_archimedean(R_min, pitch, n_turns, n_pts)
    if style == "racetrack":
        return _spiral_racetrack(R_min, pitch, n_turns, straight_L, n_pts)
    if style == "racetrack_s":
        return _spiral_racetrack_s(R_min, pitch, n_turns, straight_L, n_pts)
    if style == "archimedean_x":
        return _spiral_archimedean_cross(R_min, pitch, n_turns, n_pts)
    raise ValueError(f"Unknown spiral style: {style!r}")


# ---------------------------------------------------------------------------
# gdsfactory-backed spiral geometry
#
# For the archimedean, racetrack, and racetrack_s styles we call the
# official gdsfactory components directly and extract their polygons.
# This guarantees exact geometric fidelity with the gdsfactory library
# (see https://gdsfactory.github.io/gdsfactory/notebooks/04_components_shapes.html).
#
# Returns:
#   polygons     : list of (N_i, 2) numpy arrays in the LOCAL frame
#                  (entry of the spiral at the origin, pointing +x)
#   end_local    : (ex, ey, ea_deg) — output port in local frame
# ---------------------------------------------------------------------------
_GF_DBU = 0.001   # KLayout database unit → µm

def _spiral_gdsfactory(style, params):
    """Build a spiral with gdsfactory and return (polygons_local, end_local).
    Entry port (o1) is translated+rotated to origin with +x heading.

    Parameters understood per-style:
      archimedean:
        min_bend_radius = radius_min   (µm, default 10)
        separation      = pitch        (µm, default 3)
        number_of_loops = n_turns      (int, default 4)
      racetrack:
        length          = straight_L   (µm — target total length)
        spacing         = pitch
        n_loops         = n_turns
      racetrack_s:
        straight_length = straight_L
        spacings        = (pitch,) * n_turns  (one per nested lap)
        min_radius      = radius_min
    """
    import gdsfactory as gf
    try:
        gf.gpdk.PDK.activate()
    except Exception:
        pass   # PDK already active — fine

    R_min   = float(params.get("radius_min", 10))
    pitch   = float(params.get("pitch", 3))
    n_turns = float(params.get("n_turns", 4))
    sL_raw  = params.get("straight_L", 0)
    try:
        straight_L = float(sL_raw)
    except (TypeError, ValueError):
        straight_L = 0

    # --- build the gdsfactory component ---
    if style == "archimedean":
        comp = gf.components.spiral_double(
            min_bend_radius=R_min,
            separation=pitch,
            number_of_loops=max(1, int(round(n_turns))),
        )
    elif style == "racetrack":
        # gdsfactory's spiral() uses length (total) + spacing + n_loops.
        # Use straight_L as total length if provided, else derive from
        # R_min/pitch/n_turns heuristically.
        if straight_L > 0:
            total_len = straight_L
        else:
            total_len = max(20.0, 20 * pitch * n_turns)
        comp = gf.components.spiral(
            length=total_len,
            spacing=pitch,
            n_loops=max(1, int(round(n_turns))),
        )
    elif style == "racetrack_s":
        n_laps = max(1, int(round(n_turns)))
        sl = straight_L if straight_L > 0 else max(20.0, 15 * pitch)
        comp = gf.components.spiral_racetrack(
            straight_length=sl,
            spacings=tuple([pitch] * n_laps),
            min_radius=R_min,
        )
    else:
        raise ValueError(f"Style {style!r} not supported by gdsfactory backend")

    # --- extract polygons (in µm) ---
    polys = []
    for poly_list in comp.get_polygons(by="tuple").values():
        for p in poly_list:
            hull = np.array(
                [(pt.x * _GF_DBU, pt.y * _GF_DBU) for pt in p.each_point_hull()]
            )
            polys.append(hull)

    # --- find entry port (o1) and output port (o2) ---
    port_map = {p.name: p for p in comp.ports}
    p_in  = port_map.get("o1") or port_map.get("in")
    p_out = port_map.get("o2") or port_map.get("out")
    if p_in is None or p_out is None:
        # fall back: use bbox corners
        bb = comp.bbox_np()
        p_in_center  = (float(bb[0][0]), 0.0)
        p_in_angle   = 180.0
        p_out_center = (float(bb[1][0]), 0.0)
        p_out_angle  = 0.0
    else:
        p_in_center  = (float(p_in.center[0]),  float(p_in.center[1]))
        p_in_angle   = float(p_in.orientation)
        p_out_center = (float(p_out.center[0]), float(p_out.center[1]))
        p_out_angle  = float(p_out.orientation)

    # gdsfactory port orientation points OUTWARD (away from component).
    # An external waveguide entering this port travels in the direction
    # (orientation + 180°).  In our local frame, entry heading should be +x,
    # so rotate by: theta = 180° - p_in_angle  (mapping p_in_angle+180° → 0°).
    theta_deg = 180.0 - p_in_angle
    theta = math.radians(theta_deg)
    c, s = math.cos(theta), math.sin(theta)
    tx, ty = -p_in_center[0], -p_in_center[1]   # translate so p_in → origin

    def _xform(x, y):
        xt = x + tx; yt = y + ty
        xr = c * xt - s * yt
        yr = s * xt + c * yt
        return xr, yr

    # transform all polygons
    out_polys = []
    for P in polys:
        xs = P[:, 0] + tx; ys = P[:, 1] + ty
        xr = c * xs - s * ys
        yr = s * xs + c * ys
        out_polys.append(np.column_stack([xr, yr]))

    # transform output port
    ex, ey = _xform(*p_out_center)
    # Output port orientation in LOCAL frame: the waveguide LEAVES via
    # the orientation direction (same convention for external continuation).
    # Our convention: end_angle is the heading of waveguide going OUT.
    ea_local = p_out_angle + theta_deg
    # Normalize to [-180, 180]
    while ea_local > 180:  ea_local -= 360
    while ea_local <= -180: ea_local += 360

    return out_polys, (float(ex), float(ey), float(ea_local))


def _polys_length(polys):
    """Estimate total waveguide centerline length from polygon list.
    We approximate by half the total perimeter (a thin rectangle's
    perimeter ≈ 2·L + 2·w ≈ 2·L, so centerline length ≈ perim/2)."""
    total_perim = 0.0
    for P in polys:
        d = np.diff(P, axis=0, append=P[:1])
        total_perim += float(np.sum(np.hypot(d[:, 0], d[:, 1])))
    return total_perim / 2.0


def segment_geometry(
    kind: str,
    start: Port,
    params: Dict,
    n_pts: int = 160,
    globals_dict: Optional[Dict] = None,
) -> Tuple[np.ndarray, np.ndarray, Tuple[float, float, float, float]]:
    """
    Build a waveguide segment starting from a port.

    `params` values may be either numbers (float / int) OR strings
    containing Python expressions that reference project globals
    (e.g. "R", "width*2", "(a+b)/2"). Strings are evaluated against
    ``globals_dict`` every time this function runs, which lets the
    user tweak globals in the UI and have the geometry regenerate
    automatically without overwriting the original expressions.

    `params["width"]` (optional, >0) overrides the start-port width.
    For TAPER this is the *starting* width; `params["end_width"]` remains
    the ending width. If omitted or 0, the segment inherits the start
    port's width (back-compat with earlier project files).

    `params["xoffset"]`, `params["yoffset"]`, `params["angleoffset"]`
    (all optional, default 0) displace the segment's *effective start*
    relative to the start port's LOCAL frame:
        xoffset → along the start-port heading
        yoffset → 90° CCW of the heading (perpendicular, "left")
        angleoffset → added to the heading [deg]
    This lets you branch several segments off the same port without
    needing extra dummy base ports.

    Returns
    -------
    path_xy : (N,2) ndarray       centreline in global coordinates
    widths  : (N,) ndarray        width at each sample
    end     : (x, y, angle, width) of the endpoint port (global frame)
    """
    # --- resolve any string-valued params against globals -------------
    # Some param keys (like 'style') are string selectors, not numeric
    # expressions — those are passed through unchanged.
    _STRING_KEYS = {"style"}
    G = globals_dict or {}
    resolved: Dict = {}
    for k, v in params.items():
        if k in _STRING_KEYS:
            resolved[k] = v
            continue
        if isinstance(v, str):
            try:
                resolved[k] = _parse_num(v, G)
            except ValueError as e:
                raise ValueError(
                    f"cannot evaluate params[{k!r}]={v!r}: {e}") from None
        else:
            resolved[k] = v
    params = resolved   # shadow with all-numeric version below

    # --- offsets in the START port's local frame ----------------------
    def _p(key):
        try:
            return float(params.get(key, 0.0))
        except (TypeError, ValueError):
            return 0.0
    xoff = _p("xoffset")
    yoff = _p("yoffset")
    aoff = _p("angleoffset")

    ca_s, sa_s = _rot(start.angle)
    # translate (xoff, yoff) in the start port's local frame, then rotate
    x0 = start.x + xoff * ca_s - yoff * sa_s
    y0 = start.y + xoff * sa_s + yoff * ca_s
    a0 = start.angle + aoff
    ca, sa = _rot(a0)

    # width override (0 / missing / invalid → inherit start port)
    try:
        w_over = float(params.get("width", 0.0))
    except (TypeError, ValueError):
        w_over = 0.0
    w0 = w_over if w_over > 0 else start.width

    # --- local coordinates (x along start direction, y perpendicular) ----
    if kind == "line":
        L = float(params["length"])
        t = np.linspace(0.0, L, n_pts)
        xl = t
        yl = np.zeros_like(t)
        widths = np.full_like(t, w0)
        end_xl, end_yl = L, 0.0
        end_a = a0
        end_w = w0

    elif kind == "taper":
        L = float(params["length"])
        w1 = float(params["end_width"])
        t = np.linspace(0.0, L, n_pts)
        xl = t
        yl = np.zeros_like(t)
        widths = np.linspace(w0, w1, n_pts)
        end_xl, end_yl = L, 0.0
        end_a = a0
        end_w = w1

    elif kind == "arc":
        R = float(params["radius"])
        dth = float(params["angle"])                     # signed deg, +CCW/left
        s = 1.0 if dth >= 0 else -1.0
        phi = np.linspace(0.0, math.radians(abs(dth)), n_pts)
        # centre is perpendicular to start heading, on the +y side when s=+1
        xl = R * np.sin(phi)
        yl = s * R * (1.0 - np.cos(phi))
        widths = np.full_like(phi, w0)
        end_xl = R * math.sin(math.radians(abs(dth)))
        end_yl = s * R * (1.0 - math.cos(math.radians(abs(dth))))
        end_a = a0 + dth
        end_w = w0

    elif kind == "tapered_arc":
        # Circular bend whose width varies linearly from `w0` (start) to
        # `end_width` (end) along the arc. Centreline is a pure circular
        # arc; the polygon is obtained by offsetting ±w(s)/2 in
        # path_to_polygon, which naturally gives different inner/outer
        # radii at each end. Common in ring-coupler tapers and adiabatic
        # mode converters.
        R   = float(params["radius"])
        dth = float(params["angle"])
        w1  = float(params["end_width"])
        s   = 1.0 if dth >= 0 else -1.0
        phi = np.linspace(0.0, math.radians(abs(dth)), n_pts)
        xl  = R * np.sin(phi)
        yl  = s * R * (1.0 - np.cos(phi))
        widths = np.linspace(w0, w1, n_pts)       # linear in arc-length
        end_xl = R * math.sin(math.radians(abs(dth)))
        end_yl = s * R * (1.0 - math.cos(math.radians(abs(dth))))
        end_a = a0 + dth
        end_w = w1

    elif kind == "sbend":
        # Cosine S-bend — C¹ continuous, zero end-slope, low loss.
        dx = float(params["dx"])
        dy = float(params["dy"])
        t = np.linspace(0.0, 1.0, n_pts)
        xl = dx * t
        yl = (dy / 2.0) * (1.0 - np.cos(np.pi * t))
        widths = np.full_like(t, w0)
        end_xl, end_yl = dx, dy
        end_a = a0
        end_w = w0

    elif kind == "rel_bent":
        # Relative bend ("offset bend"): two equal circular arcs of fixed
        # radius R that reach a target y while the heading returns to the
        # start heading (net angle change = 0, so the bend always ENDS at
        # 0° relative to the start port).
        #
        #   first arc  : curvature +s/R, turns by +s·θ
        #   second arc : curvature −s/R, turns by −s·θ  → back to start angle
        #
        # `dy` is interpreted as the ABSOLUTE target y coordinate (not a
        # relative shift). The lateral displacement actually applied is
        #       h_signed = dy − y0        (y0 = effective start y)
        # so e.g. dy = yoffset(5) routes this bend to the same absolute y
        # as segment 5's endpoint. (Absolute-y targeting assumes a roughly
        # horizontal start heading, the usual rel_bent use case.)
        #
        # With h = |h_signed|:
        #       h = 2 R (1 − cos θ)   →   θ = arccos(1 − h/(2R))
        # The per-arc angle θ is fully determined by R and h.
        # Reachable range: h ≤ 4R (θ ≤ 180°); larger values are clamped.
        R  = float(params["radius"])
        dy_target = float(params["dy"])     # absolute target y
        h_signed  = dy_target - y0          # lateral displacement to apply
        s  = 1.0 if h_signed >= 0 else -1.0
        h  = abs(h_signed)
        if R <= 0 or h == 0.0:
            # degenerate → straight stub along the start heading
            xl = np.array([0.0, 0.0])
            yl = np.array([0.0, h_signed])
            widths = np.array([w0, w0])
            end_xl, end_yl = 0.0, h_signed
        else:
            c = 1.0 - min(h, 4.0 * R) / (2.0 * R)
            c = max(-1.0, min(1.0, c))
            theta = math.acos(c)
            half = max(2, n_pts // 2)
            phi = np.linspace(0.0, theta, half)
            # first arc (curving toward +s side)
            x1 = R * np.sin(phi)
            y1 = s * R * (1.0 - np.cos(phi))
            # junction state (end of first arc)
            xm, ym = float(x1[-1]), float(y1[-1])
            am = s * theta
            ca2, sa2 = math.cos(am), math.sin(am)
            # second arc: opposite curvature, built in the junction's local
            # frame then rotated/translated into the segment frame
            xl2 = R * np.sin(phi)
            yl2 = -s * R * (1.0 - np.cos(phi))
            x2 = xm + xl2 * ca2 - yl2 * sa2
            y2 = ym + xl2 * sa2 + yl2 * ca2
            xl = np.concatenate([x1, x2[1:]])
            yl = np.concatenate([y1, y2[1:]])
            widths = np.full_like(xl, w0)
            end_xl, end_yl = float(xl[-1]), float(yl[-1])
        end_a = a0                 # heading returns to the start heading
        end_w = w0

    elif kind == "euler":
        # Half-clothoid: κ(s) ramps linearly from 0 to 1/R over length L=R·θ.
        R = float(params["radius"])
        dth = float(params["angle"])                     # total bend [deg]
        s_sign = 1.0 if dth >= 0 else -1.0
        th_rad = math.radians(abs(dth))
        L = R * th_rad                                    # half-clothoid
        if L <= 0 or R <= 0:
            xl = np.zeros(n_pts); yl = np.zeros(n_pts)
            widths = np.full(n_pts, w0); end_xl = 0.0; end_yl = 0.0
        else:
            s_arr = np.linspace(0.0, L, n_pts)
            theta = s_sign * s_arr ** 2 / (2.0 * R * L) * th_rad * 2.0  # θ(L)=θ
            # ^^ ensures θ at s=L equals s_sign*th_rad
            dxl = np.cos(theta); dyl = np.sin(theta)
            xl = np.concatenate(([0.0], np.cumsum(
                (dxl[:-1] + dxl[1:]) * 0.5 * np.diff(s_arr))))
            yl = np.concatenate(([0.0], np.cumsum(
                (dyl[:-1] + dyl[1:]) * 0.5 * np.diff(s_arr))))
            widths = np.full_like(s_arr, w0)
            end_xl, end_yl = float(xl[-1]), float(yl[-1])
        end_a = a0 + dth
        end_w = w0

    elif kind == "spiral":
        # Parametric spiral delay line. Four styles are supported via
        # the ``style`` param:
        #   "archimedean"   — Image 1: gdsfactory.components.spiral_double (yin-yang)
        #   "racetrack"     — Image 2: gdsfactory.components.spiral (rectangular expanding)
        #   "racetrack_s"   — Image 3: gdsfactory.components.spiral_racetrack
        #                     (racetrack + central S-bend)
        #   "archimedean_x" — Image 4: self-implemented Archimedean with
        #                     radial crossing (no gdsfactory equivalent)
        #
        # Styles 1-3 call the official gdsfactory component via
        # _spiral_gdsfactory(), extract the polygons, and place them
        # into the `extra_polygons` slot — bypassing the centreline
        # model entirely (since gdsfactory components are multi-polygon
        # and not trivially reducible to a single centreline).
        R_min   = float(params["radius_min"])
        pitch   = float(params["pitch"])
        n_turns = float(params["n_turns"])
        style   = params.get("style", "archimedean")
        try:
            straight_L = float(params.get("straight_L", 0))
        except (TypeError, ValueError):
            straight_L = 0
        if straight_L <= 0:
            straight_L = None

        if R_min <= 0 or pitch <= 0 or n_turns <= 0:
            xl = np.zeros(n_pts); yl = np.zeros(n_pts)
            widths = np.full(n_pts, w0); end_xl = end_yl = 0.0
            end_a = a0; end_w = w0
            extra_polys = None
        elif style == "archimedean_x":
            # Self-implemented (gdsfactory has no crossing spiral)
            xl, yl, (end_xl, end_yl, end_a_local) = _spiral_archimedean_cross(
                R_min, pitch, n_turns, n_pts)
            widths = np.full_like(xl, w0)
            end_a = a0 + end_a_local
            end_w = w0
            extra_polys = None
        else:
            # Delegate to gdsfactory. Returns polygons in local frame
            # (entry at origin, heading +x) + output port location.
            local_polys, (end_xl, end_yl, end_a_local) = \
                _spiral_gdsfactory(style, {**params, "straight_L": straight_L or 0})
            # No centreline (drawn as polygons only). Keep a placeholder
            # 2-point "path" so that downstream code that inspects path[]
            # (e.g. annotation midpoint) still works.
            xl = np.array([0.0, end_xl])
            yl = np.array([0.0, end_yl])
            widths = np.array([w0, w0])
            end_a = a0 + end_a_local
            end_w = w0
            extra_polys = local_polys

    else:
        raise ValueError(f"Unknown segment kind: {kind!r}")

    # --- rotate local → global -------------------------------------------
    X = x0 + xl * ca - yl * sa
    Y = y0 + xl * sa + yl * ca
    path_xy = np.column_stack([X, Y])

    end_x = x0 + end_xl * ca - end_yl * sa
    end_y = y0 + end_xl * sa + end_yl * ca

    # Transform extra_polygons (if any) from local to global frame
    extra_polys_global = None
    if kind == "spiral" and 'extra_polys' in locals() and extra_polys is not None:
        extra_polys_global = []
        for P in extra_polys:
            Xg = x0 + P[:, 0] * ca - P[:, 1] * sa
            Yg = y0 + P[:, 0] * sa + P[:, 1] * ca
            extra_polys_global.append(np.column_stack([Xg, Yg]))

    return path_xy, widths, (end_x, end_y, end_a, end_w), extra_polys_global


def path_to_polygon(path_xy: np.ndarray, widths: np.ndarray) -> np.ndarray:
    """Offset centreline ±w/2 perpendicular to tangent → closed polygon."""
    dx = np.gradient(path_xy[:, 0])
    dy = np.gradient(path_xy[:, 1])
    mag = np.hypot(dx, dy)
    mag[mag == 0] = 1.0
    tx, ty = dx / mag, dy / mag
    # normal (rotated +90°)
    nx, ny = -ty, tx
    half = (widths / 2.0)[:, None]
    left = path_xy + np.column_stack([nx, ny]) * half
    right = path_xy - np.column_stack([nx, ny]) * half
    return np.vstack([left, right[::-1]])


# ----------------------------------------------------------------------------
#  Hit-testing helpers  (used by canvas click → segment pick in rev25)
# ----------------------------------------------------------------------------
def _point_in_polygon(x: float, y: float, poly: np.ndarray) -> bool:
    """Ray-casting point-in-polygon test. `poly` is (N, 2)."""
    if poly is None or len(poly) < 3:
        return False
    n = len(poly)
    inside = False
    x0, y0 = poly[-1]
    for i in range(n):
        x1, y1 = poly[i]
        if ((y1 > y) != (y0 > y)):
            xint = x1 + (y - y1) * (x0 - x1) / (y0 - y1 + 1e-300)
            if x < xint:
                inside = not inside
        x0, y0 = x1, y1
    return inside


def _min_dist2_to_polyline(x: float, y: float, path: np.ndarray) -> float:
    """Squared perpendicular distance from (x, y) to the closest segment
    of a polyline `path` (N, 2). Returns a large value for empty input.
    """
    if path is None or len(path) < 2:
        return float("inf")
    P = path
    x0 = P[:-1, 0]; y0 = P[:-1, 1]
    x1 = P[1:,  0]; y1 = P[1:,  1]
    dx = x1 - x0
    dy = y1 - y0
    seg_len2 = dx*dx + dy*dy
    seg_len2 = np.where(seg_len2 == 0, 1e-30, seg_len2)
    t = ((x - x0)*dx + (y - y0)*dy) / seg_len2
    t = np.clip(t, 0.0, 1.0)
    px = x0 + t*dx
    py = y0 + t*dy
    d2 = (px - x)**2 + (py - y)**2
    return float(d2.min())


# ============================================================================
#  Project (model holding all layers/ports/segments)
# ============================================================================

class Project:
    def __init__(self):
        self.layers: List[Layer] = [
            Layer("Core", 1, 0, "#1f77b4", 2.0),
        ]
        self.ports: List[Port] = []
        self.segments: List[Segment] = []
        self._next_pid = 1
        self._next_sid = 1

        # Seed one default base port at the origin so a brand-new project
        # has a valid start port (#1) for the very first "Add Segment".
        # NOTE: Project.from_dict() replaces self.ports wholesale, so this
        # seed never duplicates ports when loading a saved project.
        self.add_base_port(0.0, 0.0, 0.0, self.layers[0].width, 0, label="in")

        # Global variables — set by the user in the UI, auto-injected
        # into the script namespace on every run.
        self.globals: Dict[str, object] = {}

        # Groups — reusable "subcircuits" defined as Python callables in
        # the script (via the `@group("name")` decorator). Rebuilt every
        # script run; NOT serialised as callables.
        self._groups: Dict[str, object] = {}

        # (rev25) Auto bend-offset configuration.
        # When platform.enabled is True, rebuild() injects a mode-center
        # lateral offset at every straight↔bend and bend↔bend junction.
        self.platform: PlatformConfig = PlatformConfig()

        # Populated by rebuild() when auto-offset is enabled: list of dicts
        # { sid, start_pid, prev_sid, prev_kind, curr_kind, offset_um }.
        # Used by the UI to display junction offsets in a table.
        self.offset_log: List[Dict] = []

    # ---- ids ---------------------------------------------------------------
    def new_pid(self) -> int:
        pid = self._next_pid; self._next_pid += 1; return pid

    def new_sid(self) -> int:
        sid = self._next_sid; self._next_sid += 1; return sid

    # ---- layers ------------------------------------------------------------
    def add_layer(self, layer: Layer):
        self.layers.append(layer)

    def delete_layer(self, idx: int):
        if len(self.layers) <= 1:
            raise RuntimeError("At least one layer must remain.")
        if any(p.layer_idx == idx for p in self.ports) or \
           any(s.layer_idx == idx for s in self.segments):
            raise RuntimeError("Layer is in use by ports/segments.")
        self.layers.pop(idx)
        # re-index references above idx
        for p in self.ports:
            if p.layer_idx > idx: p.layer_idx -= 1
        for s in self.segments:
            if s.layer_idx > idx: s.layer_idx -= 1

    # ---- ports -------------------------------------------------------------
    def get_port(self, pid: int) -> Optional[Port]:
        for p in self.ports:
            if p.pid == pid: return p
        return None

    def add_base_port(self, x, y, angle, width, layer_idx, label="") -> Port:
        p = Port(self.new_pid(), x, y, angle, width, layer_idx, label, True)
        self.ports.append(p)
        return p

    def delete_port(self, pid: int):
        p = self.get_port(pid)
        if p is None: return
        if not p.is_base:
            raise RuntimeError("Auto-generated ports are deleted with their "
                               "segment.")
        if any(s.start_pid == pid for s in self.segments):
            raise RuntimeError("Port is used as the start of a segment.")
        self.ports = [q for q in self.ports if q.pid != pid]

    # ---- segments ----------------------------------------------------------
    def add_segment(self, kind: str, start_pid: int, layer_idx: int,
                    params: Dict, end_pid: Optional[int] = None) -> Segment:
        start = self.get_port(start_pid)
        if start is None:
            raise RuntimeError(f"Start port {start_pid} does not exist.")
        # endpoint port id: caller may request a specific number; otherwise
        # one is auto-assigned. A requested number must be free.
        if end_pid is None:
            end_pid = self.new_pid()
        else:
            if end_pid <= 0:
                raise RuntimeError("End port # must be a positive integer.")
            if self.get_port(end_pid) is not None:
                raise RuntimeError(
                    f"Port {end_pid} already exists; choose a free number "
                    f"for the end port.")
            # keep auto-numbering ahead of any user-chosen id
            self._next_pid = max(self._next_pid, end_pid + 1)
        seg = Segment(self.new_sid(), kind, start_pid, 0, layer_idx, dict(params))
        # create an endpoint port *now* and keep its pid stable forever
        _res = segment_geometry(
            kind, start, seg.params, globals_dict=self.eval_globals())
        ex, ey, ea, ew = _res[2]
        end_port = Port(end_pid, ex, ey, ea, ew, layer_idx,
                        label=f"S{seg.sid}", is_base=False)
        seg.end_pid = end_port.pid
        self.ports.append(end_port)
        self.segments.append(seg)
        self.rebuild()
        return seg

    def delete_segment(self, sid: int):
        seg = next((s for s in self.segments if s.sid == sid), None)
        if seg is None: return
        # orphan any segment that started from this one's endpoint
        for other in self.segments:
            if other.start_pid == seg.end_pid:
                other.start_pid = 0
        self.ports   = [p for p in self.ports   if p.pid != seg.end_pid]
        self.segments = [s for s in self.segments if s.sid != sid]
        self.rebuild()

    def move_segment(self, sid: int, delta: int):
        idx = next((i for i, s in enumerate(self.segments) if s.sid == sid), -1)
        if idx < 0: return
        new = max(0, min(len(self.segments) - 1, idx + delta))
        if new == idx: return
        self.segments[idx], self.segments[new] = self.segments[new], self.segments[idx]
        self.rebuild()

    # ---- rebuild: update positions of auto ports in chain order -----------
    def _prev_by_start(self) -> Dict[int, "Segment"]:
        """{start_pid → previous segment chained into that port}."""
        d: Dict[int, "Segment"] = {}
        for s in self.segments:
            if s.end_pid:
                d[s.end_pid] = s
        return d

    def _yoffset(self, sid) -> float:
        """Absolute global y of segment ``sid``'s endpoint port.

        Exposed to param expressions as ``yoffset(n)`` so a segment can
        reference another segment's absolute y — e.g. a rel_bent with
        ``dy = yoffset(5)`` routes to the same y as segment 5's end.
        The referenced segment must already be built (i.e. appear earlier
        in the sequence) for its endpoint y to be current.
        """
        try:
            sid = int(sid)
        except (TypeError, ValueError):
            raise ValueError(f"yoffset(): invalid segment id {sid!r}")
        seg = next((s for s in self.segments if s.sid == sid), None)
        if seg is None:
            raise ValueError(f"yoffset(): segment {sid} not found")
        port = self.get_port(seg.end_pid)
        if port is None:
            raise ValueError(f"yoffset(): segment {sid} has no endpoint port")
        return float(port.y)

    def eval_globals(self) -> Dict:
        """Globals namespace for param-expression evaluation: the user's
        globals plus built-in helper callables (currently ``yoffset``)."""
        ns = dict(self.globals)
        ns["yoffset"] = self._yoffset
        return ns

    def effective_params(self, seg: "Segment",
                          prev_map: Optional[Dict[int, "Segment"]] = None
                          ) -> Dict:
        """Return ``seg.params`` with the auto bend-offset merged into
        ``yoffset`` when ``self.platform.enabled`` is True. The original
        ``seg.params`` dict is never mutated.
        """
        eff = dict(seg.params)
        if not self.platform.enabled:
            return eff
        if prev_map is None:
            prev_map = self._prev_by_start()
        prev_seg = prev_map.get(seg.start_pid)
        G = self.eval_globals()
        auto_off = junction_offset(prev_seg, seg, self.platform,
                                    globals_dict=G)
        if auto_off != 0.0:
            user_yoff = _resolve_num(eff.get("yoffset", 0.0), G) or 0.0
            eff["yoffset"] = user_yoff + auto_off
        return eff

    def rebuild(self):
        """Recompute auto-port positions; NEVER reassign existing pids.

        (rev25) When ``self.platform.enabled`` is True, a lateral
        mode-center offset δ is injected at every junction where the
        previous segment (ending at the current segment's start port) is
        a straight↔bend or bend↔bend transition. The offset is added to
        ``params['yoffset']`` on-the-fly for the rebuild pass only — the
        user's stored params are NOT mutated.
        """
        by_pid = {p.pid: p for p in self.ports}
        prev_map = self._prev_by_start()
        G = self.eval_globals()

        self.offset_log = []

        for seg in self.segments:
            start = by_pid.get(seg.start_pid)
            if start is None:
                # dangling segment; leave its endpoint untouched
                continue

            # --- effective params with auto-offset injection -------------
            eff_params = self.effective_params(seg, prev_map)
            if self.platform.enabled and eff_params.get("yoffset") not in (
                    None, seg.params.get("yoffset")):
                # The effective yoffset differs from the stored one →
                # that means auto-offset was injected. Record for the UI.
                prev_seg = prev_map.get(seg.start_pid)
                user_yoff = _resolve_num(seg.params.get("yoffset", 0.0),
                                          G) or 0.0
                off_val = float(eff_params["yoffset"]) - user_yoff
                if abs(off_val) > 1e-15:
                    self.offset_log.append({
                        "sid":        seg.sid,
                        "start_pid":  seg.start_pid,
                        "prev_sid":   prev_seg.sid if prev_seg else None,
                        "prev_kind":  prev_seg.kind if prev_seg else "—",
                        "curr_kind":  seg.kind,
                        "offset_um":  off_val,
                    })

            try:
                _res = segment_geometry(
                    seg.kind, start, eff_params, globals_dict=G)
                ex, ey, ea, ew = _res[2]
            except Exception:
                continue
            end = by_pid.get(seg.end_pid)
            if end is None:
                # shouldn't normally happen; create one if missing
                end = Port(self.new_pid(), ex, ey, ea, ew, seg.layer_idx,
                           label=f"S{seg.sid}", is_base=False)
                self.ports.append(end)
                by_pid[end.pid] = end
                seg.end_pid = end.pid
            else:
                end.x = ex; end.y = ey
                end.angle = ea; end.width = ew
                end.layer_idx = seg.layer_idx
                end.label = f"S{seg.sid}"

    # ---- groups (reusable sub-circuits) -----------------------------------
    def register_group(self, name: str, func):
        """Register a Python callable as a named group."""
        self._groups[name] = func

    def call_group(self, name: str, *args, **kwargs):
        """Invoke a registered group. The group function is responsible
        for calling proj.add_segment(...) / proj.add_base_port(...) to
        actually build the geometry, using the arguments it receives."""
        if name not in self._groups:
            raise RuntimeError(
                f"Group {name!r} not defined. "
                f"Available: {sorted(self._groups.keys())}")
        return self._groups[name](*args, **kwargs)

    def list_groups(self) -> List[Tuple[str, str]]:
        """Return (name, signature) pairs for the defined groups."""
        import inspect
        out = []
        for name, fn in self._groups.items():
            try:
                out.append((name, str(inspect.signature(fn))))
            except (TypeError, ValueError):
                out.append((name, "(...)"))
        return out

    # ---- serialization -----------------------------------------------------
    def to_dict(self) -> Dict:
        # Only literal globals are serialisable via JSON → filter here.
        serialisable = {}
        for k, v in self.globals.items():
            try:
                import json as _json
                _json.dumps(v)
                serialisable[k] = v
            except TypeError:
                serialisable[k] = repr(v)
        return {
            "layers":   [asdict(l) for l in self.layers],
            "ports":    [asdict(p) for p in self.ports],   # include auto too
            "segments": [asdict(s) for s in self.segments],
            "globals":  serialisable,
            "platform": asdict(self.platform),              # (rev25)
            "next_pid": self._next_pid,
            "next_sid": self._next_sid,
        }

    @classmethod
    def from_dict(cls, d: Dict) -> "Project":
        proj = cls()
        proj.layers = [Layer(**l) for l in d.get("layers", [])]
        if not proj.layers:
            proj.layers = [Layer("Core", 1, 0, "#1f77b4", 2.0)]
        proj.ports = [Port(**p) for p in d.get("ports", [])]
        proj.segments = [Segment(**s) for s in d.get("segments", [])]
        proj.globals = dict(d.get("globals", {}))
        # (rev25) restore platform config if present
        plat_d = d.get("platform")
        if isinstance(plat_d, dict):
            try:
                proj.platform = PlatformConfig(**plat_d)
            except TypeError:
                # forward-compat: ignore unknown keys
                known = {k: plat_d[k] for k in plat_d
                         if k in PlatformConfig.__dataclass_fields__}
                proj.platform = PlatformConfig(**known)
        proj._next_pid = max(d.get("next_pid", 1),
                             max((p.pid for p in proj.ports), default=0) + 1)
        proj._next_sid = max(d.get("next_sid", 1),
                             max((s.sid for s in proj.segments), default=0) + 1)
        proj.rebuild()
        return proj


# ============================================================================
#  Segment parameter definitions (used to build the dynamic form)
# ============================================================================

SEG_PARAMS: Dict[str, List[Tuple[str, str, float]]] = {
    # kind    : [(param_key, label, default), ...]
    "line":        [("length",    "Length [µm]",     100.0)],
    "arc":         [("radius",    "Radius [µm]",     500.0),
                    ("angle",     "Angle [deg] (±)",  90.0)],
    "tapered_arc": [("radius",    "Radius [µm]",     500.0),
                    ("angle",     "Angle [deg] (±)",  90.0),
                    ("end_width", "End width [µm]",    3.0)],
    "sbend":       [("dx",        "dx [µm]",         300.0),
                    ("dy",        "dy [µm]  (±)",    100.0)],
    "rel_bent":    [("radius",    "Bend radius [µm]", 500.0),
                    ("dy",        "Target y [µm] (abs; e.g. yoffset(5))", 100.0)],
    "euler":       [("radius",    "End radius [µm]", 500.0),
                    ("angle",     "Angle [deg] (±)",  90.0)],
    "taper":       [("length",    "Length [µm]",     200.0),
                    ("end_width", "End width [µm]",    3.0)],
    "spiral":      [("radius_min","Min bend R [µm]", 100.0),
                    ("pitch",     "Spacing [µm]",     10.0),
                    ("n_turns",   "Turns / side",      3.0)],
    # Note: 'style' and 'straight_L' are supported too but not exposed
    # in the generic Add Segment form — use the Special Device menu.
}

SEG_DISPLAY = {
    "line":        "LINE          (straight)",
    "arc":         "ARC           (circular bend, uniform width)",
    "tapered_arc": "TAPERED ARC   (bend + linear width change)",
    "sbend":       "SBEND         (cosine S-bend)",
    "rel_bent":    "REL_BENT      (offset bend → absolute y target, ends at 0°)",
    "euler":       "EULER         (clothoid bend)",
    "taper":       "TAPER         (width transition, straight)",
    "spiral":      "SPIRAL        (Archimedean double spiral, compact delay line)",
}


DEFAULT_SCRIPT = """\
# Build the waveguide as a Python script.
# Available names : proj, Project, Layer, Port, Segment, group
# Global variables defined in the 'Globals' tab are injected as-is.
# Click  ▶  Redraw from script   (or press Ctrl+Enter) to apply.

# ----- (1) Globals — pull from the UI, with safe defaults ------------
CORE_W = globals().get('CORE_W', 2.0)   # core waveguide width  [µm]
BUS_L  = globals().get('BUS_L',  600.0) # bus straight length   [µm]
N_ARMS = int(globals().get('N_ARMS', 4))# number of parallel arms
PITCH  = globals().get('PITCH',  40.0)  # arm-to-arm pitch      [µm]

# ----- (2) Define a group (reusable sub-circuit) ---------------------
@group("dc_coupler")
def dc_coupler(start_pid, length=300.0, width=2.0, gap=3.0, layer=0):
    '''Simple directional coupler: two parallel straights with a gap.'''
    proj.add_segment("line", start_pid, layer,
                     {"length": length, "width": width,
                      "yoffset":  gap / 2.0})
    proj.add_segment("line", start_pid, layer,
                     {"length": length, "width": width,
                      "yoffset": -gap / 2.0})

@group("spiral_stub")
def spiral_stub(start_pid, loops=2, R=120.0, width=1.5, layer=0):
    '''Demo group — a snaking path that uses a local `for` loop.'''
    for i in range(loops):
        sign = 1 if i % 2 == 0 else -1
        proj.add_segment("arc", start_pid, layer,
                         {"radius": R, "angle": 180*sign, "width": width,
                          "yoffset": i * 2 * R})

# ----- (3) Main chain with a top-level for loop ----------------------
proj.add_base_port(0.0, 0.0, 0.0, CORE_W, 0, "in")                 # pid 1
proj.add_segment("line",  1, 0, {"length": 200.0, "width": CORE_W})
proj.add_segment("sbend", 2, 0, {"dx": 250.0, "dy": 80.0, "width": CORE_W})
proj.add_segment("taper", 3, 0, {"length": 100.0, "width": CORE_W,
                                 "end_width": 4.0})

# ----- (4) for-loop: N_ARMS parallel bus waveguides below pid 1 ------
for i in range(N_ARMS):
    proj.add_segment("line", 1, 0,
                     {"length": BUS_L, "width": 1.5,
                      "yoffset": -100 - i * PITCH})

# ----- (5) Call groups -----------------------------------------------
proj.add_base_port(0.0, -600.0, 0.0, CORE_W, 0, "dc_in")           # pid base
dc_in_pid = proj.ports[-1].pid
proj.call_group("dc_coupler", start_pid=dc_in_pid,
                length=350.0, gap=4.0)

# Use the spiral_stub group, shifted via a base port
proj.add_base_port(900.0, -600.0, 0.0, 1.5, 0, "spiral_in")
spiral_pid = proj.ports[-1].pid
proj.call_group("spiral_stub", start_pid=spiral_pid, loops=3, R=80.0)
"""


# ============================================================================
#  GUI
# ============================================================================

def _store_param(text, evaluated=None):
    """Decide whether a user-entered param is stored as a number or as
    a raw expression string. This is the key to parametric layouts:
    values entered as plain numbers ("2.0") are stored as floats, while
    anything symbolic ("R", "width*2", "sin(pi/4)") is kept as a string
    so it can be re-evaluated against the latest globals every time
    segment_geometry runs.
    """
    s = str(text).strip()
    try:
        return float(s)            # pure number literal → store as float
    except ValueError:
        return s                   # expression / variable name → keep raw


def _parse_num(text, globals_dict=None):
    """Parse a number from a UI text field.

    Accepts either a plain Python float literal ("2.0", "-1.5e-3") or a
    Python expression that may reference project globals and a small set
    of math helpers ("WIDTH*2", "PITCH+10", "sin(pi/4)").

    Globals whose stored *value* is itself a string are treated as
    expressions and resolved first — this lets a user define
    B = '(width+gap)/2'  in the Globals tab and then write  B  anywhere
    a numeric field is expected. Resolution iterates up to N rounds so
    chains like C = 'B*2', B = '(width+gap)/2' work regardless of the
    insertion order. Unresolvable strings are left as-is (they may or
    may not be needed by the caller's expression).
    """
    s = str(text).strip()
    if not s:
        raise ValueError("empty value")
    # Fast path: plain number
    try:
        return float(s)
    except ValueError:
        pass
    # Expression path
    import math
    safe_builtins = {
        "abs": abs, "min": min, "max": max, "round": round,
        "int": int, "float": float, "pow": pow,
    }
    env = {"__builtins__": safe_builtins}
    for name in ("sin", "cos", "tan", "asin", "acos", "atan", "atan2",
                 "sqrt", "exp", "log", "log10", "log2", "pi", "e",
                 "ceil", "floor", "degrees", "radians", "hypot"):
        env[name] = getattr(math, name, None)
    if globals_dict:
        # Pass 1 — numeric globals go in directly
        for k, v in globals_dict.items():
            if not isinstance(v, str):
                env[k] = v
        # Pass 2+ — resolve string globals as expressions, to a fixed
        # point (handles chains like C→B→width). At most len+2 rounds.
        pending = {k: v for k, v in globals_dict.items()
                   if isinstance(v, str)}
        for _ in range(len(pending) + 2):
            if not pending:
                break
            progressed = False
            for k, v in list(pending.items()):
                try:
                    env[k] = float(eval(v, dict(env)))
                    pending.pop(k)
                    progressed = True
                except Exception:
                    pass
            if not progressed:
                break
        # anything left unresolved stays as its raw string in env; if
        # the user's top-level expression doesn't reference it, no harm.
        for k, v in pending.items():
            env[k] = v
    try:
        result = eval(s, env)
    except Exception as e:
        raise ValueError(f"cannot evaluate {s!r}: {e}")
    try:
        return float(result)
    except Exception:
        raise ValueError(f"result {result!r} is not a number")


def _make_scrolled_tree(parent, columns, col_configs, height=10,
                        pack_in=True, padx=4, pady=4):
    """Wrap a Treeview in a Frame with a vertical scrollbar.

    Parameters
    ----------
    parent      : tk container
    columns     : tuple of column ids
    col_configs : list of (col_id, heading_text, pixel_width) triples
    height      : Treeview height in rows
    pack_in     : True → pack the wrapper; False → caller places it

    Returns
    -------
    (wrapper_frame, tree)
    """
    wrap = ttk.Frame(parent)
    if pack_in:
        wrap.pack(fill=tk.BOTH, expand=True, padx=padx, pady=pady)
    tree = ttk.Treeview(wrap, columns=columns, show="headings",
                        height=height)
    vsb = ttk.Scrollbar(wrap, orient=tk.VERTICAL, command=tree.yview)
    tree.configure(yscrollcommand=vsb.set)
    for c, t, w in col_configs:
        tree.heading(c, text=t)
        tree.column(c, width=w, anchor=tk.CENTER)
    tree.grid(row=0, column=0, sticky="nsew")
    vsb.grid(row=0, column=1, sticky="ns")
    wrap.rowconfigure(0, weight=1)
    wrap.columnconfigure(0, weight=1)
    return wrap, tree


class WaveguideApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("Waveguide Layout Designer — OlympiOS-style")
        root.geometry("1600x1000")

        self.project = Project()

        # middle-mouse-drag panning state (populated on button_press_event)
        self._pan_start: Optional[dict] = None

        # (rev25) Canvas pick/selection state.
        # When the user left-clicks on a segment the matching sid is
        # stored here and highlighted on redraw. Double-click opens the
        # SegmentDialog for that segment.
        self.selected_sid: Optional[int] = None
        # track a press position so a tiny drag doesn't count as a pick
        self._click_press: Optional[dict] = None

        # (rev27) Undo history for the Segment Sequence editor. Each entry is
        # a full project snapshot (Project.to_dict) captured *before* a
        # mutating action; "Undo" pops the most recent and restores it.
        self._undo_stack: List[Dict] = []
        self._UNDO_LIMIT = 50

        # (rev26) Canvas click behaviour mode.
        #   False (default): single-click selects, double-click edits
        #   True           : single-click selects AND opens editor
        # Toggled via the View menu "Edit on single click".
        self.click_to_edit: bool = False

        # three-pane layout
        paned = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)
        self._paned = paned

        self.left  = ttk.Frame(paned, width=340)
        self.mid   = ttk.Frame(paned)
        self.right = ttk.Frame(paned, width=480)
        # Give weights to all three so resizing distributes space rather
        # than letting the right pane get squeezed to zero width.
        paned.add(self.left,  weight=1)
        paned.add(self.mid,   weight=3)
        paned.add(self.right, weight=2)

        self._build_left_pane()
        self._build_center_pane()
        self._build_right_pane()
        self._build_menu()

        self.refresh_all()

        # Force initial sash positions — ttk.PanedWindow's default
        # allocation sometimes collapses the right pane below visibility
        # on first launch (especially on macOS where 'width=' hints are
        # advisory). Setting sashes after the geometry manager has run
        # guarantees both the left (340 px) and right (~480 px) panes
        # are visible.
        def _fix_sashes():
            try:
                total = paned.winfo_width()
                if total < 500:
                    # geometry not yet settled — retry shortly
                    root.after(40, _fix_sashes); return
                # Left sash fixed near 340 px; right sash ~480 px from
                # the right edge (i.e. total - 480).
                paned.sashpos(0, 340)
                paned.sashpos(1, max(total - 480, 360 + 400))
            except Exception:
                pass
        root.after(60, _fix_sashes)

    # ------------------------------------------------------------------
    #  Menu
    # ------------------------------------------------------------------
    def _build_menu(self):
        mbar = tk.Menu(self.root)
        fmenu = tk.Menu(mbar, tearoff=0)
        fmenu.add_command(label="New Project",    command=self.new_project)
        fmenu.add_command(label="Open Project…",  command=self.load_project)
        fmenu.add_command(label="Save Project…",  command=self.save_project)
        fmenu.add_separator()
        fmenu.add_command(label="Export GDS…",    command=self.export_gds)
        fmenu.add_separator()
        fmenu.add_command(label="Exit",           command=self.root.quit)
        mbar.add_cascade(label="File", menu=fmenu)

        vmenu = tk.Menu(mbar, tearoff=0)
        vmenu.add_command(label="Fit view",       command=self.fit_view)
        vmenu.add_command(label="Redraw",         command=self.redraw_canvas)
        vmenu.add_separator()
        # (rev26) Toggle: click-to-edit. When on, a single left-click on
        # a segment opens the editor dialog immediately. When off (the
        # default), single-click selects and double-click edits.
        self._click_to_edit_var = tk.BooleanVar(value=self.click_to_edit)
        vmenu.add_checkbutton(
            label="Edit on single click",
            variable=self._click_to_edit_var,
            command=self._toggle_click_to_edit,
        )
        mbar.add_cascade(label="View", menu=vmenu)

        # Special-device menu — one entry per compound/parametric
        # device. Each opens its own design dialog and writes a ready-
        # to-run snippet (usually an @group definition + call) into the
        # script editor.
        smenu = tk.Menu(mbar, tearoff=0)
        smenu.add_command(label="Spiral waveguide…",
                          command=self.add_spiral_device)
        mbar.add_cascade(label="Special device", menu=smenu)

        hmenu = tk.Menu(mbar, tearoff=0)
        hmenu.add_command(label="About", command=self._show_about)
        mbar.add_cascade(label="Help", menu=hmenu)

        self.root.config(menu=mbar)

    def _toggle_click_to_edit(self):
        """Flip click-to-edit mode from the View menu."""
        self.click_to_edit = bool(self._click_to_edit_var.get())
        # Update the canvas title hint so the user sees the change.
        self.redraw_canvas()

    def _show_about(self):
        messagebox.showinfo(
            "About",
            "Waveguide Layout Designer (OlympiOS-style) — rev27\n\n"
            "Build a photonic waveguide as a sequence of primitives\n"
            "(LINE / ARC / SBEND / REL_BENT / EULER / TAPER) starting from\n"
            "user-defined base ports. Each segment's endpoint is\n"
            "added to the port list so the next segment can chain\n"
            "off it automatically.\n\n"
            "Layout renders live; GDS export uses gdsfactory.")

    # ------------------------------------------------------------------
    #  Left pane — Notebook with Layers / Ports / Globals / Groups tabs
    # ------------------------------------------------------------------
    def _build_left_pane(self):
        nb = ttk.Notebook(self.left)
        nb.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        tab_layers  = ttk.Frame(nb); nb.add(tab_layers,  text="Layers")
        tab_ports   = ttk.Frame(nb); nb.add(tab_ports,   text="Ports")
        tab_globals = ttk.Frame(nb); nb.add(tab_globals, text="Globals")
        tab_groups  = ttk.Frame(nb); nb.add(tab_groups,  text="Groups")
        tab_offset  = ttk.Frame(nb); nb.add(tab_offset,  text="Auto-Offset")

        self._build_layers_tab(tab_layers)
        self._build_ports_tab(tab_ports)
        self._build_globals_tab(tab_globals)
        self._build_groups_tab(tab_groups)
        self._build_offset_tab(tab_offset)

    # --------- Tab 1: Layers ------------------------------------------
    def _build_layers_tab(self, parent):
        lbox = ttk.LabelFrame(parent, text="Layers")
        lbox.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        _, self.layer_tree = _make_scrolled_tree(
            lbox,
            columns=("name", "layer", "dt", "width"),
            col_configs=[("name", "Name", 90), ("layer", "GDS", 40),
                         ("dt", "DT", 30), ("width", "W [µm]", 60)],
            height=10,
        )
        self.layer_tree.bind("<Double-1>", lambda e: self.edit_layer())

        bf = ttk.Frame(lbox); bf.pack(fill=tk.X, padx=4, pady=(0, 4))
        ttk.Button(bf, text="Add",    command=self.add_layer).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Edit",   command=self.edit_layer).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Delete", command=self.delete_layer).pack(side=tk.LEFT, padx=2)

    # --------- Tab 2: Ports ------------------------------------------
    def _build_ports_tab(self, parent):
        pbox = ttk.LabelFrame(parent, text="Base Ports (anchors)")
        pbox.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        _, self.port_tree = _make_scrolled_tree(
            pbox,
            columns=("pid", "x", "y", "ang", "w", "ly", "lbl"),
            col_configs=[("pid","#",30),("x","x",55),("y","y",55),
                         ("ang","θ°",40),("w","w",40),("ly","Layer",50),
                         ("lbl","label",60)],
            height=10,
        )
        self.port_tree.bind("<Double-1>", lambda e: self.edit_port())

        bf = ttk.Frame(pbox); bf.pack(fill=tk.X, padx=4, pady=(0, 4))
        ttk.Button(bf, text="Add",    command=self.add_port).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Edit",   command=self.edit_port).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Delete", command=self.delete_port).pack(side=tk.LEFT, padx=2)

        abox = ttk.LabelFrame(parent, text="All Ports  (base + auto-generated)")
        abox.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
        _, self.all_port_tree = _make_scrolled_tree(
            abox,
            columns=("pid","x","y","ang","w","type"),
            col_configs=[("pid","#",30),("x","x",65),("y","y",65),
                         ("ang","θ°",45),("w","w",45),("type","type",55)],
            height=12,
        )

    # --------- Tab 3: Global variables -------------------------------
    def _build_globals_tab(self, parent):
        gbox = ttk.LabelFrame(parent, text="Global Variables")
        gbox.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        hint = ttk.Label(
            gbox, foreground="#666",
            text=("These names are auto-injected into the script namespace\n"
                  "before every run. Values accept Python literals:\n"
                  "  numbers, 'strings', [lists], (tuples), {dicts}, True/False/None"),
            justify=tk.LEFT, wraplength=300)
        hint.pack(fill=tk.X, padx=4, pady=(4, 4))

        _, self.glob_tree = _make_scrolled_tree(
            gbox,
            columns=("name", "value", "type"),
            col_configs=[("name","Name",90),("value","Value",160),
                         ("type","Type",60)],
            height=14,
        )
        # left-align rather than centred for readability
        for c in ("name", "value", "type"):
            self.glob_tree.column(c, anchor=tk.W)
        self.glob_tree.bind("<Double-1>", lambda e: self.edit_global())

        bf = ttk.Frame(gbox); bf.pack(fill=tk.X, padx=4, pady=(0, 4))
        ttk.Button(bf, text="Add",    command=self.add_global).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Edit",   command=self.edit_global).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Delete", command=self.delete_global).pack(side=tk.LEFT, padx=2)

    # --------- Tab 4: Groups (registered in script) ------------------
    def _build_groups_tab(self, parent):
        gbox = ttk.LabelFrame(parent, text="Defined Groups")
        gbox.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        hint = ttk.Label(
            gbox, foreground="#666",
            text=("Groups are Python callables registered with the\n"
                  "@group('name') decorator (or proj.register_group).\n"
                  "Invoke them with  proj.call_group('name', ...).\n"
                  "This list refreshes after every script run."),
            justify=tk.LEFT, wraplength=300)
        hint.pack(fill=tk.X, padx=4, pady=(4, 4))

        _, self.group_tree = _make_scrolled_tree(
            gbox,
            columns=("name", "signature"),
            col_configs=[("name", "Name", 110),
                         ("signature", "Signature", 200)],
            height=14,
        )
        for c in ("name", "signature"):
            self.group_tree.column(c, anchor=tk.W)
        self.group_tree.bind("<Double-1>",
                              lambda e: self.call_selected_group())

        bf = ttk.Frame(gbox); bf.pack(fill=tk.X, padx=4, pady=(0, 4))
        ttk.Button(bf, text="▶  Call selected group…",
                   command=self.call_selected_group
                   ).pack(side=tk.LEFT, padx=2)
        ttk.Label(bf, foreground="#666",
                  text="  (double-click도 가능)").pack(side=tk.LEFT, padx=2)

    # --------- Tab 5: Auto bend-offset  (rev25) ----------------------
    def _build_offset_tab(self, parent):
        """UI for PlatformConfig (enable + W + Δ + coefficient + presets)
        and a live table of the junction offsets that rebuild() injected.
        """
        # --- Platform box -----------------------------------------------
        pbox = ttk.LabelFrame(parent, text="Platform (for auto bend-offset)")
        pbox.pack(fill=tk.X, padx=6, pady=(6, 3))

        hint = ttk.Label(
            pbox, foreground="#666",
            text=("δ(R) = coeff · W² / (R · Δ)\n"
                  "bent waveguide의 mode center가 곡률중심 반대쪽으로 δ 만큼\n"
                  "어긋나는 현상을 straight↔bend, bend↔bend 접합에서 자동 보정\n"
                  "(+ : 진행방향 기준 왼쪽으로 shift)"),
            justify=tk.LEFT)
        hint.pack(fill=tk.X, padx=6, pady=(4, 4))

        # Enable checkbox
        self.off_enabled_var = tk.BooleanVar(value=self.project.platform.enabled)
        ttk.Checkbutton(pbox, text="Enable automatic bend-offset",
                        variable=self.off_enabled_var,
                        command=self._offset_apply_fields
                        ).pack(anchor="w", padx=6, pady=(0, 4))

        # Numeric entries
        ef = ttk.Frame(pbox); ef.pack(fill=tk.X, padx=6, pady=2)
        self.off_w_var = tk.StringVar(value=f"{self.project.platform.width:g}")
        self.off_d_var = tk.StringVar(value=f"{self.project.platform.delta:g}")
        self.off_c_var = tk.StringVar(value=f"{self.project.platform.offset_coeff:g}")

        def _row(r, label, var, col=0):
            ttk.Label(ef, text=label).grid(row=r, column=col*2,
                                            sticky="e", padx=(0, 4), pady=2)
            e = ttk.Entry(ef, textvariable=var, width=10)
            e.grid(row=r, column=col*2 + 1, sticky="w", padx=(0, 8), pady=2)
            e.bind("<Return>",   lambda _e: self._offset_apply_fields())
            e.bind("<FocusOut>", lambda _e: self._offset_apply_fields())
            return e

        _row(0, "W  [μm]:",  self.off_w_var)
        _row(1, "Δ  (rel):", self.off_d_var)
        _row(2, "coeff:",    self.off_c_var)

        # Preset buttons
        prebox = ttk.Frame(pbox); prebox.pack(fill=tk.X, padx=6, pady=(6, 4))
        ttk.Label(prebox, text="Presets:").pack(side=tk.LEFT, padx=(0, 4))
        for label, (w, d) in [
            ("Silica 0.45%", (8.0,  0.0045)),
            ("Silica 0.75%", (6.0,  0.0075)),
            ("Silica 1.5%",  (4.5,  0.015)),
            ("Silica 2.0%",  (4.0,  0.020)),
            ("Ti:LiNbO₃",    (7.0,  0.003)),
            ("LNOI",         (1.2,  0.35)),
        ]:
            ttk.Button(prebox, text=label, width=12,
                       command=lambda w=w, d=d: self._offset_load_preset(w, d)
                       ).pack(side=tk.LEFT, padx=1)

        # --- Current δ(R) preview ---------------------------------------
        self.off_preview_var = tk.StringVar(value="δ(R=1000)=…  δ(R=5000)=…")
        ttk.Label(pbox, textvariable=self.off_preview_var,
                  foreground="#064").pack(anchor="w", padx=6, pady=(0, 6))

        # --- Junction log table -----------------------------------------
        lbox = ttk.LabelFrame(parent,
                              text="Junction Offsets  (populated after rebuild)")
        lbox.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        _, self.offset_tree = _make_scrolled_tree(
            lbox,
            columns=("sid", "port", "prev", "curr", "off"),
            col_configs=[("sid",  "S#",       40),
                         ("port", "@port",    50),
                         ("prev", "prev",     90),
                         ("curr", "curr",     90),
                         ("off",  "δ [nm]",   80)],
            height=14,
        )
        for c in ("sid", "port", "prev", "curr"):
            self.offset_tree.column(c, anchor=tk.W)
        self.offset_tree.column("off", anchor=tk.E)

    # ---- helpers for the Auto-Offset tab --------------------------------
    def _offset_load_preset(self, w: float, delta: float):
        """Apply a platform preset (W, Δ) from the preset buttons."""
        self.off_w_var.set(f"{w:g}")
        self.off_d_var.set(f"{delta:g}")
        self._offset_apply_fields()

    def _offset_apply_fields(self):
        """Read the entry fields → push into proj.platform → rebuild + redraw."""
        plat = self.project.platform
        plat.enabled = bool(self.off_enabled_var.get())
        def _f(var, fallback):
            try:
                return float(var.get())
            except (tk.TclError, ValueError):
                return fallback
        plat.width        = _f(self.off_w_var, plat.width)
        plat.delta        = _f(self.off_d_var, plat.delta)
        plat.offset_coeff = _f(self.off_c_var, plat.offset_coeff)

        # preview δ at typical radii
        try:
            d1 = plat.bend_offset(1000.0) * 1000   # nm
            d5 = plat.bend_offset(5000.0) * 1000
            self.off_preview_var.set(
                f"δ(R=1000 μm) = {d1:6.1f} nm   "
                f"δ(R=5000 μm) = {d5:6.1f} nm")
        except Exception:
            pass

        # regenerate geometry + table
        self.project.rebuild()
        self._refresh_offset_table()
        self.redraw_canvas()
        # other tabs may show ports that moved
        self._refresh_port_lists()

    def _refresh_offset_table(self):
        """Push proj.offset_log into the tree widget + sync entry fields
        from proj.platform (needed after project load / script rebuild)."""
        tree = getattr(self, "offset_tree", None)
        if tree is None:
            return
        # sync entry fields with current proj.platform (tab might have been
        # built before project was loaded, or a script rewrote platform)
        plat = self.project.platform
        if hasattr(self, "off_enabled_var"):
            self.off_enabled_var.set(plat.enabled)
            self.off_w_var.set(f"{plat.width:g}")
            self.off_d_var.set(f"{plat.delta:g}")
            self.off_c_var.set(f"{plat.offset_coeff:g}")
            try:
                d1 = plat.bend_offset(1000.0) * 1000
                d5 = plat.bend_offset(5000.0) * 1000
                self.off_preview_var.set(
                    f"δ(R=1000 μm) = {d1:6.1f} nm   "
                    f"δ(R=5000 μm) = {d5:6.1f} nm")
            except Exception:
                pass
        tree.delete(*tree.get_children())
        for rec in self.project.offset_log:
            tree.insert("", "end", values=(
                f"S{rec['sid']}",
                rec["start_pid"],
                (f"S{rec['prev_sid']} ({rec['prev_kind']})"
                 if rec["prev_sid"] else "—"),
                rec["curr_kind"],
                f"{rec['offset_um']*1000:+.1f}",
            ))

    # ------------------------------------------------------------------
    #  Center pane — matplotlib canvas (top) + script editor (bottom)
    # ------------------------------------------------------------------
    def _build_center_pane(self):
        # Vertical split inside the mid pane
        self.mid_split = ttk.PanedWindow(self.mid, orient=tk.VERTICAL)
        self.mid_split.pack(fill=tk.BOTH, expand=True)

        canvas_frame = ttk.Frame(self.mid_split)
        script_frame = ttk.Frame(self.mid_split)
        # canvas : script  ≈  1 : 1 — drag the horizontal sash between
        # them to resize; the script editor starts roughly half the
        # vertical space of the centre column.
        self.mid_split.add(canvas_frame, weight=1)
        self.mid_split.add(script_frame, weight=1)

        # --- Canvas ---
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_aspect("equal", adjustable="datalim")
        self.ax.grid(True, alpha=0.3)
        self.ax.set_xlabel("x [µm]")
        self.ax.set_ylabel("y [µm]")
        self.ax.set_title("Waveguide Layout")

        self.canvas = FigureCanvasTkAgg(self.fig, master=canvas_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, canvas_frame)
        self.toolbar.update()

        # Make the canvas widget itself focusable so Tk-level clicks
        # (used as a reliable fallback on macOS where matplotlib's
        # button_press_event is eaten by the toolbar Pan/Zoom modes)
        # can be captured directly.
        try:
            self.canvas.get_tk_widget().configure(takefocus=1)
        except tk.TclError:
            pass

        # Mouse-wheel zoom (cursor-centred)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)

        # Mouse button handlers — both left-click pick and middle/right pan
        self.canvas.mpl_connect("button_press_event",   self._on_mouse_press)
        self.canvas.mpl_connect("motion_notify_event",  self._on_mouse_motion)
        self.canvas.mpl_connect("button_release_event", self._on_mouse_release)

        # (rev25) Tk-level fallback for left-click segment pick.
        # On macOS the matplotlib toolbar's pan/zoom modes often consume
        # button_press_event so the mpl handler never fires. We bind the
        # raw Tk events directly to the canvas widget as a safety net —
        # these fire regardless of toolbar mode, and we auto-disable
        # pan/zoom on click so normal picking feels natural.
        wid = self.canvas.get_tk_widget()
        wid.bind("<Button-1>",         self._on_tk_click, add="+")
        wid.bind("<Double-Button-1>",  self._on_tk_double_click, add="+")

        # --- Script editor ---
        hdr = ttk.Frame(script_frame)
        hdr.pack(fill=tk.X, padx=4, pady=(4, 0))
        ttk.Label(hdr, text="Script",
                  font=("TkDefaultFont", 9, "bold")).pack(side=tk.LEFT)
        ttk.Button(hdr, text="▶  Redraw from script  (Ctrl+Enter)",
                   command=self.run_script).pack(side=tk.LEFT, padx=8)
        ttk.Button(hdr, text="⤓  Dump current → script",
                   command=self.dump_to_script).pack(side=tk.LEFT, padx=4)

        # Explicit edit buttons — work on macOS even when Cmd keys
        # don't get forwarded to the Text widget.
        ttk.Separator(hdr, orient=tk.VERTICAL).pack(side=tk.LEFT,
                                                   fill=tk.Y, padx=6)
        for lbl, cmd in [("Sel all", self._script_select_all),
                         ("Copy",    self._script_copy),
                         ("Cut",     self._script_cut),
                         ("Paste",   self._script_paste)]:
            ttk.Button(hdr, text=lbl, width=7,
                       command=cmd).pack(side=tk.LEFT, padx=1)

        self.script_status = ttk.Label(hdr, text="", foreground="#2a8a2a")
        self.script_status.pack(side=tk.LEFT, padx=10)

        # Right-click context menu for the script editor
        self._script_menu = tk.Menu(self.root, tearoff=0)
        self._script_menu.add_command(label="Copy",        command=self._script_copy)
        self._script_menu.add_command(label="Cut",         command=self._script_cut)
        self._script_menu.add_command(label="Paste",       command=self._script_paste)
        self._script_menu.add_separator()
        self._script_menu.add_command(label="Select all",  command=self._script_select_all)

        # scrolled text area — dark theme so text stays readable under
        # any system light/dark theme (macOS dark mode otherwise inherits
        # a white foreground that disappears on the light background).
        txt_frame = ttk.Frame(script_frame)
        txt_frame.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.script_text = tk.Text(
            txt_frame,
            height=24,                         # ↑ from 18 (starts much taller)
            wrap=tk.NONE,
            font=("TkFixedFont", 11),          # slightly bigger
            undo=True,
            background="#1e1e1e",              # VS Code-style dark bg
            foreground="#e6e6e6",              # bright text
            insertbackground="#ffffff",        # bright cursor
            selectbackground="#264f78",
            selectforeground="#ffffff",
            borderwidth=0,
            highlightthickness=1,
            highlightbackground="#3c3c3c",
            highlightcolor="#5c5c5c",
        )
        sy = ttk.Scrollbar(txt_frame, orient=tk.VERTICAL,
                           command=self.script_text.yview)
        sx = ttk.Scrollbar(txt_frame, orient=tk.HORIZONTAL,
                           command=self.script_text.xview)
        self.script_text.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        self.script_text.grid(row=0, column=0, sticky="nsew")
        sy.grid(row=0, column=1, sticky="ns")
        sx.grid(row=1, column=0, sticky="ew")
        txt_frame.rowconfigure(0, weight=1)
        txt_frame.columnconfigure(0, weight=1)

        # Ctrl+Enter / Cmd+Enter runs the script
        self.script_text.bind("<Control-Return>", lambda e: (self.run_script(),
                                                             "break"))
        self.script_text.bind("<Command-Return>", lambda e: (self.run_script(),
                                                             "break"))

        # Keyboard shortcuts — use our direct-clipboard handlers
        # (event_generate("<<Paste>>") doesn't work on macOS Tkinter).
        for key, handler in [("c", lambda e: self._script_copy()),
                             ("v", lambda e: self._script_paste()),
                             ("x", lambda e: self._script_cut()),
                             ("a", lambda e: self._script_select_all())]:
            self.script_text.bind(f"<Command-{key}>", handler)
            self.script_text.bind(f"<Command-{key.upper()}>", handler)
            self.script_text.bind(f"<Control-{key}>", handler)
            self.script_text.bind(f"<Control-{key.upper()}>", handler)

        def _do_undo(e):
            try: self.script_text.edit_undo()
            except Exception: pass
            return "break"
        def _do_redo(e):
            try: self.script_text.edit_redo()
            except Exception: pass
            return "break"
        self.script_text.bind("<Command-z>", _do_undo)
        self.script_text.bind("<Command-Z>", _do_undo)
        self.script_text.bind("<Control-z>", _do_undo)
        self.script_text.bind("<Control-Z>", _do_undo)
        self.script_text.bind("<Command-Shift-Z>", _do_redo)
        self.script_text.bind("<Command-Shift-z>", _do_redo)
        self.script_text.bind("<Control-Shift-Z>", _do_redo)
        self.script_text.bind("<Control-y>", _do_redo)
        self.script_text.bind("<Control-Y>", _do_redo)

        # Right-click popup menu
        self.script_text.bind("<Button-3>", self._script_popup)      # win/linux
        self.script_text.bind("<Button-2>", self._script_popup)      # macOS
        self.script_text.bind("<Control-Button-1>", self._script_popup)  # mac ctrl+click

        self.script_text.insert("1.0", DEFAULT_SCRIPT)

    # ------------------------------------------------------------------
    #  Right pane — Segment sequence editor
    # ------------------------------------------------------------------
    def _build_right_pane(self):
        # --- sequence list ---
        sbox = ttk.LabelFrame(self.right, text="Segment Sequence")
        sbox.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        _, self.seg_tree = _make_scrolled_tree(
            sbox,
            columns=("sid","kind","start","end","len","w","layer","params"),
            col_configs=[("sid","#",30),("kind","kind",55),
                         ("start","start→",50),("end","end",45),
                         ("len","len [µm]",60),
                         ("w","w [µm]",50),
                         ("layer","layer",55),("params","params",120)],
            height=14,
        )
        self.seg_tree.configure(selectmode="extended")   # multi-select
        self.seg_tree.bind("<Double-1>", lambda e: self.edit_segment())
        # (rev25) tree → canvas: mirror current selection as canvas highlight
        self.seg_tree.bind("<<TreeviewSelect>>",
                           lambda e: self._on_seg_tree_select())

        bf = ttk.Frame(sbox); bf.pack(fill=tk.X, padx=4, pady=(0, 4))
        ttk.Button(bf, text="↶ Undo", command=self.undo).pack(side=tk.LEFT, padx=(1, 6))
        ttk.Button(bf, text="↑",      command=lambda: self.move_segment(-1)).pack(side=tk.LEFT, padx=1)
        ttk.Button(bf, text="↓",      command=lambda: self.move_segment(+1)).pack(side=tk.LEFT, padx=1)
        ttk.Button(bf, text="Edit",   command=self.edit_segment).pack(side=tk.LEFT, padx=4)
        ttk.Button(bf, text="Delete", command=self.delete_segment).pack(side=tk.LEFT, padx=2)
        ttk.Separator(bf, orient=tk.VERTICAL).pack(side=tk.LEFT,
                                                   fill=tk.Y, padx=6)
        ttk.Button(bf, text="Select all",
                   command=self.select_all_segments).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Delete all",
                   command=self.delete_all_segments).pack(side=tk.LEFT, padx=2)
        ttk.Separator(bf, orient=tk.VERTICAL).pack(side=tk.LEFT,
                                                   fill=tk.Y, padx=6)
        ttk.Button(bf, text="📏 Measure",
                   command=self.measure_segments).pack(side=tk.LEFT, padx=2)
        ttk.Separator(bf, orient=tk.VERTICAL).pack(side=tk.LEFT,
                                                   fill=tk.Y, padx=6)
        ttk.Button(bf, text="⚙  Group selected…",
                   command=self.group_selected_segments
                   ).pack(side=tk.LEFT, padx=2)

        # Keyboard: Ctrl/Cmd+Z undo, Ctrl/Cmd+A select-all (when tree focused)
        self.seg_tree.bind("<Control-z>", lambda e: (self.undo(), "break")[1])
        self.seg_tree.bind("<Command-z>", lambda e: (self.undo(), "break")[1])
        self.seg_tree.bind("<Control-a>",
                           lambda e: (self.select_all_segments(), "break")[1])
        self.seg_tree.bind("<Command-a>",
                           lambda e: (self.select_all_segments(), "break")[1])

        # --- add segment form ---
        fbox = ttk.LabelFrame(self.right, text="Add Segment")
        fbox.pack(fill=tk.X, padx=6, pady=6)

        row = ttk.Frame(fbox); row.pack(fill=tk.X, padx=4, pady=3)
        ttk.Label(row, text="Type:").pack(side=tk.LEFT)
        self.kind_var = tk.StringVar(value="line")
        kind_cb = ttk.Combobox(
            row, textvariable=self.kind_var, state="readonly", width=28,
            values=[f"{k}: {SEG_DISPLAY[k]}" for k in SEG_PARAMS])
        kind_cb.current(0)
        kind_cb.pack(side=tk.LEFT, padx=4)
        kind_cb.bind("<<ComboboxSelected>>", lambda e: self._rebuild_param_form())

        row = ttk.Frame(fbox); row.pack(fill=tk.X, padx=4, pady=3)
        ttk.Label(row, text="Start port #:").pack(side=tk.LEFT)
        self.start_pid_var = tk.StringVar(value="1")
        ttk.Entry(row, textvariable=self.start_pid_var, width=6
                  ).pack(side=tk.LEFT, padx=4)
        ttk.Label(row, text="End port #:").pack(side=tk.LEFT, padx=(8, 0))
        self.end_pid_var = tk.StringVar(value="")
        ttk.Entry(row, textvariable=self.end_pid_var, width=6
                  ).pack(side=tk.LEFT, padx=4)
        ttk.Label(row, text="(blank = auto)",
                  foreground="#666").pack(side=tk.LEFT)

        row = ttk.Frame(fbox); row.pack(fill=tk.X, padx=4, pady=3)
        ttk.Label(row, text="Layer:").pack(side=tk.LEFT)
        self.layer_var = tk.StringVar()
        self.layer_cb = ttk.Combobox(row, textvariable=self.layer_var,
                                     state="readonly", width=22)
        self.layer_cb.pack(side=tk.LEFT, padx=4)
        self.layer_cb.bind("<<ComboboxSelected>>",
                           lambda e: self._on_layer_change())

        # Width is a common field (applies to every segment type).
        row = ttk.Frame(fbox); row.pack(fill=tk.X, padx=4, pady=3)
        ttk.Label(row, text="Width [µm]:").pack(side=tk.LEFT)
        self.width_var = tk.StringVar(value="2.0")
        ttk.Entry(row, textvariable=self.width_var, width=10
                  ).pack(side=tk.LEFT, padx=4)
        ttk.Label(row, text="(taper: start width)",
                  foreground="#666").pack(side=tk.LEFT)

        # Offset fields — displace the segment's *effective start* from
        # the start port's position, in the port's local frame.
        row = ttk.Frame(fbox); row.pack(fill=tk.X, padx=4, pady=3)
        ttk.Label(row, text="Offset  x:").pack(side=tk.LEFT)
        self.xoff_var = tk.StringVar(value="0.0")
        ttk.Entry(row, textvariable=self.xoff_var, width=7
                  ).pack(side=tk.LEFT, padx=(2, 6))
        ttk.Label(row, text="y:").pack(side=tk.LEFT)
        self.yoff_var = tk.StringVar(value="0.0")
        ttk.Entry(row, textvariable=self.yoff_var, width=7
                  ).pack(side=tk.LEFT, padx=(2, 6))
        ttk.Label(row, text="θ°:").pack(side=tk.LEFT)
        self.aoff_var = tk.StringVar(value="0.0")
        ttk.Entry(row, textvariable=self.aoff_var, width=6
                  ).pack(side=tk.LEFT, padx=(2, 6))
        ttk.Label(row, text="(local frame)",
                  foreground="#666").pack(side=tk.LEFT)

        self.param_frame = ttk.Frame(fbox)
        self.param_frame.pack(fill=tk.X, padx=4, pady=4)
        self.param_vars: Dict[str, tk.StringVar] = {}
        self._rebuild_param_form()

        bf = ttk.Frame(fbox); bf.pack(fill=tk.X, padx=4, pady=4)
        ttk.Button(bf, text="＋ Append segment",
                   command=self.add_segment).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="↻ Redraw",
                   command=self.redraw_canvas).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Export GDS…",
                   command=self.export_gds).pack(side=tk.RIGHT, padx=2)

    # --- dynamic parameter form --------------------------------------
    def _selected_kind(self) -> str:
        raw = self.kind_var.get()
        # stored as "line: LINE   (straight)"
        return raw.split(":")[0].strip() if ":" in raw else raw

    def _rebuild_param_form(self):
        for ch in self.param_frame.winfo_children():
            ch.destroy()
        self.param_vars.clear()
        kind = self._selected_kind()
        for key, label, default in SEG_PARAMS[kind]:
            row = ttk.Frame(self.param_frame)
            row.pack(fill=tk.X, pady=2)
            ttk.Label(row, text=label, width=18, anchor=tk.W
                      ).pack(side=tk.LEFT)
            v = tk.StringVar(value=str(default))
            ttk.Entry(row, textvariable=v, width=12).pack(side=tk.LEFT)
            self.param_vars[key] = v

    # ------------------------------------------------------------------
    #  Refresh
    # ------------------------------------------------------------------
    def refresh_all(self):
        self._refresh_layer_list()
        self._refresh_port_lists()
        self._refresh_segment_list()
        self._refresh_layer_combobox()
        self._refresh_globals_list()
        self._refresh_groups_list()
        self._refresh_offset_table()        # (rev25)
        self.redraw_canvas()

    def _refresh_globals_list(self):
        if not hasattr(self, "glob_tree"): return
        self.glob_tree.delete(*self.glob_tree.get_children())
        for name, value in sorted(self.project.globals.items()):
            self.glob_tree.insert(
                "", tk.END, iid=name,
                values=(name, repr(value), type(value).__name__))

    def _refresh_groups_list(self):
        if not hasattr(self, "group_tree"): return
        self.group_tree.delete(*self.group_tree.get_children())
        for name, sig in self.project.list_groups():
            self.group_tree.insert("", tk.END, iid=name, values=(name, sig))

    def _refresh_layer_list(self):
        self.layer_tree.delete(*self.layer_tree.get_children())
        for i, L in enumerate(self.project.layers):
            self.layer_tree.insert("", tk.END, iid=str(i),
                                   values=(L.name, L.number, L.datatype, L.width))
            self.layer_tree.item(str(i), tags=(f"color_{i}",))
            self.layer_tree.tag_configure(f"color_{i}", background=L.color)

    def _on_layer_change(self):
        """When the user picks a new layer, pre-fill the width entry with
        that layer's default width (but only if the user hasn't typed
        something custom — we simply overwrite for simplicity)."""
        idx = self._selected_layer_idx()
        if 0 <= idx < len(self.project.layers):
            self.width_var.set(f"{self.project.layers[idx].width:g}")

    def _refresh_layer_combobox(self):
        vals = [f"{i}: {L.name}  ({L.number}/{L.datatype})"
                for i, L in enumerate(self.project.layers)]
        self.layer_cb["values"] = vals
        if vals and self.layer_var.get() not in vals:
            self.layer_cb.current(0)
        # keep the width entry in sync on first load
        if hasattr(self, "width_var"):
            idx = self._selected_layer_idx()
            if 0 <= idx < len(self.project.layers):
                # only overwrite if the field is still at its default
                cur = self.width_var.get().strip()
                if cur in ("", "2.0", "0", "0.0"):
                    self.width_var.set(f"{self.project.layers[idx].width:g}")

    def _refresh_port_lists(self):
        # base ports only in editable view
        self.port_tree.delete(*self.port_tree.get_children())
        for p in self.project.ports:
            if not p.is_base: continue
            layer = self.project.layers[p.layer_idx].name \
                    if 0 <= p.layer_idx < len(self.project.layers) else "?"
            self.port_tree.insert(
                "", tk.END, iid=str(p.pid),
                values=(p.pid, f"{p.x:.2f}", f"{p.y:.2f}",
                        f"{p.angle:.2f}", f"{p.width:.2f}", layer, p.label))
        # all ports including auto-generated
        self.all_port_tree.delete(*self.all_port_tree.get_children())
        for p in sorted(self.project.ports, key=lambda q: q.pid):
            self.all_port_tree.insert(
                "", tk.END, iid=f"a_{p.pid}",
                values=(p.pid, f"{p.x:.2f}", f"{p.y:.2f}",
                        f"{p.angle:.2f}", f"{p.width:.2f}",
                        "base" if p.is_base else "auto"))

    def _refresh_segment_list(self):
        self.seg_tree.delete(*self.seg_tree.get_children())
        prev_map = self.project._prev_by_start()
        for s in self.project.segments:
            ly = self.project.layers[s.layer_idx].name \
                 if 0 <= s.layer_idx < len(self.project.layers) else "?"
            # pull width out of params for its own column
            w_val = s.params.get("width", 0.0)
            try:
                w_str = f"{float(w_val):g}" if float(w_val) > 0 else "(port)"
            except Exception:
                w_str = "(port)"
            L = self._segment_length(s, prev_map)
            len_str = f"{L:.2f}" if L is not None else "—"
            pstr = ", ".join(f"{k}={v}" for k, v in s.params.items()
                             if k != "width")
            self.seg_tree.insert(
                "", tk.END, iid=str(s.sid),
                values=(s.sid, s.kind, s.start_pid, s.end_pid,
                        len_str, w_str, ly, pstr))

    # ------------------------------------------------------------------
    #  (rev27) Geometry helpers — segment length & spacing measurement
    # ------------------------------------------------------------------
    def _segment_geom(self, seg, prev_map=None):
        """Return (path_xy, widths, extra_polys) for a segment, or None."""
        start = self.project.get_port(seg.start_pid)
        if start is None:
            return None
        if prev_map is None:
            prev_map = self.project._prev_by_start()
        try:
            _res = segment_geometry(
                seg.kind, start,
                self.project.effective_params(seg, prev_map),
                globals_dict=self.project.eval_globals())
        except Exception:
            return None
        path = _res[0]
        widths = _res[1]
        extra = _res[3] if len(_res) > 3 else None
        return path, widths, extra

    def _segment_length(self, seg, prev_map=None):
        """Centreline arc length [µm] (polygon-perimeter/2 for spirals)."""
        g = self._segment_geom(seg, prev_map)
        if g is None:
            return None
        path, _w, extra = g
        if extra is not None:           # gdsfactory spiral → use polygons
            try:
                return _polys_length(extra)
            except Exception:
                return None
        if path is None or len(path) < 2:
            return 0.0
        return float(np.sum(np.hypot(np.diff(path[:, 0]),
                                     np.diff(path[:, 1]))))

    def measure_segments(self):
        sel = self.seg_tree.selection()
        if not sel:
            messagebox.showinfo(
                "Measure",
                "측정할 세그먼트를 표에서 선택하세요.\n"
                "  • 1개 선택 → 길이\n"
                "  • 2개 선택 → 각 길이 + 둘 사이 간격(gap)\n"
                "  • 3개 이상 → 각 길이 + 합계")
            return
        prev_map = self.project._prev_by_start()
        segs = [s for s in self.project.segments if str(s.sid) in set(sel)]
        if not segs:
            return

        lines = []
        total = 0.0
        for s in segs:
            L = self._segment_length(s, prev_map)
            if L is None:
                lines.append(f"S{s.sid} ({s.kind}): 길이 계산 불가")
            else:
                total += L
                lines.append(f"S{s.sid} ({s.kind}): {L:.3f} µm")
        lines.append(f"\n합계 길이: {total:.3f} µm")

        # exactly two → also report the gap between them
        if len(segs) == 2:
            gap = self._segment_gap(segs[0], segs[1], prev_map)
            if gap is None:
                lines.append("\n간격: 계산 불가")
            else:
                c_gap, e_gap, wa, wb = gap
                lines.append(
                    f"\n— S{segs[0].sid} ↔ S{segs[1].sid} 간격 —\n"
                    f"중심선 간 최소거리: {c_gap:.3f} µm\n"
                    f"가장자리 간 최소거리(gap): {e_gap:.3f} µm\n"
                    f"  (폭 wA≈{wa:.3g}, wB≈{wb:.3g} 적용)")
                if e_gap < 0:
                    lines.append("  ⚠ 두 도파로가 겹칩니다 (gap < 0).")

        messagebox.showinfo("Measure", "\n".join(lines))

    def _segment_gap(self, sa, sb, prev_map=None):
        """Min centre/edge distance between two segment centrelines.

        Returns (center_gap, edge_gap, w_a_at_min, w_b_at_min) or None.
        Point-to-point over the sampled centrelines — adequate for a
        layout readout.
        """
        ga = self._segment_geom(sa, prev_map)
        gb = self._segment_geom(sb, prev_map)
        if ga is None or gb is None:
            return None
        pa, wa, _ = ga
        pb, wb, _ = gb
        if pa is None or pb is None or len(pa) == 0 or len(pb) == 0:
            return None
        A = np.asarray(pa, dtype=float)
        B = np.asarray(pb, dtype=float)
        d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(axis=-1)
        ia, ib = np.unravel_index(int(np.argmin(d2)), d2.shape)
        c_gap = float(np.sqrt(d2[ia, ib]))
        try:
            wa_i = float(wa[ia]); wb_i = float(wb[ib])
        except Exception:
            wa_i = wb_i = 0.0
        e_gap = c_gap - wa_i / 2.0 - wb_i / 2.0
        return c_gap, e_gap, wa_i, wb_i

    # ------------------------------------------------------------------
    #  Canvas redraw
    # ------------------------------------------------------------------
    def redraw_canvas(self):
        self.ax.clear()
        self.ax.set_aspect("equal", adjustable="datalim")
        self.ax.grid(True, alpha=0.3)
        self.ax.set_xlabel("x [µm]")
        self.ax.set_ylabel("y [µm]")
        _mode_hint = ("click to edit"
                      if self.click_to_edit
                      else "click=select · double-click=edit")
        if self.selected_sid is not None:
            self.ax.set_title(f"Waveguide Layout   —   selected: S{self.selected_sid}   "
                              f"({_mode_hint} · middle/right-drag to pan)")
        else:
            self.ax.set_title(f"Waveguide Layout   —   "
                              f"{_mode_hint} · middle/right-drag to pan")

        # segments
        prev_map = self.project._prev_by_start()
        for seg in self.project.segments:
            start = self.project.get_port(seg.start_pid)
            if start is None: continue
            try:
                _res = segment_geometry(
                    seg.kind, start,
                    self.project.effective_params(seg, prev_map),
                    globals_dict=self.project.eval_globals())
                path   = _res[0]; widths = _res[1]
                extra_polys = _res[3] if len(_res) > 3 else None
            except Exception as e:
                print(f"[seg {seg.sid}] geometry error: {e}")
                continue
            color = self.project.layers[seg.layer_idx].color \
                    if 0 <= seg.layer_idx < len(self.project.layers) else "#888"

            # (rev25) highlight if this segment is the one picked on the canvas
            is_sel = (seg.sid == self.selected_sid)
            edge_c = "#e67e22" if is_sel else "black"   # orange when selected
            lw     = 1.8       if is_sel else 0.4

            if extra_polys is not None:
                # Multi-polygon geometry (gdsfactory-backed spiral)
                for P in extra_polys:
                    self.ax.fill(P[:, 0], P[:, 1],
                                 facecolor=color, edgecolor=edge_c,
                                 linewidth=lw, alpha=0.85)
                # Use centerline's midpoint as label anchor
                if len(path) >= 1:
                    mid = path[len(path)//2]
                    self.ax.annotate(f"S{seg.sid}", xy=mid,
                                     xytext=(4, 4), textcoords="offset points",
                                     fontsize=7,
                                     color="#e67e22" if is_sel else "black",
                                     weight="bold" if is_sel else "normal")
            else:
                poly = path_to_polygon(path, widths)
                self.ax.fill(poly[:, 0], poly[:, 1],
                             facecolor=color, edgecolor=edge_c,
                             linewidth=lw, alpha=0.85)
                # segment number at midpoint
                mid = path[len(path)//2]
                self.ax.annotate(f"S{seg.sid}", xy=mid,
                                 xytext=(4, 4), textcoords="offset points",
                                 fontsize=7,
                                 color="#e67e22" if is_sel else "black",
                                 weight="bold" if is_sel else "normal")

        # ports
        for p in self.project.ports:
            if p.is_base:
                self.ax.plot(p.x, p.y, marker="s", markersize=8,
                             markerfacecolor="red",
                             markeredgecolor="black", linestyle="None")
            else:
                self.ax.plot(p.x, p.y, marker="o", markersize=5,
                             markerfacecolor="yellow",
                             markeredgecolor="black", linestyle="None")
            # heading arrow
            dx, dy = _rot(p.angle)
            L = 8
            self.ax.annotate(
                "", xy=(p.x + L*dx, p.y + L*dy), xytext=(p.x, p.y),
                arrowprops=dict(arrowstyle="->", color="gray", lw=0.8))
            self.ax.annotate(
                f"{p.pid}", xy=(p.x, p.y),
                xytext=(6, -8), textcoords="offset points",
                fontsize=8, color="darkblue", weight="bold")

        # autoscale if empty
        if not self.project.segments and not self.project.ports:
            self.ax.set_xlim(-100, 100); self.ax.set_ylim(-100, 100)

        self.canvas.draw_idle()

    def fit_view(self):
        if not self.project.ports and not self.project.segments:
            return
        xs, ys = [], []
        for p in self.project.ports:
            xs.append(p.x); ys.append(p.y)
        prev_map = self.project._prev_by_start()
        for seg in self.project.segments:
            start = self.project.get_port(seg.start_pid)
            if start is None: continue
            try:
                _res = segment_geometry(
                    seg.kind, start,
                    self.project.effective_params(seg, prev_map),
                    globals_dict=self.project.eval_globals())
                path = _res[0]
                extra_polys = _res[3] if len(_res) > 3 else None
            except Exception:
                continue
            if extra_polys:
                for P in extra_polys:
                    xs.extend(P[:, 0].tolist()); ys.extend(P[:, 1].tolist())
            else:
                xs.extend(path[:, 0].tolist()); ys.extend(path[:, 1].tolist())
        if not xs: return
        pad = max(20.0, 0.05 * max(max(xs)-min(xs), max(ys)-min(ys)))
        self.ax.set_xlim(min(xs)-pad, max(xs)+pad)
        self.ax.set_ylim(min(ys)-pad, max(ys)+pad)
        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    #  Mouse-wheel zoom (cursor-centred)
    # ------------------------------------------------------------------
    def _on_scroll(self, event):
        """Zoom in/out of the matplotlib canvas around the cursor."""
        if event.inaxes is not self.ax:
            return
        base = 1.25
        if event.button == "up":
            scale = 1.0 / base            # zoom IN
        elif event.button == "down":
            scale = base                  # zoom OUT
        else:
            return
        x, y = event.xdata, event.ydata
        xl, xr = self.ax.get_xlim()
        yb, yt = self.ax.get_ylim()
        # Keep (x,y) fixed in data coords while rescaling
        self.ax.set_xlim(x - (x - xl) * scale, x + (xr - x) * scale)
        self.ax.set_ylim(y - (y - yb) * scale, y + (yt - y) * scale)
        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    #  Middle-mouse-button drag → pan the canvas
    # ------------------------------------------------------------------
    def _tk_event_to_data(self, event) -> Optional[Tuple[float, float]]:
        """Convert a Tk mouse event's widget coordinates to axes data
        coordinates. Returns None if the click is outside the axes.

        Tk gives (event.x, event.y) in widget pixels with (0,0) at the
        top-left; matplotlib wants figure pixels with (0,0) at the
        bottom-left. We flip y via the canvas height.
        """
        try:
            w = self.canvas.get_tk_widget()
            h = w.winfo_height()
            fx, fy = float(event.x), float(h - event.y)   # figure pixels
            inv = self.ax.transData.inverted()
            xdata, ydata = inv.transform((fx, fy))
            # Is the point within the axes bbox?
            bbox = self.ax.bbox
            if not (bbox.x0 <= fx <= bbox.x1 and bbox.y0 <= fy <= bbox.y1):
                return None
            return float(xdata), float(ydata)
        except Exception:
            return None

    def _deactivate_toolbar_modes(self):
        """Force matplotlib's toolbar out of Pan/Zoom mode so that
        subsequent clicks behave as picks. A no-op if toolbar is idle.
        """
        tb = getattr(self, "toolbar", None)
        if tb is None:
            return
        try:
            mode = getattr(tb, "mode", "")
            # matplotlib >= 3.3 stores mode as a string like "pan/zoom"
            mode_str = str(mode).lower() if mode else ""
            if "pan" in mode_str:
                tb.pan()        # toggles pan OFF
            elif "zoom" in mode_str:
                tb.zoom()       # toggles zoom OFF
        except Exception:
            pass

    def _on_tk_click(self, event):
        """Tk-level left-click handler — runs BEFORE matplotlib's toolbar
        can swallow the event. Picks the segment under the cursor and
        updates canvas highlight.

        When ``self.click_to_edit`` is True (View menu → "Edit on single
        click"), a successful pick also opens the segment editor, so the
        user can tweak segment parameters by simply clicking the shape
        on the canvas. Works unchanged at any zoom level — hit-testing
        is done in data coordinates and tolerance scales with the view.
        """
        # If toolbar is in pan/zoom mode, turn it off so this click —
        # and subsequent ones — are treated as picks, not pan/zoom drags.
        self._deactivate_toolbar_modes()

        pt = self._tk_event_to_data(event)
        if pt is None:
            return
        xdata, ydata = pt
        sid = self._segment_hit_test(xdata, ydata)
        if sid is None:
            if self.selected_sid is not None:
                self.selected_sid = None
                self.redraw_canvas()
            return
        self.selected_sid = sid
        try:
            self.seg_tree.selection_set(str(sid))
            self.seg_tree.see(str(sid))
        except tk.TclError:
            pass
        self.redraw_canvas()

        # (rev26) If click-to-edit mode is on, open the editor immediately.
        if self.click_to_edit:
            # Defer the dialog to after the current event finishes, so
            # the highlight repaint lands first (and any double-click
            # that may follow on slow trackpads doesn't re-enter).
            self.root.after_idle(self.edit_segment)

    def _on_tk_double_click(self, event):
        """Tk-level double-click → pick + immediately open the segment
        editor. Robust against toolbar mode eating the event."""
        self._deactivate_toolbar_modes()
        pt = self._tk_event_to_data(event)
        if pt is None:
            return
        xdata, ydata = pt
        sid = self._segment_hit_test(xdata, ydata)
        if sid is None:
            return
        self.selected_sid = sid
        try:
            self.seg_tree.selection_set(str(sid))
            self.seg_tree.see(str(sid))
        except tk.TclError:
            pass
        self.edit_segment()

    def _segment_hit_test(self, x: float, y: float) -> Optional[int]:
        """Return the sid of the segment whose rendered polygon (or centreline
        within a small pixel tolerance) contains the data point (x, y).

        Hit-test works identically at any zoom level because:
         1. Polygon containment uses the full rendered geometry in data
            coordinates — unaffected by view scale.
         2. Centreline fallback tolerance combines (a) the segment's own
            half-width, so clicks on the painted stroke always land, and
            (b) a screen-pixel buffer (8 px → data units) so fine
            waveguides remain clickable even when zoomed out.
        """
        # Convert an 8-pixel half-range to data units along both axes.
        # The larger of the two is used so anisotropic views still work.
        try:
            inv = self.ax.transData.inverted()
            p0 = inv.transform((0, 0))
            px_along_x = abs(inv.transform((8, 0))[0] - p0[0])
            px_along_y = abs(inv.transform((0, 8))[1] - p0[1])
            px_tol = max(px_along_x, px_along_y)
        except Exception:
            px_tol = 1.0

        prev_map = self.project._prev_by_start()
        best_sid: Optional[int] = None
        best_d2  = float("inf")
        for seg in self.project.segments:
            start = self.project.get_port(seg.start_pid)
            if start is None:
                continue
            try:
                _res = segment_geometry(
                    seg.kind, start,
                    self.project.effective_params(seg, prev_map),
                    n_pts=200,
                    globals_dict=self.project.eval_globals())
                path    = _res[0]; widths = _res[1]
                extra_p = _res[3] if len(_res) > 3 else None
            except Exception:
                continue

            # 1) Exact polygon containment — the reliable path.
            polys_to_check = []
            if extra_p is not None:
                polys_to_check.extend(extra_p)
            else:
                polys_to_check.append(path_to_polygon(path, widths))
            for P in polys_to_check:
                if _point_in_polygon(x, y, P):
                    return seg.sid

            # 2) Centreline proximity fallback — half_w + 8 px buffer.
            # When zoomed IN  : px_tol is tiny, half_w dominates (clicks
            #                   just outside a thin WG still register).
            # When zoomed OUT : px_tol dominates, so the WG is clickable
            #                   even if its painted width is sub-pixel.
            half_w  = float(widths[0]) / 2.0 if widths.size else 0.0
            eff_tol = half_w + px_tol
            d2 = _min_dist2_to_polyline(x, y, path)
            if d2 < best_d2 and d2 <= eff_tol * eff_tol:
                best_d2  = d2
                best_sid = seg.sid
        return best_sid

    def _on_mouse_press(self, event):
        """Middle/right button → pan. Left-click pick is handled by
        the Tk-level bindings (_on_tk_click / _on_tk_double_click) so
        it works even when the matplotlib toolbar's pan/zoom modes
        would otherwise consume the event.
        """
        if event.inaxes is not self.ax:
            return

        # middle OR right button → pan
        if event.button in (2, 3):
            self._pan_start = dict(
                px=event.x, py=event.y,
                xlim=self.ax.get_xlim(),
                ylim=self.ax.get_ylim(),
                inv=self.ax.transData.inverted().frozen(),
            )
            try:
                self.canvas.get_tk_widget().configure(cursor="fleur")
            except tk.TclError:
                pass
            return
        # (left button intentionally ignored here — see _on_tk_click)

    def _on_mouse_motion(self, event):
        """Shift the view by the pixel delta since press, converted to
        data coordinates via the frozen inverse transform.

        Left-click pick is handled at the Tk level. Here we only
        implement middle/right-drag panning.
        """
        if self._pan_start is None:
            return
        if event.x is None or event.y is None:
            return
        s = self._pan_start
        p0 = s["inv"].transform((s["px"], s["py"]))
        p1 = s["inv"].transform((event.x, event.y))
        ddx = p1[0] - p0[0]
        ddy = p1[1] - p0[1]
        xl, xr = s["xlim"]
        yb, yt = s["ylim"]
        self.ax.set_xlim(xl - ddx, xr - ddx)
        self.ax.set_ylim(yb - ddy, yt - ddy)
        self.canvas.draw_idle()

    def _on_mouse_release(self, event):
        """End panning on middle/right release, or handle a left-click
        pick if the user didn't drag."""
        # pan end (middle / right)
        if event.button in (2, 3):
            if self._pan_start is not None:
                self._pan_start = None
                try:
                    self.canvas.get_tk_widget().configure(cursor="")
                except tk.TclError:
                    pass
            return

        # left-click release is handled by the Tk-level bindings; nothing
        # to do here. (See _on_tk_click / _on_tk_double_click.)
        return

    def _on_seg_tree_select(self):
        """Sync canvas highlight when the user clicks a row in the
        segment list on the right pane."""
        sel = self.seg_tree.selection()
        if not sel:
            new_sid = None
        else:
            try:
                new_sid = int(sel[0])
            except ValueError:
                new_sid = None
        if new_sid != self.selected_sid:
            self.selected_sid = new_sid
            self.redraw_canvas()

    # ------------------------------------------------------------------
    #  Script editor — run / dump
    # ------------------------------------------------------------------
    def _script_status(self, msg: str, ok: bool = True):
        self.script_status.configure(
            text=msg, foreground="#2a8a2a" if ok else "#c72a2a")

    # ==================================================================
    #  Script editor clipboard helpers — direct clipboard API, so they
    #  work even when macOS Tkinter refuses to forward Cmd+C/V/X/A.
    # ==================================================================
    def _script_select_all(self):
        self.script_text.tag_remove("sel", "1.0", "end")
        self.script_text.tag_add("sel", "1.0", "end-1c")
        self.script_text.mark_set("insert", "1.0")
        self.script_text.focus_set()
        return "break"

    def _script_copy(self):
        try:
            text = self.script_text.get("sel.first", "sel.last")
        except tk.TclError:
            return "break"          # no selection
        self.root.clipboard_clear()
        self.root.clipboard_append(text)
        return "break"

    def _script_cut(self):
        try:
            text = self.script_text.get("sel.first", "sel.last")
        except tk.TclError:
            return "break"
        self.root.clipboard_clear()
        self.root.clipboard_append(text)
        self.script_text.delete("sel.first", "sel.last")
        return "break"

    def _script_paste(self):
        try:
            text = self.root.clipboard_get()
        except tk.TclError:
            return "break"          # clipboard empty / non-text
        # replace current selection, if any
        try:
            self.script_text.delete("sel.first", "sel.last")
        except tk.TclError:
            pass
        self.script_text.insert("insert", text)
        self.script_text.see("insert")
        return "break"

    def _script_popup(self, event):
        """Right-click context menu on the script editor."""
        try:
            self._script_menu.tk_popup(event.x_root, event.y_root)
        finally:
            self._script_menu.grab_release()
        return "break"

    # ==================================================================
    #  Turn a selection of segments into a @group(...) definition
    # ==================================================================
    def group_selected_segments(self):
        sel = self.seg_tree.selection()
        if not sel:
            messagebox.showinfo(
                "Group segments",
                "먼저 Segment Sequence 표에서 그룹으로 묶을 세그먼트들을\n"
                "Shift/Ctrl+클릭으로 선택한 뒤 이 버튼을 눌러 주세요.",
                parent=self.root)
            return
        sids = {int(s) for s in sel}
        # preserve original ordering
        segs = [s for s in self.project.segments if s.sid in sids]
        if not segs:
            return

        from tkinter import simpledialog
        name = simpledialog.askstring(
            "New group",
            "Group 이름 (영문 식별자, 예: 'y_splitter', 'dc'):",
            parent=self.root)
        if not name:
            return
        name = name.strip()
        if not name.isidentifier():
            messagebox.showerror(
                "Group name",
                f"'{name}'은(는) Python 식별자가 아닙니다.\n"
                "영문/밑줄/숫자만, 첫 글자는 숫자가 아니어야 합니다.",
                parent=self.root)
            return
        if name in self.project._groups:
            if not messagebox.askyesno(
                    "Overwrite?",
                    f"이미 '{name}' group이 있습니다. 덮어쓰시겠습니까?",
                    parent=self.root):
                return

        # --- build the group body ---
        first_start = segs[0].start_pid
        end_to_local: Dict[int, str] = {}   # seg.end_pid -> "s0"/"s1"/...
        lines = [
            "",
            f'# ---- Group auto-generated from {len(segs)} selected segments ----',
            f'@group("{name}")',
            f'def {name}(start_pid, layer=0):',
        ]
        multiple_entry = False
        for i, s in enumerate(segs):
            # resolve the start_pid reference
            if s.start_pid == first_start:
                start_ref = "start_pid"
            elif s.start_pid in end_to_local:
                start_ref = f"{end_to_local[s.start_pid]}.end_pid"
            else:
                # segment starts from a port outside the selection —
                # keep literal but flag as likely-broken
                start_ref = str(s.start_pid)
                multiple_entry = True
            pstr = ", ".join(f'"{k}": {v!r}' for k, v in s.params.items())
            lines.append(
                f'    s{i} = proj.add_segment("{s.kind}", {start_ref}, '
                f'layer, {{{pstr}}})')
            end_to_local[s.end_pid] = f"s{i}"
        # return last end_pid (most common pattern) — user can edit
        lines.append(f'    return s{len(segs)-1}.end_pid')
        lines.append("")

        code = "\n".join(lines)

        # Insert the @group block at the TOP of the script — before the
        # first `proj.add_*` or `proj.call_group` or any non-comment line
        current = self.script_text.get("1.0", "end")
        import re as _re
        m = _re.search(
            r"^(proj\.add_base_port|proj\.add_segment|proj\.call_group)",
            current, _re.MULTILINE)
        if m:
            line_no = current.count("\n", 0, m.start()) + 1
            insert_idx = f"{line_no}.0"
        else:
            insert_idx = "end"
        self.script_text.insert(insert_idx, code + "\n")

        msg = f"Group '{name}' 정의를 스크립트 상단에 삽입했습니다."
        if multiple_entry:
            msg += "  (주의: 선택에 외부 포트 참조가 있어 일부 start_pid는 리터럴로 남음)"
        self._script_status(msg, ok=True)

    # ==================================================================
    #  Insert a  proj.call_group(...)  line from the Groups tab
    # ==================================================================
    def call_selected_group(self):
        sel = self.group_tree.selection()
        if not sel:
            messagebox.showinfo(
                "Call group",
                "Groups 탭에서 호출할 group을 먼저 선택하세요.\n"
                "(group은 스크립트의 @group 데코레이터로 정의됩니다.)",
                parent=self.root)
            return
        name = sel[0]          # iid is the group name
        fn = self.project._groups.get(name)
        if fn is None:
            messagebox.showerror(
                "Call group", f"Group '{name}'을(를) 찾을 수 없습니다.",
                parent=self.root)
            return

        # If there's no port at all, the user typically wants to start a
        # fresh chain from a new base port. Offer to add one now so the
        # Call dialog immediately has a valid start_pid to suggest.
        if not self.project.ports:
            ans = messagebox.askyesno(
                "No base port",
                f"프로젝트에 port가 없습니다.\n\n"
                f"Group '{name}'을 호출하려면 최소 1개의 port에서 뻗어나가야 합니다.\n"
                f"지금 base port를 하나 추가하시겠습니까?\n\n"
                f"(예를 누르면 PortDialog가 열리고, OK를 누른 뒤 자동으로\n"
                f"스크립트에  proj.add_base_port(...)  줄이 삽입되고\n"
                f"▶ Redraw가 실행된 후 Call 대화창으로 넘어갑니다.)",
                parent=self.root)
            if not ans:
                return
            if not self.project.layers:
                messagebox.showerror(
                    "No layer",
                    "프로젝트에 layer가 하나도 없습니다. 먼저 Layers 탭에서 "
                    "layer를 추가하거나 스크립트에 Layer를 정의하세요.",
                    parent=self.root)
                return
            port_dlg = PortDialog(
                self.root, "Add base port",
                self.project.layers,
                initial={"x": 0.0, "y": 0.0, "angle": 0.0,
                         "width": 2.0, "label": "in", "layer_idx": 0},
                globals_dict=self.project.globals)
            if port_dlg.result is None:
                return
            pr = port_dlg.result
            port_line = (
                f'proj.add_base_port({pr["x"]:g}, {pr["y"]:g}, '
                f'{pr["angle"]:g}, {pr["width"]:g}, '
                f'{pr["layer_idx"]}, "{pr.get("label", "in")}")')
            # Append to the script then rebuild so the port is visible
            current_end = self.script_text.get("end-2c", "end-1c")
            prefix = "" if current_end == "\n" else "\n"
            self.script_text.insert(
                "end", f"{prefix}\n# Base port for group call\n{port_line}\n")
            self.script_text.see("end")
            self.run_script()
            if not self.project.ports:
                messagebox.showerror(
                    "Error",
                    "Base port를 만들지 못했습니다. 스크립트 상단의 에러를 "
                    "확인해 주세요.",
                    parent=self.root)
                return

        dlg = GroupCallDialog(self.root, name, fn, self.project)
        if dlg.result is None:
            return
        kwargs_str = ", ".join(f"{k}={v}" for k, v in dlg.result.items())
        line = f'proj.call_group("{name}"' + \
               (f", {kwargs_str}" if kwargs_str else "") + ")"
        # append at the end of the script with a separator comment
        current_end = self.script_text.get("end-2c", "end-1c")
        prefix = "" if current_end == "\n" else "\n"
        self.script_text.insert(
            "end", f"{prefix}\n# Call '{name}'\n{line}\n")
        self.script_text.see("end")
        if dlg.run_after:
            self.run_script()
            self._script_status(
                f"'{name}(...)' 호출 추가 + 자동 실행 완료.", ok=True)
        else:
            self._script_status(
                f"Call '{name}(...)' 구문을 스크립트 끝에 추가했습니다. "
                "▶ Redraw 눌러 실행하세요.", ok=True)

    # ==================================================================
    #  Special-device menu actions
    # ==================================================================
    def add_spiral_device(self):
        """Open the spiral-design dialog and insert a ready-to-run
        spiral delay line into the script."""
        if not self.project.layers:
            messagebox.showerror(
                "No layer",
                "프로젝트에 layer가 하나도 없습니다. 먼저 Layers 탭에서 "
                "layer를 추가해 주세요.",
                parent=self.root)
            return
        dlg = SpiralDesignDialog(self.root, self.project)
        if dlg.result is None:
            return
        r = dlg.result

        # Build the snippet — optionally wrap as a group for reuse.
        lines = []
        if r["ensure_base_port"]:
            lines.append(
                f'# Base port for spiral')
            lines.append(
                f'proj.add_base_port({r["x0"]:g}, {r["y0"]:g}, '
                f'{r["angle0"]:g}, {r["width"]:g}, '
                f'{r["layer_idx"]}, "{r["label"]}")')
            lines.append("")
            start_pid_expr = "proj.ports[-1].pid"
        else:
            start_pid_expr = str(r["start_pid"])

        if r["wrap_as_group"]:
            # Include straight_L only if set (non-zero) so that the
            # dispatched style's default is used otherwise.
            has_sL = r.get("straight_L") not in (None, 0, 0.0)
            lines.extend([
                f'@group("{r["group_name"]}")',
                f'def {r["group_name"]}(start_pid, layer=0,',
                f'          radius_min={r["radius_min"]:g},',
                f'          pitch={r["pitch"]:g},',
                f'          n_turns={r["n_turns"]:g},',
                f'          width={r["width"]:g},',
                f'          style={r["style"]!r}' +
                (f',\n          straight_L={r["straight_L"]:g}):'
                 if has_sL else '):'),
                f'    """Parametric spiral delay line (style={r["style"]!r})."""',
            ])
            lines.append('    params = {"radius_min": radius_min,')
            lines.append('              "pitch": pitch,')
            lines.append('              "n_turns": n_turns,')
            lines.append('              "width": width,')
            lines.append('              "style": style}')
            if has_sL:
                lines.append('    params["straight_L"] = straight_L')
            lines.append(
                '    s = proj.add_segment("spiral", start_pid, layer, params)')
            lines.append('    return s.end_pid')
            lines.append("")
            lines.append(f'# Call spiral group')
            lines.append(
                f'proj.call_group("{r["group_name"]}", '
                f'start_pid={start_pid_expr}, layer={r["layer_idx"]})')
        else:
            params_bits = [
                f'"radius_min": {r["radius_min"]:g}',
                f'"pitch": {r["pitch"]:g}',
                f'"n_turns": {r["n_turns"]:g}',
                f'"width": {r["width"]:g}',
                f'"style": {r["style"]!r}',
            ]
            if r.get("straight_L"):
                params_bits.append(f'"straight_L": {r["straight_L"]:g}')
            lines.append('# Spiral delay line')
            lines.append(
                f'proj.add_segment("spiral", {start_pid_expr}, '
                f'{r["layer_idx"]}, {{' + ", ".join(params_bits) + '})')

        snippet = "\n".join(lines) + "\n"
        # Append to the script with a blank-line separator
        current_end = self.script_text.get("end-2c", "end-1c")
        prefix = "" if current_end == "\n" else "\n"
        self.script_text.insert("end", prefix + "\n" + snippet)
        self.script_text.see("end")

        if r["run_after"]:
            self.run_script()
            if r.get("target_length"):
                self._script_status(
                    f"Spiral 추가 완료. 총 길이 = {r['target_length']:.1f} µm.",
                    ok=True)
            else:
                self._script_status(
                    "Spiral waveguide 추가 + 자동 실행 완료.", ok=True)
        else:
            msg = "Spiral snippet을 스크립트에 추가했습니다. "
            if r.get("target_length"):
                msg += f"(설계 길이 ≈ {r['target_length']:.1f} µm) "
            msg += "▶ Redraw 눌러 실행하세요."
            self._script_status(msg, ok=True)

    def run_script(self):
        """Rebuild the project from scratch by executing the script text.

        The script executes in a namespace where `proj` is a fresh Project
        (preserving the current layer definitions AND global variables).
        After the script runs, whatever `proj` points to becomes the new
        project. Global variables defined in the UI are auto-injected by
        name, and a convenience @group('name') decorator registers the
        decorated function as a reusable group on the project."""
        src = self.script_text.get("1.0", "end")
        new_proj = Project()
        # Preserve user's layer definitions and global variables
        new_proj.layers = [Layer(**asdict(l)) for l in self.project.layers]
        new_proj.globals = dict(self.project.globals)

        ns = {
            "proj":    new_proj,
            "Project": Project,
            "Layer":   Layer,
            "Port":    Port,
            "Segment": Segment,
            "__builtins__": __builtins__,
        }
        # --- inject globals by name ---------------------------------
        # Skip names that would collide with the built-in injections
        reserved = set(ns.keys())
        for k, v in new_proj.globals.items():
            if k in reserved:
                continue
            ns[k] = v
        # --- @group('name') decorator -------------------------------
        def group(name):
            def wrap(func):
                new_proj.register_group(name, func)
                return func
            return wrap
        ns["group"] = group

        import traceback
        try:
            exec(compile(src, "<script>", "exec"), ns)
        except SyntaxError as e:
            self._script_status(
                f"Syntax error (line {e.lineno}): {e.msg}", ok=False)
            return
        except Exception as e:
            tb = traceback.format_exc()
            self._script_status(f"{type(e).__name__}: {e}", ok=False)
            messagebox.showerror("Script error", tb)
            return

        resulting = ns.get("proj")
        if not isinstance(resulting, Project):
            self._script_status("Script left no valid `proj`.", ok=False)
            return
        self.project = resulting
        self.refresh_all()
        self.fit_view()
        self._script_status(
            f"OK — {len(self.project.ports)} ports, "
            f"{len(self.project.segments)} segments, "
            f"{len(self.project._groups)} groups.", ok=True)

    def dump_to_script(self):
        """Serialise the current project back into an editable script."""
        lines = [
            "# Auto-generated from the current project.",
            "# Edit and click  ▶  Redraw from script  to apply.",
            "",
        ]

        # --- pid remapping ------------------------------------------------
        # When this dumped script runs, each add_base_port and each
        # add_segment produces one new port whose pid starts from 1 and
        # counts up. We simulate that sequence here and remap all stored
        # pids (start_pid) into the "new" numbering the script will see.
        old_to_new: Dict[int, int] = {}
        next_new = 1

        # Globals first — replay them so variables referenced in params
        # (e.g. "R", "(width+gap)/2") actually exist at redraw time.
        if self.project.globals:
            lines.append("# ---- Globals ----")
            for k, v in self.project.globals.items():
                lines.append(f"proj.globals[{k!r}] = {v!r}")
            lines.append("")

        # Base ports, sorted by original pid for stable ordering
        base_ports = sorted(
            (p for p in self.project.ports if p.is_base),
            key=lambda q: q.pid)
        if base_ports:
            lines.append("# ---- Base ports ----")
        for p in base_ports:
            old_to_new[p.pid] = next_new
            next_new += 1
            lines.append(
                f'proj.add_base_port({p.x:g}, {p.y:g}, {p.angle:g}, '
                f'{p.width:g}, {p.layer_idx}, "{p.label}")'
            )

        if self.project.segments:
            lines.append("")
            lines.append("# ---- Segments ----")
        for s in self.project.segments:
            mapped_start = old_to_new.get(s.start_pid, s.start_pid)
            # each add_segment also generates one new auto-endpoint pid
            old_to_new[s.end_pid] = next_new
            next_new += 1
            pstr = ", ".join(f'"{k}": {v!r}' for k, v in s.params.items())
            lines.append(
                f'proj.add_segment("{s.kind}", {mapped_start}, '
                f'{s.layer_idx}, {{{pstr}}})'
            )
        self.script_text.delete("1.0", "end")
        self.script_text.insert("1.0", "\n".join(lines) + "\n")
        self._script_status("Dumped current project to script.", ok=True)

    # ==================================================================
    #  Layer actions
    # ==================================================================
    def add_layer(self):
        dlg = LayerDialog(self.root, "Add Layer")
        if dlg.result is None: return
        self.project.add_layer(Layer(**dlg.result))
        self.refresh_all()

    def edit_layer(self):
        sel = self.layer_tree.selection()
        if not sel: return
        idx = int(sel[0])
        L = self.project.layers[idx]
        dlg = LayerDialog(self.root, "Edit Layer", initial=asdict(L))
        if dlg.result is None: return
        self.project.layers[idx] = Layer(**dlg.result)
        self.refresh_all()

    def delete_layer(self):
        sel = self.layer_tree.selection()
        if not sel: return
        idx = int(sel[0])
        try:
            self.project.delete_layer(idx)
        except RuntimeError as e:
            messagebox.showerror("Delete layer", str(e)); return
        self.refresh_all()

    # ==================================================================
    #  Port actions
    # ==================================================================
    def add_port(self):
        dlg = PortDialog(self.root, "Add Base Port", self.project.layers,
                         globals_dict=self.project.globals)
        if dlg.result is None: return
        r = dlg.result
        self.project.add_base_port(r["x"], r["y"], r["angle"],
                                   r["width"], r["layer_idx"], r["label"])
        self.project.rebuild()
        self.refresh_all()

    def edit_port(self):
        sel = self.port_tree.selection()
        if not sel: return
        pid = int(sel[0])
        p = self.project.get_port(pid)
        if p is None or not p.is_base: return
        dlg = PortDialog(self.root, f"Edit Port {pid}",
                         self.project.layers, initial=asdict(p),
                         globals_dict=self.project.globals)
        if dlg.result is None: return
        r = dlg.result
        p.x = r["x"]; p.y = r["y"]; p.angle = r["angle"]
        p.width = r["width"]; p.layer_idx = r["layer_idx"]; p.label = r["label"]
        self.project.rebuild()
        self.refresh_all()

    def delete_port(self):
        sel = self.port_tree.selection()
        if not sel: return
        pid = int(sel[0])
        try:
            self.project.delete_port(pid)
        except RuntimeError as e:
            messagebox.showerror("Delete port", str(e)); return
        self.project.rebuild()
        self.refresh_all()

    # ==================================================================
    #  Global-variable actions
    # ==================================================================
    def add_global(self):
        dlg = GlobalDialog(self.root, "Add Global Variable",
                           existing=set(self.project.globals.keys()))
        if dlg.result is None: return
        name, value = dlg.result
        self.project.globals[name] = value
        self.project.rebuild()        # propagate to auto-port positions
        self.refresh_all()

    def edit_global(self):
        sel = self.glob_tree.selection()
        if not sel: return
        name = sel[0]
        dlg = GlobalDialog(self.root, f"Edit Global '{name}'",
                           initial={"name": name,
                                    "value": self.project.globals[name]},
                           existing=set(self.project.globals.keys()) - {name})
        if dlg.result is None: return
        new_name, new_val = dlg.result
        if new_name != name:
            self.project.globals.pop(name, None)
        self.project.globals[new_name] = new_val
        self.project.rebuild()        # propagate to auto-port positions
        self.refresh_all()

    def delete_global(self):
        sel = self.glob_tree.selection()
        if not sel: return
        name = sel[0]
        if messagebox.askyesno("Delete global",
                               f"Remove global variable '{name}'?"):
            self.project.globals.pop(name, None)
            self.project.rebuild()    # propagate (may raise for symbolic refs)
            self.refresh_all()

    # ==================================================================
    #  Segment actions
    # ==================================================================
    def _selected_layer_idx(self) -> int:
        raw = self.layer_var.get()
        try:
            return int(raw.split(":")[0])
        except Exception:
            return 0

    def add_segment(self):
        kind = self._selected_kind()
        try:
            start_pid = int(self.start_pid_var.get())
        except ValueError:
            messagebox.showerror("Add segment", "Start port # must be integer.")
            return
        if self.project.get_port(start_pid) is None:
            messagebox.showerror("Add segment",
                                 f"Port {start_pid} does not exist.")
            return

        # End port # — blank means auto-assign (default). If the user types
        # a number it must be a positive, currently-free port id.
        end_txt = self.end_pid_var.get().strip()
        end_pid = None
        if end_txt:
            try:
                end_pid = int(end_txt)
            except ValueError:
                messagebox.showerror(
                    "Add segment",
                    "End port # must be an integer (or blank for auto).")
                return
            if end_pid <= 0:
                messagebox.showerror("Add segment",
                                     "End port # must be a positive integer.")
                return
            if self.project.get_port(end_pid) is not None:
                messagebox.showerror(
                    "Add segment",
                    f"Port {end_pid} already exists. Choose a free number "
                    f"(or leave blank for auto).")
                return
        # so users can type global-variable names AND helper calls like
        # yoffset(n) (segment n's absolute end-y) in numeric fields
        G = self.project.eval_globals()

        # Width — validate that it evaluates NOW; store the raw expression
        # so that subsequent changes to globals keep flowing through.
        width_txt = self.width_var.get()
        try:
            width_val = _parse_num(width_txt, G)
        except ValueError as e:
            messagebox.showerror("Add segment", f"Invalid width:\n{e}")
            return
        if width_val <= 0:
            messagebox.showerror("Add segment", "Width must be > 0.")
            return
        params = {"width": _store_param(width_txt, width_val)}

        # Offsets — same pattern; also omit trivially-zero literals
        for k, var in [("xoffset", self.xoff_var),
                       ("yoffset", self.yoff_var),
                       ("angleoffset", self.aoff_var)]:
            txt = var.get()
            try:
                v = _parse_num(txt, G)
            except ValueError as e:
                messagebox.showerror("Add segment", f"Invalid {k}:\n{e}")
                return
            stored = _store_param(txt, v)
            # keep symbolic values always; drop numeric zeros for compactness
            if isinstance(stored, str) or stored != 0.0:
                params[k] = stored

        # Type-specific params — same parametric policy
        for key, _, _ in SEG_PARAMS[kind]:
            txt = self.param_vars[key].get()
            try:
                v = _parse_num(txt, G)
            except ValueError as e:
                messagebox.showerror("Add segment", f"Invalid {key}:\n{e}")
                return
            params[key] = _store_param(txt, v)

        self._push_undo()
        try:
            seg = self.project.add_segment(
                kind, start_pid, self._selected_layer_idx(), params,
                end_pid=end_pid)
        except RuntimeError as e:
            self._undo_stack.pop()      # add failed → discard the snapshot
            messagebox.showerror("Add segment", str(e))
            return
        self.refresh_all()
        # auto-chain: next segment starts at the new endpoint, and the
        # end-port box reverts to auto so the chosen id isn't reused.
        self.start_pid_var.set(str(seg.end_pid))
        self.end_pid_var.set("")

    def edit_segment(self):
        sel = self.seg_tree.selection()
        if not sel: return
        sid = int(sel[0])
        seg = next((s for s in self.project.segments if s.sid == sid), None)
        if seg is None: return
        dlg = SegmentDialog(self.root, f"Edit Segment {sid}",
                            self.project.layers, seg,
                            globals_dict=self.project.eval_globals())
        if dlg.result is None: return
        self._push_undo()
        seg.kind = dlg.result["kind"]
        seg.start_pid = dlg.result["start_pid"]
        seg.layer_idx = dlg.result["layer_idx"]
        seg.params = dlg.result["params"]
        self.project.rebuild()
        self.refresh_all()

    def delete_segment(self):
        sel = self.seg_tree.selection()
        if not sel: return
        self._push_undo()
        # delete every selected row (selectmode is "extended")
        for iid in sel:
            sid = int(iid)
            self.project.delete_segment(sid)
            if self.selected_sid == sid:
                self.selected_sid = None
        self.refresh_all()

    def move_segment(self, delta: int):
        sel = self.seg_tree.selection()
        if not sel: return
        self._push_undo()
        sid = int(sel[0])
        self.project.move_segment(sid, delta)
        self.refresh_all()

    # ------------------------------------------------------------------
    #  (rev27) Undo / Select-all / Delete-all  for the Segment Sequence
    # ------------------------------------------------------------------
    def _push_undo(self):
        """Snapshot the current project before a mutating action."""
        try:
            self._undo_stack.append(self.project.to_dict())
        except Exception:
            return
        if len(self._undo_stack) > self._UNDO_LIMIT:
            self._undo_stack.pop(0)

    def undo(self):
        """Restore the most recent snapshot pushed by _push_undo()."""
        if not self._undo_stack:
            messagebox.showinfo("Undo", "되돌릴 작업이 없습니다.")
            return
        snap = self._undo_stack.pop()
        groups = self.project._groups   # preserve script-defined groups
        self.project = Project.from_dict(snap)
        self.project._groups = groups
        self.selected_sid = None
        self.refresh_all()
        self.fit_view()

    def select_all_segments(self):
        items = self.seg_tree.get_children("")
        if not items: return
        self.seg_tree.selection_set(items)
        self.seg_tree.focus(items[0])

    def delete_all_segments(self):
        if not self.project.segments:
            messagebox.showinfo("Delete all", "삭제할 세그먼트가 없습니다.")
            return
        if not messagebox.askyesno(
                "Delete all",
                f"세그먼트 {len(self.project.segments)}개를 모두 삭제할까요?\n"
                f"(Undo 로 복원할 수 있습니다.)"):
            return
        self._push_undo()
        # delete by sid; copy the list since delete_segment mutates it
        for sid in [s.sid for s in list(self.project.segments)]:
            self.project.delete_segment(sid)
        self.selected_sid = None
        self.refresh_all()

    # ==================================================================
    #  Project IO
    # ==================================================================
    def new_project(self):
        if not messagebox.askyesno("New", "Discard current project?"): return
        self.project = Project()
        self.refresh_all()

    def save_project(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON project", "*.json"), ("All", "*.*")])
        if not path: return
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.project.to_dict(), f, indent=2, ensure_ascii=False)
        messagebox.showinfo("Saved", f"Project written to {path}")

    def load_project(self):
        path = filedialog.askopenfilename(
            filetypes=[("JSON project", "*.json"), ("All", "*.*")])
        if not path: return
        try:
            with open(path, "r", encoding="utf-8") as f:
                d = json.load(f)
            self.project = Project.from_dict(d)
        except Exception as e:
            messagebox.showerror("Open", f"Failed to load: {e}"); return
        self.refresh_all()
        self.fit_view()

    # ==================================================================
    #  GDS export
    # ==================================================================
    def export_gds(self):
        if not _HAS_GDSFACTORY:
            messagebox.showerror(
                "gdsfactory missing",
                "gdsfactory is not installed.\n"
                "    pip install gdsfactory")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".gds",
            filetypes=[("GDS", "*.gds"), ("All", "*.*")])
        if not path: return
        try:
            self._ensure_pdk_active()          # gdsfactory 9.x requires this
            c = self._make_unique_component(path)
            prev_map = self.project._prev_by_start()
            for seg in self.project.segments:
                start = self.project.get_port(seg.start_pid)
                if start is None: continue
                _res = segment_geometry(
                    seg.kind, start,
                    self.project.effective_params(seg, prev_map),
                    n_pts=300,
                    globals_dict=self.project.eval_globals())
                pxy = _res[0]; widths = _res[1]
                extra_polys = _res[3] if len(_res) > 3 else None
                L = self.project.layers[seg.layer_idx]
                if extra_polys is not None:
                    for P in extra_polys:
                        c.add_polygon(P, layer=(L.number, L.datatype))
                else:
                    poly = path_to_polygon(pxy, widths)
                    c.add_polygon(poly, layer=(L.number, L.datatype))
            # also add base ports for downstream tools
            for p in self.project.ports:
                if not p.is_base: continue
                L = self.project.layers[p.layer_idx]
                try:
                    c.add_port(name=f"P{p.pid}",
                               center=(p.x, p.y),
                               width=p.width,
                               orientation=p.angle,
                               layer=(L.number, L.datatype))
                except Exception:
                    pass
            c.write_gds(path)
        except Exception as e:
            messagebox.showerror("Export GDS", f"Failed: {e}"); return
        messagebox.showinfo("Exported",
                            f"GDS written to {path}\n(cell name: {c.name})")

    def _ensure_pdk_active(self):
        """gdsfactory 9.x refuses to create Components when no PDK is
        active. This helper tries several known activation paths (the
        API has shifted across versions) and silently swallows the
        failure if nothing is available — in which case the Component
        constructor itself will raise a clearer error."""
        # If there's already an active PDK, do nothing.
        try:
            gf.pdk.get_active_pdk()
            return
        except Exception:
            pass
        # Path 1 — gdsfactory.gpdk.PDK (most common on 9.x)
        try:
            gf.gpdk.PDK.activate()
            return
        except AttributeError:
            pass
        except Exception:
            pass
        # Path 2 — generic_tech.get_generic_pdk()
        try:
            from gdsfactory.generic_tech import get_generic_pdk
            get_generic_pdk().activate()
            return
        except Exception:
            pass
        # Path 3 — submodule import
        try:
            import gdsfactory.gpdk as _gpdk
            _gpdk.PDK.activate()
            return
        except Exception:
            pass
        # Path 4 — construct a bare PDK and activate it
        try:
            from gdsfactory.pdk import Pdk
            Pdk(name="generic").activate()
            return
        except Exception:
            pass
        # last-ditch: try the 8.x style
        try:
            gf.pdk.ACTIVE_PDK = gf.pdk.Pdk(name="generic")
            return
        except Exception:
            pass

    def _make_unique_component(self, path: str):
        """Create a gdsfactory Component whose cell name is unique in the
        current KCLayout. gdsfactory 9.x refuses duplicate cell names in
        the same process, which breaks second-and-later exports when
        reusing a fixed name like "waveguide_layout". Strategy:
          1. Base the cell name on the filename stem (sanitised).
          2. Try `allow_duplicate=True` first if the installed version
             supports it (gdsfactory/kfactory newer builds).
          3. Fall back to appending "_1", "_2", … until creation succeeds.
        """
        import re
        from pathlib import Path
        stem = Path(path).stem or "waveguide_layout"
        stem = re.sub(r"[^A-Za-z0-9_]", "_", stem)
        if stem and stem[0].isdigit():
            stem = "wg_" + stem

        last_err: Optional[Exception] = None
        for i in range(500):
            candidate = stem if i == 0 else f"{stem}_{i}"
            # try with allow_duplicate=True first (newer gdsfactory/kfactory)
            try:
                return gf.Component(candidate, allow_duplicate=True)
            except TypeError:
                pass        # kwarg not supported by this version
            except Exception as e:
                last_err = e
            try:
                return gf.Component(candidate)
            except Exception as e:
                msg = str(e).lower()
                if "already exist" in msg or "duplicate" in msg:
                    last_err = e
                    continue
                raise
        raise RuntimeError(
            f"Could not create a unique cell name based on '{stem}' "
            f"after 500 attempts. Last error: {last_err}")


# ============================================================================
#  Dialog boxes
# ============================================================================

class LayerDialog(tk.Toplevel):
    def __init__(self, parent, title, initial=None):
        super().__init__(parent); self.title(title)
        self.result = None
        d = initial or {}
        self.vars = {
            "name":     tk.StringVar(value=d.get("name", "Core")),
            "number":   tk.StringVar(value=str(d.get("number", 1))),
            "datatype": tk.StringVar(value=str(d.get("datatype", 0))),
            "width":    tk.StringVar(value=str(d.get("width", 2.0))),
        }
        self.color = d.get("color", "#1f77b4")

        r = 0
        for k, lbl in [("name", "Name"), ("number", "GDS layer"),
                       ("datatype", "Datatype"), ("width", "Default width [µm]")]:
            ttk.Label(self, text=lbl).grid(row=r, column=0, sticky=tk.W,
                                           padx=6, pady=4)
            ttk.Entry(self, textvariable=self.vars[k], width=18
                      ).grid(row=r, column=1, padx=6, pady=4)
            r += 1
        ttk.Label(self, text="Color").grid(row=r, column=0, sticky=tk.W,
                                           padx=6, pady=4)
        self.swatch = tk.Label(self, text="    ", bg=self.color,
                               width=8, relief=tk.SUNKEN)
        self.swatch.grid(row=r, column=1, sticky=tk.W, padx=6, pady=4)
        ttk.Button(self, text="…", command=self._pick).grid(
            row=r, column=1, sticky=tk.E, padx=6)
        r += 1
        ttk.Button(self, text="OK",     command=self._ok).grid(row=r, column=0, pady=8)
        ttk.Button(self, text="Cancel", command=self.destroy).grid(row=r, column=1, pady=8)
        self.transient(parent); self.grab_set(); self.wait_window(self)

    def _pick(self):
        c = colorchooser.askcolor(initialcolor=self.color, parent=self)
        if c and c[1]:
            self.color = c[1]; self.swatch.configure(bg=self.color)

    def _ok(self):
        try:
            self.result = dict(
                name=self.vars["name"].get().strip() or "Core",
                number=int(self.vars["number"].get()),
                datatype=int(self.vars["datatype"].get()),
                color=self.color,
                width=float(self.vars["width"].get()),
            )
        except ValueError as e:
            messagebox.showerror("Invalid", str(e), parent=self); return
        self.destroy()


class PortDialog(tk.Toplevel):
    def __init__(self, parent, title, layers: List[Layer], initial=None,
                 globals_dict=None):
        super().__init__(parent); self.title(title)
        self.result = None
        self._G = globals_dict or {}
        d = initial or {}
        self.vars = {
            "x":      tk.StringVar(value=str(d.get("x", 0.0))),
            "y":      tk.StringVar(value=str(d.get("y", 0.0))),
            "angle":  tk.StringVar(value=str(d.get("angle", 0.0))),
            "width":  tk.StringVar(value=str(d.get("width",
                                                   layers[0].width if layers else 2.0))),
            "label":  tk.StringVar(value=d.get("label", "")),
        }
        self.layer_var = tk.StringVar()
        r = 0
        for k, lbl in [("x","x [µm]"),("y","y [µm]"),("angle","angle [deg]"),
                       ("width","width [µm]"),("label","label")]:
            ttk.Label(self, text=lbl).grid(row=r, column=0, sticky=tk.W, padx=6, pady=4)
            ttk.Entry(self, textvariable=self.vars[k], width=18
                      ).grid(row=r, column=1, padx=6, pady=4)
            r += 1
        ttk.Label(self, text="Layer").grid(row=r, column=0, sticky=tk.W, padx=6, pady=4)
        vals = [f"{i}: {L.name}" for i, L in enumerate(layers)]
        cb = ttk.Combobox(self, textvariable=self.layer_var, values=vals,
                          state="readonly", width=16)
        cb.grid(row=r, column=1, padx=6, pady=4)
        cb.current(d.get("layer_idx", 0) if 0 <= d.get("layer_idx", 0) < len(vals) else 0)
        r += 1
        if self._G:
            ttk.Label(self, foreground="#666",
                      text=f"Globals available: {', '.join(sorted(self._G))}"
                      ).grid(row=r, column=0, columnspan=2,
                             padx=6, pady=(0, 4), sticky=tk.W)
            r += 1
        ttk.Button(self, text="OK",     command=self._ok).grid(row=r, column=0, pady=8)
        ttk.Button(self, text="Cancel", command=self.destroy).grid(row=r, column=1, pady=8)
        self.transient(parent); self.grab_set(); self.wait_window(self)

    def _ok(self):
        try:
            G = self._G
            self.result = dict(
                x=_parse_num(self.vars["x"].get(), G),
                y=_parse_num(self.vars["y"].get(), G),
                angle=_parse_num(self.vars["angle"].get(), G),
                width=_parse_num(self.vars["width"].get(), G),
                label=self.vars["label"].get(),
                layer_idx=int(self.layer_var.get().split(":")[0]),
            )
        except Exception as e:
            messagebox.showerror("Invalid", str(e), parent=self); return
        self.destroy()


class SegmentDialog(tk.Toplevel):
    def __init__(self, parent, title, layers: List[Layer], seg: Segment,
                 globals_dict=None):
        super().__init__(parent); self.title(title)
        self.result = None
        self.layers = layers
        self._G = globals_dict or {}

        self.kind_var = tk.StringVar(value=seg.kind)
        self.start_var = tk.StringVar(value=str(seg.start_pid))
        self.layer_var = tk.StringVar()

        # Initial display values — params may be numbers OR symbolic
        # strings ("width", "B", "(width+gap)/2"). Keep strings as-is,
        # format numbers with :g. Previously this crashed with
        # ValueError when float() hit a symbolic param.
        def _initial_text(key, numeric_fallback=0.0):
            v = seg.params.get(key, None)
            if isinstance(v, str):
                return v                              # symbolic → keep
            if v is None:
                return f"{numeric_fallback:g}"
            try:
                return f"{float(v):g}"
            except (ValueError, TypeError):
                return str(v)

        # width — special case: blank/zero falls back to layer's default
        w_val = seg.params.get("width", None)
        if isinstance(w_val, str) and w_val.strip():
            self.width_var = tk.StringVar(value=w_val)
        else:
            try:
                init_w = float(w_val or 0.0)
            except (TypeError, ValueError):
                init_w = 0.0
            if init_w <= 0:
                init_w = layers[seg.layer_idx].width if layers else 2.0
            self.width_var = tk.StringVar(value=f"{init_w:g}")

        # offsets (may be symbolic too, e.g. yoffset="B")
        self.xoff_var = tk.StringVar(value=_initial_text("xoffset"))
        self.yoff_var = tk.StringVar(value=_initial_text("yoffset"))
        self.aoff_var = tk.StringVar(value=_initial_text("angleoffset"))

        r = 0
        ttk.Label(self, text="Kind").grid(row=r, column=0, sticky=tk.W, padx=6, pady=4)
        kinds = list(SEG_PARAMS.keys())
        cb = ttk.Combobox(self, textvariable=self.kind_var, values=kinds,
                          state="readonly", width=16)
        cb.grid(row=r, column=1, padx=6, pady=4)
        cb.bind("<<ComboboxSelected>>", lambda e: self._rebuild(seg.params))
        r += 1
        ttk.Label(self, text="Start port #").grid(row=r, column=0, sticky=tk.W, padx=6, pady=4)
        ttk.Entry(self, textvariable=self.start_var, width=8
                  ).grid(row=r, column=1, padx=6, pady=4, sticky=tk.W); r += 1

        ttk.Label(self, text="Layer").grid(row=r, column=0, sticky=tk.W, padx=6, pady=4)
        vals = [f"{i}: {L.name}" for i, L in enumerate(layers)]
        cbl = ttk.Combobox(self, textvariable=self.layer_var, values=vals,
                           state="readonly", width=16)
        cbl.grid(row=r, column=1, padx=6, pady=4)
        cbl.current(min(seg.layer_idx, len(vals)-1))
        r += 1

        ttk.Label(self, text="Width [µm]").grid(row=r, column=0,
                                                 sticky=tk.W, padx=6, pady=4)
        ttk.Entry(self, textvariable=self.width_var, width=12
                  ).grid(row=r, column=1, padx=6, pady=4, sticky=tk.W); r += 1

        # offsets — three side-by-side entries
        ttk.Label(self, text="Offset  x, y, θ°").grid(
            row=r, column=0, sticky=tk.W, padx=6, pady=4)
        off_row = ttk.Frame(self)
        off_row.grid(row=r, column=1, padx=6, pady=4, sticky=tk.W)
        ttk.Entry(off_row, textvariable=self.xoff_var, width=6
                  ).pack(side=tk.LEFT, padx=(0, 3))
        ttk.Entry(off_row, textvariable=self.yoff_var, width=6
                  ).pack(side=tk.LEFT, padx=3)
        ttk.Entry(off_row, textvariable=self.aoff_var, width=5
                  ).pack(side=tk.LEFT, padx=3)
        r += 1

        self.params_frame = ttk.Frame(self)
        self.params_frame.grid(row=r, column=0, columnspan=2, padx=6, pady=4); r += 1
        self.param_vars: Dict[str, tk.StringVar] = {}
        self._rebuild(seg.params)

        if self._G:
            ttk.Label(self, foreground="#666",
                      text=f"Globals available: {', '.join(sorted(self._G))}"
                      ).grid(row=r, column=0, columnspan=2, padx=6, pady=(2, 0),
                             sticky=tk.W)
            r += 1

        bf = ttk.Frame(self); bf.grid(row=r, column=0, columnspan=2, pady=8)
        ttk.Button(bf, text="OK", command=self._ok).pack(side=tk.LEFT, padx=6)
        ttk.Button(bf, text="Cancel", command=self.destroy).pack(side=tk.LEFT, padx=6)
        self.transient(parent); self.grab_set(); self.wait_window(self)

    def _rebuild(self, current_params: Dict):
        for ch in self.params_frame.winfo_children():
            ch.destroy()
        self.param_vars.clear()
        kind = self.kind_var.get()
        for key, label, default in SEG_PARAMS[kind]:
            row = ttk.Frame(self.params_frame); row.pack(fill=tk.X, pady=2)
            ttk.Label(row, text=label, width=18, anchor=tk.W
                      ).pack(side=tk.LEFT)
            val = current_params.get(key, default)
            # strings (symbolic params) pass through; numbers use :g
            if isinstance(val, str):
                text = val
            else:
                try:
                    text = f"{float(val):g}"
                except (ValueError, TypeError):
                    text = str(val)
            v = tk.StringVar(value=text)
            ttk.Entry(row, textvariable=v, width=16).pack(side=tk.LEFT)
            self.param_vars[key] = v

    def _ok(self):
        try:
            G = self._G
            kind = self.kind_var.get()

            # Width
            width_txt = self.width_var.get()
            width_val = _parse_num(width_txt, G)
            if width_val <= 0:
                raise ValueError("Width must be > 0.")
            params = {"width": _store_param(width_txt, width_val)}

            # Type-specific params — keep symbolic
            for k, _, _ in SEG_PARAMS[kind]:
                txt = self.param_vars[k].get()
                v = _parse_num(txt, G)    # validates
                params[k] = _store_param(txt, v)

            # Offsets — only include when meaningful
            for k, var in [("xoffset", self.xoff_var),
                           ("yoffset", self.yoff_var),
                           ("angleoffset", self.aoff_var)]:
                txt = var.get()
                v = _parse_num(txt, G)
                stored = _store_param(txt, v)
                if isinstance(stored, str) or stored != 0.0:
                    params[k] = stored

            self.result = dict(
                kind=kind,
                start_pid=int(self.start_var.get()),
                layer_idx=int(self.layer_var.get().split(":")[0]),
                params=params,
            )
        except Exception as e:
            messagebox.showerror("Invalid", str(e), parent=self); return
        self.destroy()


class GlobalDialog(tk.Toplevel):
    """Edit a single global-variable name/value pair. Value accepts any
    Python literal (number, 'string', [list], (tuple), {dict}, True/False,
    None). Unparsable text is stored as-is as a string."""
    def __init__(self, parent, title, initial=None, existing=None):
        super().__init__(parent); self.title(title)
        self.result = None
        self._existing = existing or set()

        d = initial or {}
        self.name_var = tk.StringVar(value=d.get("name", ""))
        self.value_var = tk.StringVar(value=repr(d.get("value", 0.0)))

        r = 0
        ttk.Label(self, text="Name").grid(row=r, column=0, sticky=tk.W,
                                          padx=8, pady=4)
        ttk.Entry(self, textvariable=self.name_var, width=22
                  ).grid(row=r, column=1, padx=8, pady=4); r += 1

        ttk.Label(self, text="Value").grid(row=r, column=0, sticky=tk.W,
                                           padx=8, pady=4)
        ttk.Entry(self, textvariable=self.value_var, width=28
                  ).grid(row=r, column=1, padx=8, pady=4); r += 1

        hint = ttk.Label(
            self, foreground="#666", justify=tk.LEFT, wraplength=280,
            text=("Examples:   2.0  |  'Core'  |  [1, 2, 3]  |  True\n"
                  "Name must be a valid Python identifier (a-z, 0-9, _)."))
        hint.grid(row=r, column=0, columnspan=2, padx=8, pady=(0, 4))
        r += 1

        bf = ttk.Frame(self); bf.grid(row=r, column=0, columnspan=2, pady=8)
        ttk.Button(bf, text="OK",     command=self._ok
                   ).pack(side=tk.LEFT, padx=6)
        ttk.Button(bf, text="Cancel", command=self.destroy
                   ).pack(side=tk.LEFT, padx=6)

        self.transient(parent); self.grab_set(); self.wait_window(self)

    def _ok(self):
        import ast, keyword
        name = self.name_var.get().strip()
        if not name:
            messagebox.showerror("Invalid", "Name is required.",
                                 parent=self); return
        if not name.isidentifier():
            messagebox.showerror(
                "Invalid",
                f"{name!r} is not a valid Python identifier.", parent=self)
            return
        if keyword.iskeyword(name):
            messagebox.showerror("Invalid",
                                 f"{name!r} is a reserved keyword.",
                                 parent=self); return
        if name in self._existing:
            messagebox.showerror("Invalid",
                                 f"{name!r} is already defined.",
                                 parent=self); return
        value_str = self.value_var.get().strip()
        try:
            value = ast.literal_eval(value_str) if value_str else None
        except (ValueError, SyntaxError):
            # fall back to raw string if it couldn't be parsed as a literal
            value = value_str
        self.result = (name, value)
        self.destroy()


class GroupCallDialog(tk.Toplevel):
    """Prompt for the parameters of a registered group and build the
    corresponding ``proj.call_group(...)`` invocation.

    The ``start_pid`` parameter gets special treatment: the field is
    pre-filled with the most recently generated port (most common
    chaining target), a list of available ports is shown for reference,
    and OK-time validation rejects non-existent or invalid pid values.
    """

    def __init__(self, parent, name: str, fn, project):
        super().__init__(parent)
        self.title(f"Call group — {name}")
        self.result = None
        self.run_after = False
        self._project = project

        import inspect
        try:
            sig = inspect.signature(fn)
            params = list(sig.parameters.values())
        except (TypeError, ValueError):
            params = []
        self._params = params

        info = ttk.Label(
            self, foreground="#555",
            text=("각 파라미터의 값을 Python 표현식으로 입력하세요.\n"
                  "예) 1, 6.0, '+1*B', width*2, proj.ports[-1].pid\n"
                  "빈칸으로 두면 해당 인자는 생략됩니다 (default 사용)."),
            justify=tk.LEFT, wraplength=440)
        info.grid(row=0, column=0, columnspan=3, padx=8, pady=(8, 4),
                  sticky=tk.W)

        # Available-ports reference label (truncated if many)
        if project.ports:
            def _short(pr):
                tag = ("★" if pr.is_base else "")
                lbl = f"({pr.label})" if pr.label else ""
                return f"{pr.pid}{tag}{lbl}"
            port_names = [_short(p) for p in project.ports]
            MAX = 30
            if len(port_names) <= MAX:
                ports_txt = ", ".join(port_names)
            else:
                ports_txt = (", ".join(port_names[:MAX])
                             + f" … +{len(port_names)-MAX} more")
            ttk.Label(self, foreground="#0a5",
                      text=f"Available ports (★=base): {ports_txt}",
                      wraplength=440, justify=tk.LEFT
                      ).grid(row=1, column=0, columnspan=3,
                             padx=8, pady=(0, 6), sticky=tk.W)
        else:
            ttk.Label(self, foreground="#a55",
                      text=("프로젝트에 port가 없습니다. "
                            "이 group을 호출하기 전에 base port를 추가해 주세요."),
                      wraplength=440, justify=tk.LEFT
                      ).grid(row=1, column=0, columnspan=3,
                             padx=8, pady=(0, 6), sticky=tk.W)

        # Compute a smart default for start_pid — the most recent port
        # (typically the endpoint you want to chain off of)
        suggested_start_pid = (project.ports[-1].pid
                               if project.ports else None)

        self.vars: Dict[str, tk.StringVar] = {}
        self._order: List[str] = []
        base_row = 2
        for i, p in enumerate(params):
            required = (p.default is inspect.Parameter.empty)
            if required:
                default_txt, hint = "", "(필수)"
            else:
                default_txt = repr(p.default)
                hint = f"(default: {default_txt})"

            # Smart suggestion for start_pid
            if p.name == "start_pid" and suggested_start_pid is not None:
                default_txt = str(suggested_start_pid)
                hint = (f"(제안: port #{suggested_start_pid})"
                        if required
                        else f"(default: {p.default}, 제안: #{suggested_start_pid})")

            r = base_row + i
            ttk.Label(self, text=p.name).grid(
                row=r, column=0, sticky=tk.W, padx=8, pady=2)
            var = tk.StringVar(value=default_txt)
            self.vars[p.name] = var
            self._order.append(p.name)
            ttk.Entry(self, textvariable=var, width=26).grid(
                row=r, column=1, padx=6, pady=2, sticky=tk.W)
            ttk.Label(self, text=hint, foreground="#666").grid(
                row=r, column=2, padx=4, pady=2, sticky=tk.W)

        if not params:
            ttk.Label(self, text="(이 group은 인자가 없습니다.)",
                      foreground="#666").grid(row=base_row, column=0,
                                              columnspan=3, padx=8, pady=8)

        self.run_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            self, text="호출 삽입 후 ▶ Redraw 자동 실행",
            variable=self.run_var
        ).grid(row=base_row+len(params), column=0, columnspan=3,
               padx=8, pady=(8, 2), sticky=tk.W)

        bf = ttk.Frame(self)
        bf.grid(row=base_row+len(params)+1, column=0, columnspan=3,
                pady=10)
        ttk.Button(bf, text="OK",     command=self._ok
                   ).pack(side=tk.LEFT, padx=6)
        ttk.Button(bf, text="Cancel", command=self.destroy
                   ).pack(side=tk.LEFT, padx=6)

        self.transient(parent); self.grab_set(); self.wait_window(self)

    # -------------------- validation & build --------------------------
    def _ok(self):
        result = {}
        for name in self._order:
            v = self.vars[name].get().strip()
            if v:
                result[name] = v

        # Special validation for start_pid: must be integer AND exist
        if "start_pid" in result:
            import math
            env = {"__builtins__": {"abs": abs, "min": min, "max": max,
                                    "round": round, "int": int, "float": float,
                                    "pow": pow}}
            for n in ("sin","cos","tan","sqrt","exp","log","pi","e"):
                env[n] = getattr(math, n, None)
            env.update(self._project.globals)
            env["proj"] = self._project
            try:
                val = eval(result["start_pid"], env)
                pid = int(val)
            except Exception as e:
                messagebox.showerror(
                    "Invalid start_pid",
                    f"start_pid 값 '{result['start_pid']}'을(를) "
                    f"정수로 해석할 수 없습니다.\n\n{type(e).__name__}: {e}",
                    parent=self)
                return
            if self._project.get_port(pid) is None:
                available = [p.pid for p in self._project.ports]
                messagebox.showerror(
                    "Invalid start_pid",
                    f"Port #{pid} 이(가) 프로젝트에 없습니다.\n\n"
                    f"사용 가능한 port 번호: {available[:30]}"
                    + (" …" if len(available) > 30 else ""),
                    parent=self)
                return
        else:
            # start_pid is missing — is it required by the function?
            import inspect
            for p in self._params:
                if (p.name == "start_pid"
                        and p.default is inspect.Parameter.empty):
                    messagebox.showerror(
                        "Missing start_pid",
                        "start_pid는 필수 인자입니다. 유효한 port 번호를 입력하세요.",
                        parent=self)
                    return

        self.run_after = self.run_var.get()
        self.result = result
        self.destroy()


class SpiralDesignDialog(tk.Toplevel):
    """Parametric designer for an Archimedean double-spiral delay line.

    Lets the user pick the minimum bend radius, turn pitch, number of
    turns and core width, with live feedback on the resulting path
    length, footprint, and outer radius. Optionally creates a base
    port if the project is empty, and optionally wraps the result in
    a @group("spiral_delay") definition for reuse.
    """

    def __init__(self, parent, project):
        super().__init__(parent)
        self.title("Spiral waveguide — design")
        self.result = None
        self._project = project

        # --- form vars with sensible photonic defaults ---------------
        # Match typical silica-PLC / Si₃N₄ specs.
        self.style      = tk.StringVar(value="archimedean")
        self.radius_min = tk.StringVar(value="250.0")
        self.pitch      = tk.StringVar(value="20.0")
        self.n_turns    = tk.StringVar(value="4")
        self.width      = tk.StringVar(value="2.0")
        self.straight_L = tk.StringVar(value="")   # blank → style default
        # Target-length mode: if target_len > 0, solver adjusts one
        # parameter automatically.
        self.target_len     = tk.StringVar(value="")
        self.target_tune_by = tk.StringVar(value="auto")

        self.layer_var = tk.StringVar()
        self.ensure_base_port = tk.BooleanVar(value=(not project.ports))
        self.x0      = tk.StringVar(value="0.0")
        self.y0      = tk.StringVar(value="0.0")
        self.angle0  = tk.StringVar(value="0.0")
        self.label   = tk.StringVar(value="spiral_in")
        self.start_pid = tk.StringVar(
            value=str(project.ports[-1].pid) if project.ports else "1")

        self.wrap_as_group = tk.BooleanVar(value=True)
        self.group_name    = tk.StringVar(value="spiral_delay")
        self.run_after     = tk.BooleanVar(value=True)

        # --- intro ----------------------------------------------------
        ttk.Label(self, foreground="#555",
                  text=("나선 지연선을 네 가지 topology 중 하나로 설계합니다.\n"
                        "각 스타일마다 입출력 포트의 기하학적 배치가 다릅니다.\n"
                        "모든 스타일에서 곡률 ≤ 1/R_min 을 유지합니다."),
                  justify=tk.LEFT, wraplength=520
                  ).grid(row=0, column=0, columnspan=3,
                         padx=10, pady=(10, 6), sticky=tk.W)

        # --- style selector --------------------------------------------
        sf = ttk.LabelFrame(self, text="Spiral topology")
        sf.grid(row=1, column=0, columnspan=3, padx=10, pady=4, sticky="ew")
        styles = [
            ("archimedean",
             "Archimedean (Image 1) — yin-yang. Backend: gf.components.spiral_double"),
            ("racetrack",
             "Racetrack (Image 2) — rectangular expanding. Backend: gf.components.spiral"),
            ("racetrack_s",
             "Racetrack + S (Image 3) — oval + central S. Backend: gf.components.spiral_racetrack"),
            ("archimedean_x",
             "Archimedean + Crossing (Image 4) — single-arm + radial crossing (self-impl), 입력 -90° 자동"),
        ]
        for i, (val, lbl) in enumerate(styles):
            ttk.Radiobutton(sf, text=lbl, variable=self.style,
                            value=val
                            ).grid(row=i, column=0, columnspan=3,
                                   padx=8, pady=2, sticky=tk.W)

        # --- physical parameters -------------------------------------
        fr = ttk.LabelFrame(self, text="Physical parameters")
        fr.grid(row=2, column=0, columnspan=3, padx=10, pady=4, sticky="ew")

        rows = [
            ("Minimum bend radius  R_min [µm]",  self.radius_min,
             "곡률 한계. gdsfactory 파라미터: spiral_double.min_bend_radius / "
             "spiral_racetrack.min_radius (spiral에는 직접 없음)"),
            ("Turn pitch  Λ [µm]",               self.pitch,
             "인접 loop 간격. gdsfactory 파라미터: spiral_double.separation / "
             "spiral.spacing / spiral_racetrack.spacings[i]"),
            ("Number of turns  N",               self.n_turns,
             "감는 수. gdsfactory 파라미터: spiral_double.number_of_loops / "
             "spiral.n_loops / spiral_racetrack.len(spacings)"),
            ("Core width  w [µm]",               self.width,
             "단일 모드 조건에 맞춰 플랫폼 별로 설계. gdsfactory는 cross_section "
             "의 width 를 사용 (강제 0.5).  이 값은 시각화/메타데이터용."),
            ("Straight length [µm] (opt.)",      self.straight_L,
             "Racetrack 전용. Style 2: spiral.length (total); "
             "Style 3: spiral_racetrack.straight_length. 빈칸=기본값."),
        ]
        for i, (lbl, var, hint) in enumerate(rows):
            ttk.Label(fr, text=lbl).grid(
                row=i, column=0, sticky=tk.W, padx=8, pady=3)
            ttk.Entry(fr, textvariable=var, width=10).grid(
                row=i, column=1, padx=6, pady=3, sticky=tk.W)
            ttk.Label(fr, foreground="#888", text=hint, wraplength=300,
                      justify=tk.LEFT).grid(
                row=i, column=2, padx=6, pady=3, sticky=tk.W)

        # Live derived metrics — updated whenever any var changes
        self.info_lbl = ttk.Label(fr, foreground="#0a5",
                                   text="", justify=tk.LEFT)
        self.info_lbl.grid(row=len(rows), column=0, columnspan=3,
                           padx=8, pady=(6, 8), sticky=tk.W)
        for v in (self.radius_min, self.pitch, self.n_turns,
                   self.width, self.style, self.straight_L):
            v.trace_add("write", lambda *a: self._update_info())
        self._update_info()

        # --- target length (optional solver) --------------------------
        tl = ttk.LabelFrame(self, text="Target total length (optional)")
        tl.grid(row=3, column=0, columnspan=3, padx=10, pady=4, sticky="ew")
        ttk.Label(tl, foreground="#555",
                  text=("여기 값을 입력하면 위에 지정한 n_turns "
                        "또는 straight_L 을 자동으로 조정해 총 "
                        "길이를 정확히 맞춥니다. 빈 칸이면 수동 모드."),
                  justify=tk.LEFT, wraplength=520).grid(
                      row=0, column=0, columnspan=3,
                      padx=8, pady=(4, 4), sticky=tk.W)
        ttk.Label(tl, text="Target length [µm]").grid(
            row=1, column=0, sticky=tk.W, padx=8, pady=3)
        ttk.Entry(tl, textvariable=self.target_len, width=12).grid(
            row=1, column=1, sticky=tk.W, padx=6, pady=3)
        ttk.Label(tl, foreground="#888",
                  text="예) 200000  (= 20 cm),   0 = 미사용",
                  wraplength=260).grid(
            row=1, column=2, sticky=tk.W, padx=6, pady=3)

        ttk.Label(tl, text="Tune parameter").grid(
            row=2, column=0, sticky=tk.W, padx=8, pady=3)
        tune_cb = ttk.Combobox(tl, textvariable=self.target_tune_by,
                               values=["auto", "n_turns", "straight_L"],
                               state="readonly", width=14)
        tune_cb.grid(row=2, column=1, sticky=tk.W, padx=6, pady=3)
        ttk.Label(tl, foreground="#888",
                  text=("auto: Archimedean → n_turns, Racetrack → straight_L"),
                  wraplength=260).grid(
            row=2, column=2, sticky=tk.W, padx=6, pady=3)

        # --- placement ------------------------------------------------
        pl = ttk.LabelFrame(self, text="Placement")
        pl.grid(row=4, column=0, columnspan=3, padx=10, pady=4, sticky="ew")

        cb = ttk.Checkbutton(pl,
                             text="새 base port 만들기 (선택 해제 시 기존 port 사용)",
                             variable=self.ensure_base_port,
                             command=self._toggle_port_fields)
        cb.grid(row=0, column=0, columnspan=3, padx=8, pady=3, sticky=tk.W)

        port_rows = [
            ("x [µm]",            self.x0),
            ("y [µm]",            self.y0),
            ("angle θ [deg]",     self.angle0),
            ("label",             self.label),
        ]
        self.port_entries = []
        for i, (lbl, var) in enumerate(port_rows):
            ttk.Label(pl, text=lbl).grid(
                row=1+i, column=0, sticky=tk.W, padx=(24, 4), pady=2)
            e = ttk.Entry(pl, textvariable=var, width=14)
            e.grid(row=1+i, column=1, padx=4, pady=2, sticky=tk.W)
            self.port_entries.append(e)

        # existing-port fallback
        ttk.Label(pl, text="또는 기존 start_pid:").grid(
            row=len(port_rows)+1, column=0, sticky=tk.W,
            padx=(24, 4), pady=(8, 2))
        self.start_pid_entry = ttk.Entry(
            pl, textvariable=self.start_pid, width=14)
        self.start_pid_entry.grid(
            row=len(port_rows)+1, column=1, padx=4, pady=(8, 2),
            sticky=tk.W)
        if project.ports:
            ttk.Label(pl, foreground="#888",
                      text=f"(현재 ports: {[p.pid for p in project.ports][:12]}"
                           + (" …" if len(project.ports) > 12 else "") + ")"
                      ).grid(row=len(port_rows)+1, column=2,
                             padx=4, pady=(8, 2), sticky=tk.W)

        ttk.Label(pl, text="Layer").grid(
            row=len(port_rows)+2, column=0, sticky=tk.W,
            padx=(24, 4), pady=2)
        vals = [f"{i}: {L.name}" for i, L in enumerate(project.layers)]
        lb = ttk.Combobox(pl, textvariable=self.layer_var, values=vals,
                          state="readonly", width=16)
        lb.grid(row=len(port_rows)+2, column=1, padx=4, pady=2,
                sticky=tk.W)
        lb.current(0)

        self._toggle_port_fields()

        # Auto-set angle0 to -90° for Style 4 (archimedean_x) so the
        # spiral hangs DOWNWARD from the base port, matching Image 4.
        def _on_style_change(*_):
            if self.style.get() == "archimedean_x":
                self.angle0.set("-90.0")
            else:
                # Revert to 0 only if it was set to -90 by us (to not
                # override user's manual edit)
                if self.angle0.get().strip() == "-90.0":
                    self.angle0.set("0.0")
        self.style.trace_add("write", _on_style_change)
        _on_style_change()

        # --- options --------------------------------------------------
        opt = ttk.LabelFrame(self, text="Code generation")
        opt.grid(row=5, column=0, columnspan=3, padx=10, pady=4, sticky="ew")
        ttk.Checkbutton(opt, text="Group 으로 정의 (@group 래핑, 재사용 쉬움)",
                        variable=self.wrap_as_group,
                        command=self._toggle_group_name).grid(
            row=0, column=0, columnspan=2, padx=8, pady=3, sticky=tk.W)
        ttk.Label(opt, text="group 이름").grid(
            row=1, column=0, padx=(24, 4), sticky=tk.W)
        self.group_name_entry = ttk.Entry(
            opt, textvariable=self.group_name, width=20)
        self.group_name_entry.grid(row=1, column=1, padx=4,
                                   sticky=tk.W)
        ttk.Checkbutton(opt, text="스크립트에 삽입 후 ▶ Redraw 자동 실행",
                        variable=self.run_after).grid(
            row=2, column=0, columnspan=2, padx=8, pady=3, sticky=tk.W)

        # --- OK/Cancel ------------------------------------------------
        bf = ttk.Frame(self)
        bf.grid(row=6, column=0, columnspan=3, pady=10)
        ttk.Button(bf, text="OK",     command=self._ok
                   ).pack(side=tk.LEFT, padx=6)
        ttk.Button(bf, text="Cancel", command=self.destroy
                   ).pack(side=tk.LEFT, padx=6)

        self.transient(parent); self.grab_set(); self.wait_window(self)

    def _update_info(self, *_):
        try:
            R = float(self.radius_min.get())
            p = float(self.pitch.get())
            n = float(self.n_turns.get())
        except (ValueError, TypeError):
            self.info_lbl.configure(text="(숫자 입력 시 계산됨)")
            return
        if R <= 0 or p <= 0 or n <= 0:
            self.info_lbl.configure(text="(모든 값은 양수여야 함)")
            return
        r_outer = R + p * n
        # length ≈ ∫₀^{θmax} √(r² + (dr/dθ)²) dθ with r = R + aθ, a = p/(2π)
        a = p / (2 * math.pi)
        th_max = 2 * math.pi * n
        # analytical arc length of Archimedean spiral:
        import math as _m
        def _arch_len(th):
            t = th
            s = 0.5 * a * (
                t * _m.sqrt(1 + t**2)
                + _m.log(t + _m.sqrt(1 + t**2))
            ) if a > 0 else 0.0
            # But our parametrisation has offset R_min, so integrate
            # numerically for accuracy — use simpson
            return s
        # Numerical length
        import numpy as _np
        ts = _np.linspace(0, th_max, 2000)
        rs = R + a * ts
        drdt = a
        ds = _np.sqrt(rs**2 + drdt**2)
        L_one = float(_np.trapz(ds, ts))
        L_total = 2 * L_one + math.pi * (p/4)  # rough + smooth connector
        self.info_lbl.configure(
            text=(f"Path length ≈ {L_total:.0f} µm    "
                  f"Outer radius ≈ {r_outer:.0f} µm    "
                  f"Footprint ≈ {2*r_outer:.0f} × {2*r_outer:.0f} µm"))

    def _toggle_port_fields(self):
        state = "normal" if self.ensure_base_port.get() else "disabled"
        inv   = "disabled" if self.ensure_base_port.get() else "normal"
        for e in self.port_entries:
            e.configure(state=state)
        self.start_pid_entry.configure(state=inv)

    def _toggle_group_name(self):
        self.group_name_entry.configure(
            state="normal" if self.wrap_as_group.get() else "disabled")

    def _ok(self):
        try:
            R  = float(self.radius_min.get())
            p  = float(self.pitch.get())
            n  = float(self.n_turns.get())
            w  = float(self.width.get())
            if min(R, p, n, w) <= 0:
                raise ValueError("모든 값은 0보다 커야 합니다.")
            layer_idx = int(self.layer_var.get().split(":")[0])
            ensure = bool(self.ensure_base_port.get())
            if ensure:
                x0 = float(self.x0.get()); y0 = float(self.y0.get())
                a0 = float(self.angle0.get())
                start_pid = None
                label = self.label.get().strip() or "spiral_in"
            else:
                start_pid = int(self.start_pid.get())
                if self._project.get_port(start_pid) is None:
                    raise ValueError(
                        f"Port #{start_pid} 이(가) 존재하지 않습니다.")
                x0 = y0 = a0 = 0.0; label = ""

            wrap = bool(self.wrap_as_group.get())
            gname = self.group_name.get().strip()
            if wrap and not gname.isidentifier():
                raise ValueError(
                    f"Group 이름 '{gname}'은(는) Python 식별자가 아닙니다.")

            style = self.style.get()
            sL_txt = self.straight_L.get().strip()
            straight_L = float(sL_txt) if sL_txt else 0.0
            if straight_L < 0:
                raise ValueError("Straight length는 0 이상이어야 합니다.")

            # Target length solver: adjust n_turns or straight_L to
            # hit the specified total path length.
            t_txt = self.target_len.get().strip()
            achieved_L = None
            if t_txt and float(t_txt) > 0:
                target_L = float(t_txt)
                tune = self.target_tune_by.get()
                if tune == "auto":
                    tune = "n_turns" if "archimedean" in style else "straight_L"
                init_sL = straight_L if straight_L > 0 else None
                params_solved, achieved_L = _spiral_solve_length(
                    style, R, p, target_L,
                    n_turns_init=n, straight_L_init=init_sL,
                    tune=tune)
                if params_solved is None:
                    raise ValueError(
                        f"Target length {target_L:g} µm 가 R_min={R}, "
                        f"pitch={p} 조건에서 도달 불가능합니다.\n"
                        f"n_turns 나 pitch 를 바꿔 보세요.")
                # Apply solved value
                if "n_turns" in params_solved:
                    n = params_solved["n_turns"]
                if "straight_L" in params_solved:
                    straight_L = params_solved["straight_L"]
        except Exception as e:
            messagebox.showerror("Invalid", str(e), parent=self); return

        self.result = dict(
            radius_min=R, pitch=p, n_turns=n, width=w,
            style=style, straight_L=straight_L,
            target_length=achieved_L,    # None if no target was set
            layer_idx=layer_idx,
            ensure_base_port=ensure,
            x0=x0, y0=y0, angle0=a0, label=label,
            start_pid=start_pid,
            wrap_as_group=wrap, group_name=gname,
            run_after=bool(self.run_after.get()),
        )
        self.destroy()


# ============================================================================
#  Bootstrap
# ============================================================================

def main():
    root = tk.Tk()
    app = WaveguideApp(root)
    # Run the default script once so the canvas isn't empty on first launch
    app.run_script()
    root.mainloop()


if __name__ == "__main__":
    main()