#!/usr/bin/env pvpython
"""
render_drug.py — render the drug-diffusion field (B-hybrid viz) as a volume with
the CORRECT transfer functions, so you never fight the ParaView GUI.

It sets up, for a drug run's exported frames:
  * Volume representation of the drug concentration field
  * a blue(low) -> white -> red(high) color map
  * an opacity ramp that is 0 at concentration 0 (cleared tissue = invisible)
  * the range LOCKED over ALL timesteps, so the temporal decay reads as
    RED on injection days -> BLUE as it decays -> invisible after treatment
  * the cells overlaid as light spheres for context
and writes one PNG per timestep (and optionally an .ogv movie).

Run with BioDynaMo's ParaView (source thisbdm.sh first so the libs resolve):

  pvpython scripts/paraview/render_drug.py <viz_dir> [options]

<viz_dir> is a folder containing the exported frames, e.g.
  runs/treatment_prod/viz_abr/abr_s42/paraview     (Abraxane-*.pvti + *Cell-*.pvtu)
  output/pancreatic_tumor_new

Options:
  --out DIR        output folder for PNGs (default: <viz_dir>/render_<drug>)
  --drug NAME      Gemcitabine | Abraxane | auto   (default: auto-detect)
  --no-cells       don't overlay the cells
  --size WxH       image size (default: 1100x950)
  --stride N       render every Nth frame (default: 1)
  --movie          also write movie.ogv from the frames
"""
from paraview.simple import *
import argparse
import glob
import os
import re
import sys

CELL_TYPES = ["TumorCell", "StellateCell", "EffectorTCell",
              "NKCell", "HelperTCell", "TRegCell"]


def frame_series(viz_dir, prefix, ext):
    """Sorted list of prefix-<frame>.<ext> (the per-timestep group files)."""
    out = []
    for f in glob.glob(os.path.join(viz_dir, "%s-*.%s" % (prefix, ext))):
        m = re.search(r"%s-(\d+)\.%s$" % (re.escape(prefix), ext), os.path.basename(f))
        if m:
            out.append((int(m.group(1)), f))
    out.sort()
    return [f for _, f in out]


def find_drug(viz_dir, requested):
    if requested and requested != "auto":
        return requested
    for name in ("Gemcitabine", "Abraxane"):
        if frame_series(viz_dir, name, "pvti"):
            return name
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("viz_dir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--drug", default="auto")
    ap.add_argument("--no-cells", dest="cells", action="store_false", default=True)
    ap.add_argument("--size", default="1100x950")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--movie", action="store_true")
    args = ap.parse_args()

    viz = os.path.abspath(args.viz_dir)
    drug = find_drug(viz, args.drug)
    if not drug:
        print("ERROR: no Gemcitabine/Abraxane -*.pvti frames found in", viz)
        return 1
    dfiles = frame_series(viz, drug, "pvti")
    out = args.out or os.path.join(viz, "render_%s" % drug.lower())
    os.makedirs(out, exist_ok=True)
    W, H = (int(x) for x in args.size.lower().split("x"))
    print("Drug: %s  (%d frames)\nOutput: %s" % (drug, len(dfiles), out))

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [W, H]
    view.Background = [0.32, 0.34, 0.43]
    view.OrientationAxesVisibility = 0

    # ---- drug field: volume ----
    drugsrc = XMLPartitionedImageDataReader(FileName=dfiles)
    drugsrc.UpdatePipeline()
    disp = Show(drugsrc, view)
    disp.SetRepresentationType("Volume")
    ColorBy(disp, ("POINTS", "Substance Concentration"))
    ctf = GetColorTransferFunction("SubstanceConcentration")
    otf = GetOpacityTransferFunction("SubstanceConcentration")

    # Lock the range to the global data range over ALL timesteps — this is what
    # makes the decay read as red -> blue -> invisible over time.
    disp.RescaleTransferFunctionToDataRangeOverTime()
    lo, hi = ctf.RGBPoints[0], ctf.RGBPoints[-4]
    if hi <= lo:
        hi = lo + 1e-6
    span = hi - lo

    # blue (low) -> purple (mid) -> red (high); saturated all the way (no washed-
    # out white middle), so mid concentrations stay vivid.
    ctf.ColorSpace = "RGB"
    ctf.RGBPoints = [
        lo,              0.15, 0.35, 0.95,
        lo + 0.5 * span, 0.65, 0.10, 0.75,
        hi,              0.95, 0.05, 0.10,
    ]
    # opacity 0 at value 0 (cleared tissue invisible) -> visible at low -> strong
    # at high, so the cloud reads clearly and still fades to nothing after tx.
    otf.Points = [
        lo,              0.00, 0.5, 0.0,
        lo + 0.05 * span, 0.15, 0.5, 0.0,
        hi,              0.85, 0.5, 0.0,
    ]
    # Denser accumulation so the dilute field is actually visible (smaller = more
    # opaque). Scale to the domain size.
    try:
        b = drugsrc.GetDataInformation().GetBounds()
        extent = max(b[1] - b[0], b[3] - b[2], b[5] - b[4])
        disp.ScalarOpacityUnitDistance = extent / 80.0
    except Exception:
        pass
    disp.SetScalarBarVisibility(view, True)

    # ---- cells overlay ----
    if args.cells:
        for ct in CELL_TYPES:
            cfiles = frame_series(viz, ct, "pvtu")
            if not cfiles:
                continue
            cr = XMLPartitionedUnstructuredGridReader(FileName=cfiles)
            cd = Show(cr, view)
            cd.SetRepresentationType("Point Gaussian")
            try:
                cd.ShaderPreset = "Sphere"
                cd.GaussianRadius = 3.0          # small so they don't bury the cloud
            except Exception:
                pass
            cd.ColorArrayName = ["POINTS", ""]  # solid color (ColorBy(None) errors here)
            cd.AmbientColor = [0.9, 0.9, 0.9]
            cd.DiffuseColor = [0.9, 0.9, 0.9]
            cd.Opacity = 0.35                    # translucent so the drug shows through

    ResetCamera(view)
    cam = GetActiveCamera()
    cam.Azimuth(30)
    cam.Elevation(18)
    view.CenterOfRotation = view.CameraFocalPoint

    # ---- render one PNG per timestep ----
    tk = GetTimeKeeper()
    times = list(tk.TimestepValues) if tk.TimestepValues else [0.0]
    n = 0
    for i, t in enumerate(times):
        if i % args.stride:
            continue
        view.ViewTime = t
        UpdatePipeline(t, drugsrc)
        Render(view)
        SaveScreenshot(os.path.join(out, "frame_%04d.png" % i), view,
                       ImageResolution=[W, H])
        n += 1
    print("Wrote %d PNG frames to %s" % (n, out))

    if args.movie:
        SaveAnimation(os.path.join(out, "movie.ogv"), view,
                      ImageResolution=[W, H], FrameRate=10)
        print("Wrote movie.ogv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
