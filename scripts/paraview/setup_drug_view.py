#!/usr/bin/env pvpython
"""
setup_drug_view.py — configure an INTERACTIVE ParaView view of the drug cloud
with a LOCKED color range, so the diffusion volume renders reliably (no more
"visible one run, blank the next"). Unlike render_drug.py this does not save
PNGs — it leaves the pipeline ready for you to scrub/rotate in the GUI.

Run it INSIDE the ParaView GUI so you get an interactive window:
  paraview --script=scripts/paraview/setup_drug_view.py
or in the GUI: View -> Python Shell -> "Run Script..." -> pick this file.

Edit the two settings just below, then (re)run.
"""
# ======================= EDIT THESE =======================
# Point at a PER-RUN copy, e.g. runs/treatment_prod/abr_S1e4_dt24h/abr_s42/paraview
# (Do NOT use the live output/pancreatic_tumor_new — reruns overwrite it while
#  ParaView has it open, which causes the errors.)
VIZ_DIR = "output/pancreatic_tumor_new"
TRANSPARENCY = 1.0   # 1 = default; >1 more transparent, <1 more solid
# ==========================================================

from paraview.simple import *
import glob
import os
import re

CELL_TYPES = ["TumorCell", "StellateCell", "EffectorTCell",
              "NKCell", "HelperTCell", "TRegCell"]


def frame_series(viz_dir, prefix, ext):
    out = []
    for f in glob.glob(os.path.join(viz_dir, "%s-*.%s" % (prefix, ext))):
        m = re.search(r"%s-(\d+)\.%s$" % (re.escape(prefix), ext), os.path.basename(f))
        if m:
            out.append((int(m.group(1)), f))
    out.sort()
    return [f for _, f in out]


def main():
    viz = os.path.abspath(VIZ_DIR)
    drug = next((n for n in ("Gemcitabine", "Abraxane")
                 if frame_series(viz, n, "pvti")), None)
    if not drug:
        print("No Gemcitabine/Abraxane -*.pvti in", viz)
        return
    print("Drug:", drug, "  dir:", viz)

    view = GetActiveViewOrCreate("RenderView")
    view.Background = [0.32, 0.34, 0.43]
    view.OrientationAxesVisibility = 0

    drugsrc = XMLPartitionedImageDataReader(FileName=frame_series(viz, drug, "pvti"))
    drugsrc.UpdatePipeline()
    disp = Show(drugsrc, view)
    disp.SetRepresentationType("Volume")
    ColorBy(disp, ("POINTS", "Substance Concentration"))
    ctf = GetColorTransferFunction("SubstanceConcentration")
    otf = GetOpacityTransferFunction("SubstanceConcentration")

    # --- global max over ALL timesteps, then LOCK the range ---
    tk = GetTimeKeeper()
    times = list(tk.TimestepValues) if tk.TimestepValues else [0.0]
    hi, best_t = 0.0, times[0]
    for t in times:
        UpdatePipeline(time=t, proxy=drugsrc)
        ai = (drugsrc.GetDataInformation().GetPointDataInformation()
              .GetArrayInformation("Substance Concentration"))
        v = ai.GetComponentRange(0)[1] if ai else 0.0
        if v > hi:
            hi, best_t = v, t
    if hi <= 0.0:
        hi = 1e-6
    lo, span = 0.0, hi

    ctf.ColorSpace = "RGB"
    ctf.RGBPoints = [lo, 0.15, 0.35, 0.95,
                     lo + 0.5 * span, 0.65, 0.10, 0.75,
                     hi, 0.95, 0.05, 0.10]
    otf.Points = [lo, 0.00, 0.5, 0.0,
                  lo + 0.05 * span, 0.15, 0.5, 0.0,
                  hi, 0.85, 0.5, 0.0]
    # Never auto-rescale — this is the fix for the field vanishing on zero frames.
    for tf in (ctf, otf):
        try:
            tf.AutomaticRescaleRangeMode = "Never"
        except Exception:
            pass
    try:
        b = drugsrc.GetDataInformation().GetBounds()
        extent = max(b[1] - b[0], b[3] - b[2], b[5] - b[4])
        disp.ScalarOpacityUnitDistance = (extent / 80.0) * TRANSPARENCY
    except Exception:
        pass
    disp.SetScalarBarVisibility(view, True)

    # --- cells overlay (small, translucent) ---
    for ct in CELL_TYPES:
        cf = frame_series(viz, ct, "pvtu")
        if not cf:
            continue
        cr = XMLPartitionedUnstructuredGridReader(FileName=cf)
        cd = Show(cr, view)
        cd.SetRepresentationType("Point Gaussian")
        try:
            cd.ShaderPreset = "Sphere"
            cd.GaussianRadius = 1.5
        except Exception:
            pass
        cd.ColorArrayName = ["POINTS", ""]
        cd.DiffuseColor = [0.9, 0.9, 0.9]
        cd.AmbientColor = [0.9, 0.9, 0.9]
        cd.Opacity = 0.35

    # jump to the frame where the drug is strongest, so it's visible immediately
    view.ViewTime = best_t
    UpdatePipeline(time=best_t, proxy=drugsrc)
    ResetCamera(view)
    try:
        GetActiveCamera().Azimuth(30)
        GetActiveCamera().Elevation(18)
    except Exception:
        pass
    Render(view)
    print("Ready. Range locked to [0, %.5g]; jumped to peak time %.3g." % (hi, best_t))
    print("Scrub the time slider to watch it pulse; the range will NOT reset.")


main()
