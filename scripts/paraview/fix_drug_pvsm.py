#!/usr/bin/env pvbatch
"""
fix_drug_pvsm.py — repair ONLY the drug colour range in BioDynaMo's own ParaView
state file, so `bdm view` renders the diffusion field correctly.

BioDynaMo's generate_pv_state.py colours the substance by
RescaleTransferFunctionToDataRange at frame 0, where the treatment drug is 0
(dosing starts day 14). ParaView clamps that empty range to [0, FLT_MIN]
(1.1757813367477812e-38), so the field renders as a flat, unchanging haze.

This script does NOT rebuild the state — it leaves BioDynaMo's rendering exactly
as it is (sphere cells, its default colour map) and only rescales the drug's
degenerate transfer functions from [0, FLT_MIN] to [0, <real 99th-percentile>].
The true range is read from the actual .vti data via VTK/numpy (BioDynaMo writes
placeholder RangeMin/RangeMax in the files, so ParaView's own queries are
unreliable).

Run (BioDynaMo already has pvbatch on PATH via $ParaView_DIR):
  pvbatch scripts/paraview/fix_drug_pvsm.py <viz_dir>
<viz_dir> defaults to output/pancreatic_tumor_new. Called automatically at the
end of scripts/hpc/run_direct.sh for treatment (viz) runs.
"""
from paraview.simple import *
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy
import glob
import numpy as np
import os
import re
import sys

FLT_MIN = 1.1757813367477812e-38   # ParaView's clamp for an empty [0,0] range


def frame_series(viz_dir, prefix, ext):
    out = []
    for f in glob.glob(os.path.join(viz_dir, "%s-*.%s" % (prefix, ext))):
        m = re.search(r"%s-(\d+)\.%s$" % (re.escape(prefix), ext), os.path.basename(f))
        if m:
            out.append((int(m.group(1)), f))
    out.sort()
    return [f for _, f in out]


def real_max(viz, drug):
    """99th percentile of the busiest frame (robust vs source-box outliers)."""
    reader = XMLPartitionedImageDataReader(FileName=frame_series(viz, drug, "pvti"))
    reader.UpdatePipeline()
    times = list(GetTimeKeeper().TimestepValues) or [0.0]
    hi = 0.0
    for t in times:
        UpdatePipeline(time=t, proxy=reader)
        pd = sm.Fetch(reader).GetPointData()
        for i in range(pd.GetNumberOfArrays()):
            if "Concentration" in pd.GetArrayName(i):
                arr = vtk_to_numpy(pd.GetArray(i))
                if arr.size and arr.max() > 0:
                    hi = max(hi, float(np.percentile(arr, 99)))
                break
    return hi if hi > 0 else 1e-6


def rescale_substance_tfs(pvsm_path, new_max):
    """Rescale ONLY the drug's colour + opacity transfer functions to [0, new_max].

    They are identified by their registered names 'SubstanceConcentration.*' (the
    volume display's LookupTable + ScalarOpacityFunction), NOT by their range —
    BioDynaMo bakes a degenerate range that flips between [0,1] and [0,FLT_MIN]
    across runs, so a value-based test is unreliable. Every other transfer
    function (cells, gradient) is left completely untouched.
    """
    txt = open(pvsm_path).read()

    drug_ids = set()
    for coll in ("lookup_tables", "piecewise_functions"):
        m = re.search(r'<ProxyCollection name="%s">(.*?)</ProxyCollection>' % coll,
                      txt, re.S)
        if m:
            for it in re.finditer(r'<Item id="(\d+)" name="SubstanceConcentration\.',
                                  m.group(1)):
                drug_ids.add(it.group(1))

    def rescale_points(block):
        vals = [float(v) for v in re.findall(r'value="([-0-9.eE+]+)"', block)]
        xs = [vals[i] for i in range(0, len(vals), 4)]   # x is the 1st of each 4-tuple
        mx = max(xs) if xs else 0.0
        if mx <= 0.0:
            return block
        factor = new_max / mx
        new_vals = list(vals)
        for i in range(0, len(vals), 4):
            new_vals[i] = vals[i] * factor
        it = iter(new_vals)
        return re.sub(r'value="[-0-9.eE+]+"',
                      lambda _m: 'value="%s"' % ("%.10g" % next(it)), block)

    def fix_body(body):
        body = re.sub(r'<Property name="RGBPoints".*?</Property>',
                      lambda m: rescale_points(m.group(0)), body, flags=re.S)
        body = re.sub(r'<Property name="Points".*?</Property>',
                      lambda m: rescale_points(m.group(0)), body, flags=re.S)
        return body

    n = 0
    for pid in drug_ids:
        pat = re.compile(r'(<Proxy\b[^>]*\bid="%s"[^>]*>)(.*?)(</Proxy>)' % re.escape(pid),
                         re.S)
        txt, c = pat.subn(lambda m: m.group(1) + fix_body(m.group(2)) + m.group(3), txt)
        n += c
    open(pvsm_path, "w").write(txt)
    return n


def main():
    viz = os.path.abspath(sys.argv[1] if len(sys.argv) > 1
                          else "output/pancreatic_tumor_new")
    sim_name = os.path.basename(viz.rstrip("/")) or "pancreatic_tumor_new"
    pvsm = os.path.join(viz, "%s.pvsm" % sim_name)
    drug = next((n for n in ("Gemcitabine", "Abraxane")
                 if frame_series(viz, n, "pvti")), None)
    if not drug:
        print("[fix_drug_pvsm] no drug field in %s — nothing to do" % viz)
        return 0
    if not os.path.exists(pvsm):
        print("[fix_drug_pvsm] %s not found (BioDynaMo state not generated?)" % pvsm)
        return 1

    hi = real_max(viz, drug)
    n = rescale_substance_tfs(pvsm, hi)
    print("[fix_drug_pvsm] rescaled %d %s transfer function(s) to 0..%.5g"
          % (n, drug, hi))
    sys.stdout.flush()   # pvbatch may segfault on exit; flush so the line survives
    return 0


if __name__ == "__main__":
    sys.exit(main())
