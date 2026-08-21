#!/usr/bin/env pvbatch
"""
save_movie.py — export an animation video from a BioDynaMo ParaView state file,
without touching the GUI. Loads the .pvsm exactly as `bdm view` would (so it uses
the auto-fixed drug range, spheroid, colours you set) and writes a video and/or a
PNG frame sequence.

Run with BioDynaMo's ParaView (source thisbdm.sh first):
  pvbatch scripts/paraview/save_movie.py [pvsm] [out] [options]

  pvsm   state file (default: output/pancreatic_tumor_new/pancreatic_tumor_new.pvsm)
  out    output path; extension picks the format:
           .ogv / .avi  -> a single video file
           .png         -> a frame sequence out.0000.png, out.0001.png, ...
         (default: <pvsm dir>/movie.ogv)

Options:
  --size WxH     image size (default 1280x1024)
  --fps N        frame rate for video (default 10)
"""
from paraview.simple import *
import os
import sys


def main():
    args = [a for a in sys.argv[1:]]
    size = "1280x1024"
    fps = 10
    pos = []
    i = 0
    while i < len(args):
        if args[i] == "--size":
            size = args[i + 1]; i += 2
        elif args[i] == "--fps":
            fps = int(args[i + 1]); i += 2
        else:
            pos.append(args[i]); i += 1

    pvsm = pos[0] if len(pos) > 0 else "output/pancreatic_tumor_new/pancreatic_tumor_new.pvsm"
    pvsm = os.path.abspath(pvsm)
    out = pos[1] if len(pos) > 1 else os.path.join(os.path.dirname(pvsm), "movie.ogv")
    out = os.path.abspath(out)
    W, H = (int(x) for x in size.lower().split("x"))

    if not os.path.exists(pvsm):
        print("[save_movie] state file not found:", pvsm)
        return 1

    LoadState(pvsm)
    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [W, H]
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()   # one frame per exported timestep

    os.makedirs(os.path.dirname(out), exist_ok=True)
    SaveAnimation(out, view, ImageResolution=[W, H], FrameRate=fps)
    n = len(list(GetTimeKeeper().TimestepValues) or [0])
    print("[save_movie] wrote %s  (%d frames @ %d fps, %dx%d)" % (out, n, fps, W, H))
    sys.stdout.flush()   # pvbatch may segfault on exit after writing — harmless
    return 0


if __name__ == "__main__":
    sys.exit(main())
