// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey
// BioDynaMo collaboration. Apache-2.0 license.
//
// -----------------------------------------------------------------------------
#ifndef PANCREATIC_TUMOR_MODEL_H_
#define PANCREATIC_TUMOR_MODEL_H_

#include "biodynamo.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <mutex>
#include <vector>
#include "include/ptm/params.h"
#include "include/ptm/helpers.h"
#include "include/ptm/cell_types.h"
#include "include/ptm/global_census.h"
#include "include/ptm/behaviors/tumor_behavior.h"
#include "include/ptm/behaviors/psc_behavior.h"
#include "include/ptm/behaviors/eff_behavior.h"
#include "include/ptm/behaviors/nk_behavior.h"
#include "include/ptm/behaviors/helper_behavior.h"
#include "include/ptm/behaviors/treg_behavior.h"
#include "include/ptm/behaviors/reporter_behavior.h"

namespace {
  std::mutex cout_mtx;
}

namespace bdm {
namespace pancreatic_tumor {

// All components are now included from modular headers:
// - Params struct and P() function from ptm/params.h (with file loading)
// - Helper functions from ptm/helpers.h
// - Cell types from ptm/cell_types.h  
// - GlobalCensus from ptm/global_census.h
// - Behaviors from ptm/behaviors/*.h

// -------------------------------------------------------
// Simulate
// -------------------------------------------------------
inline int Simulate(int argc, const char** argv) {
  // Parse command line arguments for parameter file and output directory
  std::string param_file = "";
  std::string output_dir = ".";
  
  // Create filtered argument list for BioDynaMo (removing our custom args)
  std::vector<const char*> bdm_argv;
  bdm_argv.push_back(argv[0]); // program name
  
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--param-file" && i + 1 < argc) {
      param_file = argv[i + 1];
      setenv("PARAM_FILE", param_file.c_str(), 1);
      i++; // skip next argument
    } else if (arg == "--output-dir" && i + 1 < argc) {
      output_dir = argv[i + 1];
      setenv("OUTPUT_DIR", output_dir.c_str(), 1);
      i++; // skip next argument
    } else {
      // Pass other arguments to BioDynaMo
      bdm_argv.push_back(argv[i]);
    }
  }
  
  // Create output directory if it doesn't exist
  system(("mkdir -p " + output_dir).c_str());
  
  std::cout << "Loaded parameters from " << (param_file.empty() ? "default values" : param_file) << std::endl;
  std::cout << "Output directory: " << output_dir << std::endl;
  
  auto setp = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = P()->min_bound;
    param->max_bound = P()->max_bound;
  };

  Simulation sim(bdm_argv.size(), bdm_argv.data(), setp);
  auto* ctxt = sim.GetExecutionContext();
  auto* rng  = sim.GetRandom();

  auto rand_pos = [&]() -> Real3 {
    real_t x = rng->Uniform(P()->min_bound, P()->max_bound);
    real_t y = rng->Uniform(P()->min_bound, P()->max_bound);
    real_t z = rng->Uniform(P()->min_bound, P()->max_bound);
    return Real3{x, y, z};
  };

  // Seed populations
  for (size_t i=0; i<P()->C0; ++i) {
    auto* c = new TumorCell(ClampPoint(rand_pos()));
    c->SetCellColor(P()->color.tumor);
    c->AddBehavior(new TumorBehavior());
    ctxt->AddAgent(c);
  }
  for (size_t i=0; i<P()->P0; ++i) {
    auto* p = new StellateCell(ClampPoint(rand_pos()));
    p->SetCellColor(P()->color.psc);
    p->AddBehavior(new PSCBehavior());
    ctxt->AddAgent(p);
  }
  for (size_t i=0; i<P()->E0; ++i) {
    auto* e = new EffectorTCell(ClampPoint(rand_pos()));
    e->SetCellColor(P()->color.eff);
    e->AddBehavior(new EffectorBehavior());
    ctxt->AddAgent(e);
  }
  for (size_t i=0; i<P()->N0; ++i) {
    auto* n = new NKCell(ClampPoint(rand_pos()));
    n->SetCellColor(P()->color.nk);
    n->AddBehavior(new NKBehavior());
    ctxt->AddAgent(n);
  }
  for (size_t i=0; i<P()->H0; ++i) {
    auto* h = new HelperTCell(ClampPoint(rand_pos()));
    h->SetCellColor(P()->color.helper);
    h->AddBehavior(new HelperBehavior());
    ctxt->AddAgent(h);
  }
  for (size_t i=0; i<P()->R0; ++i) {
    auto* r = new TRegCell(ClampPoint(rand_pos()));
    r->SetCellColor(P()->color.treg);
    r->AddBehavior(new TRegBehavior());
    ctxt->AddAgent(r);
  }

  // Reporter
  auto* rep = new ReporterCell();
  rep->AddBehavior(new ReportPopCounts());
  sim.GetExecutionContext()->AddAgent(rep);

  // Run ~100 days (1 min step → 1440 steps per day)
  sim.GetScheduler()->Simulate(1440 * 2);
  std::cout << "Pancreatic tumor (C,P,E,N,H,R) ABM (global interactions) completed.\n";
  return 0;
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_TUMOR_MODEL_H_