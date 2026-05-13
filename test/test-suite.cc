// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <cmath>
#include "biodynamo.h"
#include "pancreatic_tumor_model.h"
#include "params/sim_param.h"

#define TEST_NAME typeid(*this).name()

namespace bdm {
namespace pancreatic_tumor {

// ============================================================================
// Math helper unit tests — no simulation required
// ============================================================================

TEST(MathHelpers, SatZeroOrNegativeInput) {
  EXPECT_DOUBLE_EQ(Sat(0.0, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(Sat(-1.0, 1.0), 0.0);
}

TEST(MathHelpers, SatHalfSaturation) {
  EXPECT_NEAR(Sat(100.0, 100.0), 0.5, 1e-12);
  EXPECT_NEAR(Sat(3.0, 3.0), 0.5, 1e-12);
}

TEST(MathHelpers, SatApproachesOne) {
  EXPECT_NEAR(Sat(1e9, 1.0), 1.0, 1e-6);
}

TEST(MathHelpers, SatLinearApproximation) {
  // x << K => Sat(x, K) ≈ x/K
  EXPECT_NEAR(Sat(0.001, 100.0), 0.001 / 100.0, 1e-8);
}

TEST(MathHelpers, ProbFromRateNonPositiveRate) {
  EXPECT_DOUBLE_EQ(ProbFromRate(0.0, 1.0 / 1440.0), 0.0);
  EXPECT_DOUBLE_EQ(ProbFromRate(-1.0, 1.0 / 1440.0), 0.0);
}

TEST(MathHelpers, ProbFromRateIsValidProbability) {
  double p = ProbFromRate(0.1, 1.0 / 1440.0);
  EXPECT_GT(p, 0.0);
  EXPECT_LT(p, 1.0);
}

TEST(MathHelpers, ProbFromRatePoissonExact) {
  // P = 1 - exp(-rate * dt)
  EXPECT_NEAR(ProbFromRate(1.0, 1.0), 1.0 - std::exp(-1.0), 1e-12);
}

TEST(MathHelpers, ProbFromRateSmallDtLinear) {
  // For small dt: P ≈ rate * dt
  double rate = 0.1;
  double dt   = 1.0 / 1440.0;
  EXPECT_NEAR(ProbFromRate(rate, dt), rate * dt, 1e-8);
}

TEST(MathHelpers, ProbFromRateHighRateStaysBelow1) {
  double p = ProbFromRate(1000.0, 1.0);
  EXPECT_LT(p, 1.0);
  EXPECT_GT(p, 0.999);
}

// ============================================================================
// SimParam defaults — no simulation needed beyond basic construction
// ============================================================================

TEST(SimParamTest, DefaultsAreReasonable) {
  Simulation sim(TEST_NAME);
  const auto* sp = sim.GetParam()->Get<SimParam>();

  EXPECT_GT(sp->total_days, 0);
  EXPECT_GT(sp->dt_minutes, 0.0);
  EXPECT_GT(sp->max_bound, sp->min_bound);
  EXPECT_GT(sp->cell_radius_um, 0.0);
  EXPECT_GT(sp->K_C, 0.0);
  EXPECT_GE(sp->c_base_div, 0.0);
  EXPECT_GE(sp->e_base_birth, 0.0);
  EXPECT_GE(sp->e_base_death, 0.0);
  EXPECT_GE(sp->r_decay, 0.0);
}

TEST(SimParamTest, InitialCountsPositive) {
  Simulation sim(TEST_NAME);
  const auto* sp = sim.GetParam()->Get<SimParam>();

  EXPECT_GT(sp->C0, 0u);
  EXPECT_GT(sp->E0, 0u);
  EXPECT_GT(sp->H0, 0u);
}

TEST(SimParamTest, KConstantsPositive) {
  Simulation sim(TEST_NAME);
  const auto* sp = sim.GetParam()->Get<SimParam>();

  EXPECT_GT(sp->K_C, 0.0);
  EXPECT_GT(sp->K_P, 0.0);
  EXPECT_GT(sp->K_E, 0.0);
  EXPECT_GT(sp->K_N, 0.0);
  EXPECT_GT(sp->K_H, 0.0);
  EXPECT_GT(sp->K_R, 0.0);
}

// ============================================================================
// Agent construction and GlobalCensus
// ============================================================================

TEST(AgentTest, AllCellTypesCanBeCreated) {
  Simulation sim(TEST_NAME);
  auto* rm = sim.GetResourceManager();

  rm->AddAgent(new TumorCell({10, 10, 10}));
  rm->AddAgent(new StellateCell({20, 20, 20}));
  rm->AddAgent(new EffectorTCell({30, 30, 30}));
  rm->AddAgent(new NKCell({40, 40, 40}));
  rm->AddAgent(new HelperTCell({50, 50, 50}));
  rm->AddAgent(new TRegCell({60, 60, 60}));

  EXPECT_EQ(rm->GetNumAgents(), 6u);
}

TEST(AgentTest, TumorCellColorMatchesParam) {
  Simulation sim(TEST_NAME);
  const auto* sp = sim.GetParam()->Get<SimParam>();
  auto* cell = new TumorCell({50, 50, 50});
  sim.GetResourceManager()->AddAgent(cell);
  EXPECT_EQ(cell->GetCellColor(), sp->color_tumor);
}

TEST(AgentTest, GlobalCensusCountsCorrectTypes) {
  Simulation sim(TEST_NAME);
  auto* rm = sim.GetResourceManager();

  rm->AddAgent(new TumorCell({10, 10, 10}));
  rm->AddAgent(new TumorCell({20, 20, 20}));
  rm->AddAgent(new EffectorTCell({30, 30, 30}));
  rm->AddAgent(new NKCell({40, 40, 40}));

  // Reset the step cache so RefreshIfNeeded actually runs
  auto& gc = GlobalCensus::Instance();
  gc.step_cached.store(std::numeric_limits<size_t>::max(), std::memory_order_relaxed);
  gc.RefreshIfNeeded();

  EXPECT_EQ(gc.Get().C, 2u);
  EXPECT_EQ(gc.Get().P, 0u);
  EXPECT_EQ(gc.Get().E, 1u);
  EXPECT_EQ(gc.Get().N, 1u);
  EXPECT_EQ(gc.Get().H, 0u);
  EXPECT_EQ(gc.Get().R, 0u);
}

// ============================================================================
// Smoke test: 2-day run — model must not crash and tumor must remain viable
// ============================================================================

TEST(SmokeTest, TwoDayRunStaysAlive) {
  Simulation sim(TEST_NAME);
  auto* sp_mut = const_cast<SimParam*>(sim.GetParam()->Get<SimParam>());

  sp_mut->total_days       = 2;
  sp_mut->dt_minutes       = 1.0;
  sp_mut->min_bound        = 0.0;
  sp_mut->max_bound        = 200.0;  // large domain to avoid overcrowding
  sp_mut->cell_radius_um   = 6.0;
  sp_mut->use_local_counts = false;
  sp_mut->C0 = 10; sp_mut->P0 = 5;  sp_mut->E0 = 10;
  sp_mut->N0 = 5;  sp_mut->H0 = 10; sp_mut->R0 = 3;
  sp_mut->output_dir = "/tmp/pancreatic_smoke_test";

  auto* rm     = sim.GetResourceManager();
  auto* rng    = sim.GetRandom();
  double lo    = sp_mut->min_bound + sp_mut->cell_radius_um;
  double hi    = sp_mut->max_bound - sp_mut->cell_radius_um;

  auto rpos = [&]() -> Real3 {
    return {rng->UniformDouble(lo, hi),
            rng->UniformDouble(lo, hi),
            rng->UniformDouble(lo, hi)};
  };

  for (size_t i = 0; i < sp_mut->C0; ++i) rm->AddAgent(new TumorCell(rpos()));
  for (size_t i = 0; i < sp_mut->P0; ++i) rm->AddAgent(new StellateCell(rpos()));
  for (size_t i = 0; i < sp_mut->E0; ++i) rm->AddAgent(new EffectorTCell(rpos()));
  for (size_t i = 0; i < sp_mut->N0; ++i) rm->AddAgent(new NKCell(rpos()));
  for (size_t i = 0; i < sp_mut->H0; ++i) rm->AddAgent(new HelperTCell(rpos()));
  for (size_t i = 0; i < sp_mut->R0; ++i) rm->AddAgent(new TRegCell(rpos()));

  ASSERT_EQ(rm->GetNumAgents(), 43u);

  size_t steps = static_cast<size_t>(sp_mut->total_days * 1440.0 / sp_mut->dt_minutes);
  sim.GetScheduler()->Simulate(steps);

  // Model must not go completely extinct in 2 days
  EXPECT_GT(rm->GetNumAgents(), 0u);

  // Tumor cell count should still be positive (c_base_div > 0)
  auto& gc = GlobalCensus::Instance();
  gc.step_cached.store(std::numeric_limits<size_t>::max(), std::memory_order_relaxed);
  gc.RefreshIfNeeded();
  EXPECT_GT(gc.Get().C, 0u);
}

}  // namespace pancreatic_tumor
}  // namespace bdm
