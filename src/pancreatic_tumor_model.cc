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
#include "pancreatic_tumor_model.h"

int main(int argc, const char** argv) {
  const int rc = bdm::pancreatic_tumor::Simulate(argc, argv);
  // Simulate() has returned, so BioDynaMo has written its ParaView state file;
  // repair the drug colour range so `bdm view` renders the diffusion correctly,
  // no matter how the binary was launched. No-ops for non-treatment runs.
  bdm::pancreatic_tumor::FixDrugParaviewState();
  return rc;
}
