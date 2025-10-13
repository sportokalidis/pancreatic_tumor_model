#include <iostream>
#include "include/ptm/params.h"

using namespace bdm::pancreatic_tumor;

int main() {
  // Get parameters (will load from file if available)
  Params* p = P();
  
  std::cout << "=== Parameter Loading Test ===" << std::endl;
  std::cout << "dt_minutes: " << p->dt_minutes << std::endl;
  std::cout << "C0 (initial tumor cells): " << p->C0 << std::endl;
  std::cout << "P0 (initial PSC): " << p->P0 << std::endl;
  std::cout << "c_base_div: " << p->c_base_div << std::endl;
  std::cout << "gate_C_K: " << p->gate_C_K << std::endl;
  
  return 0;
}