#ifndef PTM_REPORTER_BEHAVIOR_H_
#define PTM_REPORTER_BEHAVIOR_H_

#include "biodynamo.h"
#include "ptm/params.h"
#include "ptm/cell_types.h"
#include "ptm/global_census.h"
#include <fstream>
#include <iostream>
#include <cstdlib>

namespace bdm {
namespace pancreatic_tumor {

class ReportPopCounts : public Behavior {
  BDM_BEHAVIOR_HEADER(ReportPopCounts, Behavior, 1);
 public:
  void Run(Agent* /* reporter */) override {
    auto* sim  = Simulation::GetActive();
    auto steps = sim->GetScheduler()->GetSimulatedSteps();

    // report once per day (1440 steps)
    if (steps % 1440 == 0) {
      auto& gc = GlobalCensus::Instance();
      gc.RefreshIfNeeded();

      size_t C  = gc.C;
      size_t Pn = gc.Pn;
      size_t E  = gc.E;
      size_t N  = gc.N;
      size_t H  = gc.H;
      size_t R  = gc.R;

      real_t t_day = steps * P()->dt_minutes / 1440.0;

      // Print to console
      std::cout << "[day " << t_day << "] C=" << C << " P=" << Pn
                << " E=" << E << " N=" << N << " H=" << H << " R=" << R << "\n";

      // Get output directory from environment variable or use current directory
      std::string output_dir = ".";
      const char* env_output_dir = std::getenv("OUTPUT_DIR");
      if (env_output_dir) {
        output_dir = std::string(env_output_dir);
      }
      
      // append to CSV in the specified directory
      std::string csv_path = output_dir + "/populations.csv";
      std::ofstream csv(csv_path, std::ios::app);
      if (steps == 0) {
        csv << "step,day,C,P,E,N,H,R,total\n";
      }
      csv << steps << "," << t_day << ","
          << C << "," << Pn << "," << E << "," << N << "," << H << "," << R << ","
          << (C+Pn+E+N+H+R) << "\n";
      csv.flush();
    }
  }
};

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_REPORTER_BEHAVIOR_H_