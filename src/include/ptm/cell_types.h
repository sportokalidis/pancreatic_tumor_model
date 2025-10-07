#ifndef PTM_CELL_TYPES_H_
#define PTM_CELL_TYPES_H_

#include "biodynamo.h"
#include "ptm/params.h"

namespace bdm {
namespace pancreatic_tumor {

// -----------------------------------------------------------------------------
// Agent types
// -----------------------------------------------------------------------------
class TumorCell : public Cell {
  BDM_AGENT_HEADER(TumorCell, Cell, 1);
 public:
  TumorCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.tumor; }
  explicit TumorCell(const Real3& p) : TumorCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class StellateCell : public Cell {
  BDM_AGENT_HEADER(StellateCell, Cell, 1);
 public:
  StellateCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.psc; }
  explicit StellateCell(const Real3& p) : StellateCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class EffectorTCell : public Cell {
  BDM_AGENT_HEADER(EffectorTCell, Cell, 1);
 public:
  EffectorTCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.eff; }
  explicit EffectorTCell(const Real3& p) : EffectorTCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class NKCell : public Cell {
  BDM_AGENT_HEADER(NKCell, Cell, 1);
 public:
  NKCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.nk; }
  explicit NKCell(const Real3& p) : NKCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class HelperTCell : public Cell {
  BDM_AGENT_HEADER(HelperTCell, Cell, 1);
 public:
  HelperTCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.helper; }
  explicit HelperTCell(const Real3& p) : HelperTCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class TRegCell : public Cell {
  BDM_AGENT_HEADER(TRegCell, Cell, 1);
 public:
  TRegCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.treg; }
  explicit TRegCell(const Real3& p) : TRegCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class ReporterCell : public Cell {
  BDM_AGENT_HEADER(ReporterCell, Cell, 1);
 public:
  ReporterCell() { SetDiameter(0.1); }
};

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_CELL_TYPES_H_