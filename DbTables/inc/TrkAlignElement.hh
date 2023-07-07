//
// Alignment tables for tracker elements
//

#ifndef DbTables_TrkAlignElement_hh
#define DbTables_TrkAlignElement_hh

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DbTables/inc/TrkAlignParams.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class TrkAlignElement : public DbTable {
 public:
  typedef std::shared_ptr<TrkAlignElement> ptr_t;
  typedef std::shared_ptr<const TrkAlignElement> cptr_t;

  TrkAlignElement(const char* Name, const char* DbName, size_t nrows) :
      DbTable(Name, DbName, "index,strawid,dx,dy,dz,rx,ry,rz"), _nrows(nrows) {}
  const TrkAlignParams& rowAt(const std::size_t index) const {
    return _rows.at(index);
  }
  std::vector<TrkAlignParams> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  std::size_t nrowFix() const override { return _nrows; };
  size_t size() const override {
    return baseSize() + nrow() * sizeof(TrkAlignParams);
  };
  const std::string orderBy() const { return std::string("index"); }

  void addRow(const std::vector<std::string>& columns) override {
    int index = std::stoi(columns[0]);
    // enforce a strict sequential order
    if (index != int(_rows.size())) {
      throw cet::exception("TRKALIGNELEMENT_BAD_INDEX")
          << "TrkAlignElement::addRow found index out of order: " << index
          << " != " << _rows.size() << "\n";
    }

    _rows.emplace_back(index, StrawId(columns[1]), std::stof(columns[2]),
                       std::stof(columns[3]), std::stof(columns[4]),
                       std::stof(columns[5]), std::stof(columns[6]),
                       std::stof(columns[7]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    TrkAlignParams const& r = _rows.at(irow);
    sstream << r.index() << ",";
    sstream << r.id().plane() << "_" << r.id().panel() << "_" << r.id().straw()
            << ",";
    sstream << std::fixed
            << std::setprecision(4);  // 0.1 um displacement precision
    sstream << r.dx() << ",";
    sstream << r.dy() << ",";
    sstream << r.dz() << ",";
    sstream << std::fixed << std::setprecision(6);  // 1 urad rotation precision
    sstream << r.rx() << ",";
    sstream << r.ry() << ",";
    sstream << r.rz() << ",";
  }

  void clear() override {
    baseClear();
    _rows.clear();
  }

 private:
  std::vector<TrkAlignParams> _rows;
  size_t _nrows;
};

// unique classes for Db usage: not sure why this is needed
class TrkAlignPanel : public TrkAlignElement {
 public:
  constexpr static const char* cxname = "TrkAlignPanel";
  TrkAlignPanel() :
      TrkAlignElement(cxname, "trk.alignpanel", StrawId::_nupanels) {}
};

class TrkAlignPlane : public TrkAlignElement {
 public:
  constexpr static const char* cxname = "TrkAlignPlane";
  TrkAlignPlane() :
      TrkAlignElement(cxname, "trk.alignplane", StrawId::_nplanes) {}
};

class TrkAlignTracker : public TrkAlignElement {
 public:
  constexpr static const char* cxname = "TrkAlignTracker";
  TrkAlignTracker() : TrkAlignElement(cxname, "trk.aligntracker", size_t(1)) {}
};

}  // namespace mu2e
#endif
