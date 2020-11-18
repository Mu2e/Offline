//
// Alignment tables for tracker elements
//

#ifndef DbTables_TrkAlignElement_hh
#define DbTables_TrkAlignElement_hh

#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DataProducts/inc/StrawId.hh"
#include "DbTables/inc/DbTable.hh"
#include "DbTables/inc/AlignParams.hh"

namespace mu2e {

  class TrkAlignElement : public DbTable {
  public:

    typedef std::shared_ptr<TrkAlignElement> ptr_t;
    typedef std::shared_ptr<const TrkAlignElement> cptr_t;


    TrkAlignElement(const char* Name, const char* DbName,size_t nrows ):DbTable(Name,DbName, "index,dx,dy,dz,rx,ry,rz"), _nrows(nrows) {}
    const AlignParams& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<AlignParams> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    std::size_t nrowFix() const override { return _nrows; };
    size_t size() const override { return _csv.capacity() + nrow()*sizeof(AlignParams); };

    void addRow(const std::vector<std::string>& columns) override {
      _rows.emplace_back(std::stoi(columns[0]),
			 std::stof(columns[1]),
			 std::stof(columns[2]),
			 std::stof(columns[3]),
			 std::stof(columns[4]),
			 std::stof(columns[5]),
			 std::stof(columns[6])  );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      AlignParams const& r = _rows.at(irow);
      sstream << r.index()<<",";
      sstream << std::fixed << std::setprecision(4);  // 0.1 um displacement precision
      sstream << r.dx()<<",";
      sstream << r.dy()<<",";
      sstream << r.dz()<<",";
      sstream << std::fixed << std::setprecision(6); // 1 urad rotation precision
      sstream << r.rx()<<",";
      sstream << r.ry()<<",";
      sstream << r.rz()<<",";
    }

    void clear() override { _csv.clear(); _rows.clear(); }

  private:
    std::vector<AlignParams> _rows;
    size_t _nrows;
  };

// unique classes for Db usage: not sure why this is needed
  class TrkAlignPanel : public TrkAlignElement {
    public:
      TrkAlignPanel() : TrkAlignElement("TrkAlignPanel","trk.alignpanel",StrawId::_nupanels) {}
  };

  class TrkAlignPlane : public TrkAlignElement {
    public:
      TrkAlignPlane() : TrkAlignElement("TrkAlignPlane","trk.alignplane",StrawId::_nplanes) {}
  };
  
  class TrkAlignTracker : public TrkAlignElement {
    public:
      TrkAlignTracker() : TrkAlignElement("TrkAlignTracker","trk.aligntracker",size_t(1)) {}
  };
  
};
#endif
