#ifndef DbTables_TrkAlignStraw_hh
#define DbTables_TrkAlignStraw_hh

#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "DbTables/inc/DbTable.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TrkAlignStraw : public DbTable {
  public:
    using xyzVec = CLHEP::Hep3Vector; // switch to XYZVec TODO

    typedef std::shared_ptr<TrkAlignStraw> ptr_t;
    typedef std::shared_ptr<const TrkAlignStraw> cptr_t;

    struct Row {
      Row(int index, StrawId const& id,
	  float wire_cal_dV, float wire_cal_dW,
	  float wire_hv_dV, float wire_hv_dW,
	  float straw_cal_dV, float straw_cal_dW,
	  float straw_hv_dV, float straw_hv_dW) :
	_wire_cal_dV(wire_cal_dV), _wire_cal_dW(wire_cal_dW),
	_wire_hv_dV(wire_hv_dV), _wire_hv_dW(wire_hv_dW),
	_straw_cal_dV(straw_cal_dV), _straw_cal_dW(straw_cal_dW),
	_straw_hv_dV(straw_hv_dV), _straw_hv_dW(straw_hv_dW) {}

      StrawId const& id() const { return _id; }
      xyzVec wireDeltaUVW(StrawEnd::End end) const { return end == StrawEnd::cal ? xyzVec(0.0,_wire_cal_dV,_wire_cal_dW) : xyzVec(0.0,_wire_hv_dV,_wire_hv_dW); }
      xyzVec wireDeltaUVW(StrawEnd const& end) const { return wireDeltaUVW(end.end()); }
      xyzVec strawDeltaUVW(StrawEnd::End end) const { return end == StrawEnd::cal ? xyzVec(0.0,_straw_cal_dV,_straw_cal_dW) : xyzVec(0.0,_straw_hv_dV,_straw_hv_dW); }
      xyzVec strawDeltaUVW(StrawEnd const& end) const { return strawDeltaUVW(end.end()); }

      int _index;
      StrawId _id;
      float _wire_cal_dV, _wire_cal_dW; // wire cal end displacements from nominal
      float _wire_hv_dV, _wire_hv_dW; // wire hv end displacements from nominal
      float _straw_cal_dV, _straw_cal_dW; // straw cal end displacements from nominal
      float _straw_hv_dV, _straw_hv_dW; // straw hv end displacements from nominal
    };


    TrkAlignStraw():DbTable("TrkAlignStraw","trk.alignstraw",
	"index,StrawId,wire_cal_dV,wire_cal_dW,wire_hv_dV,wire_hv_dW,straw_cal_dV,straw_cal_dW,straw_hv_dV,straw_hv_dW") {}
    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    size_t nrow() const override { return _rows.size(); };
    size_t nrowFix() const override { return StrawId::_nustraws; };
    size_t size() const override { return _csv.capacity() + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      _rows.emplace_back(
	  std::stoi(columns[0]),
	  StrawId(columns[1]),
	  std::stof(columns[2]),
	  std::stof(columns[3]),
	  std::stof(columns[4]),
	  std::stof(columns[5]),
	  std::stof(columns[6]),
	  std::stof(columns[7]),
	  std::stof(columns[8]),
	  std::stof(columns[9]) );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << r._index <<",";
      sstream << r.id().plane() << "_" << r.id().panel() << "_" << r.id().straw() << ",";
      sstream << std::fixed << std::setprecision(4); 
      sstream << r._wire_cal_dV <<",";
      sstream << r._wire_cal_dW <<",";
      sstream << r._wire_hv_dV <<",";
      sstream << r._wire_hv_dW <<",";
      sstream << r._straw_cal_dV <<",";
      sstream << r._straw_cal_dW <<",";
      sstream << r._straw_hv_dV <<",";
      sstream << r._straw_hv_dW;
    }

    virtual void clear() { _csv.clear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
