#ifndef DbTables_TrkAlignPlane_hh
#define DbTables_TrkAlignPlane_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"
#include "GeneralUtilities/inc/HepTransform.hh"

namespace mu2e {

  class TrkAlignPlane : public DbTable {
  public:

    typedef std::shared_ptr<TrkAlignPlane> ptr_t;
    typedef std::shared_ptr<const TrkAlignPlane> cptr_t;

    class Row {
    public:
      Row(int index, float dx, float dy, float dz, 
	  float rx, float ry, float rz):
	_index(index),_dx(dx),_dy(dy),_dz(dz),
	_rx(rx),_ry(ry),_rz(rz),_transform(dx,dy,dz,rx,ry,rz) {}
      int index() const { return _index; }
      float dx() const {return _dx;}
      float dy() const {return _dy;}
      float dz() const {return _dz;}
      float rx() const {return _rx;}
      float ry() const {return _ry;}
      float rz() const {return _rz;}
      HepTransform const& transform() const {return _transform;}
    private:
      int _index;
      float _dx;
      float _dy;
      float _dz;
      float _rx;
      float _ry;
      float _rz;
      HepTransform _transform;
    };


    TrkAlignPlane():DbTable("TrkAlignPlane","trk.alignplane",
	 "index,dx,dy,dz,rx,ry,rz") {}
    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    virtual std::size_t nrowFix() const { return 36; };
    size_t size() const { return _csv.capacity() + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),
			 std::stof(columns[1]),
			 std::stof(columns[2]),
			 std::stof(columns[3]),
			 std::stof(columns[4]),
			 std::stof(columns[5]),
			 std::stof(columns[6])  );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.index()<<",";
      sstream << std::fixed << std::setprecision(4); 
      sstream << r.dx()<<",";
      sstream << r.dy()<<",";
      sstream << r.dz()<<",";
      sstream << std::fixed << std::setprecision(6); 
      sstream << r.rx()<<",";
      sstream << r.ry()<<",";
      sstream << r.rz()<<",";
    }

    virtual void clear() { _csv.clear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
