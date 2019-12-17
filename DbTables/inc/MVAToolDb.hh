#ifndef DbTables_MVAToolDb_hh
#define DbTables_MVAToolDb_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Mu2eUtilities/inc/MVATools.hh"
#include "cetlib_except/exception.h"
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class MVAToolDb : public DbTable {
  public:

    class Row {
    public:
      Row(int idx, std::string mvaname, std::string xmlfilename, int calibrated):
	_idx(idx),_mvaname(mvaname),_xmlfilename(xmlfilename),_calibrated(calibrated) {}
      int    idx() const { return _idx; }
      
      std::string  mvaname() const { return _mvaname;}
      std::string  xmlfilename() const { return _xmlfilename; }
      int calibrated() const { return _calibrated; }

      void setCalib() {
	MVATools* mva = new MVATools(xmlfilename());
	mva->getCalib(_effCalib);
      }

      const std::map<float, float>& effCalib() const { return _effCalib; }

    private:
      int _idx;
      std::string _mvaname;
      std::string _xmlfilename;
      int _calibrated;

      std::map<float, float> _effCalib;
    };


    MVAToolDb():DbTable("MVAToolDb","MVATool.db","idx,mvaname,xmlfilename,calibrated") {
      throw cet::exception("MVAToolDb") << "Shouldn't be creating a bare MVAToolDb table" << std::endl;
    }

    MVAToolDb(std::string mva):DbTable(mva.c_str(),mva.c_str(),"idx,mvaname,xmlfilename,calibrated") { }

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int idx) const { 
                return _rows.at(_chanIndex.at(idx)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    //this table should always be 3 rows
    //    virtual std::size_t nrowFix() const { return 3; }; 
    size_t size() const { return _csv.capacity() + 
	+ nrow()*nrow()/2 + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) {
      int idx = std::stoi(columns[0]);
      // enforce a strict sequential order - optional
      if(idx!=int(_rows.size())) {
	throw cet::exception("MVATOOLDB_BAD_INDEX") 
	  << "MVAToolDb::addRow found index out of order: " 
	  <<idx << " != " << _rows.back().idx()+1 <<"\n";
      }
      int calibrated = std::stoi(columns[3]);
      if (calibrated < 0 || calibrated > 1) {
	throw cet::exception("MVATOOLDB_BAD_CALIBRATED_INT") << "MVAToolDb::addRow found calibrated wasn't 0 or 1 but " << calibrated << std::endl;
      }

      _rows.emplace_back(idx,
			 columns[1],
			 columns[2],
			 calibrated);
      // add this idx to the map index - optional
      _chanIndex[_rows.back().idx()] = _rows.size()-1;

      // calibrate if we can
      if (calibrated == 1) {
	_rows.back().setCalib();
      }
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.idx()<<",";
      sstream << r.mvaname()<<",";
      sstream << r.xmlfilename();
    }


    const std::map<float, float>& getCalib(const int idx) const {
      if (row(idx).calibrated()!=1) {
	throw cet::exception("MVATOOLDB_BAD_GETCALIB") << "MVAToolDb::getCalib tried to get calibration when \"calibrated = " << row(idx).calibrated() << "\"" << std::endl;
      }
      return row(idx).effCalib();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); _chanIndex.clear();}

  private:
    std::vector<Row> _rows;
    std::map<int,std::size_t> _chanIndex;
  };
};
#endif
