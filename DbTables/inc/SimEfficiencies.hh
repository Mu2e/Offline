#ifndef DbTables_SimEfficiencies_hh
#define DbTables_SimEfficiencies_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "cetlib_except/exception.h"
#include "DbTables/inc/DbTable.hh"

namespace mu2e {
  
  class SimEfficiencies : public DbTable {
  public:
    typedef std::shared_ptr<SimEfficiencies> ptr_t;
    typedef std::shared_ptr<const SimEfficiencies> cptr_t;

    class Row {
    public:
      Row(int idx, std::string effname, double eff):
	_idx(idx),_effname(effname),_eff(eff) {}
      int    idx() const { return _idx; }
      
      std::string  effname() const { return _effname;}
      double  eff() const { return _eff; }

    private:
      int _idx;
      std::string _effname;
      double _eff;
    };


    SimEfficiencies():DbTable("SimEfficiencies","sim.efficiencies","idx,effname,eff") { }

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int idx) const { return _rows.at(idx); }
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
	throw cet::exception("SIMEFFICIENCIES_BAD_INDEX") 
	  << "SimEfficiencies::addRow found index out of order: " 
	  <<idx << " != " << _rows.back().idx()+1 <<"\n";
      }

      _rows.emplace_back(idx,
			 columns[1],
			 std::stof(columns[2]));
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.idx()<<",";
      sstream << r.effname()<<",";
      sstream << r.eff();
    }

    void findEff(std::string name, double& eff) const { 
      for (const auto& i_row : _rows) {
	if (i_row.effname() == name) {
	  eff =  i_row.eff();
	}
      }
      throw cet::exception("SIMEFFICIENCIES_BAD_EFFNAME") << "Efficiency with name " << name << " not found in database" << std::endl;
    }
    
    virtual void clear() { _csv.clear(); _rows.clear();}
    
  private:
    std::vector<Row> _rows;
  };
};
#endif
