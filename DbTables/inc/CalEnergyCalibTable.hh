#ifndef DbTables_CalEnergyCalibTable_hh
#define DbTables_CalEnergyCalibTable_hh
/*

FIXME - placeholder for the energy calib reco table
*/

#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"

namespace mu2e {

  class CalEnergyCalibTable : public DbTable {
  public:

    class Row {
    public:
      Row(int roid): _roid(roid) {}
      int       roid() const { return _roid;} 


    private:
      int   _roid; 
    //TODO - decide what goes in the combined class
    };

    constexpr static const char* cxname = "CalEnergyCalibTable";

    CalEnergyCalibTable():DbTable(cxname,"cal.energycalib",
			"roid") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int roid) const { 
                return _rows.at(_chanIndex.at(roid)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*nrow()/2 + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      int roid = std::stoi(columns[0]);
      // TODO - do we need this:
      if(roid!=int(_rows.size())) {
	      throw cet::exception("CALENERGYCALIB_BAD_INDEX") 
	        << "CalEnergyCalibTable::addRow found index out of order: " 
	        <<roid << " != " << _rows.back().roid()+1 <<"\n";
      }
       _rows.emplace_back(roid ); //TODO - add other columns here
      // add this channel to the map index - optional
      _chanIndex[_rows.back().roid()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid(); //TODO add other columns
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;
    std::map<int,std::size_t> _chanIndex;
  };
  
};
#endif
