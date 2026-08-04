#ifndef DbTables_CalAlignElement_hh
#define DbTables_CalAlignElement_hh
//
// Alignment tables for calorimeter elements
//
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DbTables/inc/CalAlignParams.hh"
#include <iomanip>
#include <sstream>
#include <string>

namespace mu2e {

  class CalAlignElement : public DbTable {
    public:
      typedef std::shared_ptr<CalAlignElement>       ptr_t;
      typedef std::shared_ptr<const CalAlignElement> cptr_t;

      CalAlignElement(const char* Name, const char* DbName, size_t nrows) :
          DbTable(Name, DbName, "index,dx,dy,dz,rx,ry,rz"), _nrows(nrows)
      {}

      const CalAlignParams& rowAt(const std::size_t index) const {
        return _rows.at(index);
      }
      const std::vector<CalAlignParams>& rows() const {return _rows;}
      std::size_t nrow()    const override {return _rows.size();};
      std::size_t nrowFix() const override {return _nrows;};
      std::size_t size()    const override {
        return baseSize() + nrow() * sizeof(CalAlignParams);
      };
      const std::string orderBy() const {return std::string("index");}

      void addRow(const std::vector<std::string>& columns) override {
        unsigned index = std::stoul(columns[0]);
        // enforce a strict sequential order
        if (index >= _nrows || index != _rows.size()) {
          throw cet::exception("CALALIGNELEMENT_BAD_INDEX")
              << "CalAlignElement::addRow found index out of order: " << index
              << " != " << _rows.size() << "\n";
        }
        _rows.emplace_back(index,std::stof(columns[1]),std::stof(columns[2]),
                           std::stof(columns[3]),std::stof(columns[4]),
                           std::stof(columns[5]),std::stof(columns[6]));
      }

      void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
        CalAlignParams const& r = _rows.at(irow);
        sstream << r.index() << ",";
        sstream << std::fixed
                << std::setprecision(3);
        sstream << r.dx() << ",";
        sstream << r.dy() << ",";
        sstream << r.dz() << ",";
        sstream << std::fixed << std::setprecision(3);
        sstream << r.rx() << ",";
        sstream << r.ry() << ",";
        sstream << r.rz() << ",";
      }

      void clear() override {
        baseClear();
        _rows.clear();
      }

    private:
      std::vector<CalAlignParams> _rows;
      size_t _nrows;
  };


  class CalAlignDisk : public CalAlignElement {
    public:
      constexpr static const char* cxname = "CalAlignDisk";
      CalAlignDisk() :
          CalAlignElement(cxname, "cal.aligndisk", CaloConst::_nDisk)
      {}
      CalAlignDisk(const char* name, const char* dbname) :
        CalAlignElement(name, dbname,  CaloConst::_nDisk)
      {}
  };

  class CalAlignCrystal : public CalAlignElement {
    public:
      constexpr static const char* cxname = "CalAlignCrystal";
      CalAlignCrystal() :
          CalAlignElement(cxname, "cal.aligncrystal", CaloConst::_nCrystal)
      {}
      CalAlignCrystal(const char* name, const char* dbname) :
        CalAlignElement(name, dbname,  CaloConst::_nCrystal)
      {}
  };


}
#endif
