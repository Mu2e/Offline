#ifndef DbTables_CalSimCrystals_hh
#define DbTables_CalSimCrystals_hh

#include <string>
#include <iomanip>
#include <sstream>
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CrystalId.hh"

namespace mu2e {

  class CalSimCrystals : public DbTable {
    public:
      typedef std::shared_ptr<CalSimCrystals> ptr_t;
      typedef std::shared_ptr<const CalSimCrystals> cptr_t;

      class Row {
        public:
          Row(CrystalId crid, float LRU, float pePerMeV_0, float pePerMeV_1):
            _crid(crid), _LRU(LRU), _pePerMeV_0(pePerMeV_0),_pePerMeV_1(pePerMeV_1) {}

          CrystalId  crid()       const {return _crid;}
          float      LRU()        const {return _LRU;}
          float      pePerMeV_0() const {return _pePerMeV_0;}
          float      pePerMeV_1() const {return _pePerMeV_1;}

        private:
          CrystalId _crid;
          float     _LRU;
          float     _pePerMeV_0;
          float     _pePerMeV_1;
      };

      constexpr static const char* cxname = "CalSimCrystal";
      CalSimCrystals():DbTable(cxname,"cal.simcrystal","crid,LRU,pePerMeV0,pePerMeV1"){}

      const Row& row(CrystalId crid) const {return _rows.at(crid.id()); }

      const std::vector<Row>& rows() const {return _rows;}
      std::size_t nrow() const override { return _rows.size(); };
      size_t size() const override { return baseSize()  + nrow()*sizeof(Row); };
      virtual std::size_t nrowFix() const override { return CaloConst::_nCrystal; };
      const std::string orderBy() const { return std::string("crid"); }
      void addRow(const std::vector<std::string>& columns) override {
        std::uint16_t index = std::stoul(columns[0]);
        // enforce order, so channels can be looked up by index
        if (index != _rows.size() || index >= CaloConst::_nCrystal) {
          throw cet::exception("CALOROCALIB_BAD_INDEX")
          << "CalSimCrystals::addRow found index out of order: "
          <<index<< " != " << _rows.size() <<"\n";

        }
        _rows.emplace_back(CrystalId(index),std::stof(columns[1]),std::stof(columns[2]),std::stof(columns[3]));
      }

      void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
        Row const& r = _rows.at(irow);
        sstream << r.crid() <<",";
        sstream << r.LRU()<<",";
        sstream << r.pePerMeV_0()<<",";
        sstream << r.pePerMeV_1()<<"\n";
      }

      virtual void clear() override { baseClear(); _rows.clear();}

    private:
      std::vector<Row> _rows;
  };
}
#endif
