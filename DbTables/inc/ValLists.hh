#ifndef DbTables_ValLists_hh
#define DbTables_ValLists_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class ValLists : public DbTable {
 public:
  class Row {
   public:
    Row(int lid, std::string name, std::string comment, std::string create_time,
        std::string create_user) :
        _lid(lid),
        _name(name), _comment(comment), _create_time(create_time),
        _create_user(create_user) {}
    int lid() const { return _lid; }
    std::string const& name() const { return _name; }
    std::string const& comment() const { return _comment; }
    std::string const& create_time() const { return _create_time; }
    std::string const& create_user() const { return _create_user; }

   private:
    int _lid;
    std::string _name;
    std::string _comment;
    std::string _create_time;
    std::string _create_user;
  };

  constexpr static const char* cxname = "ValLists";

  ValLists() :
      DbTable(cxname, "val.lists", "lid,name,comment,create_time,create_user") {
  }

  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  const Row& row(const int lid) const { return _rows.at(_index.at(lid)); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  size_t size() const override {
    size_t b = baseSize() + sizeof(this) + nrow() * 44 + nrow() * nrow() / 2;
    for (auto const& r : _rows)
      b += r.name().capacity() + r.comment().capacity();
    return b;
  };
  const std::string orderBy() const { return std::string("lid"); }

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoi(columns[0]), columns[1], columns[2],
                       columns[3], columns[4]);
    _index[_rows.back().lid()] = _rows.size() - 1;
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.lid() << ",";
    sstream << r.name() << ",";
    sstream << r.comment() << ",";
    sstream << r.create_time() << ",";
    sstream << r.create_user();
  }

  virtual void clear() override {
    baseClear();
    _rows.clear();
  }

 private:
  std::vector<Row> _rows;
  std::map<int, std::size_t> _index;
};

}  // namespace mu2e
#endif
