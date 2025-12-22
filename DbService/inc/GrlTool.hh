#ifndef DbService_GrlTool_hh
#define DbService_GrlTool_hh

// This class provides an interface to the goodrun database
// There are two main interfaces: bits and list

#include "Offline/DbService/inc/DbReader.hh"
#include "Offline/DbService/inc/DbSql.hh"
#include "Offline/DbService/inc/GrlBit.hh"
#include "Offline/DbService/inc/GrlEntry.hh"
#include "Offline/DbService/inc/GrlList.hh"
#include "Offline/DbService/inc/GrlWord.hh"
#include "Offline/DbTables/inc/DbId.hh"
#include "Offline/GeneralUtilities/inc/StringVec.hh"
#include <string>

namespace mu2e {

class GrlTool {
 public:
  GrlTool();

  // define a quality word
  int defineWord(const std::string& name, const std::string& desc = "");
  // define bits within the word
  int defineBit(const std::string& name, const std::string& bitname,
                const std::string& desc = "");
  // set quality bits for a word, if bits is empty, then no bits set, all good
  int setBits(const std::string& name, const DbIoV& iov,
              const StringVec& bitnames);

  // load the word definitions and copy into results
  int words();
  // load the bit definitions and copy into results
  int bits(const std::string& name = "");
  // load the bit entries definitions and copy into results
  int entries(const std::string& name = "", const DbIoV& iov = DbIoV());

  // apply a file of bit selections commands to create a GRL
  int genList(const std::string& comFile, const std::string& base = "");

  // load the list headers and copy into results
  int lists();
  // fetch a GRL list
  GrlList list(const std::string& name);  // getResult will be text of IoV's
  int createList(const GrlList& list);
  int extendList(const std::string& name, const std::string& siov);
  int lockList(const std::string& name);

  std::string result() const { return _result; }

  void setVerbose(int verbose);

 private:
  // read the database via query_engine
  DbId _id;
  DbReader _reader;
  DbSql _sql;
  int _verbose;
  GrlWord::WordVec _words;
  GrlBit::BitVec _bits;
  GrlEntry::EntryVec _entries;
  GrlHeader::HeaderVec _headers;
  std::string _result;
};

}  // namespace mu2e

#endif
