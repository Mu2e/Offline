//
// Main class in a primitive runtime parameter utility.
//
// $Id: SimpleConfig.cc,v 1.3 2012/07/27 19:39:39 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/27 19:39:39 $
//
// Contact person Rob Kutschke
//
// Parses a file in the format of:
// type name = value;
// vector<type> name = { list, of, values};
//
// Provides accessors by parameter name and some error checking.
// Throws when errors occur.
//

// C++ includes
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <boost/functional/hash.hpp>

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/src/SimpleConfigRecord.hh"
#include "GeneralUtilities/inc/trimInPlace.hh"

using namespace std;

namespace mu2e {
  //
  // Work list:
  // 1) The following record will fail to parse correctly.
  //    string name = "//This is not a comment";
  // 2) Check presumed eof conditions to make sure that they are eof, not errors.

  /**
   * The constructors.
   *
   */
  SimpleConfig::SimpleConfig( const string& filename,
                              bool allowReplacement,
                              bool messageOnReplacement,
                              bool messageOnDefault):
    _inputFileString(filename),
    _inputFileHash(0),
    _inputFileLines(0),
    _allowReplacement(allowReplacement),
    _messageOnReplacement(messageOnReplacement),
    _messageOnDefault(messageOnDefault)     {

    ConfigFileLookupPolicy configFile;
    _inputfile = configFile(filename);
    ReadFile();
  }

  /**
   * Return a vector<string> containing the names of all variables found in
   * the input file.
   *
   * @return a vector<string> containing all variable names.
   */
  void SimpleConfig::getNames(vector<string>& V) const{
    Rmap_type::const_iterator b = _rmap.begin();
    Rmap_type::const_iterator e = _rmap.end();
    for ( ; b!=e; ++b){
      V.push_back(b->first);
    }
  }

  bool SimpleConfig::hasName( const string& name ) const{
    return _rmap.find(name) != _rmap.end();
  }


  /**
   * Return the name of the input file.
   *
   * @return name of the input file.
   */
  string SimpleConfig::inputFile() const{
    return _inputfile;
  }

  // Accessors to named parameters, separated by data type.

  /**
   * Get a specified parameter as a string.  Works for all record types.
   *
   * @return the value of the parameter.
   */
  string SimpleConfig::getString ( const string& name ) const{
    return getRecord(name).getString();
  }

  /**
   * A helper used by markDefault.
   *
   * @return the value of the parameter.
   */
  template <class T>
  inline std::ostream& operator<<(std::ostream& ost,
                                  typename std::vector<T> const& v ){
    for ( typename std::vector<T>::const_iterator i=v.begin(), e=v.end();
          i != e; ++i ){
      ost << *i;
    }
    return ost;
  }

  /**
   * Record the name of records for which the default was taken; keep a count
   * of how often each name was asked for.
   *
   * @return the value of the parameter.
   */
  template <class T>
  void markDefault ( std::string const& rtype,
                     std::string const& name,
                     T const& t,
                     SimpleConfig::DefaultCounter_type& counter,
                     bool print ){
    int& count(counter[name]);
    ++count;
    if ( print ){
      ostringstream os;
      os << t;
      mf::LogPrint("GEOM")
        << "SimpleConfig: default value used: ("
        << setw(5) << count << ") "
        << rtype   << " "
        << name    << " "
        << os.str()
        << endl;
    }
  }

  /**
   * Get a specified parameter as a string, if not present in the file
   * return the value specified by the second argument.
   *
   * @return the value of the parameter as a string.
   */
  string SimpleConfig::getString ( const string& name,
                                   const string& def ) const {
    Record_sptr b;
    if ( getSharedPointer(name,b) ){
      return b->getString();
    }
    markDefault ( "string ", name, def, _defaultCounter, _messageOnDefault );
    return def;
  }

  /**
   * Get a specified parameter as a int.
   *
   * @return the value of the parameter as an int.
   */
  int SimpleConfig::getInt ( const string& name ) const{
    return getRecord(name).getInt();
  }

  /**
   * Get a specified parameter as a int, if not present in the file
   * return the value specified by the second argument.
   *
   * @return the value of the parameter as an int.
   */
  int SimpleConfig::getInt ( const string& name, int def )const{
    Record_sptr b;
    if ( getSharedPointer(name,b) ){
      return b->getInt();
    }
    markDefault ( "int ", name, def, _defaultCounter, _messageOnDefault );
    return def;
  }

  /**
   * Get a specified parameter as a double.
   *
   * @return the value of the parameter as an double.
   */
  double SimpleConfig::getDouble ( const string& name ) const{
    return getRecord(name).getDouble();
  }

  /**
   * Get a specified parameter as a double, if not present in the file
   * return the value specified by the second argument.
   *
   * @return the value of the parameter as an double.
   */
  double SimpleConfig::getDouble ( const string& name, double def ) const{
    Record_sptr b;
    if ( getSharedPointer(name,b) ){
      return b->getDouble();
    }
    markDefault ( "double ", name, def, _defaultCounter, _messageOnDefault );
    return def;
  }


  /**
   * Get a specified parameter as a bool.
   *
   * @return the value of the parameter as an bool.
   */
  bool SimpleConfig::getBool ( const string& name ) const{
    return getRecord(name).getBool();
  }

  /**
   * Get a specified parameter as a bool, if not present in the file
   * return the value specified by the second argument.
   *
   * @return the value of the parameter as an bool.
   */
  bool SimpleConfig::getBool ( const string& name, bool def ) const{
    Record_sptr b;
    if ( getSharedPointer(name,b) ){
      return b->getBool();
    }
    markDefault ( "bool ", name, def, _defaultCounter, _messageOnDefault );
    return def;
  }


  /**
   * Get a specified parameter as a vector<string>. Works for all parameter types.
   *
   * @return the value of the parameter as a vector<string>.
   */
  void SimpleConfig::getVectorString ( const string& name, vector<string>& v, int nRequired ) const{
    getRecord(name).getVectorString(v);
    if ( nRequired < 0 ) return;
    if ( v.size() != static_cast<size_t>(nRequired) ){
      throw cet::exception("SimpleConfig")
        << "SimpleConfig: Wrong number of elements in vector<string> "
        << name
        << " in file "
        << _inputfile
        << " Required: "
        << nRequired
        << " Found: "
        << v.size();
    }
  }

  /**
   * Get a specified parameter as a vector<string>.
   * If the parameter is absent, return the default value.
   *
   * @return the value of the parameter as a vector<string>.
   */
  void SimpleConfig::getVectorString ( const string&         name,
                                       vector<string>&       v,
                                       const vector<string>& vdefault,
                                       int                   nRequired ) const{

    // If value is present, extract it from the configuration.
    if ( hasName(name) ){
      getVectorString( name, v, nRequired);
      return;
    }

    // If asked to, check that the default value has the right length.
    if ( nRequired > -1 ) {
      if ( vdefault.size() != static_cast<size_t>(nRequired) ){
        throw cet::exception("SimpleConfig")
          << "SimpleConfig: Wrong number of elements in vector<string> "
          << name
          << " in file "
          << _inputfile
          << " Required: "
          << nRequired
          << " Found: "
          << vdefault.size();
      }
    }

    // Assign the default value;
    v = vdefault;

    markDefault ( "vector<string> ", name, v, _defaultCounter, _messageOnDefault );

  }


  /**
   * Get a specified parameter as a vector<int>.
   *
   * @return the value of the parameter as a vector<int>.
   */
  void SimpleConfig::getVectorInt ( const string& name, vector<int>& v, int nRequired ) const{
    getRecord(name).getVectorInt(v);
    if ( nRequired < 0 ) return;
    if ( v.size() != static_cast<size_t>(nRequired) ){
      throw cet::exception("SimpleConfig")
        << "SimpleConfig: Wrong number of elements in vector<int> "
        << name
        << " in file "
        << _inputfile
        << " Required: "
        << nRequired
        << " Found: "
        << v.size();
    }
  }

  /**
   * Get a specified parameter as a vector<int>.
   * If the parameter is absent, return the default value.
   *
   * @return the value of the parameter as a vector<int>.
   */
  void SimpleConfig::getVectorInt ( const string&      name,
                                    vector<int>&       v,
                                    const vector<int>& vdefault,
                                    int                nRequired ) const{

    // If value is present, extract it from the configuration.
    if ( hasName(name) ){
      getVectorInt( name, v, nRequired);
      return;
    }

    // If asked to, check that the default value has the right length.
    if ( nRequired > -1 ) {
      if ( vdefault.size() != static_cast<size_t>(nRequired) ){
        throw cet::exception("SimpleConfig")
          << "SimpleConfig: Wrong number of elements in vector<int> "
          << name
          << " in file "
          << _inputfile
          << " Required: "
          << nRequired
          << " Found: "
          << vdefault.size();
      }
    }

    // Assign the default value;
    v = vdefault;

    markDefault ( "vector<int> ", name, v, _defaultCounter, _messageOnDefault );
  }



  /**
   * Get a specified parameter as a vector<double>.
   *
   * @return the value of the parameter as a vector<double>.
   */
  void SimpleConfig::getVectorDouble ( const string& name,
                                       vector<double>& v,
                                       int nRequired ) const{
    getRecord(name).getVectorDouble(v);
    if ( nRequired < 0 ) return;
    if ( v.size() != static_cast<size_t>(nRequired) ){
      throw cet::exception("SimpleConfig")
        << "SimpleConfig: Wrong number of elements in vector<double> "
        << name
        << " in file "
        << _inputfile
        << " Required: "
        << nRequired
        << " Found: "
        << v.size();
    }

  }

  /**
   * Get a specified parameter as a vector<double>.
   * If the parameter is absent, return the default value.
   *
   * @return the value of the parameter as a vector<double>.
   */
  void SimpleConfig::getVectorDouble ( const string&         name,
                                       vector<double>&       v,
                                       const vector<double>& vdefault,
                                       int                   nRequired ) const{

    // If value is present, extract it from the configuration.
    if ( hasName(name) ){
      getVectorDouble( name, v, nRequired);
      return;
    }

    // If asked to, check that the default value has the right length.
    if ( nRequired > -1 ) {
      if ( vdefault.size() != static_cast<size_t>(nRequired) ){
        throw cet::exception("SimpleConfig")
          << "SimpleConfig: Wrong number of elements in vector<double> "
          << name
          << " in file "
          << _inputfile
          << " Required: "
          << nRequired
          << " Found: "
          << vdefault.size();
      }
    }

    // Assign the default value;
    v = vdefault;

    markDefault ( "vector<double> ", name, v, _defaultCounter, _messageOnDefault );
  }

  CLHEP::Hep3Vector SimpleConfig::getHep3Vector ( const std::string& name ) const{
    vector<double> tmp;
    getVectorDouble(name,tmp,3);
    CLHEP::Hep3Vector val( tmp[0], tmp[1], tmp[2]);
    return val;
  }

  CLHEP::Hep3Vector SimpleConfig::getHep3Vector ( const std::string& name,
                                                  const CLHEP::Hep3Vector& def ) const{
    if ( hasName(name) ) {
      vector<double> tmp;
      getVectorDouble(name,tmp,3);
      CLHEP::Hep3Vector val( tmp[0], tmp[1], tmp[2]);
      return val;
    }
    return def;
  }

  CLHEP::Hep2Vector SimpleConfig::getHep2Vector ( const std::string& name ) const{
    vector<double> tmp;
    getVectorDouble(name,tmp,2);
    CLHEP::Hep2Vector val( tmp[0], tmp[1]);
    return val;
  }

  CLHEP::Hep2Vector SimpleConfig::getHep2Vector ( const std::string& name,
                                                  const CLHEP::Hep2Vector& def ) const{
    if ( hasName(name) ) {
      vector<double> tmp;
      getVectorDouble(name,tmp,2);
      CLHEP::Hep2Vector val( tmp[0], tmp[1]);
      return val;
    }
    return def;
  }


  /**
   * Return the record as a formatted string.
   *
   * @return a formatted copy of the requested parameter.
   */
  string SimpleConfig::toString ( const string& name )const{
    return getRecord(name).toString();
  }


  /**
   * Print the config information to the specfied stream.
   * Print in input order, supressing comments, blank
   * lines and superceded items.
   *
   * @return
   */
  void SimpleConfig::print( std::ostream& ost, std::string tag ) const{

    for ( vector<SimpleConfigRecord>::size_type i=0;
          i<_image.size(); ++i ){
      if ( !_image[i]->isCommentOrBlank() && !_image[i]->isSuperceded() ) {
        ost << tag;
        _image[i]->print(ost);
        ost << endl;
      }
    }
  }

  /**
   * Print the config information to the specfied stream.
   * Print in input order, supressing comments and blank
   * lines.  But include superceded items.
   *
   * @return
   */
  void SimpleConfig::printFullImage( std::ostream& ost ) const{

    for ( vector<SimpleConfigRecord>::size_type i=0;
          i<_image.size(); ++i ){
      if ( !_image[i]->isCommentOrBlank() ) {
        _image[i]->print(ost);
        ost << endl;
      }
    }
  }

  // Private methods.

  /**
   *
   * Return true (false) if the record with the requested name does (does not) exist.
   * If the record exists, then fill the shared pointer passed as the second argument.
   *
   * @return a pointer to the requested record, or 0 if there is no such record.
   */
  bool SimpleConfig::getSharedPointer( const string& name,
                                       Record_sptr& ptr ) const{
    Rmap_type::const_iterator b = _rmap.find(name);
    Rmap_type::const_iterator e = _rmap.end();
    if (b == e) {
      return false;
    }
    ptr=b->second;
    return true;
  }

  /**
   * Return a reference to the record.  Throw if there is no such record.
   *
   * @return a formatted copy of the record for the requested parameter.
   */
  const SimpleConfigRecord& SimpleConfig::getRecord ( const string& name ) const {
    Record_sptr b;

    if ( !getSharedPointer(name,b) ) {

      // Test: fail21.conf
      throw cet::exception("SimpleConfig")
        << "SimpleConfig: No such parameter "
        << name
        << " in file "
        << _inputfile;

    }
    return *b;
  }




  /**
   * Read the input file, break it into records.
   * Keep a copy of the input file in the orginal record order.
   * Create a map to access records by parameter name.
   *
   */
  void SimpleConfig::ReadFile(){

    ifstream in(_inputfile.c_str());
    if ( !in ) {

      // No conf file for this test.
      throw cet::exception("SimpleConfig")
        << "SimpleConfig: Cannot open input file: "
        << _inputfile
        << endl;
    }

    string line;
    while ( in ){

      getline(in,line);
      if ( !in ){
        // Normal end of file.
        // ?? Check eofbit to be sure.  Could also be failbit or badbit.
        break;
      }

      // collect contribution to the hash
      boost::hash_combine<string>(_inputFileHash,line);
      _inputFileLines++;

      // Remove comments.
      line = StripComment(line);

      // If line is not semi-colon terminated, look for the line
      // to be continued.  Exception for a line starting with #include.
      if ( WantsExtension(line) && !hasInclude(line) ){

        string all(line);
        string nextline;
        while ( in ){
          getline(in,nextline);
          if ( !in ){
            // Normal end of file.
            // ?? Check eofbit to be sure.  Could also be failbit or badbit.
            break;
          }

	  // collect contribution to the hash
	  boost::hash_combine<string>(_inputFileHash,line);
	  _inputFileLines++;

          nextline = StripComment(nextline);
          if ( WantsExtension(nextline) ){
            if ( nextline.size() != 0 ) {
              all.append(nextline);
              trimInPlace(all);
            }
            line = all;
          } else{
            all.append(nextline);
            trimInPlace(all);
            line = all;
            break;
          }
        }
      }

      if ( hasInclude(line) ){
        processInclude(line);
      }
      else{
        // Convert to a record and store an image of all records that we find.
        _image.push_back(Record_sptr(new SimpleConfigRecord(line)));
      }

    }

    // Form the map after the image has been made.

    // Loop over all records in the image.
    Image_type::const_iterator b0 = _image.begin();
    Image_type::const_iterator e0 = _image.end();
    for ( ; b0 != e0; ++b0 ){

      // For readability.
      const SimpleConfigRecord& r  = **b0;

      // Add only interesting records to the map.
      if ( !r.isCommentOrBlank() ) {

        const string& key = r.getName();

        // Check for duplicate parameter names.
        Rmap_type::const_iterator b = _rmap.find(key);
        if ( b != _rmap.end() ){

          // Duplicate name.  Either replace previous value or declare an error.
          if ( _allowReplacement ){
            Record_sptr old = _rmap[key];
            old->setSuperceded();
            _rmap[key] = *b0;

            if ( _messageOnReplacement ){

              // Replace with mf::LogPrint when we learn how to
              // disable exponential backoff.
              cerr << "SimpleConfig replacing parameter in: "
                   << _inputfile << "\n"
                   << "old record: " << *old << "\n"
                   << "new record: " << **b0
                   << endl;
            }



          } else{
            // Test: fail03.conf
            throw cet::exception("SimpleConfig")
              << "Duplicate parameter name found in the input file: "
              << endl
              << "This record:     "
              << r.toString()
              << endl
              << "Previous record: "
              << b->second->toString()
              << endl;
          }

        } else{

          // Not a duplicate parameter. Add it to the map.
          _rmap[key] = *b0;
        }

      }

    }
  }

  /**
   * Test to see if this record is complete.
   *
   * A valid end of line indicator is a semi-colon as the last non-blank character
   * before any comments.  Otherwise this record needs an extension.
   *
   * @return true if this record is incomplete and false if it is complete.
   */
  bool SimpleConfig::WantsExtension( string s){
    string::size_type icomment = s.find("//");
    string line = ( icomment == string::npos ) ? s : s.substr(0,icomment);
    trimInPlace(line);

    if ( line.size() == 0 ) return false;
    // The above might be a logic bug; this should not change the state of the
    // the previous line.

    return ( line[line.size()-1] != ';' );
  }

  /**
   * Remove, comments, trailing white space and leading whitespace input string.
   *  - a comment begins with // and continues for the rest of the line.
   *
   * This will give a wrong result on a line like:
   * string name = "//This is not supposed to be a comment";
   *
   * @return a copy of the input with comments and insignificant whitespace removed.
   */
  string SimpleConfig::StripComment( string s){

    // Find comment delimiter if present.
    int islash = s.find("//");

    if ( islash < 0 ){
      trimInPlace(s);
      return s;
    }

    string s1(s.substr(0,islash));
    trimInPlace(s1);
    return s1;

  }

  /**
   * Does this line start with an include?
   * By the time we use this, any leading white space has been stripped.
   */
  bool SimpleConfig::hasInclude( const std::string& line ){
    string::size_type idx = line.find("#include");
    return idx == 0;
  }


  /**
   *  Insert contents of an included file.
   *  The required syntax is:
   *  #include "filename"
   *  where there may be any amount of white space between or after the
   *  two tokens. The filename must be delimited with double quotes.
   *
   */
  void SimpleConfig::processInclude( const std::string& line){
    string::size_type idx = line.find("#include");

    // Find the double quote that opens filename.
    string::size_type j0 = line.find("\"",idx+8);
    if ( j0 == string::npos ){
      throw cet::exception("SimpleConfig")
        << "Cannot find first double quote on include line: "
        << line
        << endl;
    }

    // Ensure that all charcaters between #include and opening " are
    // legal whitespace.  This can be triggered by a missing leading quote.
    for ( string::size_type i=idx+8; i<j0; ++i ){
      if ( line[i] != ' ' && line[i] != '\t' ){
        throw cet::exception("SimpleConfig")
          << "Unexpected characters after include and before first double quote.\n"
          << "Maybe a missing leading quote?\n"
          << line
          << endl;

      }
    }

    // Find the double quote that ends the filename.
    string::size_type j1 = line.find("\"",j0+1);
    if ( j1 == string::npos ){
      throw cet::exception("SimpleConfig")
        << "Cannot find trailing double quote on include line: \n"
        << line
        << endl;
    }

    // Check for a non-empty file name.
    if ( j1-j0-1 == 0 ){
      throw cet::exception("SimpleConfig")
        << "Empty filename in include line: \n"
        << line
        << endl;
    }

    // The filename of the file to be included.
    string fname = line.substr(j0+1,j1-j0-1);

    // Read the included file; copy verbosity control from this object.
    SimpleConfig nestedFile(fname,_allowReplacement,_messageOnReplacement,_messageOnDefault);
    // collect contribution to the hash
    boost::hash_combine<std::size_t>(_inputFileHash,nestedFile.inputFileHash());
    _inputFileLines += nestedFile.inputFileLines();

    // Copy the contents of the included file into this one.
    for ( Image_type::const_iterator i=nestedFile._image.begin();
          i != nestedFile._image.end(); ++i ){
      if ( !(**i).isCommentOrBlank() ){
        _image.push_back(*i);
      }
    }

  }

  // Some types used in printStatistics.
  struct RecordTypeStats{
    int count;           // Number of records of this type.
    int nValues;         // Number of distinct values of this type
                         //   vector<Type> objects are counted by .size() of their values.
    int accessedCount;   // As above but only if the record has been accessed.
    int accessednValues; // As above but only if the record has been accessed.
    RecordTypeStats():
      count(0),
      nValues(0),
      accessedCount(0),
      accessednValues(0){
    }
  };

  typedef std::map<string,RecordTypeStats> StatsByType;

  // A helper function used in printStatistics.
  // Find the length of the longest key in the map.
  int maxKeySize ( const StatsByType& m ){

    size_t maxSize(0);

    for ( StatsByType::const_iterator i=m.begin();
          i != m.end(); ++i ){

      const string& key = i->first;
      maxSize = ( key.size() > maxSize ) ? key.size() : maxSize;
    }

    return maxSize;

  } // end maxKeySize


  void SimpleConfig::printOpen( std::ostream& ost, std::string tag ) const{
    std::cout << tag << " file: "<< _inputfile << std::endl;
    std::cout << tag << " lines: "<< _inputFileLines 
	      <<"  hash: " << _inputFileHash << std::endl;
  }

  void SimpleConfig::printStatisticsByType ( std::ostream& ost, std::string tag ) const{

    // This will hold the accumlated information.
    StatsByType byType;

    // Loop over all records and accumulate stats.
    for ( vector<SimpleConfigRecord>::size_type i=0;
          i<_image.size(); ++i ){

      // Only look at live records.
      if ( !_image[i]->isCommentOrBlank() &&
           !_image[i]->isSuperceded() ) {

        const SimpleConfigRecord& record = *_image[i];
        RecordTypeStats& stats = byType[record.getType()];

        // Accumulate stats.
        ++stats.count;
        stats.nValues += record.size();

        if ( _image[i]->accessCount() > 0 ){
          ++stats.accessedCount;
          stats.accessednValues += record.size();
        }

      }
    }

    // Done computation; start printout.
    ost << "\nStatistics for SimpleConfig file: " << _inputfile << endl;

    if ( byType.size() == 0 ) {
      ost << "There were no live records in this file." << endl;
      return;
    }

    // Count of records that are "live", not comments, blank or superceded.
    int totalCount(0);

    // Count total values in the file.
    int totalValues(0);

    // Repeat for accessed.
    int totalAccessedCount(0);
    int totalAccessedValues(0);

    // Compute the width of field that holds the name of the type.
    int typeFieldWidth(maxKeySize(byType));

    // Make sure there is room for "Total" (in case we only have bool and int).
    typeFieldWidth = ( typeFieldWidth > 5 ) ? typeFieldWidth : 5;

    // Width of the other fields.
    int countFieldWidth(8);
    int nValuesFieldWidth(8);

    ost << tag
        << setw(typeFieldWidth+countFieldWidth+4) << " All"
        << setw(nValuesFieldWidth+countFieldWidth+7) << " Accessed\n";

    ost << tag
        << setw(typeFieldWidth)    << "Type"     << " "
        << setw(countFieldWidth)   << "   Count" << " "
        << setw(nValuesFieldWidth) << "  Values" << " | "
        << setw(countFieldWidth)   << "   Count" << " "
        << setw(nValuesFieldWidth) << "  Values" << " "
        << endl;

    // Print the information for each data type and accumulate totals.
    for ( StatsByType::const_iterator i=byType.begin();
            i != byType.end(); ++i ){

      const string&           type = i->first;
      const RecordTypeStats& stats = i->second;

      totalCount          += stats.count;
      totalValues         += stats.nValues;
      totalAccessedCount  += stats.accessedCount;
      totalAccessedValues += stats.accessednValues;

      ost << tag
          << setw(typeFieldWidth)     <<  type                 << " "
          << setw(countFieldWidth)    << stats.count           << " "
          << setw(nValuesFieldWidth)  << stats.nValues         << "   "
          << setw(countFieldWidth)    << stats.accessedCount   << " "
          << setw(nValuesFieldWidth)  << stats.accessednValues << " "
          << endl;

    }

    // Print totals.

    ost << tag
        << setw(typeFieldWidth)    << "Total"     << " "
        << setw(countFieldWidth)   << totalCount  << " "
        << setw(nValuesFieldWidth) << totalValues << "   "
        << setw(countFieldWidth)   << totalAccessedCount  << " "
        << setw(nValuesFieldWidth) << totalAccessedValues << " "
        << endl;

  } // end printStatisticsByType

  void SimpleConfig::printAccessCounts( ostream& ost, std::string const& tag) const{
    ost << "\nList of access counts:\n"
        << "From file: " << _inputfile << endl;

    int n(0);
    int nValues(0);
    for ( Image_type::const_iterator i=_image.begin(), e=_image.end();
          i != e; ++i ){

      const SimpleConfigRecord& record = **i;
      if ( !record.isCommentOrBlank() &&
           !record.isSuperceded() ) {

        ost << "Access count: "
            << tag
            << setw(4) << record.accessCount() << " "
            << record.getName()
            << endl;
        ++n;
        nValues += record.size();
      }
    }
    ost << tag << "Number of records:                " << n       << endl;
    ost << tag << "Number of values in all records : " << nValues << endl;
  } // end printAccessCounts

  /**
   * Print the names of all of the records that have not yet been accessed.
   *
   * @return
   */
  void SimpleConfig::printNeverAccessed( ostream& ost, std::string const& tag) const{

    ost << "\nList of records that have never been accessed.\n"
        << "From file: " << _inputfile << endl;

    int n(0);
    int nValues(0);
    for ( Image_type::const_iterator i=_image.begin(), e=_image.end();
          i != e; ++i ){

      const SimpleConfigRecord& record = **i;
      if ( !record.isCommentOrBlank() &&
           !record.isSuperceded() ) {

        if ( record.accessCount() < 1 ) {
          ost << "Never accessed: "
              << tag
              << record.getName()
              << endl;
          ++n;
          nValues += record.size();
        }
      }
    }
    ost << "Number of records never accessed:           " << n       << endl;
    ost << "Number of values in records never accessed: " << nValues << endl;

  } // end printNeverAccessed

  /**
   * Print the names of all of the records that have been accessed more than once.
   *
   * @return
   */
  void SimpleConfig::printAccessedMultiple( ostream& ost, std::string const& tag) const{

    ost << "\nList of records that have been accessed more than once.\n"
        << "From file: " << _inputfile << endl;

    const string header("Multiple access: ");
    size_t countFieldWidth(6);
    size_t headerWidth(countFieldWidth+header.size()+tag.size()+7);

    ost << setw(headerWidth) << "Count  Name" << endl;

    int n(0);
    int nValues(0);
    for ( Image_type::const_iterator i=_image.begin(), e=_image.end();
          i != e; ++i ){

      const SimpleConfigRecord& record = **i;
      if ( !record.isCommentOrBlank() &&
           !record.isSuperceded() ) {

        if ( record.accessCount() > 1 ) {
          ost << header
              << tag
              << setw(countFieldWidth) << record.accessCount() << " "
              << record.getName()
              << endl;
          ++n;
          nValues += record.size();
        }
      }
    }
    ost << "Number of record accessed multiple times:           " << n       << endl;
    ost << "Number of values in records accessed multiple time: " << nValues << endl;
  } // end printAccessedMultiple

  /**
   * Print the names of all of the records that were not present in the file
   * and for which the default values were used.
   *
   * @return
   */
  void SimpleConfig::printDefaultCounts( ostream& ost, std::string const& tag) const{

    ost << "\nList of names for which no record was found and for which the default value was used:\n"
        << "From file: " << _inputfile << endl;

    const std::string header("Default used: ");
    size_t countFieldWidth(6);
    size_t headerWidth(countFieldWidth+header.size()+tag.size()+7);

    ost << setw(headerWidth) << "Count  Name" << endl;

    for ( DefaultCounter_type::const_iterator i=_defaultCounter.begin(), e=_defaultCounter.end();
          i != e; ++i ){
      ost << header
          << tag
          << setw(countFieldWidth) << i->second << "  "
          << i->first
          << endl;
    }
  }

  /**
   * Print all of the summaries; this produces redundant information.
   *
   * @return
   */
  void SimpleConfig::printAllSummaries( std::ostream& ost, int verbosity, std::string const& tag) const{
    if ( verbosity <= 0 ) return;

    printStatisticsByType(ost, tag);

    if ( verbosity >= 2 ){
      printAccessCounts ( ost, tag );
    }

    if ( verbosity >= 3 ){
      printNeverAccessed    ( ost, tag );
      printAccessedMultiple ( ost, tag );
      printDefaultCounts    ( ost, tag );
    }

  }


} // end namespace mu2e
