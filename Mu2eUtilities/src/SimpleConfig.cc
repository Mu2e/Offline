/*
 *
 * Main class in a primitive runtime parameter utility.
 *
 * $Id: SimpleConfig.cc,v 1.8 2011/05/02 18:26:51 kutschke Exp $
 * $Author: kutschke $ 
 * $Date: 2011/05/02 18:26:51 $
 *
 * Original author Rob Kutschke
 *
 * Parses a file in the format of:
 * type name = value;
 * vector<type> name = { list, of, values};
 *
 * Provides accessors by parameter name and some error checking.
 * Throws when errors occur.
 *
 */


// C++ includes
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/EDMException.h"

// Mu2e includes
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/TrimInPlace.hh"
#include "Mu2eUtilities/src/SimpleConfigRecord.hh"

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
                              bool messageOnReplacement ):
    _allowReplacement(allowReplacement),
    _messageOnReplacement(messageOnReplacement){

    edm::FileInPath fip(filename);
    _inputfile = fip.fullPath();
    ReadFile();

  }

  SimpleConfig::SimpleConfig( const edm::FileInPath& fileInPath,
                              bool allowReplacement,
                              bool messageOnReplacement):
    _inputfile(fileInPath.fullPath()),
    _allowReplacement(allowReplacement),
    _messageOnReplacement(messageOnReplacement){
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
      throw edm::Exception(edm::errors::Unknown)
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
        throw edm::Exception(edm::errors::Unknown)
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
      throw edm::Exception(edm::errors::Unknown)
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
        throw edm::Exception(edm::errors::Unknown)
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
      throw edm::Exception(edm::errors::Unknown)
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
        throw edm::Exception(edm::errors::Unknown)
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
  void SimpleConfig::print( std::ostream& ost ) const{

    Image_type::const_iterator b = _image.begin();
    Image_type::const_iterator e = _image.end();
    for ( vector<SimpleConfigRecord>::size_type i=0;
          i<_image.size(); ++i ){
      if ( !_image[i]->isCommentOrBlank() && !_image[i]->isSuperceded() ) {
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

    Image_type::const_iterator b = _image.begin();
    Image_type::const_iterator e = _image.end();
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
      throw edm::Exception(edm::errors::Unknown)
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
      throw edm::Exception(edm::errors::Unknown)
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

          nextline = StripComment(nextline);
          if ( WantsExtension(nextline) ){
            if ( nextline.size() != 0 ) {
              all.append(nextline);
              TrimInPlace(all);
            }
            line = all;
          } else{
            all.append(nextline);
            TrimInPlace(all);
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
              edm::LogWarning("CONFIG")
                << "SimpleConfig replacing parameter in: "
                << _inputfile << "\n"
                << "old record: " << *old << "\n"
                << "new record: " << **b0
                << "\n";
            }



          } else{
            // Test: fail03.conf
            throw edm::Exception(edm::errors::Unknown)
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
    TrimInPlace(line);

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
      TrimInPlace(s);
      return s;
    }
  
    string s1(s.substr(0,islash));
    TrimInPlace(s1);
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
      throw edm::Exception(edm::errors::Unknown)
        << "Cannot find first double quote on include line: "
        << line
        << endl;
    }

    // Ensure that all charcaters between #include and opening " are
    // legal whitespace.  This can be triggered by a missing leading quote.
    for ( string::size_type i=idx+8; i<j0; ++i ){
      if ( line[i] != ' ' && line[i] != '\t' ){
        throw edm::Exception(edm::errors::Unknown)
          << "Unexpected characters after include and before first double quote.\n"
          << "Maybe a missing leading quote?\n"
          << line
          << endl;

      }
    }

    // Find the double quote that ends the filename.
    string::size_type j1 = line.find("\"",j0+1);
    if ( j1 == string::npos ){
      throw edm::Exception(edm::errors::Unknown)
        << "Cannot find trailing double quote on include line: \n"
        << line
        << endl;
    }

    // Check for a non-empty file name.
    if ( j1-j0-1 == 0 ){
      throw edm::Exception(edm::errors::Unknown)
        << "Empty filename in include line: \n"
        << line
        << endl;
    }

    // The filename of the file to be included.
    string fname = line.substr(j0+1,j1-j0-1);

    // Read the included file.
    SimpleConfig nestedFile(fname);

    // Copy the contents of the included file into this one.
    for ( Image_type::const_iterator i=nestedFile._image.begin();
          i != nestedFile._image.end(); ++i ){
      if ( !(**i).isCommentOrBlank() ){
        _image.push_back(*i);
      }
    }

  }

  // Some types used in printStatistics.
  struct Stats{
    int count;
    int nValues;
    Stats():
      count(0),
      nValues(0){
    }
  };
  typedef std::map<string,Stats> StatMap;

  // A helper function used in printStatistics.
  // Find the length of the longest key in the map.
  int maxKeySize ( const StatMap& m ){

    size_t maxSize(0);

    for ( StatMap::const_iterator i=m.begin();
          i != m.end(); ++i ){

      const string& key = i->first;
      maxSize = ( key.size() > maxSize ) ? key.size() : maxSize;
    }

    return maxSize;

  } // end maxKeySize


  void SimpleConfig::printStatistics ( std::ostream& ost ){

    // This will hold the accumlated information.
    StatMap statMap;

    // Loop over all records and accumulate stats.
    for ( vector<SimpleConfigRecord>::size_type i=0;
          i<_image.size(); ++i ){

      // Only look at live records.
      if ( !_image[i]->isCommentOrBlank() && 
           !_image[i]->isSuperceded() ) {

        const SimpleConfigRecord& record = *_image[i];
        Stats& stats = statMap[record.getType()];

        // Accumulate stats.
        ++stats.count;
        stats.nValues += record.size();

      }
    }

    // Done computation; start printout.
    ost << "\nStatistics for SimpleConfig file: " << _inputfile << endl;

    if ( statMap.size() == 0 ) {
      ost << "There were no live records in this file." << endl;
      return;
    }

    // Count of records that are "live", not comments, blank or superceded.
    int totalCount(0);

    // Count total values in the file.
    int totalValues(0);

    // Compute the width of field that holds the name of the type.
    int typeFieldWidth(maxKeySize(statMap));

    // Make sure there is room for "Total" (in case we only have bool and int).
    typeFieldWidth = ( typeFieldWidth > 5 ) ? typeFieldWidth : 5;

    // Width of the other fields.
    int countFieldWidth(8);
    int nValuesFieldWidth(8);

    ost << setw(typeFieldWidth)    << "Type"     << " "
        << setw(countFieldWidth)   << "   Count" << " "
        << setw(nValuesFieldWidth) << "  Values"
        << endl;

    // Print the information for each data type and accumulate totals.
    for ( StatMap::const_iterator i=statMap.begin();
            i != statMap.end(); ++i ){

      const string& type = i->first;
      const Stats& stats = i->second;

      totalCount += stats.count;
      totalValues += stats.nValues;

      ost << setw(typeFieldWidth)    <<  type << " "
          << setw(countFieldWidth)   << stats.count << " "
          << setw(nValuesFieldWidth) << stats.nValues
          << endl;
    }

    // Print totals.
    ost << setw(typeFieldWidth)    << "Total"     << " "
        << setw(countFieldWidth)   << totalCount  << " "
        << setw(nValuesFieldWidth) << totalValues
        << endl;

  } // end printStatistics


} // end namespace mu2e
