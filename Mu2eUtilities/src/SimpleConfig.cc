/*
 *
 * Main class in a primitive runtime parameter utility.
 *
 * $Id: SimpleConfig.cc,v 1.2 2010/04/22 18:30:01 kutschke Exp $
 * $Author: kutschke $ 
 * $Date: 2010/04/22 18:30:01 $
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

// Framework includes
#include "FWCore/Utilities/interface/EDMException.h"

// Mu2e includes
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/TrimInPlace.hh"
#include "Mu2eUtilities/src/SimpleConfigRecord.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {
  //
  // Work list:
  // 1) The following record will fail to parse correctly.
  //    string name = "//This is not a comment";
  // 2) Check presumed eof conditions to make sure that they are eof, not errors.

  /**
   * The only constructor.
   *
   * @return The singleton instance of this class.
   */
  SimpleConfig::SimpleConfig( const string& filename ):
    inputfile(filename)
  {
    ReadFile();
  }

  /**
   * Return a vector<string> containing the names of all variables found in
   * the input file.
   *
   * @return a vector<string> containing all variable names.
   */
  void SimpleConfig::getNames(vector<string>& V) const{
    Rmap_type::const_iterator b = rmap.begin();
    Rmap_type::const_iterator e = rmap.end();
    for ( ; b!=e; ++b){
      V.push_back(b->first);
    }
  }

  bool SimpleConfig::hasName( const string& name ) const{
    return rmap.find(name) != rmap.end();
  }


  /**
   * Return the name of the input file.
   *
   * @return name of the input file.
   */
  string SimpleConfig::inputFile() const{
    return inputfile;
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
    if ( v.size() != nRequired){
      throw edm::Exception(edm::errors::Unknown)
        << "SimpleConfig: Wrong number of elements in vector<string> "
        << name 
        << " in file " 
        << inputfile 
        << " Required: "
        << nRequired
        << " Found: "
        << v.size();
    }
  }
  
  /**
   * Get a specified parameter as a vector<int>.
   *
   * @return the value of the parameter as a vector<int>.
   */
  void SimpleConfig::getVectorInt ( const string& name, vector<int>& v, int nRequired ) const{
    getRecord(name).getVectorInt(v);
    if ( nRequired < 0 ) return;
    if ( v.size() != nRequired){
      throw edm::Exception(edm::errors::Unknown)
        << "SimpleConfig: Wrong number of elements in vector<int> "
        << name 
        << " in file " 
        << inputfile 
        << " Required: "
        << nRequired
        << " Found: "
        << v.size();
    }
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
    if ( v.size() != nRequired){
      throw edm::Exception(edm::errors::Unknown)
        << "SimpleConfig: Wrong number of elements in vector<String> "
        << name 
        << " in file " 
        << inputfile 
        << " Required: "
        << nRequired
        << " Found: "
        << v.size();
    }

  }

  Hep3Vector SimpleConfig::getHep3Vector ( const std::string& name ) const{
    vector<double> tmp;
    getVectorDouble(name,tmp,3);
    Hep3Vector val( tmp[0], tmp[1], tmp[2]);
    return val;
  }

  Hep3Vector SimpleConfig::getHep3Vector ( const std::string& name,
                                           const Hep3Vector& def ){
    if ( hasName(name) ) {
      vector<double> tmp;
      getVectorDouble(name,tmp,3);
      Hep3Vector val( tmp[0], tmp[1], tmp[2]);
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
   * Print the image of the input file to the specfied stream.
   *
   * @return
   */
  void SimpleConfig::print( std::ostream& ost ) const{

    Image_type::const_iterator b = image.begin();
    Image_type::const_iterator e = image.end();
    for ( vector<SimpleConfigRecord>::size_type i=0;
          i<image.size(); ++i ){
      if ( !image[i]->isCommentOrBlank() ) {
        image[i]->print(ost);
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
    Rmap_type::const_iterator b = rmap.find(name);
    Rmap_type::const_iterator e = rmap.end();
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
        << inputfile;

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

    ifstream in(inputfile.c_str());
    if ( !in ) {
    
      // No conf file for this test.
      throw edm::Exception(edm::errors::Unknown)
        << "SimpleConfig: Cannot open input file: "
        << inputfile
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
      
      // Add extension lines if needed.
      if ( WantsExtension(line) ){
      
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
        image.push_back(Record_sptr(new SimpleConfigRecord(line)));
      }
    
    }

    // Form the map after the image has been made.

    // Loop over all records in the image.
    Image_type::const_iterator b0 = image.begin();
    Image_type::const_iterator e0 = image.end();
    for ( ; b0 != e0; ++b0 ){

      // For readability.
      const SimpleConfigRecord& r  = **b0;

      // Add only interesting records to the map.
      if ( !r.isCommentOrBlank() ) {

        const string& key = r.getName();

        // Check for duplicate parameter names.
        Rmap_type::const_iterator b = rmap.find(key);
        Rmap_type::const_iterator e = rmap.end();
        if ( b != e ){

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

        // All Ok. Add it to the map.
        rmap[key] = *b0;
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
    int icomment = s.find("//");
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
   * Does this line contain an include?
   * Right now demands #include start in column 0.
   * Should we weaken this to first non-blank?
   */
  bool SimpleConfig::hasInclude( const std::string& line ){
    string::size_type idx = line.find("#include");
    return idx==0;
  }


  /**
   *  Insert contents of an included file.
   *  
   */
  void SimpleConfig::processInclude( const std::string& line){
  }

} // end namespace mu2e
