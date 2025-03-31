#ifndef ConfigTools_SimpleConfig_hh
#define ConfigTools_SimpleConfig_hh
//
// Main class in a primitive runtime parameter utility.
//
//
// Contact person Rob Kutschke
//
// Parses the input file into records, holds a collection of
// records and provides access by parameter name.  Provides
// error checking of the global structure of the file and
// for attempts to access non-existent parameters.
//
// Limitations:
// 1) The following record will be incorrectly indentified as containing
//    a comment.  It will then fail to parse: the closing " and ; are missing.
//     string name = "//This is not a comment";
// 2) Does not escape new lines within a string properly.
// 3) The system is designed to read an input file and thereafter become forever const.
//    The internal pointers can become invalidated if someone adds a method to
//    add new parameters.
//

// C++ includes
#include <cstdlib>
#include <map>
#include <ostream>
#include <string>
#include <vector>

// Boost includes
#include <boost/shared_ptr.hpp>

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes

namespace mu2e {

  // Forward declaration.
  class SimpleConfigRecord;

  class SimpleConfig{

  public:

    /**
     * Constructor taking the control parameters as arguments.
     *
     * filename             - The name of the input file.  It is given as a path relative one of the
     *                        roots found in MU2E_SEARCH_PATH.
     *
     * allowReplacement     - If true, a later instance of a parameter overrides the
     *                        all previous instances ( a last one wins policy ).
     *                      - If false, it is an error for a parameter to be repeated
     *                        and the code will throw an exception.
     *
     * messageOnReplacement - Print a message when a parameter is replaced with an
     *                        later instance.
     *
     * messageOnDefault     - Print a message when a parameter is not found in the file
     *                        and takes on a default value.
     *
     */
    SimpleConfig( const std::string& filename = "runtime.conf",
                  bool allowReplacement       = true,
                  bool messageOnReplacement   = false,
                  bool messageOnDefault       = false );
    ~SimpleConfig() = default;

    // This class is not copyable. See note 3.
    SimpleConfig( const SimpleConfig& ) = delete;
    SimpleConfig( SimpleConfig&& )      = delete;
    SimpleConfig& operator=(SimpleConfig const&) = delete;
    SimpleConfig& operator=(SimpleConfig&&     ) = delete;


  public:

    /**
     * Fill a vector<string> with the names of all variable names found in
     * the input file.
     *
     * @return void
     */
    void getNames( std::vector<std::string>& V) const;

    /**
     * Return true if the input file contains the named parameter.
     *
     * @return true if the input file contains the named parameter.
     */
    bool hasName( const std::string& name ) const;

    // Accessors to named parameters, separated by data type.
    // All checking is done in the SimpleConfigRecord class.

    /**
     * Get a specified parameter as a string.  Works for all record types.
     *
     * @return the value of the parameter.
     */
    std::string getString ( const std::string& name ) const;

    /**
     * Get a specified parameter as a string, if not present in the file
     * return the value specified by the second argument.
     *
     * @return the value of the parameter as a string.
     */
    std::string getString ( const std::string& name,
                            const std::string& def ) const;

    /**
     * Get a specified parameter as a int.
     *
     * @return the value of the parameter as an int.
     */
    int getInt ( const std::string& name ) const;

    /**
     * Get a specified parameter as a int, if not present in the file
     * return the value specified by the second argument.
     *
     * @return the value of the parameter as an int.
     */
    int getInt ( const std::string& name, int def ) const;

    /**
     * Get a specified parameter as a double.
     *
     * @return the value of the parameter as an double.
     */
    double getDouble ( const std::string& name ) const;


    /**
     * Get a specified parameter as a double, if not present in the file
     * return the value specified by the second argument.
     *
     * @return the value of the parameter as an double.
     */
    double getDouble ( const std::string& name, double def ) const;


    /**
     * Get a specified parameter as a bool.
     *
     * @return the value of the parameter as an bool.
     */
    bool getBool ( const std::string& name ) const;

    /**
     * Get a specified parameter as a bool, if not present in the file
     * return the value specified by the second argument.
     *
     * @return the value of the parameter as an bool.
     */
    bool getBool ( const std::string& name, bool def ) const;

    /**
     * Get a specified parameter as a vector<std::string>.
     * Works for all parameter types.
     *
     * @return the value of the parameter as a vector<std::string>.
     */
    void getVectorString ( const std::string& name,
                           std::vector<std::string>&,
                           int nRequired=-1) const;

    /**
     * Get a specified parameter as a vector<std::string>.
     * If the parameter is absent, return the default value.
     * Works for all parameter types.
     *
     * @return the value of the parameter as a vector<std::string>.
     */
    void getVectorString ( const std::string&              name,
                           std::vector<std::string>&       v,
                           const std::vector<std::string>& vdefault,
                           int                             nRequired=-1) const;

    /**
     * Get a specified parameter as a vector<int>.
     *
     * @return the value of the parameter as a vector<int>.
     */
    void getVectorInt ( const std::string& name,
                        std::vector<int>& v,
                        int nRequired=-1) const;

    /**
     * Get a specified parameter as a vector<int>.
     * If the parameter is absent, return the default value.
     *
     * @return the value of the parameter as a vector<int>.
     */
    void getVectorInt ( const std::string& name,
                        std::vector<int>& v,
                        const std::vector<int>& vdefault,
                        int nRequired=-1) const;

    /**
     * Get a specified parameter as a vector<double>.
     *
     * @return the value of the parameter as a vector<double>.
     */
    void getVectorDouble ( const std::string& name,
                           std::vector<double>& v,
                           int nRequired=-1) const;

    /**
     * Get a specified parameter as a vector<double>.
     * If the parameter is absent, return the default value.
     *
     * @return the value of the parameter as a vector<double>.
     */
    void getVectorDouble ( const std::string& name,
                           std::vector<double>& v,
                           const std::vector<double>& vdefault,
                           int nRequired=-1) const;

    /**
     * Get a specified parameter as a CLHEP::Hep3Vector.
     * It is specified in the input files as vector<double> of length 3.
     *
     * @return the value of the parameter as CLHEP::Hep3Vector.
     */
    CLHEP::Hep3Vector getHep3Vector ( const std::string& name ) const;

    /**
     * Get a specified parameter as a CLHEP::Hep3Vector.
     * It is specified in the input files as vector<double> of length 3.
     *
     * @return the value of the parameter as a CLHEP::Hep3Vector.
     */
    CLHEP::Hep3Vector getHep3Vector ( const std::string& name,
                                      const CLHEP::Hep3Vector& def ) const;

    /**
     * Get a specified parameter as a CLHEP::Hep2Vector.
     * It is specified in the input files as vector<double> of length 2.
     *
     * @return the value of the parameter as CLHEP::Hep2Vector.
     */
    CLHEP::Hep2Vector getHep2Vector ( const std::string& name ) const;

    /**
     * Get a specified parameter as a CLHEP::Hep2Vector.
     * It is specified in the input files as vector<double> of length 2.
     *
     * @return the value of the parameter as a CLHEP::Hep2Vector.
     */
    CLHEP::Hep2Vector getHep2Vector ( const std::string& name,
                                      const CLHEP::Hep2Vector& def ) const;

    /**
     * Return the requested record as a formatted string.
     *
     * @return a formatted copy of the requested parameter.
     */
    std::string toString ( const std::string& name ) const;

    /**
     * Print the config information to the specfied stream.
     * Print in input order, supressing comments, blank
     * lines and superceded items.
     *
     * @return
     */
    void print( std::ostream& ost, std::string tag = "") const;

    /**
     * Print the config information to std::cout.
     * Print in input order, supressing comments, blank
     * lines and superceded items.
     *
     * @return
     */
    void print() const{
      print(std::cout);
    }

    /**
     * Print the config information to the specfied stream.
     * Print in input order, do not supressing superceded items.
     *
     * @return
     */
    void printFullImage( std::ostream& ost) const;


    // print a summary of what was read
    void printOpen( std::ostream& ost, std::string tag="SimpleConfig" ) const;

    /**
     * Print statistics about the configuration information.
     *
     * @return
     */

    void printStatisticsByType ( std::ostream& ost, std::string tag="" ) const;

    /**
     * Return the name of the input file.
     *
     * @return name of the input file.
     */
    std::string inputFile() const;

    // hash of input file contents, available after construction
    std::size_t inputFileHash() const {return _inputFileHash;}

    // count of input file lines, available after construction
    std::size_t inputFileLines() const {return _inputFileLines;}

    /**
     * Print access counts for each record.
     *
     * @return
     */
    void printAccessCounts( std::ostream& ost, std::string const& tag) const;

    /**
     * Print the names of all of the records that have not yet been accessed.
     *
     * @return
     */
    void printNeverAccessed( std::ostream& out, std::string const& tag = "" ) const;

    /**
     * Print the names of all of the records that have been accessed more than once.
     *
     * @return
     */
    void printAccessedMultiple( std::ostream& ost, std::string const& tag) const;

    /**
     * Print the names of all of the records that were not present in the file
     * and for which the default values were used.
     *
     * @return
     */
    void printDefaultCounts( std::ostream& ost, std::string const& tag) const;

    /**
     * Print all of the summaries; this produces redundant information.
     *
     * @return
     */
    void printAllSummaries( std::ostream& ost, int verbosity, std::string const& tag) const;

    // A type used for some internal record keeping.
    // In gcc 13.1.0 and c++20 it has to be public so that it is visible to
    // the function template, markDefault, defined in SimpleConfig.cc .
    // The alternative is to friend it to all of the instantiations of that template.
    typedef std::map<std::string, int> DefaultCounter_type;

    // Private instance data.

  private:

    // input file string
    std::string _inputFileString;
    // Nname of the input file after lookup
    std::string _inputfile;
    // hash of input file computed as it is read
    std::size_t _inputFileHash;
    // count of input lines
    std::size_t _inputFileLines;

    // If a parameter is repeated in the input file, one two things can happen:
    // true  - the latest value over writes the previous value.
    // false - it is defined to be an error so throw.
    bool _allowReplacement;

    // If true, give a message when a parameter is replaced.
    bool _messageOnReplacement;

    // If true, give a message when an accessor returns a default value.
    bool _messageOnDefault;

    typedef boost::shared_ptr<SimpleConfigRecord> Record_sptr;
    typedef std::map<std::string, Record_sptr > Rmap_type;

    // Access to the input data, keyed by parameter name.
    // Records are owned by the member datum image so rmap does not need to
    // worry about deleting records.
    Rmap_type _rmap;

    // Image of the input file, with comments and empty lines stripped.
    // Records that span multiple lines in the input file are concatenated
    // into a single line in this file.
    typedef std::vector<Record_sptr> Image_type;
    Image_type _image;

    // Record the names of records that were not found in the file and
    // for which the default value was returned, count the number of times
    // each was requested.
    mutable DefaultCounter_type _defaultCounter;

    // Private methods.

    /**
     * Find a record by its name.  This is private so that users cannot
     * directly access a record - an anticipated usage pattern is that
     * people will instantiate a SimpleConfig object, extract information
     * and then let it go out of scope.  We do not want them to be able
     * hold a reference to an internal object after the object has gone
     * out of scope.
     *
     * @return a formatted copy of the record for the requested parameter.
     */
    const SimpleConfigRecord& getRecord ( const std::string& name ) const;

    /**
     *
     * Return true (false) if the record with the requested name does
     * (does not) exist. If the record exists, then fill the shared pointer
     * passed as the second argument.
     *
     * @return a pointer to the requested record, or 0 if there is no such record.
     */
    bool getSharedPointer( const std::string& name,
                           boost::shared_ptr<SimpleConfigRecord>& ) const;

    /**
     * Read the input file and break it into records.
     * Keep a copy of the input file in the orginal record order.
     * Create a map to access records by parameter name.
     *
     */
    void ReadFile();

    /**
     * Test to see if a record is complete.
     *
     * A valid end of line indicator is a semi-colon as the last non-blank character,
     * excluding trailing comments.  Otherwise this record needs an extension.
     *
     * @return true if this record is incomplete and false if it is complete.
     */
    bool WantsExtension( std::string s);

    /**
     * Remove, comments, trailing white space and leading whitespace input string.
     *  - a comment begins with // and continues for the rest of the line.
     *
     * This will give a wrong result on a line like:
     * string name = "//This is not supposed to be a comment";
     *
     * @return a copy of the input with comments and insignificant whitespace removed.
     */
    std::string StripComment( std::string s);

    /**
     * Does this line contain an include?
     */
    bool hasInclude( const std::string& line );

    /**
     *  Insert contents of an included file.
     *
     */
    void processInclude( const std::string& line);

  };  // end class SimpleConfig

  // Function to allow printing using operator<<.
  inline std::ostream& operator<<(std::ostream& ost, const SimpleConfig& c){
    c.print(ost);
    return ost;
  }

} // end namespace mu2e
#endif /* ConfigTools_SimpleConfig_hh */
