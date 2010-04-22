#ifndef ExampleExtrasSimpleConfig_HH
#define ExampleExtrasSimpleConfig_HH

/**
 *
 * Main class in a primitive runtime parameter utility.
 *
 * $Id: SimpleConfig.hh,v 1.3 2010/04/22 18:30:01 kutschke Exp $
 * $Author: kutschke $ 
 * $Date: 2010/04/22 18:30:01 $
 *
 * Original author Rob Kutschke
 *
 * Parses the input file into records, holds a collection of
 * records and provides access by parameter name.  Provides
 * error checking of the global structure of the file and
 * for attempts to access non-existent parameters.
 *
 * Limitatiopns:
 * 1) The following record will be incorrectly indentified as containing
 *    a comment.  It will then fail to parse: the closing " and ; are missing.
 *     string name = "//This is not a comment"; 
 * 2) Does not escape new lines within a string properly.
 *
 *@author $Author: kutschke $
 *@version $Id: SimpleConfig.hh,v 1.3 2010/04/22 18:30:01 kutschke Exp $
 *
 * Date $Date%
 *
 */

// C++ includes
#include <ostream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

// Boost includes
#include <boost/shared_ptr.hpp>

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Mu2e includes

namespace mu2e {

  // Forward declaration.
  class SimpleConfigRecord;

  class SimpleConfig{

  public:
    
    /**
     * Constructor.  Takes the name of an input file as an argument.
     *
     * @return The singleton instance of this class.
     */
    SimpleConfig( const std::string& filename = "runtime.conf");
    
    ~SimpleConfig(){}
    
  private:
    // This class is not copyable.  These methods are not implemented.
    // Need to deal with the auto_ptr and bare pointers to make it copyable.
    SimpleConfig( const SimpleConfig& );
    SimpleConfig& operator=(const SimpleConfig&);

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
     * Get a specified parameter as a vector<int>.
     *
     * @return the value of the parameter as a vector<int>.
     */
    void getVectorInt ( const std::string& name, 
                        std::vector<int>&,
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
     * Get a specified parameter as a Hep3Vector.
     *
     * @return the value of the parameter as Hep3Vector.
     */
    CLHEP::Hep3Vector getHep3Vector ( const std::string& name ) const;

    /**
     * Get a specified parameter as a Hep3Vector.
     *
     * @return the value of the parameter as a Hep3Vector.
     */
    CLHEP::Hep3Vector getHep3Vector ( const std::string& name,
                                      const CLHEP::Hep3Vector& def );
  
    /**
     * Return the requested record as a formatted string.
     *
     * @return a formatted copy of the requested parameter.
     */
    std::string toString ( const std::string& name ) const;
    
    /**
     * Print the image of the input file to the specfied stream.
     *
     * @return
     */
    void print( std::ostream& ost) const;
    
    /**
     * Print the image of the input file to std::cout.
     *
     * @return
     */
    void print() const{
      print(std::cout);
    }

    /**
     * Return the name of the input file.
     *
     * @return name of the input file.
     */
    std::string inputFile() const;
    
    // Private instance data.

  private:

    // Name of the input file.
    std::string inputfile;

    typedef boost::shared_ptr<SimpleConfigRecord> Record_sptr;
    typedef std::map<std::string, Record_sptr > Rmap_type;
    
    // Access to the input data, keyed by parameter name.
    // Records are owned by the variable image so rmap does not need to 
    // worry about deleting records.
    Rmap_type rmap;
  
    // Image of the input file, with comments and empty lines stripped.
    // Records that span multiple lines in the input file are concatenated
    // into a single line in this file.
    typedef std::vector<Record_sptr> Image_type;
    Image_type image;  

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


  // Function to allow printing using the shift-in operator.
  inline std::ostream& operator<<(std::ostream& ost, const SimpleConfig& c){
    c.print(ost);
    return ost;
  }

} // end namespace mu2e
#endif
