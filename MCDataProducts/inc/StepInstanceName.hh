#ifndef MCDataProducts_StepInstanceName_hh
#define MCDataProducts_StepInstanceName_hh
//
// An enum-matched-to-names class for the names of the StepPointMC
// collections produced inside G4_module.
//
// $Id: StepInstanceName.hh,v 1.10 2014/06/05 21:06:24 genser Exp $
// $Author: genser $
// $Date: 2014/06/05 21:06:24 $
//
// Contact person Rob Kutschke
//
// Notes:
// 1) The enum_type must be in the class to avoid clashes with
//    similar classes that are modelled on this.  In particular
//    unknown and lastEnum will clash.
//
// 2) I do want both the name() operator and the << operator.
//    If only <<, then users will need to make temp ostringsteams
//    if they want to do their own formatting.
//
// 3) Root always stores enum types as 32 bit ints, regardless of
//    native word length on a machine.
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

  class StepInstanceName {

  public:

    // Need to keep the enum and the following MACRO in sync.
    enum enum_type {
      unknown,
      tracker,        virtualdetector,   timeVD,            stoppingtarget,    CRV, 
      calorimeter,    calorimeterRO,     calorimeterROCard, calorimeterCrate,  ExtMonUCITof,     
      trackerDS,     protonabsorber,    PSVacuum,          stepper,           trackerSWires,     
      itrackerFWires, trackerWalls,      STMDet,            panelEBKey, DSCableRun,
      ProductionTargetCoreSection, ProductionTargetStartingCoreSection, 
      ProductionTargetFinStartingSection,ProductionTargetNegativeEndRing, ProductionTargetPositiveEndRing,
      ProductionTargetFinSection,ProductionTargetFinTopSection,ProductionTargetFinTopStartingSection,lastEnum
    };

    // Keep this in sync with the enum. Used in StepInstanceName.cc
#define STEPINSTANCENAME_NAMES \
      "unknown", \
      "tracker",       "virtualdetector", "timeVD",            "stoppingtarget",   "CRV",            \
      "calorimeter",   "calorimeterRO",   "calorimeterROCard", "calorimeterCrate", "ExtMonUCITof",   \
      "trackerDS",    "protonabsorber",  "PSVacuum",          "stepper",          "trackerSWires",  \
	"itrackerFWires", "trackerWalls",   "STMDet",            "panelEBKey", "DSCableRun", \
      "ProductionTargetCoreSection", "ProductionTargetStartingCoreSection", \
      "ProductionTargetFinStartingSection","ProductionTargetNegativeEndRing","ProductionTargetPositiveEndRing",\
	"ProductionTargetFinSection","ProductionTargetFinTopSection","ProductionTargetFinTopStartingSection"

  public:

    // The most important c'tor and accessor methods are first.
    explicit StepInstanceName( enum_type id):
      _id(id)
    {}

    enum_type id() const { return _id;}

    // Member function version of the name function.
    std::string const& name() const {
      return name(_id);
    }

    // Find the Id corresponding to the given name, or throw.
    // If the second argument is true, also throw if the result is the "undefined" value.
    static StepInstanceName findByName( std::string const& name, bool throwIfUndefined=true );

    // The number of names, including "undefined".
    static size_t size(){
      return names().size();
    }

    // ROOT requires a default c'tor.
    StepInstanceName():
      _id(unknown){
    }

    // Need this to interface with legacy code; prefer to remove it if possible.
    explicit StepInstanceName( int id):
      _id(static_cast<enum_type>(id)){
      isValidorThrow(_id);
    }

    virtual ~StepInstanceName(){}

    // This operator implements:
    //   StepInstanceName a;
    //   enum_type b;
    // a = b;
    StepInstanceName& operator=(StepInstanceName::enum_type const& c){
      _id = c;
      return *this;
    }

    // This operator implements:
    //   StepInstanceName a;
    //   enum_type b = a;
    operator StepInstanceName::enum_type ()const{
      return _id;
    }

    // Comparison operators
    bool operator==(const StepInstanceName g) const{
      return ( _id == g._id );
    }

    bool operator==(const StepInstanceName::enum_type g) const{
      return ( _id == g );
    }

    // Static version of the name method.
    static std::string const& name( enum_type id );

    // Get all names, as strings.
    static std::vector<std::string> const& names();

    // Get all values, as class instances.
    static std::vector<StepInstanceName> const& allValues();

    // Check validity of an Id.
    static bool isValid( enum_type id){
      if ( id <   unknown ) return false;
      if ( id >= lastEnum ) return false;
      return true;
    }

    // Check validity and throw if invalid.
    static bool isValidorThrow( enum_type id);

    // Print all known values and their names.
    static void printAll( std::ostream& ost = std::cout);

    // Member function version of static functions.
    bool isValid() const{
      return isValid(_id);
    }

  private:

    // The one and only per-instance member datum.
    enum_type _id;

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const StepInstanceName& id ){
    ost << "( "
        << id.id() << ": "
        << id.name()
        << " )";
    return ost;
  }

}

#endif /* MCDataProducts_StepInstanceName_hh */
