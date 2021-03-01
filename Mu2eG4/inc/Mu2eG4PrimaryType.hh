// Mu2eG4 can create G4 primaries from different Mu2e data types.
// We want to use enums to represent those types in the code, but also
// need a string representation to read them from fhicl config.
//
// An instance of the Mu2eG4PrimaryType class can be created from
// either the enum or a string, and has both enum and string-typed
// accessors.
//
// Implementation features:
//     - thread safe
//     - does not use heap
//     - a single list of entities, no other lists to keep in sync
//
// Andrei Gaponenko, 2021


#ifndef Mu2eG4_inc_Mu2eG4PrimaryInput_hh
#define Mu2eG4_inc_Mu2eG4PrimaryInput_hh

#include <string>
#include <vector>


// Master list of all Mu2eG4 primary input types
//
// We use the X-macro technique to have a single point
// of maintenance for enum values and their names as strings,
// see https://en.wikipedia.org/wiki/X_Macro or (for more
// advanced examples) https://digitalmars.com/articles/b51.html
//
#define MU2EG4_PRIMARY_TYPES                         \
  X( GenParticles )                                  \
  X( StepPoints )                                    \
  X( SimParticleLeaves )                             \
  /**/

namespace mu2e {

  class Mu2eG4PrimaryType {
  public:
   enum enum_type {
#define X(x) x,
      X( unknown)
      MU2EG4_PRIMARY_TYPES
#undef X
    };

    enum_type id() const { return id_; }
    std::string name() const { return name_; }

    explicit Mu2eG4PrimaryType(std::string name);
    Mu2eG4PrimaryType(enum_type id); // allow implicit conversion

    static const std::vector<std::string>& all_names();

  private:
    enum_type id_;
    std::string name_;
  };

} // namespace mu2e

#endif /*Mu2eG4_inc_Mu2eG4PrimaryInput_hh*/
