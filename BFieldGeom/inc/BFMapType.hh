#ifndef BFieldGeom_BFMapType_hh
#define BFieldGeom_BFMapType_hh
//
// An enum-matched-to-names class for magnetic field types.
//
//
// Original author Rob Kutschke
//
// Notes:
// 1) I put the enum_type in the class.  If I do not, and if I
//    repeat this model in another enum class, then there will
//    be ambiguities about lastEnum.
// 2) I do want both the name() operator and the << operator.
//    If only <<, then users will need to make temp ostringsteams
//    if they want to do their own formatting.
// 3) There are some notes on alternate implementations
//    at the end of MCDataProducts/inc/GenId.cc file.
// 4) Root stores enum types as 32 bit ints.
// 5) I think it is impossible to construct an invalid object.
//    So there is no isValid() member function.
//

#include <iostream>
#include <string>

namespace mu2e {

    class BFMapType {
       public:
        enum enum_type { unknown, GMC, G4BL, PARAM, lastEnum };

        // Must keep this in sync with the enum above.  Used in BFMapType.cc
        // There is no lastEnum in this list.
#define BFMAPTYPE_NAMES "unknown", "GMC", "G4BL", "PARAM"

       public:
        // The most important c'tor and accessor methods are first.
        explicit BFMapType(enum_type id) : _id(id) {}

        // Id code.
        enum_type id() const { return _id; }

        // A synonym for id().
        operator enum_type() const { return _id; }

        // Name corresponding to this id.
        std::string const& name() const { return name(_id); }

        // ROOT requires a default c'tor.
        BFMapType() : _id(unknown){};

        // This version checks the range of its argument.
        explicit BFMapType(int id);

        // Accept the compiler generator d'tor and copy c'tor.

        // Assignment from an enum_type.
        BFMapType& operator=(BFMapType::enum_type rhs) {
            _id = rhs;
            return *this;
        }

        bool operator==(BFMapType const g) const { return (_id == g._id); }

        bool operator==(BFMapType::enum_type const g) const { return (_id == g); }

        // Accessor for the version.
        static int version() { return _version; }

        // Static version of the name method.
        static std::string const& name(enum_type id);

        // Static version of the name method with an int argument.
        static std::string const& name(int id) {
            // Range safety enforced inside this c'tor.
            BFMapType bft(id);
            return bft.name();
        }

        // Check validity of an Id.  Goal posts are both invalid.
        static bool isValid(enum_type id) { return (id > unknown && id < lastEnum); }

        static void printAll(std::ostream& ost);

       private:
        // The one and only per-instance member datum.
        enum_type _id;

        // Can this make sense?  What happens if I read in two different
        // files that have different versions?  Should I use cvs version instead?
        // This is really an edm question not a question for the class itself.
        static const unsigned _version = 1000;
    };

    // Shift left (printing) operator.
    inline std::ostream& operator<<(std::ostream& ost, const BFMapType& id) {
        ost << "( " << id.id() << ": " << id.name() << " )";
        return ost;
    }

}  // end namespace mu2e

#endif /* BFieldGeom_BFMapType_hh */
