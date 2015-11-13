#ifndef GeneralUtilities_RootNameTitleHelper_hh
#define GeneralUtilities_RootNameTitleHelper_hh

//
// Helper class to construct a series of similar
// names and titles for root objects.
//
// This class helps to make a series of root objects that
// have names that differ only by a trailing integer, such as:
//
// xxx0001, xxx0002, xxx0003
//
// with corresponding titles that differ by a trailing integer,
//
// "Distribution of xxx for event = 1",
// "Distribution of xxx for event = 2",
//
// and so on, where the numbers in the names match the numbers in
// the titles.
//
//
// The arguments are:
// 1) The base string for forming the name (xxx in the above example).
// 2) The base string for the title: "Distribution of xxx for event ="
//    in the above example.
// 3) The number to be appended to the two base strings.
// 4) If positive, used to set the width of the field holding the
//    number within the name object.  If pad=3, then the names appear as:
//    xxx000, xxx001, xxx010, xxx999.
//    If non-positve, non padding is added,
//    are zero filled.  If zero or negative, no leading zeros are included
//    and the above names will appear as:
//    xxx0, xxx1, xxx10, xxx999.
//

#include <string>

class RootNameTitleHelper {

public:
  RootNameTitleHelper ( std::string const& name_base,
                        std::string const& title_base,
                        int id,
                        int pad=-1 );

  // Compiler generated versions of d'tor, copy c'tor
  // and assignment operator are OK.

  const char* name()  { return  _name.data(); }
  const char* title() { return _title.data(); }

private:
  std::string _name;
  std::string _title;

};

#endif /* GeneralUtilities_RootNameTitleHelper_hh */
