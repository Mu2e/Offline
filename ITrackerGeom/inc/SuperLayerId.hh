#ifndef ITrackerGeom_SuperLayerId_hh
#define ITrackerGeom_SuperLayerId_hh

#include <iostream>

namespace mu2e {
  class SuperLayerId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const SuperLayerId& s );

class SuperLayerId{

public:

  SuperLayerId():
    _id(-1){
  }

  explicit SuperLayerId( int id ):
    _id(id)
  {
  }

  // use compiler-generated copy c'tor, copy assignment, and d'tor

  int getSuperLayer() const{
    return _id;
  }

  bool operator==(const SuperLayerId s) const{
    return ( _id == s._id );
  }

  friend std::ostream& operator<<(std::ostream& ost,
                                  const SuperLayerId& s ){
    ost << "SuperLayer Id: "<<s._id << " )";
    return ost;
  }

private:

  int _id;

};

} //namespace mu2e

#endif /* ITrackerGeom_SuperLayerId_hh */
