#ifndef ITrackerGeom_SuperLayerId_hh
#define ITrackerGeom_SuperLayerId_hh

#include <iostream>

namespace mu2e {

struct SuperLayerId{

public:

  SuperLayerId():
    _id(-1){
  }

  SuperLayerId( int &id ):
    _id(id)
  {
  }

  ~SuperLayerId  (){
  }

  const int getSuperLayer() const{
    return _id;
  }

  bool operator==(const SuperLayerId s) const{
    return ( _id == s._id );
  }

  int _id;

};

inline std::ostream& operator<<(std::ostream& ost,
                                const SuperLayerId& s ){
  ost << "SuperLayer Id: "<<s._id << " )";
  return ost;
}

} //namespace mu2e

#endif /* ITrackerGeom_SuperLayerId_hh */
