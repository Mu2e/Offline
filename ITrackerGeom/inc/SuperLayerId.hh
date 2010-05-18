#ifndef SUPERLAYERID_HH
#define SUPERLAYERID_HH

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

#endif /*SUPERLAYERID_HH*/
