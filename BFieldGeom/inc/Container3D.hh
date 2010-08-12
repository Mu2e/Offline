#ifndef Container3D_hh
#define Container3D_hh

//
// A templated class to hold a collection of objects defined on a
// 3D grid.
//
// $Id: Container3D.hh,v 1.2 2010/08/12 16:29:27 genser Exp $
// $Author: genser $
// $Date: 2010/08/12 16:29:27 $
//

#include <vector>
#include <stdexcept>
#include <sstream>

namespace mu2e {

  template <typename OBJ>
  class Container3D
  {
  public:

    // Default constructor.
    // Need this for ROOT IO.  But normally should never use it.
    Container3D():
      _nx(),
      _ny(),
      _nz(),
      _vec(){
    }

    // Normal constructor.
    Container3D( unsigned int nx, unsigned int ny, unsigned int nz):
      _nx(nx),
      _ny(ny),
      _nz(nz),
      _vec(_nx*_ny*_nz,OBJ()){
    }

    // Copy c'tor.
    Container3D(const Container3D& rhs):
      _nx(rhs._nx),
      _ny(rhs._ny),
      _nz(rhs._nz),
      _vec(rhs._vec){
    }

    // Assignment operator.
    Container3D& operator=(const Container3D& rhs){
      _nx = rhs._nx;
      _ny = rhs._ny;
      _nz = rhs._nz;
      _vec  = rhs._vec;
    }

    ~Container3D() { }

    // Set element, without safety features.  Use if the calller has 
    // already ensured the validity of the arguments.
    void set(unsigned int ix, unsigned int iy, unsigned int iz, OBJ const& obj ){
      _vec[index(ix,iy,iz)] = obj;
    }

    // Get element, without safety features. Use if the calller has 
    // already ensured the validity of the arguments.
    OBJ const& get( unsigned int ix, unsigned int iy, unsigned int iz) const {
      return _vec[index(ix,iy,iz)];
    }

    // Synonym for get, without safety features.
    OBJ const& operator()( unsigned int ix, unsigned int iy, unsigned int iz) const {
      return _vec[index(ix,iy,iz)];
    }

    // Set, with safety features.
    void setSafe(unsigned int ix, unsigned int iy, unsigned int iz, OBJ const& obj ){
      isValidOrThrow(ix,iy,iz);
      _vec.at(index(ix,iy,iz)) = obj;
    }

    // Get, with safety features.
    OBJ const& getSafe( unsigned int ix, unsigned int iy, unsigned int iz) const {
      isValidOrThrow(ix,iy,iz);
      return _vec.at(index(ix,iy,iz));
    }

    // Check for a valid index
    bool isValid(int ix, int iy, int iz){
      if ( ix > -1 && ix < _nx &&
           iy > -1 && iy < _ny &&
           iz > -1 && iz < _nz ) return true;
      return false;
    }

    // Throw if the point is not valid.
    void isValidOrThrow( unsigned int ix, unsigned int iy, unsigned int iz){
      if ( !isValid(ix, iy, iz) ){
        std::ostringstream os;
        os << "Invalid index into Container3D: " 
           << ix << " "
           << iy << " "
           << iz << " | Limits are: "
           << _nx << " "
           << _ny << " "
           << _nz;          
        throw std::domain_error(os.str());
      }
    }

  private:

    // Dimensions of the grid.
    unsigned int _nx, _ny, _nz;

    // Container to hold everything.
    std::vector<OBJ> _vec;

    // Compute the index into the array.
    typename std::vector<OBJ>::size_type index(unsigned int ix, unsigned int iy, unsigned int iz) const {
      return ix*_ny*_nz + iy*_nz + iz;
    }

  };
}

#endif
