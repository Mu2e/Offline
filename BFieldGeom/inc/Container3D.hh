#ifndef BFieldGeom_Container3D_hh
#define BFieldGeom_Container3D_hh

//
// A templated class to hold a collection of objects defined on a
// 3D grid.
//
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
      _nx(0u),
      _ny(0u),
      _nz(0u),
      _vec(){
    }

    // Normal constructor.
    Container3D( unsigned int nx, unsigned int ny, unsigned int nz):
      _nx(nx),
      _ny(ny),
      _nz(nz),
      _vec(_nx*_ny*_nz,OBJ()){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    // Set element, without safety features.  Use if the caller has
    // already ensured the validity of the arguments.
    void set(unsigned int ix, unsigned int iy, unsigned int iz, OBJ const& obj ){
      _vec[index(ix,iy,iz)] = obj;
    }

    // Get element, without safety features. Use if the caller has
    // already ensured the validity of the arguments.
    OBJ const& get( unsigned int ix, unsigned int iy, unsigned int iz) const {
      return _vec[index(ix,iy,iz)];
    }

    // Get element, without safety features. Use if the caller has
    // already ensured the validity of the arguments.
    OBJ& get( unsigned int ix, unsigned int iy, unsigned int iz){
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
    bool isValid(unsigned int ix, unsigned int iy, unsigned int iz) const{
      return ix >= 0 && ix < _nx &&
             iy >= 0 && iy < _ny &&
             iz >= 0 && iz < _nz;
    }

    // Throw if the point is not valid.
    void isValidOrThrow( unsigned int ix, unsigned int iy, unsigned int iz) const {
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

    // Recover the memory. T
    void cleart(){
      _nx = 0;
      _ny = 0;
      _nz = 0;
      std::vector<OBJ>().swap(_vec);
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

  template <>
  class Container3D<bool>
  {
  public:

    // Default constructor.
    // Need this for ROOT IO.  But normally should never use it.
    Container3D():
      _nx(0u),
      _ny(0u),
      _nz(0u),
      _vec(){
    }

    // Normal constructor.
    Container3D( unsigned int nx, unsigned int ny, unsigned int nz):
      _nx(nx),
      _ny(ny),
      _nz(nz),
      _vec(_nx*_ny*_nz,false){
    }

    // "Normal" constructor for bool
    Container3D( unsigned int nx, unsigned int ny, unsigned int nz, bool bv):
      _nx(nx),
      _ny(ny),
      _nz(nz),
      _vec(_nx*_ny*_nz,bv){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    // Set element, without safety features.  Use if the caller has
    // already ensured the validity of the arguments.
    void set(unsigned int ix, unsigned int iy, unsigned int iz, bool obj ){
      _vec[index(ix,iy,iz)] = obj;
    }

    // Get element, without safety features. Use if the caller has
    // already ensured the validity of the arguments.
    bool get( unsigned int ix, unsigned int iy, unsigned int iz) const {
      return _vec[index(ix,iy,iz)];
    }

    // This would be needed to exactly parallel the interface of the general case
    // but this is not meaningful for bools.
    //bool& get( unsigned int ix, unsigned int iy, unsigned int iz);

    // Synonym for get, without safety features.
    bool operator()( unsigned int ix, unsigned int iy, unsigned int iz) const {
      return _vec[index(ix,iy,iz)];
    }

    // Set, with safety features.
    void setSafe(unsigned int ix, unsigned int iy, unsigned int iz, bool obj ){
      isValidOrThrow(ix,iy,iz);
      _vec.at(index(ix,iy,iz)) = obj;
    }

    // Get, with safety features.
    bool getSafe( unsigned int ix, unsigned int iy, unsigned int iz) const {
      isValidOrThrow(ix,iy,iz);
      return _vec.at(index(ix,iy,iz));
    }

    // Check for a valid index
    bool isValid(unsigned int ix, unsigned int iy, unsigned int iz) const {
      return ix >= 0 && ix < _nx &&
             iy >= 0 && iy < _ny &&
             iz >= 0 && iz < _nz;
    }

    // Throw if the point is not valid.
    void isValidOrThrow( unsigned int ix, unsigned int iy, unsigned int iz) const {
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
    std::vector<bool> _vec;

    // Compute the index into the array.
    std::vector<bool>::size_type index(unsigned int ix, unsigned int iy, unsigned int iz) const {
      return ix*_ny*_nz + iy*_nz + iz;
    }

  };
}


#endif /* BFieldGeom_Container3D_hh */
