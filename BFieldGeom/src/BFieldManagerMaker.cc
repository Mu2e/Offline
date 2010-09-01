//
// Build a BFieldManager.
//
// $Id: BFieldManagerMaker.cc,v 1.3 2010/09/01 20:29:02 genser Exp $
// $Author: genser $ 
// $Date: 2010/09/01 20:29:02 $
//

// Includes from C++
#include <iostream>
#include <set>

// Includes from C ( needed for block IO ).
#include <fcntl.h>
#include <sys/stat.h>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Includes from Mu2e
#include "BFieldGeom/inc/BFieldManagerMaker.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "BFieldGeom/inc/DiskRecord.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeneralUtilities/inc/MinMax.hh"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  BFieldManagerMaker::BFieldManagerMaker( const SimpleConfig& config ):
    _config(config){

    // Instantiate an empty BFieldManager.
    _bfmgr = auto_ptr<BFieldManager>(new BFieldManager() );

    // Add the DS, TS and PS field maps.
    loadGMC( "DS", "bfield.dsFile", "bfield.dsDimensions" );
    loadGMC( "TS", "bfield.tsFile", "bfield.tsDimensions" );
    loadGMC( "PS", "bfield.psFile", "bfield.psDimensions" );

    // Special case for the uniform DS field.
    loadUniformDS();
  }

  BFieldManagerMaker::~BFieldManagerMaker(){
  }

  // Parse the config file to learn about one magnetic field map.
  // Create an empty map and call the code to load the map from the file.
  void BFieldManagerMaker::loadGMC( const std::string& key,  
                                    const std::string& fileKey,
                                    const std::string& dimensionKey ){

    if ( ! _config.hasName(fileKey) ){
      cout << "No magnetic field file specified for: "
           << fileKey
           << "   Hope that's OK." << endl;
      return;
    }

    // Get filename and expected dimensions.
    string filename = _config.getString(fileKey);
    vector<int> dim;
    _config.getVectorInt(dimensionKey,dim, 3);

    // Create an empty map.
    BFMap& dsmap = _bfmgr->addBFMap( key,
                                     CLHEP::Hep3Vector(),
                                     dim[0],
                                     dim[1],
                                     dim[2] );

    // Fill the map.
    readGMCMap( filename, dsmap );

  }

  //
  // Read one magnetic field map file in MECO GMC format.
  //
  // This does a 2 pass operation"
  // 1) Pass 1:
  //      Read the input file into a temporary image in memory.
  //      Find the min and max values of the grid points.
  //      A the end of this pass, compute the grid spacing.
  // 2) Pass 2:
  //      Fill the 3D array from the image in memory.
  //

  void BFieldManagerMaker::readGMCMap( const string& filename,
                                       BFMap& bfmap ){

    // Open the input file.
    int fd = open( filename.c_str(), O_RDONLY );
    if ( !fd ) {
      throw cms::Exception("GEOM")
        << "Could not open file containing the magnetic filed map for: "
        << bfmap.getKey() << "\n"
        << "Filename: " 
        << filename
        << "\n";
    }

    // Compute number of records in the input file.
    const int nrecords = computeArraySize(fd,filename);

    // Image of the file in memory.
    vector<DiskRecord> data(nrecords, DiskRecord());

    // Read file into memory.
    const int nbytes = nrecords*sizeof(DiskRecord);
    ssize_t s = read( fd, &data[0], nbytes );
    if ( s != nbytes ) {
      if ( s == -1 ){
        throw cms::Exception("GEOM")
          << "Error reading magnetic field map: " 
          << bfmap.getKey() << "\n"
          << "Filename: " 
          << filename
          << "\n";
      } else{
        throw cms::Exception("GEOM")
          << "Wrong number of bytes read from magnetic field map: " 
          << bfmap.getKey() << "\n"
          << "Filename: " 
          << filename
          << "\n";
      }
    }

    // Tool to find min and max values of grid points.
    MinMax mmX, mmY, mmZ;

    // Offset needed to put this map into the Mu2e coordinate system.
    // ( Origin at center of TS ).
    const CLHEP::Hep3Vector& offset = bfmap._origin;

    // Collect distinct values of (X,Y,Z) on the grid points.
    set<float> X, Y, Z;

    // Multiply by this factor to convert from kilogauss to tesla.
    double ratio = CLHEP::kilogauss/CLHEP::tesla;

    // For the image in memory:
    //   1) Transform into the correct set of units.
    //   2) Find min/max of each dimension.
    //   3) Collect unique values of (X,Y,Z) of the grid points.
    for ( vector<DiskRecord>::iterator i=data.begin();
          i != data.end(); ++i ){

      // Modify in place.
      DiskRecord& r = *i;

      // Unit conversion: from (cm, kG) to (mm,T).
      r.x  *= CLHEP::cm; 
      r.y  *= CLHEP::cm; 
      r.z  *= CLHEP::cm;
      r.bx *= ratio; 
      r.by *= ratio; 
      r.bz *= ratio; 

      // Re-centering. - not sure if we really want this here?
      r.x += offset.x(); 
      r.y += offset.y(); 
      r.z += offset.z();
 
      // The one check I can do.
      if ( r.head != r.tail ){
        throw cms::Exception("GEOM")
          << "Error reading magnetic field map.  "
          << "Mismatched head and tail byte counts at record: " << data.size() << "\n"
          << "Could not open file containing the magnetic filed map for: "
          << bfmap.getKey() << "\n"
          << "Filename: " 
          << filename
          << "\n";
      }

      // Update min/max information.
      mmX.compare(r.x);
      mmY.compare(r.y);
      mmZ.compare(r.z);

      // Populate the set of all unique grid values.
      X.insert(r.x);
      Y.insert(r.y);
      Z.insert(r.z);
    }

    // Expected grid dimentsions.
    const int nx = bfmap._nx;
    const int ny = bfmap._ny;
    const int nz = bfmap._nz;

    // Cross-check that the grid read from the file has the size we expected.
    // This is not really a perfect check since there could be round off error
    // in the computation of the grid points.  But the MECO GMC files were written 
    // in a way that this test works.
    if ( X.size() != nx ||
         Y.size() != ny ||
         Z.size() != nz     ){
      throw cms::Exception("GEOM")
        << "Mismatch in expected and observed number of grid points for BField map: " 
        << bfmap.getKey() << "\n"
        << "From file: " 
        << filename
        << "\n"
        << "Expected/Observed x: " << nx << " " << X.size() << "\n"
        << "Expected/Observed y: " << ny << " " << Y.size() << "\n"
        << "Expected/Observed z: " << nz << " " << Z.size() << "\n";
    }

    // Cross-check that we did not find more data than we have room for.
    if ( data.size() > nx*ny*nz-1){
      throw cms::Exception("GEOM")
        << "Too many values read into the field map: " 
        << bfmap.getKey() << "\n"
        << "From file: " 
        << filename
        << "\n"
        << "Expected/Observed size: " << nx*ny*nz << " " << data.size() << "\n";
    }

    // Tell the magnetic field map what its limits are.
    bfmap.setLimits( mmX.min(), mmX.max(),
                     mmY.min(), mmY.max(),
                     mmZ.min(), mmZ.max() );

    // Store grid points and field values into 3D arrays
    for (vector<DiskRecord>::const_iterator i=data.begin(), e=data.end();
         i != e; ++i){

      DiskRecord const& r(*i);
      
      // Find indices corresponding to this grid point.
      // By construction the indices must be in bounds ( we set the limits above ).
      std::size_t ix = bfmap.iX(r.x);
      std::size_t iy = bfmap.iY(r.y);
      std::size_t iz = bfmap.iZ(r.z);
      
      // Store the information into the 3d arrays.
      bfmap._grid.set (ix, iy, iz, CLHEP::Hep3Vector(r.x,r.y,r.z));
      bfmap._field.set(ix, iy, iz, CLHEP::Hep3Vector(r.bx,r.by,r.bz));
      bfmap._isDefined.set(ix, iy, iz, true);

    }

    return;

  }

  // Compute the size of the array needed to hold the raw data of the field map.
  int BFieldManagerMaker::computeArraySize( int fd, const string& filename ){

    // Get the file size, in bytes, ( info.st_size ).
    struct stat info;
    int stat = fstat( fd, &info);

    // Check that an integral number of records fits in the file.
    int remainder = info.st_size % sizeof(DiskRecord);
    if ( remainder != 0 ){
      throw cms::Exception("GEOM")
        << "Field map file does not hold an integral number of records: \n" 
        << "Filename:  " << filename << "\n"
        << "Size:      " << info.st_size << "\n"
        << "Remainder: " << remainder << "\n";
    }

    // Compute number of records.
    int nrecords  = info.st_size / sizeof(DiskRecord);
    return nrecords;
  }

  void BFieldManagerMaker::loadUniformDS(){
    double bz = _config.getDouble("toyDS.bz", 0.) * CLHEP::tesla;
    _bfmgr->_dsUniformValue = CLHEP::Hep3Vector( 0., 0., bz);
  }

} // end namespace mu2e
