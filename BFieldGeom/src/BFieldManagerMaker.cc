//
// Build a BFieldManager.
//
// $Id: BFieldManagerMaker.cc,v 1.8 2010/11/24 22:47:26 logash Exp $
// $Author: logash $ 
// $Date: 2010/11/24 22:47:26 $
//

// Includes from C++
#include <iostream>
#include <fstream>
#include <set>

// Includes from C ( needed for block IO ).
#include <fcntl.h>
#include <sys/stat.h>

// Includes from boost
#include <boost/regex.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
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

  //
  // Utility function which adds filter to decompress file if neccessary
  //
  void decompressFile(const string& filename, boost::iostreams::filtering_istream & in) {
    boost::regex re_gz("^.*\\.gz\\s*$");
    boost::regex re_bz("^.*\\.bz2\\s*$");
    if( boost::regex_match(filename.c_str(),re_gz) ) {
      in.push (boost::iostreams::gzip_decompressor());
    } else if( boost::regex_match(filename.c_str(),re_bz) ) {
      in.push (boost::iostreams::bzip2_decompressor());
    }
  }

  BFieldManagerMaker::BFieldManagerMaker( const SimpleConfig& config ):
    _config(config){

    // Instantiate an empty BFieldManager.
    _bfmgr = auto_ptr<BFieldManager>(new BFieldManager() );

    string format = _config.getString("bfield.format","GMC");

    if( format=="GMC" ) {

      // Add the DS, TS and PS field maps.
      loadGMC( "DS", "bfield.dsFile", "bfield.dsDimensions" );
      loadGMC( "TS", "bfield.tsFile", "bfield.tsDimensions" );
      loadGMC( "PS", "bfield.psFile", "bfield.psDimensions" );
      
      //throw cms::Exception("GEOM") << "Temporal end." << "\n";

    } else if( format=="G4BL" ) {

      // These maps require torus radius of 2929 mm
      std::string torusName("toyTS.rTorus");
      const double torusRadius = 2929.0; // Required number
      if( fabs(_config.getDouble(torusName,0.0)-torusRadius)>0.1 ){
	throw cms::Exception("GEOM")
	  << "The G4BL magnetic field files require torus radius of 2929 mm."
	  << " Check " << torusName << " value in the config file." 
	  << " Maps are not loaded.\n";
      }

      // Add the DS, TSu, TSd and PS field maps.
      loadG4BL( "DS",  "bfield.dsFile"  );
      loadG4BL( "PS",  "bfield.psFile"  );
      loadG4BL( "TSu", "bfield.tsuFile" );
      loadG4BL( "TSd", "bfield.tsdFile" );
      
      //throw cms::Exception("GEOM") << "Temporal end." << "\n";

    } else {

      throw cms::Exception("GEOM")
        << "Unknown format of file with magnetic field maps: " << format
        << "\n";

    }

    // For debug purposes: print the field in the target region
    CLHEP::Hep3Vector b = _bfmgr->getBField(CLHEP::Hep3Vector(3900.0,0.0,-6550.0));
    cout << "B-field at the proton target: ("
	 << b.x() << ","
	 << b.y() << ","
	 << b.z() << ")"
	 << endl;

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
                                     dim[0],
                                     dim[1],
                                     dim[2] );

    // Fill the map.
    readGMCMap( filename, dsmap );

    //dsmap.print(std::cout);
  }

  // Parse the config file to learn about one magnetic field map.
  // Create an empty map and call the code to load the map from the file.
  void BFieldManagerMaker::loadG4BL( const std::string& key,  
				     const std::string& fileKey ) {

    if ( ! _config.hasName(fileKey) ){
      cout << "No magnetic field file specified for: "
           << fileKey
           << "   Hope that's OK." << endl;
      return;
    }

    // Get filename
    string filename = _config.getString(fileKey);
    cout << "Read " << filename << endl;

    // Open the input file.
    ifstream fin;
    fin.open(edm::FileInPath(filename).fullPath().c_str());
    if ( !fin.is_open() ) {
      throw cms::Exception("GEOM")
        << "Could not open file containing the magnetic field data. "
        << "Filename: " 
        << filename
        << "\n";
    }
    boost::iostreams::filtering_istream in;
    decompressFile(filename,in);
    in.push(fin);

    // Parse the string with parameters
    char cbuf[128];
    boost::regex re("^\\s*grid"
		    "\\s+X0=([eE\\d\\-\\+\\.]+)\\s+Y0=([eE\\d\\-\\+\\.]+)\\s+Z0=([eE\\d\\-\\+\\.]+)"
		    "\\s+nX=([eE\\d\\-\\+\\.]+)\\s+nY=([eE\\d\\-\\+\\.]+)\\s+nZ=([eE\\d\\-\\+\\.]+)"
		    "\\s+dX=([eE\\d\\-\\+\\.]+)\\s+dY=([eE\\d\\-\\+\\.]+)\\s+dZ=([eE\\d\\-\\+\\.]+)"
		    ".*$");
    boost::cmatch matches;
    bool paramFound = false;
    int nread=100; // Optimization - don't read more than 100 lines
    while( (!in.eof()) && (--nread>0) ) {
      in.getline(cbuf,128);
      if( boost::regex_match(cbuf,matches,re) ) {
	paramFound = true;
	break;
      }
    }
    in.pop();
    fin.close();

    if( ! paramFound ) {
      throw cms::Exception("GEOM")
        << "Could not find param string in magnetic firld map. "
        << "Filename: " << filename
	<< ", found " << matches.size() << " items."
        << "\n";
    }

    vector<double> X0;
    vector<int> dim;
    vector<double> dX;

    for( int i=1; i<=3; i++ ) {
      double value;
      istringstream sin(string(matches[i].first, matches[i].second));
      sin >> value; 
      X0.push_back(value);
    }

    for( int i=4; i<=6; i++ ) {
      int value;
      istringstream sin(string(matches[i].first, matches[i].second));
      sin >> value; 
      dim.push_back(value);
    }

    for( int i=7; i<=9; i++ ) {
      double value;
      istringstream sin(string(matches[i].first, matches[i].second));
      sin >> value; 
      dX.push_back(value);
    }

    // Set the offset
    CLHEP::Hep3Vector G4BL_offset(-3904.0,0.0,7929.0);
    X0[0] = X0[0] - G4BL_offset.x();
    X0[1] = X0[1] - G4BL_offset.y();
    X0[2] = X0[2] - G4BL_offset.z();

    // Create an empty map.
    BFMap& dsmap = _bfmgr->addBFMap( key,
                                     dim[0],
                                     dim[1],
                                     dim[2] );

    // Set defined region for the map
    dsmap.setLimits(X0[0],X0[0]+(dim[0]-1)*dX[0],
		    X0[1],X0[1]+(dim[1]-1)*dX[1],
		    X0[2],X0[2]+(dim[2]-1)*dX[2]);

    // Fill the map.
    readG4BLMap( filename, dsmap, G4BL_offset );
    
    //dsmap.print(std::cout);

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
    int fd = open( edm::FileInPath(filename).fullPath().c_str(), O_RDONLY );
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
    const size_t nx = bfmap._nx;
    const size_t ny = bfmap._ny;
    const size_t nz = bfmap._nz;

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

  //
  // Read one magnetic field map file in G4BL (TD) format.
  //

  void BFieldManagerMaker::readG4BLMap( const string& filename,
					BFMap& bfmap,
					CLHEP::Hep3Vector G4BL_offset ){


    // Debug print
    /*
    ifstream fin1;
    string fname("/tmp/logash/Mu2e_Rotated_Coils_DSMap.txt.bz2");
    fin1.open(fname.c_str());
    boost::iostreams::filtering_istream in;
    decompressFile(fname,in);
    in.push(fin1);
    char cbuf1[128];
    cout << "Debug file: " << endl;
    for( int i=0; i<5; i++ ) { in.getline(cbuf1,128); cout << cbuf1 << endl; }
    in.pop();
    fin1.close();
    */

    // Open the input file.
    ifstream fin;
    fin.open(edm::FileInPath(filename).fullPath().c_str());
    if ( !fin.is_open() ) throw cms::Exception("GEOM")<<"Could not open file "<<filename<<"\n";

    // Skip lines until "data" keyword
    boost::iostreams::filtering_istream in;
    decompressFile(filename,in);
    in.push(fin);
    char cbuf[128];
    boost::regex re("^\\s*data.*$");
    while( ! in.eof() ) {
      in.getline(cbuf,128);
      if( boost::regex_match(cbuf,re) ) break;
    }
    if( in.eof() ) throw cms::Exception("GEOM")<<"Can't find data keyword in "<<filename<<"\n";

    // Expected grid dimentsions.
    const int nx = bfmap._nx;
    const int ny = bfmap._ny;
    const int nz = bfmap._nz;

    // Calculate expected number of lines to read
    const int nrecord = nx*ny*nz;

    // Read data
    double x[3], b[3];
    int nread = 0;
    while( (!in.eof()) && (nread<nrecord) ) {

      // Calculate indexes
      int ix = (nread/(ny*nz));
      int iy = (nread/nz)%ny;
      int iz = nread%nz;

      // Read and parse line
      nread++;
      in.getline(cbuf,128);
      istringstream sin(cbuf);
      if( (sin>>x[0]>>x[1]>>x[2]>>b[0]>>b[1]>>b[2]).fail() ) break;

      // Store the information into the 3d arrays.
      CLHEP::Hep3Vector pos = CLHEP::Hep3Vector(x[0],x[1],x[2])-G4BL_offset;
      bfmap._grid.set (ix, iy, iz, pos);
      bfmap._field.set(ix, iy, iz, CLHEP::Hep3Vector(b[0],b[1],b[2]));
      bfmap._isDefined.set(ix, iy, iz, true);

    }

    if( nread!=nrecord ) {
      throw cms::Exception("GEOM")
	<<"Error while reading "<<filename<<"\n"
	<<"Read "<<nread<<" out of expected "<<nrecord<<" lines.\n"
	<<"Last line:\n"<<cbuf<<"\n";
    }

    return;

  }

  // Compute the size of the array needed to hold the raw data of the field map.
  int BFieldManagerMaker::computeArraySize( int fd, const string& filename ){

    // Get the file size, in bytes, ( info.st_size ).
    struct stat info;
    fstat( fd, &info);

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
