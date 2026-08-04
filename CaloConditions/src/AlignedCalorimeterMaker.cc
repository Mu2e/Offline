#include "cetlib_except/exception.h"
#include "Offline/CaloConditions/inc/AlignedCalorimeterMaker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>


namespace mu2e {

  //-----------------------------------------------------------------------------------------
  void AlignedCalorimeterMaker::alignCalorimeter(AlignedCalorimeterMaker::ptr_t ptr,
                                                 const std::vector<CalAlignParams>& disk_align_params)
  {
    DiskCalorimeter& cal = *ptr;

    for (auto& disk : cal.disks()) {
      // Use the const_cast bazooka here - this is the only place where
      // Disk parameters should be modified by an external agent
      auto& diskcc = const_cast<Disk&>(disk);
      const auto& disk_align = disk_align_params.at(diskcc.id());

      CLHEP::Hep3Vector  shift(disk_align.dx(),disk_align.dy(),disk_align.dz());
      CLHEP::HepRotation rotation(CLHEP::HepRotation::IDENTITY);
      rotation.rotateX(disk_align.rx());
      rotation.rotateY(disk_align.ry());
      rotation.rotateZ(disk_align.rz());

      diskcc.moveDisk(shift,rotation);
    }
  }



  //-----------------------------------------------------------------------------------------
  std::vector<CalAlignParams> AlignedCalorimeterMaker::readFile(const std::string& fileName, uint16_t nRowMax)
  {
    ConfigFileLookupPolicy configFile;
    std::string fullFileName = configFile(fileName);
    if (_config.verbose()>0) {
      std::cout << "AlignedCalorimeterMaker::fromFcl reading from " << fullFileName << "\n";
    }

    std::ifstream ordFile;
    ordFile.open(fullFileName);
    if (!ordFile.is_open()) {
      throw cet::exception("ALIGNEDCAL_OPEN_FAILED")
        << " failed to open file " << fullFileName << "\n";
    }

    std::vector<CalAlignParams> align_params;
    uint16_t nRead(0);
    std::string line;
    while (std::getline(ordFile, line)) {
      uint16_t index;
      float dx,dy,dz,rx,ry,rz;

      std::istringstream iss(line);
      if (!(iss >> index >> dx >> dy >> dz >> rx >> ry >> rz)) {
        throw cet::exception("ALIGNEDCAL_FORMAT")
        << "invalid format at line "<<nRead+1<<"\n";
      }

      // Check that there is nothing left on the line
      iss >> std::ws;
      if (!iss.eof()) {
        throw cet::exception("ALIGNEDCAL_FORMAT")
        << "invalid format at line "<<nRead+1<<"\n";
      }

      if (index >= nRowMax) {
        throw cet::exception("ALIGNEDCAL_RANGE") << "AlignedCalorimeterMaker read invalid offlineId "
        << index << "\n";
      }
      align_params.push_back(CalAlignParams(index,dx,dy,dz,rx,ry,rz));
      ++nRead;
    }

    if(nRead != nRowMax) {
      throw cet::exception("ALIGNEDCAL_COUNT")
        << "AlignedCalorimeterMaker read the wrong number of id's "
        << nRead << ", expected " << nRowMax << "\n";
    }

    return align_params;
  }



  //-----------------------------------------------------------------------------------------
  std::vector<CalAlignParams> AlignedCalorimeterMaker::readDb(CalAlignElement::cptr_t cptr, uint16_t nRowMax)
  {
    std::vector<CalAlignParams> align_params;
    uint16_t nRead(0);

    for (auto const& row : cptr->rows()) {
      auto index = row.index();
      auto dx    = row.dx();
      auto dy    = row.dy();
      auto dz    = row.dz();
      auto rx    = row.rx();
      auto ry    = row.ry();
      auto rz    = row.rz();

      if (index >= nRowMax) {
        throw cet::exception("AlignedCalorimeterMAKER_RANGE") << "AlignedCalorimeterMaker read invalid offlineId " << index << "\n";
      }
      align_params.push_back(CalAlignParams(index,dx,dy,dz,rx,ry,rz));
      ++nRead;
    }

    if(nRead != nRowMax) {
      throw cet::exception("AlignedCalorimeterMAKER_COUNT")
        << "AlignedCalorimeterMaker read the wrong number of id's "
        << nRead << ", expected " << nRowMax << "\n";
    }
    return align_params;
  }


  //-----------------------------------------------------------------------------------------
  AlignedCalorimeterMaker::ptr_t AlignedCalorimeterMaker::fromFcl()
  {
    //Deep copy of the calorimeter to keep the original untouched
    GeomHandle<DiskCalorimeter> cal_h;
    auto ptr = std::make_shared<DiskCalorimeter>(*cal_h);

    if ( _config.verbose() > 0 )
      std::cout << "AlignedCalorimeterMaker::fromFcl now zero aligning Calorimeter\n";

    auto disk_align_params = readFile(_config.filenameDisk(), CaloConst::_nDisk);
    alignCalorimeter(ptr, disk_align_params);

    return ptr;
  }

  //-----------------------------------------------------------------------------------------
  AlignedCalorimeterMaker::ptr_t AlignedCalorimeterMaker::fromDb(CalAlignDisk::cptr_t cad_p)
  {
    GeomHandle<DiskCalorimeter> cal_h;
    auto  ptr = std::make_shared<DiskCalorimeter>(*cal_h);

    if ( _config.verbose() > 0 )
      std::cout << "AlignedCalorimeterMaker::fromDb now aligning Calorimeter \n";

    auto disk_align_params = readDb(cad_p, CaloConst::_nDisk);
    alignCalorimeter(ptr, disk_align_params);

    return ptr;
  }
}
