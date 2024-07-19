#ifndef CalorimeterGeom_DiskCalorimeter_hh
#define CalorimeterGeom_DiskCalorimeter_hh
//
// Hold informations about the disk calorimeter
//
// Original author B. Echenard
//
// Note 1: conversion of crystal <-> readout id
//         readout_id = crystal_id*nRoPerCrystal ... crystal_id*nRoPerCrystal + nRoPerCrystal-1

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/CaloInfo.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <vector>
#include <memory>



namespace mu2e {


    class DiskCalorimeter: public Calorimeter {


        friend class DiskCalorimeterMaker;

        public:
            DiskCalorimeter();

            size_t                    nDisks()          const override {return disks_.size();}
            const Disk&               disk(size_t i)    const override {return *disks_.at(i);}
            const DiskPtrs&           diskPtrs()        const override {return disks_;}

            size_t                    nCrystals()       const override {return fullCrystalList_.size();}
            const Crystal&            crystal(size_t i) const override {return *fullCrystalList_.at(i);}
            const CrystalPtrs&        crystalPtrs()     const override {return  fullCrystalList_;}

            const CaloInfo&           caloInfo()        const override {return caloInfo_;}
            const CaloGeomUtil&       geomUtil()        const override {return geomUtil_;}

            const std::vector<int>&  neighbors(int crystalId, bool rawMap=false)             const override;
            const std::vector<int>&  nextNeighbors(int crystalId, bool rawMap=false)         const override;
                  std::vector<int>   neighborsByLevel(int crystalId, int level, bool rawMap) const override;
            int                      crystalIdxFromPosition(const CLHEP::Hep3Vector& pos)    const override;
            int                      nearestIdxFromPosition(const CLHEP::Hep3Vector& pos)    const override;

            void                     print(std::ostream &os = std::cout) const override;


        private:
            double deltaZ    (const CLHEP::Hep3Vector& p1, const CLHEP::Hep3Vector& p2)  const;
            double deltaPerp (int ic,                      const CLHEP::Hep3Vector& pos) const;

            std::vector<DiskPtr>          disks_;
            std::vector<const Crystal*>   fullCrystalList_; //non-owning crystal pointers
            CaloInfo                      caloInfo_;
            CaloGeomUtil                  geomUtil_;
     };

}

#endif
