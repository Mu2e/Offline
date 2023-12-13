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

            // calo section
            virtual unsigned                  nDisk()     const override {return nDisks_;}
            virtual const Disk&               disk(int i) const override {return *disks_.at(i);}


              // crystal section
            virtual int                       nCrystal()     const override {return fullCrystalList_.size();}
            virtual const Crystal&            crystal(int i) const override{ return *fullCrystalList_.at(i);}


            // calorimeter geometry information
            virtual const CaloInfo&           caloInfo()     const override {return caloInfo_;}
            virtual const CaloGeomUtil&       geomUtil()     const override {return geomUtil_;}


            // neighbors, indexing
            virtual const std::vector<int>&  neighbors(int crystalId, bool rawMap=false)             const override;
            virtual const std::vector<int>&  nextNeighbors(int crystalId, bool rawMap=false)         const override;
            virtual       std::vector<int>   neighborsByLevel(int crystalId, int level, bool rawMap) const override;
            virtual int                      crystalIdxFromPosition(const CLHEP::Hep3Vector& pos)    const override;
            virtual int                      nearestIdxFromPosition(const CLHEP::Hep3Vector& pos)    const override;

            // get to know me!
            virtual void                     print(std::ostream &os = std::cout) const override;


        private:
            using DiskPtr = std::shared_ptr<Disk>;

            double deltaZ   (const CLHEP::Hep3Vector& p1, const CLHEP::Hep3Vector& p2)  const;
            double deltaPerp(int ic,                      const CLHEP::Hep3Vector& pos) const;

            int                           nDisks_;
            int                           nCrates_;
            int                           nBoards_;
            std::vector<DiskPtr>          disks_;
            std::vector<const Crystal*>   fullCrystalList_; //non-owning crystal pointers
            CaloInfo                      caloInfo_;
            CaloGeomUtil                  geomUtil_;
     };

}

#endif
