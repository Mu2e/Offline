#ifndef CALORIMETER_HH
#define CALORIMETER_HH
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from 
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: Calorimeter.hh,v 1.7 2010/05/25 17:34:55 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2010/05/25 17:34:55 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include <vector>

//
// Mu2e includes
#include "GeometryService/inc/Detector.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "CalorimeterGeom/inc/VaneId.hh"

namespace mu2e {
  namespace calorimeter{
    class Calorimeter: public Detector{

      friend class CalorimeterMaker;

    public:
      Calorimeter(){}
      ~Calorimeter(){};
    
      // Compiler generated copy and assignment constructors
      // should be OK.

      virtual std::string name() const { return "Calorimeter";}


    
      // Check for legal identifiers.
      bool isLegal(VaneId v) const{
         return ( v >-1 && 
                   std::vector<Vane>::size_type(v) <_vanes.size() 
                   );
      };

    
      bool isLegal(const ZSliceId& zid) const{
         return (isLegal(zid._vid) && 
                  zid._zslice >-1   &&
                  std::vector<ZSlice>::size_type(zid._zslice) < getVane(zid._vid).getZSlices().size()
                  );
      }

      typedef std::vector<ZSlice>::size_type stypeRSlice;
      bool isLegal(const RSliceId& rslid ) const{
         return ( isLegal(rslid._zid) &&
                   rslid._rslice > -1   &&
                   std::vector<RSlice>::size_type(rslid._rslice) < getZSlice(rslid._zid).getRSlices().size()
                   );
      }

      bool isLegal(const CrystalId& cid) const{
         return ( isLegal(cid._rslid) &&
                   cid._n > -1       &&
                   std::vector<Crystal>::size_type(cid._n) < getRSlice(cid._rslid).getCrystals().size()
                   );
      }

      // Accessors
      const std::vector<Vane>& getVanes() const{ 
         return _vanes;
      }
    
      const Vane& getVane ( VaneId id) const{ 
         return _vanes.at(id);
      }
    

      const ZSlice& getZSlice ( const ZSliceId& zid ) const{
         return _vanes.at(zid.getVane()).getZSlice(zid);
      }

      const RSlice& getRSlice ( const RSliceId& rslid ) const{
         return _vanes.at(rslid.getVane()).getRSlice(rslid);
      }

      const Crystal& getCrystal ( const CrystalId& cid ) const{
         return _vanes.at(cid.getVane()).getCrystal(cid);
      }

      const Crystal& getCrystal ( CrystalIndex i ) const{
         return _allCrystals.at(i.asInt());
      }

      const std::vector<Crystal>& getAllCrystals() const {return _allCrystals;}

      const std::vector<CrystalDetail>& getCrystalDetails() const{
         return _crystalDetail;
      }


      // 
      // no point in passing reference args for ints and doubles or short strings, but do it for vectors

      uint32_t getNumberOfVanes() const{
         return _numberOfVanes;
      }

      double getROut() const{
         return _rOut;
      }

      double getHalfLength() const{
         return _halfLength;
      }

      uint32_t getNumCrystalRSlices() const{
         return _numberOfRSlices;
      }

      uint32_t getNumCrystalZSlices() const{
         return _numberOfZSlices;
      }

      const CLHEP::Hep3Vector& getCalorimeterCenter() const{
         return _calorimeterCenter;
      }


      const CLHEP::Hep3Vector& getCalorimeterCenterOffset() const{
         return _calorimeterCenterOffset;
      }



      const std::vector<double>& getCalorimeterVaneRotationsX() const{
         return _calorimeterVaneRotationsX;
      }

      const std::vector<double>& getCalorimeterVaneRotationsY() const{
         return _calorimeterVaneRotationsY;
      }

      const std::vector<double>& getCalorimeterVaneRotationsZ() const{
         return _calorimeterVaneRotationsZ;
      }

      const std::vector<double>& getCalorimeterVaneOffsetsX() const{
         return _calorimeterVaneOffsetsX;
      }

      const std::vector<double>& getCalorimeterVaneOffsetsY() const{
         return _calorimeterVaneOffsetsY;
      }

      const std::vector<double>& getCalorimeterVaneOffsetsZ() const{
         return _calorimeterVaneOffsetsZ;
      }



#ifndef __CINT__

      // Loop over all straws and call F.
      // F can be a class with an operator() or a free function.
      template <class F>
      inline void Calorimeter::forAllCrystals ( F& f) const{
         for ( std::vector<Vane>::const_iterator i=_vanes.begin(), e=_vanes.end();
               i !=e; ++i){
           i->forAllCrystals(f);
         }
      }

      template <class F>
      inline void Calorimeter::forAllRSlices ( F& f) const{
         for ( std::vector<Vane>::const_iterator i=_vanes.begin(), e=_vanes.end();
               i !=e; ++i){
           i->forAllRSlices(f);
         }
      }

      template <class F>
      inline void Calorimeter::forAllZSlices ( F& f) const{
         for ( std::vector<Vane>::const_iterator i=_vanes.begin(), e=_vanes.end();
               i !=e; ++i){
           i->forAllZSlices(f);
         }
      }
    
      template <class F>
      inline void Calorimeter::forAllVanes ( F& f) const{
         for ( std::vector<Vane>::const_iterator i=_vanes.begin(), e=_vanes.end();
               i !=e; ++i){
           f(*i);
         }
      }

#endif


    protected:


      // Outer radius and half length ( in z ) of a logical volume that will 
      // just contain the entire tracker.  Use to make the mother volume.
      double _rOut;
      double _halfLength;

      //
      // number of vanes
      uint32_t _numberOfVanes;

      //
      //half length of small side face; crystals assumed square
      double _crystalHalfTrans;

      //
      //half length of long crystal axis
      double _crystalHalfLong;

      //
      //number of rows, rows defined as starting from center and going out
      uint32_t _numberOfRSlices;

      //
      //number of columns, orthogonal to rows
      uint32_t _numberOfZSlices;

      //
      // center of calorimeter
      CLHEP::Hep3Vector _calorimeterCenter;

      //
      // offset of calorimeter wrt above
      CLHEP::Hep3Vector _calorimeterCenterOffset;

      // 
      // inner inscribed circle
      double _rInscribed;


      // Overall azimuthal rotation
      double _phi0;

      //
      // Overall tilt of calorimeter system wrt z-axis
      double _theta0;



      //
      // Individual rotations in x,y,z of vanes from central axis (so in system with theta_0 = 0)
      std::vector<double> _calorimeterVaneRotationsX;
      std::vector<double> _calorimeterVaneRotationsY;
      std::vector<double> _calorimeterVaneRotationsZ;

      //
      // Individual offsets of vanes from ideal locations
      std::vector<double> _calorimeterVaneOffsetsX;
      std::vector<double> _calorimeterVaneOffsetsY;
      std::vector<double> _calorimeterVaneOffsetsZ;


      // Detailed info about each crystal.
      std::vector<CrystalDetail> _crystalDetail;

      // A Calorimeter is made of vanes.
      std::vector<Vane> _vanes;

      // There will be pointers to the objects in this container.
      std::vector<Crystal>  _allCrystals;

      // Needed to complete the second phase of construction.
      void FillPointers1();
      void FillPointers2();

      // On readback from persistency, recursively recompute mutable members.
      void fillPointers() const;

    };
  } // namespace calorimeter
} //namespace mu2e

#endif
