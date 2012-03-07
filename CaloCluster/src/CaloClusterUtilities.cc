//
// General utilities for the calorimeter's studies
//
// $Id: CaloClusterUtilities.cc,v 1.2 2012/03/07 18:00:38 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/03/07 18:00:38 $
//
// Original author G. Pezzullo & G. Tassielli & G. Onorato
//

#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CLHEP/Vector/ThreeVector.h"

#include "CLHEP/Vector/Rotation.h"

//-----------------------------------------
// C++ includes
#include<iostream>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

// Mu2e includes
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/StatusG4.hh"

//-----------------------------------------
using namespace std;

namespace mu2e {


MCCaloUtilities::MCCaloUtilities()
{
  _localRO = 0;
  _localCrystal = 0;
  _localVane = 0;
}

MCCaloUtilities::~MCCaloUtilities()
{
}


void MCCaloUtilities::printOutCaloInfo() {

  art::ServiceHandle<GeometryService> geom;
  GeomHandle<Calorimeter> cg;

  double Hsize = cg->crystalHalfSize();
  double Hleng = cg->crystalHalfLength();
  double ROsize = cg->roHalfSize();

  cout << "Crystal HSize " << Hsize
       << "\nCrystal HLeng " << Hleng
       << "\nRO size " << ROsize << endl;

  for (int i=0; i<4; ++i) {

    Vane const& thevane = cg->getVane(i);
    cout << "Vane " << i << " : "
         << "\nOrigin: " << thevane.getOrigin()
         << "\nLocal origin: " << thevane.getOriginLocal()
         << "\nSize: " << thevane.getSize()
         << "\nRotation: " << *(thevane.getRotation()) << endl;

  }

  int nRO = cg->nRO();

  for (int j=0; j< nRO/2; ++j) {

    int thevane = cg->getVaneByRO(2*j);

    cout << "Crystal n. " << cg->getCrystalByRO(2*j);
    cout << "\tVane " << thevane;


    CLHEP::Hep3Vector cntr = cg->getCrystalOriginByRO(2*j);
    CLHEP::Hep3Vector Xaxis = cg->getCrystalAxisByRO(2*j);
    CLHEP::Hep3Vector Yaxis = cg->getCrystalAxisByRO(2*j).orthogonal();
    CLHEP::Hep3Vector Zaxis(0,0,1);

    cout << "\tcenter " << cntr
         << "\tXaxis " << Xaxis
         << "\tYaxis " << Yaxis;

    CLHEP::Hep3Vector toLeftC = ( (-Hsize) * Yaxis ) + ( (-Hleng) * Xaxis ) + ( (-Hsize) * Zaxis );

    CLHEP::Hep3Vector toRightC = ( (Hsize) * Yaxis ) + ( (Hleng+(2*ROsize)) * Xaxis ) + ( (Hsize) * Zaxis );

    //cout << "toleft" << toLeftC
    //     << "\ttoright" << toRightC;


    CLHEP::Hep3Vector lcorner = cntr + toLeftC;
    CLHEP::Hep3Vector rcorner = cntr + toRightC;

    cout << "\tleft corner " << lcorner << "\tright corner " << rcorner << endl;

  }

}

void MCCaloUtilities::setTrackAndRO(const art::Event & event,
                                    std::string const &_g4ModuleLabel,
                                    SimParticleCollection::key_type track,
                                    unsigned RO){

  _localRO = RO;
  art::ServiceHandle<GeometryService> geom;
  GeomHandle<Calorimeter> cg;
  _localCrystal  = cg->getCrystalByRO(_localRO);
  _localVane = cg->getVaneByRO(_localRO);

  art::Handle<SimParticleCollection> simParticles;
  event.getByLabel(_g4ModuleLabel, simParticles);

  SimParticle const& sim = simParticles->at(track);

  CLHEP::Hep3Vector origin = sim.startPosition();

  _startingVane = getStartingVane(origin);

  _fromOutside = (_startingVane == -1);

  _primary = !sim.hasParent();

  _generated = sim.fromGenerator();

  //    art::Handle<PhysicalVolumeInfoCollection> volumes;
  //event.getRun().getByLabel(_g4ModuleLabel, volumes);

  //PhysicalVolumeInfo const& volInfob = volumes->at(sim.startVolumeIndex());
  //PhysicalVolumeInfo const& volInfoe = volumes->at(sim.endVolumeIndex());
  //cout << "start: " << sim.startVolumeIndex() << "   " << volInfob.name() << "  " << volInfob.copyNo() << endl;
  //cout << "end:   " << sim.endVolumeIndex()   << "   " << volInfoe.name() << "  " << volInfoe.copyNo() << endl;
  //cout << "Start position: " << sim.startPosition() << '\n'
  //     << "End position:   " << sim.endPosition() << endl;
  //cout << "Particle process code " << sim.creationCode().name() << endl;
  //bool ID(false);
  // if ( sim.hasParent()) {
  //   cout  << "Parent id " << sim.parentId() << endl;
  //   ID = true;
  // }
  // if ( sim.fromGenerator() ) {
  //   cout << "Is from generator" << endl;
  //   ID =  true;
  // }
  // if  (!ID) {
  //   cout << "dunno where it come from" << endl;
  // }

}

bool MCCaloUtilities::fromOutside() {
  return _fromOutside;
}

bool MCCaloUtilities::primary() {
  return _primary;
}

bool MCCaloUtilities::generated() {
  return _generated;
}

int MCCaloUtilities::startingVane() {
  return _startingVane;
}

int MCCaloUtilities::localVane() {
  return _localVane;
}

int MCCaloUtilities::getStartingVane(CLHEP::Hep3Vector origin) {

  art::ServiceHandle<GeometryService> geom;
  GeomHandle<Calorimeter> cg;

  for (size_t i=0; i<cg->nVane(); ++i) {

    Vane const & vane = cg->getVane(i);
    CLHEP::Hep3Vector rsize = *(vane.getRotation()) * vane.getSize();
    //cout << "size " << vane.getSize() << " and rotated is " << rsize << endl;
    CLHEP::Hep3Vector vaneOr = vane.getOrigin();


    if (fabs(origin.getZ() - vaneOr.getZ()) <= fabs(rsize.getZ())) {
      if (fabs(origin.getY() - vaneOr.getY()) <= fabs(rsize.getY())) {
        if (fabs(origin.getX() - vaneOr.getX()) <= fabs(rsize.getX())) {
          return i;
        }
      }
    }
  }

  return -1;

}





//------------------------------------------------
std::string & TOUpper(std::string &in) {
  std::string::iterator i = in.begin();
  std::string::iterator end = in.end();
  while (i != end) {
    *i = std::toupper((unsigned char)*i);
    ++i;
  }
  return in;
}

void cog(CaloCluster &cluster){
        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        if(! geom->hasElement<Calorimeter>() ) return;
        GeomHandle<Calorimeter> cg;
        CLHEP::Hep3Vector res(1., 1., 1.);

        CLHEP::Hep3Vector  resError(1e-1, 1e-1, 1e-1);
        //int tmpZ = 0, tmpR = 0;//, vane = -1;
        //int timeZ = 0, timeR = 0;

        //get crystal's geometrical information
        double cryHalfSize      =      cg->crystalHalfSize();
        //double cryHalfLength    =    cg->crystalHalfLength();

        //double tmpTime = 1e10;

        //using the expression for the RMS of a flat distribution, I calculate the value of the RMS for the face of the crystals
        //double RMScryHalfSize   =   2.*cryHalfSize / TMath::Sqrt(12.);
        //double RMScryHalfLenght = 2.*cryHalfLength / TMath::Sqrt(12.);
        //double defaultError     = 1e-06;//default from geant4

//        double  deltaZ   = 0.;
//        double  deltaR   = 0.;
//        double  deltaXc  = 0.;
//        double  deltaYc  = 0.;
//        double  deltaZc  = 0.;
        double  tmpEq    = 0.;
        double  showerDepth = 1.5;

        bool isfirstCrystal=true;
        float R =0.;
        float Z =0.;

        //int vaneid = -1;
        //int counter(0);
        //cout<< "caloClusterSize = "<< cluster.caloCrystalHitsPtrVector.size() <<endl;

        for( CaloCrystalHitPtrVector::iterator itCD = cluster.caloCrystalHitsPtrVector.begin(); itCD != cluster.caloCrystalHitsPtrVector.end(); ++itCD){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);
                //now I take the first crystal of the cluster as a point of reference for calculating the cog
                if(isfirstCrystal){
                        //vaneid = cg->getVaneByRO(thehit.id() );
                        res = cg->getCrystalOriginByRO(thehit.id());
                        isfirstCrystal = false;
                        //tmpZ = cg->getCrystalZByRO(thehit.id());
//                        tmpR = cg->getCrystalRByRO(thehit.id());
//                        cout <<"vaneId = "<<vaneid<<endl;
//                        cout << "tmpZ = "<< tmpZ <<endl;
//                        cout << "tmpR = "<< tmpR <<endl;
//                        cout<< "cryref.X = " << res.getX()<<endl;
//                        cout<< "cryref.Y = " << res.getY()<<endl;
//                        cout<< "cryref.Z = " << res.getZ()<<endl;
//                        ++tmpZ;
//                        ++tmpR;
                        //vane = cg->getVaneByRO(thehit.id());
                }


                //cout<< "ctrystals read = " << ++counter<<endl;

                //Get Z and R from readout
                double tZ = cg->getCrystalZByRO(thehit.id());
                double tR = cg->getCrystalRByRO(thehit.id());


//                if( (*itCD)->time() < tmpTime){
//                        tmpTime = (*itCD)->time();
//                        timeR = tR;
//                        timeZ = tZ;
//                }

                //cout << "tZ = "<<tZ<<", tR = "<<tR<<"energy_i = "<< (*itCD)->energyDep()<<endl;

//                tZ += 1.0;
//                tR += 1.0;
                //Multiply Z and R for thhe crystal's energy, which is the weight we use in that algorithm
                tZ *=(*itCD)->energyDep();
                tR*=(*itCD)->energyDep();

                Z += tZ;
                R += tR;

                //Calculate the sum of the square of the weight for derive the error of the weighted mean
                tmpEq += std::pow((*itCD)->energyDep(), 2.);
        }

        R/=cluster.energyDep;
        Z/=cluster.energyDep;
       // cout<< "R = "<< R<<", Z = "<<Z<< ", cluster.energyDep = "<< cluster.energyDep<<endl;
//        R/=cluster.caloCrystalHitsPtrVector.size();
//        Z/=cluster.caloCrystalHitsPtrVector.size();

        //-----------------

        //const Vane & vane = cg->getVane(vaneid);

        // Crystal center in vane coordinates
//        CLHEP::Hep3Vector vlocal(0/*-1.*cg->roHalfThickness()  */,
//                        (2.*R-(double)cg->nCrystalR()+1.0)*cryHalfSize,
//                        (2.*Z-(double)cg->nCrystalZ()+1.0)*cryHalfSize );

        //res =  vane.getOrigin() + (vane.getRotation()->inverse())*vlocal;




        //-------------


        //Get the relative distance of the cog from "res" (crystal took as point of reference)
//        deltaZ = Z - tmpZ;
//        deltaR = R - tmpR;

//        CLHEP::Hep3Vector vlocal(1.*cryHalfLength /*+cg->roHalfThickness()*/,
//                                (2.*deltaR+1.0)*cryHalfSize,
//                                (2.*deltaZ)*cryHalfSize );
//        const CLHEP::HepRotation ro;
//        ro = CLHEP::HepRotation(res, 1.);

        //result = res + (ro.inverse())*vlocal;
//        res =  res + (vane.getRotation()->inverse())*vlocal;



        //deltaZc = (deltaZ - 1.0)*cryHalfSize*2.;

        //On following we calculate the correct position of the cog, which depends on the vane's orientation
//        if( vaneid == 3){
//                deltaYc = (deltaR )*cryHalfSize*2.;
//                deltaZc = (deltaZ )*cryHalfSize*2.;
//                deltaXc = cryHalfLength + 2.0*cg->roHalfThickness();
//
//                resError.setX(defaultError);
//                resError.setY( RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
//        }
//        else if(vaneid == 2){
//                deltaXc = (deltaR )*cryHalfSize*2.;
//                deltaZc = (deltaZ)*cryHalfSize*2.;
//                deltaYc = -1.*cryHalfLength - 2.0*cg->roHalfThickness();
//
//                resError.setY(defaultError);
//                resError.setX(RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
//
//        }else if(vaneid == 1){
//                deltaYc = -1.*(deltaR)*cryHalfSize*2. ;
//                deltaZc = (deltaZ)*cryHalfSize*2.;
//                deltaXc = -1.*cryHalfLength - 2.0*cg->roHalfThickness();
//
//                resError.setX(defaultError);
//                resError.setY(RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
//        }else if(vaneid == 0){
//                deltaXc = -1.0*(deltaR)*cryHalfSize*2.;
//                deltaZc = (deltaZ)*cryHalfSize*2.;
//                deltaYc = cryHalfLength + 2.0*cg->roHalfThickness();
//
//                resError.setY(defaultError);
//                resError.setX(RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
//        }



//        res.setX(res.getX() + deltaXc/* + 3904.*/);//value used to shift in tracker coordinate system
//        res.setZ(res.getZ() + deltaZc /*- 10200*/);//value used to shift in tracker coordinate system
//        res.setY(res.getY() + deltaYc);

//        resError.setZ( RMScryHalfSize /TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);

        res.setY( (R*2.+ 1.0)*cryHalfSize);//cg->nCrystalR()*cryHalfSize);//
        res.setZ( (Z*2.+1.0)*cryHalfSize);//cg->nCrystalZ()*cryHalfSize);//
        res.setX(- showerDepth);
//        Save cog#Vector and its errorVector
//        cluster.cog3VectorError = resError;
//        cout<<"inserisco row e column COG..."<<endl;
        cluster.cogRow = (int) R;
        cluster.cogColumn = (int) Z;
//        cout<<"fatto inserimento!!!! row e column COG..."<<endl;
        cluster.cog3Vector = res;

}


CLHEP::Hep3Vector cog_depth(CaloCluster &cluster, double depth){

        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        GeomHandle<Calorimeter> cg;
        CLHEP::Hep3Vector res(1., 1., 1.);

        CLHEP::Hep3Vector  resError(1e-1, 1e-1, 1e-1);
//        int tmpZ = 0, tmpR = 0;//, vane = -1;
//        int timeZ = 0, timeR = 0;

        //get crystal's geometrical information
        double cryHalfSize      =      cg->crystalHalfSize();
        //double cryHalfLength    =    cg->crystalHalfLength();

       // double tmpTime = 1e10;

        //using the expression for the RMS of a flat distribution, I calculate the value of the RMS for the face of the crystals
        //double RMScryHalfSize   =   2.*cryHalfSize / TMath::Sqrt(12.);
        //double RMScryHalfLenght = 2.*cryHalfLength / TMath::Sqrt(12.);
       // double defaultError     = 1e-06;//default from geant4

//        double  deltaZ   = 0.;
//        double  deltaR   = 0.;
        double  tmpEq    = 0.;
        double  showerDepth = depth;

        bool isfirstCrystal=true;
        float R =0.;
        float Z =0.;


        for( CaloCrystalHitPtrVector::iterator itCD = cluster.caloCrystalHitsPtrVector.begin(); itCD != cluster.caloCrystalHitsPtrVector.end(); ++itCD){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);
                //now I take the first crystal of the cluster as a point of reference for calculating the cog
                if(isfirstCrystal){
                        res = cg->getCrystalOriginByRO(thehit.id());
                        isfirstCrystal = false;
                      //  tmpZ = cg->getCrystalZByRO(thehit.id());
                    //    tmpR = cg->getCrystalRByRO(thehit.id());
                }


                //cout<< "ctrystals read = " << ++counter<<endl;

                //Get Z and R from readout
                double tZ = cg->getCrystalZByRO(thehit.id());
                double tR = cg->getCrystalRByRO(thehit.id());


//                if( (*itCD)->time() < tmpTime){
//                        tmpTime = (*itCD)->time();
//                        timeR = tR;
//                        timeZ = tZ;
//                }

               //Multiply Z and R for thhe crystal's energy, which is the weight we use in that algorithm
                tZ *=(*itCD)->energyDep();
                tR*=(*itCD)->energyDep();

                Z += tZ;
                R += tR;

                //Calculate the sum of the square of the weight for derive the error of the weighted mean
                tmpEq += std::pow((*itCD)->energyDep(), 2.);
        }

        R/=cluster.energyDep;
        Z/=cluster.energyDep;


        //Get the relative distance of the cog from "res" (crystal took as point of reference)
//        deltaZ = Z - tmpZ;
//        deltaR = R - tmpR;


//        resError.setZ( RMScryHalfSize /TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);

        res.setY( (R*2.+ 1.0)*cryHalfSize);//cg->nCrystalR()*cryHalfSize);//
        res.setZ( (Z*2.+1.0)*cryHalfSize);//cg->nCrystalZ()*cryHalfSize);//
        res.setX(- showerDepth);
//        Save cog#Vector and its errorVector
//        cluster.cog3VectorError = resError;
       // cluster._impactPoint = res;

        return res;

}


//on the following we implement an algorithm for the cog which uses the logarithm of the energy as weight (w_{i}). This is not correct, as references show, we need to calculate from simulation an offset to add at each w_{i}
//void LOGcog(CaloCluster &cluster){
CLHEP::Hep3Vector LOGcog(CaloCluster &cluster, double w){
        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;

        //if(! geom->hasElement<Calorimeter>() ) return;
        GeomHandle<Calorimeter> cg;
        CLHEP::Hep3Vector res, resError;
        int tmpZ = 0, tmpR = 0, vane = -1;

        //get crystal's geometrical information
        double cryHalfSize      =      cg->crystalHalfSize();
        double cryHalfLength    =    cg->crystalHalfLength();

        //using the expression for the RMS of a flat distribution, I calculate the value of the RMS for the face of the crystals
        double RMScryHalfSize   =   2.*cryHalfSize / TMath::Sqrt(12.);
        //double RMScryHalfLenght = 2.*cryHalfLength / TMath::Sqrt(12.);
        double defaultError     = 1e-06;//default from geant4

        double  deltaZ   = 0.;
        double  deltaR   = 0.;
        double  deltaXc  = 0.;
        double  deltaYc  = 0.;
        double  deltaZc  = 0.;
        double  tmpEq    = 0.;
        double  sumEi    = 0.;
        double offSet    = w;//0.;
        double weight    = 0.;
        double sumW      = 0.;

        bool isfirstCrystal=true;
        float R =0.;
        float Z =0.;


        for( CaloCrystalHitPtrVector::iterator itCD = cluster.caloCrystalHitsPtrVector.begin(); itCD != cluster.caloCrystalHitsPtrVector.end(); ++itCD){

                sumEi += (*itCD)->energyDep();
        }

        for( CaloCrystalHitPtrVector::iterator itCD = cluster.caloCrystalHitsPtrVector.begin(); itCD != cluster.caloCrystalHitsPtrVector.end(); ++itCD){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);

                //now I take the first crystal of the cluster as a point of reference for calculating the cog
                if(isfirstCrystal){
                        res = cg->getCrystalOriginByRO(thehit.id());
                        isfirstCrystal = false;
                        tmpZ = cg->getCrystalZByRO(thehit.id());
                        tmpR = cg->getCrystalRByRO(thehit.id());
                        vane = cg->getVaneByRO(thehit.id());
                }

                //Get Z and R from readout
                double tZ = cg->getCrystalZByRO(thehit.id());
                double tR = cg->getCrystalRByRO(thehit.id());

                //Multiply Z and R for thhe crystal's energy, which is the weight we use in that algorithm
                weight = offSet + TMath::Log((*itCD)->energyDep()/sumEi);
                if(weight < 0.0){
                        weight = 0.0;
                }
                tZ *= weight;
                tR *= weight;

                Z += tZ;
                R += tR;
                sumW += weight;

                //Calculate the sum of the square of the weight for derive the error of the weighted mean
                tmpEq += std::pow(weight, 2.);
        }

        R /= sumW;
        Z /= sumW;

        //Get the relative distance of the cog from "res" (crystal took as point of reference)
        deltaZ = Z - tmpZ;
        deltaR = R - tmpR;


        deltaZc = (deltaZ )*cryHalfSize*2.;

        //On following we calculate the correct position of the cog, which depends on the vane's orientation
        if( vane == 3){
                       deltaYc = (deltaR )*cryHalfSize*2.;
                       //deltaZc = (deltaZ )*cryHalfSize*2.;
                       deltaXc = cryHalfLength + 2.0*cg->roHalfThickness();

                       resError.setX(defaultError);
                       resError.setY( RMScryHalfSize/TMath::Sqrt(cluster.clusterSize) * TMath::Sqrt(tmpEq)/cluster.energyDep);
               }
               else if(vane == 2){
                       deltaXc = (deltaR )*cryHalfSize*2.;
                       //deltaZc = (deltaZ)*cryHalfSize*2.;
                       deltaYc = -1.*cryHalfLength - 2.0*cg->roHalfThickness();

                       resError.setY(defaultError);
                       resError.setX(RMScryHalfSize/TMath::Sqrt(cluster.clusterSize) * TMath::Sqrt(tmpEq)/cluster.energyDep);

               }else if(vane == 1){
                       deltaYc = -1.*(deltaR)*cryHalfSize*2. ;
                       //deltaZc = (deltaZ)*cryHalfSize*2.;
                       deltaXc = -1.*cryHalfLength - 2.0*cg->roHalfThickness();

                       resError.setX(defaultError);
                       resError.setY(RMScryHalfSize/TMath::Sqrt(cluster.clusterSize) * TMath::Sqrt(tmpEq)/cluster.energyDep);
               }else if(vane == 0){
                       deltaXc = -1.0*(deltaR)*cryHalfSize*2.;
                       //deltaZc = (deltaZ)*cryHalfSize*2.;
                       deltaYc = cryHalfLength + 2.0*cg->roHalfThickness();

                       resError.setY(defaultError);
                       resError.setX(RMScryHalfSize/TMath::Sqrt(cluster.clusterSize) * TMath::Sqrt(tmpEq)/cluster.energyDep);
               }



        res.setX(res.getX() + deltaXc /*+ 3904.*/);//value used to shift in tracker coordinate system
        res.setZ(res.getZ() + deltaZc /*- 10200*/);//value used to shift in tracker coordinate system
        res.setY(res.getY() + deltaYc);

        resError.setZ( RMScryHalfSize /TMath::Sqrt(cluster.clusterSize) * TMath::Sqrt(tmpEq)/cluster.energyDep);

        //Save cog#Vector and its errorVector
        //        cluster.cog3VectorError = resError;
        //        cluster._impactPoint = res;
        return res;
}


}
