//
// General utilities for the calorimeter's studies
//
// $Id: CaloClusterUtilities.cc,v 1.4 2012/04/10 20:28:57 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/04/10 20:28:57 $
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


double cry(double val){
        double size = 30.0000;
        int index = (int) (val / size);
        double result =  val - index*size;

        return result;
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

        double  tmpEq    = 0.;
        double  showerDepth = 1.5;

        bool isfirstCrystal=true;
        float R =0.;
        float Z =0., RZ = 0., RQ = 0.;

        //int vaneid = -1;
        //int counter(0);
        //cout<< "calosize() = "<< cluster.caloCrystalHitsPtrVector.size() <<endl;

//        cout<<"---------cog 1 --------"<<endl;
        int c = 0;
        for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){
//                cout<<"----------cog1."<<c<<" --------"<<endl;
                ++c;
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
//                cout<<"----------cog1."<<c<<".1"<<" --------"<<endl;
                CaloHit const& thehit = *ROIds.at(0);
//                cout<<"----------cog1."<<c<<".2"<<" --------"<<endl;
                //now I take the first crystal of the cluster as a point of reference for calculating the cog
                if(isfirstCrystal){
                        //vaneid = cg->getVaneByRO(thehit.id() );
//                        cout<<"----------cog1."<<c<<".2@"<<" --------"<<endl;

                        res = cg->getCrystalOriginByRO(thehit.id());
//                        cout<<"----------cog1."<<c<<".2@@"<<" --------"<<endl;

                        isfirstCrystal = false;
                }


                //cout<< "ctrystals read = " << ++counter<<endl;

                //Get Z and R from readout
//                cout<<"----------cog1."<<c<<".2@@@"<<" --------"<<endl;

                double tZ = cg->getCrystalZByRO(thehit.id());
//                cout<<"----------cog1."<<c<<".2@@@@"<<" --------"<<endl;

                double tR = cg->getCrystalRByRO(thehit.id());
//                cout<<"----------cog1."<<c<<".3"<<" --------"<<endl;
                Z += tZ;
                RQ += pow(tR, 2);
                R += tR;
                RZ += tZ*tR;
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
//                cout<<"----------cog1."<<c<<".4"<<" --------"<<endl;

//                Z += tZ;
//                RQ += pow(tR, 2);
//                R += tR;
//                RZ += tZ*tR;

                //Calculate the sum of the square of the weight for derive the error of the weighted mean
                tmpEq += std::pow((*itCD)->energyDep(), 2.);
        }
//        cout<<"-------------cog 2--------------"<<endl;
        float m = R*Z - RZ;
        m /= (pow(R,2) - RQ);

        R/=cluster.energyDep();
        Z/=cluster.energyDep();
//        cout<<"-------------cog 3--------------"<<endl;


        //        res.setX(res.getX() + deltaXc/* + 3904.*/);//value used to shift in tracker coordinate system
        //        res.setZ(res.getZ() + deltaZc /*- 10200*/);//value used to shift in tracker coordinate system
        //        res.setY(res.getY() + deltaYc);

        //        resError.setZ( RMScryHalfSize /TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);

        res.setY( (R*2.+ 1.0)*cryHalfSize);//cg->nCrystalR()*cryHalfSize);//
        res.setZ( (Z*2.+1.0)*cryHalfSize);//cg->nCrystalZ()*cryHalfSize);//
        res.setX(- showerDepth);

        cluster.SetShowerDir(m);
        cluster.SetCogRow((int) R);
        cluster.SetCogColumn((int) Z);
        cluster.SetCog3Vector(res);
}


void cog_depth(CaloCluster &cluster, double depth, ClusterMap &clusterMap){
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

        double  tmpEq    = 0.;
        double  showerDepth = depth;

        bool isfirstCrystal=true;
        float R =0.;
        float Z =0., V= 0.0, W=0.0, VW=0.0, VQ = 0.0, RQ=0.0, RZ=0.0;

        clusterMap._cluSize = cluster.size();

        for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);
                //now I take the first crystal of the cluster as a point of reference for calculating the cog
                if(isfirstCrystal){
                        res = cg->getCrystalOriginByRO(thehit.id());
                        isfirstCrystal = false;
                        //  tmpZ = cg->getCrystalZByRO(thehit.id());
                        //    tmpR = cg->getCrystalRByRO(thehit.id());
                        int vane = cg->getVaneByRO(thehit.id());
                        clusterMap._vane = vane;

                }


                //cout<< "ctrystals read = " << ++counter<<endl;

                //Get Z and R from readout
                double tZ = cg->getCrystalZByRO(thehit.id());
                double tR = cg->getCrystalRByRO(thehit.id());

                V += tR;
                W += tZ;
                VW += tR*tZ;
                VQ += pow(tR,2);

                clusterMap._rowVec.push_back(tR);
                clusterMap._columnVec.push_back(tZ);
                clusterMap._cryEdepVec.push_back((*itCD)->energyDep());

                clusterMap._COGrowVec.push_back(tR);
                clusterMap._COGcolumnVec.push_back(tZ);


                //                if( (*itCD)->time() < tmpTime){
                //                        tmpTime = (*itCD)->time();
                //                        timeR = tR;
                //                        timeZ = tZ;
                //                }

                //Multiply Z and R for thhe crystal's energy, which is the weight we use in that algorithm
                tZ *=(*itCD)->energyDep();
                tR*=(*itCD)->energyDep();

                Z += tZ;
                RQ += pow(tR, 2);
                R += tR;
                RZ += tZ*tR;
                //Calculate the sum of the square of the weight for derive the error of the weighted mean
                tmpEq += std::pow((*itCD)->energyDep(), 2.);
        }

        int size = cluster.size();
        clusterMap._COGcrySize = size;
        float m = V*W - size*VW;//R*Z - RZ;//
        m /= (pow(V,2) - size*VQ);//(pow(R,2) - RQ);


        float D = RQ - pow(R,2);
        float errM = (1- pow(m,2)*D);
        if(errM>=0.0){
                errM = sqrt(errM);
        }else{
                errM = 1e-4;
        }

        R/=cluster.energyDep();
        Z/=cluster.energyDep();

        clusterMap._cluCogRow = R;
        clusterMap._cluCogColumn = Z;
        //        resError.setZ( RMScryHalfSize /TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);

        res.setY( (R*2.+ 1.0)*cryHalfSize);//cg->nCrystalR()*cryHalfSize);//
        res.setZ( (Z*2.+ 1.0)*cryHalfSize);//cg->nCrystalZ()*cryHalfSize);//
        res.setX(- showerDepth);
        //        Save cog#Vector and its errorVector
        //        cluster.cog3VectorError = resError;
        // cluster._impactPoint = res;

        clusterMap._cluCOG = res;
        clusterMap._showerDir = m;
        clusterMap._errShowerDir = errM;
        //  return res;

}


//on the following we implement an algorithm for the cog which uses the logarithm of the energy as weight (w_{i}). This is not correct, as references show, we need to calculate from simulation an offset to add at each w_{i}
//void LOGcog(CaloCluster &cluster){
//CLHEP::Hep3Vector LOGcog(CaloCluster &cluster, double w, double depth){
void LOGcogMap(CaloCluster &cluster, double w, double depth, ClusterMap &clusterMap ){

        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;

        //if(! geom->hasElement<Calorimeter>() ) return;
        GeomHandle<Calorimeter> cg;
        CLHEP::Hep3Vector res, resError;
        //int tmpZ = 0, tmpR = 0, vane = -1;

        //get crystal's geometrical information
        double cryHalfSize      =      cg->crystalHalfSize();
        //double cryHalfLength    =    cg->crystalHalfLength();

        //using the expression for the RMS of a flat distribution, I calculate the value of the RMS for the face of the crystals
        //double RMScryHalfSize   =   2.*cryHalfSize / TMath::Sqrt(12.);
        //double RMScryHalfLenght = 2.*cryHalfLength / TMath::Sqrt(12.);
        //double defaultError     = 1e-06;//default from geant4

        double  tmpEq    = 0.;
        double  sumEi    = 0.;
        double offSet    = w;
        double weight    = 0.;
        double sumW      = 0.;

        bool isfirstCrystal=true;
        float R =0.;
        float Z =0.;

        int count = 0;
        for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){

                sumEi += (*itCD)->energyDep();
        }

        clusterMap._cluSize = cluster.size();

        for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);

                //now I take the first crystal of the cluster as a point of reference for calculating the cog
                if(isfirstCrystal){
                        res = cg->getCrystalOriginByRO(thehit.id());
                        isfirstCrystal = false;
                        //tmpZ = cg->getCrystalZByRO(thehit.id());
                        //tmpR = cg->getCrystalRByRO(thehit.id());
                        int vane = cg->getVaneByRO(thehit.id());
                        clusterMap._vane = vane;
                }

                //Get Z and R from readout
                double tZ = cg->getCrystalZByRO(thehit.id());
                double tR = cg->getCrystalRByRO(thehit.id());

                clusterMap._rowVec.push_back(tR);
                clusterMap._columnVec.push_back(tZ);

                //Multiply Z and R for thhe crystal's energy, which is the weight we use in that algorithm
                weight = offSet + TMath::Log((*itCD)->energyDep()/sumEi);
                if(weight < 0.0){
                        weight = 0.0;
                }else{
                        count++;
                        clusterMap._COGrowVec.push_back(tR);
                        clusterMap._COGcolumnVec.push_back(tZ);
                        clusterMap._cryEdepVec.push_back((*itCD)->energyDep());
                }
                tZ *= weight;
                tR *= weight;

                Z += tZ;
                R += tR;
                sumW += weight;

                //Calculate the sum of the square of the weight for derive the error of the weighted mean
                tmpEq += std::pow(weight, 2.);
        }
//        if(count <=1){
//                cout<< "--> ALLERT! It was used only one crystal for the cog computation..."<<endl;
//        }

        clusterMap._COGcrySize = count;

        R /= sumW;
        Z /= sumW;

        clusterMap._cluCogRow = R;
        clusterMap._cluCogColumn = Z;


        res.setY( (R*2.+ 1.0)*cryHalfSize);//cg->nCrystalR()*cryHalfSize);//
        res.setZ( (Z*2.+1.0)*cryHalfSize);//cg->nCrystalZ()*cryHalfSize);//
        res.setX(-depth/*res.getX() + deltaXc + 3904.*/);//value used to shift in tracker coordinate system


        //        resError.setZ( RMScryHalfSize /TMath::Sqrt(cluster.size()) * TMath::Sqrt(tmpEq)/cluster.energyDep);

        //Save cog#Vector and its errorVector
        //        cluster.cog3VectorError = resError;
        //        cluster._impactPoint = res;


        clusterMap._cluCOG = res;
        //        return res;
}


CLHEP::Hep3Vector LOGcog(CaloCluster &cluster, double w, double depth){
        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;

        //if(! geom->hasElement<Calorimeter>() ) return;
        GeomHandle<Calorimeter> cg;
        CLHEP::Hep3Vector res, resError;
        //int tmpZ = 0, tmpR = 0, vane = -1;

        //get crystal's geometrical information
        double cryHalfSize      =      cg->crystalHalfSize();
        //double cryHalfLength    =    cg->crystalHalfLength();

        //using the expression for the RMS of a flat distribution, I calculate the value of the RMS for the face of the crystals
        //double RMScryHalfSize   =   2.*cryHalfSize / TMath::Sqrt(12.);
        //double RMScryHalfLenght = 2.*cryHalfLength / TMath::Sqrt(12.);
        //double defaultError     = 1e-06;//default from geant4

        double  tmpEq    = 0.;
        double  sumEi    = 0.;
        double offSet    = w;
        double weight    = 0.;
        double sumW      = 0.;

        bool isfirstCrystal=true;
        float R =0.;
        float Z =0.;

        int count = 0;
        for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){

                sumEi += (*itCD)->energyDep();
        }

        //clusterMap._cluSize = cluster.size();

        for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);

                //now I take the first crystal of the cluster as a point of reference for calculating the cog
                if(isfirstCrystal){
                        res = cg->getCrystalOriginByRO(thehit.id());
                        isfirstCrystal = false;
                        //tmpZ = cg->getCrystalZByRO(thehit.id());
                        //tmpR = cg->getCrystalRByRO(thehit.id());
                        //int vane = cg->getVaneByRO(thehit.id());
                        //clusterMap._vane = vane;
                }

                //Get Z and R from readout
                double tZ = cg->getCrystalZByRO(thehit.id());
                double tR = cg->getCrystalRByRO(thehit.id());

                //        clusterMap._rowVec.push_back(tR);
                //        clusterMap._columnVec.push_back(tZ);

                //Multiply Z and R for thhe crystal's energy, which is the weight we use in that algorithm
                weight = offSet + TMath::Log((*itCD)->energyDep()/sumEi);
                if(weight < 0.0){
                        weight = 0.0;
                }else{
                        count++;
                        //                clusterMap._COGrowVec.push_back(tR);
                        //                clusterMap._COGcolumnVec.push_back(tZ);
                }
                tZ *= weight;
                tR *= weight;

                Z += tZ;
                R += tR;
                sumW += weight;

                //Calculate the sum of the square of the weight for derive the error of the weighted mean
                tmpEq += std::pow(weight, 2.);
        }
//        if(count <=1){
//                cout<< "--> ALLERT! It was used only one crystal for the cog computation..."<<endl;
//        }

        //clusterMap._COGcrySize = count;

        R /= sumW;
        Z /= sumW;

//        clusterMap._cluCogRow = R;
//        clusterMap._cluCoglumn = Z;


        res.setY( (R*2.+ 1.0)*cryHalfSize);//cg->nCrystalR()*cryHalfSize);//
        res.setZ( (Z*2.+1.0)*cryHalfSize);//cg->nCrystalZ()*cryHalfSize);//
        res.setX(-depth/*res.getX() + deltaXc + 3904.*/);//value used to shift in tracker coordinate system


        //        resError.setZ( RMScryHalfSize /TMath::Sqrt(cluster.size()) * TMath::Sqrt(tmpEq)/cluster.energyDep);

        //Save cog#Vector and its errorVector
        //        cluster.cog3VectorError = resError;
        //        cluster._impactPoint = res;


        //        clusterMap._cluCOG = res;
        return res;
}

}
