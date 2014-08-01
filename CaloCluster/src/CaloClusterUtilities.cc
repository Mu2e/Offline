//
// General utilities for the calorimeter's studies
//
// $Id: CaloClusterUtilities.cc,v 1.10 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Original author G. Pezzullo & G. Tassielli & G. Onorato
//

#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "BaBar/BaBar/include/Constants.hh"

// C++ includes
#include<iostream>


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
        GeomHandle<VaneCalorimeter> cg;

        double Hsize = cg->caloGeomInfo().crystalHalfTrans();
        double Hleng = cg->caloGeomInfo().crystalHalfLength();
        double ROsize = cg->caloGeomInfo().roHalfTrans();

        cout << "Crystal HSize " << Hsize
                        << "\nCrystal HLeng " << Hleng
                        << "\nRO size " << ROsize << endl;

        for (int i=0; i<4; ++i) {

                Vane const& thevane = cg->vane(i);
                cout << "Vane " << i << " : "
                                << "\nOrigin: " << thevane.origin()
                                << "\nLocal origin: " << thevane.originLocal()
                                << "\nSize: " << thevane.size()
                                << "\nRotation: " << (thevane.rotation()) << endl;

        }

        int nRO = cg->nRO();

        for (int j=0; j< nRO/2; ++j) {

                int thevane = cg->vaneByRO(2*j);

                cout << "Crystal n. " << cg->crystalByRO(2*j);
                cout << "\tVane " << thevane;


                CLHEP::Hep3Vector cntr = cg->crystalOriginByRO(2*j);
                CLHEP::Hep3Vector Xaxis = cg->crystalAxisByRO(2*j);
                CLHEP::Hep3Vector Yaxis = cg->crystalAxisByRO(2*j).orthogonal();
                CLHEP::Hep3Vector Zaxis(0,0,1);

                cout << "\tcenter " << cntr
                                << "\tXaxis " << Xaxis
                                << "\tYaxis " << Yaxis;

                CLHEP::Hep3Vector toLeftC = ( (-Hsize) * Yaxis ) + ( (-Hleng) * Xaxis ) + ( (-Hsize) * Zaxis );

                CLHEP::Hep3Vector toRightC = ( (Hsize) * Yaxis ) + ( (Hleng+(2*ROsize)) * Xaxis ) + ( (Hsize) * Zaxis );

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

        GeomHandle<VaneCalorimeter> cg;
        _localCrystal  = cg->crystalByRO(_localRO);
        _localVane = cg->vaneByRO(_localRO);

        art::Handle<SimParticleCollection> simParticles;
        event.getByLabel(_g4ModuleLabel, simParticles);

        SimParticle const& sim = simParticles->at(track);

        CLHEP::Hep3Vector origin = sim.startPosition();

        _startingVane = getStartingVane(origin);

        _fromOutside = (_startingVane == -1);

        _primary = !sim.hasParent();

        _generated = sim.fromGenerator();

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
        GeomHandle<VaneCalorimeter> cg;

        for (size_t i=0; i<cg->nVane(); ++i) {

                Vane const & vane = cg->vane(i);
                CLHEP::Hep3Vector rsize = (vane.rotation()) * vane.size();
                CLHEP::Hep3Vector vaneOr = vane.origin();


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
        GeomHandle<VaneCalorimeter> cg;
        double size = 2.0*cg->caloGeomInfo().crystalHalfTrans();//[mm]
        int index = (int) (val / size);
        double result =  val - index*size;

        return result;
}
double indexToCoor(double ind){
        GeomHandle<VaneCalorimeter> cg;
        double res = ind*2.0*cg->caloGeomInfo().crystalHalfTrans();
        res += cg->caloGeomInfo().crystalHalfTrans();
        return res;
}
double indexToCoor(int ind){
        GeomHandle<VaneCalorimeter> cg;
        double res = (double)(ind*2.0*cg->caloGeomInfo().crystalHalfTrans());
        res += cg->caloGeomInfo().crystalHalfTrans();
        return res;
}

//the following method set as the cog of the Cluster the cog which has only the w coordinate
//corrected using the information of the position of the most energetic crystal
void w_correction_0(double& clCOGw,double& clCOGwErr, int& clCryEnergyMaxColumn){
        GeomHandle<VaneCalorimeter> cg;
        double x = clCOGw - indexToCoor(clCryEnergyMaxColumn), xErr = clCOGwErr;
        xErr -= (cg->caloGeomInfo().crystalHalfTrans()/sqrt(12.0));

        double val = 0.0, valErr = 0.0;

        if(x<-50.) x=-50.;
        if(x>50.) x=50.;

        double p[10] = {-19.48, 1.178, -0.0223, -0.0001045, 3.044e-6, 4.378e-08, 4.862e-10, -1.462e-11, -1.134e-13, 1.403e-15};
        for(int j=0; j<9; ++j){
                val += p[j]*pow(x, j);
        }
        int exp = 0;
        for(int k=1; k<9; ++k){
                exp = k-1;
                valErr += k*p[k]*pow(xErr,exp);
        }
        valErr +=(cg->caloGeomInfo().crystalHalfTrans()/sqrt(12.0));
        clCOGwErr = valErr;

        clCOGw = val + indexToCoor(clCryEnergyMaxColumn);

}

void v_correction_0( float& extrapolThetaV,  double& clCOGv,  double& clCOGvErr){// the first argument is in [deg]
        GeomHandle<VaneCalorimeter> cg;
        double val = 0.0, valErr = 0.0;
        double x = extrapolThetaV;
        if(x > 100.)  x = 100.;
        if(x < -100.) x = -100.;

        const int d = 5;

        double p[d] = {0.9204, 0.4834, -0.0003745, -4.387e-5, -4.722e-7};

        for(int j=0; j<d; ++j){
                val += p[j]*pow(x, j);
        }

        for(int j=1; j<d; ++j){
                valErr += j*p[j]*pow(clCOGvErr, j-1);
        }
        valErr += (cg->caloGeomInfo().crystalHalfTrans()/sqrt(12.0));
        val += clCOGv;
        clCOGv = val;
        clCOGvErr = valErr;

}

void v_correction_0( double& extrapolThetaV,  double& clCOGv,  double& clCOGvErr){// the first argument is in [deg]
        GeomHandle<VaneCalorimeter> cg;
        double val = 0.0, valErr = 0.0;
        double x = extrapolThetaV;
        if(x > 60.)  x = 60.;
        if(x < -80.) x = -80.;

        const int d = 5
                        , de = 10;

        double p[d] = {0.9204, 0.4834, -0.0003745, -4.387e-5, -4.722e-7};

        double pe[de] = {6.926, -0.003834, 0.0005851, 1.217e-5, -1.422e-7, -2.188e-8, -1.26e-10, 6.477e-12, 8.272e-14, 2.32e-16 };

        for(int j=0; j<d; ++j){
                val += p[j]*pow(x, j);
        }

        for(int j=0; j<de; ++j){
                valErr += pe[j]*pow(x, j);
        }

        //valErr += (cg->caloGeomInfo().crystalHalfTrans()/sqrt(12.0));
        valErr *= valErr;
        valErr += clCOGvErr*clCOGvErr;
        valErr = std::sqrt(valErr);

        val += clCOGv;

        clCOGv = val;

        clCOGvErr = valErr;

}

void w_correction_1(double& clCOGw,double& clCOGwErr, int& wSize){
        GeomHandle<VaneCalorimeter> cg;
        double p[10] = {-36.49, 34.92, -16.15, 3.546, -0.1358, -0.06925, 0.006722, 0.001008, -0.0001886, 8.151e-6};
        double spo = 0.0;
        if(wSize < 1.0) wSize = 1.0;
        if(wSize > 10.0) wSize = 10.0;

        for(int i=0; i<10; ++i){
                spo += p[i]*pow(wSize, i);
        }
        clCOGw -= spo;//the minus is not an error!

        double xErr = (cg->caloGeomInfo().crystalHalfTrans()/sqrt(12.0))/wSize;

        double  valErr = 0.0;

        int exp = 0;
        for(int k=1; k<9; ++k){
                exp = k-1;
                valErr += k*p[k]*pow(xErr,exp);
        }
        clCOGwErr += valErr;
}

void v_correction_1(int& clCryEnergyMaxRow,  double& clCOGv,  double& clCOGvErr){
        GeomHandle<VaneCalorimeter> cg;
        double val = 0.0, valErr = 0.0;
        double x = clCOGv - indexToCoor(clCryEnergyMaxRow);
        if(x > 100.)  x = 100.;
        if(x < -100.) x = -100.;

        if( x < -50.) x=-50.;
        if( x > 30.) x = 30.;

        double p[4] = {-0.2246, 0.9762, 0.001116, 7.389e-5};

        for(int j=0; j<4; ++j){
                val += p[j]*pow(x, j);
        }

        for(int j=1; j<4; ++j){
                valErr += j*p[j]*pow(clCOGvErr, j-1);
        }
        valErr += (cg->caloGeomInfo().crystalHalfTrans()/sqrt(12.0));
        val += indexToCoor(clCryEnergyMaxRow);
        clCOGv = val;
        clCOGvErr = valErr;

}


void cog_correction_0(CaloCluster &cluster){
        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        if(! geom->hasElement<VaneCalorimeter>() ) return;
        GeomHandle<VaneCalorimeter> cg;
        CLHEP::Hep3Vector res(1., 1., 1.);

        CLHEP::Hep3Vector  resError(1e-1, 1e-1, 1e-1);

        //get crystal's geometric information
        double cryHalfSize      =      cg->caloGeomInfo().crystalHalfTrans();
        double  tmpEq    = 0.;

        float R =0.;
        float Z =0., RZ = 0., RQ = 0.;

        int c = 0;
        for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){
                ++c;
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);
                //now I take the first crystal of the cluster as a point of reference for calculating the cog

                //Get Z and R from readout
                double tZ = cg->crystalZByRO(thehit.id());
                double tR = cg->crystalRByRO(thehit.id());

                Z += tZ;
                RQ += pow(tR, 2);
                R += tR;
                RZ += tZ*tR;

                //Multiply Z and R for the crystal's energy, which is the weight we use in that algorithm
                tZ *=(*itCD)->energyDep();
                tR*=(*itCD)->energyDep();

                //Calculate the sum of the square of the weight for derive the error of the weighted mean
                tmpEq += std::pow((*itCD)->energyDep(), 2.);
        }
        float m = R*Z - RZ;
        m /= (pow(R,2) - RQ);

        R/=cluster.energyDep();
        Z/=cluster.energyDep();

        res.setY( (R*2.+ 1.0)*cryHalfSize);
        res.setZ( (Z*2.+1.0)*cryHalfSize);
        res.setX(0.0);

}

//the following procedure fill the CaloCluster object with the COGVector, and its error vector
// with the respective values obtained using an energy weighted algorithm
void cog(CaloCluster &cluster){
        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        if(! geom->hasElement<VaneCalorimeter>() ) return;
        GeomHandle<VaneCalorimeter> cg;
        CLHEP::Hep3Vector res(1., 1., 1.);

        CLHEP::Hep3Vector  resError(1e-1, 1e-1, 1e-1);

        //get crystal's geometrical information
        double cryHalfSize      =      cg->caloGeomInfo().crystalHalfTrans();

        double  tmpEq    = 0.;
        double  showerDepth = 0.0;

        bool isfirstCrystal=true;
        float R =0.;
        float Z =0., RZ = 0., RQ = 0.;

        int c = 0;
        for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){
                ++c;
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);
                //now I take the first crystal of the cluster as a point of reference for calculating the cog
                if(isfirstCrystal){
                        res = cg->crystalOriginByRO(thehit.id());
                        isfirstCrystal = false;
                }

                //Get Z and R from readout
                double tZ = cg->crystalZByRO(thehit.id());
                double tR = cg->crystalRByRO(thehit.id());

                tZ *=(*itCD)->energyDep();
                tR*=(*itCD)->energyDep();

                Z += tZ;
                RQ += pow(tR, 2);
                R += tR;
                RZ += tZ*tR;

                //Calculate the sum of the square of the weight for derive the error of the weighted mean
                tmpEq += std::pow((*itCD)->energyDep(), 2.);
        }
        float m = R*Z - RZ;
        m /= (pow(R,2) - RQ);

        R/=cluster.energyDep();
        Z/=cluster.energyDep();

        res.setY( (R*2.+ 1.0)*cryHalfSize);
        res.setZ( (Z*2.+1.0)*cryHalfSize);
        res.setX(- showerDepth);

        double error = cryHalfSize*2.0;
        error /= cluster.energyDep();
        error /= sqrt(12.0);

        resError.setX(error);
        resError.setY(error);
        resError.setZ(error);

        cluster.SetCogRow((int) R);
        cluster.SetCogColumn((int) Z);
        cluster.SetCog3Vector(res);
        cluster.SetCog3VectorError(resError);
}


void cog_depth(CaloCluster &cluster, double depth, ClusterMap &clusterMap){
        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        GeomHandle<VaneCalorimeter> cg;
        CLHEP::Hep3Vector res(1., 1., 1.);

        CLHEP::Hep3Vector  resError(1e-1, 1e-1, 1e-1);

        //get crystal's geometrical information
        double cryHalfSize      =      cg->caloGeomInfo().crystalHalfTrans();

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
                        res = cg->crystalOriginByRO(thehit.id());
                        isfirstCrystal = false;
                        int vane = cg->vaneByRO(thehit.id());
                        clusterMap._vaneId = vane;

                }

               //Get Z and R from readout
                double tZ = cg->crystalZByRO(thehit.id());
                double tR = cg->crystalRByRO(thehit.id());

                V += tR;
                W += tZ;
                VW += tR*tZ;
                VQ += pow(tR,2);

                clusterMap._rowVec.push_back(tR);
                clusterMap._columnVec.push_back(tZ);
                clusterMap._cryEdepVec.push_back((*itCD)->energyDep());
                clusterMap._timeVec.push_back((*itCD)->time());

                clusterMap._COGrowVec.push_back(tR);
                clusterMap._COGcolumnVec.push_back(tZ);

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
        float m = V*W - size*VW;
        m /= (pow(V,2) - size*VQ);


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

        res.setY( (R*2.+ 1.0)*cryHalfSize);
        res.setZ( (Z*2.+ 1.0)*cryHalfSize);
        res.setX(- showerDepth);

        clusterMap._cluCOG = res;
        clusterMap._showerDir = m;
        clusterMap._errShowerDir = errM;

}


//on the following we implement an algorithm for the cog which uses the logarithm of the energy as weight (w_{i}). This is not correct, as references show, we need to calculate from simulation an offset to add at each w_{i}
//void LOGcog(CaloCluster &cluster){
//CLHEP::Hep3Vector LOGcog(CaloCluster &cluster, double w, double depth){
void LOGcogMap(CaloCluster &cluster, double w, double depth, ClusterMap &clusterMap ){

        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;

        GeomHandle<VaneCalorimeter> cg;
        CLHEP::Hep3Vector res, resError;

        //get crystal's geometric information
        double cryHalfSize      =      cg->caloGeomInfo().crystalHalfTrans();

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
                        res = cg->crystalOriginByRO(thehit.id());
                        isfirstCrystal = false;
                        int vane = cg->vaneByRO(thehit.id());
                        clusterMap._vaneId = vane;
                }

                //Get Z and R from readout
                double tZ = cg->crystalZByRO(thehit.id());
                double tR = cg->crystalRByRO(thehit.id());

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

        clusterMap._COGcrySize = count;

        R /= sumW;
        Z /= sumW;

        clusterMap._cluCogRow = R;
        clusterMap._cluCogColumn = Z;


        res.setY( (R*2.+ 1.0)*cryHalfSize);
        res.setZ( (Z*2.+1.0)*cryHalfSize);
        res.setX(-depth);

        clusterMap._cluCOG = res;
}


CLHEP::Hep3Vector LOGcog(CaloCluster &cluster, double w, double depth){
        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;

        GeomHandle<VaneCalorimeter> cg;
        CLHEP::Hep3Vector res, resError;

        //get crystal's geometrical information
        double cryHalfSize      =      cg->caloGeomInfo().crystalHalfTrans();

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

        for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);

                //now I take the first crystal of the cluster as a point of reference for calculating the cog
                if(isfirstCrystal){
                        res = cg->crystalOriginByRO(thehit.id());
                        isfirstCrystal = false;

                }

                //Get Z and R from readout
                double tZ = cg->crystalZByRO(thehit.id());
                double tR = cg->crystalRByRO(thehit.id());

                //Multiply Z and R for thhe crystal's energy, which is the weight we use in that algorithm
                weight = offSet + TMath::Log((*itCD)->energyDep()/sumEi);
                if(weight < 0.0){
                        weight = 0.0;
                }else{
                        count++;
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

        res.setY( (R*2.+ 1.0)*cryHalfSize);
        res.setZ( (Z*2.+1.0)*cryHalfSize);
        res.setX(-depth);

        return res;
}

std::vector<int> vecRoId(CaloCluster &cluster){

    std::vector<int> vec;
      
     for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
                CaloHit const& thehit = *ROIds.at(0);

		vec.push_back(thehit.id());
     }
     
     return vec;
  }

}
