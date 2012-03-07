//
// implementation of different algorithm to reconstruct the impact position on the electrons on the calorimeter
//
// $Id: CaloClusterCog_module.cc,v 1.3 2012/03/07 18:00:38 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/03/07 18:00:38 $
//
// Original author G. Pezzullo & G. Tassielli
//


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "Mu2eUtilities/inc/LinePointPCA.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//calorimeter packages
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
//#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
//#include "CaloCluster/inc/CaloClusterFinder.hh"
//#include "CaloCluster/inc/ClosestCaloClusterFinder.hh"

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"


// Other includes.
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib/exception.h"

//#include "Analyses/inc/MCCaloUtilities.hh"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"

//root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"

// From the art tool-chain
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "TMath.h"

//#include "CaloCluster/inc/CaloClusterUtilities.hh"

using namespace std;

namespace mu2e {


struct elecData{
        double _cluEnergy;
        double _cluTime;
        double _impTime;
        double _impEnergy;
        CLHEP::Hep3Vector _cluCog;
        CLHEP::Hep3Vector _impPos;
        CLHEP::Hep3Vector _impPosCryFrame;
        unsigned int _vane;
        unsigned int _row;
        unsigned int _colum;
        CLHEP::Hep3Vector _impMom3Vec;


        bool operator<( const elecData other) const{
                return ( _impTime< other._impTime);
        }
        elecData & operator=(const elecData& other) {
                _cluEnergy = other._cluEnergy;
                _cluTime   = other._cluTime;
                _impTime = other._impTime;
                _impEnergy = other._impEnergy;
                _cluCog    = other._cluCog;
                _impPos = other._impPos;
                _impPosCryFrame = other._impPosCryFrame;
                _vane     = other._vane;
                _row      = other._row;
                _colum    = other._colum;
                _impMom3Vec =other._impMom3Vec;

                return *this;
        }
        elecData():
                _impTime(1e10),
                _impEnergy(0.0){
        }
};

//the key is the the vane
typedef std::map<unsigned int,std::map<unsigned int, elecData > > ElecMap;


static int ncalls(0);

class CaloClusterCog : public art::EDAnalyzer {
public:
        explicit CaloClusterCog(fhicl::ParameterSet const& pset):
        _diagLevel(pset.get<int>("diagLevel",0)),
        _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
        _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
        _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel", "makeCaloCluster")),
        _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
        _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "time")),
        _producerName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
        _elextractModuleLabel(pset.get<std::string>("elextractModuleLabel", "extractElData")),
        _extractElectronsData(pset.get<string>("elextractModuleLabel")),
        _rowToCanc(pset.get<int>("rowToCancel",-1)),//row index running from 1, 2, ... nCryR
        _columnToCanc(pset.get<int>("columnToCancel",-1)),//column index running from 1, 2, ... nCryZ
        _depthPitch(pset.get<double>("depthPitch",0.5)),//[mm]
        _thetaVimpact(pset.get<double>("thetaVimpact",0.785)),//45 [grad]
        _thetaWimpact(pset.get<double>("thetaWimpact",0.785)),//45 [grad]
        _nAnalyzed(0),
        //_hTHistGlobalEffNorm(0),
        _hTHistDeltaEnergy(0),
        _hTHistDeltaEnergVRec(0),
        //_hTHistGlobalEffRec(0),
        _hTHistDeltaUquality(0),
        _hTHistDeltaURec(0),
        _hTHistDeltaVquality(0),
        _hTHistDeltaVRec(0),
        _hTHistDeltaWquality(0),
        _hTHistDeltaWRec(0),
        _hTHistImpactParam(0),
        _hTHistDistanceToCog(0),
        _hTHistMomDotVaneNorm(0),
        _hTHistRecImpactParam(0),
        _hTHistRecDistanceToCog(0),
        _hTHistMomRecDotVaneNorm(0),
        _CogOffSet(pset.get<double>("cogOffset",3.)),
        _MaxDepth(pset.get<double>("maxDepth",11.)),
        _application(0),
        _directory(0)
        {
        }
        virtual ~CaloClusterCog() {
        }
        virtual void beginJob();
        virtual void endJob();

        void analyze(art::Event const& e );

private:

        void doCalorimeter(art::Event const& evt, bool skip);

        // Diagnostic level
        int _diagLevel;

        // Label of the generator.
        std::string _generatorModuleLabel;

        // Label of the G4 module
        std::string _g4ModuleLabel;

        // Label of the calo readout hits maker
        std::string _caloReadoutModuleLabel;

        // Label of the calo crystal hists maker
        std::string _caloCrystalModuleLabel;

        // Label of the calo clusters  maker
        std::string _caloClusterModuleLabel;

        string _caloClusterAlgorithm;
        string _caloClusterSeeding;
        string _producerName;
        string _elextractModuleLabel;
        string _extractElectronsData;

        double _minimumEnergy; //minimum energy deposition of hits
        int _rowToCanc, _columnToCanc;
        double _depthPitch;
        double _thetaVimpact;
        double _thetaWimpact;
        //number of analyzed events
        int _nAnalyzed;

        //TH1D* _hTHistGlobalEffNorm;
        TH1D* _hTHistDeltaEnergy;
        TH1D* _hTHistDeltaEnergVRec;
        //TH1D* _hTHistGlobalEffRec;

        TH2D* _hTHistDeltaUquality;
        TH2D* _hTHistDeltaURec;
        TH2D* _hTHistDeltaVquality;
        TH2D* _hTHistDeltaVRec;
        TH2D* _hTHistDeltaWquality;
        TH2D* _hTHistDeltaWRec;
        TH2D* _hTHistImpactParam;
        TH2D* _hTHistDistanceToCog;
        TH2D* _hTHistMomDotVaneNorm;
        TH2D*_hTHistRecImpactParam;
        TH2D*_hTHistRecDistanceToCog;
        TH2D*_hTHistMomRecDotVaneNorm;

        double _CogOffSet;
        double _MaxDepth;

        double globalNtrkCut;
        double globalNtrkTot;
        double *globalCaloCut;
        double *globalRecCaloCut;

        std::auto_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;

        // The job needs exactly one instance of TApplication.  See note 1.
        auto_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;


};

bool findTrkId(std::vector<unsigned int> vec, unsigned int t){
        bool res = false;

        unsigned int size = vec.size();
        if(size!=0){
                unsigned int cont = 0;
                while(/*!res ||*/ cont!=size){
                        if(vec[cont] == t) {
                                res = true;
                        }
                        ++cont;

                }
        }
        return res;
}


void CaloClusterCog::beginJob( ) {

        cout << "start CaloClusterCog..."<<endl;

        CaloManager = auto_ptr<MCCaloUtilities>(new MCCaloUtilities());

        // If needed, create the ROOT interactive environment. See note 1.
        if ( !gApplication ){
                int    tmp_argc(0);
                char** tmp_argv(0);
                _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
        }

        gStyle->SetPalette(1);
        gROOT->SetStyle("Plain");

        _directory = gDirectory;


}

void CaloClusterCog::analyze(art::Event const & evt ) {

        ++_nAnalyzed;
        ++ncalls;

        CaloClusterer c;

        art::ServiceHandle<GeometryService> geom;
        if (ncalls == 1) {

                // cout << "This should be done only in the first event" << endl;

                art::ServiceHandle<art::TFileService> tfs;
                art::TFileDirectory cogResol = tfs->mkdir("CogResol");
                art::TFileDirectory energyResol = tfs->mkdir("EnergyResol");
                art::TFileDirectory dcaCog = tfs->mkdir("DcaCog");

                _hTHistDeltaEnergy   = energyResol.make<TH1D>( "DeltaEnergy", "DeltaEnergy;Eseed-Eclu [MeV]", 100, -20., 110.);
                _hTHistDeltaEnergVRec= energyResol.make<TH1D>( "DeltaEnergVRec", "DeltaEnergVRec;Eseed-Eclu [MeV]", 100, -20., 110.);

                _hTHistDeltaUquality = cogResol.make<TH2D>( "DeltaUquality", "DeltaUquality; depth_{MaxShower} [mm];Xseed-Xclu [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 200, -100., 100. );
                _hTHistDeltaURec     = cogResol.make<TH2D>( "DeltaURec", "DeltaURec; depth_{MaxShower} [mm];Xseed-Xclu [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 200, -100., 100.);
                _hTHistDeltaVquality = cogResol.make<TH2D>( "DeltaVquality", "DeltaVquality; depth_{MaxShower} [mm];Yseed-Yclu [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 200, -100., 100.);
                _hTHistDeltaVRec     = cogResol.make<TH2D>( "DeltaVRec", "DeltaVRec; depth_{MaxShower} [mm]; Yseed-Yclu [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 200, -100., 100.);
                _hTHistDeltaWquality = cogResol.make<TH2D>( "DeltaWquality", "DeltaWquality; depth_{MaxShower} [mm];vZseed-Zclu [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 200, -100., 100.);
                _hTHistDeltaWRec     = cogResol.make<TH2D>( "DeltaWRec", "DeltaWRec; depth_{MaxShower} [mm];Zseed-Zclu [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth,200, -100., 100.);

                _hTHistImpactParam   = dcaCog.make<TH2D>( "ImpactParam", "ImpactParam;  depth_{MaxShower} [mm];ImpactParam [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 600,  0., 150.);
                _hTHistDistanceToCog = dcaCog.make<TH2D>( "DistanceToCog", "DistanceToCog; depth_{MaxShower} [mm]; |#vec{x}_{seed} - #vec{x}_{cog}| [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth,  600,  0., 150.);
                _hTHistMomDotVaneNorm= dcaCog.make<TH2D>( "MomDotVaneNorm","MomDotVaneNorm ; depth_{MaxShower} [mm] ;ArcCos(#vec{p}, #hat{n}) / |#vec{p}| ", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 720 , -180., 180.0);

                _hTHistRecImpactParam   = dcaCog.make<TH2D>( "RecImpactParam", "RecImpactParam;  depth_{MaxShower} [mm];ImpactParam [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 400,  0., 100.);
                _hTHistRecDistanceToCog = dcaCog.make<TH2D>( "RecDistanceToCog", "RecDistanceToCog; depth_{MaxShower} [mm]; |#vec{x}_{seed} - #vec{x}_{cog}| [mm]", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 400,  0., 100.);
                _hTHistMomRecDotVaneNorm= dcaCog.make<TH2D>( "MomRecDotVaneNorm","MomRecDotVaneNorm ; depth_{MaxShower} [mm]; ArcCos(#vec{p}, #hat{n}) / |#vec{p}| ", (_MaxDepth - 0.) / _depthPitch, 0., _MaxDepth, 720 , -180., 180.0);
        }

        doCalorimeter(evt, _skipEvent);

} // end of analyze

void CaloClusterCog::endJob() {

}



void CaloClusterCog::doCalorimeter(art::Event const& evt, bool skip){


        if ( _diagLevel > 0 ) cout << "MakeCaloCluster: produce() begin" << endl;

        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        if(! geom->hasElement<Calorimeter>() ) return;
        GeomHandle<Calorimeter> cg;

        // Get handles to calorimeter collections
        art::Handle<CaloHitCollection> caloHits;
        evt.getByLabel(_caloReadoutModuleLabel, caloHits);

        art::Handle<CaloCrystalHitCollection>  caloCrystalHits;
        evt.getByLabel(_caloCrystalModuleLabel, caloCrystalHits);

        // Get the persistent data about pointers to StepPointMCs
        art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
        evt.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);

        PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();
        if (!( caloHits.isValid())) {
                return;
        }

        if (!caloCrystalHits.isValid()) {
                cout << "NO CaloCrystalHits" << endl;
                return;
        }

        // Get handles to the generated and simulated particles.
        art::Handle<GenParticleCollection> genParticles;
        evt.getByLabel(_generatorModuleLabel, genParticles);

        art::Handle<SimParticleCollection> simParticles;
        evt.getByLabel(_g4ModuleLabel, simParticles);

        art::Handle<CaloClusterCollection> caloClusters;
        evt.getByLabel(_caloClusterModuleLabel,_producerName,caloClusters );

        auto_ptr<CaloClusterCollection>             caloClustersPointer(new CaloClusterCollection);

        art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
        evt.getByLabel(_extractElectronsData,genEltrksHandle);

        VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
        std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

        double trkMomCut = 100.0;//MeV
        std::vector<unsigned int> tmpV, tmpVTot;
        int NtrkCut =0;
        int NtrkTot = 0;

        //mapping &counting the electrons with quality cuts in the TRK
        for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
                ++NtrkTot;
                VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
                GenElHitData& hdil = iEltrk.getithLoopHit(0);

                if(!findTrkId(tmpVTot, iEltrk.getTrkID().asUint() ) ){

                        //cout<<"faccio il puschback in tmpV..."<<endl;
                        tmpVTot.push_back( iEltrk.getTrkID().asUint() );
                        //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                }

                if (iEltrk.getNumOfHit()>=20){



                        if(hdil._hitMomentum.mag() >= trkMomCut){
                                NtrkCut++;// # of electrons at the entrance of the TRK which have trkMomCut MeV of kinetic energy and which will light at least 20 straws
                                if(!findTrkId(tmpV, iEltrk.getTrkID().asUint() ) ){

                                        //cout<<"faccio il puschback in tmpV..."<<endl;
                                        tmpV.push_back( iEltrk.getTrkID().asUint() );
                                        //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                                }

                                //                                for(unsigned int b=1; b <=(unsigned int) _hTHistEff->GetNbinsX(); ++b){
                                //
                                //                                        if (hdil._hitMomentum.mag() >= (_hTHistEff->GetXaxis()->GetBinCenter(b) -  0.5*_hTHistEff->GetXaxis()->GetBinWidth(b) )){
                                //                                               // NGoodElec[b-1]++;
                                //                                                if ( iEltrk.isConversionEl() ) {
                                //                                                     //   NGoodConvEl[b-1]++;
                                //                                                }
                                //                                        }
                                //
                                //                                }
                        }

                }

        }//end TRK mapping
        //cout <<"NtrkCut = " <<NtrkCut<<endl;

        if (NtrkCut==0) return;

        globalNtrkCut += NtrkCut;
        globalNtrkTot += NtrkTot;

        for(int bin =1; bin <= _hTHistDeltaURec->GetNbinsX(); ++bin){
                std::vector<unsigned int> trkVec, trkVecTot;
                for(unsigned int f = 0; f != tmpV.size(); ++f){
                        trkVec.push_back(tmpV[f]);
                }
                for(unsigned int f = 0; f != tmpVTot.size(); ++f){
                        trkVecTot.push_back(tmpVTot[f]);
                }
                ElecMap elecMap;

                if(caloClusters->size()>0 ){
                        int iVane;
                        for(size_t icl=0; icl<caloClusters->size(); ++icl){
                                CaloCluster const& clu = (*caloClusters).at(icl);

                                caloClustersPointer->push_back(clu);
                                CaloClusterCollection::iterator tmpCluster = caloClustersPointer->end();
                                tmpCluster--;


                                CLHEP::Hep3Vector cogDepth = cog_depth(*tmpCluster,_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) );//- 0.5*_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin) );

                                iVane = clu.vaneId;

                                CaloCrystalHitPtrVector caloClusterHits = clu.caloCrystalHitsPtrVector;

                                for(size_t i=0; i<caloClusterHits.size(); ++i){
                                        CaloCrystalHit const& hit = *(caloClusterHits.at(i));

                                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                                        if(ROIds.size()<1 ) continue;

                                        CaloHit const& thehit = *ROIds.at(0);
                                        size_t collectionPosition = ROIds.at(0).key();

                                        PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
                                        size_t nHitsPerCrystal = mcptr.size();

                                        for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

                                                //cout<< "Start loop..."<< "j2 = "<< j2<<endl;

                                                StepPointMC const& mchit = *mcptr[j2];

                                                // The simulated particle that made this hit.
                                                SimParticleCollection::key_type trackId(mchit.trackId());
                                                CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame(thehit.id(), mchit.position());

                                                if( (cg->getCrystalRByRO(thehit.id()) != ( _rowToCanc - 1 ) && cg->getCrystalZByRO(thehit.id()) != ( _columnToCanc - 1 ) ) ){

                                                        if(elecMap[iVane][trackId.asUint()]._impTime > mchit.time() ){
                                                                elecMap[iVane][trackId.asUint()]._cluEnergy = clu.energyDep;
                                                                elecMap[iVane][trackId.asUint()]._cluTime = clu.time;
                                                                elecMap[iVane][trackId.asUint()]._impTime = mchit.time();
                                                                elecMap[iVane][trackId.asUint()]._impEnergy = mchit.momentum().mag();
                                                                elecMap[iVane][trackId.asUint()]._cluCog = cogDepth;// clu._impactPoint;//logCog;
                                                                elecMap[iVane][trackId.asUint()]._impPos = mchit.position();
                                                                elecMap[iVane][trackId.asUint()]._impPosCryFrame = cryFrame;
                                                                elecMap[iVane][trackId.asUint()]._vane = iVane;
                                                                elecMap[iVane][trackId.asUint()]._row = cg->getCrystalRByRO(thehit.id());
                                                                elecMap[iVane][trackId.asUint()]._colum = cg->getCrystalZByRO(thehit.id());
                                                                elecMap[iVane][trackId.asUint()]._impMom3Vec = mchit.momentum();

                                                                if(_diagLevel < 0){
                                                                cout<< "###################"<<endl;
                                                                cout<< "idVande = "<< iVane<<endl;
                                                                cout << "cluU = "<<elecMap[iVane][trackId.asUint()]._impPos.getX()<<endl;
                                                                cout << "cluV = "<<elecMap[iVane][trackId.asUint()]._impPos.getY()<<endl;
                                                                cout << "cluW = "<<elecMap[iVane][trackId.asUint()]._impPos.getZ()<<endl;
                                                                cout<< "###################"<<endl;
                                                                }

                                                        }
                                                }
                                        }

                                }


                        }//end for(caloClsuters.size)
                }//end if(caloClsuters.size() > 0)

                //counting how many reconstructed clusters match with the generated particles
                unsigned int canc2 = 0;

                //cout <<"-------------------> numero di elettroni generati = "<< trkVecTot.size() <<endl;
                unsigned int size2 = trkVecTot.size();
                unsigned int it2=0;
                while( it2 < size2){
                        ElecMap::iterator ite = elecMap.begin();
                        bool trovato = false;
                        while(!trovato && ite!=elecMap.end() ){
                                if(ite->second.find(trkVecTot[it2]) != ite->second.end()){

                                        //                                        cout << "delta X = " << ite->second[trkVec[it2]]._impPos.getX() - ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                        //                                        cout << "delta Y = " << ite->second[trkVec[it2]]._impPos.getY() - ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                        //                                        cout << "delta Z = " << ite->second[trkVec[it2]]._impPos.getZ() - ite->second[trkVec[it2]]._cluCog.getZ()<<endl;



                                        CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, ite->second[trkVecTot[it2]]._impPos);


                                        CLHEP::Hep3Vector dirMom = ite->second[trkVecTot[it2]]._impMom3Vec.unit();
                                        if(_diagLevel < 0){
                                                cout<<"vaneId = "<<ite->first<<endl;
                                                cout<<"posX = "<<ite->second[trkVecTot[it2]]._impPos.getX()<<", posY = "<<ite->second[trkVecTot[it2]]._impPos.getY()<<", posZ = "<< ite->second[trkVecTot[it2]]._impPos.getZ()<<endl;
                                                cout<<"on vane ref..."<<endl;
                                                cout<<"posU = "<<vaneFrame.getX()<<", posV = "<<vaneFrame.getY()<<", posW = "<< vaneFrame.getZ()<<endl;
                                                cout << "dirMomX = "<< dirMom.getX()<<", dirMomY = "<< dirMom.getY()<<", dirMomZ = "<<dirMom.getZ() <<endl;
                                        }
                                        Vane const &vane = cg->getVane(ite->first);
                                        CLHEP::Hep3Vector dirMom_rotated = *(vane.getRotation())*dirMom;

                                        if(_diagLevel < 0){
                                        cout<<"after rotation..."<<endl;
                                        cout << "dirMomU = "<< dirMom_rotated.getX()<<", dirMomV = "<< dirMom_rotated.getY()<<", dirMomW = "<<dirMom_rotated.getZ() <<endl;
                                        cout<<"cog coordinates... on vane ref!"<<endl;
                                        cout<<"cluU = "<< ite->second[trkVecTot[it2]]._cluCog.getX()<< ", cluV = "<<ite->second[trkVecTot[it2]]._cluCog.getY() <<", cluW = "<<ite->second[trkVecTot[it2]]._cluCog.getZ()<<endl;
                                        }

                                        LinePointPCA lppca(vaneFrame,dirMom_rotated,  ite->second[trkVecTot[it2]]._cluCog);
//                                        LinePointPCA lppca(ite->second[trkVecTot[it2]]._impPos,dirMom,  ite->second[trkVecTot[it2]]._cluCog);

                                        CLHEP::Hep3Vector dcaV = lppca.unit();
                                        double impactParam = lppca.dca();


                                        double distanceToCog = sqrt( pow(vaneFrame.getX() -ite->second[trkVecTot[it2]]._cluCog.getX() ,2) + pow(vaneFrame.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() ,2) + pow(vaneFrame.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() ,2) );
//                                        double distanceToCog = sqrt( pow(ite->second[trkVecTot[it2]]._impPos.getX() -ite->second[trkVecTot[it2]]._cluCog.getX() ,2) + pow(ite->second[trkVecTot[it2]]._impPos.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() ,2) + pow(ite->second[trkVecTot[it2]]._impPos.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() ,2) );
                                        double tmpDepth = _hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin)*0.5;

                                        _hTHistRecDistanceToCog->Fill(tmpDepth, distanceToCog);
                                        _hTHistRecImpactParam->Fill(tmpDepth, impactParam);
                                        _hTHistMomRecDotVaneNorm->Fill(tmpDepth, asin(impactParam / distanceToCog)*180./TMath::Pi());
                                        if(_diagLevel < 0){
                                                cout<<"------------------------------------------------------"<<endl;
                                                cout<< "dcaV_U = "<<dcaV.getX()<< ", dcaV_V = "<<dcaV.getY()<<", dcaV_W = "<<dcaV.getZ()<<endl;
                                                cout << "impactParam = "<< impactParam<<endl;
                                                cout<< "distance ImpPointToCog = "<< distanceToCog << endl;
                                                cout<< "-----------------------------> vane = "<< ite->first<<endl;
                                                cout << "impX = "<< vaneFrame.getX()<<endl;
                                                cout << "impY = "<< vaneFrame.getY()<<endl;
                                                cout << "impZ = "<< vaneFrame.getZ()<<endl;
                                                cout<<"------------------------------------------------------"<<endl;
                                                cout << "cluX = "<<ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                                cout << "cluY = "<<ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                                cout << "cluZ = "<<ite->second[trkVec[it2]]._cluCog.getZ()<<endl;
                                                cout<<"------------------------------------------------------"<<endl;
                                        }
                                                                                        if(bin == 1){
                                                _hTHistDeltaEnergVRec->Fill( ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);


                                        }
                                        //CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, vane.getOrigin() );

                                        double deltaZ = (tmpDepth - impactParam/cos(_thetaWimpact) )*tan(_thetaWimpact);
                                        int sign = -1;
                                        if(dcaV.getY() < 0){
                                                sign = 1;
                                        }
                                        double deltaY = sign*(tmpDepth - impactParam/cos(_thetaVimpact) )*tan(_thetaVimpact);

                                        _hTHistDeltaURec->Fill(tmpDepth, vaneFrame.getX() - ite->second[trkVecTot[it2]]._cluCog.getX() );
                                        _hTHistDeltaWRec->Fill(tmpDepth, vaneFrame.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() + deltaZ );
                                        _hTHistDeltaVRec->Fill(tmpDepth, vaneFrame.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() + deltaY );
//                                        if(ite->first == 0 || ite->first == 2){
//
//                                                _hTHistDeltaURec->Fill(_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impPos.getX() - ite->second[trkVecTot[it2]]._cluCog.getX() );
//                                        }else{
//                                                _hTHistDeltaVRec->Fill(_hTHistDeltaVRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impPos.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() );
//                                        }
//                                        _hTHistDeltaWRec->Fill(_hTHistDeltaWRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impPos.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() );

                                        trovato = true;
                                        canc2 = it2;
                                }
                                ++ite;
                        }

                        if(trovato){
                                std::vector<unsigned int>::iterator er = trkVecTot.begin();
                                er +=canc2;
                                trkVecTot.erase(er);
                                size2 = trkVecTot.size();
                        }else{
                                ++it2;
                        }

                }
                if(_diagLevel < 0){
                        cout <<"after.... trkVecTot.size = "<< trkVecTot.size() <<endl;
                }

                //counting how many reconstructed clusters match with the generated particles which pass the quality cuts of the tracker
                unsigned int canc = 0;
                if(_diagLevel < 0){
                        cout <<"-------------------> number of quality generated electrons = "<< trkVec.size() <<endl;
                }
                unsigned int size = trkVec.size();
                unsigned int it=0;
                while( it < size){
                        ElecMap::iterator ite = elecMap.begin();
                        bool trovato = false;
                        while(!trovato && ite!=elecMap.end() ){
                                if(ite->second.find(trkVec[it]) != ite->second.end()){

                                        if(_diagLevel < 0){
                                                cout << "impX = "<<ite->second[trkVec[it]]._impPosCryFrame.getX()<<endl;
                                                cout << "impY = "<<ite->second[trkVec[it]]._impPosCryFrame.getY()<<endl;
                                                cout << "impZ = "<<ite->second[trkVec[it]]._impPosCryFrame.getZ()<<endl;
                                        }

                                        //globalCaloCut[bin - 1] += 1.0;

                                        trovato = true;
                                        canc = it;

                                        CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, ite->second[trkVec[it]]._impPos);

                                        CLHEP::Hep3Vector dirMom = ite->second[trkVec[it]]._impMom3Vec.unit();

                                        Vane const &vane = cg->getVane(ite->first);
                                        CLHEP::Hep3Vector dirMom_rotated = *(vane.getRotation())*dirMom;

                                        LinePointPCA lppca(vaneFrame,dirMom_rotated,  ite->second[trkVec[it]]._cluCog);

                                        CLHEP::Hep3Vector dcaV = lppca.unit();
                                        double impactParam = lppca.dca();

                                        double distanceToCog = sqrt( pow(vaneFrame.getX() -ite->second[trkVec[it]]._cluCog.getX() ,2) + pow(vaneFrame.getY() - ite->second[trkVec[it]]._cluCog.getY() ,2) + pow(vaneFrame.getZ() - ite->second[trkVec[it]]._cluCog.getZ() ,2) );

                                        _hTHistDistanceToCog->Fill(_hTHistDistanceToCog->GetXaxis()->GetBinCenter(bin) -_hTHistDistanceToCog->GetXaxis()->GetBinWidth(bin)*0.5, distanceToCog);
                                        _hTHistImpactParam->Fill(_hTHistImpactParam->GetXaxis()->GetBinCenter(bin) -_hTHistImpactParam->GetXaxis()->GetBinWidth(bin)*0.5, impactParam);
                                        _hTHistMomDotVaneNorm->Fill(_hTHistMomDotVaneNorm->GetXaxis()->GetBinCenter(bin) -_hTHistMomDotVaneNorm->GetXaxis()->GetBinWidth(bin)*0.5, asin(impactParam / distanceToCog)*180./TMath::Pi());

                                        if(bin ==1){
                                                _hTHistDeltaEnergy->Fill(ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);
                                                //_hTHistDeltaWquality->Fill(vaneFrame.getY() - ite->second[trkVec[it]]._cluCog.getY() );
                                        }
                                        //CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, vane.getOrigin() );
                                        double tmpDepth = _hTHistDeltaUquality->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaUquality->GetXaxis()->GetBinWidth(bin)*0.5;

                                        double impactParamPrjW = impactParam*(pow(dcaV.getX(), 2) + pow(dcaV.getZ(), 2));
                                        double impactParamPrjV = impactParam*(pow(dcaV.getX(), 2) + pow(dcaV.getY(), 2));

                                        double deltaW = (tmpDepth - impactParamPrjW/cos(_thetaWimpact) )*tan(_thetaWimpact);
                                        int sign = -1;
                                        if(dcaV.getY() < 0){
                                                sign = 1;
                                        }
                                        double deltaV = sign*(tmpDepth - impactParamPrjV/cos(_thetaVimpact) )*tan(_thetaVimpact);

                                        _hTHistDeltaUquality->Fill(tmpDepth, vaneFrame.getX() - ite->second[trkVec[it]]._cluCog.getX() );
                                        _hTHistDeltaWquality->Fill(tmpDepth, vaneFrame.getZ() - ite->second[trkVec[it]]._cluCog.getZ() + deltaW );
                                        _hTHistDeltaVquality->Fill(tmpDepth, vaneFrame.getY() - ite->second[trkVec[it]]._cluCog.getY() + deltaV );
                                      //  _hTHistDeltaUquality->Fill(_hTHistDeltaUquality->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaUquality->GetXaxis()->GetBinWidth(bin)*0.5, vaneFrame.getX() - ite->second[trkVec[it]]._cluCog.getX() );
                                       // _hTHistDeltaWquality->Fill(_hTHistDeltaWquality->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaWquality->GetXaxis()->GetBinWidth(bin)*0.5, vaneFrame.getZ() - ite->second[trkVec[it]]._cluCog.getZ() );

                                }
                                ++ite;
                        }

                        if(trovato){
                                std::vector<unsigned int>::iterator er = trkVec.begin();
                                er +=canc;
                                trkVec.erase(er);
                                size = trkVec.size();
                        }else{
                                ++it;
                        }

                }

                if(_diagLevel < 0){
                        cout <<"after.... trkVec.size = "<< trkVec.size() <<endl;
                        cout<<"...\n...\n "<<endl;
                }

        }//end for(_hTHistEff->GetNbinsX())
        cout << "Event "<<evt.id().event()<<" CaloClusterCog done..."<<endl;

}

}


using mu2e::CaloClusterCog;
DEFINE_ART_MODULE(CaloClusterCog);

