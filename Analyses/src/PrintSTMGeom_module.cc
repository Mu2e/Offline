//
// Print some information about the STM geometry
//
//
// Original author Rob Kutschke
// Modified by A. Palladino
//


#include "GeometryService/inc/GeomHandle.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "STMGeom/inc/STM.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace mu2e {

  class PrintSTMGeom : public art::EDAnalyzer {
  public:

    explicit PrintSTMGeom(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& e) override;

    void beginRun ( const art::Run& r) override;

  private:

  };

  PrintSTMGeom::PrintSTMGeom(fhicl::ParameterSet const& pset ):
    EDAnalyzer(pset){
  }

  void PrintSTMGeom::analyze(const art::Event& ){}

  void PrintSTMGeom::beginRun(const art::Run& run){
    // Get access to the master geometry system and its run time config.
    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const * _config = &(geom->config());

    int const verbosityLevel = _config->getInt("stm.verbosityLevel",0);
    if ( verbosityLevel > 0 ) {
       
       
       
      //std::cout << __func__ << " mstmMotherPositionInMu2e = " << mstmMotherPositionInMu2e << endl;       
      //std::cout << __func__ << " mstmReferencePositionInMu2e = " << mstmReferencePositionInMu2e << endl; 
      //It's critical for STM studies to know the radii and thicknesses of materials along the
      //beamline between the Stopping Target and the STM detector, print a summary:
      GeomHandle<StoppingTarget> target;
      std::cout << __func__ << " Stopping Target n_foils   = " << target->nFoils() << std::endl;       
      std::cout << __func__ << " Stopping Target z_min     = " << target->foil(0).centerInMu2e().z()-target->foil(0).halfThickness() << std::endl;       
      std::cout << __func__ << " Stopping Target z_max     = " << target->foil(target->nFoils()-1).centerInMu2e().z()+target->foil(target->nFoils()-1).halfThickness() << std::endl;
      
      std::vector<double> mbs_plug_rIn; _config->getVectorDouble("mbs.CLV2InnerRadii", mbs_plug_rIn);
      std::cout << __func__ << " MBS endplug inner radius     = " << mbs_plug_rIn[0] << endl;       
      std::cout << __func__ << " MBS endplug z_center         = " << _config->getDouble("mbs.MBSCZ") + _config->getDouble("mbs.CLV2ZrelCntr") << endl;       
      std::vector<double> mbs_plug_length; _config->getVectorDouble("mbs.CLV2Lengths", mbs_plug_length);
      double mbs_halfZ = mbs_plug_length[0]/2.0 + mbs_plug_length[1]/2.0;
      std::cout << __func__ << " MBS endplug z_halflength     = " << mbs_halfZ << endl;
      std::cout << __func__ << " MBS endplug z extent         = [" << _config->getDouble("mbs.MBSCZ") + _config->getDouble("mbs.CLV2ZrelCntr") - mbs_halfZ <<","
                                                                   << _config->getDouble("mbs.MBSCZ") + _config->getDouble("mbs.CLV2ZrelCntr") + mbs_halfZ << std::endl;
      
      GeomHandle<DetectorSolenoid> ds;
      std::cout << __func__ << " IFB window material          = " << _config->getString("ifb.endwindow.material") << endl;       
      std::cout << __func__ << " IFB window outer radius      = " << _config->getDouble("ifb.endwindow.rOut") << endl;       
      std::cout << __func__ << " IFB window z_center          = " << ds->cryoZMax()+_config->getDouble("ifb.endwindow.z") << endl;       
      std::cout << __func__ << " IFB window z_halflength      = " << _config->getDouble("ifb.endwindow.halfLength") << endl;
      std::cout << __func__ << " IFB window z extent          = [" << ds->cryoZMax()+_config->getDouble("ifb.endwindow.z")-_config->getDouble("ifb.endwindow.halfLength") <<","
                                                                   << ds->cryoZMax()+_config->getDouble("ifb.endwindow.z")+_config->getDouble("ifb.endwindow.halfLength")<<"]"<<std::endl;
      std::cout << __func__ << " IFB flange material          = " << _config->getString("ifb.material") << endl;       
      std::cout << __func__ << " IFB flange inner radius      = " << _config->getDouble("ifb.endplug.rIn") << endl;       
      std::cout << __func__ << " IFB flange z_center          = " << ds->cryoZMax()+_config->getDouble("ifb.endplug.z") << endl;       
      std::cout << __func__ << " IFB flange z_halflength      = " << _config->getDouble("ifb.endplug.halfLength") << endl;
      std::cout << __func__ << " IFB flange z extent          = [" << ds->cryoZMax()+_config->getDouble("ifb.endplug.z")-_config->getDouble("ifb.endplug.halfLength") <<","
                                                                   << ds->cryoZMax()+_config->getDouble("ifb.endplug.z")+_config->getDouble("ifb.endplug.halfLength")<<"]"<<std::endl;

      
      std::cout << __func__ << " Ext. Neutron Shielding hole rOut    = " << _config->getDouble("ExtShieldDownstream.holeRadiusType11Box4Hole1") << endl;
      std::vector<double> endcap_center; _config->getVectorDouble("ExtShieldDownstream.centerType11Box4",endcap_center);
      std::cout << __func__ << " Ext. Neutron (Endcap) Shielding z_center     = " << endcap_center[2] << endl;
      std::cout << __func__ << " Ext. Neutron (Endcap) Shielding z_fulllength = " << _config->getDouble("ExtShieldDownstream.holeLengthType11Box4Hole1") << endl;
      std::cout << __func__ << " Ext. Neutron (Endcap) Shielding z_halflength = " << _config->getDouble("ExtShieldDownstream.holeLengthType11Box4Hole1")/2.0 << endl;
      std::cout << __func__ << " Ext. Neutron (Endcap) Shielding z extent     = [" << endcap_center[2]-_config->getDouble("ExtShieldDownstream.holeLengthType11Box4Hole1")/2.0 <<","
                                                                                   << endcap_center[2]+_config->getDouble("ExtShieldDownstream.holeLengthType11Box4Hole1")/2.0<<"]"<<std::endl;


      GeomHandle<CosmicRayShield> CRS;
      const CRSScintillatorShield &sectionCRVD = CRS->getCRSScintillatorShield(13); //for CRV-D or 18, 19, 20 for CRV-D2,3,4
      std::cout<<"CRV section: "<<sectionCRVD.getName()<<":"<<std::endl;
      double y_min_CRVD = 10000.0;
      double z_min_CRVD = 1000000.0; 
      double z_max_CRVD = -1000000.0;
      size_t nModules = sectionCRVD.nModules();
      for(size_t m=0; m<nModules; m++)
	{
	  std::cout<<"  Module "<<m<<":"<<std::endl;
	  const CRSScintillatorModule &module = sectionCRVD.getModule(m);
	  size_t nLayers = module.nLayers();
	  for(size_t l=0; l<nLayers; l++)
	    {
	      const CRSScintillatorLayer &layer = module.getLayer(l);
	      const CLHEP::Hep3Vector &position = layer.getPosition();
	      double halfThickness = layer.getHalfThickness();
	      double halfWidth = layer.getHalfWidth();
	      double halfLength = layer.getHalfLength();
              if (position[1]-halfWidth < y_min_CRVD){
                y_min_CRVD = position[1]-halfWidth;
	      }
              if (position[2]-halfThickness < z_min_CRVD){
                z_min_CRVD = position[2]-halfThickness;
              }
              if (position[2]+halfThickness > z_max_CRVD){
                z_max_CRVD = position[2]+halfThickness;
              }
	      std::cout<<"    Layer "<<l<<" position: "<<position<<"   half thickness (z): "<<halfThickness<<"   half width (y): "<<halfWidth<<"   half length (x): "<<halfLength<<std::endl;
	    }
	}
      const CRSScintillatorShield &sectionCRVD2 = CRS->getCRSScintillatorShield(18); //for CRV-D or 18, 19, 20 for CRV-D2,3,4                                                                                       
      std::cout<<"CRV section: "<<sectionCRVD2.getName()<<":"<<std::endl;
      double x_max_CRVD2 = -10000.0;
      size_t nModulesCRVD2 = sectionCRVD2.nModules();
      for(size_t m=0; m<nModulesCRVD2; m++)
        {
	  std::cout<<"  Module "<<m<<":"<<std::endl;
          const CRSScintillatorModule &module = sectionCRVD2.getModule(m);
          size_t nLayers = module.nLayers();
          for(size_t l=0; l<nLayers; l++)
            {
              const CRSScintillatorLayer &layer = module.getLayer(l);
              const CLHEP::Hep3Vector &position = layer.getPosition();
              double halfThickness = layer.getHalfThickness();
              double halfWidth = layer.getHalfWidth();
              double halfLength = layer.getHalfLength();
              if (position[0]+halfLength > x_max_CRVD2){
                x_max_CRVD2 = position[0]+halfLength;
              }
	      std::cout<<"    Layer "<<l<<" position: "<<position<<"   half thickness (z): "<<halfThickness<<"   half width (y): "<<halfWidth<<"   half length (x): "<<halfLength<<std::endl;
            }
        }

      const CRSScintillatorShield &sectionCRVD3 = CRS->getCRSScintillatorShield(19); //for CRV-D or 18, 19, 20 for CRV-D2,3,4                                                                                       
      std::cout<<"CRV section: "<<sectionCRVD3.getName()<<":"<<std::endl;
      double x_min_CRVD3 = 10000.0;
      size_t nModulesCRVD3 = sectionCRVD3.nModules();
      for(size_t m=0; m<nModulesCRVD3; m++)
        {
	  std::cout<<"  Module "<<m<<":"<<std::endl;
          const CRSScintillatorModule &module = sectionCRVD3.getModule(m);
          size_t nLayers = module.nLayers();
          for(size_t l=0; l<nLayers; l++)
            {
              const CRSScintillatorLayer &layer = module.getLayer(l);
              const CLHEP::Hep3Vector &position = layer.getPosition();
              double halfThickness = layer.getHalfThickness();
              double halfWidth = layer.getHalfWidth();
              double halfLength = layer.getHalfLength();
              if (position[0]-halfLength < x_min_CRVD3){
                x_min_CRVD3 = position[0]-halfLength;
              }
	      std::cout<<"    Layer "<<l<<" position: "<<position<<"   half thickness (z): "<<halfThickness<<"   half width (y): "<<halfWidth<<"   half length (x): "<<halfLength<<std::endl;
            }
        }

      const CRSScintillatorShield &sectionCRVD4 = CRS->getCRSScintillatorShield(20); //for CRV-D or 18, 19, 20 for CRV-D2,3,4                                                                                       
      std::cout<<"CRV section: "<<sectionCRVD4.getName()<<":"<<std::endl;
      double y_max_CRVD4 = -10000.0;
      size_t nModulesCRVD4 = sectionCRVD4.nModules();
      for(size_t m=0; m<nModulesCRVD4; m++)
        {
	  std::cout<<"  Module "<<m<<":"<<std::endl;
          const CRSScintillatorModule &module = sectionCRVD4.getModule(m);
          size_t nLayers = module.nLayers();
          for(size_t l=0; l<nLayers; l++)
            {
              const CRSScintillatorLayer &layer = module.getLayer(l);
              const CLHEP::Hep3Vector &position = layer.getPosition();
              double halfThickness = layer.getHalfThickness();
              double halfWidth = layer.getHalfWidth();
              double halfLength = layer.getHalfLength();
              if (position[1]+halfWidth > y_max_CRVD4){
                y_max_CRVD4 = position[1]+halfWidth;
              }
	      std::cout<<"    Layer "<<l<<" position: "<<position<<"   half thickness (z): "<<halfThickness<<"   half width (y): "<<halfWidth<<"   half length (x): "<<halfLength<<std::endl;
            }
        }
      std::cout << __func__ << " CRV z_min = "<<z_min_CRVD<<" mm"<<std::endl;
      std::cout << __func__ << " CRV z_max = "<<z_max_CRVD<<" mm"<<std::endl;
      float x_center = _config->getDouble("ExtShieldDownstream.detecHoleX");
      float y_center = 0.0;//_config->getDouble("ExtShieldDownstream.detecHoleY");
      std::cout << __func__ << " CRV axial opening +y = "<<y_min_CRVD -y_center<<" mm"<<std::endl;
      std::cout << __func__ << " CRV axial opening +x = "<<x_min_CRVD3-x_center<<" mm"<<std::endl;
      std::cout << __func__ << " CRV axial opening -x = "<<x_center-x_max_CRVD2<<" mm"<<std::endl;
      std::cout << __func__ << " CRV axial opening -y = "<<y_center-y_max_CRVD4<<" mm"<<std::endl;
      //std::cout << __func__ << " STM UpStr Concrete Shielding center       = " << mstmUpStreamWallPositionInMu2e << endl;
      //std::cout << __func__ << " Ext. Neutron Shielding hole rOut    = " << _config->getDouble("ExtShieldDownstream.holeRadiusType11Box4Hole1")<<endl;
      //double ExtShieldDownstream.holeLengthType11Box4Hole1 =  915.0;

      GeomHandle<STM> STM;
      
      const double crvShieldPipeZcenter  =   STM->getSTMMagnetPtr()->originInMu2e().z() 
                                          - STM->getSTMMagnetPtr()->zHalfLength() 
                                          - STM->getSTMShieldPipePtr()->dnStrSpace()
                                          - 2.0*STM->getSTMShieldPipePtr()->dnStrWallHalflength()
                                          - STM->getSTMShieldPipePtr()->pipeHalfLength();
      const double crvShieldPipeZmin  =  crvShieldPipeZcenter - STM->getSTMShieldPipePtr()->pipeHalfLength();
      const double crvShieldPipeZmax  =  crvShieldPipeZcenter + STM->getSTMShieldPipePtr()->pipeHalfLength();
      std::cout << __func__ << " STM CRVshield cylinder z center = "<<crvShieldPipeZcenter<<std::endl;
      std::cout << __func__ << " STM CRVshield cylinder z extent = ["<<crvShieldPipeZmin<<","<<crvShieldPipeZmax<<"]"<<std::endl;
      std::cout << __func__ << " STM CRVshield cylinder r inner  = "<<STM->getSTMShieldPipePtr()->radiusIn()<<std::endl;
      std::cout << __func__ << " STM CRVshield cylinder r outer  = "<<STM->getSTMShieldPipePtr()->radiusOut()<<std::endl;
      
      const double crvShieldWallZcenter  =   STM->getSTMMagnetPtr()->originInMu2e().z() 
                                          - STM->getSTMMagnetPtr()->zHalfLength() 
                                          - STM->getSTMShieldPipePtr()->dnStrSpace()
                                          - STM->getSTMShieldPipePtr()->dnStrWallHalflength();
      const double crvShieldWallZmin  =  crvShieldWallZcenter - STM->getSTMShieldPipePtr()->dnStrWallHalflength();
      const double crvShieldWallZmax  =  crvShieldWallZcenter + STM->getSTMShieldPipePtr()->dnStrWallHalflength();
      std::cout << __func__ << " STM CRVshield wall z center     = "<<crvShieldWallZcenter<<std::endl;
      std::cout << __func__ << " STM CRVshield wall z extent     = ["<<crvShieldWallZmin<<","<<crvShieldWallZmax<<"]"<<std::endl;
      
      std::cout << __func__ << " STM magnet center               = "<<STM->getSTMMagnetPtr()->originInMu2e()<< std::endl;
      std::cout << __func__ << " STM magnet z halflength         = "<<STM->getSTMMagnetPtr()->zHalfLength()<< std::endl;
      std::cout << __func__ << " STM magnet z extent             = ["<<STM->getSTMMagnetPtr()->originInMu2e().z()-STM->getSTMMagnetPtr()->zHalfLength()<<","<<STM->getSTMMagnetPtr()->originInMu2e().z()+STM->getSTMMagnetPtr()->zHalfLength()<<"]"<<std::endl;
      std::cout << __func__ << " STM magnet x opening halflength = "<<STM->getSTMMagnetPtr()->xHoleHalfLength()<<std::endl;
      std::cout << __func__ << " STM magnet y opening halflength = "<<STM->getSTMMagnetPtr()->yHoleHalfLength()<<std::endl;
      
     
      std::cout<<__func__<<" STM FOV Coll (lead) z_center              = "<< STM->getSTMFOVCollimatorPtr()->originInMu2e().z() <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (lead) z_halflength          = "<< STM->getSTMFOVCollimatorPtr()->halfLength()<<std::endl;
      std::cout<<__func__<<" STM FOV Coll (lead) z_min                 = "<< STM->getSTMFOVCollimatorPtr()->originInMu2e().z()-STM->getSTMFOVCollimatorPtr()->halfLength()<<std::endl;
      std::cout<<__func__<<" STM FOV Coll (lead) z_max                 = "<< STM->getSTMFOVCollimatorPtr()->originInMu2e().z()+STM->getSTMFOVCollimatorPtr()->halfLength()<<std::endl;
      std::cout<<__func__<<" STM FOV Coll (poly) z_center              = "<< STM->getSTMFOVCollimatorPtr()->originInMu2e().z() <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (poly) z_halflength          = "<< STM->getSTMFOVCollimatorPtr()->linerHalfLength() <<std::endl;    
      std::cout<<__func__<<" STM FOV Coll (poly) z_min                 = "<< STM->getSTMFOVCollimatorPtr()->originInMu2e().z()-STM->getSTMFOVCollimatorPtr()->linerHalfLength() <<std::endl;    
      std::cout<<__func__<<" STM FOV Coll (poly) z_max                 = "<< STM->getSTMFOVCollimatorPtr()->originInMu2e().z()+STM->getSTMFOVCollimatorPtr()->linerHalfLength() <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (poly) r                     = "<< STM->getSTMFOVCollimatorPtr()->hole1RadiusUpStr() <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (poly) absorber z_halflength = "<< _config->getDouble("stm.FOVcollimator.absorber.halfLength") <<std::endl;
       
      std::cout<<__func__<<" STM SS Coll (lead)     z_center     = "<< STM->getSTMSSCollimatorPtr()->originInMu2e().z() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (lead)     z_halflength = "<< STM->getSTMSSCollimatorPtr()->halfLength() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (lead)     z_min        = "<< STM->getSTMSSCollimatorPtr()->originInMu2e().z()-STM->getSTMSSCollimatorPtr()->halfLength() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (lead)     z_max        = "<< STM->getSTMSSCollimatorPtr()->originInMu2e().z()+STM->getSTMSSCollimatorPtr()->halfLength() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (liner material)        = "<< STM->getSTMSSCollimatorPtr()->linerMaterial()<<std::endl;    
      std::cout<<__func__<<" STM SS Coll (tungsten) r_DnStr left = "<< STM->getSTMSSCollimatorPtr()->hole1RadiusDnStr() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (tungsten) r_DnStr right= "<< STM->getSTMSSCollimatorPtr()->hole2RadiusDnStr() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (tungsten) x_halfwidth  = "<< STM->getSTMSSCollimatorPtr()->linerHalfWidth() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (tungsten) y_halfheight = "<< STM->getSTMSSCollimatorPtr()->linerHalfHeight()<<std::endl;    
      
      
      //std::cout<<__func__<<" STM SS Coll (tungsten) z_min        = "<< stmSSCollPositionInMu2e2.z()-stmSSCollHalfLength2 <<std::endl;  
      
      
//       std::cout << __func__ << " STM UpStr Concrete Shielding z_halflength = " << _config->getDouble("mstm.wallUpStr.halfLength") << endl;
//       std::cout << __func__ << " STM UpStr Concrete Shielding hole radius  = " << _config->getDouble("mstm.wallUpStr.holeRadius") << endl;
//       //std::cout << __func__ << " STM Magnet center            = " << mstmMagnetPositionInMu2e << endl;
//       std::cout << __func__ << " STM Magnet z_halflength      = " << _config->getDouble("mstm.magnet.halfLength") << endl;
//       std::cout << __func__ << " STM Magnet hole halfHeight   = " << _config->getDouble("mstm.magnet.holeHalfHeight") << endl;
//       //std::cout << __func__ << " STM Collimator1 center       = " << mstmColl1PositionInMu2e << endl;
//       std::cout << __func__ << " STM Collimator1 z_halflength = " << _config->getDouble("mstm.collimator1.halfLength") << endl;
//       std::cout << __func__ << " STM Collimator1 hole radius  = " << _config->getDouble("mstm.collimator1.holeRadius") << endl;
//       //std::cout << __func__ << " STM Collimator2 center       = " << mstmColl2PositionInMu2e << endl;
//       std::cout << __func__ << " STM Collimator2 z_halflength = " << _config->getDouble("mstm.collimator2.halfLength") << endl;
//       std::cout << __func__ << " STM Collimator2 hole radius  = " << _config->getDouble("mstm.collimator2.holeRadius") << endl;
//       //std::cout << __func__ << " STM absorber material        = " << _config->getString("mstm.absorber.material") << endl;       
//       //std::cout << __func__ << " STM absorber hole radius     = " << _config->getDouble("mstm.absorber.rOut") << endl;
//       //std::cout << __func__ << " STM absorber center          = " << mstmAbsorberPositionInMu2e << endl;
//       //std::cout << __func__ << " STM absorber z_halflength    = " << _config->getDouble("mstm.absorber.halfLength") << endl;
//       //std::cout << __func__ << " STM Collimator3 center       = " << mstmColl3PositionInMu2e << endl;
//       std::cout << __func__ << " STM Collimator3 z_halflength = " << _config->getDouble("mstm.collimator3.halfLength") << endl;
//       //std::cout << __func__ << " STM Collimator3 hole rad UpStr=" << _config->getDouble("mstm.collimator3.holeRadiusUpStr") << endl;
//       //std::cout << __func__ << " STM Collimator3 hole rad DnStr=" << _config->getDouble("mstm.collimator3.holeRadiusDnStr") << endl;
//       //std::cout << __func__ << " STM Detector center          = " << mstmCrystalPositionInMu2e << endl;
//       std::cout << __func__ << " STM Detector z_halflength    = " << _config->getDouble("mstm.crystal.halfLength") << endl;
//       std::cout << __func__ << " STM Detector radius          = " << _config->getDouble("mstm.crystal.rOut") << endl;

    }

    // Tracker const& tracker(*GeomHandle<Tracker>());
    // cout << "Tracker: " << tracker.nDevices() << endl;
    // for ( auto const& dev : tracker.getDevices() ){
    //   for ( auto const& sec : dev.getSectors() ){
    //     StrawId sid( sec.id(), 0, 0 );
    //     Straw const& straw = sec.getStraw(sid);
    //     double phi  = straw.direction().phi();
    //     double z    = straw.getMidPoint().z() - dev.origin().z();
    //     double phi1 = phi/M_PI*180.;
    //     cout << "sector: "
    //          << sec.id()      << " "
    //          << sid           << " "
    //          << straw.index() << " : "
    //          << z             << " "
    //          << phi1
    //          << endl;
    //   }
    // }

  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::PrintSTMGeom);
