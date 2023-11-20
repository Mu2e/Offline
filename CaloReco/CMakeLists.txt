cet_make_library(
    SOURCE
      src/CaloRawWFProcessor.cc
      src/CaloTemplateWFProcessor.cc
      src/CaloTemplateWFUtil.cc
    LIBRARIES PUBLIC
      
      Offline::ConditionsService
      Offline::Mu2eUtilities
)

cet_build_plugin(CaloHitMakerFast art::module
    REG_SOURCE src/CaloHitMakerFast_module.cc
    LIBRARIES REG
      Offline::CaloReco
      
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(CaloHitMaker art::module
    REG_SOURCE src/CaloHitMaker_module.cc
    LIBRARIES REG
      Offline::CaloReco
      
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloRecoDigiMaker art::module
    REG_SOURCE src/CaloRecoDigiMaker_module.cc
    LIBRARIES REG
      Offline::CaloReco
      
      Offline::ConditionsService
      Offline::DAQConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloReco/fcl)
