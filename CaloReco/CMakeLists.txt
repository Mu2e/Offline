cet_make_library(
    SOURCE
      src/CaloRawWFProcessor.cc
      src/CaloTemplateWFProcessor.cc
      src/CaloTemplateWFUtil.cc
    LIBRARIES PUBLIC
      
      offline::ConditionsService
      offline::Mu2eUtilities
)

cet_build_plugin(CaloHitMakerFast art::module
    REG_SOURCE src/CaloHitMakerFast_module.cc
    LIBRARIES REG
      offline::CaloReco
      
      offline::ConditionsService
      offline::DataProducts
      offline::GeneralUtilities
      offline::RecoDataProducts
)

cet_build_plugin(CaloHitMaker art::module
    REG_SOURCE src/CaloHitMaker_module.cc
    LIBRARIES REG
      offline::CaloReco
      
      offline::CalorimeterGeom
      offline::DataProducts
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(CaloRecoDigiMaker art::module
    REG_SOURCE src/CaloRecoDigiMaker_module.cc
    LIBRARIES REG
      offline::CaloReco
      
      offline::ConditionsService
      offline::DAQConditions
      offline::ProditionsService
      offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloReco)
