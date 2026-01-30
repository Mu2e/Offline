#include "Offline/BTrkLegacy/inc/ExternalInfo.hh"

#include <stdexcept>

FileFinderInterface   const* ExternalInfo::_findFile     = nullptr;
ParticleInfoInterface const* ExternalInfo::_particleInfo = nullptr;

FileFinderInterface const* ExternalInfo::fileFinderInstance(){
  if ( _findFile != nullptr ) return _findFile;
  throw std::runtime_error("BaBar ExternalInfo::fileFinderInstance: instance pointer is null");
}

ParticleInfoInterface const* ExternalInfo::particleInfoInstance(){
  if ( _particleInfo !=  nullptr ) return _particleInfo;
  throw std::runtime_error("BaBar ExternalInfo::particleInfoInstance: instance pointer is null");
}
