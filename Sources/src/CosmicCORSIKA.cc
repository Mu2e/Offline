////////////////////////////////////////////////////////////////////////
/// \file  CosmicCORSIKA_module.cc
/// \brief Generator for cosmic-ray secondaries based on pre-generated CORSIKA binary outputs.
///
/// \author Stefano Roberto Soleti roberto@lbl.gov
////////////////////////////////////////////////////////////////////////

#include "Sources/inc/CosmicCORSIKA.hh"

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

namespace mu2e {

  CosmicCORSIKA::CosmicCORSIKA(const Config &conf, SeedService::seed_t seed)
      : _fluxConstant(conf.fluxConstant()),
        _tOffset(conf.tOffset()),
        _projectToTargetBox(conf.projectToTargetBox()),
        _showerAreaExtension(conf.showerAreaExtension()), // mm
        _targetBoxXmin(conf.targetBoxXmin()), // mm
        _targetBoxXmax(conf.targetBoxXmax()),  // mm
        _targetBoxYmin(conf.targetBoxYmin()), // mm
        _targetBoxYmax(conf.targetBoxYmax()),  // mm
        _targetBoxZmin(conf.targetBoxZmin()), // mm
        _targetBoxZmax(conf.targetBoxZmax()),  // mm
        _engine(seed),
        _randFlatX(_engine, -(_targetBoxXmax-_targetBoxXmin+_showerAreaExtension)/2, +(_targetBoxXmax-_targetBoxXmin+_showerAreaExtension)/2),
        _randFlatZ(_engine, -(_targetBoxZmax-_targetBoxZmin+_showerAreaExtension)/2, +(_targetBoxZmax-_targetBoxZmin+_showerAreaExtension)/2)
  {
  }

  void CosmicCORSIKA::openFile(ifstream *f, unsigned &runNumber, float &lowE, float &highE)
  {
    input = f;
    _current_event_number = -1;
    _event_count = 0;
    _run_number = -1;
    _infmt = Format::UNDEFINED;

    while( input->read(_buf.ch, 4)) {
      unsigned reclen = _buf.in[0];
      // CORSIKA records are in units of 4 bytes
      if(reclen % 4) {
        throw std::runtime_error("Error: record size not a multiple of 4");
      }

      // We will be looking at at least 8 bytes to determine the
      // input file format, and all real CORSIKA records are longer
      // than that.
      if(reclen < 2*4) {
        throw std::runtime_error("Error: reclen too small");
      }

      // Read the full record
      if(!input->read(_buf.ch, reclen)) {
        break;
      }

      // Determine the format and and store the decision for future blocks.
      // We are starting file read, so should see the RUNH marker
      // In COMPACT format each block is preceded by 4 bytes
      // giving the size of the block in words.

      if(!strncmp(_buf.ch+0, "RUNH", 4)) {
        std::cout<<"Reading NORMAL format"<<std::endl;
        _infmt = Format::NORMAL;
      }
      else if(!strncmp(_buf.ch+4, "RUNH", 4)) {
        std::cout<<"Reading COMPACT format"<<std::endl;
        _infmt = Format::COMPACT;
      }
      else {
        throw std::runtime_error("Error: did not find the RUNH record to determine COMPACT flag");
      }

      unsigned iword = 0;
      if(_infmt == Format::COMPACT) {
        // Move to the beginning of the actual block
        ++iword;
      }

      if(!strncmp(_buf.ch+4*iword, "RUNH", 4)) {
        runNumber = lrint(_buf.fl[1+iword]);
        lowE = float(_buf.fl[16+iword]);
        highE = float(_buf.fl[17+iword]);
      }

      break;

    }
    input->clear();
    input->seekg(0, ios::beg);
  }

  CosmicCORSIKA::~CosmicCORSIKA(){
  }


  float CosmicCORSIKA::wrapvarBoxNo(const float var, const float low, const float high, int &boxno)
  {
    //wrap variable so that it's always between low and high
    boxno = int(floor(var / (high - low)));
    return (var - (high - low) * floor(var / (high - low))) + low;
  }

  bool CosmicCORSIKA::genEvent(std::map<std::pair<int,int>, GenParticleCollection> &particles_map) {

      const float xOffset = _randFlatX.fire();
      const float zOffset = _randFlatZ.fire();

      // FORTRAN sequential records are prefixed with their length
      // in a 4-byte word
      while( input->read(_buf.ch, 4)) {

        unsigned reclen = _buf.in[0];
        // CORSIKA records are in units of 4 bytes
        if(reclen % 4) {
          throw std::runtime_error("Error: record size not a multiple of 4");
        }

        // We will be looking at at least 8 bytes to determine the
        // input file format, and all real CORSIKA records are longer
        // than that.
        if(reclen < 2*4) {
          throw std::runtime_error("Error: reclen too small");
        }

        if(reclen > 4*_fbsize_words) {
          throw std::runtime_error("Error: reclen too big");
        }

        // Read the full record
        if(!input->read(_buf.ch, reclen)) {
          break;
        }

        unsigned n_part = 0;
        //================================================================
        // Go over blocks in the record
        for(unsigned iword = 0; iword < reclen/4; ) {

          unsigned block_words = (_infmt == Format::COMPACT) ?
            _buf.in[iword] : 273;

          if(!block_words) {
            throw std::runtime_error("Got block_words = 0\n");
          }

          if(_infmt == Format::COMPACT) {
            // Move to the beginning of the actual block
            ++iword;
          }

          std::string event_marker =
            (_infmt == Format::NORMAL || !_event_count) ? "EVTH" : "EVHW";

          // Determine the type of the data block
          if(!strncmp(_buf.ch+4*iword, "RUNH", 4)) {
            _run_number = lrint(_buf.fl[1+iword]);
          }
          else if(!strncmp(_buf.ch+4*iword, "RUNE", 4)) {
            unsigned end_run_number = lrint(_buf.fl[1+iword]);
            unsigned end_event_count = lrint(_buf.fl[2+iword]);
            if(end_run_number != _run_number) {
              throw std::runtime_error("Error: run number mismatch in end of run record\n");
            }
            if(_event_count != end_event_count) {
              std::cerr<<"RUNE: _event_count = "<<_event_count<<" end record = "<<end_event_count<<std::endl;
              throw std::runtime_error("Error: event count mismatch in end of run record\n");
            }
            // Exit the read loop at the end of run
            _primaries = 0;
            return false;
          }
          else if(!strncmp(_buf.ch+4*iword, event_marker.data(), 4)) {
            ++_event_count;
            _current_event_number = lrint(_buf.fl[1+iword]);
            ++_primaries;
          }
          else if(!strncmp(_buf.ch+4*iword, "EVTE", 4)) {
            unsigned end_event_number = lrint(_buf.fl[1+iword]);
            if(end_event_number != _current_event_number) {
              throw std::runtime_error("Error: event number mismatch in end of event record\n");
            }
          }
          else {
            for (unsigned i_part = 0; i_part < block_words; i_part+=7) {
              unsigned id = _buf.fl[iword + i_part] / 1000;
              if (id == 0)
                continue;
              n_part++;
              const int pdgId = corsikaToPdgId.at(id);
              const float P_x = _buf.fl[iword + i_part + 2] * _GeV2MeV;
              const float P_y = -_buf.fl[iword + i_part + 3] * _GeV2MeV;
              const float P_z = _buf.fl[iword + i_part + 1] * _GeV2MeV;

              int boxnox = 0, boxnoz = 0;

              const float x = wrapvarBoxNo(_buf.fl[iword + i_part + 5] * _cm2mm + xOffset, _targetBoxXmin - _showerAreaExtension, _targetBoxXmax + _showerAreaExtension, boxnox);
              const float z = wrapvarBoxNo(-_buf.fl[iword + i_part + 4] * _cm2mm + zOffset, _targetBoxZmin - _showerAreaExtension, _targetBoxZmax + _showerAreaExtension, boxnoz);
              std::pair xz(boxnox, boxnoz);
              const float m = pdt->particle(pdgId).ref().mass(); // to MeV

              const float energy = safeSqrt(P_x * P_x + P_y * P_y + P_z * P_z + m * m);

              const Hep3Vector position(x, _targetBoxYmax, z);
              const HepLorentzVector mom4(P_x, P_y, P_z, energy);

              const float particleTime = _buf.fl[iword + i_part + 6] * _ns2s;
              GenParticle part(static_cast<PDGCode::type>(pdgId),
                    GenId::cosmicCORSIKA, position, mom4,
                    particleTime);

              if (particles_map.count(xz) == 0) {
                GenParticleCollection parts;
                parts.push_back(part);
                particles_map.insert({xz, parts});
              } else {
                particles_map[xz].push_back(part);
              }
            }

          }

          // Move to the next block
          iword += block_words;

        } // loop over blocks in a record

        // Here we expect the FORTRAN end of record padding,
        // read and verify its value.
        if(!input->read(_buf.ch, 4)) {
          break;
        }
        if(_buf.in[0] != reclen) {
          throw std::runtime_error("Error: unexpected FORTRAN record end padding");
        }

        if (n_part > 0) {
          return true;
        }

      } // loop over records
      return true;
  }

  bool CosmicCORSIKA::generate( GenParticleCollection& genParts, unsigned int &primaries)
  {
    // loop over particles in the truth object
    bool passed = false;
    while (!passed) {
      if (_particles_map.size() == 0)
      {
        if (!genEvent(_particles_map)) {
          return false;
        }
      }

      GenParticleCollection particles = _particles_map.begin()->second;
      GenParticleCollection crossingParticles;

      float timeOffset = std::numeric_limits<float>::max();
      primaries = _primaries;

      for (unsigned int i = 0; i < particles.size(); i++) {
        GenParticle particle = particles[i];

        _targetBoxIntersections.clear();
        VectorVolume particleTarget(particle.position(), particle.momentum().vect(),
                                    _targetBoxXmin, _targetBoxXmax,
                                    _targetBoxYmin, _targetBoxYmax,
                                    _targetBoxZmin, _targetBoxZmax);

        particleTarget.calIntersections(_targetBoxIntersections);

        if (_targetBoxIntersections.size() > 0 || !_projectToTargetBox)
          crossingParticles.push_back(particle);

        if (particle.time() < timeOffset)
          timeOffset = particle.time();
      }

      for (unsigned int i = 0; i < crossingParticles.size(); i++) {
          GenParticle part = crossingParticles[i];
          genParts.push_back(GenParticle(part.pdgId(), part.generatorId(), part.position(), part.momentum(), part.time()+_tOffset-timeOffset));
      }
      _particles_map.erase(_particles_map.begin()->first);

      if (genParts.size() != 0) {
        passed = true;
        _primaries = 0;
      }

    }

    return true;

  }

}// end namespace
