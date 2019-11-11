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

  CosmicCORSIKA::CosmicCORSIKA(const Config &conf)
      : _fluxConstant(conf.fluxConstant()),
        _tOffset(conf.tOffset()),
        _projectToTargetBox(conf.projectToTargetBox()),
        _showerAreaExtension(conf.showerAreaExtension()), // mm
        _targetBoxXmin(conf.targetBoxXmin()), // mm
        _targetBoxXmax(conf.targetBoxXmax()),  // mm
        _targetBoxYmin(conf.targetBoxYmin()), // mm
        _targetBoxYmax(conf.targetBoxYmax()),  // mm
        _targetBoxZmin(conf.targetBoxZmin()), // mm
        _targetBoxZmax(conf.targetBoxZmax())  // mm
  {
  }

  const unsigned int CosmicCORSIKA::getNumShowers() {
    return _primaries;
  }

  void CosmicCORSIKA::openFile(FILE *f) {
    in = f;
    fread(&_garbage, 4, 1, in);
  }

  const float CosmicCORSIKA::getLiveTime()
  {
    const float area = (_targetBoxXmax + 2 * _showerAreaExtension - _targetBoxXmin) * (_targetBoxZmax + 2 * _showerAreaExtension - _targetBoxZmin) * 1e-6; //m^2
    const float eslope = -2.7;
    const float lowE = 1.3;     // GeV
    const float highE = 1e6; // GeV
    const float EiToOneMinusGamma = pow(lowE, 1 + eslope);
    const float EfToOneMinusGamma = pow(highE, 1 + eslope);
    // http://pdg.lbl.gov/2018/reviews/rpp2018-rev-cosmic-rays.pdf eq. 29.2
    return _primaries / (M_PI * area * _fluxConstant * (EfToOneMinusGamma - EiToOneMinusGamma) / (1. + eslope));
  }

  CosmicCORSIKA::~CosmicCORSIKA(){
    std::cout << "Total number of primaries: " << getNumShowers() << std::endl;
    std::cout << "Simulated live-time: " << getLiveTime() << std::endl;
  }


  float CosmicCORSIKA::wrapvarBoxNo(const float var, const float low, const float high, int &boxno)
  {
    //wrap variable so that it's always between low and high
    boxno = int(floor(var / (high - low)));
    return (var - (high - low) * floor(var / (high - low))) + low;
  }

  float CosmicCORSIKA::wrapvar( const float var, const float low, const float high){
    //wrap variable so that it's always between low and high
    return (var - (high - low) * floor(var/(high-low))) + low;
  }


 bool CosmicCORSIKA::genEvent(std::map<std::pair<int,int>, GenParticleCollection> &particles_map) {
    float block[273]; // all data blocks have this size
    // static TDatabasePDG *pdgt = TDatabasePDG::Instance();

    bool running = true; // end run condition
    bool validParticleSubBlock = true;
    int particleDataSubBlockCounter = 0;

    if (feof(in))
      return false;

    while (running)
    {
      std::vector<GenParticle> showerParticles;
      fread(block, sizeof(block), 1, in); // blocks have 273 informations
      char blockName[5];
      memcpy(blockName, block, 4);
      blockName[4] = '\0';
      _loops++;

      if (_loops % 21 == 0)
      { // 2 _garbage data floats between every 21 blocks

        fread(&_garbage, 4, 1, in);
        fread(&_garbage, 4, 1, in);
      }

      // std::cout << blockName << " " << sizeof(block) << " " << (int)block[0] << std::endl;

      //__________RUN HEADER

      if (strcmp(blockName, "RUNH") == 0)
      {
        continue;
      }

      //__________EVENT HEADER
      else if (strcmp(blockName, "EVTH") == 0)
      {
        _primaries++;
        validParticleSubBlock = true;
      }
      //__________EVENT END

      else if (strcmp(blockName, "EVTE") == 0)
      {
        if (particleDataSubBlockCounter > 0)
        {
          running = false;
        }
        particleDataSubBlockCounter = 0;
      }
      //__________END OF RUN
      else if (strcmp(blockName, "RUNE") == 0)
      {
        _loops = 0;
        return false;
        // end run condition
      }
      else if (strcmp(blockName, "LONG")==0)
      {
        continue;
      }
      //__________PARTICLE DATA SUB-BLOCKS

      else
      {

        if (validParticleSubBlock == true)
        {
          for (int l = 1; l < 40; l++)
          { // each sub-block has up to 39 particles. If a sub-block is not fulfilled, trailing zeros are added

            int k = 7 * (l - 1);
            float id = block[k];
            //cout << "	" << k << " " << block[k] << " " << block[k+1] << " " << block[k+2] << " " << block[k+3] << " " << block[k+4] << " " << block[k+5] << " " << block[k+6] << endl;	// for verification purposes only

            if ((int)(id / 1000) != 0)
            {
              if (id < 75000)
              {
                // float dxyz[3] = {-block[k + 1], -block[k + 3], block[k + 2]};
                // float xyz[3] = {-block[k + 4], FDHalfHeight, block[k + 5]};
                // if (CheckTPCIntersection(xyz, dxyz)) {
                if (block[k + 1] != 0)
                {
                  id = (int)(id / 1000);
                  const int pdgId = corsikaToPdgId.at(id);
                  const float P_x = block[k + 2] * _GeV2MeV;
                  const float P_y = -block[k + 3] * _GeV2MeV;
                  const float P_z = block[k + 1] * _GeV2MeV;

                  int boxnox = 0, boxnoz = 0;

                  const float x = wrapvarBoxNo(block[k + 5] * _cm2mm, _targetBoxXmin - _showerAreaExtension, _targetBoxXmax + _showerAreaExtension, boxnox);
                  const float z = wrapvarBoxNo(-block[k + 4] * _cm2mm, _targetBoxZmin - _showerAreaExtension, _targetBoxZmax + _showerAreaExtension, boxnoz);
                  std::pair xz(boxnox, boxnoz);
                  const float m = pdt->particle(pdgId).ref().mass(); // to MeV

                  const float energy = safeSqrt(P_x * P_x + P_y * P_y + P_z * P_z + m * m);

                  const Hep3Vector position(x, 0, z);
                  const HepLorentzVector mom4(P_x, P_y, P_z, energy);

                  const float particleTime = block[k + 6];

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

                  particleDataSubBlockCounter++;
                }
              }
            }
            else
            {
              validParticleSubBlock = false;
            }
          }
        }
      }
    }

    return true;
  }


  bool CosmicCORSIKA::generate( GenParticleCollection& genParts)
  {

    // loop over particles in the truth object

    bool passed = false;
    while (!passed) {

      if (_particles_map.size() == 0)
      {
        if (!genEvent(_particles_map))
        {
          return false;
        }
      }

      GenParticleCollection particles = _particles_map.begin()->second;
      GenParticleCollection crossingParticles;

      float timeOffset = std::numeric_limits<float>::max();

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
          genParts.push_back(GenParticle(part.pdgId(), part.generatorId(), part.position(), part.momentum(), part.time()-timeOffset+_tOffset));
      }
      _particles_map.erase(_particles_map.begin()->first);

      if (genParts.size() != 0)
        passed = true;

    }

    return true;

  }

}// end namespace
