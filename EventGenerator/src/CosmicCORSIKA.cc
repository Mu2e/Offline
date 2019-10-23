////////////////////////////////////////////////////////////////////////
/// \file  CosmicCORSIKA_module.cc
/// \brief Generator for cosmic-ray secondaries based on pre-generated CORSIKA binary outputs.
///
/// \author Stefano Roberto Soleti roberto@lbl.gov
////////////////////////////////////////////////////////////////////////

#include "EventGenerator/inc/CosmicCORSIKA.hh"

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

namespace mu2e {

CosmicCORSIKA::CosmicCORSIKA(art::Run &run,
                             const Config &conf,
                             CLHEP::HepRandomEngine &engine)
    : _showerInputFiles(conf.showerInputFiles()),
      _fluxConstant(conf.fluxConstant()),
      _refY0(conf.refY0()), // mm
      _projectToTargetBox(conf.projectToTargetBox()),
      _tOffset(conf.tOffset()),
      _buffer(conf.buffer()), // mm
      _showerAreaExtension(conf.showerAreaExtension()), // mm
      _refPointChoice(conf.refPointChoice()),
      _targetBoxXmin(conf.targetBoxXmin()), // mm
      _targetBoxXmax(conf.targetBoxXmax()),  // mm
      _targetBoxYmin(conf.targetBoxYmin()), // mm
      _targetBoxYmax(conf.targetBoxYmax()),  // mm
      _targetBoxZmin(conf.targetBoxZmin()), // mm
      _targetBoxZmax(conf.targetBoxZmax()),  // mm
      _randomXZshift(conf.randomXZshift()), // mm
      _flat(engine, -_randomXZshift, _randomXZshift)
{

  _showerInputs = _showerInputFiles.size();

  if (_refY0 == 0.)
    mf::LogInfo("CosmicCORSIKA") << "Using 0. for _refY0!";

  _fileIndex = 0;
  mf::LogInfo("CosmicCORSIKA") << "Opening: " << _showerInputFiles[_fileIndex] << std::endl;
  in = fopen(_showerInputFiles[_fileIndex].c_str(), "r");
  fread(&_garbage, 4, 1, in); // skipping the first (and useless) information

  }

  const unsigned int CosmicCORSIKA::getNumShowers() {
    return _primaries;
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
    fclose(in);
  }

  bool CosmicCORSIKA::pointInBox(float x, float y, float x0, float y0,
                              float x1, float y1)
  {
    bool ret = false;
    if ((x >= x0) && (x <= x1) && (y >= y0) && (y <= y1))
    {
      ret = true;
    }
    return ret;
  }

  float CosmicCORSIKA::distance(const CLHEP::Hep3Vector &u, const CLHEP::Hep3Vector &v)
  {
    return safeSqrt((u.x() - v.x()) * (u.x() - v.x()) +
                    (u.y() - v.y()) * (u.y() - v.y()) +
                    (u.z() - v.z()) * (u.z() - v.z()));
  }

  void CosmicCORSIKA::calIntersections(Hep3Vector orig, Hep3Vector dir,
                                  std::vector<CLHEP::Hep3Vector> &intersections, float xMin, float xMax,
                                  float yMin, float yMax, float zMin, float zMax)
  {
    // roof: _targetBoxYmax, _targetBoxXmin, _targetBoxXmax, _targetBoxZmin, _targetBoxZmax
    // skip projection if the particle goes parallely to the plane
    if (dir.y() != 0.)
    {
      const float t = (yMax - orig.y()) / dir.y();
      const float x1 = dir.x() * t + orig.x();
      const float z1 = dir.z() * t + orig.z();
      // std::cout << "x1 " << x1 << " " << z1 << std::endl;

      if (pointInBox(x1, z1, xMin, zMin, xMax, zMax))
      {
        intersections.push_back(Hep3Vector(x1, yMax, z1));
      }
    }

    //east: zMin, xMin, xMax, yMin, yMax
    if (dir.z() != 0.)
    {
      const float t = (zMin - orig.z()) / dir.z();
      const float x1 = dir.x() * t + orig.x();
      const float y1 = dir.y() * t + orig.y();
      if (pointInBox(x1, y1, xMin, yMin, xMax, yMax))
      {
        intersections.push_back(Hep3Vector(x1, y1, zMin));
      }
    }

    //west: zMax, xMin, xMax, yMin, yMax
    if (dir.z() != 0.)
    {
      const float t = (zMax - orig.z()) / dir.z();
      const float x1 = dir.x() * t + orig.x();
      const float y1 = dir.y() * t + orig.y();
      if (pointInBox(x1, y1, xMin, yMin, xMax, yMax))
      {
        intersections.push_back(Hep3Vector(x1, y1, zMax));
      }
    }

    //south: xMin, yMin, yMax, zMin, zMax
    if (dir.x() != 0.)
    {
      const float t = (xMin - orig.x()) / dir.x();
      const float z1 = dir.z() * t + orig.z();
      const float y1 = dir.y() * t + orig.y();
      if (pointInBox(z1, y1, zMin, yMin, zMax, yMax))
      {
        intersections.push_back(Hep3Vector(xMin, y1, z1));
      }
    }

    //north: xMax, yMin, yMax, zMin, zMax
    if (dir.x() != 0.)
    {
      const float t = (xMax - orig.x()) / dir.x();
      const float z1 = dir.z() * t + orig.z();
      const float y1 = dir.y() * t + orig.y();
      if (pointInBox(z1, y1, zMin, yMin, zMax, yMax))
      {
        intersections.push_back(Hep3Vector(xMax, y1, z1));
      }
    }
  }

  float CosmicCORSIKA::wrapvar( const float var, const float low, const float high){
    //wrap variable so that it's always between low and high
    return (var - (high - low) * floor(var/(high-low))) + low;
  }

  void CosmicCORSIKA::genEvent(GenParticleCollection &genParts) {
    float block[273]; // all data blocks have this size
    // static TDatabasePDG *pdgt = TDatabasePDG::Instance();

    bool running = true; // end run condition
    bool validParticleSubBlock = true;
    int particleDataSubBlockCounter = 0;

    while (running && !feof(in))
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
        fclose(in);

        _fileIndex++;
        if (_fileIndex >= _showerInputs) {
          _fileIndex = 0;
          mf::LogWarning("CosmicCORSIKA") << "RE-USING FIRST CORSIKA FILE" << std::endl;
        }

        mf::LogInfo("CosmicCORSIKA") << "Opening: " << _showerInputFiles[_fileIndex] << std::endl;
        in = fopen(_showerInputFiles[_fileIndex].c_str(), "r");
        fread(&_garbage, 4, 1, in);
        _loops = 0;
        running = false; // end run condition
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
                  const float x = block[k + 5];
                  const float z = -block[k + 4];

                  GlobalConstantsHandle<ParticleDataTable> pdt;
                  ParticleDataTable const &pdt_ = *pdt;
                  const float m = pdt_.particle(pdgId).ref().mass(); // to MeV

                  const float energy = safeSqrt(P_x * P_x + P_y * P_y + P_z * P_z + m * m);

                  const Hep3Vector position(x, 0, z);
                  const HepLorentzVector mom4(P_x, P_y, P_z, energy);

                  const float particleTime = block[k + 6];

                  genParts.push_back(GenParticle(static_cast<PDGCode::type>(pdgId),
                                                  GenId::cosmicCORSIKA, position, mom4,
                                                  particleTime));
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
  }


  void CosmicCORSIKA::generate( GenParticleCollection& genParts )
  {

    // loop over particles in the truth object

    GeomHandle<Mu2eEnvelope> env;
    GeomHandle<WorldG4> worldGeom;
    GeomHandle<DetectorSystem> detsys;

    const float deltaX = 1; // mm
    _worldXmin = worldGeom->mu2eOriginInWorld().x() - worldGeom->halfLengths()[0] + deltaX;
    _worldXmax = worldGeom->mu2eOriginInWorld().x() + worldGeom->halfLengths()[0] - deltaX;
    _worldYmin = env->ymin() + deltaX;
    _worldYmax = env->ymax() - deltaX;
    _worldZmin = worldGeom->mu2eOriginInWorld().z() - worldGeom->halfLengths()[2] + deltaX;
    _worldZmax = worldGeom->mu2eOriginInWorld().z() + worldGeom->halfLengths()[2] - deltaX;

    if (_refPointChoice == "TRACKER")
    {
      _cosmicReferencePointInMu2e = Hep3Vector(detsys->getOrigin().x(),
                                               _refY0, detsys->getOrigin().z());
    }
    else if (_refPointChoice == "EXTMONFNAL")
    {
      GeomHandle<ExtMonFNAL::ExtMon> extMonFNAL;
      _cosmicReferencePointInMu2e =
          Hep3Vector(extMonFNAL->detectorCenterInMu2e().x(),
                      _refY0, extMonFNAL->detectorCenterInMu2e().z());
    }
    else if (_refPointChoice == "CALO")
    {
      GeomHandle<Calorimeter> calorimeter;
      _cosmicReferencePointInMu2e = Hep3Vector(detsys->getOrigin().x(),
                                                _refY0, calorimeter->disk(0).geomInfo().origin().z());
    }
    else if (_refPointChoice == "UNDEFINED")
      _cosmicReferencePointInMu2e = Hep3Vector(0., _refY0, 0.);

    bool passed = false;
    while (!passed)
    {

      GenParticleCollection pretruth;
      genEvent(pretruth);

      for (unsigned int i = 0; i < pretruth.size(); ++i)
      {
        GenParticle particle = pretruth[i];
        const HepLorentzVector mom4 = particle.momentum();

        const Hep3Vector position(
            wrapvar(particle.position().x() + _flat(), _targetBoxXmin - _showerAreaExtension, _targetBoxXmax + _showerAreaExtension) + _cosmicReferencePointInMu2e.x(),
            particle.position().y(),
            wrapvar(particle.position().z() + _flat(), _targetBoxZmin - _showerAreaExtension, _targetBoxZmax + _showerAreaExtension) + _cosmicReferencePointInMu2e.z()
        );

        if (_projectToTargetBox)
        {

          _targetBoxIntersections.clear();
          // std::cout << "Position " << position.x() << " " << position.y() << " " << position.z() << std::endl;
          // std::cout << "Direction " << mom4.vect().x() << " " << mom4.vect().y() << " " << mom4.vect().z() << std::endl;

          // std::cout << "Target box " << _targetBoxXmin << " " << _targetBoxXmax << " " << _targetBoxYmin << " " << _targetBoxYmax << " " << _targetBoxZmin << " " << _targetBoxZmax << std::endl;

          calIntersections(position, mom4.vect(),
                            _targetBoxIntersections, _targetBoxXmin - _buffer, _targetBoxXmax + _buffer,
                            _targetBoxYmin - _buffer, _targetBoxYmax + _buffer, _targetBoxZmin - _buffer, _targetBoxZmax + _buffer);

          if (_targetBoxIntersections.size() > 0)
          {
            _worldIntersections.clear();
            calIntersections(position, mom4.vect(), _worldIntersections,
                _worldXmin, _worldXmax, _worldYmin, _worldYmax, _worldZmin, _worldZmax);

            if (_worldIntersections.size() > 0) {
                int idx = 0;
                float closestDistance = distance(_worldIntersections.at(0), position);
                for (unsigned i = 0; i < _worldIntersections.size(); ++i) {
                    if (distance(_worldIntersections.at(i), position) < closestDistance) {
                        idx = i;
                        closestDistance = _targetBoxIntersections.at(idx).y();
                    }
                }

              const Hep3Vector projectedPos = _worldIntersections.at(idx);

              genParts.push_back(GenParticle(static_cast<PDGCode::type>(particle.pdgId()),
                          GenId::cosmicCORSIKA, projectedPos, mom4,
                          particle.time() + _tOffset));
            }
        }
        } else {
          genParts.push_back(GenParticle(static_cast<PDGCode::type>(particle.pdgId()),
                      GenId::cosmicCORSIKA, position, mom4,
                      particle.time() + _tOffset));
        }
      }

      if (genParts.size() != 0) {
        passed = true;
      }

    }

  }

}// end namespace
