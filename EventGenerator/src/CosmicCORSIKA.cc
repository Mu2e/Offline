////////////////////////////////////////////////////////////////////////
/// \file  CosmicCORSIKA_module.cc
/// \brief Generator for cosmic-ray secondaries based on pre-generated CORSIKA shower databases.
///
/// \author  Matthew.Bass@physics.ox.ac.uk
////////////////////////////////////////////////////////////////////////

#include "EventGenerator/inc/CosmicCORSIKA.hh"

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

namespace mu2e {

  CosmicCORSIKA::CosmicCORSIKA(art::Run &run,
                              const SimpleConfig &config,
                              CLHEP::HepRandomEngine &engine)
      : _refY0(config.getDouble("cosmicCORSIKA.refY0", 0.)), // mm
        _projectToTargetBox(config.getBool("cosmicCORSIKA.projectToTargetBox", false)),
        _tOffset(config.getDouble("cosmicCORSIKA.timeOffset", 0.)),
        _buffer(config.getDouble("cosmicCORSIKA.buffer", 0.)),
        _showerAreaExtension(config.getDouble("cosmicCORSIKA.showerAreaExtension", 2000.)), // mm
        _refPointChoice(config.getString("cosmicCORSIKA.refPoint", "UNDEFINED")),
        _targetBoxXmin(config.getDouble("cosmicCORSIKA.targetBoxXmin", -1000)), // mm
        _targetBoxXmax(config.getDouble("cosmicCORSIKA.targetBoxXmax", 1000)),  // mm
        _targetBoxYmin(config.getDouble("cosmicCORSIKA.targetBoxYmin", -1000)), // mm
        _targetBoxYmax(config.getDouble("cosmicCORSIKA.targetBoxYmax", 1000)),  // mm
        _targetBoxZmin(config.getDouble("cosmicCORSIKA.targetBoxZmin", -1000)), // mm
        _targetBoxZmax(config.getDouble("cosmicCORSIKA.targetBoxZmax", 1000)),  // mm
        _randomXZshift(config.getDouble("cosmicCORSIKA.randomXZshift", 500)), // mm
        _flat(engine, -_randomXZshift, _randomXZshift),
        _fluxConstant(config.getDouble("cosmicCORSIKA.fluxConstant", 1.8e4))
  {
    config.getVectorString("cosmicCORSIKA.showerInputFiles", _showerInputFiles);

    _showerInputs = _showerInputFiles.size();

    if (_refY0 == 0.)
      mf::LogInfo("CosmicCORSIKA") << "Using 0. for _refY0!";

    _fileIndex = 0;
    mf::LogInfo("CosmicCORSIKA") << "Opening: " << _showerInputFiles[_fileIndex] << std::endl;
    in = fopen(_showerInputFiles[_fileIndex].c_str(), "r");
    fread(&_garbage, 4, 1, in); // skipping the first (and useless) information

  }

  int CosmicCORSIKA::corsikaToHepevtID(const int corsikaID)
  {

    switch (corsikaID)
    {

    case 1:
      return 22; // gamma
    case 2:
      return -11; // e+
    case 3:
      return 11; // e-
    case 5:
      return -13; // mu+
    case 6:
      return 13; // mu-
    case 7:
      return 111; // pi0
    case 8:
      return 211; // pi+
    case 9:
      return -211; // pi-
    case 10:
      return 130; // K0_L
    case 11:
      return 321; // K+
    case 12:
      return -321; // K-
    case 13:
      return 2112; // n
    case 14:
      return 2212; // p
    case 15:
      return -2212; // pbar
    case 16:
      return 310; // K0_S
    case 17:
      return 221; // eta
    case 18:
      return 3122; // Lambda
    case 19:
      return 3222; // Sigma+
    case 20:
      return 3212; // Sigma0
    case 21:
      return 3112; // Sigma-
    case 22:
      return 3322; // Cascade0
    case 23:
      return 3312; // Cascade-
    case 24:
      return 3334; // Omega-
    case 25:
      return -2112; // nbar
    case 26:
      return -3122; // Lambdabar
    case 27:
      return -3112; // Sigma-bar
    case 28:
      return -3212; // Sigma0bar
    case 29:
      return -3222; // Sigma+bar
    case 30:
      return -3322; // Cascade0bar
    case 31:
      return -3312; // Cascade+bar
    case 32:
      return -3334; // Omega+bar

    case 50:
      return 223; // omega
    case 51:
      return 113; // rho0
    case 52:
      return 213; // rho+
    case 53:
      return -213; // rho-
    case 54:
      return 2224; // Delta++
    case 55:
      return 2214; // Delta+
    case 56:
      return 2114; // Delta0
    case 57:
      return 1114; // Delta-
    case 58:
      return -2224; // Delta--bar
    case 59:
      return -2214; // Delta-bar
    case 60:
      return -2114; // Delta0bar
    case 61:
      return -1114; // Delta+bar
    case 62:
      return 10311; // K*0
    case 63:
      return 10321; // K*+
    case 64:
      return -10321; // K*-
    case 65:
      return -10311; // K*0bar
    case 66:
      return 12; // nu_e
    case 67:
      return -12; // nu_ebar
    case 68:
      return 14; // nu_mu
    case 69:
      return -14; // nu_mubar

    case 116:
      return 421; // D0
    case 117:
      return 411; // D+
    case 118:
      return -411; // D-bar
    case 119:
      return -421; // D0bar
    case 120:
      return 431; // D+_s
    case 121:
      return -431; // D-_sbar
    case 122:
      return 441; // eta_c
    case 123:
      return 423; // D*0
    case 124:
      return 413; // D*+
    case 125:
      return -413; // D*-bar
    case 126:
      return -423; // D*0bar
    case 127:
      return 433; // D*+_s
    case 128:
      return -433; // D*-_s

    case 130:
      return 443; // J/Psi
    case 131:
      return -15; // tau+
    case 132:
      return 15; // tau-
    case 133:
      return 16; // nu_tau
    case 134:
      return -16; // nu_taubar

    case 137:
      return 4122; // Lambda+_c
    case 138:
      return 4232; // Cascade+_c
    case 139:
      return 4132; // Cascade0_c
    case 140:
      return 4222; // Sigma++_c
    case 141:
      return 4212; // Sigma+_c
    case 142:
      return 4112; // Sigma0_c
    case 143:
      return 4322; // Cascade'+_c
    case 144:
      return 4312; // Cascade'0_c
    case 145:
      return 4332; // Omega0_c
    case 149:
      return -4122; // Lambda-_cbar
    case 150:
      return -4232; // Cascade-_cbar
    case 151:
      return -4132; // Cascade0_cbar
    case 152:
      return -4222; // Sigma--_cbar
    case 153:
      return -4212; // Sigma-_cbar
    case 154:
      return -4112; // Sigma0_cbar
    case 155:
      return -4322; // Cascade'-_cbar
    case 156:
      return -4312; // Cascade'0_cbar
    case 157:
      return -4332; // Omega0_cbar
    case 161:
      return 4224; // Sigma*++_c
    case 162:
      return 1214; // Sigma*+_c
    case 163:
      return 4114; // Sigma*0_c

    case 171:
      return -4224; // Sigma*--_cbar
    case 172:
      return -1214; // Sigma*-_cbar
    case 173:
      return -4114; // Sigma*0_cbar
    case 176:
      return 511; // B0
    case 177:
      return 521; // B+
    case 178:
      return -521; // B-bar
    case 179:
      return -511; // B0bar
    case 180:
      return 531; // B0_s
    case 181:
      return -531; // B0_sbar
    case 182:
      return 541; // B+_c
    case 183:
      return -541; // B-_cbar
    case 184:
      return 5122; // Lambda0_b
    case 185:
      return 5112; // Sigma-_b
    case 186:
      return 5222; // Sigma+_b
    case 187:
      return 5232; // Cascade0_b
    case 188:
      return 5132; // Cascade-_b
    case 189:
      return 5332; // Omega-_b
    case 190:
      return -5112; // Lambda0_bbar
    case 191:
      return -5222; // Sigma+_bbar
    case 192:
      return -5112; // Sigma-_bbar
    case 193:
      return -5232; // Cascade0_bbar
    case 194:
      return -5132; // Cascade+_bbar
    case 195:
      return -5332; // Omega+_bbar
    }

    return 0;
  }

  const unsigned int CosmicCORSIKA::getNumShowers() {
    return _primaries;
  }

  const double CosmicCORSIKA::getLiveTime()
  {
    const float area = (_targetBoxXmax + 2 * _showerAreaExtension - _targetBoxXmin) * (_targetBoxZmax + 2 * _showerAreaExtension - _targetBoxZmin) * 1e-6; //m^2
    const float eslope = -2.7;
    const float lowE = 1.3;     // GeV
    const float highE = 100000; // GeV
    const double EiToOneMinusGamma = pow(lowE, 1 + eslope);
    const double EfToOneMinusGamma = pow(highE, 1 + eslope);
    return _primaries / (M_PI * area * _fluxConstant * (EfToOneMinusGamma - EiToOneMinusGamma) / (1. + eslope));
  }

  CosmicCORSIKA::~CosmicCORSIKA(){
    fclose(in);
  }

  bool CosmicCORSIKA::pointInBox(double x, double y, double x0, double y0,
                              double x1, double y1)
  {
    bool ret = false;
    if ((x >= x0) && (x <= x1) && (y >= y0) && (y <= y1))
    {
      ret = true;
    }
    return ret;
  }

  double CosmicCORSIKA::distance(const CLHEP::Hep3Vector &u, const CLHEP::Hep3Vector &v)
  {
    return safeSqrt((u.x() - v.x()) * (u.x() - v.x()) +
                    (u.y() - v.y()) * (u.y() - v.y()) +
                    (u.z() - v.z()) * (u.z() - v.z()));
  }

  void CosmicCORSIKA::calIntersections(Hep3Vector orig, Hep3Vector dir,
                                  std::vector<CLHEP::Hep3Vector> &intersections, double xMin, double xMax,
                                  double yMin, double yMax, double zMin, double zMax)
  {
    // roof: _targetBoxYmax, _targetBoxXmin, _targetBoxXmax, _targetBoxZmin, _targetBoxZmax
    // skip projection if the particle goes parallely to the plane
    if (dir.y() != 0.)
    {
      const double t = (yMax - orig.y()) / dir.y();
      const double x1 = dir.x() * t + orig.x();
      const double z1 = dir.z() * t + orig.z();
      // std::cout << "x1 " << x1 << " " << z1 << std::endl;

      if (pointInBox(x1, z1, xMin, zMin, xMax, zMax))
      {
        intersections.push_back(Hep3Vector(x1, yMax, z1));
      }
    }

    //east: zMin, xMin, xMax, yMin, yMax
    if (dir.z() != 0.)
    {
      const double t = (zMin - orig.z()) / dir.z();
      const double x1 = dir.x() * t + orig.x();
      const double y1 = dir.y() * t + orig.y();
      if (pointInBox(x1, y1, xMin, yMin, xMax, yMax))
      {
        intersections.push_back(Hep3Vector(x1, y1, zMin));
      }
    }

    //west: zMax, xMin, xMax, yMin, yMax
    if (dir.z() != 0.)
    {
      const double t = (zMax - orig.z()) / dir.z();
      const double x1 = dir.x() * t + orig.x();
      const double y1 = dir.y() * t + orig.y();
      if (pointInBox(x1, y1, xMin, yMin, xMax, yMax))
      {
        intersections.push_back(Hep3Vector(x1, y1, zMax));
      }
    }

    //south: xMin, yMin, yMax, zMin, zMax
    if (dir.x() != 0.)
    {
      const double t = (xMin - orig.x()) / dir.x();
      const double z1 = dir.z() * t + orig.z();
      const double y1 = dir.y() * t + orig.y();
      if (pointInBox(z1, y1, zMin, yMin, zMax, yMax))
      {
        intersections.push_back(Hep3Vector(xMin, y1, z1));
      }
    }

    //north: xMax, yMin, yMax, zMin, zMax
    if (dir.x() != 0.)
    {
      const double t = (xMax - orig.x()) / dir.x();
      const double z1 = dir.z() * t + orig.z();
      const double y1 = dir.y() * t + orig.y();
      if (pointInBox(z1, y1, zMin, yMin, zMax, yMax))
      {
        intersections.push_back(Hep3Vector(xMax, y1, z1));
      }
    }
  }

  double CosmicCORSIKA::wrapvar( const double var, const double low, const double high){
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
            double id = block[k];
            //cout << "	" << k << " " << block[k] << " " << block[k+1] << " " << block[k+2] << " " << block[k+3] << " " << block[k+4] << " " << block[k+5] << " " << block[k+6] << endl;	// for verification purposes only

            if ((int)(id / 1000) != 0)
            {
              if (id < 75000)
              {
                // double dxyz[3] = {-block[k + 1], -block[k + 3], block[k + 2]};
                // double xyz[3] = {-block[k + 4], FDHalfHeight, block[k + 5]};
                // if (CheckTPCIntersection(xyz, dxyz)) {
                if (block[k + 1] != 0)
                {
                  id = (int)(id / 1000);
                  const int pdgId = corsikaToHepevtID(id);

                  const double P_x = block[k + 2] * _GeV2MeV;
                  const double P_y = -block[k + 3] * _GeV2MeV;
                  const double P_z = block[k + 1] * _GeV2MeV;
                  const double x = block[k + 5];
                  const double z = -block[k + 4];

                  GlobalConstantsHandle<ParticleDataTable> pdt;
                  ParticleDataTable const &pdt_ = *pdt;
                  const double m = pdt_.particle(pdgId).ref().mass(); // to MeV

                  const double energy = safeSqrt(P_x * P_x + P_y * P_y + P_z * P_z + m * m);

                  const Hep3Vector position(x, 0, z);
                  const HepLorentzVector mom4(P_x, P_y, P_z, energy);

                  const double particleTime = block[k + 6];

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
                double closestDistance = distance(_worldIntersections.at(0), position);
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
