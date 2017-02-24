// Original author: Andrei Gaponenko, 2012

#include "Sources/inc/ExtMonFNALMARSUtils.hh"

#include <map>
#include <cmath>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/LorentzVector.h"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    std::ostream& operator<<(std::ostream& os, const MARSParticle& mp) {
      return os<<"MARSParticle("<<mp.protonNumber
               <<", pid="<<mp.pid
               <<", kE="<<mp.kineticEnergy
               <<", w="<<mp.weight
               <<", tof="<<mp.tof
               <<")";
    }

    //================================================================
    bool readMARSLine(std::istream& file, MARSParticle& res) {
      std::string line;
      std::getline(file, line);

      if(line.empty()) {
        res.protonNumber = -1U;
        return false;
      }
      else {
        std::istringstream is(line);

        if(!(is>>res.protonNumber
             >>res.pid>>res.kineticEnergy>>res.weight
             >>res.x>>res.y>>res.z
             >>res.dcx>>res.dcy>>res.dcz
             >>res.tof)
           ) {

          throw cet::exception("BADINPUT")<<" parseMARSLine(): error processing line: "<<line<<"\n";
        }
        return true;
      }
    }

    //================================================================
    GenParticle MARSMu2eConverter::marsToMu2eParticle(const MARSParticle& mp) {
      PDGCode::type pdgId = marsToMu2eParticleCode(mp.pid);

      //try {
      const double mass = pdt_->particle(pdgId).ref().mass().value();
      // }
      //// ParticleDataTable throws a string?!
      // catch(cet::exception& e) { throw; }

      const double energy = mass + marsToMu2eEnergy(mp.kineticEnergy);
      const double p3mag = sqrt((energy-mass)*(energy+mass));

      const CLHEP::HepLorentzVector p4(mp.dcx * p3mag,
                                       mp.dcy * p3mag,
                                       mp.dcz * p3mag,
                                       energy
                                       );

      return GenParticle(pdgId,
                         GenId::MARS,
                         marsToMu2ePosition(mp.x, mp.y, mp.z),
                         p4,
                         marsToMu2eTime(mp.tof)
                         );
    }

    //================================================================
    // returns pdgId
    PDGCode::type MARSMu2eConverter::marsToMu2eParticleCode(int marsPID) {
      static std::map<int,int> table;
      if(table.empty()) {
        table[1] = +2212; // proton
        table[2] = +2112; // neutron
        table[3] = +211; // pi+
        table[4] = -211; // pi-
        table[5] = +321; // K+
        table[6] = -321; // K-
        table[7] = -13; // mu+
        table[8] = +13; // mu-
        table[9] =  22; // gamma
        table[10] = -11; // e-
        table[11] = +11; // e+
        table[12] = -2212; // p- (antiproton)
        table[13] = 111; // pi0
        table[14] = 1000010020; // deutron
        table[15] = 1000010030; // H3
        table[16] = 1000020030; // He3
        table[17] = 1000020040; // He4
        table[18] = +14; // numu
        table[19] = -14; // numubar
        table[20] = +12; // nue
        table[21] = -12; // nuebar
        table[22] = 130; // K0L
        table[23] = 310; // K0S
        //table[24] = ; // K0
        //table[25] = ; // K0bar
        table[26] = +3122; // Lambda
        table[27] = -3122; // Lambdabar
        table[28] = +3222; // Sigma+
        table[29] = +3212; // Sigma0
        table[30] = +3112; // Sigma-
        table[31] = -2112; // nbar
        table[32] =  3322; // Xi0
        table[33] =  3312; // Xi-
        table[34] =  3334; // Omega-
        table[35] = -3222; // anti-Sigma-
        table[36] = -3212; // anti-Sigma0
        table[37] = -3112; // anti-Sigma+
        table[38] = -3322; // anti-Xi0
        table[39] = -3312; // anti-Xi+
        table[40] = -3334; // anti-Omega+
      }
      std::map<int,int>::const_iterator ip = table.find(marsPID);
      if(ip == table.end()) {
        throw cet::exception("BADINPUT")<<" marsToMu2eParticleCode(): unknonw MARS particle code "<<marsPID<<"\n";
      }
      return PDGCode::type(ip->second);
    }
  } // namespace ExtMonFNAL
} // namespace mu2e
