//
// Create a disk and fills it with crystals. 
// We assume that the real crystal position is close to the ideal one. If not, 
//     this module needs to be completely rewritten
//
// Original author B Echenard
//

// Notes: CrystalMap tesselates a plane with square crystals and let you know which 
//        ones are fully contained inside an annulus

#include "CalorimeterGeom/inc/CrystalCondReader.hh"
#include "cetlib_except/exception.h"

#include <sstream>
#include <fstream>

namespace mu2e {

      
      CrystalCondReader::CrystalCondReader(const std::string& filename) : dataPosition_(),dataSize_()
      { 
	  read(filename);
      }
      

      int CrystalCondReader::nCrystal() const {return dataPosition_.size();}
     
      void CrystalCondReader::read(const std::string& filename)
      {   
          std::ifstream myfile(filename.c_str());
          
          //put error message here instead
          if (! myfile.is_open()) throw cet::exception("CrystalCondReader") << " unknown file "<<filename<<"\n";

          int crId;
          std::vector<double> pos,size;
          std::string line,substring;
          while(std::getline(myfile,line))
          {                 
             pos.clear(); 
             size.clear();
             std::istringstream iss(line);
             iss >> substring; crId  = atoi(substring.c_str()); 
             iss >> substring; pos.push_back(atof(substring.c_str())); 
             iss >> substring; pos.push_back(atof(substring.c_str())); 
             iss >> substring; size.push_back(atof(substring.c_str())); 
             iss >> substring; size.push_back(atof(substring.c_str())); 
             dataPosition_[crId] = pos;
             dataSize_[crId] = size;           
          }

          myfile.close();
      }


      //-----------------------------------------------------------------------------
      const CLHEP::Hep3Vector CrystalCondReader::position(int crId) const
      {
          auto iter = dataPosition_.find(crId);
          if (iter == dataPosition_.end())
             throw cet::exception("CrystalCondReader") << " unknown crystal id "<<crId<<"\n";

          return CLHEP::Hep3Vector(iter->second.at(0),iter->second.at(1),0.0);          
      }

      //-----------------------------------------------------------------------------
      const CLHEP::Hep3Vector CrystalCondReader::size(int crId) const
      {
          auto iter = dataSize_.find(crId);
          if (iter == dataSize_.end())
              throw cet::exception("CrystalCondReader") << " unknown crystal id "<<crId<<"\n";
          return CLHEP::Hep3Vector(iter->second.at(0),iter->second.at(1),0.0);          
      }


      const void CrystalCondReader::print(std::ostream &os) const
      {
           os<<"Crystal position"<<std::endl;
           for (const auto iter : dataPosition_) 
             os<<iter.first<<" "<<iter.second.at(0)<<" "<<iter.second.at(1)<<std::endl;
             
           os<<"Crystal size"<<std::endl;
           for (const auto iter : dataSize_) 
             os<<iter.first<<" "<<iter.second.at(0)<<" "<<iter.second.at(1)<<std::endl;
      }  

}

