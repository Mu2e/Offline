// NotchManager.cc
// This is the definitions file for the NotchManager class, a hopefully
// unobtrusive way to add notches/holes to Mu2e building geometry objects.
// David Norvil Brown, U. Louisville, December 2017

#include "GeometryService/inc/NotchManager.hh"
#include "GeomPrimitives/inc/Notch.hh"
#include <sstream>

void mu2e::NotchManager::loadNotches( const SimpleConfig& config ) {
  if ( hasLoaded_ ) return;

  // Get all the variables in SimpleConfig
  std::vector<std::string> variables;
  config.getNames(variables);

  // Another vector for just the names of the notch-containing parts
  // and how many notches each
  std::vector<std::string>  usedVariables;
  std::string searchString(".Notch.numberOfNotches");
  for ( auto theVar : variables ) {
    std::size_t pos = 0;
    if ( (pos = theVar.find(searchString)) != std::string::npos ) {
      // Assume variable names go like this:  
      // part_name.Notch.notch_property_name.number   OR
      // part_name.Notch.numberOfNotches  
      // We trigger on the latter here 
      std::string partName = theVar.substr(0,pos);
      usedVariables.push_back(partName);
    } // end of if Notch found
  } // end of loop over variables to find those with Notch specs

  // Have list of parts with notches now.  Build the notch map.
  for ( auto partName : usedVariables ) {
    // For each part that has a notch or notches, we have to first
    // ask how many notches for the part, then loop and for each notch,
    // find the dimensions, position, and orientation.  Dimensions are in
    // the form of a vector of 3 doubles.  Position is in the form of a
    // Hep3Vector (entered in SimpleConfig as a vector of 3 doubles).  
    // Orientation is a 3-character string.
    std::ostringstream tmpName1;
    tmpName1 << partName << ".Notch.numberOfNotches";
    int nNotches = config.getInt(tmpName1.str(),0);
    // Now create the vector of notches for this part
    std::vector<Notch> tmpVecNotch;
    tmpVecNotch.reserve(nNotches);
    for ( int iNotch = 0; iNotch < nNotches; iNotch++ ) {
      // Get the (half) dimensions
      std::ostringstream tmpName2;
      tmpName2 << partName << ".Notch.halfDims." << iNotch+1;
      std::vector<double> tmpDims;
      config.getVectorDouble(tmpName2.str(),tmpDims);
      // Get the position
      std::ostringstream tmpName3;
      tmpName3 << partName << ".Notch.position." << iNotch+1;
      CLHEP::Hep3Vector tmpPosition = config.getHep3Vector(tmpName3.str());
      // Get the orientation
      std::ostringstream tmpName4;
      tmpName4 << partName << ".Notch.orientation." << iNotch+1;
      std::string tmpOri = config.getString(tmpName4.str());

      // Make the notch from this info
      Notch tmpNotch(tmpDims, tmpPosition, tmpOri);

      // Put the notch in the vector
      tmpVecNotch.push_back(tmpNotch);
    } // Loop over notches for part
    // Now get the associate volumeName for the part
    std::ostringstream tmpName5;
    tmpName5 << partName << ".name";
    std::string volName = config.getString(tmpName5.str());
    // Now put the volume in the map
    theMap_.emplace( volName, tmpVecNotch );
  } // end of loop to build map

  // Information cached as long as this instance shall live
  hasLoaded_ = true;
} // end def of NotchManager::loadNotches

const bool mu2e::NotchManager::hasNotches( const std::string& part ) const {
  if ( theMap_.empty() ) return false;
  if ( theMap_.find( part ) == theMap_.end() ) return false;
  return true;  // Am I forgetting any scenarios?
} // end def of NotchManager::hasNotches

const std::vector<mu2e::Notch>& mu2e::NotchManager::getNotchVector ( const std::string& part ) const {

  auto whatPart = theMap_.find( part );
  if ( whatPart != theMap_.end() ) return whatPart->second;

  // User should have used above function to check for this possibility, but
  // let's assume them didn't and return an empty vector.
  return emptyVec_;

} // end def of NotchManager::getNotchVector
