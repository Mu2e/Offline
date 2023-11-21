// NotchHoleManager.cc
// This is the definitions file for the NotchManager class, a hopefully
// unobtrusive way to add notches/holes to Mu2e building geometry objects.
// David Norvil Brown, U. Louisville, December 2017

#include "Offline/GeometryService/inc/NotchHoleManager.hh"
#include "Offline/GeomPrimitives/inc/Notch.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/GeomPrimitives/inc/Hole.hh"
#include <sstream>

void mu2e::NotchHoleManager::loadNotches( const SimpleConfig& config ) {
  if ( hasLoaded_ ) return;

  // Get all the variables in SimpleConfig
  std::vector<std::string> variables;
  config.getNames(variables);

  // Another vector for just the names of the notch-containing parts
  // and how many notches each
  std::vector<std::string>  usedVariables;
  std::vector<std::string>  usedVariablesHH;
  std::string searchString(".Notch.numberOfNotches");
  std::string searchString1(".Hole.numberOfHoles");
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
    theMap_.emplace( volName, tmpVecNotch);
   // end of loop to build map //theMap has been called twice, first loadNotches and loadholes, err?
  }

  // Information cached as long as this instance shall live
  hasLoaded_ = true;
} // end def of NotchHoleManager::loadHoles



    //--------------
    //----load holes-----------
    //--------------

void mu2e::NotchHoleManager::loadHoles( const SimpleConfig& config ) {
  if ( hasLoadedH_ ) return;  //check if err

  // Get all the variables in SimpleConfig for holes
  std::vector<std::string> variablesH;
  config.getNames(variablesH);

  // Another vector for just the names of the hole-containing parts
  // and how many holes each
  std::vector<std::string>  usedVariablesH;
  std::string searchString2(".Hole.numberOfHoles");
  for ( auto theVarH : variablesH ) {
    std::size_t pos = 0; // check if err
    if ( (pos = theVarH.find(searchString2)) != std::string::npos ) {
      // Assume variable names go like this:
      // part_name.Hole.numberOfHoles
      // We trigger on the latter here
      std::string partNameH = theVarH.substr(0,pos);
      usedVariablesH.push_back(partNameH);
    } // end of if Holes found
  } // end of loop over variables to find those with Holes

  // Have list of parts with holes now.  Build the hole map.
  for ( auto partNameH : usedVariablesH ) {
    // For each part that has a hole or holes, we have to first
    // ask how many holes for the part, then loop and for each hole,
    // find the radii,lengths, position, and orientation.  Dimensions (radii & lengths) are in
    // the form of a vector of doubles.  Position is in the form of a
    // Hep3Vector (entered in SimpleConfig as a vector of 3 doubles).
    // Orientation is a 3-character string.
    std::ostringstream tmpName11;
    tmpName11 << partNameH << ".Hole.numberOfHoles";
    int nHoles = config.getInt(tmpName11.str(),0);
    // Now create the vector of holes for this part
    std::vector<Hole> tmpVecHole;
    tmpVecHole.reserve(nHoles);
    for ( int iHole = 0; iHole < nHoles; iHole++ ) {
      // Get the radius
      std::ostringstream tmpName22;
      tmpName22 << partNameH << ".Hole.radius." << iHole+1;
      double tmpRadius=config.getDouble(tmpName22.str());
      // Get the half-lengths
      std::ostringstream tmpName33;
      tmpName33 << partNameH << ".Hole.halfLength." << iHole+1;
      double tmphalfLength=config.getDouble(tmpName33.str());
      // Get the position
      std::ostringstream tmpName44;
      tmpName44 << partNameH << ".Hole.position." << iHole+1;
      CLHEP::Hep3Vector tmpPositionH = config.getHep3Vector(tmpName44.str());
      // Get the orientation
      std::ostringstream tmpName55;
      tmpName55 << partNameH << ".Hole.orientation." << iHole+1;
      std::string tmpOriH = config.getString(tmpName55.str());
      // Make the hole from this info
      Hole tmpHole(tmpRadius,tmphalfLength,tmpPositionH,tmpOriH);

      // Put the hole in the vector
      tmpVecHole.push_back(tmpHole);
    } // Loop over holes for part


    // Now get the associate volumeName for the part
    std::ostringstream tmpName66;
    tmpName66 << partNameH << ".name";
    std::string volNameH = config.getString(tmpName66.str());
    // Now put the volume in the map
    theMapH_.emplace( volNameH, tmpVecHole);
  } // end of loop to build map

  // Information cached as long as this instance shall live
  hasLoadedH_ = true;
} // end def of NotchHoleManager::loadHoles


const bool mu2e::NotchHoleManager::hasNotches( const std::string& part ) const {
  if ( theMap_.empty() ) return false;
  if ( theMap_.find( part ) == theMap_.end() ) return false;
  return true;  // Am I forgetting any scenarios?
} // end def of NotchHoleManager::hasNotches

const bool mu2e::NotchHoleManager::hasHoles( const std::string& partH ) const {
  if ( theMapH_.empty() ) return false;
  if ( theMapH_.find( partH ) == theMapH_.end() ) return false;
  return true;
} // end def of NotchHoleManager::hasHoles

const std::vector<mu2e::Notch>& mu2e::NotchHoleManager::getNotchVector ( const std::string& part ) const {

  auto whatPart = theMap_.find( part );
  if ( whatPart != theMap_.end() ) return whatPart->second;
  // User should have used above function to check for this possibility, but
  // let's assume them didn't and return an empty vector.
  return emptyVec_;

} // end def of NotchHoleManager::getNotchVector

const std::vector<mu2e::Hole>& mu2e::NotchHoleManager::getHoleVector ( const std::string& partH ) const {

  auto whatPartH = theMapH_.find( partH );
  if ( whatPartH != theMapH_.end() ) return whatPartH->second;
  return emptyVecH_;
} // end def of NotchHoleManager::getHoleVector
