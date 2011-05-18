#ifndef G4Helper_AntiLeakRegistry_hh
#define G4Helper_AntiLeakRegistry_hh
//
// An anti-leak system to aid in using G4 from the Mu2e framework.
//
// $Id: AntiLeakRegistry.hh,v 1.4 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
//
// Original author Rob Kutschke
//
// This header defines an object registry that is used to avoid memory leaks in G4.
// The API for G4 often requires that user code create a new object on the heap and
// pass a pointer to that object to G4.  In a few cases, G4 takes ownership of the
// newly created object and deletes it at the appropriate time.  In most cases,
// however, G4 does not take ownership of the object and it is the responsibility of
// the user code to delete the object when appropriate.  In most cases, the appropriate
// time is either at the end of a G4 run or at the end of the job; for Mu2e these
// cases are not distinguished because we do not intend to run multiple G4 runs within
// one framework job.
//
// One diffculty with this API is that the author of the Mu2e code just has to know whether
// or not G4 will take ownership of a particular object.  The anti-leak registry does not
// help with that problem; you still just have to know.  The anti-leak registry does
// provide a mechanism to register an object so that it will be deleted when the registry
// goes out of scope at the end of the job.  If the correct time to delete your object is
// not at the end of job, then you will need to find another way to manage its lifetime.
//
// This registry is accessed via the G4Helper service.
//
// Here is an example of a call to G4 that leaks memory:
//
//    G4LogicalVolume lvol = ...;
//    lvol->SetVisAttributes(  new G4VisAttributes(true, color) );
//
// The logical volume does not take ownership of the G4VisAttributes object.
// The recommended practice for this situation is:
//
//      AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
//      G4LogicalVolume v = ...;
//
//      lvol->SetVisAttributes( reg.add(G4VisAttributes(true, color)) );
// or
//      lvol->SetVisAttributes( reg.add( new G4VisAttributes(true, color)) );
// or
//      G4VisAttributes* visAtt = new G4VisAttributes(true, color);
//      visAtt->Set ... ; // modify some property.
//      lvol->SetVisAttributes( reg.add( visAtt ) );
//
// There are two distinct add methods, one of which takes an object by const reference
// and one of which takes an object by pointer.  The first add method will make a copy of
// its argument, placing that copy in a std::list of the correct type and will return a
// pointer to the object in the list. The type std::list was chosen ( not, for example,
// std::vector ) because subsequent insertions into the list do not invalidate the pointer
// that was returned.
//
// The second add method puts a shared pointer to the object into a list of the correct type.
// So no copy is necessary and the returned pointer points directly to the new'ed object.
//
// The first method has the advantage that the interface is cleaner but it has the
// disadvantage that it requires a copy, which has two sorts of potential problems. The big
// problem is if the type T is not copyable ( or should not be copyable because copying
// breaks smething); the solution is that must use add a pointer to the object.  The lesser
// problem is that the copy operation is wasteful of time and memory; however the add
// operations are done once at the start of the job, which makes the waste a minor consideration;
// with C++0X, which supports move aware objects,  we may be able to remove this objection too.
//

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <typeinfo>

#include "boost/shared_ptr.hpp"

namespace mu2e
{

  class AntiLeakRegistry {

  private:

  // Base class that defines an interface obeyed by the concrete list classes.
  // The concrete classes are held in the collection as pointer to base.
  class ObjectListBase{
  public:
    virtual ~ObjectListBase(){}
    virtual size_t size() const = 0;
  };

  // Template for a concrete class to hold a list of objects of templated type.
  template< class T>
  class ObjectList: public ObjectListBase{

  public:
    ObjectList():ObjectListBase(){}
    virtual ~ObjectList(){}
    size_t size() const { return list.size(); }

    std::list<T> list;

  };

  public:
    AntiLeakRegistry():
      _mapOfLists(),
      _maxNameLength(0){}

    ~AntiLeakRegistry(){}

    // Place a copy of the argument in the list and return a pointer to the
    // object in the list.
    template< class T>
    T* add( T const& p ){

      // Get the correct list
      ObjectList<T>* olt = findOrCreateList<T>();

      // Add the object to be managed; this does a copy.
      olt->list.push_back(p);

      // Return the address of the copy of the object.
      return &olt->list.back();
    }

    // Take ownership of the object pointed to by p.  Place a shared pointer to the
    // object into the list.  Return the original pointer.
    template< class T>
    T* add( T* p ){

      // Object will be stored by shared pointer so that deletion is automatic.
      typedef boost::shared_ptr<T> SPtr;

      // Get the correct list
      ObjectList<SPtr>* olt = findOrCreateList<SPtr>();

      // Add shared pointer to the object.
      olt->list.push_back(SPtr(p));

      // Return the address of the object.
      //return olt->list.back().get();
      return p;
    }

    // Clear everything; calls all destructors.
    void clear();

    // Print information about types and numbers of saved objects.
    void print( std::ostream& ost = std::cout ) const;

    // How many types of objects have we stored?
    size_t size() const { return _mapOfLists.size(); }

  private:

    // A helper function to find or create the relevant list.
    template< class T>
    ObjectList<T>* findOrCreateList(){

      // This will be the key in the map.
      std::string name(typeid(T).name());

      // Keep track of the length of the names.
      _maxNameLength = ( name.size() > _maxNameLength )
        ? name.size() : _maxNameLength;

      // Find the key in the map; if absent, create a new entry in the map.
      ListMap::iterator i = _mapOfLists.find(name);
      if ( i == _mapOfLists.end() ){
        std::pair<ListMap::iterator,bool> r =
          _mapOfLists.insert( std::make_pair(name,ListPtr(new ObjectList<T>())));
        i = r.first;
      }

      // Return a pointer to the list object of the right type.
      return dynamic_cast<ObjectList<T>*>(i->second.get());
    }

    // This object should not be copyable.
    AntiLeakRegistry( AntiLeakRegistry const&);
    AntiLeakRegistry& operator=(AntiLeakRegistry const&);

    typedef boost::shared_ptr<ObjectListBase> ListPtr;
    typedef std::map<std::string,ListPtr> ListMap;

    // The actual registry is this.
    ListMap _mapOfLists;

    // Keep track of string lengths for nice printout.
    size_t _maxNameLength;

  };

}

#endif /* G4Helper_AntiLeakRegistry_hh */
