#ifndef MU2EG4_ANTILEAKREGISTRY_H
#define MU2EG4_ANTILEAKREGISTRY_H
//
// The Mu2e anti-leak system for G4.
//
// $Id: AntiLeakRegistry.hh,v 1.1 2010/11/11 23:19:41 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/11/11 23:19:41 $
//
// Original author Rob Kutschke
//
// This header defines an object registry that is used to avoid memory leaks in G4.
// The G4 examples often contain code like:
//
//    G4LogicalVolume v = ...;
//    v->SetVisAttributes(  new G4VisAttributes(true, color) );
//
// This works fine in the limited cases that G4 takes ownership of the new'ed object.
//
// In most cases, however, G4 does not take ownership of the new'ed object, which 
// creates a memory leak.  In Mu2e, the recommended practice for this situation is:
//
//      AntiLeakRegistry& reg(Mu2eG4RunManager::AntiLeakRegistry());
//      G4LogicalVolume v = ...;
//
//      G4VisAttributes* va = v->SetVisAttributes( reg.add(G4VisAttributes(true, color)) );
// or
//      v->SetVisAttributes( reg.addByPointer( new G4VisAttributes(true, color)) );
// or
//      G4VisAttributes* visAtt = new G4VisAttributes(true, color);
//      v->SetVisAttributes( reg.addByPointer( visAtt ) );
//
// where the first form is prefered provided that the new'ed object is safely copyable.
//
// The add method will make a copy of its argument, placing that copy in
// a std::list of the correct type and will return a pointer to the object in the list.
// The type std::list was chosen ( not, for example, std::vector ) because subsequent
// insertions into the list do not invalidate the pointer that was returned.
//
// The addByPtr method puts a shared pointer to the object into a list of the correct type.
// So there is no copy and the returned pointer points directly to the new'ed object.
//
// The first method has the disadvantage that it requires a copy. This has two sorts of
// potential problems. The big problem is if the type T is not copyable ( or should
// not be copyable because copying breaks something).  In that case use one of the
// other methods.  The lesser problem is that the copy is wasteful of time and memory; 
// however each add operation is usually done once at the start of the job,
// which makes the waste a minor consideration; with C++0X, which supports move aware objects, 
// we may be able to remove this objection.
//

#include <iostream>
#include <iomanip>
#include <list>
#include <map>
#include <string>
#include <typeinfo>

#include "boost/shared_ptr.hpp"

namespace mu2e
{

  // Base class that defines an interface obeyed by the concrete classes.
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

  class AntiLeakRegistry {
    
  public:
    AntiLeakRegistry():
      _mapOfLists(),
      _maxNameLength(0){}

    ~AntiLeakRegistry(){}

    template< class T>
    T* add( T& p ){

      // This will be the key in the map.
      std::string name(typeid(T).name());

      _maxNameLength = ( name.size() > _maxNameLength ) 
        ? name.size() : _maxNameLength;
      
      // Find the key in the map; if absent, create a new entry in the map.
      ListMap::iterator i = _mapOfLists.find(name);
      if ( i == _mapOfLists.end() ){
        ListPtr ol(new ObjectList<T>());
        std::pair<ListMap::iterator,bool> x = _mapOfLists.insert( std::make_pair(name,ol));
        i = x.first;
      }

      // Get a pointer to the list object of the right type.
      ObjectList<T>* olt = dynamic_cast<ObjectList<T>*>(i->second.get());

      // Add the object to be managed; this does a copy.
      olt->list.push_back(p);

      // Return the address of the object.
      return &olt->list.back();
    }

    template< class T>
    T* addByPointer( T* p ){

      // Object will be stored by shared pointer so that deletion is automatic.
      typedef boost::shared_ptr<T> SPtr;

      // The key in the map is the full name.
      std::string name(typeid(SPtr).name());

      _maxNameLength = ( name.size() > _maxNameLength ) 
        ? name.size() : _maxNameLength;

      // Find the key in the map; if absent, create a new entry in the map.
      ListMap::iterator i = _mapOfLists.find(name);
      if ( i == _mapOfLists.end() ){
        ListPtr ol(new ObjectList<SPtr>());
        std::pair<ListMap::iterator,bool> x = _mapOfLists.insert( std::make_pair(name,ol));
        i = x.first;
      }

      // Get a pointer to the list object of the right type.
      ObjectList<SPtr>* olt = dynamic_cast<ObjectList<SPtr>*>(i->second.get());

      // Add the object to be managed; this does a copy.
      olt->list.push_back(SPtr(p));

      // Return the address of the object.
      return olt->list.back().get();
    }

    // Clear everything; calls all destructors.
    void clear(){
      _mapOfLists.clear();
      _maxNameLength = 0;
    }

    // Print information about types and numbers of saved objects.
    void print( std::ostream& ost = std::cout ) const{

      ost << "\nMu2e Registry for G4 objects " << std::endl;
      ost << "Number of entries: " << _mapOfLists.size() << std::endl;

      for ( ListMap::const_iterator i = _mapOfLists.begin();
            i != _mapOfLists.end(); ++i ){
        ost << "  Data type: "
            << std::setw(_maxNameLength)
            << i->first << " "
            << i->second->size() << " Number of entries: "
            << std::endl;
      }

    }

    // How many types of objects have we stored?
    size_t size() const { return _mapOfLists.size(); }

  private:

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

#endif
