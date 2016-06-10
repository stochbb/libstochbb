#include "object.hh"
#include <list>
#include <iostream>
#include <typeinfo>
#include "logger.hh"

using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of GC
 * ********************************************************************************************* */
GC *GC::_instance = 0;

GC::GC()
  : _objects(), _lock()
{
  // pass...
}

GC::~GC() {
  // pass...
}

void
GC::run() {
  _lock.lock();
  // Mark all boxed objects:
  std::unordered_set<Object *>::iterator item=_objects.begin();
  for (; item!=_objects.end(); item++) {
    if ((*item)->refcount()) { (*item)->mark(); }
  }
  // Delete all unmarked objects:
  for (item=_objects.begin(); item!=_objects.end(); ) {
    if ((*item)->isMarked()) {
      (*item)->unmark(); item++;
    } else {
      delete *item;
      item = _objects.erase(item);
    }
  }
  _lock.unlock();
}

void
GC::add(Object *obj) {
  _lock.lock();
  _objects.insert(obj);
  _lock.unlock();
}

GC&
GC::get() {
  if (0 == GC::_instance) {
    GC::_instance = new GC();
  }
  return *(GC::_instance);
}


/* ********************************************************************************************* *
 * Implementation of Object
 * ********************************************************************************************* */
Object::Object()
  : _marked(false), _refcount(1)
{
  GC::get().add(this);
}

Object::~Object() {
  // pass...
}

void
Object::unref() {
  if (! _refcount) {
    throw std::exception();
  }
  _refcount--;
  if (!_refcount)
    GC::get().run();
}

void
Object::mark() {
  if (_marked) { return; }
  _marked = true;
}

void
Object::unmark() {
  _marked = false;
}

void
Object::print(std::ostream &stream) const {
  stream << "<Object #" << (void *) this << ">";
}
