// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SUBTRIANGULATION_IO_BINARY_HH
#define DUNE_SUBTRIANGULATION_IO_BINARY_HH

#include <iostream>

template<class T>
class Binary
{
private:
  T* ptr;
public:
  Binary () {}
  Binary (T & d) {
    ptr = & d;
  }
  Binary & operator = (T & d) {
    ptr = & d;
    return *this;
  }
  operator T ()
  {
    return *ptr;
  }
  operator const T () const
  {
    return *ptr;
  }
  char* data() {
    return reinterpret_cast<char*>(ptr);
  }
  const char* data() const {
    return reinterpret_cast<const char*>(ptr);
  }
};

template<class T>
class Binary<const T>
{
  // Binary<const T> is not allowed, to avoid problems with temporary objects
};

template <class T>
Binary<T> binary(T & t)
{
  static_assert(! std::is_const<T>::value,
                     "binary of temporary is not allowed");

  return Binary<T>(t);
}

template <class T>
Binary<T> binary(const T & t)
{
  static_assert(! std::is_const<const T>::value,
                     "binary of temporary is not allowed");

  return Binary<const T>(t);
}

template<class T>
std::istream & operator >> (std::istream & s, Binary<T> & d)
{
  s.read(d.data(), sizeof(T));
  return s;
}

template<class T>
std::ostream & operator << (std::ostream & s, const Binary<T> & d)
{
  s.write(d.data(), sizeof(T));
  return s;
}

#endif
