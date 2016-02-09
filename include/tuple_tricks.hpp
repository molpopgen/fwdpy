#ifndef FWDPY_TUPLE_TRICK_HPP
#define FWDPY_TUPLE_TRICK_HPP

#include <type_traits>
#include <tuple>

namespace fwdpy
{
  template<std::size_t I,typename ostream_t, typename tuple_t>
  struct ser_details
  {
    inline void operator()(ostream_t & o, const tuple_t & t) const
    {
      static_assert( std::is_pod<typename std::tuple_element<I-1,tuple_t>::type>::value,
		     "tuple element types must be POD" );
      auto x = std::get<I-1>(t);
      using type = typename std::remove_reference<decltype(x)>::type;
      o.write( reinterpret_cast<const char *>(&x),sizeof(type) );
      ser_details<I-1,ostream_t,tuple_t>()(o,t);
    }
  };

  template<typename ostream_t, typename tuple_t>
  struct ser_details<0,ostream_t,tuple_t>
  {
    inline void operator()(ostream_t & o, const tuple_t & t) const
    {
      static_assert( std::is_pod<typename std::tuple_element<0,tuple_t>::type>::value,
		     "tuple element types must be POD" );
      auto x = std::get<0>(t);
      using type = typename std::remove_reference<decltype(x)>::type;
      o.write( reinterpret_cast<const char *>(&x),sizeof(type) );
    }
  };

  template<typename ostream_t,typename tuple_t>
  void serialize_tuple_POD(ostream_t & o, const tuple_t & t)
  {
    ser_details<std::tuple_size<tuple_t>::value,ostream_t,tuple_t>()(o,t);
  }

  template<std::size_t I,typename istream_t, typename tuple_t>
  struct deser_details
  {
    inline void operator()(istream_t & i,  tuple_t & t) const
    {
      static_assert( std::is_pod<typename std::tuple_element<I-1,tuple_t>::type>::value,
		     "tuple element types must be POD" );
      using type = typename std::tuple_element<I-1,tuple_t>::type;
      type x;
      i.read( reinterpret_cast<char *>(&x),sizeof(type) );
      std::get<I-1>(t)=x;
      deser_details<I-1,istream_t,tuple_t>()(i,t);
    }
  };

  template<typename istream_t, typename tuple_t>
  struct deser_details<0,istream_t,tuple_t>
  {
    inline void operator()(istream_t & i, tuple_t & t) const
    {
      static_assert( std::is_pod<typename std::tuple_element<0,tuple_t>::type>::value,
		     "tuple element types must be POD" );
      using type = typename std::tuple_element<0,tuple_t>::type;
      type x;
      i.read( reinterpret_cast<char *>(&x),sizeof(type) );
      std::get<0>(t)=x;
    }
  };

  template<typename istream_t,typename tuple_t>
  void deserialize_tuple_POD(istream_t & i, tuple_t & t)
  {
    deser_details<std::tuple_size<tuple_t>::value,istream_t,tuple_t>()(i,t);
  }
}

#endif

