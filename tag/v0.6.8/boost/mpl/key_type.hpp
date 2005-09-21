
#ifndef BOOST_MPL_KEY_TYPE_HPP_INCLUDED
#define BOOST_MPL_KEY_TYPE_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2003-2004
// Copyright David Abrahams 2003-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Source: /cvsroot/boost/boost/boost/mpl/key_type.hpp,v $
// $Date: 2004/09/02 15:40:41 $
// $Revision: 1.2 $

#include <boost/mpl/key_type_fwd.hpp>
#include <boost/mpl/sequence_tag.hpp>
#include <boost/mpl/aux_/na_spec.hpp>
#include <boost/mpl/aux_/lambda_support.hpp>

namespace boost { namespace mpl {

template<
      typename BOOST_MPL_AUX_NA_PARAM(AssociativeSequence)
    , typename BOOST_MPL_AUX_NA_PARAM(T)
    >
struct key_type
    : key_type_impl< typename sequence_tag<AssociativeSequence>::type >
        ::template apply<AssociativeSequence,T>
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(2,key_type,(AssociativeSequence,T))
};

BOOST_MPL_AUX_NA_SPEC(2, key_type)

}}

#endif // BOOST_MPL_KEY_TYPE_HPP_INCLUDED
