/**
 *  $File: boost_pch.h $
 *  $LastChangedDate: 2015-08-29 11:10:32 -0500 (Sat, 29 Aug 2015) $
 *  $Rev: 4983 $
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _BOOST_PCH_H
#define _BOOST_PCH_H

#if TR1_SUPPORT == 0
#  include <map>
#elif TR1_SUPPORT == 1
#  include <unordered_map>
#else
#  include <tr1/unordered_map>
#endif


#include <boost/format.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>

#include "boost/tuple/tuple.hpp"
#include <boost/numeric/ublas/matrix.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/version.hpp>

#include "boost/lambda/lambda.hpp"

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "boost/lexical_cast.hpp"
#include "boost/dynamic_bitset/detail/lowest_bit.hpp"

#include "boost/regex.hpp"

#endif
