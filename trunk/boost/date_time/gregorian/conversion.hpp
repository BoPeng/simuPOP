#ifndef _GREGORIAN__CONVERSION_HPP___
#define _GREGORIAN__CONVERSION_HPP___

/* Copyright (c) 2004 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst
 * $Date: 2004/10/07 22:28:19 $
 */

#include <exception>
#include "boost/date_time/gregorian/gregorian_types.hpp"
#if defined(BOOST_DATE_TIME_INCLUDE_LIMITED_HEADERS)
#include "boost/date_time/gregorian/formatters_limited.hpp"
#else
#include "boost/date_time/gregorian/formatters.hpp"
#endif // BOOST_DATE_TIME_INCLUDE_LIMITED_HEADERS
#include "boost/date_time/c_time.hpp"

namespace boost {

namespace gregorian {


  //! Converts a date to a tm struct. Throws out_of_range exception if date is a special value
  inline
  tm to_tm(const date& d) 
  {
    if(d.is_pos_infinity() || d.is_neg_infinity() || d.is_not_a_date()){
      std::string s("tm unable to handle date value of " + to_simple_string(d));
      throw std::out_of_range(s);
    }
    tm datetm;
    boost::gregorian::date::ymd_type ymd = d.year_month_day();
    datetm.tm_year = ymd.year-1900; 
    datetm.tm_mon = ymd.month-1; 
    datetm.tm_mday = ymd.day;
    datetm.tm_wday = d.day_of_week();
    datetm.tm_yday = d.day_of_year()-1;
    datetm.tm_hour = datetm.tm_min = datetm.tm_sec = 0;
    datetm.tm_isdst = -1; // negative because not enough info to set tm_isdst
    return datetm;
  }

  //! Converts a tm structure into a date dropping the any time values.
  inline
  date date_from_tm(const tm& datetm) 
  {
    return date(datetm.tm_year+1900, datetm.tm_mon+1, datetm.tm_mday);
  }
  

} } //namespace boost::gregorian




#endif

