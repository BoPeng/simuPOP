/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu
*                                                                         *
*   $LastChangedDate: 2008-06-04 10:46:34 -0500 (Wed, 04 Jun 2008) $
*   $Rev: 1587 $                                                       *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along with this program; if not, write to the                         *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/

// This file is needed because swig generated wrap file needs
// stdint.h...

#include "boost/cstdint.hpp"

// boost defines these types but I can not use boost namespace in 
// simuPOP_cfg.h because 
// 1. other systems defines their own uint#_t types
// 2. SWIG does not understand boost::uint#_t and will raise type error
//    in the interfaces.
using boost::uint8_t;
using boost::uint16_t;
using boost::int32_t;
using boost::uint32_t;
