/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006)        *
*   $Rev: 191$                                                            *
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

#include "virtualSubPop.h"

namespace simuPOP {


duplicateSplitter::duplicateSplitter(vectori const & offWeights)
	: virtualSplitter(offWeights)
{
}


infoSplitter::infoSplitter(vectori const & offWeights)
	: virtualSplitter(offWeights)
{
}


proportionSplitter::proportionSplitter(vectori const & offWeights)
	: virtualSplitter(offWeights)
{
}


rangeSplitter::rangeSplitter(vectori const & offWeights)
	: virtualSplitter(offWeights)
{
}




}
