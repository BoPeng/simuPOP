
/**
 *  $File: sandbox.cpp $
 *  $LastChangedDate: 2010-06-04 13:29:09 -0700 (Fri, 04 Jun 2010) $
 *  $Rev: 3579 $
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

#include "sandbox.h"

namespace simuPOP {


double InfSitesSelector::indFitness(Individual * ind, ULONG gen) const
{
    if (m_mode == MULTIPLICATIVE)
        return randomSelMulFitness(ind->genoBegin(), ind->genoEnd());
    else if (m_mode == ADDITIVE)
        return randomSelAddFitness(ind->genoBegin(), ind->genoEnd());
    else if (m_mode == EXPONENTIAL)
        return randomSelExpFitness(ind->genoBegin(), ind->genoEnd());
	return 0;
}

double InfSitesSelector::getFitnessValue(int mutant) const
{
    int sz = m_selDist.size();
    double s;
    if (sz == 0)
        // call a function
        s = m_selDist.func()(PyObj_As_Double, "()");
    else if (sz == 2)
        // constant
        s = m_selDist[1];
    else
        // a gamma distribution
        s = getRNG().randGamma(m_selDist[1], m_selDist[2]);
    m_selFactory[mutant] = s;
    return s;    
}

double InfSitesSelector::randomSelMulFitness(GenoIterator it, GenoIterator it_end) const 
{
    double s = 1;
    for (; it != it_end; ++it) {
        if (*it == 0)
            continue;
        intDict::iterator sit = m_selFactory.find(*it);
        if (sit == m_selFactory.end())
            s *= 1 + getFitnessValue(*it);
        else 
            s *= 1 + sit->second;
    }
    return s;
}

double InfSitesSelector::randomSelAddFitness(GenoIterator it, GenoIterator it_end) const
{
    double s = 0;
    for (; it != it_end; ++it) {
        if (*it == 0)
            continue;
        intDict::iterator sit = m_selFactory.find(*it);
        if (sit == m_selFactory.end())
            s += getFitnessValue(*it);
        else 
            s += sit->second;
    }
    return 1 + s > 0 ? 1 + s : 0;
}

double InfSitesSelector::randomSelExpFitness(GenoIterator it, GenoIterator it_end) const
{
    double s = 0;
    for (; it != it_end; ++it) {
        if (*it == 0)
            continue;
        intDict::iterator sit = m_selFactory.find(*it);
        if (sit == m_selFactory.end())
            s += getFitnessValue(*it);
        else 
            s += sit->second;
    }
    return exp(s);
}


}
