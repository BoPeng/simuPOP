/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate$
*   $Rev$
*
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

#include "terminator.h"

namespace simuPOP {

bool terminateIf::apply(population & pop)
{
    // experssion return true
    m_expr.setLocalDict(pop.dict());

    if (m_expr.valueAsBool() == true) {
        // store the generations this replicate stops
        int gen = pop.gen();                              // mainVars().getVarAsInt("gen");

        pop.setIntVar(m_var, gen);
        if (!this->noOutput() ) {
            ostream & out = this->getOstream(pop.dict());
            out << gen << endl;
            this->closeOstream();
        }
        if (!m_message.empty())
            cout << m_message << endl;
        if (m_stopAll)
            throw StopEvolution(m_message);
        return false;                                             // return false, this replicate will be stopped
    } else
        return true;
}


bool continueIf::apply(population & pop)
{
    // experssion return true
    m_expr.setLocalDict(pop.dict());

    if (m_expr.valueAsBool() == false) {
        // store the generations this replicate stops
        int gen = pop.gen();                              // mainVars().getVarAsInt("gen");

        pop.setIntVar(m_var, gen);
        if (!this->noOutput() ) {
            ostream & out = this->getOstream(pop.dict());
            out << gen << endl;
            this->closeOstream();
        }
        if (!m_message.empty())
            cout << m_message << endl;
        if (m_stopAll)
            throw StopEvolution(m_message);
        return false;                                             // return false, this replicate will be stopped
    } else
        return true;
}




}
