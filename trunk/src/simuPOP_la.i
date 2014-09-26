//
// $File: simuPOP_la.i $
// $LastChangedDate$
// $Rev$
//
// This file is part of simuPOP, a forward-time population genetics
// simulation environment. Please visit http://simupop.sourceforge.net
// for details.
//
// Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//

%module simuPOP_la

#define LONGALLELE

%include "simuPOP_common.i"
%pythoncode %{
defdict = _simuPOP_la.defdict
%}
