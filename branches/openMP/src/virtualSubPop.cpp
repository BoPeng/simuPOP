/**
 *  $File: virtualSubPop.cpp $
 *  $LastChangedDate: 2010-01-14 21:43:45 -0600 (Wed, 14 Jan 2010) $
 *  $Rev: 2335 $
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

#include "simuPOP_cfg.h"
#include "virtualSubPop.h"
#include "population.h"

#if PY_VERSION_HEX >= 0x03000000
#  define PyInt_AsLong(x) PyLong_AsLong(x)
#  define PyString_Check(x) PyUnicode_Check(x)
#endif

namespace simuPOP {

vspID::vspID(PyObject * obj) : m_subPop(InvalidSubPopID), m_virtualSubPop(InvalidSubPopID),
	m_spName(""), m_vspName(""), m_allAvailSP(false), m_allAvailVSP(false)
{
	if (PyNumber_Check(obj)) {
		// accept a number
		m_subPop = PyInt_AsLong(obj);
	} else if (PyString_Check(obj)) {
		m_spName = PyObj_AsString(obj);
#if PY_VERSION_HEX >= 0x03000000
	} else if (PyBytes_Check(obj)) {
		m_spName = PyBytes_AsString(obj);
#endif
	} else if (PySequence_Check(obj)) {
		int sz = PySequence_Size(obj);
		if (sz > 2) {
			cerr << "Invalid subpopulation ID." << endl;
			throw ValueError("Invalid subpopulation ID");
		}
		if (sz >= 1) {
			PyObject * sp = PySequence_GetItem(obj, 0);
			if (PyNumber_Check(sp))
				m_subPop = PyInt_AsLong(sp);
			else if (PyString_Check(sp)) {
				m_spName = PyObj_AsString(sp);
			}
#if PY_VERSION_HEX >= 0x03000000
			else if (PyBytes_Check(sp)) {
				m_spName = PyBytes_AsString(sp);
			}
#endif
			else {
				cerr << "Invalid vsp id for a subpopulation." << endl;
				throw ValueError("Invalid vsp id for a subpopulation");
			}
			Py_DECREF(sp);
		}
		if (sz == 2) {
			PyObject * vsp = PySequence_GetItem(obj, 1);
			if (PyNumber_Check(vsp))
				m_virtualSubPop = PyInt_AsLong(vsp);
			else if (PyString_Check(vsp)) {
				m_vspName = PyObj_AsString(vsp);
			}
#if PY_VERSION_HEX >= 0x03000000
			else if (PyBytes_Check(vsp)) {
				m_vspName = PyBytes_AsString(vsp);
			}
#endif
			else if (!PyBool_Check(vsp)) {
				cerr << "Invalid vsp id for a subpopulation." << endl;
				throw ValueError("Invalid vsp id for a subpopulation");
			}
			Py_DECREF(vsp);
		}
	} else {
		DBG_FAILIF(true, ValueError, "Invalid input for a (virtual) subpopulation.");
	}
}


vspID vspID::resolve(const Population & pop) const
{
	SubPopID sp = m_subPop;
	SubPopID vsp = m_virtualSubPop;

	if (!m_spName.empty())
		sp = pop.subPopByName(m_spName);
	if (!m_vspName.empty()) {
		DBG_ASSERT(pop.hasVirtualSubPop(), ValueError,
			"No virtual subpopulation is defined.");
		vsp = pop.virtualSplitter()->vspByName(m_vspName);
	}
	return vspID(sp, vsp);
}


subPopList::subPopList(PyObject * obj) : m_subPops(), m_allAvail(false)
{
	if (obj == NULL || obj == Py_None)
		// accept NULL
		m_allAvail = true;
	else if (PyBool_Check(obj)) {
		// accept True/False
		m_allAvail = obj == Py_True;
		return;
	} else if (PyNumber_Check(obj)) {
		// accept a number
		m_allAvail = false;
		m_subPops.push_back(vspID(PyInt_AsLong(obj)));
	} else if (PyString_Check(obj)) {
		m_allAvail = false;
		string name = PyObj_AsString(obj);
		m_subPops.push_back(vspID(InvalidSubPopID, InvalidSubPopID, false, false, name));
#if PY_VERSION_HEX >= 0x03000000
	} else if (PyBytes_Check(obj)) {
		m_allAvail = false;
		string name = PyBytes_AsString(obj);
		m_subPops.push_back(vspID(InvalidSubPopID, InvalidSubPopID, false, false, name));
#endif
	} else if (PySequence_Check(obj)) {
		m_subPops.resize(PySequence_Size(obj));
		// assign values
		for (size_t i = 0, iEnd = m_subPops.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			DBG_FAILIF(PyBool_Check(item), ValueError,
				"Invalid input for a (virtual) subpopulation in a subpopulation list.");
			if (PyNumber_Check(item))           // subpopulation
				m_subPops[i] = vspID(PyInt_AsLong(item));
			else if (PyString_Check(item)) {
				string name = PyObj_AsString(item);
				m_subPops[i] = vspID(InvalidSubPopID, InvalidSubPopID, false, false, name);
			}
#if PY_VERSION_HEX >= 0x03000000
			else if (PyBytes_Check(item)) {
				string name = PyBytes_AsString(item);
				m_subPops[i] = vspID(InvalidSubPopID, InvalidSubPopID, false, false, name);
			}
#endif
			else if (PySequence_Check(item)) {  // virtual subpopulation
				int sz = PySequence_Size(item);
				if (sz == 0) {
					m_subPops[i] = vspID();
					continue;
				}
				if (sz > 2) {
					cerr << "Invalid subpopulation ID." << endl;
					throw ValueError("Invalid subpopulation ID");
				}
				PyObject * sp = PySequence_GetItem(item, 0);
				SubPopID sp_id(InvalidSubPopID);
				string sp_name;
				if (PyNumber_Check(sp))
					sp_id = PyInt_AsLong(sp);
				else if (PyString_Check(sp)) {
					sp_name = PyObj_AsString(sp);
				}
#if PY_VERSION_HEX >= 0x03000000
				else if (PyBytes_Check(sp)) {
					sp_name = PyBytes_AsString(sp);
				}
#endif
				else if (!PyBool_Check(sp)) {
					cerr << "Invalid vsp id for a subpopulation." << endl;
					throw ValueError("Invalid vsp id for a subpopulation");
				}
				Py_DECREF(sp);
				if (sz == 1) {
					m_subPops[i] = vspID(sp_id, InvalidSubPopID, sp == Py_True, false, sp_name);
				} else {
					PyObject * vsp = PySequence_GetItem(item, 1);
					SubPopID vsp_id(InvalidSubPopID);
					string vsp_name;
					if (PyNumber_Check(vsp))
						vsp_id = PyInt_AsLong(vsp);
					else if (PyString_Check(vsp)) {
						vsp_name = PyObj_AsString(vsp);
					}
#if PY_VERSION_HEX >= 0x03000000
					else if (PyBytes_Check(vsp)) {
						vsp_name = PyBytes_AsString(vsp);
					}
#endif
					else if (!PyBool_Check(vsp)) {
						cerr << "Invalid vsp id for a subpopulation." << endl;
						throw ValueError("Invalid vsp id for a subpopulation");
					}
					Py_DECREF(vsp);
					m_subPops[i] = vspID(sp_id, vsp_id, sp == Py_True, vsp == Py_True, sp_name, vsp_name);
				}
			} else {
				cerr << "Invalid input for a list of (virtual) subpopulations." << endl;
				throw ValueError("Invalid input for a list of (virtual) subpopulations.");
			}
			Py_DECREF(item);
		}
	} else {
		cerr << "Invalid input for a list of (virtual) subpopulations." << endl;
		throw ValueError("Invalid input for a list of (virtual) subpopulations.");
	}
}


subPopList::subPopList(const vectorvsp & subPops)
	: m_subPops(subPops), m_allAvail(false)
{
	for (size_t i = 0; i < m_subPops.size(); ++i) {
		DBG_ASSERT(m_subPops[i].valid(), ValueError,
			"Invalid subpopulation ID");
	}
}


subPopList subPopList::expandFrom(const Population & pop) const
{
	DBG_FAILIF(m_allAvail && !m_subPops.empty(), SystemError,
		"Only when no subpopulation is specified can this function be called."
		"This is likely caused by the use of persistent subPops for different populations.");
	vectorvsp vsps;
	if (allAvail()) {
		for (size_t sp = 0; sp < pop.numSubPop(); ++sp)
			vsps.push_back(vspID(sp));
	} else {
		// otherwise, handle vsps such as (ALL_AVAIL, vsp)
		vectorvsp::const_iterator it = m_subPops.begin();
		vectorvsp::const_iterator it_end = m_subPops.end();
		for (; it != it_end; ++it) {
			if (it->allAvailSP()) {
				for (size_t sp = 0; sp < pop.numSubPop(); ++sp) {
					if (it->allAvailVSP()) {
						if (pop.numVirtualSubPop() == 0)
							vsps.push_back(vspID(sp));
						else
							for (size_t vsp = 0; vsp < pop.numVirtualSubPop(); ++vsp)
								vsps.push_back(vspID(sp, vsp));
					} else if (it->vspName().empty()) {
						// valid or invalid...
						vsps.push_back(vspID(sp, it->virtualSubPop()));
					} else {
						// with a name
						DBG_FAILIF(pop.numVirtualSubPop() == 0, ValueError,
							"No virtual subpopulation is defined.");
						vsps.push_back(vspID(sp, pop.virtualSplitter()->vspByName(it->vspName())));
					}
				}
			} else {
				SubPopID sp = InvalidSubPopID;
				if (it->spName().empty())
					sp = it->subPop();
				else
					sp = pop.subPopByName(it->spName());
				if (it->allAvailVSP()) {
					if (pop.numVirtualSubPop() == 0)
						vsps.push_back(vspID(sp));
					else
						for (size_t vsp = 0; vsp < pop.numVirtualSubPop(); ++vsp)
							vsps.push_back(vspID(sp, vsp));
				} else if (it->vspName().empty())
					vsps.push_back(vspID(sp, it->virtualSubPop()));
				else {
					// with a name
					DBG_FAILIF(pop.numVirtualSubPop() == 0, ValueError,
						"No virtual subpopulation is defined.");
					vsps.push_back(vspID(sp, pop.virtualSplitter()->vspByName(it->vspName())));
				}
			}
		}
	}
	return subPopList(vsps);
}


ostream & operator<<(ostream & out, const vspID & vsp)
{
	out << vsp.subPop();
	if (vsp.isVirtual())
		out << "," << vsp.virtualSubPop();
	return out;
}


ULONG BaseVspSplitter::countVisibleInds(const Population & pop, SubPopID subPop) const
{
	if (activatedSubPop() != subPop)
		return pop.subPopSize(subPop);
	ULONG count = 0;
	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it)
		if (it->visible())
			++count;
	return count;
}


UINT BaseVspSplitter::vspByName(const string & vspName) const
{
	if (!m_names.empty()) {
		vectorstr::const_iterator it = std::find(m_names.begin(), m_names.end(), vspName);
		DBG_FAILIF(it == m_names.end(), ValueError, "An invalid virtual subpopulation name is given: "
			+ vspName + ". Available names are: " + toStr(m_names));
		return it - m_names.begin();
	}
	for (size_t i = 0; i < numVirtualSubPop(); ++i) {
		if (vspName == name(i))
			return i;
	}
	string allNames;
	for (size_t i = 0; i < numVirtualSubPop(); ++i)
		allNames += name(i) + ", ";
	DBG_FAILIF(true, ValueError, "An invalid virtual subpopulation name is given: "
		+ vspName + ". Available names are: " + allNames);
	return 0;
}


CombinedSplitter::CombinedSplitter(const vectorsplitter & splitters,
	const intMatrix & vspMatrix, const stringList & names)
	: BaseVspSplitter(names), m_splitters(0), m_vspMap(0)
{
	for (size_t i = 0; i < splitters.size(); ++i)
		m_splitters.push_back(splitters[i]->clone());
	// default vsp map
	const matrixi & vspMap = vspMatrix.elems();
	if (vspMap.empty()) {
		size_t idx = 0;
		for (size_t i = 0; i < splitters.size(); ++i)
			for (size_t j = 0; j < splitters[i]->numVirtualSubPop(); ++j, ++idx)
				m_vspMap.push_back(vspList(1, vspPair(i, j)));
	} else {
		for (size_t i = 0; i < vspMap.size(); ++i) {
			vspList list;
			for (size_t j = 0; j < vspMap[i].size(); ++j) {
				// find out which splitter and which vsp
				int lower = 0;
				int higher = 0;
				bool done = false;
				for (size_t s = 0; s < splitters.size(); ++s) {
					higher += splitters[s]->numVirtualSubPop();
					if (vspMap[i][j] >= lower && vspMap[i][j] < higher) {
						list.push_back(vspPair(s, vspMap[i][j] - lower));
						done = true;
						break;
					}
					lower = higher;
				}
				DBG_ASSERT(done, IndexError,
					"Given vsp index " + toStr(vspMap[i][j]) + " is larger than the number of total VSPs");
			}
			m_vspMap.push_back(list);
		}
	}
}


CombinedSplitter::CombinedSplitter(const CombinedSplitter & rhs) :
	BaseVspSplitter(rhs), m_splitters(), m_vspMap(rhs.m_vspMap)
{
	for (size_t i = 0; i < rhs.m_splitters.size(); ++i)
		m_splitters.push_back(rhs.m_splitters[i]->clone());
}


CombinedSplitter::~CombinedSplitter()
{
	for (size_t i = 0; i < m_splitters.size(); ++i)
		delete m_splitters[i];
}


BaseVspSplitter * CombinedSplitter::clone() const
{
	return new CombinedSplitter(*this);
}


ULONG CombinedSplitter::size(const Population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_vspMap.size(), IndexError,
		"Virtual subpopulation index out of range");

	if (m_vspMap[virtualSubPop].empty())
		return 0;

	const vspList & list = m_vspMap[virtualSubPop];
	if (list.size() == 1)
		return m_splitters[list[0].first]->size(pop, subPop, list[0].second);

	size_t count = 0;
	for (size_t ind = 0; ind < pop.subPopSize(subPop); ++ind) {
		bool ok = false;
		for (size_t s = 0; s < list.size(); ++s) {
			if (m_splitters[list[s].first]->contains(pop, ind, vspID(subPop, list[s].second))) {
				ok = true;
				break;
			}
		}
		if (ok)
			++count;
	}
	return count;
}


bool CombinedSplitter::contains(const Population & pop, ULONG ind, vspID vsp) const
{
	SubPopID virtualSubPop = vsp.virtualSubPop();

	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_vspMap.size(), IndexError,
		"Virtual subpopulation index out of range");

	const vspList & list = m_vspMap[virtualSubPop];
	for (size_t s = 0; s < list.size(); ++s)
		if (m_splitters[list[s].first]->contains(pop, ind, vspID(vsp.subPop(), list[s].second)))
			return true;
	return false;
}


void CombinedSplitter::activate(const Population & pop, SubPopID subPop, SubPopID virtualSubPop)
{
	const vspList & list = m_vspMap[virtualSubPop];

	if (list.size() == 0)
		return;

	if (list.size() == 1)
		m_splitters[list[0].first]->activate(pop, subPop, list[0].second);
	else {
		const vspList & list = m_vspMap[virtualSubPop];
		for (size_t ind = 0; ind < pop.subPopSize(subPop); ++ind) {
			bool ok = false;
			for (size_t i = 0; i < list.size(); ++i) {
				if (m_splitters[list[i].first]->contains(pop, ind, vspID(subPop, list[i].second))) {
					ok = true;
					break;
				}
			}
			pop.individual(ind, subPop).setVisible(ok);
		}
	}
	m_activated = subPop;
}


string CombinedSplitter::name(SubPopID sp) const
{
	DBG_FAILIF(static_cast<UINT>(sp) >= numVirtualSubPop(), IndexError,
		"Virtual subpopulation index out of range");

	DBG_ASSERT(m_names.empty() || m_names.size() == numVirtualSubPop(), ValueError,
		"VSP names, if given, should be assigned to all VSPs");

	if (!m_names.empty())
		return m_names[sp];

	const vspList & list = m_vspMap[sp];
	string name;

	for (size_t i = 0; i < list.size(); ++i) {
		if (i != 0)
			name += " or ";
		name += m_splitters[list[i].first]->name(list[i].second);
	}
	return name;
}


ProductSplitter::ProductSplitter(const vectorsplitter & splitters, const stringList & names)
	: BaseVspSplitter(names), m_numVSP(0), m_subIndexes()
{
	for (size_t i = 0; i < splitters.size(); ++i) {
		if (m_numVSP == 0)
			m_numVSP = 1;
		m_numVSP *= splitters[i]->numVirtualSubPop();
		m_splitters.push_back(splitters[i]->clone());
	}
	for (size_t vsp = 0; vsp < m_numVSP; ++vsp) {
		vectori res(splitters.size());
		UINT tmpMod = m_numVSP;
		UINT tmpIdx = vsp;
		for (size_t i = 0; i < m_splitters.size(); ++i) {
			tmpMod /= m_splitters[i]->numVirtualSubPop();
			res[i] = tmpIdx / tmpMod;
			tmpIdx %= tmpMod;
		}
		m_subIndexes.push_back(res);
	}
}


ProductSplitter::ProductSplitter(const ProductSplitter & rhs) :
	BaseVspSplitter(rhs), m_splitters(), m_numVSP(rhs.m_numVSP), m_subIndexes(rhs.m_subIndexes)
{
	for (size_t i = 0; i < rhs.m_splitters.size(); ++i)
		m_splitters.push_back(rhs.m_splitters[i]->clone());
}


ProductSplitter::~ProductSplitter()
{
	for (size_t i = 0; i < m_splitters.size(); ++i)
		delete m_splitters[i];
}


BaseVspSplitter * ProductSplitter::clone() const
{
	return new ProductSplitter(*this);
}


ULONG ProductSplitter::size(const Population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	DBG_FAILIF(virtualSubPop >= m_numVSP, IndexError, "Subpopulation index out of range.");

	const vectori & idx = m_subIndexes[virtualSubPop];
	size_t count = 0;

	for (size_t i = 0; i < pop.subPopSize(subPop); ++i) {
		bool ok = true;
		for (size_t s = 0; s < m_splitters.size(); ++s) {
			if (!m_splitters[s]->contains(pop, i, vspID(subPop, idx[s]))) {
				ok = false;
				break;
			}
		}
		if (ok)
			++count;
	}
	return count;
}


bool ProductSplitter::contains(const Population & pop, ULONG ind, vspID vsp) const
{
	DBG_FAILIF(vsp.virtualSubPop() >= m_numVSP, IndexError, "Subpopulation index out of range.");
	const vectori & idx = m_subIndexes[vsp.virtualSubPop()];

	for (size_t s = 0; s < m_splitters.size(); ++s)
		if (!m_splitters[s]->contains(pop, ind, vspID(vsp.subPop(), idx[s])))
			return false;
	return true;
}


void ProductSplitter::activate(const Population & pop, SubPopID subPop, SubPopID virtualSubPop)
{
	DBG_FAILIF(virtualSubPop >= m_numVSP, IndexError, "Subpopulation index out of range.");
	const vectori & idx = m_subIndexes[virtualSubPop];

	for (size_t i = 0; i < pop.subPopSize(subPop); ++i) {
		bool ok = true;
		for (size_t s = 0; s < m_splitters.size(); ++s) {
			if (!m_splitters[s]->contains(pop, i, vspID(subPop, idx[s]))) {
				ok = false;
				break;
			}
		}
		pop.individual(i, subPop).setVisible(ok);
	}
	m_activated = subPop;
}


string ProductSplitter::name(SubPopID sp) const
{
	DBG_FAILIF(static_cast<UINT>(sp) >= numVirtualSubPop(), IndexError,
		"Virtual subpopulation index out of range");

	DBG_ASSERT(m_names.empty() || m_names.size() == numVirtualSubPop(), ValueError,
		"VSP names, if given, should be assigned to all VSPs");

	if (!m_names.empty())
		return m_names[sp];

	const vectori & idx = m_subIndexes[sp];
	string name;

	for (size_t i = 0; i < idx.size(); ++i) {
		if (i != 0)
			name += ", ";
		name += m_splitters[i]->name(idx[i]);
	}
	return name;
}


ULONG SexSplitter::size(const Population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	ULONG count = 0;
	Sex s = virtualSubPop == 0 ? MALE : FEMALE;

	for (; it != it_end; ++it)
		if (it->sex() == s)
			++count;
	return count;
}


bool SexSplitter::contains(const Population & pop, ULONG ind, vspID vsp) const
{
	return (vsp.virtualSubPop() == 0 ? MALE : FEMALE) == pop.individual(ind, vsp.subPop()).sex();
}


void SexSplitter::activate(const Population & pop, SubPopID subPop, SubPopID virtualSubPop)
{
	Sex s = virtualSubPop == 0 ? MALE : FEMALE;

	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);

	for (; it != it_end; ++it)
		it->setVisible(it->sex() == s);
	m_activated = subPop;
}


string SexSplitter::name(SubPopID vsp) const
{
	DBG_FAILIF(vsp > 1, IndexError, "Virtual subpopulation index out of range");

	DBG_ASSERT(m_names.empty() || m_names.size() == 2, ValueError,
		"VSP names, if given, should be assigned to all VSPs");

	if (!m_names.empty())
		return m_names[vsp];

	return vsp == 0 ? "Male" : "Female";
}


ULONG AffectionSplitter::size(const Population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	ULONG count = 0;
	// 0 is unaffected
	bool aff = virtualSubPop == 0 ? false : true;

	for (; it != it_end; ++it)
		if (it->affected() == aff)
			++count;
	return count;
}


bool AffectionSplitter::contains(const Population & pop, ULONG ind, vspID vsp) const
{
	return (vsp.virtualSubPop() == 0 ? false : true) == pop.individual(ind, vsp.subPop()).affected();
}


void AffectionSplitter::activate(const Population & pop, SubPopID subPop, SubPopID virtualSubPop)
{
	bool aff = virtualSubPop == 0 ? false : true;

	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);

	for (; it != it_end; ++it)
		it->setVisible(it->affected() == aff);
	m_activated = subPop;
}


string AffectionSplitter::name(SubPopID vsp) const
{
	DBG_FAILIF(vsp > 1, IndexError, "VSP index out of range");

	DBG_ASSERT(m_names.empty() || m_names.size() == 2, ValueError,
		"VSP names, if given, should be assigned to all VSPs");

	if (!m_names.empty())
		return m_names[vsp];

	return vsp == 0 ? "Unaffected" : "Affected";
}


InfoSplitter::InfoSplitter(string info, const floatList & values,
	const floatList & cutoff, const floatMatrix & ranges, const stringList & names)
	: BaseVspSplitter(names),
	m_info(info), m_values(values.elems()), m_cutoff(cutoff.elems()), m_ranges(ranges.elems())
{
	DBG_FAILIF(m_values.empty() && m_cutoff.empty() && m_ranges.empty(),
		ValueError, "Please specify either a list of values, a set of cutoff values or ranges");
	DBG_FAILIF(m_values.empty() + m_cutoff.empty() + m_ranges.empty() != 2,
		ValueError, "Please specify only one of parameters values, cutoff or ranges.");
	// cutoff value has to be ordered
	if (!m_cutoff.empty()) {
		for (size_t i = 1; i < m_cutoff.size(); ++i) {
			DBG_ASSERT(m_cutoff[i - 1] < m_cutoff[i],
				ValueError,
				"Cutoff values have to be in increasing order");
		}
	}
	if (!m_ranges.empty()) {
		for (size_t i = 1; i < m_ranges.size(); ++i) {
			DBG_FAILIF(m_ranges[i].size() != 2, ValueError,
				"Invalid information range.");
			DBG_ASSERT(m_ranges[i][1] > m_ranges[i][0],
				ValueError, "Invalid range.");
		}
	}
}


ULONG InfoSplitter::size(const Population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	UINT idx = pop.infoIdx(m_info);

	ULONG count = 0;

	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);

	if (!m_cutoff.empty()) {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) > m_cutoff.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_cutoff.size()));

		// using cutoff, below
		if (virtualSubPop == 0) {
			for (; it != it_end; ++it)
				if (it->info(idx) < m_cutoff[0])
					count++;
			return count;
		} else if (static_cast<UINT>(virtualSubPop) == m_cutoff.size()) {
			double v = m_cutoff.back();
			for (; it != it_end; ++it)
				if (it->info(idx) >= v)
					count++;
			return count;
		} else {         // in between
			double v1 = m_cutoff[virtualSubPop - 1];
			double v2 = m_cutoff[virtualSubPop];
			for (; it != it_end; ++it) {
				double v = it->info(idx);
				if (v >= v1 && v < v2)
					count++;
			}
			return count;
		}
	} else if (!m_values.empty()) {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_values.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_values.size() - 1));
		double v = m_values[virtualSubPop];
		for (; it != it_end; ++it)
			if (fcmp_eq(it->info(idx), v))
				count++;
		return count;
	} else {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_ranges.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_ranges.size() - 1));
		double v1 = m_ranges[virtualSubPop][0];
		double v2 = m_ranges[virtualSubPop][1];
		for (; it != it_end; ++it) {
			double v = it->info(idx);
			if (v >= v1 && v < v2)
				count++;
		}
		return count;
	}
	// should never reach here.
	return 0;
}


UINT InfoSplitter::numVirtualSubPop() const
{
	if (!m_cutoff.empty())
		return m_cutoff.size() + 1;
	else if (!m_values.empty())
		return m_values.size();
	else
		return m_ranges.size();
}


bool InfoSplitter::contains(const Population & pop, ULONG ind, vspID vsp) const
{
	SubPopID virtualSubPop = vsp.virtualSubPop();
	UINT idx = pop.infoIdx(m_info);

	if (!m_cutoff.empty()) {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) > m_cutoff.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_cutoff.size()));

		// using cutoff, below
		if (virtualSubPop == 0)
			return pop.individual(ind, vsp.subPop()).info(idx) < m_cutoff[0];
		else if (static_cast<UINT>(virtualSubPop) == m_cutoff.size())
			return pop.individual(ind, vsp.subPop()).info(idx) >= m_cutoff.back();
		else {         // in between
			double v1 = m_cutoff[virtualSubPop - 1];
			double v2 = m_cutoff[virtualSubPop];
			double v = pop.individual(ind, vsp.subPop()).info(idx);
			return v >= v1 && v < v2;
		}
	} else if (!m_values.empty()) {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_values.size(), IndexError,
			"Virtual Subpoplation index " + toStr(virtualSubPop) + " out of range of 0 ~ "
			+ toStr(m_values.size() - 1));
		return fcmp_eq(pop.individual(ind, vsp.subPop()).info(idx), m_values[virtualSubPop]);
	} else {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_ranges.size(), IndexError,
			"Virtual Subpoplation index " + toStr(virtualSubPop) + " out of range of 0 ~ "
			+ toStr(m_ranges.size() - 1));
		double v = pop.individual(ind, vsp.subPop()).info(idx);
		return v >= m_ranges[virtualSubPop][0] && v < m_ranges[virtualSubPop][1];
	}
	// this should not be reached.
	return false;
}


void InfoSplitter::activate(const Population & pop, SubPopID subPop, SubPopID virtualSubPop)
{
	UINT idx = pop.infoIdx(m_info);

	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);

	if (!m_cutoff.empty()) {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) > m_cutoff.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_cutoff.size()));

		// using cutoff, below
		if (virtualSubPop == 0) {
			for (; it != it_end; ++it)
				it->setVisible(it->info(idx) < m_cutoff[0]);
		} else if (static_cast<UINT>(virtualSubPop) == m_cutoff.size()) {
			double v = m_cutoff.back();
			for (; it != it_end; ++it)
				it->setVisible(it->info(idx) >= v);
		} else {         // in between
			double v1 = m_cutoff[virtualSubPop - 1];
			double v2 = m_cutoff[virtualSubPop];
			for (; it != it_end; ++it) {
				double v = it->info(idx);
				it->setVisible(v >= v1 && v < v2);
			}
		}
	} else if (!m_values.empty()) {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_values.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_values.size() - 1));
		double v = m_values[virtualSubPop];
		for (; it != it_end; ++it)
			it->setVisible(fcmp_eq(it->info(idx), v));
	} else {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_ranges.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_ranges.size() - 1));
		double v1 = m_ranges[virtualSubPop][0];
		double v2 = m_ranges[virtualSubPop][1];
		for (; it != it_end; ++it) {
			double v = it->info(idx);
			it->setVisible(v >= v1 && v < v2);
		}
	}
	m_activated = subPop;
}


string InfoSplitter::name(SubPopID sp) const
{
	DBG_FAILIF(static_cast<UINT>(sp) >= numVirtualSubPop(), IndexError,
		"Virtual subpopulation index out of range");

	DBG_ASSERT(m_names.empty() || m_names.size() == numVirtualSubPop(), ValueError,
		"VSP names, if given, should be assigned to all VSPs");

	if (!m_names.empty())
		return m_names[sp];

	if (!m_cutoff.empty()) {
		DBG_FAILIF(static_cast<UINT>(sp) > m_cutoff.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_cutoff.size()));
		if (sp == 0)
			return m_info + " < " + toStr(m_cutoff[0]);
		else if (static_cast<UINT>(sp) == m_cutoff.size())
			return m_info + " >= " + toStr(m_cutoff[sp - 1]);
		else
			return toStr(m_cutoff[sp - 1]) + " <= " + m_info + " < "
			       + toStr(m_cutoff[sp]);
	} else if (!m_values.empty()) {
		DBG_FAILIF(static_cast<UINT>(sp) >= m_values.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_values.size() - 1));
		return m_info + " = " + toStr(m_values[sp]);
	} else {
		DBG_FAILIF(static_cast<UINT>(sp) >= m_ranges.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_ranges.size() - 1));
		return toStr(m_ranges[sp][0]) + " <= " + m_info + " < "
		       + toStr(m_ranges[sp][1]);
	}
}


ProportionSplitter::ProportionSplitter(vectorf const & proportions, const stringList & names)
	: BaseVspSplitter(names), m_proportions(proportions)
{
	DBG_ASSERT(fcmp_eq(std::accumulate(proportions.begin(),
				proportions.end(), 0.), 1.), ValueError,
		"Passed proportions should add up to one");
}


ULONG ProportionSplitter::size(const Population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_proportions.size(), IndexError,
		"Virtual subpopulation index out of range");
	//
	if (static_cast<UINT>(virtualSubPop) < m_proportions.size() - 1)
		return static_cast<ULONG>(pop.subPopSize(subPop) * m_proportions[virtualSubPop]);
	// to avoid floating point problem, the last subpop is specially treated
	ULONG size = pop.subPopSize(subPop);
	ULONG spSize = size;
	// virtualSubPop == m_proportions.size() - 1
	for (int i = 0; i < virtualSubPop; ++i)
		spSize -= static_cast<ULONG>(size * m_proportions[i]);
	return spSize;
}


UINT ProportionSplitter::numVirtualSubPop() const
{
	return m_proportions.size();
}


bool ProportionSplitter::contains(const Population & pop, ULONG ind, vspID vsp) const
{
	DBG_FAILIF(static_cast<UINT>(vsp.virtualSubPop()) >= m_proportions.size(), IndexError,
		"Virtual subpopulation index out of range");

	ULONG size = pop.subPopSize(vsp.subPop());
	vectoru count(m_proportions.size());
	propToCount(m_proportions, size, count);

	ULONG lower = 0;
	ULONG higher = 0;
	for (UINT sp = 0; sp < m_proportions.size(); ++sp) {
		higher += count[sp];
		if (ind >= lower && ind < higher)
			return vsp.virtualSubPop() == static_cast<int>(sp);
		lower = higher;
	}
	return false;
}


void ProportionSplitter::activate(const Population & pop, SubPopID subPop, SubPopID virtualSubPop)
{
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_proportions.size(), IndexError,
		"Virtual subpopulation index out of range");

	ULONG size = pop.subPopSize(subPop);
	vectoru count(m_proportions.size());
	propToCount(m_proportions, size, count);
	// determine range
	ULONG lower = std::accumulate(count.begin(), count.begin() + virtualSubPop, 0L);
	ULONG higher = lower + count[virtualSubPop];

	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	for (size_t idx = 0; it != it_end; ++it, ++idx)
		it->setVisible(idx >= lower && idx < higher);
	m_activated = subPop;
}


string ProportionSplitter::name(SubPopID subPop) const
{
	DBG_FAILIF(static_cast<UINT>(subPop) >= numVirtualSubPop(), IndexError,
		"Virtual subpopulation index out of range");

	DBG_ASSERT(m_names.empty() || m_names.size() == numVirtualSubPop(), ValueError,
		"VSP names, if given, should be assigned to all VSPs");

	if (!m_names.empty())
		return m_names[subPop];

	return "Prop " + toStr(m_proportions[subPop]);
}


RangeSplitter::RangeSplitter(const intMatrix & ranges, const stringList & names)
	: BaseVspSplitter(names), m_ranges(ranges.elems())
{
	for (size_t i = 0; i < m_ranges.size(); ++i) {
		DBG_FAILIF(m_ranges[i].size() != 2
			|| m_ranges[i][0] > m_ranges[i][1],
			ValueError, "Wrong range [" +
			toStr(m_ranges[i][0]) + ", " +
			toStr(m_ranges[i][1]) + ")");
	}
}


ULONG RangeSplitter::size(const Population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_ranges.size(), IndexError,
		"Virtual subpopulation index out of range");
	if (static_cast<UINT>(m_ranges[virtualSubPop][0]) > pop.subPopSize(subPop))
		return 0;
	if (static_cast<UINT>(m_ranges[virtualSubPop][1]) > pop.subPopSize(subPop))
		return pop.subPopSize(subPop) - m_ranges[virtualSubPop][0];
	return m_ranges[virtualSubPop][1] - m_ranges[virtualSubPop][0];
}


UINT RangeSplitter::numVirtualSubPop() const
{
	return m_ranges.size();
}


bool RangeSplitter::contains(const Population & pop, ULONG ind, vspID vsp) const
{
	SubPopID virtualSubPop = vsp.virtualSubPop();

	return ind >= static_cast<ULONG>(m_ranges[virtualSubPop][0]) &&
	       ind < static_cast<ULONG>(m_ranges[virtualSubPop][1]);
}


void RangeSplitter::activate(const Population & pop, SubPopID subPop, SubPopID virtualSubPop)
{
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_ranges.size(), IndexError,
		"Virtual subpopulation index out of range");

	ULONG low = m_ranges[virtualSubPop][0];
	ULONG high = m_ranges[virtualSubPop][1];
	ULONG idx = 0;

	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it, ++idx)
		it->setVisible(idx >= low && idx < high);
	m_activated = subPop;
}


string RangeSplitter::name(SubPopID subPop) const
{
	DBG_FAILIF(static_cast<UINT>(subPop) >= numVirtualSubPop(), IndexError,
		"Virtual subpopulation index out of range");

	DBG_ASSERT(m_names.empty() || m_names.size() == numVirtualSubPop(), ValueError,
		"VSP names, if given, should be assigned to all VSPs");

	if (!m_names.empty())
		return m_names[subPop];

	return "Range [" + toStr(m_ranges[subPop][0]) + ", " +
	       toStr(m_ranges[subPop][1]) + ")";
}


GenotypeSplitter::GenotypeSplitter(const lociList & loci,
	const intMatrix & alleles, bool phase, const stringList & names)
	: BaseVspSplitter(names), m_loci(loci), m_alleles(alleles.elems()),
	m_phase(phase)
{
}


ULONG GenotypeSplitter::size(const Population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_alleles.size(), IndexError,
		"Virtual subpopulation index out of range");
	const vectori & alleles = m_alleles[virtualSubPop];
	ULONG count = 0;
	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it)
		if (match(&*it, alleles))
			++count;
	return count;
}


UINT GenotypeSplitter::numVirtualSubPop() const
{
	return m_alleles.size();
}


bool GenotypeSplitter::contains(const Population & pop, ULONG ind, vspID vsp) const
{
	SubPopID virtualSubPop = vsp.virtualSubPop();

	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_alleles.size(), IndexError,
		"Virtual subpopulation index out of genotype");

	const vectori & alleles = m_alleles[virtualSubPop];
	// this can be very slow if contains is used extensively.
	m_loci.locate(&pop);

	return match(&pop.individual(ind, vsp.subPop()), alleles);
}


void GenotypeSplitter::activate(const Population & pop, SubPopID subPop, SubPopID virtualSubPop)
{
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_alleles.size(), IndexError,
		"Virtual subpopulation index out of genotype");

	m_loci.locate(&pop);
	const vectori & alleles = m_alleles[virtualSubPop];
	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it)
		it->setVisible(match(&*it, alleles));
	m_activated = subPop;
}


//
// 1: 0 1
// 1: 0 1 1 1
// Two possibilities (depending on ploidy)
// 1, 2: 0 1 1 0 1 1 1 1
string GenotypeSplitter::name(SubPopID subPop) const
{
	DBG_FAILIF(static_cast<UINT>(subPop) >= numVirtualSubPop(), IndexError,
		"Virtual subpopulation index out of range");

	DBG_ASSERT(m_names.empty() || m_names.size() == numVirtualSubPop(), ValueError,
		"VSP names, if given, should be assigned to all VSPs");

	if (!m_names.empty())
		return m_names[subPop];

	string label = "Genotype ";
	if (m_loci.allAvail())
		label += "all loci";
	else if (m_loci.dynamic()) {
		for (size_t i = 0; i < m_loci.size(); ++i) {
			if (i != 0)
				label += ", ";
			label += m_loci.name(i);
		}
	} else {
		const vectoru & loci = m_loci.elems(NULL);
		for (size_t i = 0; i < m_loci.size(); ++i) {
			if (i != 0)
				label += ", ";
			label += toStr(loci[i]);
		}
	}
	label += ":";
	for (size_t i = 0; i < m_alleles[subPop].size(); ++i)
		label += " " + toStr(m_alleles[subPop][i]);
	return label;
}


bool GenotypeSplitter::match(const Individual * it, const vectori & alleles) const
{
	int ploidy = it->ploidy();
	int numLoci = m_loci.allAvail() ? m_loci.elems(NULL).size() : m_loci.size();

	unsigned choices = alleles.size() / ploidy / numLoci;

	DBG_FAILIF(alleles.size() != choices * ploidy * numLoci,
		ValueError, "Given genotype does not match population ploidy.");

	for (unsigned t = 0; t < choices; ++t) {
		vectori partial(alleles.begin() + t * ploidy * numLoci,
		                alleles.begin() + (t + 1) * ploidy * numLoci);
		if (matchSingle(it, partial))
			return true;
	}
	return false;
}


bool GenotypeSplitter::matchSingle(const Individual * it, const vectori & alleles) const
{
	int ploidy = it->ploidy();
	const vectoru & loci = m_loci.elems(NULL);

	if (m_phase || ploidy == 1) {
		// if phase=True, has to match exactly.
		UINT idx = 0;
		uintList::const_iterator loc = loci.begin();
		uintList::const_iterator loc_end = loci.end();
		for (; loc != loc_end; ++loc)
			for (int p = 0; p < ploidy; ++p)
				if (static_cast<int>(it->allele(*loc, p)) != alleles[idx++])
					return false;
		return true;
	} else if (ploidy == 2) {
		UINT idx = 0;
		uintList::const_iterator loc = loci.begin();
		uintList::const_iterator loc_end = loci.end();
		UINT numLoci = loci.size();
		for (; loc != loc_end; ++loc, ++idx) {
			int a1 = it->allele(*loc, 0);
			int a2 = it->allele(*loc, 1);
			int a3 = alleles[idx];
			int a4 = alleles[idx + numLoci];
			if (!((a1 == a3 && a2 == a4) || (a1 == a4 && a2 == a3)))
				return false;
		}
		return true;
	} else {
		UINT idx = 0;
		uintList::const_iterator loc = loci.begin();
		uintList::const_iterator loc_end = loci.end();
		UINT numLoci = loci.size();
		vectori v1(ploidy);
		vectori v2(ploidy);
		for (; loc != loc_end; ++loc, ++idx) {
			for (int p = 0; p < ploidy; ++p) {
				v1[p] = it->allele(*loc, p);
				v2[p] = alleles[idx + p * numLoci];
			}
			std::sort(v1.begin(), v1.end());
			std::sort(v2.begin(), v2.end());
			for (int p = 0; p < ploidy; ++p)
				if (v1[p] != v2[p])
					return false;
		}
		return true;
	}
	return false;
}


}     // namespce simuPOP

