/**
 *  $File: utility.h $
 *  $LastChangedDate: 2012-03-27 13:52:12 -0700 (Tue, 27 Mar 2012) $
 *  $Rev: 4509 $
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

#ifndef _MUTANT_VECTOR_H
#define _MUTANT_VECTOR_H

#ifdef MUTANTALLELE

#  include <map>
#  include <iostream>

namespace simuPOP {

class vectorm
{
public:
	typedef Allele & reference;
	typedef const Allele & const_reference;
	typedef Allele * pointer;
	typedef const Allele * const_pointer;
	typedef std::map<size_t, Allele> storage;
	typedef storage::iterator val_iterator;
	typedef storage::const_iterator const_val_iterator;

	// Construction and destruction
	//
	// size:     the supposed real size, which is the boundary of indexes in m_data
	//
	inline vectorm (size_t size = 0) : m_size(size), m_data()
	{
	}


	inline vectorm (const vectorm & v) :
		m_size(v.m_size), m_data(v.m_data)
	{
	}


	// Accessors
	inline size_t size() const
	{
		return m_size;
	}


	inline const storage & data() const
	{
		return m_data;
	}


	inline storage & data()
	{
		return m_data;
	}


	void validate() const
	{
#  ifndef OPTIMIZED
		const_val_iterator it = m_data.begin();
		const_val_iterator it_end = m_data.end();
		for (; it != it_end; ++it) {
			DBG_ASSERT(it->second != 0, RuntimeError,
				(boost::format("Mutant with zero value is detected at location %1%") % it->first).str());
		}
#  endif
	}


public:
	inline void resize(size_t size, bool preserve = true)
	{
		m_size = size;
		if (preserve)
			m_data.erase(m_data.upper_bound(size), m_data.end());
		else
			m_data.clear();
	}


	// Zeroing, but do not set size to zero
	inline void clear()
	{
		m_data.clear();
	}


	inline void clear(size_t beg, size_t end)
	{
		m_data.erase(m_data.lower_bound(beg), m_data.lower_bound(end));
	}


	// Assignment
	inline vectorm & operator =(const vectorm & v)
	{
		if (this != &v) {
			m_size = v.m_size;
			m_data = v.m_data;
		}
		return *this;
	}


	// Swapping
	inline void swap(vectorm & v)
	{
		if (this != &v) {
			std::swap(m_size, v.m_size);
			m_data.swap(v.m_data);
		}
	}


	inline friend void swap(vectorm & v1, vectorm & v2)
	{
		v1.swap(v2);
	}


	// Back element insertion and erasure
	// This function does not change size_....
	inline void push_back(size_t i, const_reference t)
	{
		DBG_ASSERT(t != 0, RuntimeError, "Cannot store zero as mutant");
		m_data.insert(m_data.end(), storage::value_type(i, t));
	}


	class iterator;
	class const_iterator;

	// insert, added by Bo
	// NOTE: This function always insert at the back. Parameter dest is NOT used, it
	// is kept to make this function compatible to std::insert. Unlike push_back,
	// This function changes the size of vectorm.
	inline void insert(const iterator &, const const_iterator & ibeg, const const_iterator iend)
	{
		const_val_iterator beg = ibeg.get_val_iterator();
		const_val_iterator end = iend.get_val_iterator();
		ssize_t shift = m_size - ibeg.index();

		m_size += iend.index() - ibeg.index();
		const_val_iterator ptr = beg;
		for (; ptr != end; ++ptr) {
			DBG_ASSERT(ptr->second != 0, RuntimeError, "Cannot store zero as mutant");
			// we are inserting to the end, which should be constant instead of log(n) time
			m_data.insert(m_data.end(), storage::value_type(ptr->first + shift, ptr->second));
		}
	}


	// copy regions, added by Bo
	inline void copy_region(const const_iterator & begin, const const_iterator & end,
	                        iterator & it)
	{
		size_t iend = it.index() + (end - begin);
		ssize_t lagging = it.index() - begin.index();

		// remove old data
		if (!m_data.empty() && it.index() <= m_data.rbegin()->first)
			m_data.erase(m_data.lower_bound(it.index()),
				iend > m_size ? m_data.end() : m_data.lower_bound(iend));
		// insert new data
		const_val_iterator vbeg = begin.get_val_iterator();
		const_val_iterator vend = (end - (iend > m_size ? iend - m_size : 0)).get_val_iterator();
		if (vbeg == vend)
			return;
		// the first element is insert to get the right location for future insertion
		DBG_ASSERT(vbeg->second != 0, RuntimeError, "Cannot store zero as mutant");
		val_iterator dest = m_data.insert(m_data.end(), val_iterator::value_type(vbeg->first + lagging, vbeg->second));
		++vbeg;
		// according to the documentation, if ip in insert(ip, val) is set to the position
		// precedes the value to be inserted, this will make a very fast insertion.
		// because we are insertion from small to large numbers, we are setting ip to the
		// location of the last item
		for (; vbeg != vend; ++vbeg) {
			DBG_ASSERT(vbeg->second != 0, RuntimeError, "Cannot store zero as mutant");
			dest = m_data.insert(dest, val_iterator::value_type(vbeg->first + lagging, vbeg->second));
		}
#  if 0
		/*
		   The following code copies elements one by one, which can be more efficient if the
		   elements overlap a lot so that we do not have to remove and insert values at the
		   same location. In practice this rarely happens because we clear genotypes of the
		   offspring before mating.
		 */
		size_t iend = it.index() + (end - begin);
		ssize_t lagging = it.index() - begin.index();

		const_val_iterator sbeg = begin.get_val_iterator();
		const_val_iterator send = (end - (iend > m_size ? iend - m_size : 0)).get_val_iterator();
		if (sbeg == send)
			return;

		val_iterator dbeg = m_data.lower_bound(it.index());
		val_iterator dend = iend > m_size ? m_data.end() : m_data.lower_bound(iend);

		for (; sbeg != send; ++sbeg) {
			if (dbeg == dend) {
				// destination empty, no comparison is needed
				for ( ; sbeg != send; ++sbeg)
					dbeg = m_data.insert(dbeg, val_iterator::value_type(sbeg->first + lagging, sbeg->second));
				return;
			}
			if (sbeg->first < dbeg->first)          // insert
				m_data.insert(dbeg, val_iterator::value_type(sbeg->first + lagging, sbeg->second));
			else if (sbeg->first == dbeg->first) {  // assign
				dbeg->second = sbeg->second;
				++dbeg;
			} else  // remove
				m_data.erase(dbeg++);
		}
#  endif
	}


	//
	class iterator
	{
protected:
		vectorm * m_container;
		size_t m_index;

public:
		typedef std::input_iterator_tag iteratorcategory;
		typedef Allele & reference;
		typedef Allele * pointer;
		typedef long int difference_type;
		typedef Allele value_type;

		iterator (const iterator & iter) :
			m_container(iter.m_container), m_index(iter.m_index)
		{
		}


		iterator () : m_container(NULL), m_index(0)
		{
		}


		vectorm & operator()() const
		{
			return *m_container;
		}


		iterator (vectorm & v, size_t index)
			: m_container(&v), m_index(index)
		{
		}


		bool operator ==(const iterator & iter) const
		{
			return m_index == iter.m_index;
		}


		bool operator !=(const iterator & iter) const
		{
			return m_index != iter.m_index;
		}


		bool operator >=(const iterator & iter) const
		{
			return m_index >= iter.m_index;
		}


		bool operator <=(const iterator & iter) const
		{
			return m_index <= iter.m_index;
		}


		bool operator >(const iterator & iter) const
		{
			return m_index > iter.m_index;
		}


		bool operator <(const iterator & iter) const
		{
			return m_index < iter.m_index;
		}


		size_t index() const
		{
			return m_index;
		}


		size_t to_next() const
		{
			const_val_iterator it = (*this)().data().lower_bound(m_index + 1);

			return it == (*this)().data().end() ? (*this)().size() - m_index : it->first - m_index;
		}


		val_iterator get_val_iterator()
		{
			return (*this)().data().lower_bound(m_index);
		}


		const_val_iterator get_val_iterator() const
		{
			return (*this)().data().lower_bound(m_index);
		}


		const_reference value() const
		{
			const_val_iterator it((*this)().data().find(m_index));

			return (it == (*this)().data().end()) ? zero_ : it->second;
		}


		/// CPPONLY pre-incrment return by-reference
		iterator & operator++()
		{
			++m_index;
			return *this;
		}


		/// CPPONLY post-incrment return by-value
		iterator operator++(int)
		{
			iterator orig = *this;

			++m_index;
			return orig;
		}


		iterator & operator+=(const size_t size)
		{
			m_index += size;
			return *this;
		}


		iterator operator +(const size_t size) const
		{
			iterator result = *this;

			result.m_index += size;
			return result;
		}


		size_t operator -(const iterator & iter) const
		{
			return m_index - iter.m_index;
		}


		iterator operator -(const size_t size) const
		{
			iterator result = *this;

			result.m_index -= size;
			return result;
		}


		void assignIfDiffer(const_reference value)
		{
			// find the lower bound
			val_iterator it((*this)().data().lower_bound(m_index));

			// if the element does not exist
			if (it == (*this)().data().end() || it->first != m_index) {
				if (value != 0)
					// use lower_bound instead of find so that we can use the insert(it, val)
					// version of insert, which should be faster (constant vs log(n))
					(*this)().data().insert(it, storage::value_type(m_index, value));
				// if the element exists, but value is zero, remove it
			} else if (value == 0)
				(*this)().data().erase(it);
			// finally, update it directly
			else
				it->second = value;
		}


		friend class const_iterator;

		operator const_iterator() const
		{
			return *reinterpret_cast<const_iterator *>(const_cast<iterator *>(this));
		}
	};


	class const_iterator
	{
protected:
		const vectorm * m_container;
		size_t m_index;

public:
		typedef std::input_iterator_tag iteratorcategory;
		typedef Allele & reference;
		typedef Allele * pointer;
		typedef long int difference_type;
		typedef Allele value_type;

		const_iterator (const const_iterator & iter) :
			m_container(iter.m_container), m_index(iter.m_index)
		{
		}


		const_iterator () : m_container(NULL), m_index(0)
		{
		}


		const vectorm & operator()() const
		{
			return *m_container;
		}


		const_iterator (const vectorm & v, size_t index)
			: m_container(&v), m_index(index)
		{
		}


		bool operator ==(const const_iterator & iter) const
		{
			return m_index == iter.m_index;
		}


		bool operator !=(const const_iterator & iter) const
		{
			return m_index != iter.m_index;
		}


		bool operator >=(const const_iterator & iter) const
		{
			return m_index >= iter.m_index;
		}


		bool operator <=(const const_iterator & iter) const
		{
			return m_index <= iter.m_index;
		}


		bool operator >(const const_iterator & iter) const
		{
			return m_index > iter.m_index;
		}


		bool operator <(const const_iterator & iter) const
		{
			return m_index < iter.m_index;
		}


		size_t index() const
		{
			return m_index;
		}


		const_val_iterator get_val_iterator() const
		{
			return (*this)().data().lower_bound(m_index);
		}


		const_reference value() const
		{
			const_val_iterator it((*this)().data().find(m_index));

			return (it == (*this)().data().end()) ? zero_ : it->second;
		}


		/// CPPONLY pre-incrment return by-reference
		const_iterator & operator++()
		{
			++m_index;
			return *this;
		}


		/// CPPONLY post-incrment return by-value
		const_iterator operator++(int)
		{
			const_iterator orig = *this;

			++(*this);
			return orig;
		}


		const_iterator & operator+=(const size_t size)
		{
			m_index += size;
			return *this;
		}


		const_iterator operator +(const size_t size) const
		{
			const_iterator result = *this;

			result.m_index += size;
			return result;
		}


		size_t operator -(const const_iterator & iter) const
		{
			return m_index - iter.m_index;
		}


		const_iterator operator -(const size_t size) const
		{
			const_iterator result = *this;

			result.m_index -= size;
			return result;
		}


	};


	iterator begin()
	{
		return iterator(*this, 0);
	}


	const_iterator begin() const
	{
		return const_iterator(*this, 0);
	}


	iterator end()
	{
		return iterator(*this, size());
	}


	const_iterator end() const
	{
		return const_iterator(*this, size());
	}


private:
	size_t m_size;

	storage m_data;

	static const Allele zero_;

	friend class iterator;
	friend class const_iterator;
};


//typedef unsigned char Allele;
typedef Allele & AlleleRef;

}
#endif

#endif
