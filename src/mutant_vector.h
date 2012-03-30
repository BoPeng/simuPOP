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

#  include "boost/numeric/ublas/vector_sparse.hpp"

using boost::numeric::ublas::sparse_bidirectional_iterator_tag;
using boost::numeric::ublas::reverse_iterator_base;
using boost::numeric::ublas::vector_container;
using boost::numeric::ublas::vector_assign;
using boost::numeric::ublas::vector_assign_scalar;
using boost::numeric::ublas::scalar_minus_assign;
using boost::numeric::ublas::scalar_plus_assign;
using boost::numeric::ublas::scalar_multiplies_assign;
using boost::numeric::ublas::scalar_divides_assign;
using boost::numeric::ublas::container_const_reference;
using boost::numeric::ublas::container_reference;
using boost::numeric::ublas::vector_expression;
using boost::numeric::ublas::scalar_assign;
using boost::numeric::ublas::bidirectional_iterator_base;
using boost::numeric::ublas::sparse_tag;
using boost::numeric::ublas::vector_reference;
using boost::numeric::ublas::unbounded_array;

namespace simuPOP {

//typedef unbounded_array<std::size_t> IndexArray;
//typedef unbounded_array<Allele> ValueArray;
typedef std::vector<std::size_t> IndexArray;
typedef std::vector<Allele> ValueArray;

template<class I, class C>
inline I _lower_bound(const I & begin, const I & end, const size_t & t, C compare)
{
	//
	// if empty              ==> return begin
	// if t <= first element ==> return begin
	// if t > last element   ==> return end
	// return the first element that insertion will not change order
	//
	// t <= *begin <=> ! (*begin < t)
	if (begin == end || !compare(*begin, t))
		return begin;
	if (compare(*(end - 1), t))
		return end;
	return std::lower_bound(begin, end, t, compare);
}


class vectorm :
	public vector_container<vectorm>
{
	typedef Allele & true_reference;
	typedef Allele * pointer;
	typedef const Allele * const_pointer;
	typedef vectorm self_type;

public:
#  ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
	using vector_container<self_type>::operator ();
#  endif
	// ISSUE require type consistency check
	// is_convertable (IndexArray::size_type, ValueArray::size_type)
	//typedef IndexArray::value_type size_type;
	typedef size_t size_type;
	typedef IndexArray::difference_type difference_type;
	typedef Allele value_type;
	typedef const Allele & const_reference;
	typedef Allele & reference;
	typedef IndexArray index_array_type;
	typedef ValueArray value_array_type;
	typedef const vector_reference<const self_type> const_closure_type;
	typedef vector_reference<self_type> closure_type;
	typedef self_type vector_temporary_type;
	typedef sparse_tag storage_category;

	// Construction and destruction
	//
	// size:     the supposed real size, which is the boundary of indexes in index_data_
	// filled:   number of elements filled, that are guaranteed to be in order in index_data_
	// capacity: the size of the array holding data.
	//
	// restrict_capacity: greater than 1, less than (or equal to size_)
	//
	inline vectorm () :
		vector_container<self_type> (),
		size_(0), capacity_(restrict_capacity(0)), filled_(0),
		index_data_(capacity_), value_data_(capacity_)
	{
		storage_invariants();
	}


	explicit inline vectorm (size_type size, size_type non_zeros = 0) :
		vector_container<self_type> (),
		size_(size), capacity_(restrict_capacity(non_zeros)), filled_(0),
		index_data_(capacity_), value_data_(capacity_)
	{
		storage_invariants();
	}


	inline vectorm (const vectorm & v) :
		vector_container<self_type> (),
		size_(v.size_), capacity_(v.capacity_), filled_(v.filled_),
		index_data_(v.index_data_), value_data_(v.value_data_)
	{
		storage_invariants();
	}


	// index and value are swapped in for best performance
	// in this function filled_ == capacity_ so there is no room to grow
	inline vectorm (size_type size, IndexArray & index, ValueArray & value) :
		vector_container<self_type> (),
		size_(size), capacity_(restrict_capacity(index.size())), filled_(index.size()),
		index_data_(), value_data_()
	{
		index_data_.swap(index);
		value_data_.swap(value);
		storage_invariants();
	}


	template<class AE>
	inline vectorm (const vector_expression<AE> & ae, size_type non_zeros = 0) :
		vector_container<self_type> (),
		size_(ae().size()), capacity_(restrict_capacity(non_zeros)), filled_(0),
		index_data_(capacity_), value_data_(capacity_)
	{
		storage_invariants();
		vector_assign<scalar_assign> (*this, ae);
	}


	// Accessors
	inline size_type size() const
	{
		return size_;
	}


	inline size_type nnz_capacity() const
	{
		return capacity_;
	}


	inline size_type nnz() const
	{
		return filled_;
	}


	inline index_array_type::size_type filled() const
	{
		return filled_;
	}


	inline const index_array_type & index_data() const
	{
		return index_data_;
	}


	inline const value_array_type & value_data() const
	{
		return value_data_;
	}


	inline void set_filled(const index_array_type::size_type & filled)
	{
		filled_ = filled;
		storage_invariants();
	}


	inline index_array_type & index_data()
	{
		return index_data_;
	}


	inline value_array_type & value_data()
	{
		return value_data_;
	}


	// Resizing

private:
	inline size_type restrict_capacity(size_type non_zeros) const
	{
		non_zeros = (std::max)(non_zeros, size_type(1));
		non_zeros = (std::min)(non_zeros, size_);
		return non_zeros;
	}


public:
	inline void resize(size_type size, bool preserve = true)
	{
		size_ = size;
		capacity_ = restrict_capacity(capacity_);
		if (preserve) {
			index_data_.resize(capacity_, size_type());
			value_data_.resize(capacity_, value_type());
			filled_ = (std::min)(capacity_, filled_);
			while ((filled_ > 0) && (index_data_[filled_ - 1] >= size)) {
				--filled_;
			}
		} else {
			index_data_.resize(capacity_);
			value_data_.resize(capacity_);
			filled_ = 0;
		}
		storage_invariants();
	}


	// Reserving
	inline void reserve(size_type non_zeros, bool preserve = true)
	{
		capacity_ = restrict_capacity(non_zeros);
		if (preserve) {
			index_data_.resize(capacity_, size_type());
			value_data_.resize(capacity_, value_type());
			filled_ = (std::min)(capacity_, filled_);
		}else {
			index_data_.resize(capacity_);
			value_data_.resize(capacity_);
			filled_ = 0;
		}
		storage_invariants();
	}


	// Element support
	inline pointer find_element(size_type i)
	{
		return const_cast<pointer> (const_cast<const self_type &>(*this).find_element(i));
	}


	inline const_pointer find_element(size_type i) const
	{
		const_subiterator_type it(_lower_bound(index_data_.begin(), index_data_.begin() + filled_, i, std::less<size_type> ()));

		if (it == index_data_.begin() + filled_ || *it != i)
			return 0;
		return &value_data_ [it - index_data_.begin()];
	}


	// Element access
	inline const_reference operator ()(size_type i) const
	{
		BOOST_UBLAS_CHECK(i < size_, bad_index());
		const_subiterator_type it(_lower_bound(index_data_.begin(), index_data_.begin() + filled_, i, std::less<size_type> ()));
		if (it == index_data_.begin() + filled_ || *it != i)
			return zero_;
		return value_data_ [it - index_data_.begin()];
	}


	inline true_reference ref(size_type i)
	{
		BOOST_UBLAS_CHECK(i < size_, bad_index());
		subiterator_type it(_lower_bound(index_data_.begin(), index_data_.begin() + filled_, i, std::less<size_type> ()));
		if (it == index_data_.begin() + filled_ || *it != i)
			return insert_element(i, value_type/*zero*/ ());
		else
			return value_data_ [it - index_data_.begin()];
	}


	inline reference operator ()(size_type i)
	{
		return ref(i) ;
	}


	inline const_reference operator [](size_type i) const
	{
		return (*this)(i);
	}


	inline reference operator [](size_type i)
	{
		return (*this)(i);
	}


	// Element assignment
	inline true_reference insert_element(size_type i, const_reference t)
	{
		BOOST_UBLAS_CHECK(!find_element(i), bad_index());               // duplicate element
		if (filled_ >= capacity_)
			reserve(2 * capacity_, true);
		subiterator_type it(_lower_bound(index_data_.begin(), index_data_.begin() + filled_, i, std::less<size_type> ()));
		// ISSUE max_capacity limit due to difference_type
		std::iterator_traits<subiterator_type>::difference_type n = it - index_data_.begin();
		BOOST_UBLAS_CHECK(filled_ == 0 || filled_ == index_array_type::size_type(n) || *it != i, internal_logic());           // duplicate found by _lower_bound
		++filled_;
		it = index_data_.begin() + n;
		std::copy_backward(it, index_data_.begin() + filled_ - 1, index_data_.begin() + filled_);
		*it = i;
		value_array_type::iterator itt(value_data_.begin() + n);
		std::copy_backward(itt, value_data_.begin() + filled_ - 1, value_data_.begin() + filled_);
		*itt = t;
		storage_invariants();
		return *itt;
	}


	inline void erase_element(size_type i)
	{
		subiterator_type it(_lower_bound(index_data_.begin(), index_data_.begin() + filled_, i, std::less<size_type> ()));

		std::iterator_traits<subiterator_type>::difference_type n = it - index_data_.begin();
		if (filled_ > index_array_type::size_type(n) && *it == i) {
			std::copy(it + 1, index_data_.begin() + filled_, it);
			value_array_type::iterator itt(value_data_.begin() + n);
			std::copy(itt + 1, value_data_.begin() + filled_, itt);
			--filled_;
		}
		storage_invariants();
	}


	// Zeroing
	inline void clear()
	{
		filled_ = 0;
		storage_invariants();
	}


	// Assignment
	inline vectorm & operator =(const vectorm & v)
	{
		if (this != &v) {
			size_ = v.size_;
			capacity_ = v.capacity_;
			filled_ = v.filled_;
			index_data_ = v.index_data_;
			value_data_ = v.value_data_;
		}
		storage_invariants();
		return *this;
	}


	// Swapping
	inline void swap(vectorm & v)
	{
		if (this != &v) {
			std::swap(size_, v.size_);
			std::swap(capacity_, v.capacity_);
			std::swap(filled_, v.filled_);
			index_data_.swap(v.index_data_);
			value_data_.swap(v.value_data_);
		}
		storage_invariants();
	}


	inline friend void swap(vectorm & v1, vectorm & v2)
	{
		v1.swap(v2);
	}


	// Back element insertion and erasure
	inline void push_back(size_type i, const_reference t)
	{
		BOOST_UBLAS_CHECK(filled_ == 0 || index_data_ [filled_ - 1] < i, external_logic());
		if (filled_ >= capacity_)
			reserve(2 * capacity_, true);
		BOOST_UBLAS_CHECK(filled_ < capacity_, internal_logic());
		index_data_ [filled_] = i;
		value_data_ [filled_] = t;
		++filled_;
		storage_invariants();
	}

	// clear some ... added by Bo
	inline void clear(size_t beg, size_t end)
	{
		if (beg == end)
			return;
		BOOST_UBLAS_CHECK(beg < end && end <= size_, external_logic());
		size_t b = _lower_bound(index_data_.begin(), index_data_.begin() + filled_, beg, std::less<size_type> ()) - index_data_.begin();
		size_t e = _lower_bound(index_data_.begin(), index_data_.begin() + filled_, end, std::less<size_type> ()) - index_data_.begin();
		for (size_t i = b; i < e; ++i)
			value_data_[i] = 0;
	}

	class const_val_iterator;

	// insert, added by Bo
	inline void insert_back(const const_val_iterator & beg, const const_val_iterator end, int shift = 0)
	{
		BOOST_UBLAS_CHECK(filled_ == 0 || index_data_ [filled_ - 1] < i, external_logic());
		if (filled_ + (end - beg) >= capacity_)
			reserve(2 * std::max(capacity_, filled_ + (end - beg)), true);
		BOOST_UBLAS_CHECK(filled_ < capacity_, internal_logic());
		const_val_iterator ptr = beg;
		BOOST_UBLAS_CHECK(beg == end || beg != (*this)().end(), bad_index());
		for (; ptr != end; ++ptr) {
			if (*ptr != 0) {
				index_data_ [filled_] = ptr.index() + shift;
				value_data_ [filled_] = *ptr;
				++filled_;
			}
		}
		storage_invariants();
	}


	inline void pop_back()
	{
		BOOST_UBLAS_CHECK(filled_ > 0, external_logic());
		--filled_;
		storage_invariants();
	}


	// Iterator types

private:
	// Use index array iterator
	typedef IndexArray::const_iterator const_subiterator_type;
	typedef IndexArray::iterator subiterator_type;

	inline true_reference at_element(size_type i)
	{
		BOOST_UBLAS_CHECK(i < size_, bad_index());
		subiterator_type it(_lower_bound(index_data_.begin(), index_data_.begin() + filled_, i, std::less<size_type> ()));
		BOOST_UBLAS_CHECK(it != index_data_.begin() + filled_ && *it == i, bad_index());
		return value_data_ [it - index_data_.begin()];
	}


public:
	class const_val_iterator;
	class val_iterator;

	// Element lookup
	// inline  This function seems to be big. So we do not let the compiler inline  it.
	const_val_iterator find(size_type i) const
	{
		return const_val_iterator(*this, _lower_bound(index_data_.begin(), index_data_.begin() + filled_, i, std::less<size_type> ()));
	}


	// inline  This function seems to be big. So we do not let the compiler inline  it.
	val_iterator find(size_type i)
	{
		return val_iterator(*this, _lower_bound(index_data_.begin(), index_data_.begin() + filled_, i, std::less<size_type> ()));
	}


	class const_val_iterator :
		public container_const_reference<vectorm>,
		public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
		                                   const_val_iterator, value_type>
	{
public:
		typedef vectorm::value_type value_type;
		typedef vectorm::difference_type difference_type;
		typedef vectorm::const_reference reference;
		typedef const vectorm::pointer pointer;

		// Construction and destruction
		inline const_val_iterator () :
			container_const_reference<self_type> (), it_() {}
		inline const_val_iterator (const self_type & v, const const_subiterator_type & it) :
			container_const_reference<self_type> (v), it_(it) {}
		inline const_val_iterator (const self_type::val_iterator & it) :          // ISSUE self_type:: stops VC8 using std::iterator here
			container_const_reference<self_type> (it()), it_(it.it_) {}

		// Arithmetic
		inline const_val_iterator & operator ++()
		{
			++it_;
			return *this;
		}


		inline const_val_iterator & operator --()
		{
			--it_;
			return *this;
		}


		// Dereference
		inline const_reference operator *() const
		{
			BOOST_UBLAS_CHECK(index() < (*this)().size(), bad_index());
			return (*this)().value_data_ [it_ - (*this)().index_data_.begin()];
		}


		// Index
		inline size_type index() const
		{
			BOOST_UBLAS_CHECK(*this != (*this)().end(), bad_index());
			BOOST_UBLAS_CHECK(*it_ < (*this)().size(), bad_index());
			return *it_;
		}


		// added by Bo to calcualte number of elements between things
		inline size_type operator -(const const_val_iterator & it) const
		{
			return it_ - it.it_;
		}


		// Assignment
		inline const_val_iterator & operator =(const const_val_iterator & it)
		{
			container_const_reference<self_type>::assign(&it());
			it_ = it.it_;
			return *this;
		}


		// Comparison
		inline bool operator ==(const const_val_iterator & it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			return it_ == it.it_;
		}


private:
		const_subiterator_type it_;
	};

	inline const_val_iterator val_begin() const
	{
		return find(0);
	}


	inline const_val_iterator val_begin(size_type idx) const
	{
		return find(idx);
	}


	inline const_val_iterator val_end() const
	{
		return find(size_);
	}


	class val_iterator :
		public container_reference<vectorm>,
		public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
		                                   val_iterator, value_type>
	{
public:
		typedef vectorm::value_type value_type;
		typedef vectorm::difference_type difference_type;
		typedef vectorm::true_reference reference;
		typedef vectorm::pointer pointer;

		// Construction and destruction
		inline val_iterator () :
			container_reference<self_type> (), it_() {}
		inline val_iterator (self_type & v, const subiterator_type & it) :
			container_reference<self_type> (v), it_(it) {}

		// Arithmetic
		inline val_iterator & operator ++()
		{
			++it_;
			return *this;
		}


		inline val_iterator & operator --()
		{
			--it_;
			return *this;
		}


		// Dereference
		inline reference operator *() const
		{
			BOOST_UBLAS_CHECK(index() < (*this)().size(), bad_index());
			return (*this)().value_data_ [it_ - (*this)().index_data_.begin()];
		}


		// Index
		inline size_type index() const
		{
			BOOST_UBLAS_CHECK(*this != (*this)().end(), bad_index());
			BOOST_UBLAS_CHECK(*it_ < (*this)().size(), bad_index());
			return *it_;
		}


		// Assignment
		inline val_iterator & operator =(const val_iterator & it)
		{
			container_reference<self_type>::assign(&it());
			it_ = it.it_;
			return *this;
		}


		// Comparison
		inline bool operator ==(const val_iterator & it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			return it_ == it.it_;
		}


private:
		subiterator_type it_;

		friend class const_val_iterator;
	};

	inline val_iterator val_begin(size_type idx)
	{
		return find(idx);
	}


	inline val_iterator val_begin()
	{
		return find(0);
	}


	inline val_iterator val_end()
	{
		return find(size_);
	}


	// Reverse iterator
	typedef reverse_iterator_base<const_val_iterator> const_reverse_val_iterator;
	typedef reverse_iterator_base<val_iterator> reverse_val_iterator;

	inline const_reverse_val_iterator val_rbegin() const
	{
		return const_reverse_val_iterator(val_end());
	}


	inline const_reverse_val_iterator val_rend() const
	{
		return const_reverse_val_iterator(val_begin());
	}


	inline reverse_val_iterator val_rbegin()
	{
		return reverse_val_iterator(val_end());
	}


	inline reverse_val_iterator val_rend()
	{
		return reverse_val_iterator(val_begin());
	}


	class const_iterator;

	class iterator :
		public container_reference<self_type>
	{
protected:
		mutable size_t m_index;
		mutable ssize_t m_com_index;

public:
		typedef std::input_iterator_tag iteratorcategory;
		typedef Allele & reference;
		typedef Allele * pointer;
		typedef long int difference_type;
		typedef Allele value_type;

		iterator (const iterator & iter) : container_reference<self_type>(const_cast<self_type &>(iter())),
			m_index(iter.m_index), m_com_index(iter.m_com_index)
		{
		}


		iterator () : container_reference<self_type>(), m_index(0), m_com_index(-1)
		{
		}


		iterator (self_type & v) : container_reference<self_type>(const_cast<self_type &>(v)), m_index(0), m_com_index(-1)
		{
			// m_com_index is the smallest index of mutants (or -1 is nothing is there)
		}


		iterator (const self_type & v, size_t index)
			: container_reference<self_type>(const_cast<self_type &>(v)), m_index(index)
		{
			if (index == 0) {
				if ((*this)().index_data().size() > 0)
					if ((*this)().index_data()[0] > 0)
						m_com_index = -1;
					else
						m_com_index = 0;
				else
					m_com_index = -1;
			} else {
				index_array_type::iterator lower =
				    std::lower_bound((*this)().index_data().begin(),
						(*this)().index_data().begin() + (*this)().filled(), index);
				m_com_index = lower - (*this)().index_data().begin();
			}
			// find the smallest m_com_index  that is larger than index
		}


		iterator & operator =(const iterator & iter)
		{
			container_reference<self_type>::assign(&iter());
			m_index = iter.m_index;
			m_com_index = iter.m_com_index;
			return *this;
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

		val_iterator getValIterator() const
		{
			return (*this)().val_begin(m_index);
		}

		const_reference value() const
		{
			static const Allele zero = 0;

			if (m_com_index < (int)(*this)().index_data().size()) {
				if (m_com_index == -1)
					return zero;
				else if (m_index < (*this)().index_data()[m_com_index])
					return zero;
				else
					return (*this)().value_data()[m_com_index];
			} else
				return zero;
		}


		const_reference operator *() const
		{
			static const Allele zero = 0;

			if (m_com_index < (int)(*this)().index_data().size()) {
				if (m_com_index == -1)
					return zero;
				else if (m_index < (*this)().index_data()[m_com_index])
					return zero;
				else
					return (*this)().value_data()[m_com_index];
			} else
				return zero;
		}


		reference operator *()
		{
			return (*this)()[m_index];
		}


		reference operator [](const size_t i)
		{
			return (*this)()[m_index + i];
		}


		const_reference operator [](const size_t i) const
		{
			return (*this)()[m_index + i];
		}


		/// CPPONLY pre-incrment return by-reference
		iterator & operator++()
		{
			++m_index;
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			if (m_com_index < (int)(*this)().index_data().size() && m_index > (*this)().index_data()[m_com_index])
				++m_com_index;
			return *this;
		}


		const iterator & operator++() const
		{
			++m_index;
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			if (m_com_index < (int)(*this)().index_data().size() && m_index > (*this)().index_data()[m_com_index])
				++m_com_index;
			return *this;
		}


		/// CPPONLY post-incrment return by-value
		iterator operator++(int)
		{
			iterator orig = *this;

			++(*this);
			return orig;
		}


		const iterator operator++(int)  const
		{
			const iterator orig = *this;

			++(*this);
			return orig;
		}


		iterator & operator+=(const iterator & iter)
		{
			m_index += iter.m_index;
			if (m_index > (*this)().size())
				m_index = (*this)().size();
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin() + m_com_index,
					(*this)().index_data().begin() + (*this)().filled(), m_index);
			m_com_index = lower - (*this)().index_data().begin();
			return *this;
		}


		const iterator & operator+=(const iterator & iter)  const
		{
			m_index += iter.m_index;
			if (m_index > (*this)().size())
				m_index = (*this)().size();
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin() + m_com_index, (*this)().index_data().begin() + (*this)().filled(), m_index);
			m_com_index = lower - (*this)().index_data().begin();
			return *this;
		}


		iterator & operator+=(const size_t size)
		{
			m_index += size;
			if (m_index > (*this)().size())
				m_index = (*this)().size();
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin() + m_com_index, (*this)().index_data().begin() + (*this)().filled(), m_index);
			m_com_index = lower - (*this)().index_data().begin();
			return *this;
		}


		const iterator & operator+=(const size_t size)  const
		{
			m_index += size;
			if (m_index > (*this)().size())
				m_index = (*this)().size();
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin() + m_com_index, (*this)().index_data().begin() + (*this)().filled(), m_index);
			m_com_index = lower - (*this)().index_data().begin();
			return *this;
		}


		/*size_t operator + (iterator & iter)
		   {
		    return (*this).m_index + iter.m_index;
		   }
		 */

		iterator operator +(const iterator & iter) const
		{
			iterator result = *this;

			result.m_index += iter.m_index;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin(), (*this)().index_data().begin() + (*this)().filled(), result.m_index);
			result.m_com_index = lower - (*this)().index_data().begin();

			return result;
		}


		iterator operator +(const size_t size) const
		{
			iterator result = *this;

			result.m_index += size;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin(), (*this)().index_data().begin() + (*this)().filled(), result.m_index);
			result.m_com_index = lower - (*this)().index_data().begin();
			return result;
		}


		size_t operator -(const iterator & iter) const
		{

			return m_index - iter.m_index;
		}


		/*
		   iterator operator -(const iterator & iter) const
		   {
		    iterator result = *this;

		    result.m_index -= iter.m_index;
		    if (m_com_index == -1)
		        return result;
		    index_array_type::const_reverse_iterator upper =
		        std::upper_bound((*this)().index_data().rbegin() + ((*this)().index_data().size() - m_com_index), (*this)().index_data().rend(), result.m_index);
		    result.m_com_index = (*this)().index_data().rend() - upper;
		    return result;
		   }
		 */

		iterator operator -(const size_t size) const
		{
			iterator result = *this;

			result.m_index -= size;
			if (m_com_index == -1)
				return result;
			index_array_type::const_reverse_iterator upper =
			    std::upper_bound((*this)().index_data().rbegin() + ((*this)().index_data().size() - m_com_index), (*this)().index_data().rend(), result.m_index);
			result.m_com_index = (*this)().index_data().rend() - upper;
			return result;
		}


		size_t findPositionIndexData()
		{
			for (size_t i = 0; i < (*this)().filled(); i++) {
				if (i < (*this)().filled() - 1 && m_index > (*this)().index_data()[i]) {
					continue;
				} else if (i == (*this)().filled() - 1 && m_index > (*this)().index_data()[i]) {
					return i + 1;
				}else {
					return i;
				}
			}

			return 0;
		}


		/*
		        size_t findPositionIndexData(iterator & it, size_t idxStart = 0) const
		        {
		            for (size_t i = idxStart; i < it.(*this).filled(); i++) {
		                if (i < it.(*this).filled() - 1 && it.m_index > it.(*this).index_data()[i]) {
		                    continue;
		                } else if (i == it.(*this).filled() - 1 && it.m_index > it.(*this).index_data()[i]) {
		                    return i + 1;
		                }else {
		                    return i;
		                }
		            }

		            return 0;
		        }
		 */

		index_array_type::iterator getIndexIterator() const
		{
			//return (*this).index_data().begin() + findPositionIndexData();
			//return std::lower_bound((*this).index_data().begin(), (*this).index_data().begin() + (*this).filled(), m_index);
			if (m_com_index == -1)
				return (*this)().index_data().begin();
			else if (m_com_index > (ssize_t)(*this)().filled())
				return (*this)().index_data().begin() + (*this)().filled();
			else
				return (*this)().index_data().begin() + m_com_index;
		}


		value_array_type::iterator getValueIterator() const
		{
			//return (*this)().value_data().begin() + findPositionIndexData();
			//size_t com_index = getIndexIterator() - (*this)().index_data().begin();
			//return (*this)().value_data().begin() + com_index;
			if (m_com_index == -1)
				return (*this)().value_data().begin();
			else if (m_com_index > (ssize_t)(*this)().filled())
				return (*this)().value_data().begin() + (*this)().filled();
			else
				return (*this)().value_data().begin() + m_com_index;
		}


		val_iterator getCompressedVectorIterator()
		{
			return (*this)().find(m_index);
		}


		size_t getIndex() const
		{
			return m_index;
		}


		size_t getComIndex() const
		{
			return m_com_index;
		}


		vectorm * getContainer()
		{
			return &((*this)());
		}


		/*
		   void deleted()
		   {
		    (*this)().erase_element(m_index);
		   }
		 */

		void assignIfDiffer(const_reference value)
		{
			if (value != this->value())
				(*this)()[m_index] = value;
		}


		friend class const_iterator;

		operator const_iterator() const
		{
			return *reinterpret_cast<const_iterator *>(const_cast<iterator *>(this));
		}
	};


	class const_iterator :
		public container_const_reference<self_type>
	{
protected:
		mutable size_t m_index;
		mutable ssize_t m_com_index;

public:
		typedef std::input_iterator_tag iteratorcategory;
		typedef Allele & reference;
		typedef Allele * pointer;
		typedef long int difference_type;
		typedef Allele value_type;

		const_iterator (const const_iterator & iter) : container_const_reference<self_type>(const_cast<self_type &>(iter())),
			m_index(iter.m_index), m_com_index(iter.m_com_index)
		{
		}


		const_iterator () : container_const_reference<self_type>(), m_index(0), m_com_index(-1)
		{
		}


		const_iterator (self_type & v) : container_const_reference<self_type>(const_cast<self_type &>(v)), m_index(0), m_com_index(-1)
		{
			// m_com_index is the smallest index of mutants (or -1 is nothing is there)
		}


		const_iterator (const self_type & v, size_t index)
			: container_const_reference<self_type>(const_cast<self_type &>(v)), m_index(index)
		{
			if (index == 0) {
				if ((*this)().index_data().size() > 0)
					if ((*this)().index_data()[0] > 0)
						m_com_index = -1;
					else
						m_com_index = 0;
				else
					m_com_index = -1;
			} else {
				index_array_type::const_iterator lower =
				    std::lower_bound((*this)().index_data().begin(),
						(*this)().index_data().begin() + (*this)().filled(), index);
				m_com_index = lower - (*this)().index_data().begin();
			}
			// find the smallest m_com_index  that is larger than index
		}


		const_iterator & operator =(const const_iterator & iter)
		{
			container_const_reference<self_type>::assign(&iter());
			m_index = iter.m_index;
			m_com_index = iter.m_com_index;
			return *this;
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

		const_val_iterator getValIterator() const
		{
			return (*this)().val_begin(m_index);
		}

		const_reference value() const
		{
			static const Allele zero = 0;

			if (m_com_index < (int)(*this)().index_data().size()) {
				if (m_com_index == -1)
					return zero;
				else if (m_index < (*this)().index_data()[m_com_index])
					return zero;
				else
					return (*this)().value_data()[m_com_index];
			} else
				return zero;
		}


		const_reference operator *() const
		{
			static const Allele zero = 0;

			if (m_com_index < (int)(*this)().index_data().size()) {
				if (m_com_index == -1)
					return zero;
				else if (m_index < (*this)().index_data()[m_com_index])
					return zero;
				else
					return (*this)().value_data()[m_com_index];
			} else
				return zero;
		}


		const_reference operator [](const size_t i) const
		{
			return (*this)()[m_index + i];
		}


		/// CPPONLY pre-incrment return by-reference
		const_iterator & operator++()
		{
			++m_index;
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			if (m_com_index < (int)(*this)().index_data().size() && m_index > (*this)().index_data()[m_com_index])
				++m_com_index;
			return *this;
		}


		const const_iterator & operator++() const
		{
			++m_index;
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			if (m_com_index < (int)(*this)().index_data().size() && m_index > (*this)().index_data()[m_com_index])
				++m_com_index;
			return *this;
		}


		/// CPPONLY post-incrment return by-value
		const_iterator operator++(int)
		{
			const_iterator orig = *this;

			++(*this);
			return orig;
		}


		const const_iterator operator++(int)  const
		{
			const const_iterator orig = *this;

			++(*this);
			return orig;
		}


		const_iterator & operator+=(const const_iterator & iter)
		{
			m_index += iter.m_index;
			if (m_index > (*this)().size())
				m_index = (*this)().size();
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin() + m_com_index,
					(*this)().index_data().begin() + (*this)().filled(), m_index);
			m_com_index = lower - (*this)().index_data().begin();
			return *this;
		}


		const const_iterator & operator+=(const const_iterator & iter)  const
		{
			m_index += iter.m_index;
			if (m_index > (*this)().size())
				m_index = (*this)().size();
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin() + m_com_index, (*this)().index_data().begin() + (*this)().filled(), m_index);
			m_com_index = lower - (*this)().index_data().begin();
			return *this;
		}


		const_iterator & operator+=(const size_t size)
		{
			m_index += size;
			if (m_index > (*this)().size())
				m_index = (*this)().size();
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin() + m_com_index, (*this)().index_data().begin() + (*this)().filled(), m_index);
			m_com_index = lower - (*this)().index_data().begin();
			return *this;
		}


		const const_iterator & operator+=(const size_t size)  const
		{
			m_index += size;
			if (m_index > (*this)().size())
				m_index = (*this)().size();
			if ((*this)().index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin() + m_com_index, (*this)().index_data().begin() + (*this)().filled(), m_index);
			m_com_index = lower - (*this)().index_data().begin();
			return *this;
		}


		/*size_t operator + (iterator & iter)
		   {
		    return (*this).m_index + iter.m_index;
		   }
		 */

		const_iterator operator +(const const_iterator & iter) const
		{
			const_iterator result = *this;

			result.m_index += iter.m_index;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin(), (*this)().index_data().begin() + (*this)().filled(), result.m_index);
			result.m_com_index = lower - (*this)().index_data().begin();

			return result;
		}


		const_iterator operator +(const size_t size) const
		{
			const_iterator result = *this;

			result.m_index += size;
			index_array_type::const_iterator lower =
			    std::lower_bound((*this)().index_data().begin(), (*this)().index_data().begin() + (*this)().filled(), result.m_index);
			result.m_com_index = lower - (*this)().index_data().begin();
			return result;
		}


		size_t operator -(const const_iterator & iter) const
		{

			return m_index - iter.m_index;
		}


		/*
		   iterator operator -(const iterator & iter) const
		   {
		    iterator result = *this;

		    result.m_index -= iter.m_index;
		    if (m_com_index == -1)
		        return result;
		    index_array_type::const_reverse_iterator upper =
		        std::upper_bound((*this)().index_data().rbegin() + ((*this)().index_data().size() - m_com_index), (*this)().index_data().rend(), result.m_index);
		    result.m_com_index = (*this)().index_data().rend() - upper;
		    return result;
		   }
		 */

		const_iterator operator -(const size_t size) const
		{
			const_iterator result = *this;

			result.m_index -= size;
			if (m_com_index == -1)
				return result;
			index_array_type::const_reverse_iterator upper =
			    std::upper_bound((*this)().index_data().rbegin() + ((*this)().index_data().size() - m_com_index), (*this)().index_data().rend(), result.m_index);
			result.m_com_index = (*this)().index_data().rend() - upper;
			return result;
		}


		size_t findPositionIndexData()
		{
			for (size_t i = 0; i < (*this)().filled(); i++) {
				if (i < (*this)().filled() - 1 && m_index > (*this)().index_data()[i]) {
					continue;
				} else if (i == (*this)().filled() - 1 && m_index > (*this)().index_data()[i]) {
					return i + 1;
				}else {
					return i;
				}
			}

			return 0;
		}


		/*
		        size_t findPositionIndexData(iterator & it, size_t idxStart = 0) const
		        {
		            for (size_t i = idxStart; i < it.(*this).filled(); i++) {
		                if (i < it.(*this).filled() - 1 && it.m_index > it.(*this).index_data()[i]) {
		                    continue;
		                } else if (i == it.(*this).filled() - 1 && it.m_index > it.(*this).index_data()[i]) {
		                    return i + 1;
		                }else {
		                    return i;
		                }
		            }

		            return 0;
		        }
		 */

		index_array_type::const_iterator getIndexIterator() const
		{
			//return (*this).index_data().begin() + findPositionIndexData();
			//return std::lower_bound((*this).index_data().begin(), (*this).index_data().begin() + (*this).filled(), m_index);
			if (m_com_index == -1)
				return (*this)().index_data().begin();
			else if (m_com_index > (ssize_t)(*this)().filled())
				return (*this)().index_data().begin() + (*this)().filled();
			else
				return (*this)().index_data().begin() + m_com_index;
		}


		value_array_type::const_iterator getValueIterator() const
		{
			//return (*this)().value_data().begin() + findPositionIndexData();
			//size_t com_index = getIndexIterator() - (*this)().index_data().begin();
			//return (*this)().value_data().begin() + com_index;
			if (m_com_index == -1)
				return (*this)().value_data().begin();
			else if (m_com_index > (ssize_t)(*this)().filled())
				return (*this)().value_data().begin() + (*this)().filled();
			else
				return (*this)().value_data().begin() + m_com_index;
		}


		const_val_iterator getCompressedVectorIterator()
		{
			return (*this)().find(m_index);
		}


		size_t getIndex() const
		{
			return m_index;
		}


		size_t getComIndex() const
		{
			return m_com_index;
		}


		const vectorm * getContainer()
		{
			return &((*this)());
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


	void erase(iterator begin, iterator end)
	{
		if (end != this->end()) {
			std::copy(end.getValueIterator(), this->end().getValueIterator(), begin.getValueIterator());
			index_array_type::iterator it = end.getIndexIterator();
			index_array_type::iterator index_end = end.getIndexIterator();
			for (; it != this->end().getIndexIterator(); ++it) {
				*it -= (end - begin);
			}
			std::copy(index_end, this->end().getIndexIterator(), begin.getIndexIterator());
		}
		resize(size() - (end.getIndex() - begin.getIndex()));
	}


	void insert(iterator it, const_iterator begin, const_iterator end)
	{
		index_array_type::const_iterator src_index_begin = begin.getIndexIterator();
		index_array_type::const_iterator src_index_end = end.getIndexIterator();
		value_array_type::const_iterator src_value_begin = begin.getValueIterator();
		value_array_type::const_iterator src_value_end = end.getValueIterator();
		size_t insert_size = end - begin;
		size_t index_insert_size = src_index_end - src_index_begin;
		size_t filled_size = filled();

		resize(size() + insert_size);
		if (filled_size + index_insert_size >= it.getContainer()->nnz_capacity()) {
			reserve(2 * filled_size + index_insert_size, true);
		}

		index_array_type::iterator dest_index_begin = it.getIndexIterator();
		index_array_type::iterator dest_index_end = this->end().getIndexIterator();
		value_array_type::iterator dest_value_begin = it.getValueIterator();

		std::copy_backward(dest_index_begin, it.getContainer()->index_data().begin() + filled_size, it.getContainer()->index_data().begin() + filled_size + index_insert_size);
		for (size_t i = 0; i < (size_t)(dest_index_end - dest_index_begin); i++) {
			*(dest_index_begin + index_insert_size + i) += index_insert_size;
		}
		for (size_t i = 0; i < index_insert_size; i++) {
			size_t range = *(src_index_begin + i) - begin.getIndex();
			*(dest_index_begin + i) = it.getIndex() + range;
		}
		std::copy_backward(dest_value_begin, it.getContainer()->value_data().begin() + filled_size, it.getContainer()->value_data().begin() + filled_size + index_insert_size);
		std::copy(src_value_begin, src_value_end, dest_value_begin);
		set_filled(filled_size + index_insert_size);
	}


private:
	void storage_invariants() const
	{
		BOOST_UBLAS_CHECK(capacity_ == index_data_.size(), internal_logic());
		BOOST_UBLAS_CHECK(capacity_ == value_data_.size(), internal_logic());
		BOOST_UBLAS_CHECK(filled_ <= capacity_, internal_logic());
		BOOST_UBLAS_CHECK((0 == filled_) || (index_data_[filled_ - 1] < size_), internal_logic());
	}


	size_type size_;
	index_array_type::size_type capacity_;
	index_array_type::size_type filled_;
	index_array_type index_data_;
	value_array_type value_data_;
	static const value_type zero_;

	friend class val_iterator;
	friend class const_val_iterator;
};


//typedef unsigned char Allele;
typedef vectorm::reference AlleleRef;

inline vectorm::iterator find(vectorm::iterator begin, vectorm::iterator end, Allele value)
{
	vectorm::iterator it = begin;

	for (; it != end; ++it) {
		if (*it == value)
			return it;
	}
	return end;

}


}
#endif

#endif
