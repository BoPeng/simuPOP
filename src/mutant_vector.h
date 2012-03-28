#include "boost/numeric/ublas/vector_sparse.hpp"

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

#ifndef _MUTANT_VECTOR_H
#  define _MUTANT_VECTOR_H

namespace simuPOP {
template<class T, class IA = unbounded_array<std::size_t>, class TA = unbounded_array<T> >
class compressed_vector;

template<class I, class T, class C>
inline I _lower_bound(const I & begin, const I & end, const T & t, C compare)
{
	// t <= *begin <=> ! (*begin < t)
	if (begin == end || !compare(*begin, t))
		return begin;
	if (compare(*(end - 1), t))
		return end;
	return std::lower_bound(begin, end, t, compare);
}


template<class T, class IA, class TA>
class compressed_vector :
	public vector_container<compressed_vector<T, IA, TA> >
{

	typedef T & true_reference;
	typedef T * pointer;
	typedef const T * const_pointer;
	typedef compressed_vector<T, IA, TA> self_type;

public:
#  ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
	using vector_container<self_type>::operator ();
#  endif
	// ISSUE require type consistency check
	// is_convertable (IA::size_type, TA::size_type)
	typedef typename IA::value_type size_type;
	typedef typename IA::difference_type difference_type;
	typedef T value_type;
	typedef const T & const_reference;
	typedef T & reference;
	typedef IA index_array_type;
	typedef TA value_array_type;
	typedef const vector_reference<const self_type> const_closure_type;
	typedef vector_reference<self_type> closure_type;
	typedef self_type vector_temporary_type;
	typedef sparse_tag storage_category;

	// Construction and destruction
	inline compressed_vector () :
		vector_container<self_type> (),
		size_(0), capacity_(restrict_capacity(0)), filled_(0),
		index_data_(capacity_), value_data_(capacity_)
	{
		storage_invariants();
	}


	explicit inline compressed_vector (size_type size, size_type non_zeros = 0) :
		vector_container<self_type> (),
		size_(size), capacity_(restrict_capacity(non_zeros)), filled_(0),
		index_data_(capacity_), value_data_(capacity_)
	{
		storage_invariants();
	}


	inline compressed_vector (const compressed_vector & v) :
		vector_container<self_type> (),
		size_(v.size_), capacity_(v.capacity_), filled_(v.filled_),
		index_data_(v.index_data_), value_data_(v.value_data_)
	{
		storage_invariants();
	}


	template<class AE>
	inline compressed_vector (const vector_expression<AE> & ae, size_type non_zeros = 0) :
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


	inline typename index_array_type::size_type filled() const
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


	inline void set_filled(const typename index_array_type::size_type & filled)
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
		}else {
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
		typename std::iterator_traits<subiterator_type>::difference_type n = it - index_data_.begin();
		BOOST_UBLAS_CHECK(filled_ == 0 || filled_ == typename index_array_type::size_type(n) || *it != i, internal_logic());           // duplicate found by _lower_bound
		++filled_;
		it = index_data_.begin() + n;
		std::copy_backward(it, index_data_.begin() + filled_ - 1, index_data_.begin() + filled_);
		*it = i;
		typename value_array_type::iterator itt(value_data_.begin() + n);
		std::copy_backward(itt, value_data_.begin() + filled_ - 1, value_data_.begin() + filled_);
		*itt = t;
		storage_invariants();
		return *itt;
	}


	inline void erase_element(size_type i)
	{
		subiterator_type it(_lower_bound(index_data_.begin(), index_data_.begin() + filled_, i, std::less<size_type> ()));

		typename std::iterator_traits<subiterator_type>::difference_type n = it - index_data_.begin();
		if (filled_ > typename index_array_type::size_type(n) && *it == i) {
			std::copy(it + 1, index_data_.begin() + filled_, it);
			typename value_array_type::iterator itt(value_data_.begin() + n);
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
	inline compressed_vector & operator =(const compressed_vector & v)
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


	template<class C>              // Container assignment without temporary
	inline compressed_vector & operator =(const vector_container<C> & v)
	{
		resize(v().size(), false);
		assign(v);
		return *this;
	}


	inline compressed_vector & assign_temporary(compressed_vector & v)
	{
		swap(v);
		return *this;
	}


	template<class AE>
	inline compressed_vector & operator =(const vector_expression<AE> & ae)
	{
		self_type temporary(ae, capacity_);

		return assign_temporary(temporary);
	}


	template<class AE>
	inline compressed_vector & assign(const vector_expression<AE> & ae)
	{
		vector_assign<scalar_assign> (*this, ae);
		return *this;
	}


	// Computed assignment
	template<class AE>
	inline compressed_vector & operator +=(const vector_expression<AE> & ae)
	{
		self_type temporary(*this + ae, capacity_);

		return assign_temporary(temporary);
	}


	template<class C>              // Container assignment without temporary
	inline compressed_vector & operator +=(const vector_container<C> & v)
	{
		plus_assign(v);
		return *this;
	}


	template<class AE>
	inline compressed_vector & plus_assign(const vector_expression<AE> & ae)
	{
		vector_assign<scalar_plus_assign> (*this, ae);
		return *this;
	}


	template<class AE>
	inline compressed_vector & operator -=(const vector_expression<AE> & ae)
	{
		self_type temporary(*this - ae, capacity_);

		return assign_temporary(temporary);
	}


	template<class C>              // Container assignment without temporary
	inline compressed_vector & operator -=(const vector_container<C> & v)
	{
		minus_assign(v);
		return *this;
	}


	template<class AE>
	inline compressed_vector & minus_assign(const vector_expression<AE> & ae)
	{
		vector_assign<scalar_minus_assign> (*this, ae);
		return *this;
	}


	template<class AT>
	inline compressed_vector & operator *=(const AT & at)
	{
		vector_assign_scalar<scalar_multiplies_assign> (*this, at);
		return *this;
	}


	template<class AT>
	inline compressed_vector & operator /=(const AT & at)
	{
		vector_assign_scalar<scalar_divides_assign> (*this, at);
		return *this;
	}


	// Swapping
	inline void swap(compressed_vector & v)
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


	inline friend void swap(compressed_vector & v1, compressed_vector & v2)
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


	inline void pop_back()
	{
		BOOST_UBLAS_CHECK(filled_ > 0, external_logic());
		--filled_;
		storage_invariants();
	}


	// Iterator types

private:
	// Use index array iterator
	typedef typename IA::const_iterator const_subiterator_type;
	typedef typename IA::iterator subiterator_type;

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
		public container_const_reference<compressed_vector>,
		public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
		                                   const_val_iterator, value_type>
	{
public:
		typedef typename compressed_vector::value_type value_type;
		typedef typename compressed_vector::difference_type difference_type;
		typedef typename compressed_vector::const_reference reference;
		typedef const typename compressed_vector::pointer pointer;

		// Construction and destruction
		inline const_val_iterator () :
			container_const_reference<self_type> (), it_() {}
		inline const_val_iterator (const self_type & v, const const_subiterator_type & it) :
			container_const_reference<self_type> (v), it_(it) {}
		inline const_val_iterator (const typename self_type::iterator & it) :          // ISSUE self_type:: stops VC8 using std::iterator here
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


	inline const_val_iterator val_end() const
	{
		return find(size_);
	}


	class iterator :
		public container_reference<compressed_vector>,
		public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
		                                   iterator, value_type>
	{
public:
		typedef typename compressed_vector::value_type value_type;
		typedef typename compressed_vector::difference_type difference_type;
		typedef typename compressed_vector::true_reference reference;
		typedef typename compressed_vector::pointer pointer;

		// Construction and destruction
		inline iterator () :
			container_reference<self_type> (), it_() {}
		inline iterator (self_type & v, const subiterator_type & it) :
			container_reference<self_type> (v), it_(it) {}

		// Arithmetic
		inline iterator & operator ++()
		{
			++it_;
			return *this;
		}


		inline iterator & operator --()
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
		inline iterator & operator =(const iterator & it)
		{
			container_reference<self_type>::assign(&it());
			it_ = it.it_;
			return *this;
		}


		// Comparison
		inline bool operator ==(const iterator & it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			return it_ == it.it_;
		}


private:
		subiterator_type it_;

		friend class const_val_iterator;
	};

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


private:
	void storage_invariants() const
	{
		BOOST_UBLAS_CHECK(capacity_ == index_data_.size(), internal_logic());
		BOOST_UBLAS_CHECK(capacity_ == value_data_.size(), internal_logic());
		BOOST_UBLAS_CHECK(filled_ <= capacity_, internal_logic());
		BOOST_UBLAS_CHECK((0 == filled_) || (index_data_[filled_ - 1] < size_), internal_logic());
	}


	size_type size_;
	typename index_array_type::size_type capacity_;
	typename index_array_type::size_type filled_;
	index_array_type index_data_;
	value_array_type value_data_;
	static const value_type zero_;

	friend class val_iterator;
	friend class const_val_iterator;
};

template<class T, class IA, class TA>
const typename compressed_vector<T, IA, TA>::value_type compressed_vector<T, IA, TA>::zero_ = value_type/*zero*/ ();


/// CPPONLY
template <class T>
class mutant_vector
{
public:
	class iterator;
	typedef const iterator const_iterator;

	mutant_vector()
	{
	}


	mutant_vector(size_t size)
	{
		m_container.resize(size);
	}


	mutant_vector(size_t size, const T & value)
	{
		if (value == 0)
			m_container.resize(size);
		else {
			m_container.resize(size);
			m_container.reserve(size);
			for (size_t i = 0; i < size; i++)
				m_container.push_back(i, value);
		}
	}


	void resize(size_t size)
	{
		m_container.resize(size);
	}


	void reserve(size_t size)
	{
		m_container.reserve(size);
	}


	size_t size() const
	{
		return m_container.size();
	}


	typename mutant_vector<T>::iterator  begin()
	{
		return typename mutant_vector<T>::iterator(&m_container, 0);
	}


	const typename mutant_vector<T>::iterator  begin() const
	{
		return typename mutant_vector<T>::iterator(const_cast<compressed_vector<T> *>(&m_container), 0);
	}


	typename mutant_vector<T>::iterator  end()
	{
		return typename mutant_vector<T>::iterator(&m_container, m_container.size());
	}


	const typename mutant_vector<T>::iterator  end() const
	{
		return typename mutant_vector<T>::iterator(const_cast<compressed_vector<T> *>(&m_container), m_container.size());
	}


	typename compressed_vector<T>::reference operator [](size_t i)
	{
		return m_container[i];
	}


	void push_back(typename compressed_vector<T>::size_type i, typename compressed_vector<T>::const_reference t)
	{
		m_container.push_back(i, t);
	}


	void swap(mutant_vector<T> & vec)
	{
		m_container.swap(vec.getContainer());

	}


	compressed_vector<T> & getContainer()
	{
		return m_container;
	}


	const compressed_vector<T> & getContainer() const
	{
		return m_container;
	}


	class iterator
	{
protected:
		compressed_vector<T> * m_container;
		mutable size_t m_index;
		mutable ssize_t m_com_index;

public:
		typedef std::input_iterator_tag iterator_category;
		typedef typename compressed_vector<T>::reference reference;
		typedef T * pointer;
		typedef long int difference_type;
		typedef T value_type;

		iterator (const iterator & iter)
		{
			m_container = iter.m_container;
			m_index = iter.m_index;
			m_com_index = iter.m_com_index;
			/*	if (m_index == 0) {
			        if (m_container->index_data().size() > 0)
			            m_com_index = 0;
			        else
			            m_com_index = -1;
			    } else {
			        typename compressed_vector<size_t>::index_array_type::iterator lower =
			            std::lower_bound(m_container->index_data().begin(), m_container->index_data().begin() + m_container->filled(), m_index);
			        m_com_index = lower - m_container->index_data().begin();
			    }
			 */
		}


		iterator (compressed_vector<T> * c) : m_container(c), m_index(0), m_com_index(-1)
		{
			// m_com_index is the smallest index of mutants (or -1 is nothing is there)
		}


		iterator (const compressed_vector<T> * c) : m_container(c), m_index(0), m_com_index(-1)
		{
		}


		iterator (compressed_vector<T> * c, size_t index) : m_container(c), m_index(index)
		{
			if (index == 0) {
				if (m_container->index_data().size() > 0)
					if (m_container->index_data()[0] > 0)
						m_com_index = -1;
					else
						m_com_index = 0;
				else
					m_com_index = -1;
			} else {
				typename compressed_vector<size_t>::index_array_type::iterator lower =
				    std::lower_bound(m_container->index_data().begin(), m_container->index_data().begin() + m_container->filled(), index);
				m_com_index = lower - m_container->index_data().begin();
			}
			// find the smallest m_com_index  that is larger than index
		}


		iterator (const compressed_vector<T> * c, size_t index) : m_container(c), m_index(index)
		{
			if (index == 0) {
				if (m_container->index_data().size() > 0)
					if (m_container->index_data()[0] > 0)
						m_com_index = -1;
					else
						m_com_index = 0;
				else
					m_com_index = -1;
			} else {
				typename compressed_vector<size_t>::index_array_type::iterator lower =
				    std::lower_bound(m_container->index_data().begin(), m_container->index_data().begin() + m_container->filled(), index);
				m_com_index = lower - m_container->index_data().begin();
			}
		}


		iterator () : m_container(NULL), m_index(0), m_com_index(-1)
		{

		}


		iterator & operator =(const iterator & iter)
		{
			m_container = iter.m_container;
			m_index = iter.m_index;
			m_com_index = iter.m_com_index;
			return *this;
		}


		bool operator ==(const iterator & iter) const
		{
			if (m_index == iter.m_index)
				return true;
			return false;
		}


		bool operator !=(const iterator & iter) const
		{
			if (m_index != iter.m_index)
				return true;
			return false;
		}


		bool operator >=(const iterator & iter) const
		{
			if (m_index >= iter.m_index)
				return true;
			return false;
		}


		bool operator <=(const iterator & iter) const
		{
			if (m_index <= iter.m_index)
				return true;
			return false;
		}


		bool operator >(const iterator & iter) const
		{
			if (m_index > iter.m_index)
				return true;
			return false;
		}


		bool operator <(const iterator & iter) const
		{
			if (m_index < iter.m_index)
				return true;
			return false;
		}


		typename compressed_vector<T>::const_reference operator*() const
		{
			static const T zero = 0;

			if (m_com_index < (int)m_container->index_data().size()) {
				if (m_com_index == -1)
					return zero;
				else if (m_index < m_container->index_data()[m_com_index])
					return zero;
				else
					return m_container->value_data()[m_com_index];
			} else
				return zero;

		}


		typename compressed_vector<T>::reference operator*()
		{
			return (*m_container)[m_index];
		}


		typename compressed_vector<T>::reference operator [](const size_t i)
		{
			return (*m_container)[m_index + i];
		}


		typename compressed_vector<T>::const_reference operator [](const size_t i) const
		{
			return (*m_container)[m_index + i];
		}


		/// CPPONLY pre-incrment return by-reference
		iterator & operator++()
		{
			++m_index;
			if (m_container->index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			if (m_com_index < (int)m_container->index_data().size() && m_index > m_container->index_data()[m_com_index])
				++m_com_index;
			return *this;
		}


		const iterator & operator++() const
		{
			++m_index;
			if (m_container->index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			if (m_com_index < (int)m_container->index_data().size() && m_index > m_container->index_data()[m_com_index])
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
			if (m_index > m_container->size())
				m_index = m_container->size();
			if (m_container->index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			typename compressed_vector<T>::index_array_type::const_iterator lower =
			    std::lower_bound(m_container->index_data().begin() + m_com_index, m_container->index_data().begin() + m_container->filled(), m_index);
			m_com_index = lower - m_container->index_data().begin();
			return *this;
		}


		const iterator & operator+=(const iterator & iter)  const
		{
			m_index += iter.m_index;
			if (m_index > m_container->size())
				m_index = m_container->size();
			if (m_container->index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			typename compressed_vector<T>::index_array_type::const_iterator lower =
			    std::lower_bound(m_container->index_data().begin() + m_com_index, m_container->index_data().begin() + m_container->filled(), m_index);
			m_com_index = lower - m_container->index_data().begin();
			return *this;
		}


		iterator & operator+=(const size_t size)
		{
			m_index += size;
			if (m_index > m_container->size())
				m_index = m_container->size();
			if (m_container->index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			typename compressed_vector<T>::index_array_type::const_iterator lower =
			    std::lower_bound(m_container->index_data().begin() + m_com_index, m_container->index_data().begin() + m_container->filled(), m_index);
			m_com_index = lower - m_container->index_data().begin();
			return *this;
		}


		const iterator & operator+=(const size_t size)  const
		{
			m_index += size;
			if (m_index > m_container->size())
				m_index = m_container->size();
			if (m_container->index_data().size() == 0)
				return *this;
			if (m_com_index == -1)
				m_com_index = 0;
			typename compressed_vector<T>::index_array_type::const_iterator lower =
			    std::lower_bound(m_container->index_data().begin() + m_com_index, m_container->index_data().begin() + m_container->filled(), m_index);
			m_com_index = lower - m_container->index_data().begin();
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
			typename compressed_vector<T>::index_array_type::const_iterator lower =
			    std::lower_bound(m_container->index_data().begin(), m_container->index_data().begin() + m_container->filled(), result.m_index);
			result.m_com_index = lower - m_container->index_data().begin();

			return result;
		}


		iterator operator +(const size_t size) const
		{
			iterator result = *this;

			result.m_index += size;
			typename compressed_vector<T>::index_array_type::const_iterator lower =
			    std::lower_bound(m_container->index_data().begin(), m_container->index_data().begin() + m_container->filled(), result.m_index);
			result.m_com_index = lower - m_container->index_data().begin();
			return result;
		}


		size_t operator -(iterator & iter)
		{

			return (*this).m_index - iter.m_index;
		}


		iterator operator -(const iterator & iter) const
		{
			iterator result = *this;

			result.m_index -= iter.m_index;
			if (m_com_index == -1)
				return result;
			typename compressed_vector<T>::index_array_type::const_reverse_iterator upper =
			    std::upper_bound(m_container->index_data().rbegin() + (m_container->index_data().size() - m_com_index), m_container->index_data().rend(), result.m_index);
			result.m_com_index = m_container->index_data().rend() - upper;
			return result;
		}


		iterator operator -(const size_t size) const
		{
			iterator result = *this;

			result.m_index -= size;
			if (m_com_index == -1)
				return result;
			typename compressed_vector<T>::index_array_type::const_reverse_iterator upper =
			    std::upper_bound(m_container->index_data().rbegin() + (m_container->index_data().size() - m_com_index), m_container->index_data().rend(), result.m_index);
			result.m_com_index = m_container->index_data().rend() - upper;
			return result;
		}


		size_t findPositionIndexData()
		{
			for (size_t i = 0; i < m_container->filled(); i++) {
				if (i < m_container->filled() - 1 && m_index > m_container->index_data()[i]) {
					continue;
				} else if (i == m_container->filled() - 1 && m_index > m_container->index_data()[i]) {
					return i + 1;
				}else {
					return i;
				}
			}

			return 0;
		}


		size_t findPositionIndexData(iterator & it, size_t idxStart = 0) const
		{
			for (size_t i = idxStart; i < it.m_container->filled(); i++) {
				if (i < it.m_container->filled() - 1 && it.m_index > it.m_container->index_data()[i]) {
					continue;
				} else if (i == it.m_container->filled() - 1 && it.m_index > it.m_container->index_data()[i]) {
					return i + 1;
				}else {
					return i;
				}
			}

			return 0;
		}


		typename compressed_vector<T>::index_array_type::iterator getIndexIterator()
		{
			//return m_container->index_data().begin() + findPositionIndexData();
			//return std::lower_bound(m_container->index_data().begin(), m_container->index_data().begin() + m_container->filled(), m_index);
			if (m_com_index == -1)
				return m_container->index_data().begin();
			else if (m_com_index > (ssize_t)m_container->filled())
				return m_container->index_data().begin() + m_container->filled();
			else
				return m_container->index_data().begin() + m_com_index;
		}


		typename compressed_vector<T>::value_array_type::iterator getValueIterator()
		{
			//return m_container->value_data().begin() + findPositionIndexData();
			//size_t com_index = getIndexIterator() - m_container->index_data().begin();
			//return m_container->value_data().begin() + com_index;
			if (m_com_index == -1)
				return m_container->value_data().begin();
			else if (m_com_index > (ssize_t)m_container->filled())
				return m_container->value_data().begin() + m_container->filled();
			else
				return m_container->value_data().begin() + m_com_index;
		}


		typename compressed_vector<T>::index_array_type::const_iterator getIndexIterator() const
		{
			//return m_container->index_data().begin() + findPositionIndexData();
			//return std::lower_bound(m_container->index_data().begin(), m_container->index_data().begin() + m_container->filled(), m_index);
			if (m_com_index == -1)
				return m_container->index_data().begin();
			else if (m_com_index > (ssize_t)m_container->filled())
				return m_container->index_data().begin() + m_container->filled();
			else
				return m_container->index_data().begin() + m_com_index;
		}


		typename compressed_vector<T>::value_array_type::const_iterator getValueIterator() const
		{
			//return m_container->value_data().begin() + findPositionIndexData();
			//size_t com_index = getIndexIterator() - m_container->index_data().begin();
			//return m_container->value_data().begin() + com_index;
			if (m_com_index == -1)
				return m_container->value_data().begin();
			else if (m_com_index > (ssize_t)m_container->filled())
				return m_container->value_data().begin() + m_container->filled();
			else
				return m_container->value_data().begin() + m_com_index;
		}


		typename compressed_vector<T>::iterator getCompressedVectorIterator()
		{
			return m_container->find(m_index);
		}


		size_t getIndex() const
		{
			return m_index;
		}


		size_t getComIndex() const
		{
			return m_com_index;
		}


		compressed_vector<T> * getContainer()
		{
			return m_container;
		}


		void deleted()
		{
			m_container->erase_element(m_index);
		}


		void assign(typename compressed_vector<T>::const_reference value)
		{
			if (value == 0u && (*m_container)[m_index] == 0u)
				return;
			else
				(*m_container)[m_index] = value;
		}


	};

	void erase(typename mutant_vector<T>::iterator begin, typename mutant_vector<T>::iterator end)
	{

		if (end != this->end()) {
			std::copy(end.getValueIterator(), this->end().getValueIterator(), begin.getValueIterator());
			typename compressed_vector<T>::index_array_type::iterator it = end.getIndexIterator();
			typename compressed_vector<T>::index_array_type::iterator index_end = end.getIndexIterator();
			for (; it != this->end().getIndexIterator(); ++it) {
				*it -= (end - begin);
			}
			std::copy(index_end, this->end().getIndexIterator(), begin.getIndexIterator());
		}
		m_container.resize(m_container.size() - (end.getIndex() - begin.getIndex()));
	}


	void insert(typename mutant_vector<T>::iterator it, typename mutant_vector<T>::iterator begin, typename mutant_vector<T>::iterator end)
	{
		typename compressed_vector<T>::index_array_type::iterator src_index_begin = begin.getIndexIterator();
		typename compressed_vector<T>::index_array_type::iterator src_index_end = end.getIndexIterator();
		typename compressed_vector<T>::value_array_type::iterator src_value_begin = begin.getValueIterator();
		typename compressed_vector<T>::value_array_type::iterator src_value_end = end.getValueIterator();
		size_t insert_size = end - begin;
		size_t index_insert_size = src_index_end - src_index_begin;
		size_t filled_size = m_container.filled();

		m_container.resize(m_container.size() + insert_size);
		if (filled_size + index_insert_size >= it.getContainer()->nnz_capacity()) {
			m_container.reserve(2 * filled_size + index_insert_size, true);
		}

		typename compressed_vector<T>::index_array_type::iterator dest_index_begin = it.getIndexIterator();
		typename compressed_vector<T>::index_array_type::iterator dest_index_end = this->end().getIndexIterator();
		typename compressed_vector<T>::value_array_type::iterator dest_value_begin = it.getValueIterator();

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
		m_container.set_filled(filled_size + index_insert_size);
	}


	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*file_version*/)
	{
		ar & m_container;
	}


private:
	compressed_vector<T> m_container;
};

}

#  ifdef MUTANTALLELE
typedef unsigned int Allele;
typedef simuPOP::compressed_vector<Allele>::reference AlleleRef;
typedef simuPOP::mutant_vector<Allele> mutant_vectora;

namespace simuPOP {

inline void copy(mutant_vectora::iterator begin, mutant_vectora::iterator end, mutant_vectora::iterator it)
{

	size_t size = end - begin;
	mutant_vectora::iterator it_end = it + size;

	if (it_end.getIndex() <= it_end.getContainer()->size()) {
		compressed_vector<Allele>::index_array_type::const_iterator src_index_begin = begin.getIndexIterator();
		compressed_vector<Allele>::index_array_type::const_iterator src_index_end = end.getIndexIterator();
		compressed_vector<Allele>::value_array_type::const_iterator src_value_begin = begin.getValueIterator();
		compressed_vector<Allele>::value_array_type::const_iterator src_value_end = end.getValueIterator();
		compressed_vector<Allele>::index_array_type::iterator dest_index_begin = it.getIndexIterator();
		compressed_vector<Allele>::index_array_type::iterator dest_index_end = it_end.getIndexIterator();
		compressed_vector<Allele>::value_array_type::iterator dest_value_begin = it.getValueIterator();
		size_t src_size = src_index_end - src_index_begin;
		size_t dest_size = dest_index_end - dest_index_begin;
		size_t src_idx_num_begin = begin.getIndex();
		size_t dest_idx_num_begin = it.getIndex();
		if (src_size > dest_size) {
			size_t diff_size = src_size - dest_size;
			size_t filled_size = it.getContainer()->filled();
			if (filled_size + diff_size >= it.getContainer()->nnz_capacity()) {
				it.getContainer()->reserve(2 * filled_size + diff_size, true);
				//After using reserve(), dest_value_begin no longer valid, get a new one
				dest_value_begin = it.getValueIterator();
				dest_index_begin = it.getIndexIterator();
				src_index_begin = begin.getIndexIterator();
			}
			std::copy_backward(dest_index_begin, it.getContainer()->index_data().begin() + filled_size, it.getContainer()->index_data().begin() + filled_size + diff_size);
			for (size_t i = 0; i < src_size; i++) {
				size_t range = *(src_index_begin + i) - src_idx_num_begin;
				*(dest_index_begin + i) = dest_idx_num_begin + range;
			}
			std::copy_backward(dest_value_begin, it.getContainer()->value_data().begin() + filled_size, it.getContainer()->value_data().begin() + filled_size + diff_size);
			std::copy(begin.getValueIterator(), end.getValueIterator(), dest_value_begin);
			it.getContainer()->set_filled(filled_size + diff_size);
		} else if (src_size < dest_size) {
			size_t diff_size = dest_size - src_size;
			size_t filled_size = it.getContainer()->filled();
			std::copy(dest_index_begin + diff_size, it.getContainer()->index_data().begin() + filled_size, dest_index_begin);
			for (size_t i = 0; i < src_size; i++) {
				size_t range = *(src_index_begin + i) - src_idx_num_begin;
				*(dest_index_begin + i) = dest_idx_num_begin + range;
			}
			std::copy(dest_value_begin + diff_size, it.getContainer()->value_data().begin() + filled_size, dest_value_begin);
			std::copy(src_value_begin, src_value_end, dest_value_begin);
			it.getContainer()->set_filled(filled_size - diff_size);

		} else {
			for (size_t i = 0; i < src_size; i++) {
				size_t range = *(src_index_begin + i) - src_idx_num_begin;
				*(dest_index_begin + i) = dest_idx_num_begin + range;
			}
			std::copy(src_value_begin, src_value_end, dest_value_begin);
		}
	}
}


inline void fill(mutant_vectora::iterator begin, mutant_vectora::iterator end, Allele value)
{
	compressed_vector<Allele>::index_array_type::iterator index_begin = begin.getIndexIterator();
	compressed_vector<Allele>::index_array_type::iterator index_end = end.getIndexIterator();
	compressed_vector<Allele>::value_array_type::iterator value_begin = begin.getValueIterator();
	size_t diff_size = index_end - index_begin;
	size_t filled_size = begin.getContainer()->filled();
	if (value == 0) {
		std::copy(index_begin + diff_size, begin.getContainer()->index_data().begin() + filled_size, index_begin);
		std::copy(value_begin + diff_size, begin.getContainer()->value_data().begin() + filled_size, value_begin);
		begin.getContainer()->set_filled(filled_size - diff_size);
		//  If the fill function is extensively use, reserve() can be remove in order to improve performance .
		begin.getContainer()->reserve(filled_size - diff_size);
	} else {
		size_t value_size = end - begin;
		size_t insert_size = value_size - diff_size;
		if (filled_size + insert_size >= begin.getContainer()->nnz_capacity()) {
			begin.getContainer()->reserve(filled_size + insert_size, true);
			// After using reserve(), index_begin and value_begin no longer valid, get a new one
			index_begin = begin.getIndexIterator();
			value_begin = begin.getValueIterator();
		}
		std::copy_backward(index_begin, begin.getContainer()->index_data().begin() + filled_size, begin.getContainer()->index_data().begin() + filled_size + insert_size);
		for (size_t i = 0; i < value_size; i++) {
			*(index_begin + i) = *index_begin + i;
		}
		std::copy_backward(value_begin, begin.getContainer()->value_data().begin() + filled_size, begin.getContainer()->value_data().begin() + filled_size + insert_size);
		std::fill(value_begin, value_begin + value_size, value);
		begin.getContainer()->set_filled(filled_size + insert_size);

	}
}


inline mutant_vectora::iterator find(mutant_vectora::iterator begin, mutant_vectora::iterator end, Allele value)
{
	mutant_vectora::const_iterator it = begin;

	for (; it != end; ++it) {
		if (*it == value)
			return it;
	}
	return end;

}


}
#  endif

#endif
