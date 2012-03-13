#include "boost/numeric/ublas/vector_sparse.hpp"
using boost::numeric::ublas::compressed_vector;
using boost::numeric::ublas::sparse_vector_element;

#ifndef _MUTANT_VECTOR_H
#  define _MUTANT_VECTOR_H

namespace simuPOP {

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
				}else  {
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
				}else  {
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
typedef compressed_vector<Allele>::reference AlleleRef;
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
