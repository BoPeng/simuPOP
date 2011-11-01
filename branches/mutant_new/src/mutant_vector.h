#include "boost/numeric/ublas/vector_sparse.hpp"
using boost::numeric::ublas::compressed_vector;

#ifndef _MUTANT_VECTOR_H
#define _MUTANT_VECTOR_H


namespace simuPOP {

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

		mutant_vector(size_t size, const T& value)
		{
			if (value == 0)
				m_container.resize(size, value);
			else {
				m_container.resize(size,size);
				for ( size_t i = 0; i < size; i++)
					m_container.push_back(i, value);
			}
		}

		void resize (size_t size)
		{
			m_container.resize(size);
		}

		void reserve (size_t size)
		{
			m_container.reserve(size);
		}

		size_t size() const {
			return m_container.size();
		}

		mutant_vector<T>::iterator  begin() 
		{
			return mutant_vector<T>::iterator(&m_container, 0); 
		}

		const mutant_vector<T>::iterator  begin() const
		{
			return mutant_vector<T>::iterator(const_cast<compressed_vector<T> *>(&m_container), 0); 
		}

		mutant_vector<T>::iterator  end() 
		{
			return mutant_vector<T>::iterator(&m_container, m_container.size()); 
		}
	
		const mutant_vector<T>::iterator  end() const
		{
			return mutant_vector<T>::iterator(const_cast<compressed_vector<T> *>(&m_container), m_container.size()); 
		}

		typename compressed_vector<T>::reference operator [] (size_t i)
		{
			return m_container[i];
		}

		void push_back (typename compressed_vector<T>::size_type i, typename compressed_vector<T>::const_reference t)
		{
			m_container.push_back(i, t);
		}

		void swap(mutant_vector<T> & vec)
		{
			m_container.swap(vec.getContainer());	

		}

			
		compressed_vector<T> &  getContainer()
		{
			return m_container;		
		}

		const compressed_vector<T> &  getContainer() const
		{
			return m_container;		
		}

		class iterator 
		{
			protected:
				compressed_vector<T> * m_container;
				mutable size_t m_index;
				
			public:
				
				iterator (compressed_vector<T> * c) : m_container(c), m_index(0)
				{
				}

				iterator (const compressed_vector<T> * c) : m_container(c), m_index(0)
				{
				}

				iterator (compressed_vector<T> * c, size_t index) : m_container(c), m_index(index)
				{
				}

				iterator (const compressed_vector<T> * c, size_t index) : m_container(c), m_index(index)
				{
				}

				iterator () : m_container(NULL), m_index(0)
				{
				}
/*
				iterator (const iterator & other)
				{
					std::cout << "call copy constructor: " << m_index << std::endl;
					m_container = other.m_container;
					m_index = other.m_index;

				}
*/

				iterator & operator= (const iterator & iter) 
				{
					m_container = iter.m_container;
					m_index = iter.m_index;
					return *this;
				}

				bool operator == (const iterator & iter) const 
				{
					if (m_index == iter.m_index)
						return true;
					return false;
				}

				bool operator != (const iterator & iter) const 
				{
					if (m_index != iter.m_index)
						return true;
					return false;
				}
				bool operator >= (const iterator & iter) const 
				{
					if (m_index >= iter.m_index)
						return true;
					return false;
				}

				bool operator <= (const iterator & iter) const 
				{
					if (m_index <= iter.m_index)
						return true;
					return false;
				}

				bool operator > (const iterator & iter) const 
				{
					if (m_index > iter.m_index)
						return true;
					return false;
				}

				bool operator < (const iterator & iter) const 
				{
					if (m_index < iter.m_index)
						return true;
					return false;
				}

				typename compressed_vector<T>::reference operator* () 
				{
					return (*m_container)[m_index];
				}

				typename compressed_vector<T>::const_reference operator* () const
				{
					return (*m_container)[m_index];
				}

				typename compressed_vector<T>::reference operator [] (const size_t i)
				{
					return (*m_container)[m_index + i];	
				}

				typename compressed_vector<T>::const_reference operator [] (const size_t i) const
				{
					return (*m_container)[m_index + i];	
				}

				/// CPPONLY pre-incrment return by-reference 
				iterator & operator++ () 
				{
					++m_index;
					return *this;
				}

				const iterator & operator++ () const
				{
					++m_index;
					return *this;
				}

				/// CPPONLY post-incrment return by-value 
				iterator operator++ ( int ) 
				{
					iterator orig = *this;
					++(*this);
					return orig;
				}

				const iterator operator++ ( int )  const
				{
					const iterator orig = *this;
					++(*this);
					return orig;
				}

				iterator & operator+= (const iterator & iter) 
				{
					m_index+=iter.m_index;
					if (m_index > m_container->size())
						m_index = m_container->size();
					return *this;
				}

				const iterator & operator+= (const iterator & iter)  const
				{
					m_index+=iter.m_index;
					if (m_index > m_container->size())
						m_index = m_container->size();
					return *this;
				}

				iterator & operator+= (const size_t size) 
				{
					m_index+=size;
					if (m_index > m_container->size())
						m_index = m_container->size();
					return *this;
				} 

				const iterator & operator+= (const size_t size)  const
				{
					m_index+=size;
					if (m_index > m_container->size())
						m_index = m_container->size();
					return *this;
				} 

				/*size_t operator + (iterator & iter)
				{
					return (*this).m_index + iter.m_index;
				}
				*/

				iterator operator + (const iterator & iter) const
				{
					iterator result = *this;
					result.m_index += iter.m_index;
					return result;
				}

				iterator operator + (const size_t size) const
				{
					iterator result = *this;
					result.m_index += size;
					return result;
				}

				/*	
				size_t operator - (iterator & iter)
				{
					return (*this).m_index - iter.m_index;
				}
				*/
				

				iterator operator - (const iterator & iter) const
				{
					iterator result = *this;
					result.m_index -= iter.m_index;
					return result;
				}

				iterator operator - (const size_t size) const
				{
					iterator result = *this;
					result.m_index -= size;
					return result;
				}

				size_t findPositionIndexData () {
					for (size_t i = 0; i < m_container->filled(); i++) {
						if (i < m_container->filled() - 1 && m_index > m_container->index_data()[i]) {
							continue;
						} else if (i == m_container->filled() - 1 &&  m_index > m_container->index_data()[i]) {
							return i + 1;		
						}
						else { 
							return i;
						}
					}
					return 0;
				}	

				typename compressed_vector<T>::index_array_type::iterator getIndexIterator () {
					return m_container->index_data().begin() + findPositionIndexData();
				}	

				typename compressed_vector<T>::value_array_type::iterator getValueIterator () {
					return m_container->value_data().begin() + findPositionIndexData();
				}	

				typename compressed_vector<T>::iterator getCompressedVectorIterator()
				{
					return m_container->find(m_index);
				}

				size_t getIndex()
				{
					return m_index;
				}

				compressed_vector<T> * getContainer()
				{
					return m_container;	
				}

				void deleted()
				{
					m_container->erase_element(m_index);
				}	

		};

		void erase (mutant_vector<T>::iterator begin,mutant_vector<T>::iterator end)
		{
			mutant_vector<T>::iterator it = begin;
			mutant_vector<T>::iterator it2 = end;
			for(;it != end; ++it) {
				it.deleted();
			}
			it = begin;
			if (end != this->end())
				for(it2 = end; it2 != this->end(); ++it2, ++it) {
					*it = *it2;
				}
			//  compressed vector doesn't decrease the size when it is erased.
			//  we have to resize it. 
			m_container.resize(m_container.size() - (end.getIndex() - begin.getIndex()));
		}

		void insert(mutant_vector<T>::iterator it, mutant_vector<T>::iterator begin, mutant_vector<T>::iterator end)
		{
			compressed_vector<size_t>::index_array_type::iterator it_src_begin = begin.getIndexIterator();
			compressed_vector<size_t>::index_array_type::iterator iend   = end.getIndexIterator();
			size_t src_begin = *it_src_begin;
			size_t src_index = *it_src_begin != begin.getIndex() ? *it_src_begin - begin.getIndex() : 0;	
			size_t dest_begin = it.getIndex();
			//size_t reserve_size = iend - it_src_begin;

			m_container.resize(m_container.size() + (end.getIndex() - begin.getIndex()));
			//m_container.reserve(m_container.nnz_capacity() + reserve_size);

			for (;it_src_begin  != iend; ++it_src_begin) {
				m_container.push_back(((*it_src_begin + src_index) - src_begin) + dest_begin, (*begin.getContainer())[*it_src_begin]);
			}

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

namespace simuPOP 
{

inline void copy(mutant_vectora::iterator begin, mutant_vectora::iterator end, mutant_vectora::iterator  it) 
{
	/*compressed_vector<Allele>::index_array_type::iterator it_src_begin = begin.getIndexIterator();
	compressed_vector<Allele>::index_array_type::iterator iend   = end.getIndexIterator();
	size_t src_begin = *it_src_begin;
	size_t src_index = *it_src_begin != begin.getIndex() ? *it_src_begin - begin.getIndex() : 0;	
	size_t dest_begin = it.getIndex();
	for (;it_src_begin  != iend; ++it_src_begin) {
		(*it.getContainer())[ ((*it_src_begin + src_index) - src_begin) + dest_begin] = (*begin.getContainer())[*it_src_begin];
	}*/
	mutant_vectora::iterator itt = begin;
	for (;itt != end; ++itt, ++it)
		*it = *itt;	
}

inline void fill(mutant_vectora::iterator begin, mutant_vectora::iterator end, Allele value) 
{
	if ( value == 0 )	
        	for (mutant_vectora::iterator it = begin; it != end; ++it)
			it.deleted();
	else
        	for (mutant_vectora::iterator it = begin; it != end; ++it)
			*it = value;



}

}
#  endif

#endif
