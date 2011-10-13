#ifndef _MUTANT_VECTOR_H
#define _MUTANT_VECTOR_H

#include "boost/numeric/ublas/vector_sparse.hpp"
using boost::numeric::ublas::compressed_vector;

namespace simuPOP {

#  ifdef MUTANTALLELE
template <class T>
class mutant_vector
{
	public:

		class iterator;

		mutant_vector()
		{
		}

		mutant_vector(size_t size)
		{
			m_container.resize(size);
		}

		void resize (size_t size)
		{
			m_container.resize(size);
		}

		size_t size() {
			return m_container.size();
		}

		mutant_vector<T>::iterator  begin()
		{
			return mutant_vector<T>::iterator(&m_container, 0); 
		}

		mutant_vector<T>::iterator  end()
		{
			return mutant_vector<T>::iterator(&m_container, m_container.size()); 
		}

		void swap(mutant_vector<T> & vec)
		{
			m_container.swap(vec.getContainer());	

		}

		compressed_vector<T> &  getContainer()
		{
			return m_container;		
		}

		class iterator 
		{
			protected:
				compressed_vector<T> * m_container;
				size_t m_index;
				typename compressed_vector<T>::iterator m_iter;	
				
			public:
				iterator (compressed_vector<T> * c) : m_container(c), m_index(0), m_iter(m_container->begin())
				{
				}

				iterator (compressed_vector<T> * c, size_t index) : m_container(c), m_index(index)
				{
					if (index == 0) { 
						m_iter = m_container->begin();
					}
					else if (index == m_container->size()) {
						m_iter = m_container->end();
					}
				}

				iterator () : m_container(NULL), m_index(0)
				{
				}

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


				typename compressed_vector<T>::reference operator* () 
				{
					return (*m_container)[m_index];
				}

				/// CPPONLY pre-incrment return by-reference 
				iterator & operator++ () 
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

				iterator & operator+= (const iterator & iter) 
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

				size_t operator + (const iterator & iter) 
				{
					return (*this).m_index + iter.m_index;
				}

				iterator operator + (const size_t size) const
				{
					iterator result = *this;
					result.m_index += size;
					return result;
				}

				size_t operator - (const iterator & iter) const
				{
					return (*this).m_index - iter.m_index;
				}

				iterator operator - (const size_t size) 
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
		};

	private:
		compressed_vector<T> m_container;
};

}

typedef unsigned int Allele;
typedef unsigned int & AlleleRef;
typedef simuPOP::mutant_vector<Allele> vectora;

namespace simuPOP 
{

inline void copy(vectora::iterator begin, vectora::iterator end, vectora::iterator  it) 
{
	compressed_vector<Allele>::index_array_type::iterator it_src_begin = begin.getIndexIterator();
	compressed_vector<Allele>::index_array_type::iterator iend   = end.getIndexIterator();
	compressed_vector<Allele>::index_array_type::iterator it_dest_begin = it.getIndexIterator();
	size_t src_begin = *it_src_begin;
	size_t src_index = *it_src_begin != begin.getIndex() ? *it_src_begin - begin.getIndex() : 0;	
	size_t dest_begin = *it_dest_begin;
	for (;it_src_begin  != iend; ++it_src_begin) {
		(*it.getContainer())[ ((*it_src_begin + src_index) - src_begin) + dest_begin ] = (*begin.getContainer())[*it_src_begin];
	}
}

}

#  endif
#endif
