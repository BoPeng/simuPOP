#ifndef _MUTANT_VECTOR_H
#define _MUTANT_VECTOR_H

#include "boost/numeric/ublas/vector_sparse.hpp"
using boost::numeric::ublas::compressed_vector;

namespace simuPOP {

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
			public:
				iterator (compressed_vector<T> * c) : m_container(c), m_index(0)
				{
				}

				iterator (compressed_vector<T> * c, size_t index) : m_container(c), m_index(index) 
				{
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


				void setIndex(size_t index)
				{
					m_index = index;
				}
				void setContainer(compressed_vector<T> *  c)
				{
					m_container = c;	
				}
		};

	private:
		compressed_vector<T> m_container;
};

}

#endif
