/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                    *
 *                                                                         *
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

#ifndef _OUTPUTER_H
#define _OUTPUTER_H
/**
\file
\brief head file of class outputer: public baseOperator
*/
#include "utility.h"
#include "operator.h"
#include <iostream>

#include <iomanip>
using std::setw;
using std::hex;
using std::dec;

namespace simuPOP
{
	/**
	\brief Base class of all operators that out information.
	different format.

	@author Bo Peng
	*/

	class outputer: public baseOperator
	{

		public:
			/// constructor. 
			outputer(string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			baseOperator(output, outputExpr, stage, begin, end, step, at, rep, grp, infoFields)
			{
			};

			/// destructor
			virtual ~outputer(){};

			virtual baseOperator * clone() const
			{
				return new outputer(*this);
			}

	};

	/// Output a given string.
	/**
	A common usage is <tt>pyOutpue('\n', rep=REP_LAST)</tt>
	*/
	class pyOutput: public outputer
	{

		public:
			/// Create a pyOutput operator that output a given string.
			/**
				\param str string to be outputted
			*/
			pyOutput(string str="", string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			outputer( output, outputExpr, stage, begin, end,
				step, at, rep, grp, infoFields), m_string(str)
			{
			}

			/// simply output some info
			virtual bool apply(population& pop)
			{
				ostream& out = this->getOstream(pop.dict());
				out << m_string;
				this->closeOstream();
				return true;
			}

			///// destructor
			virtual ~pyOutput()
			{
			}

			virtual baseOperator * clone() const
			{
				return new pyOutput(*this);
			}

			/// set output string.
			void setString(const string str)
			{
				m_string = str;
			}

			virtual string __repr__()
			{
				string reprStr;
				for(size_t i=0; i<10 && i < m_string.size(); ++i)
					if(m_string[i] != '\n')
						reprStr += m_string[i];
				if( m_string.size() > 10)
					reprStr += "... ";
				return "<simuPOP::output " + reprStr + "> " ;
			}

		private:
			string m_string;
	};

	/// Output a given string.
	/**
	A common usage is <tt>pyOutpue('\n', rep=REP_LAST)</tt>

	This operator, renamed to output, in the python interface
	is obsolete. It is replaced by pyOutput.
	*/
	class outputHelper: public outputer
	{

		public:
			/// Create a outputHelper operator that output a given string.
			/**
				\param str string to be outputted
			*/
			outputHelper(string str="", string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			outputer( output, outputExpr, stage, begin, end,
				step, at, rep, grp, infoFields), m_string(str)
			{
			}

			/// simply output some info
			virtual bool apply(population& pop)
			{
				ostream& out = this->getOstream(pop.dict());
				out << m_string;
				this->closeOstream();
				return true;
			}

			///// destructor
			virtual ~outputHelper()
			{
			}

			virtual baseOperator * clone() const
			{
				return new outputHelper(*this);
			}

			/// set output string.
			void setString(const string str)
			{
				m_string = str;
			}

			virtual string __repr__()
			{
				string reprStr;
				for(size_t i=0; i<10 && i < m_string.size(); ++i)
					if(m_string[i] != '\n')
						reprStr += m_string[i];
				if( m_string.size() > 10)
					reprStr += "... ";
				return "<simuPOP::output " + reprStr + "> " ;
			}

		private:
			string m_string;
	};

	/// dump the content of a population.

	class dumper: public outputer
	{
		public:
			/** \brief dump population

			\param alleleOnly only display allele
			\param infoOnly only display info
			\param dispWidth width of allele display, default to 1
			\param ancestralPops whether or not display ancestral populations, default to False
			\param chrom chromsoome(s) to display
			\param loci loci to display
			\param subPop only display subPop(s)
			\param indRange range(s) of individuals to display
			\param max max number of individuals to display, default to 100.
			This is to avoid careless dump of huge populations.
			\param output output file, default to standard output.
			\param outputExpr and other parameters: refer to help(baseOperator.__init__)

			*/
			dumper( bool alleleOnly=false, bool infoOnly=false, bool ancestralPops=false, int dispWidth=1, UINT max=100,
				const vectori& chrom=vectori(), const vectori& loci=vectori(), const vectoru& subPop=vectoru(),
				const vectorlu& indRange=vectorlu(),
				string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			outputer(output, outputExpr, stage, begin, end, step, at, rep, grp, infoFields),
				m_alleleOnly(alleleOnly), m_infoOnly(infoOnly), m_dispAncestry(ancestralPops), m_width(dispWidth),
				m_chrom(chrom), m_loci(loci), m_subPop(subPop), m_indRange(indRange), m_max(max)
				{}

			virtual baseOperator * clone() const
			{
				return new dumper(*this);
			}

			/// only show alleles (not structure, gene information?
			bool alleleOnly()
			{
				return m_alleleOnly;
			}

			///
			void setAlleleOnly(bool alleleOnly)
			{
				m_alleleOnly = alleleOnly;
			}

			/// only show info
			bool infoOnly()
			{
				return m_infoOnly;
			}

			///
			void setInfoOnly(bool infoOnly)
			{
				m_infoOnly = infoOnly;
			}

			virtual bool apply(population& pop);

			virtual ~dumper(){};

			virtual string __repr__()
			{
				return "<simuPOP::dumper>" ;
			}

		private:
			/// only output alleles, not structure info
			bool m_alleleOnly;

			/// only display info
			bool m_infoOnly;

			/// whether or not display ancestral populations
			bool m_dispAncestry;

			/// disp width when outputing alleles
			int m_width;

			///
			vectori m_chrom;

			///
			vectori m_loci;

			///
			vectoru m_subPop;

			///
			vectorlu m_indRange;

			/// only output first ... individuals. Good for large population
			UINT m_max;
	};

	/// save population to a file

	class savePopulation: public outputer
	{
		public:
			savePopulation( string output="", string outputExpr="",
				string format = "bin", bool compress=true, int stage=PostMating, int begin=0, int end=-1,
				int step=1, vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			outputer( "", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_filename(output), m_filenameParser(outputExpr), m_format(format), m_compress(compress)
			{
				if(output == "" && outputExpr == "")
					throw ValueError("Please specify one of output and outputExpr.");
			}

			~savePopulation()
			{
			}

			virtual baseOperator * clone() const
			{
				return new savePopulation(*this);
			}

			virtual bool apply(population& pop);

			virtual string __repr__()
			{
				return "<simuPOP::save population>" ;
			}

		private:
			/// filename,
			string m_filename;

			/// or an expression that will be evaluated dynamically
			Expression m_filenameParser;

			/// can specify format, default to 'auto'
			string m_format;

			/// whether or not compress population
			bool m_compress;
	};

}
#endif
