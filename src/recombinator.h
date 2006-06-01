/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                     *
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

#ifndef _RECOMBINATOR_H
#define _RECOMBINATOR_H
/**
\file
\brief head file of class recombinator:public Operator
*/
#include "operator.h"

#include <iterator>
using std::ostream;
using std::ostream_iterator;

namespace simuPOP
{
	/**
	\brief  Recombination

	- only works for diploids (and for females in haplodiploids) population.

	- Free recombination between loci. Loci behave completely independently.

	- otherwise there will be some linkage between loci, user
	  need to specify physical recombination rate between adjacent loci
	(ie between locus n and n+1)

	- The recombination rate must be comprised between 0.0 and 0.5.

	- A recombination rate of 0.0 means that the loci are completely linked
	and thus behave together as a single linked locus.

	- A recombination rate of 0.5 is equivalent to free
	recombination.

	- All values in between will represent various linkage intensities
	between adjacent pairs of loci. The recombination rate is equivalent to
	1-linkage and represents the probability that the allele at the next locus is
	randomly drawn.

	*/

	class recombinator: public Operator
	{
		public:
			/// recombine chromosomes from parents
			/**

			\param intensity recombination rate per unit of loci distance. I.e., the really recombination
			rate between two loci is determined by intensity*loci distance between them.
			\param rate recombination rate regardless of loci distance; it can also be
			  an array of recombination rates. Must be the same length
			  as afterLoci or totNumOfLoci(). If totNumLociThe last item can be ignored.
			\param afterLoci an array of loci index. If rates is also specified, they
			should have the same length. Default to all loci (but meaningless for those
			loci locate at the end of chromosome.) If given, afterLoci should
			be ordered, and can not include loci at the end of a chromosome.
			\maleRate recombination rate for male individuals. If given,
			parameter rate will be considered as female rate.
			\maleIntensity recombination intensity for male individuals. If given,
			parameter intensity will be considered as female intensity.
			\maleAfterLoci if given, male will recombine at different locations.
			This is rarely used.
			\note rate gives recombination rate PER unit, rates gives plain rates.
			\note there is no recombination between sex chromosomes of male individuals.
			(sexChrom=True).
			*/
			recombinator(double intensity=-1,
				vectorf rate=vectorf(),
				vectoru afterLoci=vectoru(),
				double maleIntensity=-1,
				vectorf maleRate=vectorf(),
				vectoru maleAfterLoci=vectoru(),
				int method=1,
				int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL)
				:
			Operator( "", "", DuringMating, begin, end, step, at,
				rep, grp)
				, m_intensity(intensity), m_maleIntensity(maleIntensity),
				m_rate(rate), m_maleRate(maleRate),
				m_afterLoci(afterLoci), m_maleAfterLoci(maleAfterLoci),
				m_recBeforeLoci(0), m_maleRecBeforeLoci(0),
				m_bt(rng()), m_maleBt(rng()), m_recCount(0), m_method(method)
			{
				// tells mating schemes that this operator will form
				// the genotype of offspring so they do not have to
				// generate default genotype for offspring
				this->setFormOffGenotype(true);
			};

			virtual ~recombinator()
			{
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new recombinator(*this);
			}

			virtual string __repr__()
			{
				return "<simuPOP::recombination>" ;
			}

			/// this function takes intensity, rate, afterLoci, ...
			/// inputs and return a bernulli trailer and a recBeforeLoci
			/// vector.
			void prepareRecRates(population& pop,
				double intensity,
				vectorf rate,
				vectoru afterLoci,				  //
				bool sexChrom,					  // whether or not recombine the last chromosome
				vectoru& recBeforeLoci,			  // return before loci vector
				vectorf& vecP);					  // return recombination rate

			/// return recombination count
			ULONG recCount(size_t locus)
			{
				DBG_ASSERT( locus < m_recCount.size(), IndexError,
					"locus index " + toStr(locus) + " is out of range");
				return m_recCount[locus];
			}

			/// return recombination counts
			vectoru recCounts()
			{
				return m_recCount;
			}

			// this function implement how to recombine
			// parental chromosomes and set one copy of offspring chromosome
			// bt contains the bernulli trailer
			void recombine(
				individual* parent,				  // one of the parent
				population::IndIterator offspring,// offspring
				int offPloidy,					  // which offspring ploidy to fill
				BernulliTrials& bt,
				const vectoru& recBeforeLoci,
				bool setSex=false);

			virtual bool applyDuringMating(population& pop,
				population::IndIterator offspring,
				individual* dad=NULL,
				individual* mom=NULL);

		private:

			/// intensity
			double m_intensity, m_maleIntensity;

			/// differnt rates
			vectorf m_rate, m_maleRate;

			/// initial parameter
			vectoru m_afterLoci, m_maleAfterLoci;

			/// position to recombine, changed to fit a special pop
			vectoru m_recBeforeLoci, m_maleRecBeforeLoci;

			/// bernulli trials
			//  vector<BernulliTrials*> m_bt;
			BernulliTrials m_bt, m_maleBt;

			/// whether or not set sex (population having sex chromosome)
			bool m_setSex;

			/// report the number of recombination events
			vectoru m_recCount;

			/// which recombine method to use
			int m_method;
	};

}
#endif
