/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__DESCENDER_HH_INCLUDED
#define CORE__DESCENDER_HH_INCLUDED

#ifndef CORE__CONFIGURATION_HH_INCLUDED
# include "configuration.hh"
#endif

namespace core
{

  class ARG;

  class Descender
  {
    public:
      Descender(const Configuration &conf) : i_conf(conf) {}
      ~Descender() {}

      // assign evolution to the ARG as specified by the configuration.
      // If the traits cannot be assigned according to specification, the
      // method returns `false' which means that a new ARG should be build
      // and processed.  If everything goes well, evolve returns `true'.
      void evolve(ARG &arg) const;

    private:
      const Configuration &i_conf;
  };

}
#endif
