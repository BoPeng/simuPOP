#!/usr/bin/env python

#
# $File: simuUtil.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit http://simupop.sourceforge.net
# for details.
#
# Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#


"""
simuPOP utilities.

This module provides some commonly used operators
and format conversion utilities.

"""

import exceptions, operator, types, os, sys, re

from simuPOP import *
from simuOpt import simuOptions

def ViewVars(var, gui=None):
    '''
    list a variable in tree format, either in text format or in a
        wxPython window.

    var
        A dictionary variable to be viewed. Dictionary wrapper objects returned
        by ``population.dvars()`` and ``simulator.dvars()`` are also acceptable.

    gui
        If gui is ``False`` or ``'Tkinter'``, a text presentation (use the
        pprint module) of the variable will be printed to the screen. If gui is
        ``'wxPython'`` and wxPython is available, a wxPython windows will be
        used. The default mode is determined by the global gui mode (see also
        ``simuOpt.setOptions``).
    '''
    if gui is None:
        gui = simuOptions['GUI']
    #
    if gui == False or gui == 'Tkinter':
        import pprint
        pprint.pprint(var)
        return

    try:
        import wx, wx.py.filling as fill
    except ImportError:
        pprint.pprint(var)
        return

    app = wx.PySimpleApp()
    wx.InitAllImageHandlers()
    if var==None:
        fillFrame = fill.FillingFrame()
    else:
        if type(var) == type( dw({}) ):
            fillFrame = fill.FillingFrame(rootObject=var.__dict__,
                rootLabel='var')
        else:
            fillFrame = fill.FillingFrame(rootObject=var,
                rootLabel='var')
    fillFrame.Show(True)
    app.MainLoop()


# migration rate matrix generators
def MigrIslandRates(r, n):
    '''migration rate matrix

::

         x m/(n-1) m/(n-1) ....
         m/(n-1) x ............
         .....
         .... m/(n-1) m/(n-1) x

    where x = 1-m
    '''
    # n==1?
    if n == 1:
        return [[1]]
    #
    m = []
    for i in range(0,n):
        m.append([r/(n-1.)]*n)
        m[-1][i] = 1-r
    return m


def MigrHierarchicalIslandRates(r1, r2, n):
    '''
    Return the migration rate matrix for a hierarchical island model
    where there are different migration rate within and across groups
    of islands.

    r1
        Within group migration rates. It can be a number or a list of numbers
        for each group of the islands.

    r2
        Across group migration rates which is the probability that someone will
        migrate to a subpopulation outside of his group. A list of r2 could be
        specified for each group of the islands.

    n
        Number of islands in each group. E.g. n=[5, 4] specifies two groups of
        islands with 5 and 4 islands each.

    For individuals in an island, the probability that it remains in the same
    island is 1-r1-r2 (r1, r2 might vary by island groups), that it migrates
    to another island in the same group is r1 and migrates to another island
    outside of the group is r2. Migrate rate to a specific island depends on
    the size of group.
    '''
    if type(n) not in [type(()), type([])]:
        raise exceptions.ValueError('A list of size of island groups is expected for parameter n')
    nIslands = sum(n)
    if type(r1) in [type(0), type(1.)]:
        r1 = [r1] * len(n)
    elif len(r1) != len(n):
        raise exceptions.ValueError('If multiple r1 is given, it should be given to all island groups.')
    #
    if type(r2) in [type(0), type(1.)]:
        r2 = [r2] * len(n)
    elif len(r2) != len(n):
        raise exceptions.ValueError('If multiple r2 is given, it should be given to all island groups.')
    #
    m = []
    for groupIdx, groupSize in enumerate(n):
        nOther = nIslands - groupSize
        groupStart = sum(n[:groupIdx])
        groupEnd = groupStart + groupSize
        for island in range(groupStart, groupEnd):
            m.append([])
            for i in range(groupStart):
                m[-1].append(r2[groupIdx] * 1.0 / nOther)
            for i in range(groupStart, groupEnd):
                if i == island:
                    m[-1].append(1 - r1[groupIdx] - r2[groupIdx])
                else:
                    m[-1].append(r1[groupIdx] * 1.0 / groupSize)
            for i in range(groupEnd, nIslands):
                m[-1].append(r2[groupIdx] * 1.0 / nOther)
    return m


def MigrSteppingStoneRates(r, n, circular=False):
    '''migration rate matrix, circular stepping stone model (X=1-m)

::

           X   m/2               m/2
           m/2 X   m/2           0
           0   m/2 x   m/2 ......0
           ...
           m/2 0 ....       m/2  X

or non-circular

::

           X   m/2               m/2
           m/2 X   m/2           0
           0   m/2 X   m/2 ......0
           ...
           ...              m   X
    '''
    if n < 2:
        raise exceptions.ValueError("Can not define step stone model for n < 2")
    elif n == 2:
        return [[1-r,r],[r,1-r]]
    # the normal case (n>2)
    m = []
    for i in range(0, n):
        m.append([0]*n)
        m[i][i] = 1-r
        m[i][(i+1)%n] = r/2.
        m[i][(i+n-1)%n] = r/2.
    if not circular:
        m[0][1] = r
        m[0][-1] = 0
        m[n-1][0] = 0
        m[n-1][n-2] = r
    return m


def SaveCSV(pop, filename='', fields=[], loci=[], header=True,
    shift=1, combine=None,
        sexCode={Male: '1', Female: '2'}, affectionCode={True: '1', False: '2'},
        **kwargs):
    '''Save a simuPOP population ``pop`` in csv format.

    pop
        A simuPOP population object. If a string is given, it will be loaded.

    filename
        Output filename.

    fileds
        Information fields to be outputted.

    loci
        If a list of loci is given, only genotype at these loci will be
        written.

    header
        Whether or not a header should be written. These headers will include
        information fields, sex (if ``sexCode`` is not ``None``), affection
        status (if ``affectionCode`` is not ``None``) and loci names. If
        genotype at a locus needs more than one column, ``_1``, ``_2`` etc will
        be appended to locus names.

    genotype
        list of loci to output, default to all.

    combine
        how to combine the markers. Default to None.
        A function can be specified, that takes the form::

             def func(markers):
                 return markers[0]+markers[1]

    shift
        since alleles in simuPOP is 0-based, shift=1 is usually needed to
        output alleles starting from allele 1. This parameter is ignored if
        combine is used.

    '''
    if loci == []:
        loci = range(0, pop.totNumLoci())
    try:
        out = open(filename, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + filename +" to write."
    # keep the content of pieces in strings first
    content = [''] * pop.numChrom()
    # write out header
    print >> out, 'id, ', ', '.join(fields), ', ',
    if combine is None:
        print >> out, ', '.join(['%s_1, %s_2' % (pop.locusName(loc), pop.locusName(loc)) for loc in loci])
    else:
        print >> out, ', '.join(['%s' % pop.locusName(loc) for loc in loci])
    # write out
    id = 1
    pldy = pop.ploidy()
    for ind in pop.individuals():
        print >> out, id,
        for f in fields:
            if f == 'sex':
                print >> out, ',', sexCode[ind.sex()],
            elif f == 'affection':
                print >> out, ',', affectionCode[ind.affected()],
            else:
                print >> out, ',', ind.info(f),
        for marker in loci:
            if combine is None:
                for p in range(pldy):
                    out.write(", %d" % (ind.allele(marker, p) + shift))
            else:
                out.write(", %d" % combine([ind.allele(marker, p) for p in range(pldy)]))
        print >> out
        id += 1
    out.close()


class _baseProgressBar:
    def __init__(self, message, totalCount):
        '''
        message
            Title of the progress bar

        totalCount
            Total expected steps.

        done
            Message displayed when the job is finished.
        '''
        self.message = message
        self.totalCount = totalCount
        self.percent = 0
        self.completed = False

    def update(self, count):
        '''
        Update the progreebar.
        '''
        count = min(count, self.totalCount)
        self.progress = int(round(100*count/self.totalCount))
        if self.progress <= self.percent:
            return False
        else:
            return True

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if self.completed:
            return False
        else:
            self.completed = True
            return True

class _textProgressBar(_baseProgressBar):
    def __init__(self, message, totalCount, progressChar='.', block=2, done=' Done.\n'):
        '''
        message
            Title of the progress bar

        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        _baseProgressBar.__init__(self, message, totalCount)
        self.percent = 0
        self.progressChar = progressChar
        self.block = block
        self.doneMsg = done
        sys.stdout.write(message)
        sys.stdout.flush()

    def update(self, count):
        '''
        Update the text progreebar.
        '''
        if not _baseProgressBar.update(self, count):
            return
        for p in range(self.percent + 1, self.progress + 1):
            if p == 100:
                self.done()
            elif p % 10 == 0:
                sys.stdout.write(str(p/10))
            elif p % self.block == 0:
                sys.stdout.write(self.progressChar)
        sys.stdout.flush()
        self.percent = self.progress
        if self.percent == 100:
            self.done()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if not _baseProgressBar.done(self):
            return
        sys.stdout.write(self.doneMsg)
        sys.stdout.flush()


class _tkProgressBar(_baseProgressBar):
    def __init__(self, message, totalCount):
        '''
        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        _baseProgressBar.__init__(self, message, totalCount)
        import Tkinter as tk
        self.width = 300
        self.height = 30
        self.max = 100
        self.fillColor = 'blue'
        self.labelColor = 'black'
        self.label = 'Progress'
        #
        self.app = tk.Tk()
        self.app.title(self.label)
        self.frame = tk.Frame(self.app, bd=0)
        self.canvas = tk.Canvas(self.frame, bd=0, width=self.width+40,
            height = self.height + 70, highlightthickness=0)
        self.label = self.canvas.create_text(20, 20, 
            text='', anchor="w", fill=self.labelColor, font=('Verdana', 10))
        self.scale = self.canvas.create_rectangle(
            20, 50, self.width + 20, 50 + self.height, fill=self.fillColor)
        self.rect = self.canvas.create_rectangle(
            20, 50, self.width + 20, 50 + self.height)
        self.canvas.pack(side='top', fill='x', expand='yes', padx=0)
        self.update(0)
        self.frame.pack(padx=0, pady=0)

    def update(self, count):
        '''
        Update the progreebar.
        '''
        if not _baseProgressBar.update(self, count):
            return
        #
        self.canvas.coords(self.scale, 20, 50,
            20 + self.progress * 1.0 / self.max * self.width, 50 + self.height)
        # Now update the colors
        self.canvas.itemconfig(self.scale, fill=self.fillColor)
        self.canvas.itemconfig(self.label, fill=self.labelColor)
        # And update the label
        if self.progress > 0:
            self.canvas.itemconfig(self.label, text=self.message + "\n%d%% completed." % self.progress)
        else:
            self.canvas.itemconfig(self.label, text=self.message)
        self.canvas.update_idletasks()
        self.app.update()
        #
        self.percent = self.progress
        if self.percent == 100:
            self.done()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if not _baseProgressBar.done(self):
            return
        self.app.destroy()
        del self.app


class _wxProgressBar(_baseProgressBar):
    def __init__(self, message, totalCount):
        '''
        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        _baseProgressBar.__init__(self, message, totalCount)
        import wx
        self.app = wx.PySimpleApp(0)
        self.dialog = wx.ProgressDialog(
            'Progress', self.message + '\n', self.totalCount,
            style = \
                # wx.PD_CAN_ABORT | \
                # wx.PD_CAN_SKIP | \
                wx.PD_ELAPSED_TIME | \
                # wx.PD_ESTIMATED_TIME | \
                wx.PD_AUTO_HIDE | \
                wx.PD_REMAINING_TIME
            )
        self.dialog.Update(0)

    def update(self, count):
        '''
        Update the progreebar.
        '''
        if not _baseProgressBar.update(self, count):
            return
        self.dialog.Update(count, self.message + "\n%d%% completed." % self.progress)
        self.percent = self.progress
        if self.percent == 100:
            self.done()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if not _baseProgressBar.done(self):
            return
        self.dialog.Destroy()
        del self.app



class simuProgress:
    '''The ``simuProgress`` class defines a progress bar. This class will use a
    text-based progress bar that outputs progressing dots (.) with intermediate
    numbers (e.g. 5 for 50%) under a non-GUI mode (``gui=False``). In the GUI
    mode, a Tkinter or wxPython progress dialog will be used (``gui=Tkinter``
    or ``gui=wxPython``). The default mode is determined by the global gui mode
    of simuPOP (see also ``simuOpt.setOptions``).

    This class is usually used as follows::

        progress = simuProgress("Start simulation", 500)
        for i in range(500):
            progress.update(i+1)
        # if you would like to make sure the done message is displayed.
        progress.done()
    '''
    def __init__(self, message, totalCount, progressChar='.', block=2, done=' Done.\n', gui=None):
        '''
        message
            Title of the progress bar.

        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress. This is only used for
            text-based progress bars.

        block
            Intervals at which progresses will be displayed. Default to 2 (percent).

        done
            Message displayed when the job is finished.
        '''
        if gui is None:
            self.gui = simuOptions['GUI']
        else:
            self.gui = gui
        if self.gui in ['wxPython', True]:
            try:
                import wx
                self.gui = 'wxPython'
            except ImportError:
                self.gui = 'Tkinter'
        if self.gui == 'Tkinter':
            try:
                import Tkinter
            except ImportError:
                self.gui = False
        if self.gui == 'wxPython':
            self.progressBar = _wxProgressBar(message, totalCount)
        elif self.gui == 'Tkinter':
            self.progressBar = _tkProgressBar(message, totalCount)
        else:
            self.progressBar = _textProgressBar(message, totalCount, progressChar, block, done)

    def update(self, count):
        '''
        Update the progreebar.
        '''
        self.progressBar.update(count)

    def done(self):
        '''
        Finish progressbar, print 'done' message if in text-mode.
        '''
        self.progressBar.done()


class trajectory:
    '''A ``trajectory`` object contains frequencies of one or more loci in one
    or more subpopulations over several generations. It is usually returned by
    member functions of class ``trajectorySimulator`` or equivalent global
    functions ``ForwardTrajectory`` and ``BackwardTrajectory``.

    The ``trajectory`` object provides several member functions to facilitate
    the use of trajectory-simulation techiniques. For example,
    ``trajectory.func()`` returns a trajectory function that can be provided
    directly to a ``controlledOffspringGenerator``; ``trajectory.mutators()``
    provides a list of ``pointMutator`` that insert mutants at the right
    generations to initialize a trajectory.

    For more information about trajectory simulation techniques and related
    controlled random mating scheme, please refer to the simuPOP user's guide,
    and Peng et al (PLoS Genetics 3(3), 2007).
    '''
    def __init__(self, endGen, nLoci):
        '''Create a ``trajectory`` object of alleles at *nLoci* loci with
        ending generation *endGen*. *endGen* is the generation when expected
        allele frequencies are reached after mating. Therefore, a trajectory
        for 1000 generations should have ``endGen=999``.
        '''
        self.traj = {}
        self.endGen = endGen
        self.nLoci = nLoci
       
    def freq(self, gen):
        '''Return frequencies of all loci at generation *gen* of the simulated
        trajectory. Allele frequencies are assumed to be zero if *gen* is out
        of range of the simulated trajectory.
        '''
        if not self.traj.has_key(gen):
            return [0.] * self.nLoci
        return self.traj[gen]

    def func(self):
        '''Return a Python function that returns allele frequencies for each
        locus at specified loci. If there are multiple subpopulations, allele
        frequencies are arranged in the order of ``loc0_sp0``, ``loc0_sp1``,
        ..., ``loc1_sp0``, ``loc1_sp1``, ... and so on. The returned function
        can be supplied directly to the ``freqFunc`` parameter of a controlled
        random mating scheme (``controlledRandomMating``) or a homogeneous
        mating scheme that uses a controlled offspring generator
        (``controlledOffspringGenerator``).
        '''
        def trajFunc(gen):
            return self.freq(gen)
        return trajFunc
    
    def mutators(self, inds=0, allele=1, *args, **kwargs):
        '''Return a list of ``pointMutator`` operators that introduce mutants
        at the beginning of simulated trajectories. These mutators should be
        added to the ``ops`` parameter of ``simulator.evolve`` function to
        introduce a mutant at the beginning of a generation with zero allele
        frequency before mating, and a positive allele frequency after mating.
        Other than default parameters ``inds=0`` and ``allele=1``, additional
        parameters could be passed to point mutator as keyward parameters.
        '''
        gens = self.traj.keys()
        gens.sort()
        if len(gens) == 0:
            return []
        mut = []
        for gen in gens[:-1]:
            # no introduction of mutants with population merge or split.
            if len(self.traj[gen]) != len(self.traj[gen + 1]):
                continue
            # we may need to introduce mutant at each subpopulation.
            nSP = len(self.traj[gen]) / self.nLoci
            for loc in range(self.nLoci):
                for sp in range(nSP):
                    if (self.traj[gen][sp + loc * nSP] == 0
                        and self.traj[gen + 1][sp + loc * nSP] > 0):
                        mut.append(pointMutator(inds=inds, loci=loc, allele=allele,
                            subPops=sp, at=gen + 1, stage=PreMating, *args, **kwargs))
        return mut

    def _setFreq(self, freq, gen, nSubPop):
        '''This function sets frequency *freq* at specified generation *gen*.
        *nSubPop* is used if frequency for multiple subpopulations are given.
        '''
        if type(freq) in [type(0), type(0.)]:
            self.traj[gen] = [freq] * self.nLoci * nSubPop
        elif len(freq) == 1:
            self.traj[gen] = freq * self.nLoci * nSubPop
        elif len(freq) == self.nLoci:
            f = []
            for idx, ind in enumerate(freq):
                for sp in range(nSubPop):
                    f.append(ind)
            self.traj[gen] = f
        elif len(freq) == self.nLoci * nSubPop:
            self.traj[gen] = freq
        else:
            raise exceptions.ValueError('please define freq as a single number or a list with specified number of elements.')

    def plot(self, filename=None, **kwargs):
        '''Plot simulated trajectory using ``R`` through a Pythonmodule ``rpy``.
        The function will return silently if module ``simuRPy`` cannot be
        imported.

        This function will use different colors to plot trajectories at
        different loci. The trajectories are plotted from generation 0 to
        *endGen* even if the trajectories are short. The y-axis ranges from 0
        to 1 and is labeled ``Allele frequency``. If a valid *filename* is
        given, the figure will be saved to *filename* in a format specified by
        file extension. Currently supported formats/extensions are eps, jpg,
        bmp, tif, png and pdf. The availability of formats may be limited by
        your version of R.

        This function makes use of the derived keyword parameter feature of
        module ``simuRPy``. Allowed prefixes are ``par``, ``plot``, ``lines``
        and ``dev_print``. Allowed repeating suffix is ``loc``. For example,
        you could use parameter ``plot_ylim`` to reset the default value of
        ``ylim`` in R function ``plot``.
        '''
        try:
            import simuRPy
        except ImportError:
            return
        #
        args = simuRPy.derivedArgs(
            defaultFuncs = ['plot', 'lines'],
            allFuncs = ['par', 'plot', 'lines', 'dev_print'],
            suffixes = ['loc'],
            defaultParams = {
                'plot_ylim': [0, 1],
                'plot_ylab': 'Allele frequency',
                'plot_xlab': 'Generation'
            },
            **kwargs
        )
        # device
        simuRPy.newDevice()
        #
        gens = self.traj.keys()
        simuRPy.rpy.r.par(**args.getArgs('par', None))
        simuRPy.rpy.r.plot(gens[0], 0,
            **args.getArgs('plot', None, xlim=(min(gens), max(gens)), type='n'))
        lines = []
        for gen in gens:
            if len(lines) != len(self.traj[gen]):
                if len(lines) != 0:
                    nSP = len(lines) / self.nLoci
                    for idx, line in enumerate(lines):
                        # remove leading zero's
                        for start in range(len(line)):
                            if line[start] != 0:
                                break
                        if start == len(line) - 1:
                            continue
                        r.lines(x=range(genBegin - len(line) + start, gen), y=line[start:],
                                col = int(idx / nSP) + 1)
                # the first point
                lines = [[x] for x in self.traj[gen]]
            else:
                for idx, x in enumerate(self.traj[gen]):
                    lines[idx].append(x)
        # plot the lines again
        for idx, line in enumerate(lines):
            nSP = len(lines) / self.nLoci
            # remove leading zero's
            for start in range(len(line)):
                if line[start] != 0:
                    break
            if start == len(line) - 1:
                continue
            genBegin = gen - len(line) + start
            simuRPy.rpy.r.lines(x=range(genBegin, gen), y=line[start:],
                    **args.getArgs('lines', None, col = int(idx / nSP) + 1))
        #
        simuRPy.saveFigure(**args.getArgs('dev_print', None, file=filename))
        

class trajectorySimulator:
    '''A trajectory simulator takes basic demographic and genetic (natural
    selection) information of an evolutionary process of a diploid population
    and allow the simulation of trajectory of allele frequencies of one or
    more loci. Trajectories could be simulated in two ways: forward-time and
    backward-time. In a forward-time simulation, the simulation starts from
    certain allele frequency and simulate the frequency at the next generation
    using given demographic and genetic information. The simulation continues
    until an ending generation is reached. A trajectory is successfully
    simulated if the allele frequency at the ending generation falls into a
    specified range. In a backward-time simulation, the simulation starts from
    the ending generation with a desired allele frequency and simulate the
    allele frequency at previous generations one by one until the allele gets
    lost (allele frequency equals zero).

    The result of a trajectory simulation is a trajectory object which can be
    used to direct the simulation of a special random mating process that
    controls the evolution of one or more disease alleles so that allele 
    frequencies are consistent across replicate simulations. For more
    information about trajectory simulation techniques and related controlled
    random mating scheme, please refer to the simuPOP user's guide, and Peng et
    al (PLoS Genetics 3(3), 2007).
    '''

    def __init__(self, N, nLoci=1, fitness=None, logger=None):
        '''Create a trajectory simulator using provided demographic and genetic
        (natural selection) parameters. Member functions *simuForward* and
        *simuBackward* can then be used to simulate trajectories within certain
        range of generations. Parameter *N* can be a constant or a demographic
        function that returns the population size (or a list of subpopulation
        sizes) at each generation. Parameter *fitness* can be a list of fitness
        values for genotype *AA*, *Aa*, and *aa, respectively; or a function
        that returns fitness values at each generation. When multiple loci are
        involved (parameter *nLoci*), *fitness* can be a list of 3 (the same
        fitness values for all loci), a list of 3*nLoci (different fitness
        values for each locus) or a list of 3**nLoci (fitness value for each
        combination of genotype). This simulator by default does not produce
        any output. However, if a valid logger (see Python module ``logging``)
        is given, detailed information about how trajectories are simulated
        will be written to this logger object.
        '''
        # a vector of subpopulation sizes is needed
        if type(N) in [type(1), type(1L)]:
            self.N = [N]
        else: # N is a list or a function
            self.N = N
        if fitness is None:
            self.fitness = [1, 1, 1]
        else:
            if type(fitness) in [type(()), type([])] and len(fitness) not in [3, 3*nLoci, 3**nLoci]:
                raise exceptions.ValueError('Invalid list of fitness.')
            self.fitness = fitness
        self.logger = logger
        self.nLoci = nLoci
        self.errorCount = {'tooLong' : 0, 'tooShort': 0, 'invalid': 0}
        self.maxMutAge = 0
        self.minMutAge = 0
    
    def _Nt(self, gen):
        'Get Nt(gen) depending on the type of N'
        # _Nt() expects parameter gen
        if callable(self.N):
            nt = self.N(gen)
            if type(nt) in [type(1), type(1L)]:
                return [nt]
            else:
                return nt
        else:
            return self.N

    def _marginalFitness(self, fitness, freq):
        '''Convert interaction fitness (3**n elements) to marginal fitness
        (3*n elements) using given allele frequency. The marginal fitnesses
        are calculated using formula:
                f(X=Aa) = Sum_g P(Y=g) * f(X=Aa, Y=g)
        where g is genotype at all other loci.
        '''
        sAll = [0] * (3 * self.nLoci)
        # each locus
        for loc in range(self.nLoci):
            # each genotype AA, Aa and aa (geno is the number of disease allele)
            for geno in range(3):
                # iterate through OTHER DSL
                allgeno = [0] * self.nLoci
                # set myself
                allgeno[loc] = geno
                # iterate through genotype at other loci
                f = 0.
                for it in range(3**(self.nLoci - 1)):
                    # assign allgeno, using it as a 3-based integer.
                    num = it
                    for l in range(self.nLoci):
                        if l == loc:
                            continue
                        allgeno[l] = num % 3
                        num /= 3
                    # calculate P(Y=g) and f(X=Aa, Y=g)
                    index = 0
                    fq = 1.
                    for i in range(len(allgeno)):
                        if i != loc:
                            if allgeno[i] == 0:
                                fq *= (1 - freq[i]) * (1 - freq[i])
                            elif allgeno[i] == 1:
                                fq *= 2 * (1 - freq[i]) * freq[i]
                            else:
                                fq *= freq[i] * freq[i]
                        # index is determined by genotype.
                        index = index * 3 + allgeno[i]
                    f += fitness[index] * fq
                # sum over other genotype
                sAll[loc * 3 + geno] = f
            # convert to form 0, s1, s2
            sAll[3 * loc + 1] = float(sAll[3 * loc + 1]) / sAll[3 * loc] - 1.
            sAll[3 * loc + 2] = float(sAll[3 * loc + 2]) / sAll[3 * loc] - 1.
            sAll[3 * loc] = 0.
        return sAll

    def _getS(self, gen, freq=None):
        '''Get fitness(gen) depending on the type of fitness. If self.fitness
        is a function, it is called to get a generation dependent fitness. The
        fitness value is then translated to 0, s1, s2. Marginal fitness is
        calculated using optional allele frequency (``freq``) if interaction is
        involved.
        '''
        # _fitness() expects parameters gen 
        if callable(self.fitness):
            fit = self.fitness(gen)
        else:
            fit = self.fitness
        s = []
        # simplest case when fitness only depends on gen if defined in fitness func:      
        # case 1: 3x self.nLoci no interaction 
        if len(fit) == 3 * self.nLoci:
            for i in range(self.nLoci):
                if fit[3 * i] == 0:
                    raise exceptions.ValueError('fitness['+ str(3 * i) + '] should be a non zero value.')
                s.append(0.)
                s.append(float(fit[3 * i + 1]) / float(fit[3 * i]) - 1.)
                s.append(float(fit[3 * i + 2]) / float(fit[3 * i]) - 1.)
        # case 2: same fitness for multiple loci
        elif len(fit) == 3 and self.nLoci > 1:
            if fit[0] == 0:
                raise exceptions.ValueError('fitness[0] should be a non zero value.')
            s.append(0.)
            s.append(float(fit[1]) / float(fit[0]) - 1.)
            s.append(float(fit[2]) / float(fit[0]) - 1.)
            s = s * self.nLoci
        # case 3: 3**self.nLoci, interaction
        elif len(fit) == 3**self.nLoci:
            # from fitness list, get s using allele frequency
            # Allele frequency for each subpopulation is passed and there will be
            # different s for each subpopulation because different allele frequencies.
            nSP = len(freq) / self.nLoci
            for sp in range(nSP):
                s.extend(self._marginalFitness(fit, [freq[sp + nSP * i] for i in range(self.nLoci)]))
        else:
            raise exceptions.ValueError('Wrong length of list of fitness: ' + str(len(fit)))
        return s

    def _getNextXt(self, curXt, Nt, s):
        '''Solve y from the formula and simulate forward trajectory.'''
        it = []
        xt = []
        nSP = len(Nt)
        for loc in range(self.nLoci):
            for sp in range(nSP):
                x = curXt[loc * nSP + sp]
                # if current allele freq in subpop sp at locus loc has already been 0 or 1,
                # set it to be 0 or 1 for next gens
                if x in [0, 1]:
                    xt.append(x)
                    continue
                # In the interaction case, s1, s2 will be different
                # from subpopulation to subpopulation.
                if len(s) == 3 * self.nLoci:
                    s1 = s[3 * loc + 1]
                    s2 = s[3 * loc + 2]
                else:
                    s1 = s[3 * loc + 1 + 3 * self.nLoci * sp]
                    s2 = s[3 * loc + 2 + 3 * self.nLoci * sp]                
                # with s1 and s2 on hand, calculate freq at the next generation
                y = x * (1 + s2 * x + s1 * (1 - x)) / (1 + s2 * x * x + 2 * s1 * x * (1 - x))
                # y is obtained, is the expected allele frequency for the next generation t+1
                it = GetRNG().randBinomial(2 * Nt[sp], y)
                xt.append(float(it) / (2 * Nt[sp]))
        return xt
        

    def _getPrevXt(self, curXt, Nt, s):
        '''Solve y from the backward formula and simulate backward trajectories
        in time. It takes parameter *curXt* in the form of a list. Nt is passed
        as a list for current population size/subpopulation sizes. The function
        returns a list, which contains allele frequencies of the previous
        generation at all loci through every subpopulation, arranged in the
        order of loc0_sp0, loc0_sp1, ..., loc1_sp0, loc1_sp1, ... and so on.
        '''
        #
        # given x(t)
        # calculate y=x(t-1)' by solving an equation
        #
        # x_t = y(1+s2 y+s1 (1-y))/(1+s2 y+2 s1 Y(1-y))
        it = []
        xt = []
        nSP = len(Nt)
        for loc in range(self.nLoci):
            for sp in range(nSP):
                x = curXt[loc * nSP + sp]
                # if current allele freq in subpop sp at locus loc has already been 0 or 1,
                # set it to be 0 or 1 for previous gens
                if x in [0, 1]:
                    xt.append(x)
                    continue
                # In the interaction case, s1, s2 will be different
                # from subpopulation to subpopulation.
                if len(s) == 3 * self.nLoci:
                    s1 = s[3 * loc + 1]
                    s2 = s[3 * loc + 2]
                else:
                    s1 = s[3 * loc + 1 + 3 * self.nLoci * sp]
                    s2 = s[3 * loc + 2 + 3 * self.nLoci * sp]
                # with s1 and s2 on hand, calculate freq at the previous generation
                if s1 == 0 and s2 == 0:
                    # special case when a = 0
                    y = x
                else:
                    a = s2 * x - 2 * s1 * x - s2 + s1
                    b = 2 * s1 * x - 1 - s1
                    c = float(x)
                    b2_4ac = b * b - 4 * a * c
                
                    if abs(a) < 1e-8:
                        y1 = float(-c) / float(b)
                        # y1 should be valid
                        y2 = 1000.
                    else:
                        y1 = (-b + b2_4ac**0.5) / (2 * a)
                        y2 = (-b - b2_4ac**0.5) / (2 * a)
                    #
                    # choose one of the solutions
                    if y1 >= 0 or y1 <= 1:
                         y = y2
                    else:
                         y = y1
                # y is obtained, is the expected allele frequency for the previous generation t-1
                it = GetRNG().randBinomial(int(2 * Nt[sp]), y)
                xt.append(float(it) / (2 * Nt[sp]))
        return xt

    def _simuForward(self, freq, freqEnd, genBegin, genEnd):
        '''Simulates a trajectory froward in time, starting from frequency
        ``freq`` at generation ``genBegin`` to frequency ranges specified in
        ``freqEnd`` at generation ``genEnd``. During the evolution, multiple
        subpopulations can be merged into one, and one population can be split
        into several subpopulations. The number of subpopulation is determined
        by the demographic function. The function raises an exception if the 
        simulated trajectory does not fall into ``destFeq``.
        '''
        # initialize a trajectory
        xt = trajectory(endGen = genEnd, nLoci = self.nLoci)
        xt._setFreq(freq, gen = genBegin, nSubPop = len(self._Nt(genBegin)))
        # begin at genBegin, go forward.
        gen = genBegin
        while 1:
            # when generation reaches the end, break.
            if gen == genEnd:
                # at all loci in all subpopulations, check if any allele frequency is outbound.
                nSP = len(Nt)
                for loc in range(self.nLoci):
                    for sp in range(nSP):
                        afq = nextXt[loc*nSP + sp]
                        if not freqEnd[loc][0] <= afq <= freqEnd[loc][1]:
                            if self.logger is not None:
                                self.logger.debug('Exception-F, restart due to:  Nsubpop= %d Nloci= %d alleleFreq= %.2f' % (sp, loc, afq))
                            raise exceptions.Exception('invalid')
                break
            # first get curXt, N(t+1), then calculate nextXt.
            curXt = xt.freq(gen)
            gen += 1
            nextNt = self._Nt(gen)
            if len(nextNt) > len(curXt) / self.nLoci:
                # split(forward sense) from one population to nSP subpopulations
                # get xt at next generation: t+1
                Nt = [sum(nextNt)]
                nextXt = self._getNextXt(curXt, Nt, self._getS(gen, xt.freq(gen-1)))
                # it nextXt to multiple subpopulations and assign... nextXt.
                NEWnextXt = []
                for idx, ind in enumerate(nextXt):
                    for sp in range(len(nextNt)):
                        NEWnextXt.append(ind) 
                xt._setFreq(NEWnextXt, gen=gen, nSubPop = len(nextNt))
            elif len(nextNt) < len(curXt) / self.nLoci:
                # check length of next Nt.
                if len(nextNt) != 1:
                    raise exceptions.ValueError('When merging event occurs between two adjacent generations in the forward sense, '
                                     + 'all subpops are required to merge into one population.')
                # merge (forward sense) from multiple subpopulations to one pop.
                Ntp = self._Nt(gen - 1)
                Nt = [int(x * 1.0 / sum(Ntp) * nextNt[0]) for x in Ntp]
                nextXt = self._getNextXt(curXt, Nt, self._getS(gen, xt.freq(gen-1)))
                # Merge nextXt and set to trajectory
                NEWnextXt = [
                    # sum over all subpopulations and calculate overall allele frequency
                    sum([x*y for x,y in zip(nextXt[i*len(Nt):(i+1)*len(Nt)], Nt)]) / sum(Nt)
                        # for all loci
                        for i in range(self.nLoci)]
                xt._setFreq(NEWnextXt, gen=gen, nSubPop = len(nextNt))
            else:
                Nt = nextNt
                nextXt = self._getNextXt(curXt, Nt, self._getS(gen, xt.freq(gen-1)))
                xt._setFreq(nextXt, gen=gen, nSubPop = len(nextNt))
        return xt
                
    def _simuBackward(self, genEnd, freq, minMutAge, maxMutAge):
        '''Simulates a trajectory backward from allele frequency ``freq`` at
        generation ``genEnd``. During the evolution, multiple subpopulations can
        be merged into one, and one population can be split into several
        subpopulations. The number of subpopulation is determined by the
        demographic function. If a simulated trajectory is shorter than
        ``minMutAge`` or is longer than ``maxMutAge``, the function will raise
        an exception.
        '''
        if genEnd <= 0:
            raise exceptions.ValueError("A positive ending generation is needed.")
        # done[i] is used to track at which generation a trajectory
        # is successfully generated at locus i.
        done = [False] * self.nLoci
        # initialize a trajectory
        xt = trajectory(endGen=genEnd, nLoci = self.nLoci)
        xt._setFreq(freq, gen=genEnd, nSubPop = len(self._Nt(genEnd)))
        # start from genEnd, go backward.
        gen = genEnd
        while 1:
            # first get curXt, N(t-1), then calculate prevXt
            curXt = xt.freq(gen)
            if self.logger is not None:
                self.logger.debug('Gen=%d, xt=%s'  % (gen, xt.freq(gen)))           
            gen -= 1
            if gen < 0 or gen + maxMutAge < genEnd:
                raise exceptions.Exception('tooLong')
            prevNt = self._Nt(gen)
            if len(prevNt) > len(curXt) / self.nLoci:
                # merge (forward sense)--> split from one population
                # to nSP subpopulations (backward sense)
                # get xt at previous generation: t-1
                Nt = [sum(prevNt)]
                prevXt = self._getPrevXt(curXt, Nt, self._getS(gen, xt.freq(gen+1)))
                # SPLIT prevXt to multiple subpopulations and assign an expanded
                # prevXt.
                assert len(prevXt) == self.nLoci
                assert len(Nt) == 1
                p = [float(x)/Nt[0] for x in prevNt]
                NEWprevXt = []
                for loc in range(self.nLoci):
                    NewIt = GetRNG().randMultinomialVal(int(prevXt[loc]*Nt[0]), p)
                    NEWprevXt.extend([float(x)/y for x,y in zip(NewIt, prevNt)])
                xt._setFreq(NEWprevXt, gen=gen, nSubPop = len(prevNt))
            elif len(prevNt) < len(curXt) / self.nLoci:
                # check length of previous Nt.
                if len(prevNt) != 1:
                    raise exceptions.ValueError('When merging event occurs between two adjacent generations in the backward sense, '
                                     + 'all subpops are required to merge into one population.') 
                # split (forward sense)--> merge from multiple subpopulations
                # to one population (backward sense)
                Ntp = self._Nt(gen + 1)
                Nt = [int(x * 1.0 / sum(Ntp) * prevNt[0]) for x in Ntp]
                prevXt = self._getPrevXt(curXt, Nt, self._getS(gen, xt.freq(gen+1)))
                # MERGE prevXt and set to trajectory
                NEWprevXt = [
                    # sum over all subpopulations and calculate overall allele frequency
                    sum([x*y for x,y in zip(prevXt[i*len(Nt):(i+1)*len(Nt)], Nt)]) / sum(Nt)
                    # for all loci
                    for i in range(self.nLoci)]
                xt._setFreq(NEWprevXt, gen=gen, nSubPop = len(prevNt))
            else:
                Nt = prevNt
                prevXt = self._getPrevXt(curXt, Nt, self._getS(gen, xt.freq(gen+1)))
                xt._setFreq(prevXt, gen=gen, nSubPop = len(prevNt))
            # set freq for previous generation
            # at all loci, check when prevIt is 0, if curIt is 1
            nSP = len(Nt)
            for loc in range(self.nLoci):
                doneNSP = [False] * nSP
                if done[loc]:
                    continue
                # loop over subpopulation
                for idx, (xtCur, xtPrev) in enumerate(zip(curXt[loc*nSP:(loc+1)*nSP],
                        prevXt[loc*nSP:(loc+1)*nSP])):
                    it = xtCur * self._Nt(gen + 1)[idx] * 2
                    if it == 0:
                        # already done in a previous generation
                        doneNSP[idx] = True
                        continue
                    if xtPrev == 0:
                        # success (judge when a trajectory is successfully generated)
                        doneNSP[idx] = gen
                        if self.logger is not None:
                            self.logger.debug('Exception-B1:  it= %d Nsubpop= %d Nloci= %d' % (it, idx, loc))
                        if genEnd - gen < minMutAge:
                            raise exceptions.Exception('tooShort')
                        break
                    elif xtPrev == 1: # fixed
                        if self.logger is not None:
                            self.logger.debug('Exception-B2:  it= %d Nsubpop= %d Nloci= %d' % (it, idx, loc))
                        raise exceptions.Exception('invalid')
                if False not in doneNSP:
                    done[loc] = True
            if False not in done:
                break
        return xt
            
    def simuForward(self, freq, freqEnd, genBegin = 0, genEnd = 0, maxAttempts = 10000):
        '''Simulate trajectories of multiple disease susceptibility loci using a
        forward time approach. A ``trajectory`` object is returned if the
        simulation succeeds. Otherwise, value ``None`` will be returned.

        This function accepts allele frequencies of alleles of multiple unlinked
        loci at the beginning generation (parameter ``freq``) at generation
        ``genBegin``, and expected *range* of allele frequencies of these
        alleles (parameter ``freqEnd``) at the end of generation ``genEnd``.
        Depending on the number of loci and subpopulations, ``freq`` can be
        a number or a list of frequencies for each locus at each subpopulation,
        and ``freqEnd`` can be a list, or a list for ranges of frequencies for
        each locus at each subpopulation. This simulator will simulate a
        trajectory generation by generation and restart if the resulting
        frequencies do not fall into specified range of frequencies. This
        simulator will return ``None`` if no valid trajectory is found after
        ``maxAttempts`` attemps.
        '''
        # freqEnd
        if type(freqEnd) not in [type(()), type([])] or len(freqEnd) == 0:
            raise exceptions.ValueError('A list of frequency range is expected')
        elif type(freqEnd[0]) not in [type(()), type([])]:
            if len(freqEnd) == 2:
                freqEnd = [freqEnd]
            else:
                raise exceptions.ValueError('A list of frequency range is expected.')
        if len(freqEnd) != self.nLoci:
            raise exceptions.ValueError('Please specify a frequency range for each locus')
        for i in range(self.nLoci):
            if len(freqEnd[i]) != 2:
                raise exceptions.ValueError('Please specify frequency range of each marker')
            if freqEnd[i][0] >= freqEnd[i][1]:
                raise exceptions.ValueError('Invalid frequency range %s' % freq[i])
        if not(genBegin <= genEnd):
            raise exceptions.ValueError('Beginning generation should be less than ending generation')
        self.errorCount['invalid'] = 0
        for failedCount in range(maxAttempts):
            try:
                return self._simuForward(freq, freqEnd, genBegin, genEnd)
            except exceptions.Exception, e:
                if e.args[0] == 'invalid':
                    self.errorCount['invalid'] += 1
                else:
                    print 'Unknown error occurs'
                    print e
                    raise
        return None
    
    def simuBackward(self, genEnd, freq, minMutAge = 0, maxMutAge = 0, maxAttempts = 1000):
        '''Simulate trajectories of multiple disease susceptibility loci using a
        forward time approach. A ``trajectory`` object is returned if the
        simulation succeeds. Otherwise, value ``None`` will be returned.

        This function accepts allele frequencies of alleles of multiple unlinked
        loci at the beginning generation (parameter ``freq``) at generation
        ``genBegin``, and expected *range* of allele frequencies of these
        alleles (parameter ``freqEnd``) at the end of generation ``genEnd``.
        Depending on the number of loci and subpopulations, ``freq`` can be
        a number or a list of frequencies for each locus at each subpopulation,
        and ``freqEnd`` can be a list, or a list for ranges of frequencies for
        each locus at each subpopulation. This simulator will simulate a
        trajectory generation by generation and restart if the resulting
        frequencies do not fall into specified range of frequencies. During the
        evolution, multiple subpopulations can be merged into one, and one
        population can be split into several subpopulations. The number of
        subpopulation is determined by the demographic function. This
        simulator will return ``None`` if no valid trajectory is found after
        ``maxAttempts`` attemps.
        '''
        self.maxMutAge = maxMutAge
        self.minMutAge = minMutAge
        
        if genEnd > 0 and minMutAge > genEnd:
            minMutAge = genEnd
        if genEnd > 0 and maxMutAge == 0:
            maxMutAge = genEnd

        if not maxMutAge >= minMutAge:
            raise exceptions.ValueError('maxMutAge should >= minMutAge')
        if genEnd == 0 and (callable(self.N) or callable(self.fitness)):
            raise exceptions.ValueError('genEnd should be > 0 if N or fitness is defined in the form of function')
        if genEnd > 0 and genEnd < maxMutAge:
            raise exceptions.ValueError('genEnd should be >= maxMutAge')
        self.errorCount['invalid'] = 0
        self.errorCount['tooLong'] = 0
        self.errorCount['tooShort'] = 0
        for failedCount in range(maxAttempts):
            try:
                return self._simuBackward(genEnd, freq, minMutAge, maxMutAge)
            except exceptions.Exception, e:
                if e.args[0] == 'tooLong':
                    self.errorCount['tooLong'] += 1
                elif e.args[0] == 'tooShort':
                    self.errorCount['tooShort'] += 1
                elif e.args[0] == 'invalid' :
                    self.errorCount['invalid'] += 1
                else:
                    print e
                    raise
        return None


def ForwardTrajectory(N, freq, freqEnd, genBegin, genEnd, nLoci = 1,
        fitness = None, maxAttempts = 10000, logger=None):
    '''Given a demographic model (parameter *N*, which can be a constant or a
    function that returns subpopulation sizes at each generation) and the
    fitness of genotype at one or more loci (parameter *fitness*, which can be
    ``None`` (no selection), a list of 3 items (for genotype *AA*, *Aa*, and
    *aa*), 3*nLoci (for multiple loci), or 3**nLoci (for each combination of
    genotype)), this function simulates a trajectory of one or more unlinked
    loci (parameter *nLoci*) from allele frequency *freq* at generation
    *genBegin* forward in time, until it reaches generation *genEnd*. A
    ``trajectory`` object will be returned if the allele frequency falls
    into specified ranges (``freqEnd``). ``None`` will be returned if no valid
    trajectory is simulated after ``maxAttempts`` attempts. Please refer to
    class ``trajectory``, ``trajectorySimulator`` and their member functions
    for more details.
    '''
    return trajectorySimulator(N, nLoci, fitness, logger).simuForward(freq, freqEnd, genBegin, genEnd,
                                                              maxAttempts)
    
def BackwardTrajectory(N, genEnd, freq, nLoci=1, fitness=None,
        minMutAge = 0, maxMutAge = 0, maxAttempts = 1000, logger=None):
    '''Given a demographic model (parameter *N*, which can be a constant or a
    function that returns subpopulation sizes at each generation) and the
    fitness of genotype at one or more loci (parameter *fitness*, which can be
    ``None`` (no selection), a list of 3 (for genotype *AA*, *Aa*, and *aa*),
    3*nLoci (for multiple loci), or 3**nLoci (for each combination of genotype)
    items), this function simulates a trajectory of one or more unlinked loci
    (parameter *nLoci*) from allele frequency *freq* at generation *genEnd*
    backward in time, until all alleles get lost. A ``trajectory`` object will
    be returned if the length of simulated trajectory is longer than
    ``minMutAge`` and shorter than ``maxMutAge`` (if specified). ``None`` will
    be returned if no valid trajectory is simulated after ``maxAttempts``
    attempts. Please refer to class ``trajectory``, ``trajectorySimulator`` and
    their member functions for more details.
    '''
    return trajectorySimulator(N, nLoci, fitness, logger).simuBackward(genEnd, freq, minMutAge, maxMutAge, 
                                                              maxAttempts)


if __name__ == "__main__":
    pass

