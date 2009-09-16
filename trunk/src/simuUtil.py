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
        raise ValueError('A list of size of island groups is expected for parameter n')
    nIslands = sum(n)
    if type(r1) in [type(0), type(1.)]:
        r1 = [r1] * len(n)
    elif len(r1) != len(n):
        raise ValueError('If multiple r1 is given, it should be given to all island groups.')
    #
    if type(r2) in [type(0), type(1.)]:
        r2 = [r2] * len(n)
    elif len(r2) != len(n):
        raise ValueError('If multiple r2 is given, it should be given to all island groups.')
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
            raise ValueError('please define freq as a single number or a list with specified number of elements.')

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
    '''
    Simulate trajectories of disease susceptibility loci using an extension of the backward
    method described in Slatkin 2001 or forward algorithms.

    Tracking allele frequencies of alleles on all loci.

    Class trajectory(...) takes four arguments, at least three (N, fitness, nLoci) need to be
    specified by the user:

    N
        population size, which may be passed in as a constant number or an array of subpop sizes
        or a user defined function, NtFunc(gen), which returns population size at each generation.
        
    fitness
        selection pressure for all loci, which can be passed in as a constant array showing
        fitness for [AA, Aa, aa, BB, Bb, bb,...] or a user defined function, fitnessFunc(gen),
        which returns selection pressure at each generation.
    
    nLoci
        number of Loci, which should be passed in as a constant integer number with its value
        equal to or larger than 1.

    logger
        Logged messages have levels of importance associated with them. The default levels
        provided are DEBUG, INFO, WARNING, ERROR and CRITICAL. User should indicate the importance
        of a logged message. Default value for logger is None.
    '''

    def __init__(self, N, fitness, nLoci, logger=None):
        '''
        Initialization of global parameters.
        
        Parameter description:

        N
            pass in population size ''N'' to global variable self.N, which may be a number, a list
            of subpop sizes or a user-defined function.

        fitness
            pass in selection pressure ''fitness'' to global variable self.fitness, which may take
            forms such as, a list of three values in the case of single locus or same fitness for any
            locus within multiple loci, a list of 3 * nLoci values in the case of multi-loci without
            interaction, a list of 3**nLoci values when multi-loci and interaction situation are both considered.

        nLoci
            pass in number of loci ''nLoci'' to global variable self.nLoci.

        logger
            pass in level of importance ''logger'' to global variable self.logger, which is a logging
            object that can be used to record warnings and error messages. 

        errorCount
            define global dictionary object self.errorCount as the variable to
            record counts of distinct forms of errors that occur during the simulation.        
        '''
        # a vector of subpopulation sizes is needed
        if type(N) in [type(1), type(1L)]:
            self.N = [N]
        else: # N is a list or a function
            self.N = N
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

    def _interFitness(self, fitness, freq):
        'Get fitness for all loci when there is interaction'
        sAll = [0] * (3 * self.nLoci)
        # each locus
        for loc in range(self.nLoci):
            # each genotype AA, Aa and aa (geno is the number of disease allele)
            # calculate
            #     f(X=Aa) = Sum_g P(Y=g) * f(X=Aa, Y=g)
            for geno in range(3):
                # iterate through OTHER DSL
                allgeno = [0] * self.nLoci
                # set myself
                allgeno[loc] = geno
                # iterate through others
                f = 0.
                for it in range(3**(self.nLoci - 1)):
                    # assign allgeno
                    num = it
                    for l in range(self.nLoci):
                        if l == loc:
                            continue
                        allgeno[l] = num % 3
                        num /= 3
                    f += self._fitOfGeno(loc, allgeno, fitness, freq)
                # sum over other genotype
                sAll[loc * 3 + geno] = f
            # convert to form 1, s1, s2
            sAll[3 * loc + 1] = float(sAll[3 * loc + 1]) / sAll[3 * loc] - 1.
            sAll[3 * loc + 2] = float(sAll[3 * loc + 2]) / sAll[3 * loc] - 1.
            sAll[3 * loc] = 0.
        return sAll

    def _fitOfGeno(self, loc, allgeno, fitness, freq):
        'Return the fitness of genotype for func _interFitness'
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
            index = index * 3 + allgeno[i]
            if fq == 0:
                return 0.
        return fitness[index] * fq

    def _getS(self, gen, freq=None):
        '''
        Get fitness(gen) depending on the type of fitness.
        Get s1, s2 from f1, f2, f3.
        Currently, selection pressure only varies among multiple loci and/or depends on gen,
        but does not rely on changes of allele frequency or different subpopulations.
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
                    raise ValueError('fitness['+ str(3 * i) + '] should be a non zero value.')
                s.append(0.)
                s.append(float(fit[3 * i + 1]) / float(fit[3 * i]) - 1.)
                s.append(float(fit[3 * i + 2]) / float(fit[3 * i]) - 1.)
        # case 2: same fitness for multiple loci
        elif len(fit) == 3 and self.nLoci > 1:
            if fit[0] == 0:
                raise ValueError('fitness[0] should be a non zero value.')
            s.append(0.)
            s.append(float(fit[1]) / float(fit[0]) - 1.)
            s.append(float(fit[2]) / float(fit[0]) - 1.)
            s = s * self.nLoci
        # case 3: 3^^self.nLoci, interaction
        elif len(fit) == 3**self.nLoci:
            # from fitness list, get s using allele frequency
            # Allele frequency for each subpopulation is passed and there will be
            # different s for each subpopulation because different allele frequencies.
            nSP = len(freq) / self.nLoci
            for sp in range(nSP):
                s.extend(self._interFitness(fit, [freq[sp + nSP * i] for i in range(self.nLoci)]))
        else:
            raise ValueError('Wrong length of list of fitness: ' + str(len(fit)))
        return s

    def _getNextXt(self, curXt, Nt, s, ploidy):
        '''
        Solve y from the formula and simulate forward trajectory.
        '''
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
                it = GetRNG().randBinomial(ploidy * Nt[sp], y)
                xt.append(float(it) / (ploidy * Nt[sp]))
        return xt
        

    def _getPrevXt(self, curXt, Nt, s, ploidy):
        '''
        Solve y from the backward formula and simulate backward trajectories in time
        It takes parameter curXt in the form of a list.
        Nt is passed in as a list for current population size/subpopulation sizes.
        The function returns a list, which contains allele frequencies of one
        previous generation at all loci through every subpopulation, it is arranged
        in the order of loc0_sp0, loc0_sp1, ..., loc1_sp0, loc1_sp1, ... and so on. 
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

                    # choose one of the solutions
                    if y1 >= 0 or y1 <= 1:
                         y = y2
                    else:
                         y = y1
                # y is obtained, is the expected allele frequency for the previous generation t-1
                it = GetRNG().randBinomial(int(ploidy * Nt[sp]), y)
                xt.append(float(it) / (ploidy * Nt[sp]))
        return xt
        
    def _simuForward(self, freq, destFreq, genBegin, genEnd, ploidy, maxAttempts):
        '''
        This function simulates forward trajectory for one time.
        During the simulation, variable number of subpopulations within different
        generations will be allowed. It's introduced by merging events and splitting
        events in forward sense. Changes of subpopulation sizes need to meet specific
        requirements. If subpopulations are merged, they must combine into *one*
        population. If a splitting event has to occur, it can only originate from
        generations where no subpopulation existed. 
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
                        if not destFreq[loc][0] <= afq <= destFreq[loc][1]:
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
                nextXt = self._getNextXt(curXt, Nt, self._getS(gen, xt.freq(gen-1)), ploidy)
                # it nextXt to multiple subpopulations and assign... nextXt.
                NEWnextXt = []
                for idx, ind in enumerate(nextXt):
                    for sp in range(len(nextNt)):
                        NEWnextXt.append(ind) 
                xt._setFreq(NEWnextXt, gen=gen, nSubPop = len(nextNt))
            elif len(nextNt) < len(curXt) / self.nLoci:
                # check length of next Nt.
                if len(nextNt) != 1:
                    raise ValueError('When merging event occurs between two adjacent generations in the forward sense, '
                                     + 'all subpops are required to merge into one population.')
                # merge (forward sense) from multiple subpopulations to one pop.
                Ntp = self._Nt(gen - 1)
                Nt = [int(x * 1.0 / sum(Ntp) * nextNt[0]) for x in Ntp]
                nextXt = self._getNextXt(curXt, Nt, self._getS(gen, xt.freq(gen-1)), ploidy)
                # Merge nextXt and set to trajectory
                NEWnextXt = [
                    # sum over all subpopulations and calculate overall allele frequency
                    sum([x*y for x,y in zip(nextXt[i*len(Nt):(i+1)*len(Nt)], Nt)]) / sum(Nt)
                        # for all loci
                        for i in range(self.nLoci)]
                xt._setFreq(NEWnextXt, gen=gen, nSubPop = len(nextNt))
            else:
                Nt = nextNt
                nextXt = self._getNextXt(curXt, Nt, self._getS(gen, xt.freq(gen-1)), ploidy)
                xt._setFreq(nextXt, gen=gen, nSubPop = len(nextNt))

        return xt
                
    def _simuBackward(self, genEnd, freq, minMutAge, maxMutAge, ploidy,
                     restartIfFail):
        '''
        This function simulates backward trajectory for one time.
        During the simulation, various number of subpopulations within different
        generations will be allowed. It's introduced by merging events and splitting
        events in forward sense, correspondingly splitting and merging in backward
        simulation. Changes of subpopulation sizes need to meet specific requirements.
        If subpopulations are merged, they must combine into *one* population.
        If a splitting event has to occur, it can only originate from generations
        where no subpopulation existed.
        '''
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
                prevXt = self._getPrevXt(curXt, Nt, self._getS(gen, xt.freq(gen+1)), ploidy)
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
                    raise ValueError('When merging event occurs between two adjacent generations in the backward sense, '
                                     + 'all subpops are required to merge into one population.') 
                # split (forward sense)--> merge from multiple subpopulations
                # to one population (backward sense)
                Ntp = self._Nt(gen + 1)
                Nt = [int(x * 1.0 / sum(Ntp) * prevNt[0]) for x in Ntp]
                prevXt = self._getPrevXt(curXt, Nt, self._getS(gen, xt.freq(gen+1)), ploidy)
                # MERGE prevXt and set to trajectory
                NEWprevXt = [
                    # sum over all subpopulations and calculate overall allele frequency
                    sum([x*y for x,y in zip(prevXt[i*len(Nt):(i+1)*len(Nt)], Nt)]) / sum(Nt)
                    # for all loci
                    for i in range(self.nLoci)]
                xt._setFreq(NEWprevXt, gen=gen, nSubPop = len(prevNt))
            else:
                Nt = prevNt
                prevXt = self._getPrevXt(curXt, Nt, self._getS(gen, xt.freq(gen+1)), ploidy)
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
                    it = xtCur * self._Nt(gen + 1)[idx] * ploidy
                    if it == 0:
                        # already done in a previous generation
                        doneNSP[idx] = True
                        continue
                    if xtPrev == 0:
                        if it > 1:
                            # When prevIt equals to 2 and curIt is 0, we consider
                            # that case to be a successful simulation the same as that
                            # when pervIt is 1.
                            if it <= 2 and self._Nt(gen + 1)[idx] > 100:
                                xt.freq(gen + 1)[loc * nSP + idx] /= it
                                it = 1
                            else:
                                if self.logger is not None:
                                    self.logger.debug('Exception-B1:  it= %d Nsubpop= %d Nloci= %d' % (it, idx, loc))
                                raise exceptions.Exception('invalid')
                        else:  # it == 1
                            assert abs(it - 1.) < 1e-5
                            # success (judge when a trajectory is successfully generated)
                            doneNSP[idx] = gen
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
            
    def simuForward(self, freq, destFreq, genBegin = 0, genEnd = 0, ploidy=2, maxAttempts = 10000):
        '''
        simulate trajectories of multiple disease susceptibility loci using a forward time approach

        Return the trajectory for each locus at each subpopulation. In the order of
        OC0: sp0, sp1, sp2, ..., LOC1: sp0, sp1, sp2,...
        Each trajectory will have length genEnd - genBegin + 1
        
        parameter description:

        genBegin
            starting generation number
        
        genEnd
            ending generation number
        
        freq
            allele frequencies of alleles of multiple unlinked loci at the beginning
            generation.
        
        destFreq
            expected *range* of allele frequencies of alleles of multiple unlinked loci,
            at generation genEnd , with all subpopulation combined. If there are two loci,
            it can be [[0.08, 0.12], [0.19, 0.21]]
        
        ploidy
            number of chromosomes will be N*ploidy

        maxAttempts
            How many times to try to get a valid path? Default 10,000
        '''
        for i in range(self.nLoci):
            if len(destFreq[i]) != 2:
                raise ValueError('Please specify frequency range of each marker')
        if not(genBegin <= genEnd):
            raise ValueError('Beginning generation should be less than ending generation')
        self.errorCount['invalid'] = 0
        for failedCount in range(maxAttempts):
            try:
                return self._simuForward(freq, destFreq, genBegin, genEnd, ploidy, maxAttempts)
            except exceptions.Exception, e:
                if e.args[0] == 'invalid':
                    self.errorCount['invalid'] += 1
                else:
                    print 'Unknown error occurs'
                    print e
                    raise          
    
    def simuBackward(self, genEnd, freq, minMutAge = 0, maxMutAge = 100000, ploidy = 2,
                     restartIfFail = False, maxAttempts = 1000):
        '''
        simulate trajectories of multiple disease susceptibility loci using an extension of the
        backward method described in Slatkin 2001.

        parameter description:
        
        genEnd
            genenration number in the end, namely the current generation number in the
            backward trajectory.
        
        freq
            expected allele frequencies of alleles of multiple unlinked loci.
            It may take forms of a single value or a list of values with number of elememts
            equal to nLoci. In the case of single value, such freq will be shared by all loci. 
            FIXME: There are multiple subpopultions and only one frequency is given,
            the same frequency will be used for all subpopulations. Users can specify
            different allele frequencies for each subpopulation using the long form...
            
        minMutAge
            minimum generation number. The process will restart if the trajectory is
            less than it. Default to 0.
        
        maxMutAge
            maximum generation number. The process will terminate or restart if it can
            not reach allele zero after T generations. Default to 100,000, roughly
            2,000,000 years which is longer than human history.
        
        ploidy
            number of chromosomes will be N * ploidy
        
        restartIfFail
            If the process can not finish after T generations, restart if
            restartIfFail = True, otherwise return. Default to False.
        
        maxAttempts
            How many times to try to get a valid path? Default to 1000.
        
        Return the trajectory for each locus at each subpopulation. In the order of
            LOC0: sp0, sp1, sp2,..., LOC1: sp0, sp1, sp2,...
        '''
        self.maxMutAge = maxMutAge
        self.minMutAge = minMutAge
        
        if genEnd > 0 and minMutAge > genEnd:
            minMutAge = genEnd
        if genEnd > 0 and maxMutAge > genEnd:
            maxMutAge = genEnd

        if not(maxMutAge >= minMutAge):
            raise ValueError('maxMutAge should >= minMutAge')
        if (genEnd == 0 and (callable(self.N) or callable(self.fitness))):
            raise ValueError('genEnd should be > 0 if N or fitness is defined in the form of function')
        if genEnd > 0 and genEnd < maxMutAge:
            raise ValueError('genEnd should be >= maxMutAge')
        self.errorCount['invalid'] = 0
        for failedCount in range(maxAttempts):
            try:
                return self._simuBackward(genEnd, freq, minMutAge, maxMutAge, ploidy, restartIfFail)
            except exceptions.Exception, e:
                if e.args[0] == 'tooLong':
                    self.errorCount['tooLong'] += 1
                    self.errorCount['invalid'] += 1
                elif e.args[0] == 'tooShort':
                    self.errorCount['tooShort'] += 1
                    self.errorCount['invalid'] += 1
                elif e.args[0] == 'invalid' :
                    self.errorCount['invalid'] += 1
                else:
                    print 'Unknown error occurs'
                    print e
                    raise
        # if no valid trajectory is successfully simulated when reaching maximum attempts,
        # return an empty object to class trajectory 
        return trajectory(endGen=genEnd, nLoci = self.nLoci)
        
                    
    def message(self):
        '''
        Report how many 'too long paths' errors, 'too short paths' errors
        and 'invalid paths errors' generated during the simulation.
        '''
        if self.errorCount['tooLong'] > 0:
            print 'Trajectories regenerated due to long', self.maxMutAge, 'path', self.errorCount['tooLong'], 'times.'
        if self.errorCount['tooShort'] > 0:
            print 'Trajectories regenerated due to short', self.minMutAge, 'path', self.errorCount['tooShort'], 'times.'
        if self.errorCount['invalid'] > 0:
            print 'Trajectories regenerated due to invalid path:', self.errorCount['invalid'], 'times.'
        else:
            print 'No invalid trajectory generated.'
        return 


def ForwardTrajectory(N, fitness, nLoci, freq, destFreq, genBegin = 0, genEnd = 0, ploidy=2, maxAttempts = 10000,
                    logger=None):
    '''
    Return an object to class trajectory in forward simulation.
    '''
    return trajectorySimulator(N, fitness, nLoci, logger).simuForward(freq, destFreq, genBegin, genEnd,
                                                              ploidy, maxAttempts)
    
def BackwardTrajectory(N, fitness, nLoci, genEnd, freq, minMutAge = 0, maxMutAge = 100000, ploidy = 2,
                     restartIfFail = False, maxAttempts = 1000, logger=None):
    '''
    Return an object to class trajectory in backward simulation.
    '''
    return trajectorySimulator(N, fitness, nLoci, logger).simuBackward(genEnd, freq, minMutAge, maxMutAge, ploidy,
                                                              restartIfFail, maxAttempts)


if __name__ == "__main__":
    pass



