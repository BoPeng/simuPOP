#!/usr/bin/env python

#
# $File: utils.py $
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

__all__ = [
    'ViewVars',
    'MigrIslandRates',
    'MigrHierarchicalIslandRates',
    'MigrSteppingStoneRates',
    'SaveCSV',
    'simuProgress',
    'trajectory',
    'trajectorySimulator',
    'BackwardTrajectory',
    'ForwardTrajectory',
]

import exceptions, operator, types, os, sys, re

from simuOpt import simuOptions

# avoid duplicated banner message
q = simuOptions['Quiet']
simuOptions['Quiet'] = True
from simuPOP import Male, Female, pointMutator, GetRNG
simuOptions['Quiet'] = q

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
        # self.traj stores a list of frequencies for each loci.
        # at each generation, the frequencies are saved as
        #   [[loc0_sp0, loc1_sp0], [loc0_sp1, loc1_sp1], ]
        # and so on. That is to say, the frequencies should be accessed as
        # self.traj[gen][sp][loc]
        self.traj = {}
        self.endGen = endGen
        self.nLoci = nLoci
    
    def _beginGen(self):
        '''Return starting generation of all trajectories'''
        return min(self.traj.keys())

    def _freq(self, gen):
        '''Return frequencies at all subpopulations at generation *gen*.'''
        if not self.traj.has_key(gen):
            # assuming no subpopulations
            return [[0.] * self.nLoci]
        assert len(self.traj[gen][0]) == self.nLoci
        return self.traj[gen]
    
    def freq(self, gen, subPop):
        '''Return frequencies of all loci in subpopulation *subPop* at
        generation *gen* of the simulated trajectory. Allele frequencies are
        assumed to be zero if *gen* is out of range of the simulated
        trajectory.
        '''
        if not self.traj.has_key(gen):
            # assuming no subpopulations
            return [0.] * self.nLoci
        assert len(self.traj[gen][subPop]) == self.nLoci
        return self.traj[gen][subPop]

    def func(self):
        '''Return a Python function that returns allele frequencies for each
        locus at specified loci. If there are multiple subpopulations, allele
        frequencies are arranged in the order of ``loc0_sp0``, ``loc1_sp0``,
        ..., ``loc0_sp1``, ``loc1_sp1``, ... and so on. The returned function
        can be supplied directly to the ``freqFunc`` parameter of a controlled
        random mating scheme (``controlledRandomMating``) or a homogeneous
        mating scheme that uses a controlled offspring generator
        (``controlledOffspringGenerator``).
        '''
        def trajFunc(gen):
            if not self.traj.has_key(gen):
                return [0.] * self.nLoci
            freq = []
            for spFreq in self.traj[gen]:
                freq.extend(spFreq)
            return freq
        return trajFunc
    
    def mutators(self, loci, inds=0, allele=1, *args, **kwargs):
        '''Return a list of ``pointMutator`` operators that introduce mutants
        at the beginning of simulated trajectories. These mutators should be
        added to the ``ops`` parameter of ``simulator.evolve`` function to
        introduce a mutant at the beginning of a generation with zero allele
        frequency before mating, and a positive allele frequency after mating.
        A parameter ``loci`` is needed to specify actual loci indexes in the
        real forward simulation. Other than default parameters ``inds=0`` and
        ``allele=1``, additional parameters could be passed to point mutator
        as keyward parameters.
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
            for sp in range(len(self.traj[gen])):
                for loc in range(self.nLoci):
                    if self.traj[gen][sp][loc] == 0 and self.traj[gen + 1][sp][loc] > 0:
                        if type(loci) in [type(()), type([])]:
                            if len(loci) != self.nLoci:
                                raise exceptions.ValueError('%d loci is expected' % self.nLoci)
                            mut.append(pointMutator(inds=inds, loci=loci[loc], allele=allele,
                                subPops=sp, at=gen + 1, *args, **kwargs))
                        elif self.nLoci == 1 and type(loci) == type(0):
                            mut.append(pointMutator(inds=inds, loci=loci, allele=allele,
                                subPops=sp, at=gen + 1, *args, **kwargs))
                        else:
                            raise exceptions.ValueError('Invalid parameter loci')
        return mut

    def _setFreq(self, freq, gen):
        '''This function sets frequency *freq* at specified generation *gen*.
        *nSubPop* is used if frequency for multiple subpopulations are given.
        '''
        assert type(freq) in [type(()), type([])]
        # deep copy to avoid trouble.
        self.traj[gen] = []
        for spFreq in freq:
            assert len(spFreq) == self.nLoci
            assert type(spFreq[0]) not in [type(()), type([])]
            self.traj[gen].append([x for x in spFreq])

    def plot(self, filename=None, **kwargs):
        '''Plot simulated trajectory using ``R`` through a Python module ``rpy``.
        The function will return silently if module ``plotter`` cannot be
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
        module ``plotter``. Allowed prefixes are ``par``, ``plot``, ``lines``
        and ``dev_print``. Allowed repeating suffix are ``loc`` and ``sp``.
        For example, you could use parameter ``plot_ylim`` to reset the default
        value of ``ylim`` in R function ``plot``.
        '''
        try:
            import plotter
        except ImportError:
            return
        #
        args = plotter.derivedArgs(
            defaultFuncs = ['plot', 'lines'],
            allFuncs = ['par', 'plot', 'lines', 'dev_print'],
            suffixes = ['loc', 'sp'],
            defaultParams = {
                'plot_ylim': [0, 1],
                'plot_ylab': 'Allele frequency',
                'plot_xlab': 'Generation'
            },
            **kwargs
        )
        # device
        plotter.newDevice()
        #
        gens = self.traj.keys()
        gens.sort()
        plotter.r.par(**args.getArgs('par', None))
        plotter.r.plot(gens[0], 0,
            **args.getArgs('plot', None, xlim=(min(gens), max(gens)), type='n'))
        allLines = []
        for gen in gens:
            for sp in range(len(self.traj[gen])):
                if len(allLines) <= sp:
                    allLines.append([])
                lines = allLines[sp]
                if len(lines) != len(self.traj[gen][0]):
                    if len(lines) != 0:
                        for loc, line in enumerate(lines):
                            # remove leading zero's
                            for start in range(len(line)):
                                if line[start] != 0:
                                    break
                            if start == len(line) - 1:
                                continue
                            r.lines(x=range(beginGen - len(line) + start, gen), y=line[start:],
                                    **args.getArgs('lines', None, sp=sp, loc=loc))
                    # the first point
                    allLines[sp] = [[x] for x in self.traj[gen][sp]]
                else:
                    for loc, x in enumerate(self.traj[gen][sp]):
                        lines[loc].append(x)
        # plot the lines again
        for sp, lines in enumerate(allLines):
            for loc, line in enumerate(lines):
                # remove leading zero's
                for start in range(len(line)):
                    if line[start] != 0:
                        break
                if start == len(line) - 1:
                    continue
                beginGen = gen - len(line) + start
                plotter.r.lines(x=range(beginGen, gen), y=line[start:],
                        **args.getArgs('lines', None, sp=sp, loc=loc))
        #
        plotter.saveFigure(**args.getArgs('dev_print', None, file=filename))
        

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
        range of generations. This class accepts the following parameters
        
        N
            Parameter *N* accepts either a constant number for population size
            (e.g. N=1000), a list of subpopulation sizes (e.g. N=[1000, 2000]),
            or a demographic function that returns population or subpopulation
            sizes at each generation. During the evolution, multiple
            subpopulations can be merged into one, and one population can be
            split into several subpopulations. The number of subpopulation is
            determined by the return value of the demographic function. Note
            that *N* should be considered as the population size at the end of
            specified generation.

        nLoci
            Number of unlinked loci for which trajectories of allele
            frequencies are simulated. We assume a diploid population with
            diallelic loci. The trajectory represents frequencies of a 

        fitness
            Parameter fitness can be ``None`` (no selection), a list of fitness
            values for genotype with 0, 1, and 2 disease alleles (*AA*, *Aa*,
            and *aa*) at one or more loci; or a function that returns fitness
            values at each generation. When multiple loci are involved
            (*nLoci*), *fitness* can be a list of 3 (the same fitness values
            for all loci), a list of 3*nLoci (different fitness values for each
            locus) or a list of 3**nLoci (fitness value for each combination of
            genotype). The fitness function should accept generation number and
            a subpopulation index. The latter parameter allows, and is the only
            way to specify different fitness in each subpopulation.
        
        logger
            A logging object (see Python module ``logging``) that can be used
            to output intermediate results with debug information.
        '''
        # a vector of subpopulation sizes is needed
        if type(N) in [type(1), type(1L)]:
            self.N = [N]
        else: # N is a list or a function
            self.N = N
        if fitness is None:
            self.fitness = [1, 1, 1]
        else:
            # fitness is a list or a function
            if type(fitness) in [type(()), type([])] and len(fitness) not in [3, 3*nLoci, 3**nLoci]:
                raise exceptions.ValueError('Invalid list of fitness.')
            self.fitness = fitness
        self.logger = logger
        self.nLoci = nLoci
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
            # a constant list
            return self.N

    def _marginalFitness(self, fitness, freq):
        '''Convert interaction fitness (3**n elements) to marginal fitness
        (3*n elements) using given allele frequency. The marginal fitnesses
        are calculated using formula:
                f(X=Aa) = Sum_g P(Y=g) * f(X=Aa, Y=g)
        where g is genotype at all other loci.
        '''
        assert len(freq) == 2
        assert len(fitness) == 3 ** self.nLoci
        s = [0] * (3 * self.nLoci)
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
                        if l != loc:
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
                s[loc * 3 + geno] = f
            # convert to form 0, s1, s2
            s[3 * loc + 1] = float(s[3 * loc + 1]) / s[3 * loc] - 1.
            s[3 * loc + 2] = float(s[3 * loc + 2]) / s[3 * loc] - 1.
            s[3 * loc] = 0.
        return s

    def _getS(self, gen, subPop, freq):
        '''Get s1, s2 for subpopulation *subPop* at generation *gen*. If
        self.fitness is a function, it is called with *gen* and *subPop* to get
        a generation and subpopulation specific fitness value. The fitness
        value is then translated to 0, s1, s2. If interactions are involved,
        marginal fitness is calculated using allele frequency (``freq``) in
        subpopulation *subPop*. 
        '''
        assert len(freq) == self.nLoci
        # _fitness() expects parameters gen and a subpopulation index
        if callable(self.fitness):
            fit = self.fitness(gen, subPop)
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
            s.extend(self._marginalFitness(fit, freq))
        else:
            raise exceptions.ValueError('Wrong length of list of fitness: ' + str(len(fit)))
        return s

    def _getNextXt(self, curXt, Nt, s):
        '''Solve y from the formula and simulate allele frequencies in the next
        generation. All parameters are assumed to be for one subpopulation. Nt
        is the population size of at the end of the current generation (or the
        next generation)'''
        assert len(curXt) == self.nLoci
        assert type(Nt) not in [type(()), type([])]
        it = []
        xt = []
        for loc in range(self.nLoci):
            # if current allele freq in subpop sp at locus loc has already been 0 or 1,
            # set it to be 0 or 1 for next gens
            x = curXt[loc]
            if x in [0, 1]:
                xt.append(x)
                continue
            s1 = s[3 * loc + 1]
            s2 = s[3 * loc + 2]
            # with s1 and s2 on hand, calculate freq at the next generation
            y = x * (1 + s2 * x + s1 * (1 - x)) / (1 + s2 * x * x + 2 * s1 * x * (1 - x))
            # y is obtained, is the expected allele frequency for the next generation t+1
            it = GetRNG().randBinomial(2 * Nt, y)
            xt.append(float(it) / (2 * Nt))
        return xt

    def _getPrevXt(self, curXt, Nt, s):
        '''Solve y from the backward formula and simulate allele frequencies in
        the previous generation. All parameters are assumed to be for one
        subpopulation. Nt is the population size at the begining of the current
        generation, e.g. the population size at the previous generation.
        '''
        assert type(Nt) not in [type(()), type([])]
        assert len(curXt) == self.nLoci
        #
        # given x(t)
        # calculate y=x(t-1)' by solving an equation
        #
        # x_t = y(1+s2 y+s1 (1-y))/(1+s2 y+2 s1 Y(1-y))
        it = []
        xt = []
        for loc in range(self.nLoci):
            x = curXt[loc]
            # if current allele freq in subpop sp at locus loc has already been 0 or 1,
            # set it to be 0 or 1 for previous gens
            if x in [0, 1]:
                xt.append(x)
                continue
            # In the interaction case, s1, s2 will be different
            # from subpopulation to subpopulation.
            s1 = s[3 * loc + 1]
            s2 = s[3 * loc + 2]
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
            it = GetRNG().randBinomial(int(2 * Nt), y)
            xt.append(float(it) / (2 * Nt))
        return xt

    def _simuForward(self, freq, endFreq, beginGen, endGen):
        '''Simulates a trajectory froward in time, starting from frequency
        ``freq`` at generation ``beginGen`` to frequency ranges specified in
        ``endFreq`` at generation ``endGen``. During the evolution, multiple
        subpopulations can be merged into one, and one population can be split
        into several subpopulations. The number of subpopulation is determined
        by the demographic function. The function raises an exception if the 
        simulated trajectory does not fall into ``destFeq``.
        '''
        # initialize a trajectory
        # freq is assumed to be at the beginning of the beginGen.
        # so we do not it does not count into the trajectory.
        xt = trajectory(endGen = endGen, nLoci = self.nLoci)
        # go through each generation
        for gen in range(beginGen, endGen + 1):
            # first get beginXt, N(t+1), then calculate nextXt.
            if gen == beginGen:
                beginXt = freq
            else:
                # current Xt is the frequency at the previous generation.
                beginXt = xt._freq(gen - 1)
            # _Ne(gen) is the population size at the end of this generation.
            Nt = self._Nt(gen)
            #
            if len(Nt) > len(beginXt):
                # split (forward sense) from one population to nSP subpopulations
                if len(beginXt) != 1:
                    raise exceptions.RuntimeError('Can only split from one subpopulation.')
                #
                # get NextXt using one subpopulation, then split...
                tmpXt = self._getNextXt(beginXt[0], sum(Nt), self._getS(gen, 0, beginXt[0]))
                # split tmpXt to multiple subpopulations and assign endingXt.
                # here we assume a multi-nomial distribution of disease alleels.
                endingXt = [[0]*self.nLoci for x in Nt]
                p = [float(x) / sum(Nt) for x in Nt]
                for loc in range(self.nLoci):
                    it = GetRNG().randMultinomial(int(tmpXt[loc]*sum(Nt)), p)
                    for sp in range(len(Nt)):
                        endingXt[sp][loc] = float(it[sp]) / Nt[sp]
            elif len(Nt) < len(beginXt):
                # check length of next Nt.
                if len(Nt) != 1:
                    raise exceptions.RuntimeError('Can only merge into one subpopulation')
                # merge (forward sense) from multiple subpopulations to one pop.
                Nt_prev = self._Nt(gen - 1)
                if len(beginXt) != len(Nt_prev):
                    raise exceptions.RuntimeError('Subpopulation size and allele frequency mismatch.')
                Nt_tmp = [int(x * 1.0 / sum(Nt_prev) * Nt[0]) for x in Nt_prev]
                #
                endingXt = [[0] * self.nLoci]
                for sp in range(len(Nt_prev)):
                    # simulate frequency in each subpopulation
                    tmpXt = self._getNextXt(beginXt[sp], Nt_tmp[sp], self._getS(gen, sp, beginXt[sp]))
                    # and accumulate alleles in the final merged frequency
                    for loc in range(self.nLoci):
                        endingXt[0][loc] += tmpXt[loc] * Nt_tmp[sp] / Nt[0]
            else:
                endingXt = [self._getNextXt(beginXt[sp], Nt[sp], self._getS(gen, sp, beginXt[sp]))
                    for sp in range(len(Nt))]
            #
            assert len(endingXt) == len(Nt)
            # set frequency at the end of this generation
            if self.logger is not None:
                self.logger.debug('Gen=%d, xt=%s'  % (gen, endingXt))           
            xt._setFreq(endingXt, gen)
            # and now we go to the next generation...
        # not we have a trajectory... is it valid?
        freq = xt._freq(endGen)
        Nt = self._Nt(endGen)
        for loc in range(self.nLoci):
            # case 1: allele frequency at each subpopulation
            if len(endFreq) == self.nLoci * len(Nt):
                for sp in range(len(Nt)):
                    if freq[sp][loc] < endFreq[sp * self.nLoci + loc][0] or \
                        freq[sp][loc] > endFreq[sp * self.nLoci + loc][1]:
                        if self.logger is not None:
                            self.logger.info('Forward trajectory restarted, hitting allele requency %s' % freq)
                        return None
            # case 2: combined allele frequency
            else:
                allFreq = 0
                for sp in range(len(Nt)):
                    allFreq += freq[sp][loc] * Nt[sp]
                allFreq /= sum(Nt)
                if allFreq < endFreq[loc][0] or allFreq > endFreq[loc][1]:
                    if self.logger is not None:
                        self.logger.info('Forward trajectory restarted, hitting allele frequency %s (combined %.3f)' \
                            % (freq, allFreq))
                    return None
        if self.logger is not None:
            self.logger.info('Forward trajectory succeed, hitting allele frequency %s' % freq)
        return xt

    def _simuBackward(self, endGen, freq, minMutAge, maxMutAge):
        '''Simulates a trajectory backward from allele frequency ``freq`` at
        generation ``endGen``. During the evolution, multiple subpopulations can
        be merged into one, and one population can be split into several
        subpopulations. The number of subpopulation is determined by the
        demographic function. If a simulated trajectory is shorter than
        ``minMutAge`` or is longer than ``maxMutAge``, the function will raise
        an exception.
        '''
        if endGen <= 0:
            raise exceptions.ValueError("A positive ending generation is needed.")
        # done[i] is used to track at which generation a trajectory
        # is successfully generated at locus i.
        done = [False] * self.nLoci
        # initialize a trajectory
        xt = trajectory(endGen=endGen, nLoci = self.nLoci)
        # because freq is the allele frequency at the end of the last generation,
        # it is part of the trajectory.
        xt._setFreq(freq, gen=endGen)
        # start from endGen, go backward.
        for gen in range(endGen, -1, -1):
            # first get curXt, N(t-1), then calculate prevXt
            endingXt = xt._freq(gen)
            # Nt is the size at the beginning of the current generation.
            Nt = self._Nt(gen - 1)
            Nt_end = self._Nt(gen)
            if len(Nt) > len(endingXt):
                if len(endingXt) != 1:
                    raise exceptions.RuntimeError('Can only merge to one subpopulation')
                # merge (forward sense)
                tmpXt = self._getPrevXt(endingXt[0], sum(Nt), self._getS(gen, 0, endingXt[0]))
                assert len(tmpXt) == self.nLoci
                # SPLIT tmpXt to multiple subpopulations and assign an expanded beginXt.
                beginXt = [[0]*self.nLoci for x in Nt]
                p = [float(x)/sum(Nt) for x in Nt]
                for loc in range(self.nLoci):
                    it = GetRNG().randMultinomial(int(tmpXt[loc]*sum(Nt)), p)
                    beginXt[sp][loc] = float(it[sp]) / Nt[sp]
            elif len(Nt) < len(endingXt):
                # check length of previous Nt.
                if len(Nt) != 1:
                    raise exceptions.ValueError('Can only split from one subpoplation.')
                # split (forward sense)
                Nt_tmp = [int(float(x) / sum(Nt_end) * Nt[0]) for x in Nt_end]
                beginXt = [[0] * self.nLoci]
                for sp in range(len(Nt)):
                    tmpXt = self._getPrevXt(endingXt[sp], Nt_tmp[sp], self._getS(gen, sp, endingXt[sp]))
                    for loc in range(self.nLoci):
                        beginXt[0][loc] += tmpXt[loc] * Nt_tmp[sp] / Nt[0]
            else:
                beginXt = [self._getPrevXt(endingXt[sp], Nt[sp], self._getS(gen, sp, endingXt[sp]))
                    for sp in range(len(Nt))]
            #
            assert len(beginXt) == len(Nt)
            # set frequency at the end of this generation
            if self.logger is not None:
                self.logger.debug('Gen=%d, xt=%s'  % (gen - 1, beginXt)) 
            #
            xt._setFreq(beginXt, gen - 1)
            # check all loci and see if beginXt is 0
            for loc in range(self.nLoci):
                doneSP = [False] * len(Nt)
                if done[loc]:
                    continue
                # loop over subpopulation
                for sp in range(len(Nt)):
                    if (len(Nt_end) == 1 and endingXt[0][loc] == 0.) or \
                        (len(Nt_end) > 1 and endingXt[sp][loc] == 0.):
                        # already done in a previous generation
                        doneSP[sp] = True
                        continue
                    if beginXt[sp][loc] == 0.:
                        # success
                        doneSP[sp] = True
                        if endGen - gen < minMutAge:
                            if self.logger is not None:
                                self.logger.info('Backward failed - trajectory too short. gen = %d subPop=%d locus = %d' \
                                    % (gen, sp, loc))
                            return None
                        if self.logger is not None:
                            self.logger.info('Backward success: gen = %d subPop=%d locus = %d' % (gen, sp, loc))
                        break
                    elif beginXt[sp][loc] == 1: # fixed
                        if self.logger is not None:
                            self.logger.info('Backward failed - allele gets fixed. gen = %d subPop=%d locus = %d' \
                                    % (gen, sp, loc))
                        return None
                if False not in doneSP:
                    done[loc] = True
            if False not in done:
                # success
                if self.logger is not None:
                    self.logger.info('Backward trajectory succeded at gen = %d' % gen)
                return xt
            # go back gen == 0 and not successful, or if the trajectory is too long
            if gen == 0 or gen + self.maxMutAge < endGen:
                if self.logger is not None:
                    self.logger.info('Backward failed - trajectory too long. gen = %d' % gen)
                return None
            
    def simuForward(self, beginGen, endGen, beginFreq, endFreq, maxAttempts = 10000):
        '''Simulate trajectories of multiple disease susceptibility loci using a
        forward time approach. This function accepts allele frequencies of
        alleles of multiple unlinked loci at the beginning generation (``freq``)
        at generation ``beginGen``, and expected *range* of allele frequencies
        of these alleles (``endFreq``) at the end of generation ``endGen``.
        Depending on the number of loci and subpopulations, these parameters
        accept the following inputs:
        
        beginGen
            Starting generation. The initial frequecies are considered as
            frequencies at the *beginning* of this generation.

        endGen
            Ending generation. The ending frequencies are considerd as
            frequencies at the *end* of this generation.
            
        beginFreq
            The initial allele frequency of involved loci in all subpopulations.
            It can be a number (same frequency for all loci in all
            subpopulations), or a list of frequencies for each locus (same
            frequency in all subpopulations), or a list of frequencies for each
            locus in each subpopulation in the order of ``loc0_sp0``,
            ``loc1_sp0``, ..., ``loc0_sp1``, ``loc1_sp1``, ... and so on.

        endFreq
            The range of acceptable allele frequencies at the ending generation.
            The ranges can be specified for all loci in all subpopulations,
            for all loci (allele frequency in the whole population is
            considered), or for all loci in all subpopulations, in the order
            of ``loc0_sp0``, ``loc1_sp0``, .... ``loc0_sp1``, ... and so on.

        This simulator will simulate a trajectory generation by generation and
        restart if the resulting frequencies do not fall into specified range
        of frequencies. This simulator will return ``None`` if no valid
        trajectory is found after ``maxAttempts`` attemps.
        '''
        #
        # This functin wraps around _simuForward. It handles parameter
        # validation and maxAttempts.
        #
        # endGen
        if not beginGen <= endGen or endGen <= 0:
            raise exceptions.ValueError('Beginning generation should be less than ending generation')
        # beginFreq
        if type(beginFreq) in [type(0), type(0.)]:
            freq = [[beginFreq] * self.nLoci for sp in self._Nt(beginGen)]
        elif type(beginFreq) in [type(()), type([])]:
            if len(beginFreq) == self.nLoci:
                freq = [beginFreq for sp in self._Nt(beginGen)]
            elif len(beginFreq) == self.nLoci * len(self._Nt(beginGen)):
                freq = []
                for sp in range(len(self._Nt(endGen))):
                    freq.append(beginFreq[self.nLoci*sp : self.nLoci * (sp+1)])
            else:
                raise exceptions.ValueError("Invalid initial frequency list")
        else:
            raise exceptions.ValueError("Invalid initial frequency list")
        #
        # endFreq
        if type(endFreq) not in [type(()), type([])] or len(endFreq) == 0:
            raise exceptions.ValueError('A list of frequency range is expected')
        elif type(endFreq[0]) not in [type(()), type([])]:
            if len(endFreq) == 2:
                endFreq = [endFreq]
            else:
                raise exceptions.ValueError('A list of frequency range is expected.')
        if len(endFreq) not in [self.nLoci, self.nLoci * len(self._Nt(endGen))]:
            raise exceptions.ValueError('Please specify a frequency range for each locus')
        for rng in endFreq:
            if len(rng) != 2:
                raise exceptions.ValueError('Please specify frequency range of each marker')
            if rng[0] >= rng[1]:
                raise exceptions.ValueError('Invalid frequency range %s' % rng)
        for failedCount in range(maxAttempts):
            xt = self._simuForward(freq, endFreq, beginGen, endGen)
            if xt is not None:
                return xt
        if self.logger is not None:
            self.logger.info('Simulation failed after %d attempts' % failedCount)
        return None
    
    def simuBackward(self, endGen, endFreq, minMutAge = None, maxMutAge = None,
        maxAttempts = 1000):
        '''Simulate trajectories of multiple disease susceptibility loci using
        a forward time approach. This function accepts allele frequencies of
        alleles of multiple unlinked loci (*endFreq*) at the end of generation
        *endGen*. Depending on the number of loci and subpopulations, parameter
        *beginFreq* can be a number (same frequency for all loci in all
        subpopulations), or a list of frequencies for each locus (same
        frequency in all subpopulations), or a list of frequencies for each
        locus in each subpopulation in the order of ``loc0_sp0``, ``loc1_sp0``,
        ..., ``loc0_sp1``, ``loc1_sp1``, ... and so on.

        This simulator will simulate a trajectory generation by generation and
        restart if the disease allele got fixed (instead of lost), or if the 
        length simulated trajectory does not fall into *minMutAge* and
        *maxMutAge* (ignored if ``None`` is given). This simulator will return
        ``None`` if no valid trajectory is found after ``maxAttempts`` attemps.
        '''
        #
        # This functin wraps around _simuBackward. It handles parameter
        # validation and maxAttempts.
        #
        if endGen <= 0:
            raise exceptions.ValueError('A positive ending generation number is expected.')

        if minMutAge is not None and minMutAge > endGen:
            raise exceptions.ValueError('Minimal mutation age is larger than ending generation.')
        #
        if minMutAge is None:
            self.minMutAge = 0
        else:
            self.minMutAge = minMutAge
        #
        if maxMutAge is None:
            self.maxMutAge = endGen
        else:
            self.maxMutAge = maxMutAge
        
        if not self.maxMutAge >= self.minMutAge:
            raise exceptions.ValueError('maxMutAge should >= minMutAge')
        if endGen == 0 and (callable(self.N) or callable(self.fitness)):
            raise exceptions.ValueError('endGen should be > 0 if N or fitness is defined in the form of function')
        if endGen > 0 and endGen < self.maxMutAge:
            raise exceptions.ValueError('endGen should be >= maxMutAge')
        #
        # endFreq
        if type(endFreq) in [type(0), type(0.)]:
            freq = [[endFreq] * self.nLoci for sp in self._Nt(endGen)]
        elif type(endFreq) in [type(()), type([])]:
            if len(endFreq) == self.nLoci:
                freq = [endFreq for sp in self._Nt(endGen)]
            elif len(endFreq) == self.nLoci * len(self._Nt(endGen)):
                freq = []
                for sp in range(len(self._Nt(endGen))):
                    freq.append(endFreq[self.nLoci*sp : self.nLoci * (sp+1)])
            else:
                raise exceptions.ValueError("Invalid ending frequency list")
        else:
            raise exceptions.ValueError("Invalid ending frequency list")
        #
        for failedCount in range(maxAttempts):
            xt = self._simuBackward(endGen, freq, minMutAge, maxMutAge)
            if xt is not None:
                return xt
        if self.logger is not None:
            self.logger.info('Simulation failed after %d attempts' % failedCount)
        return None


def ForwardTrajectory(N, beginGen, endGen, beginFreq, endFreq, nLoci = 1,
        fitness = None, maxAttempts = 10000, logger=None):
    '''Given a demographic model (*N*) and the fitness of genotype at one or
    more loci (*fitness*), this function simulates a trajectory of one or more
    unlinked loci (*nLoci*) from allele frequency *freq* at generation
    *beginGen* forward in time, until it reaches generation *endGen*. A
    ``trajectory`` object will be returned if the allele frequency falls
    into specified ranges (*endFreq*). ``None`` will be returned if no valid
    trajectory is simulated after ``maxAttempts`` attempts. Please refer to
    class ``trajectory``, ``trajectorySimulator`` and their member functions
    for more details about allowed input for these parameters.
    '''
    return trajectorySimulator(N, nLoci, fitness, logger).simuForward(
        beginGen, endGen, beginFreq, endFreq, maxAttempts)
    
def BackwardTrajectory(N, endGen, endFreq, nLoci=1, fitness=None,
        minMutAge = None, maxMutAge = None, maxAttempts = 1000, logger=None):
    '''Given a demographic model (*N*) and the fitness of genotype at one or
    more loci (*fitness*), this function simulates a trajectory of one or more
    unlinked loci (*nLoci*) from allele frequency *freq* at generation *endGen*
    backward in time, until all alleles get lost. A ``trajectory`` object will
    be returned if the length of simulated trajectory with ``minMutAge`` and
    ``maxMutAge`` (if specified). ``None`` will be returned if no valid
    trajectory is simulated after ``maxAttempts`` attempts. Please refer to
    class ``trajectory``, ``trajectorySimulator`` and their member functions
    for more details about allowed input for these parameters.
    '''
    return trajectorySimulator(N, nLoci, fitness, logger).simuBackward(
        endGen, endFreq, minMutAge, maxMutAge, maxAttempts)


if __name__ == "__main__":
    pass

