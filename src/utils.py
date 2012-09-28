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
# Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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
    'viewVars',
    'migrIslandRates',
    'migrHierarchicalIslandRates',
    'migrSteppingStoneRates',
    'saveCSV',
    'Exporter',
    'export',
    'importPopulation',
    'ProgressBar',
    'Trajectory',
    'TrajectorySimulator',
    'simulateBackwardTrajectory',
    'simulateForwardTrajectory',
]

import sys
import time

from simuOpt import simuOptions

from simuPOP import moduleInfo, MALE, FEMALE, Population, PointMutator, getRNG, ALL_AVAIL, PyOperator

def viewVars(var, gui=None):
    '''
    list a variable in tree format, either in text format or in a
        wxPython window.

    var
        A dictionary variable to be viewed. Dictionary wrapper objects returned
        by ``Population.dvars()`` and ``Simulator.dvars()`` are also acceptable.

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
    if gui in [False, 'batch', 'interactive', 'Tkinter']:
        import pprint
        try:
            # a dvars() object
            pprint.pprint(var.__dict__)
        except:
            pprint.pprint(var)
        return

    try:
        import wx, wx.py.filling as fill
    except ImportError:
        import pprint
        pprint.pprint(var)
        return

    app = wx.PySimpleApp()
    wx.InitAllImageHandlers()
    if var==None:
        fillFrame = fill.FillingFrame()
    else:
        try:
            # a dvars() object?
            fillFrame = fill.FillingFrame(rootObject=var.__dict__,
                rootLabel='var')
        except:
            fillFrame = fill.FillingFrame(rootObject=var,
                rootLabel='var')
    fillFrame.Show(True)
    app.MainLoop()


# migration rate matrix generators
def migrIslandRates(r, n):
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


def migrHierarchicalIslandRates(r1, r2, n):
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
    outside of the group is r2. migrate rate to a specific island depends on
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


def migrSteppingStoneRates(r, n, circular=False):
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

    This function returns [[1]] when there is only one subpopulation.
    '''
    if n < 2:
        return [[1]]
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


def saveCSV(pop, filename='', infoFields=[], loci=ALL_AVAIL, header=True,
        subPops=ALL_AVAIL, genoFormatter=None, infoFormatter=None,
        sexFormatter={MALE: 'M', FEMALE: 'F'},
        affectionFormatter={True: 'A', False: 'U'}, sep=', ', **kwargs):
    '''This function is deprecated. Please use ``export(format='csv')`` instead.
    Save a simuPOP population ``pop`` in csv format. Columns of this
    file is arranged in the order of information fields (``infoFields``),
    sex (if ``sexFormatter`` is not ``None``), affection status (if
    ``affectionFormatter`` is not ``None``), and genotype (if ``genoFormatter`` is
    not ``None``). This function only output individuals in the present
    generation of population ``pop``. This function accepts the following
    parameters:

    pop
        A simuPOP population object.

    filename
        Output filename. Leading '>' characters are ignored. However, if the first
        character of this filename is '!', the rest of the name will be evalulated
        in the population's local namespace. If ``filename`` is empty, the content
        will be written to the standard output.

    infoFileds
        Information fields to be outputted. Default to none.

    loci
        If a list of loci is given, only genotype at these loci will be
        written. Default to ``ALL_AVAIL``, meaning all available loci. You can
        set this parameter to ``[]`` if you do not want to output any genotype.

    header
        Whether or not a header should be written. These headers will include
        information fields, sex (if ``sexFormatter`` is not ``None``), affection
        status (if ``affectionFormatter`` is not ``None``) and loci names. If
        genotype at a locus needs more than one column, ``_1``, ``_2`` etc will
        be appended to loci names. Alternatively, a complete header (a string)
        or a list of column names could be specified directly.

    subPops
        A list of (virtual) subpopulations. If specified, only individuals
        from these subpopulations will be outputed.

    infoFormatter
        A format string that is used to format all information fields. If
        unspecified, ``str(value)`` will be used for each information field.

    genoFormatter
        How to output genotype at specified loci. Acceptable values include
        ``None`` (output allele names), a dictionary with genotype as keys,
        (e.g. ``genoFormatter={(0,0):1, (0,1):2, (1,0):2, (1,1):3}``, or a function
        with genotype (as a tuple of integers) as inputs. The dictionary value
        or the return value of this function can be a single or a list of
        number or strings.

    sexFormatter
        How to output individual sex. Acceptable values include ``None`` (no
        output) or a dictionary with keys ``MALE`` and ``FEMALE``.

    affectionFormatter
        How to output individual affection status. Acceptable values include
        ``None`` (no output) or a dictionary with keys ``True`` and ``False``.

    Parameters ``genoCode``, ``sexCode``, and ``affectionCode`` from version
    1.0.0 have been renamed to ``genoFormatter``, ``sexFormatter`` and 
    ``affectionFormatter`` but can still be used.
    '''
    if moduleInfo()['debug']['DBG_COMPATIBILITY']:
        print >> sys.stderr, 'WARNING: Function saveCSV is deprecated. Use export(format="csv") instead.'
    # handle obsolete parameters affectionCode, sexCode and genoCode
    if kwargs.has_key('genoCode'):
        if moduleInfo()['debug']['DBG_COMPATIBILITY']:
            print >> sys.stderr, 'WARNING: Parameter genoCode is obsolete. Use genoFormatter instead.'
        genoFormatter = kwargs['genoCode']
    if kwargs.has_key('sexCode'):
        if moduleInfo()['debug']['DBG_COMPATIBILITY']:
            print >> sys.stderr, 'WARNING: Parameter sexCode is obsolete. Use sexFormatter instead.'
        sexFormatter = kwargs['sexCode']
    if kwargs.has_key('affectionCode'):
        if moduleInfo()['debug']['DBG_COMPATIBILITY']:
            print >> sys.stderr, 'WARNING: Parameter genoCode is obsolete. Use sexFormatter instead.'
        affectionFormatter = kwargs['affectionCode']
    for key in kwargs.keys():
        if key not in ('genoCode', 'sexCode', 'affectionCode'):
            raise ValueError("Unrecognized keyword parameter %s" % key)
    # parameter pop
    if not isinstance(pop, Population):
        raise ValueError("Passed population should either be a population object")
    # parameter loci
    if loci == ALL_AVAIL:
        loci = range(0, pop.totNumLoci())
    elif type(loci) == type(1):
        loci = [loci]
    if not type(loci) in [type([]) or type(())]:
        raise ValueError("Passed loci should be ALL_AVAIL or a list of loci.")
    # parameter infoFields (allow single input)
    if type(infoFields) == type(''):
        infoFields = [infoFields]
    # parameter filename
    if filename.startswith('!'):
        filename = str(pop.evalulate(filename[1:]))
    if filename.startswith('>'):
        filename = filename.lstrip('>')
    #
    try:
        if filename:
            out = open(filename, "w")
        else:
            out = sys.stdout
    except IOError:
        raise IOError, "Can not open file " + filename +" to write."
    # parameter subPops
    if subPops == ALL_AVAIL:
        subPops = range(pop.numSubPop())
    #
    # figure out columns per genotype
    ploidy = pop.ploidy()
    colPerGenotype = 0
    if len(loci) > 0 and pop.totNumLoci() > 0 and pop.popSize() > 0:
        if genoFormatter is None:
            value = [0]*ploidy
        elif isinstance(genoFormatter, dict):
            if len(genoFormatter) == 0:
                raise ValueError("genoFormatter cannot be empty")
            value = genoFormatter.values()[0]
        else:
            if not callable(genoFormatter):
                raise ValueError("genoFormatter should be a None, a dictionary or a callable function")
            value = genoFormatter(tuple([pop.individual(0).allele(0, p) for p in range(ploidy)]))
        try:
            if type(value) == type(''):
                colPerGenotype = 1
            else:  # a sequece?
                colPerGenotype = len(value)
        except:
            colPerGenotype = 1
    # header
    if header is True:
        names = [x for x in infoFields]
        if sexFormatter is not None:
            names.append('sex')
        if affectionFormatter is not None:
            names.append('aff')
        if colPerGenotype == 1:
            names.extend([pop.locusName(loc) for loc in loci])
        elif colPerGenotype > 1:
            for loc in loci:
                names.extend(['%s_%d' % (pop.locusName(loc), x+1) for x in range(colPerGenotype)])
        # output header
        print >> out, sep.join(names)
    elif type(header) == type(''):
        print >> out, header
    elif type(header) in [type(()), type([])]:
        print >> out, sep.join(header)
    for subPop in subPops:
        for ind in pop.individuals(subPop):
            # information fields
            if infoFormatter is None:
                values = [str(ind.info(x)) for x in infoFields]
            elif type(infoFormatter) == type(''):
                values = [infoFormatter % tuple([ind.info(x) for x in infoFields])]
            else:
                raise ValueError('Parameter infoFormatter can only be None or a format string.')
            # sex
            if sexFormatter is not None:
                values.append(str(sexFormatter[ind.sex()]))
            # affection status
            if affectionFormatter is not None:
                values.append(str(affectionFormatter[ind.affected()]))
            # genotype
            for loc in loci:
                if genoFormatter is None:
                    values.extend([ind.alleleChar(loc, p) for p in range(ploidy)])
                else:
                    genotype = [ind.allele(loc, p) for p in range(ploidy)]
                    if isinstance(genoFormatter, dict):
                        code = genoFormatter[tuple(genotype)]
                    else:
                        code = genoFormatter(genotype)
                    if type(code) in [type([]), type(())]:
                        values.extend(['%s' % x for x in code])
                    else:
                        values.append(str(code))
            # output
            print >> out, sep.join(values)
    # clode output
    if filename:
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
        self.count = 0
        self.percent = 0
        self.completed = False

    def update(self, count=None):
        '''
        Update the progress bar with ``count`` progress. If ``count`` is ``None``,
        it updates by 1 count (not percent).
        '''
        if count is None:
            self.count += 1
        else:
            self.count = min(count, self.totalCount)
        self.progress = int(round(100*self.count/self.totalCount))
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
        ''' Update the progress bar.'''
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
        '''Update the progress bar.'''
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
        '''Update the progreebar.'''
        if not _baseProgressBar.update(self, count):
            return
        self.dialog.Update(self.count, self.message + "\n%d%% completed." % self.progress)
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



class ProgressBar:
    '''The ``ProgressBar`` class defines a progress bar. This class will use a
    text-based progress bar that outputs progressing dots (.) with intermediate
    numbers (e.g. 5 for 50%) under a non-GUI mode (``gui=False``). In the GUI
    mode, a Tkinter or wxPython progress dialog will be used (``gui=Tkinter``
    or ``gui=wxPython``). The default mode is determined by the global gui mode
    of simuPOP (see also ``simuOpt.setOptions``).

    This class is usually used as follows::

        progress = ProgressBar("Start simulation", 500)
        for i in range(500):
            # i+1 can be ignored if the progress bar is updated by 1 step
            progress.update(i+1)   
        # if you would like to make sure the done message is displayed.
        progress.done()
    '''
    def __init__(self, message, totalCount, progressChar='.', block=2, done=' Done.\n', gui=None):
        '''Create a progress bar with ``message``, which will be the title of
        a progress dialog or a message for textbased progress bar. Parameter
        ``totalCount`` specifies total expected steps. If a text-based progress
        bar is used, you could specified progress character and intervals at
        which progresses will be displayed using parameters ``progressChar``
        and ``block``. A ending message will also be displayed in text mode.
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

    def update(self, count=None):
        '''
        Update the progreebar with ``count`` steps done. The dialog or textbar
        may not be updated if it is updated by full percent(s). If ``count`` is
        ``None``, the progressbar increases by one step (not percent).
        '''
        self.progressBar.update(count)

    def done(self):
        '''
        Finish progressbar, print 'done' message if in text-mode.
        '''
        self.progressBar.done()


class Trajectory:
    '''A ``Trajectory`` object contains frequencies of one or more loci in one
    or more subpopulations over several generations. It is usually returned by
    member functions of class ``TrajectorySimulator`` or equivalent global
    functions ``simulateForwardTrajectory`` and ``simulateBackwardTrajectory``.

    The ``Trajectory`` object provides several member functions to facilitate
    the use of Trajectory-simulation techiniques. For example,
    ``Trajectory.func()`` returns a trajectory function that can be provided
    directly to a ``ControlledOffspringGenerator``; ``Trajectory.mutators()``
    provides a list of ``PointMutator`` that insert mutants at the right
    generations to initialize a trajectory.

    For more information about Trajectory simulation techniques and related
    controlled random mating scheme, please refer to the simuPOP user's guide,
    and Peng et al (PLoS Genetics 3(3), 2007).
    '''
    def __init__(self, endGen, nLoci):
        '''Create a ``Trajectory`` object of alleles at *nLoci* loci with
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
        generation *gen* of the simulated Trajectory. Allele frequencies are
        assumed to be zero if *gen* is out of range of the simulated
        Trajectory.
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
        random mating scheme (``ControlledRandomMating``) or a homogeneous
        mating scheme that uses a controlled offspring generator
        (``ControlledOffspringGenerator``).
        '''
        def trajFunc(gen):
            if not self.traj.has_key(gen):
                return [0.] * self.nLoci
            freq = []
            for spFreq in self.traj[gen]:
                freq.extend(spFreq)
            return freq
        return trajFunc

    def mutants(self):
        '''Return a list of mutants in the form of (loc, gen, subPop)'''
        gens = self.traj.keys()
        gens.sort()
        if len(gens) == 0:
            return []
        mut = []
        for gen in gens[:-1]:
            # no introduction of mutants with Population merge or split.
            if len(self.traj[gen]) != len(self.traj[gen + 1]):
                continue
            # we may need to introduce mutant at each subpopulation.
            for sp in range(len(self.traj[gen])):
                for loc in range(self.nLoci):
                    if self.traj[gen][sp][loc] == 0 and self.traj[gen + 1][sp][loc] > 0:
                        mut.append((loc, gen + 1, sp))
        return mut

    
    def mutators(self, loci, inds=0, allele=1, *args, **kwargs):
        '''Return a list of ``PointMutator`` operators that introduce mutants
        at the beginning of simulated trajectories. These mutators should be
        added to the ``preOps`` parameter of ``Simulator.evolve`` function to
        introduce a mutant at the beginning of a generation with zero allele
        frequency before mating, and a positive allele frequency after mating.
        A parameter ``loci`` is needed to specify actual loci indexes in the
        real forward simulation. Other than default parameters ``inds=0`` and
        ``allele=1``, additional parameters could be passed to point mutator
        as keyward parameters.
        '''
        ops = []
        if hasattr(loci, '__iter__') and len(loci) != self.nLoci:
            raise ValueError('%d loci is expected' % self.nLoci)
        for loc, gen, sp in self.mutants():
            if self.nLoci == 1 and type(loci) == type(0):
                ops.append(PointMutator(inds=inds, loci=loci, allele=allele,
                                subPops=sp, at=gen, *args, **kwargs))
            elif hasattr(loci, '__iter__'):
                ops.append(PointMutator(inds=inds, loci=loci[loc], allele=allele,
                                subPops=sp, at=gen, *args, **kwargs))
            else:
                raise ValueError('Invalid value for parameter loci')
        return ops

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
        '''Plot simulated Trajectory using ``R`` through a Python module ``rpy``.
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
        args = plotter.DerivedArgs(
            defaultFuncs=['plot', 'lines'],
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
        plotter.r.dev_off()
        

class TrajectorySimulator:
    '''A Trajectory Simulator takes basic demographic and genetic (natural
    selection) information of an evolutionary process of a diploid population
    and allow the simulation of Trajectory of allele frequencies of one or
    more loci. Trajectories could be simulated in two ways: forward-time and
    backward-time. In a forward-time simulation, the simulation starts from
    certain allele frequency and simulate the frequency at the next generation
    using given demographic and genetic information. The simulation continues
    until an ending generation is reached. A Trajectory is successfully
    simulated if the allele frequency at the ending generation falls into a
    specified range. In a backward-time simulation, the simulation starts from
    the ending generation with a desired allele frequency and simulate the
    allele frequency at previous generations one by one until the allele gets
    lost (allele frequency equals zero).

    The result of a trajectory simulation is a trajectory object which can be
    used to direct the simulation of a special random mating process that
    controls the evolution of one or more disease alleles so that allele 
    frequencies are consistent across replicate simulations. For more
    information about Trajectory simulation techniques and related controlled
    random mating scheme, please refer to the simuPOP user's guide, and Peng et
    al (PLoS Genetics 3(3), 2007).
    '''

    def __init__(self, N, nLoci=1, fitness=None, logger=None):
        '''Create a trajectory Simulator using provided demographic and genetic
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
            diallelic loci. The Trajectory represents frequencies of a 

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
                raise ValueError('Invalid list of fitness.')
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
            # the return value of a demographic function sometimes is not integer.
            if type(nt) in [int, long, float]:
                return [int(nt)]
            else:
                return [int(x) for x in nt]
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
        # case 3: 3**self.nLoci, interaction
        elif len(fit) == 3**self.nLoci:
            # from fitness list, get s using allele frequency
            # Allele frequency for each subpopulation is passed and there will be
            # different s for each subpopulation because different allele frequencies.
            s.extend(self._marginalFitness(fit, freq))
        else:
            raise ValueError('Wrong length of list of fitness: ' + str(len(fit)))
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
            it = getRNG().randBinomial(2 * Nt, y)
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
            # if current allele freq in subpop sp at locus loc has already been 0,
            # it to be 0 for previous gens
            if x == 0:
                xt.append(x)
                continue
            # if current allele freq in subop sp is 1, we assume that it just reaches
            # here by losing one allele
            if x == 1:
                xt.append(float(2 * Nt - 1) / (2 * Nt))
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
            it = getRNG().randBinomial(int(2 * Nt), y)
            xt.append(float(it) / (2 * Nt))
        return xt

    def _simuForward(self, freq, endFreq, beginGen, endGen):
        '''Simulates a trajectory froward in time, starting from frequency
        ``freq`` at generation ``beginGen`` to frequency ranges specified in
        ``endFreq`` at generation ``endGen``. During the evolution, multiple
        subpopulations can be merged into one, and one population can be split
        into several subpopulations. The number of subpopulation is determined
        by the demographic function. The function returns the ending allele
        frequency if the simulated Trajectory does not fall into ``destFeq``,
        and a ``Trajectory`` object otherwise.
        '''
        # initialize a trajectory
        # freq is assumed to be at the beginning of the beginGen.
        # so we do not it does not count into the Trajectory.
        xt = Trajectory(endGen = endGen, nLoci = self.nLoci)
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
                    raise RuntimeError('Can only split from one subpopulation.')
                #
                # get NextXt using one subpopulation, then split...
                tmpXt = self._getNextXt(beginXt[0], sum(Nt), self._getS(gen, 0, beginXt[0]))
                # split tmpXt to multiple subpopulations and assign endingXt.
                # here we assume a multi-nomial distribution of disease alleels.
                endingXt = [[0]*self.nLoci for x in Nt]
                p = [float(x) / sum(Nt) for x in Nt]
                for loc in range(self.nLoci):
                    it = getRNG().randMultinomial(int(tmpXt[loc]*sum(Nt)), p)
                    for sp in range(len(Nt)):
                        endingXt[sp][loc] = float(it[sp]) / Nt[sp]
            elif len(Nt) < len(beginXt):
                # check length of next Nt.
                if len(Nt) != 1:
                    raise RuntimeError('Can only merge into one subpopulation')
                # merge (forward sense) from multiple subpopulations to one pop.
                Nt_prev = self._Nt(gen - 1)
                if len(beginXt) != len(Nt_prev):
                    raise RuntimeError('Subpopulation size and allele frequency mismatch.')
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
            #if self.logger:
            #    self.logger.debug('Gen=%d, xt=%s'  % (gen, endingXt))           
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
                        if self.logger:
                            self.logger.debug('Forward Trajectory restarted, hitting allele requency ' + str(freq))
                        return freq
            # case 2: combined allele frequency
            else:
                allFreq = 0
                for sp in range(len(Nt)):
                    allFreq += freq[sp][loc] * Nt[sp]
                allFreq /= sum(Nt)
                if allFreq < endFreq[loc][0] or allFreq > endFreq[loc][1]:
                    if self.logger:
                        self.logger.debug('Forward Trajectory restarted, hitting allele frequency %s (combined %.3f)' \
                            % (freq, allFreq))
                    return allFreq
        if self.logger:
            self.logger.info('Forward Trajectory succeed, hitting allele frequency ' + str(freq))
        return xt

    def _avgOfNestedList(self, value):
        '''Take average of each element of a nested list of the same shape. For
        example, _avgOfNestedList([[1,[2,3]], [2,[3,4]]]) would return
        [1.5, [2.5, 3.5]]. This is used to return summary statistics of failed
        attempts.
        '''
        if len(value) == 0:
            return []
        if type(value[0]) in [type(()), type([])]:
            avg = []
            for i in range(len(value[0])):
                avg.append(self._avgOfNestedList([val[i] for val in value]))
        else:
            return float(sum(value)) / len(value)
        return avg

    def _simuBackward(self, endGen, freq, minMutAge, maxMutAge):
        '''Simulates a trajectory backward from allele frequency ``freq`` at
        generation ``endGen``. During the evolution, multiple subpopulations can
        be merged into one, and one population can be split into several
        subpopulations. The number of subpopulation is determined by the
        demographic function. If a simulated Trajectory is shorter than
        ``minMutAge`` or is longer than ``maxMutAge``, the function will raise
        an exception.
        '''
        if endGen <= 0:
            raise ValueError("A positive ending generation is needed.")
        # done[i] is used to track at which generation a trajectory
        # is successfully generated at locus i.
        done = [False] * self.nLoci
        # initialize a trajectory
        xt = Trajectory(endGen=endGen, nLoci = self.nLoci)
        # because freq is the allele frequency at the end of the last generation,
        # it is part of the Trajectory.
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
                    raise RuntimeError('Can only merge to one subpopulation')
                # merge (forward sense)
                tmpXt = self._getPrevXt(endingXt[0], sum(Nt), self._getS(gen, 0, endingXt[0]))
                assert len(tmpXt) == self.nLoci
                # SPLIT tmpXt to multiple subpopulations and assign an expanded beginXt.
                beginXt = [[0]*self.nLoci for x in Nt]
                p = [float(x)/sum(Nt) for x in Nt]
                for loc in range(self.nLoci):
                    it = getRNG().randMultinomial(int(tmpXt[loc]*sum(Nt)), p)
                    beginXt[sp][loc] = float(it[sp]) / Nt[sp]
            elif len(Nt) < len(endingXt):
                # check length of previous Nt.
                if len(Nt) != 1:
                    raise ValueError('Can only split from one subpoplation.')
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
            if self.logger:
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
                            if self.logger:
                                self.logger.debug('Backward failed - Trajectory too short. gen = %d subPop=%d locus = %d' \
                                    % (gen, sp, loc))
                            return (gen, beginXt)
                        if self.logger:
                            self.logger.debug('Backward success: gen = %d subPop=%d locus = %d' % (gen, sp, loc))
                        break
                    elif beginXt[sp][loc] == 1: # fixed
                        if self.logger:
                            self.logger.debug('Backward failed - allele gets fixed. gen = %d subPop=%d locus = %d' \
                                    % (gen, sp, loc))
                        return (gen, beginXt)
                if False not in doneSP:
                    done[loc] = True
            if False not in done:
                # success
                if self.logger:
                    self.logger.info('Backward Trajectory succeded at gen = %d' % gen)
                return xt
            # go back gen == 0 and not successful, or if the Trajectory is too long
            if gen == 0 or gen + self.maxMutAge < endGen:
                if self.logger:
                    self.logger.debug('Backward failed - Trajectory too long. gen = %d' % gen)
                return (gen, beginXt)
            
    def simuForward(self, beginGen, endGen, beginFreq, endFreq, maxAttempts=10000):
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
        Trajectory is found after ``maxAttempts`` attemps.
        '''
        #
        # This functin wraps around _simuForward. It handles parameter
        # validation and maxAttempts.
        #
        # endGen
        if not beginGen <= endGen or endGen <= 0:
            raise ValueError('Beginning generation should be less than ending generation')
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
                raise ValueError("Initial frequency should be provided for each locus (nLoci) or each locus at each subpopulation (nLoci * len(N)).")
        else:
            raise ValueError("Invalid initial frequency list")
        #
        # endFreq
        if type(endFreq) not in [type(()), type([])] or len(endFreq) == 0:
            raise ValueError('A list of frequency range is expected')
        elif type(endFreq[0]) not in [type(()), type([])]:
            if len(endFreq) == 2:
                endFreq = [endFreq]
            else:
                raise ValueError('A list of frequency range is expected.')
        if len(endFreq) not in [self.nLoci, self.nLoci * len(self._Nt(endGen))]:
            raise ValueError('Please specify a frequency range for each locus')
        for rng in endFreq:
            if len(rng) != 2:
                raise ValueError('Please specify frequency range of each marker')
            if rng[0] > rng[1]:
                raise ValueError('Invalid frequency range %f - %f' % (rng[0], rng[1]))
        failedFreq = []
        for failedCount in range(maxAttempts):
            xt = self._simuForward(freq, endFreq, beginGen, endGen)
            if isinstance(xt, Trajectory):
                if self.logger:
                    self.logger.info('Simulation succeed after %d attempts with average ending frequencies %s.' \
                        % (failedCount, self._avgOfNestedList(failedFreq)))
                return xt
            else:
                failedFreq.append(xt)
        if self.logger:
            self.logger.debug('Ending frequencies:')
            for freq in failedFreq:
                self.logger.debug('    ' + str(freq))
            self.logger.info(('Simulation failed after %d attempts with average frequencies ' % failedCount) \
                + str(self._avgOfNestedList(failedFreq)))
        return None
    
    def simuBackward(self, endGen, endFreq, minMutAge=None, maxMutAge=None,
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
        length simulated Trajectory does not fall into *minMutAge* and
        *maxMutAge* (ignored if ``None`` is given). This simulator will return
        ``None`` if no valid Trajectory is found after ``maxAttempts`` attemps.
        '''
        #
        # This functin wraps around _simuBackward. It handles parameter
        # validation and maxAttempts.
        #
        if endGen <= 0:
            raise ValueError('A positive ending generation number is expected.')

        if minMutAge is not None and minMutAge > endGen:
            raise ValueError('Minimal mutation age is larger than ending generation.')
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
            raise ValueError('maxMutAge should >= minMutAge')
        if endGen == 0 and (callable(self.N) or callable(self.fitness)):
            raise ValueError('endGen should be > 0 if N or fitness is defined in the form of function')
        if endGen > 0 and endGen < self.maxMutAge:
            raise ValueError('endGen should be >= maxMutAge')
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
                raise ValueError("Invalid ending frequency list")
        else:
            raise ValueError("Invalid ending frequency list")
        #
        failedFreq = []
        for failedCount in range(maxAttempts):
            xt = self._simuBackward(endGen, freq, self.minMutAge, self.maxMutAge)
            if isinstance(xt, Trajectory):
                if self.logger:
                    self.logger.info(('Simulation succeeded after %d attempts with average generation and frequencies' \
                        % failedCount) + str(self._avgOfNestedList(failedFreq)))
                return xt
            else:
                failedFreq.append(xt)
        if self.logger:
            self.logger.debug('Beginning generation and frequencies:')
            for freq in failedFreq:
                self.logger.debug('    ' + str(freq))
            self.logger.info(('Simulation failed after %d attempts with average starting generation and frequencies ' % failedCount) \
                + str(self._avgOfNestedList(failedFreq)))
        return None


def simulateForwardTrajectory(N, beginGen, endGen, beginFreq, endFreq, nLoci=1,
        fitness=None, maxAttempts=10000, logger=None):
    '''Given a demographic model (*N*) and the fitness of genotype at one or
    more loci (*fitness*), this function simulates a trajectory of one or more
    unlinked loci (*nLoci*) from allele frequency *freq* at generation
    *beginGen* forward in time, until it reaches generation *endGen*. A
    ``Trajectory`` object will be returned if the allele frequency falls
    into specified ranges (*endFreq*). ``None`` will be returned if no valid
    Trajectory is simulated after ``maxAttempts`` attempts. Please refer to
    class ``Trajectory``, ``TrajectorySimulator`` and their member functions
    for more details about allowed input for these parameters. If a *logger*
    object is given, it will send detailed debug information at ``DEBUG``
    level and ending allele frequencies at the ``INFO`` level. The latter
    can be used to adjust your fitness model and/or ending allele frequency
    if a trajectory is difficult to obtain because of parameter mismatch.
    '''
    return TrajectorySimulator(N, nLoci, fitness, logger).simuForward(
        beginGen, endGen, beginFreq, endFreq, maxAttempts)

def simulateBackwardTrajectory(N, endGen, endFreq, nLoci=1, fitness=None,
        minMutAge=None, maxMutAge=None, maxAttempts=1000, logger=None):
    '''Given a demographic model (*N*) and the fitness of genotype at one or
    more loci (*fitness*), this function simulates a trajectory of one or more
    unlinked loci (*nLoci*) from allele frequency *freq* at generation *endGen*
    backward in time, until all alleles get lost. A ``Trajectory`` object will
    be returned if the length of simulated Trajectory with ``minMutAge`` and
    ``maxMutAge`` (if specified). ``None`` will be returned if no valid
    Trajectory is simulated after ``maxAttempts`` attempts. Please refer to
    class ``Trajectory``, ``TrajectorySimulator`` and their member functions
    for more details about allowed input for these parameters. If a *logger*
    object is given, it will send detailed debug information at ``DEBUG``
    level and ending generation and frequency at the ``INFO`` level. The latter
    can be used to adjust your fitness model and/or ending allele frequency
    if a trajectory is difficult to obtain because of parameter mismatch.
    '''
    return TrajectorySimulator(N, nLoci, fitness, logger).simuBackward(
        endGen, endFreq, minMutAge, maxMutAge, maxAttempts)

#
# STRUCTURE format (no import yet)
#
class StructureExporter:
    '''An exporter to export given population in structure format'''
    def __init__(self, markerNames=True, recessiveAlleles=None, interMarkerDistances=True,
        phaseInformation=None, label=True, popData=True, popFlag=None, locData=None,
        phenotype=None):
        self.markerNames = markerNames
        self.recessiveAlleles = recessiveAlleles
        self.interMarkerDistances = interMarkerDistances
        self.phaseInformation = phaseInformation
        self.label = label
        self.popData = popData
        self.popFlag = popFlag
        self.locData = locData
        self.phenotype = phenotype
    
    def export(self, pop, filename, subPops, infoFields, gui):
        '''export in structure format '''
        # http://pritch.bsd.uchicago.edu/structure_software/release_versions/v2.3.4/structure_doc.pdf
        with open(filename, 'w') as out:
            #
            # first line: marker names
            #
            if self.markerNames is True:
                names = pop.lociNames()
                if names:
                    out.write('\t'.join(names) + '\n')
            elif hasattr(self.markerNames, '__iter__'):
                if len(self.markerNames) != pop.totNumLoci():
                    raise ValueError('%d names are provided for %d markers' % (len(self.markerNames), pop.totNumLoci()))
                out.write('\t'.join(self.markerNames) + '\n')
            else:
                raise ValueError('Please provide a list of marker names for parameter markerNames')
            #
            # second line: recessive alleles
            #
            if self.recessiveAlleles is not None:
                if self.recessiveAlleles not in [0, 1]:
                    raise ValueError('Only 0 or 1 is acceptable for parameter revessiveAlleles')
                out.write('%d\n' % self.recessiveAlleles)
            # 
            # third line: inter marker distance
            #
            if self.interMarkerDistances is True:
                loci_pos = pop.lociPos()
                # get difference
                loci_dist = [-1] + [loci_pos[i] - loci_pos[i-1] for i in range(1, len(loci_pos))]
                # set beginning of each chromosome to -1
                for ch in range(pop.numChrom()):
                    loci_dist[pop.chromBegin(ch)] = -1
                out.write('\t'.join(['%s' % x for x in loci_dist]) + '\n')
            #
            # fourth line: phase information
            #
            if self.phaseInformation is not None:
                if self.phaseInformation not in [0, 1]:
                    raise ValueError('Only 0 or 1 is acceptable for parameter revessiveAlleles')
                out.write('%d\n' % self.phaseInformation)
            # 
            # sixth line and later: genotype lines
            #
            # progress bar might be wrong with subPops parameter...
            prog = ProgressBar(filename, pop.popSize(), gui=gui)
            count = 0
            for vsp in subPops:
                sp = vsp if type(vsp) == type(0) else vsp[0]
                for idx, ind in enumerate(pop.individuals(vsp)):
                    items = []
                    #
                    # label
                    #
                    if self.label:
                        items.append(str(idx + 1))
                    #
                    # popData
                    #
                    if self.popData:
                        items.append(str(sp + 1))
                    #
                    # popFlag
                    #
                    if self.popFlag is not None:
                        if self.popFlag not in [0, 1]:
                            raise ValueError('Only 0 or 1 is acceptable for parameter popFlag')
                        items.append(str(self.popFlag))
                    #
                    # locData
                    #
                    if self.locData is not None:
                        try:
                            items.append(str(int(ind.info(self.locData))))
                        except:
                            raise ValueError('Population does not have information field %s as locData' % self.locData)
                    #
                    # phenotype
                    #
                    if self.phenotype is not None:
                        try:
                            items.append(str(int(ind.info(self.phenotype))))
                        except:
                            raise ValueError('Population does not have information field %s as phenotype' % self.locData)
                    #
                    # genotype
                    #
                    for p in range(pop.ploidy()):
                        if items:
                            out.write('%s\t%s\n' % ('\t'.join(items), '\t'.join([str(x) for x in ind.genotype(p)])))
                        else:
                            out.write('%s\n' % '\t'.join([str(x) for x in ind.genotype(p)]))
                    #
                    # update progress bar
                    #
                    count += 1
                    prog.update(count)
            prog.done()


#
# GenePop format
#
class GenePopExporter:
    '''An exporter to export given population in structure format'''
    def __init__(self, title=None, adjust=1):
        self.title = title.rstrip() if title is not None else None
        self.adjust = adjust
    
    def export(self, pop, filename, subPops, infoFields, gui):
        ''' Export in genepop format '''
        # http://genepop.curtin.edu.au/help_input.html
        if pop.ploidy() != 2:
            raise ValueError('simuPOP currently can only export diploid populations in GenePop format.')
        #
        with open(filename, 'w') as out:
            #
            # first line: title
            #
            if self.title is not None:
                out.write(self.title + '\n')
            else:
                out.write('%s in GenePop format. Outputted by simuPOP at %s\n' % (
                    filename, time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
            #
            # second line: allele names
            #
            names = pop.lociNames()
            if names:
                # if names are specified
                out.write(', '.join(names) + '\n')
            else:
                names = []
                for ch in range(pop.numChrom()):
                    for loc in range(pop.numLoci(ch)):
                        names.append('ch%d-loc%d' % (ch + 1, loc + 1))
                out.write(', '.join(names) + '\n')
            # 
            # output genotype
            #
            # progress bar might be wrong with subPops parameter...
            alleleWidth = 3 if max(pop.genotype()) >= 99 else 2
            format_string = '%%0%dd%%0%dd' % (alleleWidth, alleleWidth)
            prog = ProgressBar(filename, pop.popSize(), gui=gui)
            count = 0
            numLoci = pop.totNumLoci()
            for vsp in subPops:
                # 
                # for each subpopulation, output pop
                #
                out.write('POP\n')
                # the name might contain space etc
                name = ''.join([x for x in pop.subPopName(vsp) if x.isalnum()])
                if not name:
                    name = 'SubPop%d' % (vsp if type(vsp) == type(0) else vsp[0])
                #
                for idx, ind in enumerate(pop.individuals(vsp)):
                    #
                    # label
                    #
                    out.write('%s-%d, ' % (name, idx + 1))
                    #
                    # genotype
                    #
                    geno = ind.genotype()
                    out.write(' '.join([format_string % (geno[x] + self.adjust, geno[numLoci + x] + self.adjust) for x in range(numLoci)]) + '\n')
                    #
                    # update progress bar
                    #
                    count += 1
                    prog.update(count)
            prog.done()


class GenePopImporter:
    def __init__(self, adjust=0):
        self.adjust = adjust

    def importFrom(self, filename):
        with open(filename, 'r') as input:
            #
            # ignore the first line
            #
            input.readline()
            #
            # read all loci names
            #
            loci_names = []
            while True:
                line = input.readline()
                if not line.rstrip():
                    raise ValueError('No POP line is found. This file must not be in GenePop format')
                if line.lower().rstrip() == 'pop':
                    break
                loci_names.extend([x.strip() for x in line.split(',')])
            # 
            # read genotypes
            #
            popSize = [0]
            genotypes = []
            while True:
                line = input.readline()
                if not line.rstrip():
                    break
                # new subpopulation
                if line.lower().rstrip() == 'pop':
                    popSize.append(0)
                    continue
                # increase pop size count
                popSize[-1] = popSize[-1] + 1
                # 
                try:
                    # ignore ind name
                    name, geno = line.split(',', 1)
                    # split genotype into pieces
                    geno = [x.strip() for x in geno.split() if x.strip()]
                    # get alleles (adjusted with self.adjust)
                    alleles = [(int(x[:3]) + self.adjust, int(x[3:]) + self.adjust) if len(x) == 6 else \
                        (int(x[:2]) + self.adjust, int(x[2:]) + self.adjust) for x in geno]
                    # append alleles in simuPOP order
                    genotypes.extend([x[0] for x in alleles])
                    genotypes.extend([x[1] for x in alleles])
                    if len(geno) != len(loci_names):
                        raise ValueError('Incorrect number of genotype (%d expected)' % len(loci_names))
                except Exception as e:
                    raise ValueError('Invalid input genotype line (%s). The file must not be in GenePop format. %s' % (line, e))
            #
            # create a population
            pop = Population(size=popSize, loci = len(loci_names), lociNames=loci_names)
            pop.setGenotype(genotypes)
        return pop

#
# FSTAT format
#
# The first line contains 4 numbers: the number of samples, np , the 
# number of loci, nl, the highest number used to label an allele, nu, 
# and a 1 if the code for alleles is a one digit number (1-9), a 2 if 
# code for alleles is a 2 digit number (01-99) or a 3 if code for 
# alleles is a 3 digit number (001-999). These 4 numbers need to be
# separated by any number of spaces. 
#
# The first line is immediately followed by nl lines, each containing the 
# name of a locus, in the order they will appear in the rest of the file. 
#
# On line nl+2, a series of numbers as follow: 
# 1     0102   0103   0101  0203          0      0303
#
# The first number identifies the sample to which the individual belongs,
# the second is the genotype of the individual at the first locus, coded
# with a 2 digits number for each allele, the third is the genotype at the
# second locus, until locus nl is entered (in the example above, nl=6).
# Missing genotypes are encoded with 0. Note that 0001 or 0100 are not 
# a valid format, that is, both alleles at a locus have to be known, 
# otherwise, the genotype is considered as missing. No empty lines 
# are needed between samples.
#
class FStatExporter:
    '''An exporter to export given population in fstat format'''
    def __init__(self, lociNames=None, adjust=1):
        self.lociNames = lociNames
        self.adjust = adjust
    
    def export(self, pop, filename, subPops, infoFields, gui):
        '''Export in FSTAT format
        '''
        with open(filename, 'w') as out:
            #
            # first line: np, nl, nu and nd
            #
            np = pop.numSubPop()
            nl = pop.totNumLoci()
            nu = max(pop.genotype()) + self.adjust
            if nu < 10:
                nd = 1
            elif nu < 100:
                nd = 2
            elif nu < 1000:
                nd = 3
            else: # FSTAT can not handle this now. how many digits?
                nd = len(str(nu))
            #
            out.write( '%d %d %d %d\n' % (np, nl, nu, nd))
            #
            # loci names
            #
            if self.lociNames:
                if len(self.lociNames) != pop.totNumLoci():
                    raise ValueError('Parameter lociNames, if specified, should give all %d loci a name' % pop.totNumLoci())
                [out.write(x + '\n') for x in self.lociNames]
            else:
                names = pop.lociNames()
                if names:
                    [out.write(x + '\n') for x in names]
                else:
                    # cook up some name
                    for ch in range(pop.numChrom()):
                        for loc in range(pop.numLoci(ch)):
                            out.write('chr%d_%d\n' % (ch, loc))
            #
            #  genotype
            #
            format_string = '%%0%dd%%0%dd' % (nd, nd)
            numLoci = pop.totNumLoci()
            prog = ProgressBar(filename, pop.popSize(), gui=gui)
            count = 0
            for vsp in subPops:
                sp = vsp if type(vsp) == type(0) else vsp[0]
                for ind in pop.individuals(vsp):
                    geno = ind.genotype()
                    out.write("%d " % (sp + 1) + ' '.join([format_string % (geno[x] + self.adjust, geno[numLoci + x] + self.adjust) for x in range(numLoci)]) + '\n')
                    count += 1
                    prog.update(count)
            prog.done()


class FStatImporter:
    def __init__(self, adjust=0):
        self.adjust = adjust

    def importFrom(self, filename):
        with open(filename, 'r') as input:
            # file is opened. get basic parameters
            try:
                # get numSubPop(), totNumLoci(), maxAllele(), digit
                [np, nl, nu, nd] = map(int, input.readline().split())
            except ValueError:
                raise ValueError("The first line does not have 4 numbers. Are you sure this is a FSTAT file?")
            # now, ignore nl lines, if loci is empty try to see if we have info here
            # following lines with loci name.
            lociNames = []
            for al in range(nl):
                lociNames.append(input.readline().strip())
            #
            # get all the genotypes
            subPopIndex = []
            genotypes = []
            for line in input.readlines():
                try:
                    items = line.split()
                    if len(items) != nl + 1:
                        raise ValueError('Genotype line (%s) has incorrect number of items' % line)
                    subPopIndex.append(int(items[0]))
                    # 
                    # split genotype into pieces
                    geno = [x.strip() for x in items[1:]]
                    # get alleles (adjusted with self.adjust)
                    alleles = [(int(x[:nd]) + self.adjust, int(x[nd:]) + self.adjust) if x != '0' else (self.adjust, self.adjust) for x in geno]
                    # append alleles in simuPOP order
                    genotypes.extend([x[0] for x in alleles])
                    genotypes.extend([x[1] for x in alleles])
                    if len(geno) != nl:
                        raise ValueError('Incorrect number of genotype (%d expected)' % len(loci_names))
                except Exception as e:
                    raise ValueError('Invalid input genotype line (%s). The file must not be in FSTAT format. %s' % (line, e))
                    
        # subpop size?
        # count number of subpopulations
        subPopSize = [0] * (max(subPopIndex) + 1)
        for idx in subPopIndex:
            subPopSize[idx] += 1
        if len([x for x in subPopSize if x != 0]) != np:
            raise ValueError("Number of subpop does not match")
        # we have all the information, create a population
        pop = Population(size=[x for x in subPopSize if x != 0],
            subPopNames=[str(idx) for idx,x in enumerate(subPopSize) if x != 0],
            loci=len(lociNames), lociNames=lociNames)
        # set genotype
        pop.setGenotype(genotypes)
        return pop

#
# Format CSV
#
class CSVExporter:
    '''An exporter to export given population in csv format'''
    def __init__(self, header=True, genoFormatter=None, infoFormatter=None,
        sexFormatter={MALE: 'M', FEMALE: 'F'},
        affectionFormatter={True: 'A', False: 'U'}, delimiter=','):
        self.header = header
        self.genoFormatter = genoFormatter
        self.infoFormatter = infoFormatter
        self.sexFormatter = sexFormatter
        self.affectionFormatter = affectionFormatter
        self.delimiter = delimiter
    
    def _genoFromDict(self, geno):
        return self.genoFormatter[geno]

    def _genoDirect(self, geno):
        return geno

    def _genoCallable(self, geno):
        return self.genoFormatter(geno)

    def export(self, pop, filename, subPops, infoFields, gui):
        '''Export in FSTAT format
        '''
        ploidy = pop.ploidy()
        colPerGenotype = 0
        if pop.totNumLoci() > 0 and pop.popSize() > 0:
            if self.genoFormatter is None:
                _genoFunc = self._genoDirect
                colPerGenotype = ploidy
            elif isinstance(self.genoFormatter, dict):
                value = self.genoFormatter.values()[0]
                colPerGenotype = 1 if type(value) in [type(''), type(1), type(1L)] else len(value)
                _genoFunc = self._genoFromDict
            else:
                if not callable(self.genoFormatter):
                    raise ValueError("genoFormatter should be a None, a dictionary or a callable function")
                value = self.genoFormatter(tuple([pop.individual(0).allele(0, p) for p in range(ploidy)]))
                colPerGenotype = 1 if type(value) in [type(''), type(1), type(1L)] else len(value)
                _genoFunc = self._genoCallable
            print colPerGenotype, ploidy
        #
        with open(filename, 'w') as out:
            # header
            if self.header is True:
                names = list(infoFields)
                if self.sexFormatter is not None:
                    names.append('sex')
                if self.affectionFormatter is not None:
                    names.append('aff')
                if colPerGenotype == 1:
                    names.extend([pop.locusName(loc) for loc in range(pop.totNumLoci())])
                elif colPerGenotype > 1:
                    for loc in range(pop.totNumLoci()):
                        names.extend(['%s_%d' % (pop.locusName(loc), x+1) for x in range(colPerGenotype)])
                # output header
                out.write(self.delimiter.join(names) + '\n')
            elif type(self.header) == type(''):
                out.write(self.header + '\n')
            elif type(self.header) in [type(()), type([])]:
                out.write(self.delimiter.join([str(x) for x in self.header]) + '\n')
            # progress bar
            prog = ProgressBar(filename, pop.popSize(), gui=gui)
            count = 0
            for vsp in subPops:
                for ind in pop.individuals(vsp):
                    # information fields
                    if self.infoFormatter is None:
                        values = [str(ind.info(x)) for x in infoFields]
                    elif type(self.infoFormatter) == type(''):
                        values = [self.infoFormatter % tuple([ind.info(x) for x in infoFields])]
                    else:
                        raise ValueError('Parameter infoFormatter can only be None or a format string.')
                    # sex
                    if self.sexFormatter is not None:
                        values.append(str(self.sexFormatter[ind.sex()]))
                    # affection status
                    if self.affectionFormatter is not None:
                        values.append(str(self.affectionFormatter[ind.affected()]))
                    # genotype
                    for geno in zip(*[ind.genotype(p) for p in range(ploidy)]):
                        val = _genoFunc(geno)
                        if type(val) in [type([]), type(())]:
                            values.extend(['%s' % x for x in val])
                        else:
                            values.append(str(val))
                    # output
                    out.write(self.delimiter.join(values) + '\n')
                    count += 1
                    prog.update(count)
            # clode output
            if filename:
                out.close()
            prog.done()


#
# Wrapper of the importer and expoters.
#
# Class Exporter
# function export
# function importPopulation
#


class Exporter(PyOperator):
    '''An operator to export the current population in specified format.
    Currently supported file formats include:

    STRUCTURE (http://pritch.bsd.uchicago.edu/structure.html). This format
    accepts the following parameters:

    markerNames
        If set to True (default), output names of loci that are specified by parameter
        *lociNames* of the ``Population`` class. No names will be outputted if loci are
        anonymous. A list of loci names are acceptable which will be outputted directly.
    
    recessiveAlleles
        If specified, value of this parameter will be outputted after the marker names
        line.

    interMarkerDistances
        If set to True (default), output distances between markers. The first marker
        of each chromosome has distance -1, as required by this format.

    phaseInformation
        If specified, output the value (0 or 1) of this parameter after the inter marker
        distances line. Note that simuPOP populations always have phase information.

    label
        Output 1-based indexes of individuals if this parameter is true (default)

    popData
        Output 1-based index of subpopulation if this parameter is set to true (default).
    
    popFlag
        Output value of this parameter (0 or 1) after popData if this parameter specified.

    locData
        Name of an information field with location information of each individual. Default
        to None (no location data)

    phenotype
        Name of an information field with phenotype information of each individual. Default
        to None (no phenotype)
    
        
    Genotype information are always outputted. Alleles are coded the same way (0, 1, 2, etc)
    as they are stored in simuPOP.

    GENEPOP (http://genepop.curtin.edu.au/). The genepop format accepts the following
    parameters:

    title
        The tile line. If unspecified, a line similar to 'produced by simuPOP on XXX'
        will be outputted.

    adjust
        Adjust values of alleles by specified value (1 as default). This adjustment is
        necessary in many cases because GENEPOP treats allele 0 as missing values, and 
        simuPOP treats allele 0 as a valid allele. Exporting alleles 0 and 1 as 1 and 2
        will allow GENEPOP to analyze simuPOP-exported files correctly.

    Because 0 is reserved as missing data in this format, allele A is outputted as A+adjust.
    simuPOP will use subpopulation names (if available) and 1-based individual index
    to output individual label (e.g. SubPop2-3). If parameter subPops is used to output
    selected individuals, each subpop will be outputted as a separate subpopulation even 
    if there are multiple virtual subpopulations from the same subpopulation. simuPOP 
    currently only export diploid populations to this format.

    FSTAT (http://www2.unil.ch/popgen/softwares/fstat.htm). The fstat format accepts
    the following parameters:

    lociNames
        Names of loci that will be outputted. If unspecified, simuPOP will try to use
        names of loci that are specified by parameter *lociNames* of the ``Population``
        class, or names in the form of chrX-Y.

    adjust
        Adjust values of alleles by specified value (1 as default). This adjustment is
        necessary in many cases because FSTAT treats allele 0 as missing values, and 
        simuPOP treats allele 0 as a valid allele. Exporting alleles 0 and 1 as 1 and 2
        will allow FSTAT to analyze simuPOP-exported files correctly.
        
        
    CSV (comma separated values). This is a general format that output genotypes in
    comma (or tab etc) separated formats. The function form of this operator 
    ``export(format='csv')`` is similar to the now-deprecated ``saveCSV`` function,
    but its interface has been adjusted to match other formats supported by this
    operator. This format outputs a header (optiona), and one line for each individual
    with values of specified information fields, sex, affection status, and genotypes.
    All fields except for genotypes are optional. The output format is controlled by the
    following parameters:
    
    infoFileds
        Information fields to be outputted. Default to none.

    header
        Whether or not a header should be written. These headers will include
        information fields, sex (if ``sexFormatter`` is not ``None``), affection
        status (if ``affectionFormatter`` is not ``None``) and loci names. If
        genotype at a locus needs more than one column, ``_1``, ``_2`` etc will
        be appended to loci names. Alternatively, a complete header (a string)
        or a list of column names could be specified directly.

    infoFormatter
        A format string that is used to format all information fields. If
        unspecified, ``str(value)`` will be used for each information field.

    genoFormatter
        How to output genotype at specified loci. Acceptable values include
        ``None`` (output allele values), a dictionary with genotype as keys,
        (e.g. ``genoFormatter={(0,0):1, (0,1):2, (1,0):2, (1,1):3}``, or a function
        with genotype (as a tuple of integers) as inputs. The dictionary value
        or the return value of this function can be a single or a list of
        number or strings.

    sexFormatter
        How to output individual sex. Acceptable values include ``None`` (no
        output) or a dictionary with keys ``MALE`` and ``FEMALE``.

    affectionFormatter
        How to output individual affection status. Acceptable values include
        ``None`` (no output) or a dictionary with keys ``True`` and ``False``.

    delimiter
        Delimiter used to separate values, default to ','.

    This operator supports the usual applicability parameters such as begin,
    end, step, at, reps, and subPops. If subPops are specified, only
    individuals from specified (virtual) subPops are exported. Because this
    exporter will always overwite any existing file, leading '>' of parameter
    output is ignored, although the '!expr' format can still be used to
    generate context dependent output filename. Unless explicitly stated for
    a particular format, this operator exports individuals from the current
    generation if there are multiple ancestral generations in the population.
    
    The Exporter class will make use of a progress bar to show the progress. The
    interface of the progress bar is by default determined by the global GUI status
    but you can also set it to, for example, ``gui=False`` to forcefully use a 
    text-based progress bar.
    '''
    def __init__(self, format, output, begin=0, end=-1, step=1, at=[],
        reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[], gui=None, *args, **kwargs):
        self.output = output
        self.subPops = subPops
        self.infoFields = [infoFields] if type(infoFields) == type('') else infoFields
        self.gui = gui
        if format.lower() == 'structure':
            self.exporter = StructureExporter(*args, **kwargs)
        elif format.lower() == 'genepop':
            self.exporter = GenePopExporter(*args, **kwargs)
        elif format.lower() == 'fstat':
            self.exporter = FStatExporter(*args, **kwargs)
        elif format.lower() == 'csv':
            self.exporter = CSVExporter(*args, **kwargs)
        else:
            raise ValueError('Unrecognized fileformat: {}.'.format(format))
        PyOperator.__init__(self, func=self._export, begin=begin, end=end,
            step=step, at=at, reps=reps, subPops=ALL_AVAIL, infoFields=[])

    def _determineOutput(self, pop):
        # determine the output filename
        # the expression form
        if self.output.startswith('!'):
            # evaluate output with the population variable
            output = eval(self.output[1:], pop.vars(), pop.vars())
            return output.lstrip('>')
        else:
            return self.output.lstrip('>')
            
    def _determineSubPops(self, pop):
        # this is basically subPopList::expandFrom(pop)
        if self.subPops == ALL_AVAIL:
            return range(pop.numSubPop())
        elif type(self.subPops) == type(0):
            return [self.subPops]
        elif type(self.subPops) == type(''):
            try:
                return [pop.subPopNames().index(self.subPops)]
            except:
                raise ValueError('%s is not a valid subpop name' % self.subPops)
        # handle vsps such as (ALL_AVAIL, vsp)
        subPops = []
        for vsp in self.subPops:
            # is it a number?
            if type(vsp) == type(0):
                subPops.append(vsp)
            elif type(vsp) == type(''):
                subPops.append(pop.subPopNames().index(vsp))
            else:
                if type(vsp[0]) == type(''):
                    try:
                        vsp[0] = pop.subPopNames().index(vsp[0])
                    except:
                        raise ValueError('%s is not a valid subpop name' % vsp[0])
                if type(vsp[1]) == type(''):
                    try:
                        vsp[1] = pop.virtualSplitter().vspByName(vsp[1])
                    except:
                        raise ValueError('Population does not have any virtual subpopulation %s' % vsp[1])
                if vsp[0] == ALL_AVAIL:
                    for u in range(pop.numSubPop()):
                        if vsp[1] == ALL_AVAIL:
                            for v in range(pop.numVirtualSubPops()):
                                subPops.append([u, v])
                        else:
                            subPops.append([u, vsp[1]])
                else:
                    if vsp[1] == ALL_AVAIL:
                        for v in range(pop.numVirtualSubPops()):
                            subPops.append([vsp[0], v])
                    else:
                        subPops.append(vsp)
        return subPops
        
    def _export(self, pop):
        self.exporter.export(pop, self._determineOutput(pop),
            self._determineSubPops(pop), self.infoFields, gui=self.gui)
        return True

def export(pop, format, *args, **kwargs):
    '''Apply operator ``Exporter`` to population *pop* in format *format*.'''
    Exporter(format, *args, **kwargs).apply(pop)


def importPopulation(format, filename, *args, **kwargs):
    '''This function import and return a population from a file *filename* in
    specified *format*. Format-specific parameters can be used to define how the
    input should be interpreted and imported. This function supports the following
    file format.

    GENEPOP (http://genepop.curtin.edu.au/). For input file of this format, this
    function ignores the first title line, load the second line as loci names,
    and import genotypes of different POP sections as different subpopulations.
    This format accepts the following parameters:

    adjust
        Adjust alleles by specified value (default to 0 for no adjustment). This
        parameter is mostly used to convert alleles 1 and 2 in a GenePop file to
        alleles 0 and 1 (with adjust=-1) in simuPOP.

    FSTAT (http://www2.unil.ch/popgen/softwares/fstat.htm). This format accepts
    the following parameters:

    adjust
        Adjust alleles by specified value (default to 0 for no adjustment). This
        parameter is mostly used to convert alleles 1 and 2 in a GenePop file to
        alleles 0 and 1 (with adjust=-1) in simuPOP.
    '''
    if format.lower() == 'genepop':
        importer = GenePopImporter(*args, **kwargs)
    elif format.lower() == 'fstat':
        importer = FStatImporter(*args, **kwargs)
    else:
        raise ValueError('Importing genotypes in format %s is currently not supported' % format)
    return importer.importFrom(filename)

if __name__ == "__main__":
    pass

