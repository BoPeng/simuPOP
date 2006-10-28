#!/usr/bin/env python
#
# This scripts help users of simuPOP to submit his jobs to 
# a cluster system. To use this script
#
#
# 
# $LastChangedDate: 2006-01-02 20:59:22 -0600 (Mon, 02 Jan 2006) $
# $Rev: 110 $
# 

import os, sys, getopt, re

Usage = '''
1. create a list file (a python script), that defines

script: a script that has field $0, $1, $2, ..., $name,
    and other variables provided through command line or configuration
    file ($HOME/.simuCluster).

joblist: a multi-line string of the form
      name1: var1: var2:
      name2: var2: var2:
   The values will be used to substitute $0, $1, ...
   The job will be saved as $0.pbs. $0 is also named $name.


separator: separator to separate joblist, default to ':'

Note that if you have to use $ in the script, use $$ instead.

The list file can define script and jobs by itself, or simply code 
like:
    script = open('job.template').read()
    jobs = open('jobs').read()

2. use of variables:

  Any command line argument key=value (without leading --), and any 
  environment variables will be used to subsitute the script as well. 
  Exceptions are
    -t/--time, which can be either a number of short (4), medium (24), 
        long (48) or extended (96), and the corresponding number will be used. 
    -q/--queue, is the same as queue=xxx

3. configuration file:
  You can define machine specific variables in a configuration file 
  $HOME/.simuCluster. simuCluster.py process this python file and use the
  variables to subsititute variables in script. 


To use this script:
  python simuCluster.py [options]

where options can be:
  -h: help information
  -l: list file, default to simulation.lst
  -t|--time: walltime (number, or short, medium, long, extended)
  -q|--queue: queue to submit job to
  -f|--force: submit the job even it has $ in the (resulting) script 
  -c|--command: command used to submit jobs.

Example of one such script:

script = """
#!/bin/bash
#PBS -N $name
#PBS -l nodes=1:ppn=2:myrinet,walltime=$time:00:00
#PBS -M emailaddress
#PBS -m abe
#PBS -o /path/to/output
#PBS -e /path/to/error
#
mkdir -p /shared.scratch/jobs
cd /shared.scratch/jobs
/bin/cp -f /users/bpeng/somescript somescript
# command
python script.py --opt1=$opt1 --opt2=$opt2
"""

joblist = """
simu1: 1: 2
simu2: 2: 1
"""

Example $HOME/.simuCluster

command = 'bsub -n 1 <'
queue = 'normal'

'''

alljobs = []


def getJobs():
    ''' process jobs or joblist and create alljobs '''
    global alljobs
    alljobs = []
    # now process jobs
    if globals().has_key('separator'):
        sep = separator
    else:
        sep = ':'
    print "Getting job list from joblist"    
    for job in joblist.split('\n'):
        # get keys 0, 1, 2, ...
        keys = {}
        # remove empty lines
        if sep not in job:
            continue
        for i,key in enumerate(job.split(sep)):
            keys[str(i)] = key.strip()
        keys['name'] = keys['0']
        alljobs.append(keys)
    return


def getScript(name, options):
    ''' return the script for simulation 'name' '''
    if not globals().has_key('script'):
        print 'Vairable script is not defined'
        sys.exit(0)
    #
    for job in alljobs:
        # a dictionary
        if job['name'] == name:
            s = script
            for k,v in (job.items() + options.items() + os.environ.items()):
                # first protect double $ sign.
                s = re.sub(r'\$\$', r'%%%', s)
                # subsitute in script
                if True in [x in str(v) for x in [' ', '*', ',', '[', ']']]:
                    if '"' in str(v):
                        quote = '"'
                    else:
                        quote = "'"
                    s = re.sub(r'\$%s([^a-zA-Z])' % k, r"%s%s%s\1" % (quote, str(v), quote), s)
                    s = re.sub(r'\$\{%s\}' % k, "%s%s%s" % (quote, str(v), quote), s)
                else:
                    s = re.sub(r'\$%s([^a-zA-Z])' % k, r"%s\1" % str(v), s)
                    s = re.sub(r'\$\{%s\}' % k, str(v), s)
                # $$ is replaced with single $
                s = re.sub(r'%%%', r'$', s)
            return s


def allJobs():
    ''' list all jobs '''
    names = []
    for job in alljobs:
        names.append(job['name'])
    return names


def readConfigFile():
    ''' read variables from configuration file '''
    config = os.path.join(os.environ['HOME'], '.simuCluster')
    tmp = {}
    res = {}
    if os.path.isfile(config):
        execfile(config, tmp, res)
    return res


if __name__ == '__main__':
    # get pptions
    # options can can be used to subst fields in script
    options = {'time':96}
    simuList = 'simulation.lst'
    proc_jobs = []
    run = False
    force = False
    # default command to run the job, can be, for example sh
    options['command'] = 'qsub'
    # read configuration file
    options.update(readConfigFile())
    #
    optlist, args = getopt.gnu_getopt(sys.argv[1:], 't:l:s:ahrq:fc:', 
      ['list=', 'show=', 'time=', 'all', 'run', 'help', 'force', 'command'])
    for opt in optlist:
        if opt[0] == '-t' or opt[0] == '--time':
            try:
                time = {'short':4, 'medium':24, 'long':48, 'extended':96, 'verylong':96}[opt[1]]
            except: # if not one of them
                try:
                    time = int(opt[1])
                except:
                    print opt[1], " is not a valid time name (or number)"
                    sys.exit(1)
            options['time'] = str(time)
        elif opt[0] in ['-q', '--queue']:
            options['queue'] = opt[1]
        elif opt[0] == '-l' or opt[0] == '--list':
            simuList = opt[1]
        elif opt[0] in ['-s', '--show']:
            proc_jobs = [opt[1]]
        elif opt[0] in ['-a', '--all']:
            proc_jobs = 'all'
        elif opt[0] in ['-r', '--run']:
            run = True
        elif opt[0] in ['-f', '--force']:
            force = True
        elif opt[0] in ['-h', '--help']:
            print Usage
            sys.exit(0)
        elif opt[0] in ['-c', '--command']:
            options['command'] = command
    #
    # this is a special case
    #
    if not os.path.isfile(simuList):
        print 'Simulation list file does not exist'
        sys.exit(1)
    #
    print 'Loading simulation list ', simuList
    execfile(simuList)
    #
    getJobs()
    # 
    all_jobs = allJobs()
    if proc_jobs == 'all':
        proc_jobs = all_jobs
    #
    # submit some jobs
    scan = re.compile('(\w+[^\d]+)(\d+)-(\d+)')
    for j in args:
        if '=' in j:
            (key,val) = j.split('=')
            options[key] = val
        else:
            try:
                (name, start, end) = scan.match(j).groups()
                for i in range(int(start), int(end)+1):
                    for n in ['%s%d' % (name, i), '%s%02d' % (name, i), '%s%03d' % (name, i), '%s%04d' % (name, i)]:
                        if n not in proc_jobs and n in all_jobs:
                            proc_jobs.append(n)
            except:
                if j not in proc_jobs and j in all_jobs:
                    proc_jobs.append(j)
    # create .pbs scripts, run them if -r
    print "Creating .pbs files for job ", ' '.join(proc_jobs)
    for job in proc_jobs:
        if job in all_jobs:
            pbs_script = getScript(job, options)
            pbs = open(job + '.pbs', 'w')
            print >> pbs, pbs_script
            pbs.close()
            if '$' in pbs_script and not force:
                print pbs_script
                print
                print 'Warning: symbol $ exists in the script, indicating unsubstituted variables'
                print
                for line in pbs_script.split():
                    if '$' in line:
                        print '>> ', line
                print
                print 'Please check your pbs_script, if there is no problem, please use option -f (--force)'
                print 'to submit the job'
                sys.exit(1)
            if run:
                command = options['command'].replace('$name', job)
                print "Submitting job using command '%s %s.pbs'" % (command, job)
                os.system('%s %s.pbs' % (command, job))
        else:
            print "Job %s does not exist" % job

