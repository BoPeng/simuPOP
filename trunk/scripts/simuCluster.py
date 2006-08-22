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
1. create a list file (a python script), that define

script: a script that has field $0, $1, $2, ...
   or $name1, $name2, ...
   ${0} or ${name} can be used to separate variables.

jobs: a list of
    a list which items will be used to substitute $0, $1, ...
    The job will be saved as $0.pbs
  or
    a dictionary with keys to subsitute $name1, $name2, one
    of the keys have to be 'name'. The job will be saved as
    $name.pbs
 
  Any command line argument --key=value will be used to subsitute
  the script as well. An except is -t/--time, which can be either
  a number of short (4), medium (24), long (48) or extended (96),
  and the corresponding number will be used.


joblist: a string of the form
  name1: var1: var2:
  name2: var2: var2:

separator: separator to separate joblist, default to ':'


Either jobs or joblist should exist.



To use this script:
  python simuCluster.py [options]

where options can be:
  -h: help information
  -l: list file, default to simulation.lst
  -t|--time: walltime (number, or short, medium, long, extended)
  -s|--show: show all jobs
  -d|--dryrun: does not actually submit job
  
  
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
    jobs = [
        {'name': 'simu1', 'opt1':1, 'opt2'2},
        {'name': 'simu2', 'opt1':1, opt2':2}
    ]

'''

alljobs = []


def getJobs():
    ''' process jobs or joblist and create alljobs '''
    global alljobs
    alljobs = []
    if globals().has_key('jobs'):
        #
        for job in jobs:
            # a dictionary
            if type(job) == type({}):
                if not job.has_key('name'):
                    print 'Dictionary job entry does not have key name'
                    sys.exit(1)
                else:
                    alljobs.append(job)
            elif type(job) in [type([]), type(())]:
                keys = {}
                for i,key in enumerate(job):
                    keys[str(i)] = key
                keys['name'] = keys['0']
                alljobs.append(keys)
        return
    elif not globals().has_key('joblist'):
        print 'Either jobs or jobliss has to be defined'
        sys.exit(1)
    # now process jobs
    if globals().has_key('separator'):
        sep = separator
    else:
        sep = ':'
    print "Getting job list from joblist"    
    for job in joblist.split('\n'):
        # get keys 0, 1, 2, ...
        keys = {}
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
            for k,v in (job.items() + options.items()):
                # subsitute in script
                if True in [x in str(v) for x in [' ', '*', ',', '[', ']']]:
                    if '"' in str(v):
                        quote = '"'
                    else:
                        quote = "'"
                    s = re.sub(r'\$%s(\W)' % k, r"%s%s%s\1" % (quote, str(v), quote), s)
                    s = re.sub(r'\$\{%s\}' % k, "%s%s%s" % (quote, str(v), quote), s)
                else:
                    s = re.sub(r'\$%s(\W)' % k, r"%s\1" % str(v), s)
                    s = re.sub(r'\$\{%s\}' % k, str(v), s)
            return s


def allJobs():
    ''' list all jobs '''
    names = []
    for job in alljobs:
        names.append(job['name'])
    return names


if __name__ == '__main__':
    # get pptions
    # options can can be used to subst fields in script
    options = {'time':96}
    simuList = 'simulation.lst'
    proc_jobs = []
    run = False
    #
    optlist, args = getopt.getopt(sys.argv[1:], 't:l:s:ahr', 
      ['list=', 'show=', 'all', 'run', 'help'])
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
        elif opt[0] == '-l' or opt[0] == '--list':
            simuList = opt[1]
        elif opt[0] in ['-s', '--show']:
            proc_jobs = [opt[1]]
        elif opt[0] in ['-a', '--all']:
            proc_jobs = 'all'
        elif opt[0] in ['-r', '--run']:
            run = True
        elif opt[0] in ['-h', '--help']:
            print Usage
            sys.exit(0)
    #
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
        if j[:2] == '--' and '=' in j:
            (key,val) = j[2:].split('=')
            options[key] = val
        else:
            try:
                (name, start, end) = scan.match(j).groups()
                for i in range(int(start), int(end)+1):
                    if name+str(i) not in proc_jobs:
                        proc_jobs.append( name+str(i))
            except:
                if j not in proc_jobs:
                    proc_jobs.append(j)
    #
    if not run:
        if len(proc_jobs) > 1:
            print ' '.join(filter(lambda x: x in all_jobs, proc_jobs))
        elif len(proc_jobs) == 1:
            if proc_jobs[0] in all_jobs:
                print "\nList script", proc_jobs[0]
                print getScript(proc_jobs[0], options)
            else:
                print "\nJob %s does not exist. " % proc_jobs[0]
    else:
        # submit the jobs
        for job in proc_jobs:
            if job in all_jobs:
                pbs = open(job + '.pbs', 'w')
                print >> pbs, getScript(job, options)
                pbs.close()
                print "Submitting job via 'qsub %s.pbs'" % job 
                os.system('qsub %s.pbs' % job)
            else:
                print "Job %s does not exist" % job

