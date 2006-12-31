#!/usr/bin/env python

# make src and binary distribution on all platforms 
# (currently mac, linux and solaris) as daily snapshot
#
# This is part of simuPOP
# 
# Bo Peng (bpeng@rice.edu)
#

import os, sys, time

last_version_file = '/var/www/html/simuPOP/download/latestversion'
last_revision_file = '/var/www/html/simuPOP/download/latestrevision'
release_file = '/home/bpeng/simuPOP/simuPOP.release'
download_directory = '/var/www/html/simuPOP/download'
doc_directory = '/var/www/html/simuPOP_doc/doc'
user_tmp_directory = '/home/bpeng/tmp'


def cmdOutput(cmd):
    ''' utility function: run a command and get its output as a string
        cmd: command to run
    '''
    fout = os.popen(cmd)
    output = fout.read()
    fout.close()
    return output.strip()


def removeTempFiles():
    ''' remove unnecessary files '''
    # walk down the file hierarchy
    for root,path,files in os.walk('.'):
        # check file type
        for file in files:
            if True in [('_%s.cpp' % x in file) for x in ['std', 'op', 'la', 'laop', 'ba', 'baop']] or \
                True in [('_%s_wrap.cpp' % x in file) for x in ['std', 'op', 'la', 'laop', 'ba', 'baop']]:
                # force the use of / since miktex uses / even under windows
                print "Remove file %s..." % os.path.join(root, file)
                os.remove(os.path.join(root, file))
    

def commitModification():
    ''' if there are changes, commit it '''
    if cmdOutput('svn diff') != '':
        cmd = 'svn ci -m "automatic checkin on %s"' % time.asctime()
        print cmd
        os.system(cmd)
    os.system('svn update')


def checkRevision():
    if os.path.isfile(last_revision_file):
        last_rev = open(last_revision_file).read()
    else:
        last_rev = None
    # check revision number
    cur_rev = cmdOutput('svnversion .')
    #
    # write to the file
    file = open(last_revision_file, 'w')
    file.write(cur_rev)
    file.close()
    return (last_rev, cur_rev)


def writeVersionFile(release):
    # write to version file
    file = open(last_version_file, 'w')
    file.write(release)
    file.close()
    

def makeReleaseTag(release):
    cmd = 'svn copy https://svn.sourceforge.net/svnroot/simupop/trunk ' + \
        'https://svn.sourceforge.net/svnroot/simupop/tag/v%s' % release + \
        ' -m "Version %s released at %s"' % (release, time.asctime())
    print cmd
    os.system(cmd)


def writeReleaseFile(release, revision):
    ' write release file only in release mode '
    res = {}
    execfile(release_file, res, res)
    file = open(release_file, 'w')
    file.write('''SIMUPOP_VER = "%s"
SIMUPOP_REV = "%s"''' % (release, revision))
    file.close()
    return (res['SIMUPOP_VER'], res['SIMUPOP_REV'])


def buildDocument():
    ' build doxygen document '
    os.environ['SIMUPOP_DOC_DIR'] = doc_directory
    os.system('doxygen Doxyfile')
    os.system('python tools/doxy2swig.py %s/xml/index.xml tmp.i' % doc_directory)
    os.system('perl tools/processDocString.pl tmp.i > src/simuPOP/download.i')
    os.remove('tmp.i')
    os.chdir('doc')
    print 'Entering ', os.getcwd()
    os.system('postGuide.sh')
    os.chdir('..')
    print 'Entering ', os.getcwd()


def buildTarget(dest, script, message, wait=os.P_WAIT):
    dest_file = os.path.join(download_directory, dest % release)
    print 'Building %s...' % dest_file
    if os.path.isfile(dest_file):
        os.remove(dest_file)
    #
    try:
        print "Running command sh %s" % script
        os.spawnlp(wait, 'sh', 'sh', script)
    except:
        print "Can not start process %s " % script
    #
    #if not os.path.isfile(dest_file):
    #    print 'Can not find %s' % dest_file, message
    #    sys.exit(1)


if __name__ == '__main__':
    force_build = False
    release = 'snapshot'
    actions = []
    actions_exclude = []
    all_actions = ['svn', 'doc', 'src', 'linux', 'mac', 'mdk', 'win', 'sol', 'rhel4']

    ## Parse the command line
    for op in sys.argv[1:]:   # default shell/for list is $*, the options
        if op in [ '-help', '--help', '-h' ]:
            print '''Usage: snapshot [options]
Options:
    --help                  show this help lines
    --force-build           force building even if there is no change
    --release=version       release simupop with version number
    --actions=%s
                            actions to run, default to all
    --exclude=[]            exclude some actions.
''' % (','.join(all_actions))
            sys.exit(0)
        elif op == '--force-build':
            force_build = True
        elif '--release' in op:
            force_build = True
            release = op[10:]
        elif '--actions' in op:
            actions = op[10:].split(',')
        elif '--exclude' in op:
            actions_exclude = op[10:].split(',')
        else:
            print "Unknown option", op
            sys.exit(1)
    # decide what to do ....
    if actions == []:
        actions = all_actions
    else:
        for act in actions:
            if act not in all_actions:
                print "Wrong action %s " % act
                sys.exit(1)
    for act in actions_exclude:
        if act in actions:
            actions.remove(act)
    #
    # go to top directory
    os.chdir('..')
    removeTempFiles()
    #
    if 'svn' in actions:
        # check in unsubmitted changes
        commitModification()
    # get the old and new revision number
    (old_rev, revision) = checkRevision()
    #
    # avoid uncessary rebuild
    if old_rev == revision and not force_build:
        sys.exit(0)
    #
    if 'svn' in actions and release != 'snapshot':
        # if --release, make a tag
        makeReleaseTag(release)
        writeReleaseFile(release, revision)
        commitModification()
    #
    # start building
    print "Building simuPOP revision %s..." % revision
    writeVersionFile(release)
    os.environ['SIMUPOP_VER'] = release
    os.environ['SIMUPOP_REV'] = revision
    if 'doc' in actions:
        buildDocument()
    os.chdir('tools')
    print 'Entering ', os.getcwd()
    if 'src' in actions:
        # make release file
        # In the case of snapshot, we need to restore it.
        (old_ver, old_rev) = writeReleaseFile(release, revision)
        buildTarget('simuPOP-%s-src.tar.gz', 'make_src.sh', 'Building src fails')
        # restore it
        writeReleaseFile(old_ver, old_rev)
    #if 'sol' in actions:
    #    buildTarget('simuPOP-%s-sol-py23.tar.gz', 'make_sol.sh', 'Building solaris binary fails', os.P_NOWAIT)
    if 'mac' in actions:
        buildTarget('simuPOP-%s.darwin-7.7.0-PowerMacintosh.tar.gz', 'make_mac.sh', 'Building mac binary fails', os.P_NOWAIT)
    if 'linux' in actions:
        buildTarget('simuPOP-%s.linux-i686-py23.tar.gz', 'make_linux.sh', 'Building linux binary fails', os.P_NOWAIT)
    if 'rhel4' in actions:
        buildTarget('simuPOP-%s.linux-i686-py23.tar.gz', 'make_rhel4.sh', 'Building linux binary fails', os.P_WAIT)
    if 'win' in actions:
        buildTarget('simuPOP-%s.win32-py2.3.exe', 'make_win.sh', 'Building windows binary fails', os.P_WAIT)
    if 'mdk' in actions:
        buildTarget('simuPOP-%s-mdk-py23.tar.gz', 'make_mdk.sh', 'Building mandriva binary fails', os.P_WAIT)

    
