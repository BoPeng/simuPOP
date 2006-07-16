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
release_file = 'simuPOP.release'
download_directory = '/var/www/html/simuPOP/download'
doc_directory = '/var/www/html/simuPOP_doc/doc'


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
        os.system('svn ci -m "automatic checkin on %s"' % time.asctime())
    os.system('svn update')


def checkRevision():
    if os.path.isfile(last_revision_file):
        rev = open(last_revision_file).read()
    else:
        rev = None
    # check revision number
    cur_rev = cmdOutput('svnversion')
    #
    # write to the file
    file = open(last_revision_file, 'w')
    file.write(cur_rev)
    file.close()
    return (rev, cur_rev)


def writeVersionFile(release):
    # write to version file
    file = open(last_version_file, 'w')
    file.write(release)
    file.close()
    

def makeReleaseTag(release):
    os.system('svn copy http://bp6.stat.rice.edu:8080/svn/simuPOP/trunk ' + \
        'http://bp6.stat.rice.edu:8080/svn/simuPOP/tag/v%s' % release + \
        '-m "Version $1 released at %s"' % (release, time.asctime()))


def writeReleaseFile(release, revision):
    ' write release file only in release mode '
    file = open(release_file, 'w')
    file.write('''SIMUPOP_VER = "%s"
SIMUPOP_REV = "%s"''' % (release, revision))
    file.close()


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


def buildTarget(dest, script, message):
    dest_file = os.path.join(download_directory, dest % release)
    print 'Building %s...' % dest_file
    if os.path.isfile(dest_file):
        os.remove(dest_file)
    #
    os.system('sh ' + script)
    #
    if not os.path.isfile(dest_file):
        print 'Can not find %s' % dest_file, message
        sys.exit(1)


if __name__ == '__main__':
    force_build = False
    release = 'snapshot'
    actions = []
    all_actions = ['svn', 'doc', 'src', 'linux', 'mac', 'win', 'sol']

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
''' % (','.join(all_actions))
            sys.exit(0)
        elif op == '--force-build':
            force_build = True
        elif op[:9] == '--release':
            force_build = True
            release = op[10:]
        elif op[:9] == '--actions':
            actions = op[10:].split(',')
        else:
            print "Unknown option", op
            sys.exit(1)
    if actions == []:
        actions = all_actions
    else:
        for act in actions:
            if act not in all_actions:
                print "Wrong action %s " % act
                sys.exit(1)
    # go to top directory
    os.chdir('..')
    print 'Entering ', os.getcwd()
    #
    removeTempFiles()
    #
    if 'svn' in actions:
        # check in unsubmitted changes
        commitModification()
        #
    (old_rev, revision) = checkRevision()
    #
    if old_rev == revision and not force_build:
        sys.exit(0)
    if 'svn' in actions:
        # if --release, make a tag
        if release != 'snapshot':
            makeReleaseTag(release)
            writeReleaseFile(release, revision)
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
        buildTarget('simuPOP-%s-src.tar.gz', 'make_src.sh', 'Building src fails')
    if 'linux' in actions:
        buildTarget('simuPOP-%s.linux-i686-py23.tar.gz', 'make_linux.sh', 'Building linux binary fails')
    if 'win' in actions:
        buildTarget('simuPOP-%s.win32-py2.3.exe', 'make_win.sh', 'Building windows binary fails')
    if 'sol' in actions:
        buildTarget('simuPOP-%s-sol-py23.tar.gz', 'make_sol.sh', 'Building solaris binary fails')
    if 'mdk' in actions:
        buildTarget('simuPOP-%s-mdk-py23.tar.gz', 'make_mdk.sh', 'Building mandriva binary fails')
    if 'mac' in actions:
        buildTarget('simuPOP-%s.darwin-7.7.0-PowerMacintosh.tar.gz', 'make_mac.sh', 'Building mac binary fails')

    
