#!/usr/bin/env python

# make src and binary distribution on all platforms 
# (currently mac, linux and solaris) as daily snapshot
#
# This is part of simuPOP
# 
# Bo Peng (bpeng@mdanderson.org)
#

import os, sys, time, shutil

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
    modules = ['std', 'op', 'la', 'laop', 'ba', 'baop', \
        'mpi', 'opmpi', 'lampi', 'laopmpi', 'bampi', 'baopmpi']
    for root,path,files in os.walk('.'):
        # check file type
        for file in files:
            if True in [('_%s.cpp' % x in file) for x in modules]:
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


def writeReleaseFile(release, revision):
    ' write release file only in release mode '
    res = {}
    execfile(release_file, res, res)
    file = open(release_file, 'w')
    file.write('''SIMUPOP_VER = "%s"
SIMUPOP_REV = "%s"''' % (release, revision))
    file.close()
    return (res['SIMUPOP_VER'], res['SIMUPOP_REV'])


def setVersionRevision(release):
    ''' if release = snapshot, do not change simuPOP.release
        by returning current ver and rev numbers
        otherwise, update simuPOP.release

        last_version_file and last_revision_file will always
        be updated.
    '''
    if release == 'snapshot':
        rev = cmdOutput('svnversion .')
        if ':' in rev:
            rev = rev.split(':')[0]
    # write last_revision_file and last_version_file
    file = open(last_revision_file, 'w')
    file.write(rev)
    file.close()
    file = open(last_version_file, 'w')
    file.write(release)
    file.close()
    # replace simuPOP.release
    (old_ver, old_rev) = writeReleaseFile(release, rev)
    return (release, rev, old_ver, old_rev)
    

def makeReleaseTag(release):
    cmd = 'svn copy https://svn.sourceforge.net/svnroot/simupop/trunk ' + \
        'https://svn.sourceforge.net/svnroot/simupop/tag/v%s' % release + \
        ' -m "Version %s released at %s"' % (release, time.asctime())
    print cmd
    os.system(cmd)



def build_doc(ver, rev):
    ' build doxygen document '
    os.environ['SIMUPOP_DOC_DIR'] = doc_directory
    os.environ['SIMUPOP_VER'] = ver
    os.environ['SIMUPOP_REV'] = rev
    os.system('doxygen Doxyfile')
    os.system('python tools/doxy2swig.py %s/xml/index.xml tmp.i' % doc_directory)
    os.system('perl tools/processDocString.pl tmp.i > src/simuPOP/download.i')
    os.remove('tmp.i')
    os.chdir('doc')


def build_src():
    ver = open(last_version_file).read()

    if ver == 'snapshot':
        rev = cmdOutput('svnversion .')
        if ':' in rev:
            rev = rev.split(':')[0]
    else:
        rev = open(last_revision_file).read()
    #
    # replace simuPOP.release file
    (old_ver, old_rev) = writeReleaseFile(ver, rev)
    # build source
    os.system('python setup.py sdist --formats=gztar,zip')
    # write old release file back
    writeReleaseFile(old_ver, old_rev)
    # coppy files
    shutil.copy('dist/simuPOP-%s.tar.gz' % ver, '%s/simuPOP-%s-src.tar.gz' % (download_directory, ver))
    shutil.copy('dist/simuPOP-%s.zip' % ver, '%s/simuPOP-%s-src.zip' % (download_directory, ver))


def build_x86_64():
    if not os.path.isfile('%s/simuPOP-%s.zip' % (download_directory, ver)):
        print 'Source package does not exist. Please run build_src.py first'
        sys.exit(1)
    # 
    # build
    print 'Copying source package to user temp directory...'
    os.system('/bin/rm -rf simuPOP-%s' % ver)
    shutil.copy('%s/simuPOP-%s-src.tar.gz' % (download_directory, ver), 'simuPOP-%s-src.tar.gz' % ver)
    print 'Unpacking ...'
    os.system('tar zxf simuPOP-%s-src.tar.gz' % ver)
    os.chdir('simuPOP-%s' % ver)
    print 'Building ...'
    os.system('python setup.py bdist --formats=gztar,rpm')
    # coppy files
    shutil.copy('dist/simuPOP-%s-1.x86_64.rpm' % ver, '%s/simuPOP-%s-1.x86_64.rpm' % (download_directory, ver))
    shutil.copy('dist/simuPOP-%s-1.src.rpm' % ver, '%s/simuPOP-%s-1.src.rpm' % (download_directory, ver))


def build_mac():
    if not os.path.isfile('%s/simuPOP-%s.zip' % (download_directory, ver)):
        print 'Source package does not exist. Please run build_src.py first'
        sys.exit(1)
    #   
    print 'Copying source package to remote machine ...'
    os.system("ssh %s '/bin/rm -rf temp &&  mkdir temp && /bin/rm -rf simuPOP'" % mac_name)
    os.system('scp %s/simuPOP-%s-src.tar.gz %s:temp' % (download_directory, ver, mac_name))
    #
    print 'Building ...'
    unpack = 'tar zxf simuPOP-%s-src.tar.gz' % ver
    build = 'python setup.py bdist_dumb'
    os.system('ssh -X %s "cd temp && %s && cd simuPOP-%s && %s"' % (mac_name, unpack, ver, build))
    os.system('scp %s:temp/simuPOP-%s/dist/* %s' % (mac_name, ver, download_directory))


if __name__ == '__main__':
    config_file = 'build.cfg' 
    release = 'snapshot'
    actions = []
    all_actions = ['all', 'svn', 'src', 'doc', 'x86_64', 'rhel4', 'mac', 'win', 'mdk', 'sol', 'fedora']

    ## Parse the command line
    for op in sys.argv[1:]:   # default shell/for list is $*, the options
        if op in [ '-help', '--help', '-h' ]:
            print '''Usage: build.py [options] action1 action2 
Options:
    --help                  show this help lines
    --release=version       release simupop with version number, default to snapshot
    --config=file           configuration file
        This configuration file should define:
        last_version_file: a file that has the last version, e.g.  1.7.5
        last_revision_file: a file that has the last revision. e.g. 567
        release_file: the release file or simuPOP distribution
        download_directory: where to put resulting packages
        doc_directory: doc directory of simuPOP source
        user_tmp_directory: temp directory to build simuPOP locally.
        mac_name: name of the mac machine
actions:    
    svn:    submit changes
    src:    build source package .tar.gz and .zip
    doc:    process simuPOP document
    x96_64: build x86_64 rpm package, including -src.rpm
    rhel4:  build rpm and .tar.gz packages for rhel4/i386
    mac:    package for macintosh
    win:    packages for windows
    mdk:    packages for madravia
    sol:    packages for solaris
    fedora: packages for fedora
    all:    src + x86_64 + rhel4 + mac + win
''' % ' '.join(all_actions)
            sys.exit(0)
        elif '--release' in op:
            release = op[10:]
        elif '--config' in op:
            config_file = op[9:]
        elif op == 'all':
            actions.extend(['src', 'x86_64', 'rhel4', 'mac', 'win'])
        elif op in all_actions and op not in actions:
            actions.append(op)
        else:
            print "Unknown option", op
            sys.exit(1)
    # 
    if os.path.isfile(config_file):
        execfile(config_file)
    else:
        print 'Configuration file not found'
        sys.exit(1)
    # get revision number and update last_revision_file
    (ver, rev, old_ver, old_rev) = setVersionRevision(release)
    #
    os.chdir('..')
    removeTempFiles()
    if 'svn' in actions:
        commitModification()
        if release != 'snapshot':
            makeReleaseTag(release)
    if 'src' in actions:
        build_src()
    if 'doc' in actions:
        build_doc(ver, rev)
    if 'x86_64' in actions:
        build_x86_64()
    if 'mac' in actions:
        build_mac()
    # 
    # restore simuPOP.release
    if release == 'snapshot':
        writeReleaseFile(old_ver, old_rev)

    
