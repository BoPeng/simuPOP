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

def run(cmd):
    print cmd
    if not dryrun:
        os.system(cmd)

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
        run(cmd)
    run('svn update')


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
    if release != 'snapshot':
        rev = cmdOutput('svnversion .')
        if ':' in rev:
            rev = rev.split(':')[0]
    else:
        rev = '9999'
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
    run(cmd)



def build_doc(ver, rev):
    ' build doxygen document '
    os.environ['SIMUPOP_DOC_DIR'] = doc_directory
    os.environ['SIMUPOP_VER'] = ver
    os.environ['SIMUPOP_REV'] = rev
    run('doxygen Doxyfile')
    run('python tools/doxy2swig.py %s/xml/index.xml tmp.i' % doc_directory)
    run('perl tools/processDocString.pl tmp.i > src/simuPOP_doc.i')
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
    run('python setup.py --extra_link_args "-static" sdist --formats=gztar,zip')
    # write old release file back
    writeReleaseFile(old_ver, old_rev)
    # coppy files
    shutil.copy('dist/simuPOP-%s.tar.gz' % ver, '%s/simuPOP-%s-src.tar.gz' % (download_directory, ver))
    print 'Saving ', '%s/simuPOP-%s-src.tar.gz' % (download_directory, ver)
    shutil.copy('dist/simuPOP-%s.zip' % ver, '%s/simuPOP-%s-src.zip' % (download_directory, ver))
    print 'Saving ', '%s/simuPOP-%s-src.zip' % (download_directory, ver)


def build_x86_64(ver):
    source = '%s/simuPOP-%s-src.tar.gz' % (download_directory, ver)
    if not os.path.isfile(source):
        print 'Source package %s does not exist. Please run "build.py src" first' % source
        sys.exit(1)
    # 
    # build
    d = os.getcwd()
    os.chdir(user_tmp_directory)
    print 'Copying source package to user temp directory...'
    run('/bin/rm -rf simuPOP-%s' % ver)
    shutil.copy(source, 'simuPOP-%s-src.tar.gz' % ver)
    print 'Unpacking ...'
    run('tar zxf simuPOP-%s-src.tar.gz' % ver)
    os.chdir('simuPOP-%s' % ver)
    print 'Building ...'
    run('python setup.py --extra_link_args "-static" bdist --formats=gztar,rpm')
    # coppy files
    shutil.copy('dist/simuPOP-%s-1.x86_64.rpm' % ver, '%s/simuPOP-%s-1.rhel4.x86_64.rpm' % (download_directory, ver))
    shutil.copy('dist/simuPOP-%s.linux-x86_64.tar.gz' % ver, '%s/simuPOP-%s.rhel4.x86_64.tar.gz' % (download_directory, ver))
    shutil.copy('dist/simuPOP-%s-1.src.rpm' % ver, '%s/simuPOP-%s-1.src.rpm' % (download_directory, ver))
    os.chdir(d)


def build_vm(ver, name, pyver, vm, vm_port, vm_name):
    # set display, since this script will be launched from crontab.
    run('vmrun start %s' % vm)
    # wait for the vm to start.
    run('sleep 30')
    run('ssh -p %d %s "/bin/rm -rf simuPOP-%s.zip simuPOP-%s.tar.gz simuPOP-%s"' % \
        (vm_port, vm_name, ver, ver, ver))
    run('scp -P %s %s/simuPOP-%s-src.tar.gz %s:' % \
        (vm_port, download_directory, ver, vm_name))
    #
    run('ssh -X -p %d %s "tar zxf simuPOP-%s-src.tar.gz && cd simuPOP-%s && python setup.py --extra_link_args "-static" bdist --formats=gztar,rpm"' % \
        (vm_port, vm_name, ver, ver))
    run('scp -P %d %s:simuPOP-%s/dist/simuPOP-%s.linux-i686.tar.gz %s/simuPOP-%s-%s-py2%2d.tar.gz' % \
        (vm_port, vm_name, ver, ver, download_directory, ver, name, pyver))
    run("scp -P %d '%s:simuPOP-%s/dist/simuPOP-%s*' %s/simuPOP-%s-%s-py%2d.tar.gz" % \
        (vm_port, vm_name, ver, ver, download_directory, ver, name, pyver))
    run('vmrun suspend %s' % vm)


def build_mdk(ver):
    build_vm(ver, 'mdk', 23, mdk_vm, mdk_port, mdk_vm_name)

def build_fedora5(ver):
    build_vm(ver, 'fedora5', 23, fedora5_vm, fedora5_port, fedora5_vm_name)

def build_rhel4(ver):
    build_vm(ver, 'rhel4', 23, rhel4_vm, rhel4_port, rhel4_vm_name)

def build_suse(ver):
    build_vm(ver, 'suse', 23, suse_vm, suse_port, suse_vm_name)

def build_mac():
    source = '%s/simuPOP-%s-src.tar.gz' % (download_directory, ver)
    if not os.path.isfile(source):
        print 'Source package %s does not exist. Please run "build.py src" first' % source
        sys.exit(1)
    #   
    print 'Copying source package to remote machine ...'
    run("ssh %s '/bin/rm -rf temp &&  mkdir temp && /bin/rm -rf simuPOP'" % mac_name)
    run('scp %s/simuPOP-%s-src.tar.gz %s:temp' % (download_directory, ver, mac_name))
    #
    print 'Building ...'
    unpack = 'tar zxf simuPOP-%s-src.tar.gz' % ver
    build = 'python setup.py --extra_link_args "-static" bdist_dumb'
    run('ssh -X %s "cd temp && %s && cd simuPOP-%s && %s"' % (mac_name, unpack, ver, build))
    run('scp %s:temp/simuPOP-%s/dist/* %s' % (mac_name, ver, download_directory))

Usage = '''Usage: build.py [options] action1 action2 
Options:
    --help                  show this help lines
    --release=version       release simupop with version number, default to snapshot
    --dryrun                show command, do not run
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
    x86_64: build x86_64 rpm package, including -src.rpm
    rhel4:  build rpm and .tar.gz packages for rhel4/i386
    mac:    package for macintosh
    win:    packages for windows
    mdk:    packages for madravia
    sol:    packages for solaris
    fedora: packages for fedora
    all:    src + x86_64 + rhel4 + mac + win
'''

if __name__ == '__main__':
    config_file = 'build.cfg' 
    release = 'snapshot'
    actions = []
    dryrun = False
    all_actions = ['all', 'svn', 'src', 'doc', 'x86_64', 'rhel4', 'mac', 'win', 'mdk', 'suse', 'sol', 'fedora5']

    if len(sys.argv) == 1:
        print Usage
        sys.exit(0)
    ## Parse the command line
    for op in sys.argv[1:]:   # default shell/for list is $*, the options
        if op in [ '-help', '--help', '-h' ]:
            print Usage
            sys.exit(0)
        elif '--release' in op:
            release = op[10:]
        elif '--config' in op:
            config_file = op[9:]
        elif op == 'all':
            actions.extend(['src', 'doc', 'svn', 'x86_64', 'rhel4', 'mac', 'win', 'fedora5', 'suse'])
        elif op in all_actions and op not in actions:
            actions.append(op)
        elif op == '--dryrun':
            dryrun = True
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
    print 'New version: %s, revision: %s ' % (ver, rev)
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
        build_x86_64(ver)
    if 'mdk' in actions:
        build_mdk(ver)
    if 'fedora5' in actions:
        build_fedora5(ver)
    if 'rhel4' in actions:
        build_rhel4(ver)
    if 'suse' in actions:
        build_suse(ver)
    if 'mac' in actions:
        build_mac()
    # 
    # restore simuPOP.release
    if release == 'snapshot':
        writeReleaseFile(old_ver, old_rev)

    
