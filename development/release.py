#!/usr/bin/env python

# make src and binary distribution on all platforms 
# (currently mac, linux and solaris) as daily snapshot
#
# This is part of simuPOP
# 
# Bo Peng (bpeng@mdanderson.org)
#

import os
import sys
import time
import shutil
import platform
import subprocess
import argparse
import glob
import platform
import re

download_directory = 'download'
MODULES = ['std', 'op', 'la', 'laop', 'ba', 'baop', 'mu', 'muop', 'lin', 'linop']

#
# Utility
#
def run_command(cmd):
    print('Running: {}'.format(cmd))
    return subprocess.call(cmd, shell=True)

#
# Actions
#
def setVersionRevision(release):
    #
    # get revisionision of local tree
    revision = subprocess.check_output('svnversion .', shell=True).strip()
    if ':' in revision:
        revision = revision.split(':')[0]
    if revision.endswith('M'):
        revision = revision[:-1]
        print 'Warning: Please commit all changes before releasing a source package'
    #
    # get recorded release
    ver = {}
    execfile('simuPOP_version.py', globals(), ver)
    if release is None:
        release = ver['SIMUPOP_VER']
    # 
    # write release file
    with open('simuPOP_version.py', 'w') as version_file:
        version_file.write('SIMUPOP_VER="{}"\nSIMUPOP_REV="{}"\n'
            .format(release, revision))
    # 
    # update documents
    for filename in ['doc/userGuide.lyx', 'doc/refManual.lyx']:
        with open(filename) as manual:
            content = manual.readlines()
            for idx,line in enumerate(content):
                if line.startswith('\\setreleaseinfo'):
                    content[idx] = '\\setreleaseinfo{Release %s (\\mbox{Rev: %s})}\n' % (release, revision)
                    break
        with open(filename, 'w') as manual:
            manual.write(''.join(content))
    return release, revision


def prepareEnvironment():
    ''' remove unnecessary files '''
    # walk down the file hierarchy
    for root,path,files in os.walk('.'):
        # check file type
        for file in files:
            if True in [('_%s.cpp' % x in file) for x in MODULES]:
                # force the use of / since miktex uses / even under windows
                print "Remove file %s..." % os.path.join(root, file)
                os.remove(os.path.join(root, file))

def generateSWIGWrappers():
    '''Generate Wrapping files'''
    SWIG = 'swig'
    SWIG_CPP_FLAGS = '-O -templatereduce -shadow -python -c++ -keyword '\
        '-nodefaultctor -w-503,-312,-511,-362,-383,-384,-389,-315,-509,-525 -Ibuild'
    SWIG_CC_FLAGS = '-python -keyword'
    if sys.version_info[0] == 3:
        SWIG_CPP_FLAGS += ' -py3'
        SWIG_CC_FLAGS += ' -py3'
    SWIG_RUNTIME_FLAGS = '-python -external-runtime'
    SWIG_OUTDIR = 'src'
    WRAP_INFO = {
        'std':    ['src/simuPOP_std_wrap.cpp', 'src/simuPOP_std.i', ''],
        'op':     ['src/simuPOP_op_wrap.cpp', 'src/simuPOP_op.i', '-DOPTIMIZED'],
        'la':     ['src/simuPOP_la_wrap.cpp', 'src/simuPOP_la.i', '-DLONGALLELE'],
        'laop':   ['src/simuPOP_laop_wrap.cpp', 'src/simuPOP_laop.i', '-DLONGALLELE -DOPTIMIZED'],
        'ba':     ['src/simuPOP_ba_wrap.cpp', 'src/simuPOP_ba.i', '-DBINARYALLELE'],
        'baop':   ['src/simuPOP_baop_wrap.cpp', 'src/simuPOP_baop.i', '-DBINARYALLELE -DOPTIMIZED'],
        'mu':     ['src/simuPOP_mu_wrap.cpp', 'src/simuPOP_mu.i', '-DMUTANTALLELE'],
        'muop':   ['src/simuPOP_muop_wrap.cpp', 'src/simuPOP_muop.i', '-DMUTANTALLELE -DOPTIMIZED'],
        'lin':    ['src/simuPOP_lin_wrap.cpp', 'src/simuPOP_lin.i', '-DLINEAGE'],
        'linop':  ['src/simuPOP_linop_wrap.cpp', 'src/simuPOP_linop.i', '-DLINEAGE -DOPTIMIZED'],
    }
    fout = subprocess.Popen(SWIG + ' -version', shell=True, stdout=subprocess.PIPE).stdout
    output = fout.readlines()[1].decode('utf8')
    #
    try:
        version = re.match('SWIG Version\s*(\d+).(\d+).(\d+).*', output).groups()
    except:
        print('Can not obtain swig version, please install swig')
        sys.exit(1)
    v1, v2, v3 = [int(x) for x in version]
    if (v1, v2, v3) < (1, 3, 35):
        print('Swig >= 1.3.35 is required, please upgrade it.')
        sys.exit(1)
    if sys.version_info[0] >= 3 and sys.version_info[1] >= 2 and \
        (v1, v2, v3) < (2, 0, 4):
        print('Swig >= 2.0.4 is required for Python 3.2 or higher')
        sys.exit(1)
    # generate header file 
    print("Generating external runtime header file src/swigpyrun.h...")
    run_command('%s %s src/swigpyrun.h' % (SWIG, SWIG_RUNTIME_FLAGS))
    # try the first option set with the first library
    for lib in MODULES:
        print("Generating wrapper file " + WRAP_INFO[lib][0])
        if run_command('%s %s -outdir %s %s -o %s %s' % (SWIG, SWIG_CPP_FLAGS, \
            SWIG_OUTDIR, WRAP_INFO[lib][2], WRAP_INFO[lib][0], WRAP_INFO[lib][1])) != 0:
            print("Calling swig failed. Please check your swig version.")
            sys.exit(1)
    print("Generating wrapper file src/gsl_wrap.c")
    if run_command('%s %s -outdir %s %s -o %s %s' % (SWIG, SWIG_CC_FLAGS, \
        SWIG_OUTDIR, '', 'src/gsl_wrap.c', 'src/gsl.i')) != 0:
        print("Calling swig failed. Please check your swig version.")
        sys.exit(1)
    print("\nAll wrap files are generated successfully.\n")


def generateDocuments(ver, rev):
    ' build doxygen document '
    d = os.getcwd()
    os.environ['SIMUPOP_DOC_DIR'] = 'doc'
    os.environ['SIMUPOP_VER'] = ver
    os.environ['SIMUPOP_REV'] = rev
    run_command('doxygen')
    os.chdir('development')
    run_command('python doxy2swig.py')
    os.chdir(d)


def buildSimuPOP():
    'Build simuPOP'
    run_command('python setup.py install')


def generateExamples():
    d = os.getcwd()
    os.chdir('doc')
    run_command('python runSampleCode.py')
    os.chdir(d)

def uploadSourceDocuments(ver, rev):
    ' build doxygen document and update source document to sourceforge'
    os.environ['SIMUPOP_DOC_DIR'] = 'doc'
    os.environ['SIMUPOP_VER'] = ver
    os.environ['SIMUPOP_REV'] = rev
    run_command('doxygen Doxy_web')
    d = os.getcwd()
    os.chdir('doc')
    run_command('make src_doc')
    os.chdir(d)

def uploadDocuments(ver, rev):
    d = os.getcwd()
    os.chdir('doc')
    run_command('make pdf')
    run_command('make dist_release')
    os.chdir(d)
    


def buildSourcePackage(ver):
    # build source
    run_command('python setup.py  sdist --formats=gztar,zip')
    # coppy files
    shutil.copy('dist/simuPOP-%s.tar.gz' % ver, '%s/simuPOP-%s-src.tar.gz' % (download_directory, ver))
    print('Source package saved to {}/simuPOP-{}-src.tar.gz'.format(download_directory, ver))
    shutil.copy('dist/simuPOP-%s.zip' % ver, '%s/simuPOP-%s-src.zip' % (download_directory, ver))
    print('Source package saved to {}/simuPOP-{}-src.zip'.format(download_directory, ver))


def build_x86_64(ver):
    src_path = os.path.realpath(download_directory)
    source = os.path.join(src_path, 'simuPOP-%s-src.tar.gz' % ver)
    if not os.path.isfile(source):
        print 'Source package %s does not exist. Please run "build.py src" first' % source
        sys.exit(1)
    # 
    # build
    d = os.getcwd()
    os.chdir(user_tmp_directory)
    print 'Copying source package to user temp directory...'
    run_command('/bin/rm -rf simuPOP-%s' % ver)
    shutil.copy(source, 'simuPOP-%s-src.tar.gz' % ver)
    print 'Unpacking ...'
    run_command('tar zxf simuPOP-%s-src.tar.gz' % ver)
    os.chdir('simuPOP-%s' % ver)
    print 'Building ...'
    run_command('python setup.py  bdist --formats=gztar,rpm')
    # coppy files
    shutil.copy('dist/simuPOP-%s-1.x86_64.rpm' % ver, '%s/simuPOP-%s-1.rhel4.x86_64.rpm' % (src_path, ver))
    shutil.copy('dist/simuPOP-%s.linux-x86_64.tar.gz' % ver, '%s/simuPOP-%s.rhel4.x86_64.tar.gz' % (src_path, ver))
    shutil.copy('dist/simuPOP-%s-1.src.rpm' % ver, '%s/simuPOP-%s-1.src.rpm' % (src_path, ver))
    os.chdir(d)


def build_vm(ver, name, pyver, vm, vm_port, vm_name):
    # set display, since this script will be launched from crontab.
    run_command('vmrun start %s' % vm)
    # wait for the vm to start.
    run_command('sleep 30')
    run_command('ssh -p %d %s "/bin/rm -rf simuPOP-%s.zip simuPOP-%s.tar.gz simuPOP-%s"' % \
        (vm_port, vm_name, ver, ver, ver))
    run_command('scp -P %s %s/simuPOP-%s-src.tar.gz %s:' % \
        (vm_port, download_directory, ver, vm_name))
    #
    run_command('ssh -X -p %d %s "tar zxf simuPOP-%s-src.tar.gz && cd simuPOP-%s && python setup.py  bdist --formats=gztar,rpm"' % \
        (vm_port, vm_name, ver, ver))
    run_command('scp -P %d %s:simuPOP-%s/dist/simuPOP-%s.linux-i686.tar.gz %s/simuPOP-%s-%s-py%2d.tar.gz' % \
        (vm_port, vm_name, ver, ver, download_directory, ver, name, pyver))
    run_command("scp -P %d %s:simuPOP-%s/dist/simuPOP-%s-i386.rpm %s/simuPOP-%s-%s-py%2d.i386.rpm" % \
        (vm_port, vm_name, ver, ver, download_directory, ver, name, pyver))
    run_command('vmrun suspend %s' % vm)


def build_mdk(ver):
    build_vm(ver, 'mdk', 23, mdk_vm, mdk_port, mdk_vm_name)

def build_fedora5(ver):
    build_vm(ver, 'fedora5', 23, fedora5_vm, fedora5_port, fedora5_vm_name)

def build_rhel4(ver):
    build_vm(ver, 'rhel4', 23, rhel4_vm, rhel4_port, rhel4_vm_name)

def build_suse(ver):
    build_vm(ver, 'suse', 23, suse_vm, suse_port, suse_vm_name)

def build_remote(ver, remote_machine):
    source = os.path.realpath('%s/simuPOP-%s-src.tar.gz' % (download_directory, ver))
    if not os.path.isfile(source):
        print 'Source package %s does not exist. Please run "build.py src" first' % source
        sys.exit(1)
    #   
    print 'Copying source package to remote machine ...'
    run_command("ssh %s '/bin/rm -rf temp &&  mkdir temp && /bin/rm -rf simuPOP'" % remote_machine)
    run_command('scp %s/simuPOP-%s-src.tar.gz %s:temp' % (download_directory, ver, remote_machine))
    #
    print 'Building ...'
    unpack = 'tar zxf simuPOP-%s-src.tar.gz' % ver
    build = 'python setup.py bdist_dumb'
    run_command('ssh -X %s "cd temp && %s && cd simuPOP-%s && %s"' % (remote_machine, unpack, ver, build))
    run_command('scp %s:temp/simuPOP-%s/dist/* %s' % (remote_machine, ver, download_directory))


def createMacPackage(ver, pyver):
    #
    if not platform.platform().startswith('Darwin'):
        sys.exit('Can only build darwin binary from a mac')
    #
    source = os.path.realpath('%s/simuPOP-%s-src.tar.gz' % (download_directory, ver))
    if not os.path.isfile(source):
        print 'Source package %s does not exist. Please run "build.py src" first' % source
        build_src(ver)
    #
    # copy to a directory
    temp_dir = os.path.expanduser('~/Temp/mpkg')
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
    # 
    src_dir = os.path.join(temp_dir, 'simuPOP-%s' % ver)
    if os.path.isdir(src_dir):
        shutil.rmtree(src_dir)
    #
    # copy files
    shutil.copy(source, temp_dir)
    #
    old_dir = os.getcwd()
    os.chdir(temp_dir)
    # decompress
    run_command('tar -zxf simuPOP-%s-src.tar.gz' % ver)
    # 
    os.chdir(src_dir)
    run_command('python3 setup.py bdist')
    os.chdir(old_dir)
    # 
    # Move results to download directory
    dest_targz = os.path.join(download_directory, 'simuPOP-{}-py{}.tar.gz'.format(ver, pyver))
    if os.path.isdir(dest_targz):
        shutil.rmtree(dest_targz)
    print('Moving tar.gz package from %s/dist to %s' % (src_dir, dest_targz))
    shutil.move(glob.glob(os.path.join(src_dir, 'dist', 'simuPOP-{}*.tar.gz'.format(ver)))[0],
        dest_targz)
    # for python 2, also build mpkg
    if pyver.startswith('2'):
        os.chdir(src_dir)
        run_command('python setup.py bdist_mpkg')
        os.chdir(old_dir)
        dest_mpkg = os.path.join(download_directory, 'simuPOP-{}-py{}.mpkg'.format(ver, pyver))
        if os.path.isdir(dest_mpkg):
            shutil.rmtree(dest_mpkg)
        print('Moving mpkg package from %s/dist to %s' % (src_dir, dest_mpkg))
        shutil.move(glob.glob(os.path.join(src_dir, 'dist', 'simuPOP-{}*.mpkg'.format(ver)))[0],
            dest_mpkg)

def createMacImage(ver, pyver):
    #
    mpkg = os.path.join(download_directory, 'simuPOP-{}-py{}.mpkg'.format(ver, pyver))
    if not os.path.isdir(mpkg):
        print('mpkg is not available, building a mpkg package')
        buildMacPackage(ver, pyver)
    #
    # create a dmg file for the package
    temp_dir = os.path.expanduser('~/Temp')
    dest = os.path.join(temp_dir, 'simuPOP-{}-py{}'.format(ver, pyver))
    if os.path.isdir(dest):
        shutil.rmtree(dest)
    #
    print('Copy %s to %s' % (mpkg, dest))
    shutil.copytree(mpkg, '{0}/simuPOP-{1}-py{2}.mpkg'.format(dest, ver, pyver))
    #
    dmg = os.path.join(download_directory, 'simuPOP-{}-py{}.dmg'.format(ver, pyver))
    if os.path.isfile(dmg):
        os.remove(dmg)
    run_command('hdiutil create {0} -volname simuPOP-{1} -fs HFS+ -srcfolder {2}/simuPOP-{1}-py{3}.mpkg'
        .format(dmg, ver, dest, pyver))


def condaRelease():
    '''Build conda binary'''
    d = os.getcwd()
    os.chdir('development')
    run_command('conda build conda')
    os.chdir(d)

def tagRelease(release):
    ''' if there are changes, commit it '''
    if subprocess.check_output('svn diff', shell=True) != '':
        cmd = 'svn ci -m "automatic checkin on %s"' % time.asctime()
        run_command(cmd)
    run_command('svn update')
    cmd = ('svn copy https://sourceforge.net/p/simupop/code/HEAD/tree/trunk ' + 
         'https://sourceforge.net/p/simupop/code/HEAD/tree/trunk/tag/v%s') % release + \
        ' -m "Version %s released at %s"' % (release, time.asctime())
    print cmd
    run_command(cmd)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Create source distribution 
        and binary installers for a simuPOP release. In addition to optional
        parameters version and tag, extra parameters would be specified and 
        will be passed directly to the 'python setup.py install' process. ''')
    parser.add_argument('--version',
        help='''Modify simuPOP_version.py to the specified version string and
            make the release.''')
    parser.add_argument('actions', nargs='*', default=[
        'build', 'src', 'doc', 'mac'],
        help='Actions to take to make the release.')
    # go to the top source directory        
    os.chdir('..')
    #
    args, argv = parser.parse_known_args()
    #
    # get revision number and update last_revision_file
    ver, rev = setVersionRevision(args.version)
    print 'New version: %s, revision: %s ' % (ver, rev)
    #
    if 'build' in args.actions:
        prepareEnvironment()
        generateSWIGWrappers()
        # documents need to be generated before simuPOP is built
        generateDocuments(ver, rev)
        buildSimuPOP()
    if 'src' in args.actions:
        buildSourcePackage(ver)
    if 'doc' in args.actions:
        generateDocuments(ver, rev)
        generateExamples()
        uploadDocuments(ver, rev)
        uploadSourceDocuments(ver, rev)
    if 'src_doc' in args.actions:
        uploadSourceDocuments(ver, rev)
    if 'mac' in args.actions and platform.platform().startswith('Darwin'):
        createMacPackage(ver, pyver='2.7')
        createMacPackage(ver, pyver='3.3')
        createMacImage(ver, pyver='2.7')
    if 'tag' in args.actions:
        tagRelease(ver)
    if 'conda' in args.actions:
        condaRelease()
