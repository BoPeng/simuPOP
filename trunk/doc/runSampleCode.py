#!/usr/bin/env python
#
# $File: runSampleCode.py $
# $LastChangedDate: 2009-09-05 13:56:43 -0500 (Sat, 05 Sep 2009) $
# $Rev: 2894 $
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

#
# This script get a filename from sys.argv[1], open and run it as if in an
# interactive session. The output is separated according to special comments
# in the output file. Allowed intructions are
#
# * logging output between these two lines to a filename
#   #file filename
#   #end
#
# * execute, but do not write the output between these two lines
#   #begin_ignore
#   #end
# 
# * expect error so do not stop when an error happens
#   #expecterror

import code, sys, os, re, tempfile

#  run a script interatively
class runScriptInteractively(code.InteractiveConsole):
    def __init__(self, locals=None, filename="<console>", file = None):
        self.file = file or open(filename)
        code.InteractiveConsole.__init__(self, locals, filename)

    def raw_input(self, prompt):
        l = self.file.readline()
        if l == '':
            raise EOFError
        sys.stdout.write(prompt + l)
        return l.strip("\n")
    
    def exit(self):
        self.file.close()

# do not stop on help() function since stdin is then no longer
# a terminal.
class wrapper:
  def __init__(self, file):
    self.file = file
  def isatty(self):
    return 0
  def __getattr__(self, key):
    return getattr(self.file, key)

def runScript(inputFile, outputFile):
    '''Run a script and return its output as a list of strings'''
    # out to a file
    #
    oldIn = sys.stdin
    oldOut = sys.stdout
    oldErr = sys.stderr
    #
    # set stdin, stderr, stdout
    outFile = open(outputFile, 'w')
    sys.stdin = wrapper(sys.stdin)
    sys.stderr = outFile
    sys.stdout = outFile
    #
    b = runScriptInteractively(locals=locals(), filename = inputFile)
    b.interact(None)
    b.exit()
    #
    # reset io streams
    sys.stdin = oldIn
    sys.stdout = oldOut
    sys.stderr = oldErr
    #
    outFile.close()

def writeFile(content, srcFile, logFile=False):
    dir = os.path.split(srcFile)[0]
    if dir != '' and not os.path.isdir(dir):
        os.mkdir(dir)
    #
    src = open(srcFile, 'w')
    ignore = False
    expect_error = False
    start = not logFile
    for line in content:
        if line.startswith('#begin_file') or line.startswith('>>> #begin_file'):
            continue
        if logFile and line.startswith('>>>'):
            start = True
        if not start:
            continue
        if line.startswith('#begin_ignore') or line.startswith('>>> #begin_ignore'):
            ignore = True
        elif line.startswith('#end_ignore') or line.startswith('>>> #end_ignore'):
            ignore = False
        elif line.startswith('#expect_error') or line.startswith('>>> #expect_error'):
            expect_error = True
        elif not ignore:
            src.write(line)
    # if there is error
    if not expect_error and True in ['Error' in x for x in content]:
        print()
        print('An Error occured in log file %s ' % srcFile)
        print("If this is expected, please add '#expecterror' in your source code.")
        print()
        print(''.join(content))
        print()
        sys.exit(1)
    src.close()


def runSampleCode(srcFile, names):
    begin_re = re.compile('^#begin_file\s*(.*)')
    end_re = re.compile('^#end_file\s*$')
    #
    src = open(srcFile, 'r')
    tmpSrc = tmpSrcName = None
    filename = None
    count = 0
    skip = False
    for lineno, line in enumerate(src.readlines()):
        if begin_re.match(line):
            if tmpSrc is not None:
                print('ERROR (Unmatched file/end at line %s): %s' % (lineno, line))
                sys.exit(1)
            filename = begin_re.match(line).groups()[0].strip()
            if len(names) > 0 and not (True in [name in filename for name in names]):
                skip = True
                continue
            else:
                skip = False
            tmp, tmpSrcName = tempfile.mkstemp()
            os.close(tmp)
            tmpSrc = open(tmpSrcName, 'w')
            count += 1
        if skip:
            continue
        if end_re.match(line):
            if tmpSrc is None:
                print('ERROR (Unmatched file/end at line %d): %s' % (lineno, line))
                sys.exit(1)
            #
            print('Processing %s...' % filename)
            sys.stdout.flush()
            tmpSrc.close()
            tmpSrc = None
            tmp, tmpLogName = tempfile.mkstemp()
            os.close(tmp)
            os.system('%s %s %s %s' % (sys.executable, sys.argv[0], tmpSrcName, tmpLogName))
            #
            writeFile(open(tmpSrcName).readlines(), filename)
            logFile = filename.replace('.py', '.log')
            writeFile(open(tmpLogName).readlines(), logFile, True)
            #
            os.remove(tmpSrcName)
            os.remove(tmpLogName)
        else:
            if tmpSrc is not None:
                tmpSrc.write(line)
            elif line.strip() != '' and not line.startswith('#'):
                print('Unprocessed:', line)
    src.close()
    print('Finished processing %d examples.' % count)

if __name__ == '__main__':
    if len(sys.argv) == 3 and os.path.isfile(sys.argv[1]):
        runScript(sys.argv[1], sys.argv[2])
    else:
        runSampleCode('userGuide.py', sys.argv[1:])

