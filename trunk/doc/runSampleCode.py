#!/usr/bin/env python
#
# This script get sys.argv[1], open it and run as if in a 
# interactive session. Then, the output is separated accorindg
# to the instructions within.
#
# 
import code, sys, os, re



#  run a script interatively
class runScriptInteractively(code.InteractiveConsole):
    def __init__(self, locals=None, filename="<console>", file = None):
        self.file = file or open(filename)
        code.InteractiveConsole.__init__(self, locals, filename)

    def raw_input(self, prompt):
        l = self.file.readline()
        if l == '': raise EOFError
        sys.stdout.write(prompt + l)
        return l.strip("\n")

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
    # out to a file
    outFile = open(outputFile, 'w')
    
    oldIn = sys.stdin
    oldOut = sys.stdout
    oldErr = sys.stderr
    
    # set stdin, stderr, stdout
    sys.stdin = wrapper(sys.stdin)
    sys.stderr = outFile
    sys.stdout = outFile
    
    b = runScriptInteractively(locals=locals(), filename = sys.argv[1])
    b.interact(None)
    
    # reset io streams
    sys.stdin = oldIn
    sys.stdout = oldOut
    sys.stderr = oldErr
    
    outFile.close()
    
    print "Finished executing ", sys.argv[1]
    print sys.argv


def splitFile(outputFile, runCommand=True):
    if runCommand:
        begin_re = re.compile('^(>>>|\.\.\.)\s*#file\s*(.*)')
        end_re = re.compile('^(>>>|\.\.\.)\s*#end')
        cmd_re = re.compile('^(>>>|\.\.\.)\s*#PS\s*(.*)')
    else:
        begin_re = re.compile('^()*#file\s*(.*)')
        end_re = re.compile('^()*#end')
        cmd_re = re.compile('^()*#PS\s*(.*)')
    outFile = open(outputFile, 'r')
    out = ''
    first = False
    ignore = False
    for line in outFile.readlines():
        if line.strip().endswith('#beginignore'):
            ignore = True
        elif line.strip().endswith('#endignore'):
            ignore = False
            continue
        if ignore:
            continue
        if begin_re.match(line):
            (tmp, file) = begin_re.match(line).groups()
            file = file.strip()
            # in a special split input mode
            if runCommand: 
                print "Writing log to %s" % file
            else:
                file = file.replace('.log', '.py')
                print "Writing source to %s" % file
            dir = os.path.split(file)[0]
            if dir != '' and not os.path.isdir(dir):
                os.mkdir(dir)
            out = open(file, 'w')
            first = True
        elif end_re.match(line):
            if type(out) != type(''):
                if runCommand:
                    print >> out, '>>>'
                out.close()
                out = ""
        elif runCommand and cmd_re.match(line):
            (tmp, cmd) = cmd_re.match(line).groups()
            print "Running command %s" % cmd
            os.system(cmd)
        else:
            if type(out) == type(''):
                continue
            if first:
                if runCommand:
                    print >> out, '>>> %s' % line[4:],
                    if 'Error' in line:
                        print >> sys.stderr, " * * * * * Error is encounted in this file (this is perhaps intentional)"
                else:
                    print >> out, '#!/usr/bin/env python'
                    print >> out, 'from simuPOP import *'
                    print >> out, 'rng().setSeed(12345)'
                    print >> out
                    print >> out, line,
                first = False
            else:
                print >> out, line,
                if 'Error' in line:
                    print >> sys.stderr, " * * * * * Error is encounted in this file (this is perhaps intentional)"
    outFile.close()

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] == '-h':
        print 'Usage: runSampleCode scriptToRun'
        print '    -h: view this help information'
        sys.exit(0)

    inputFile = sys.argv[-1]
    splitFile(inputFile, False)
    outputFile = inputFile.split('.')[0] + '.out'
    runScript(inputFile, outputFile)
    splitFile(outputFile, True)
    os.remove(outputFile)

