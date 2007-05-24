#!/usr/bin/env python
#
# This script get sys.argv[1], open it and run as if in a 
# interactive session. Then, the output is separated accorindg
# to the instructions within.
#
# 
import code, sys, os, re

if len(sys.argv) < 2 or sys.argv[1] == '-h':
    print 'Usage: runSampleCode scriptToRun'
    sys.exit(0)

outputFile = sys.argv[1].split('.')[0] + '.out'
perlFile = sys.argv[1].split('.')[0] + '.pl'

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

# out to a file
outFile = open(outputFile, 'w')

oldIn = sys.stdin
oldOut = sys.stdout
oldErr = sys.stderr

# set stdin, stderr, stdout
sys.stdin = wrapper(sys.stdin)
#sys.stderr = sys.stdout
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

begin_re = re.compile('^(>>>|\.\.\.)\s*#file\s*(.*)')
end_re = re.compile('^>>>\s*#end')
cmd_re = re.compile('^(>>>|\.\.\.)\s*#PS\s*(.*)')
outFile = open(outputFile, 'r')
out = ''
first = False
for line in outFile.readlines():
    if begin_re.match(line):
        (tmp, file) = begin_re.match(line).groups()
        file = file.strip()
        print "Writing to %s" % file
        out = open(file, 'w')
        first = True
    elif end_re.match(line):
        if type(out) != type(''):
            print >> out, '>>>'
            out.close()
            out = ""
    elif cmd_re.match(line):
        (tmp, cmd) = cmd_re.match(line).groups()
        print "Running command %s" % cmd
        os.system(cmd)
    else:
        if type(out) == type(''):
            continue
        if first:
            print >> out, '>>> %s' % line[4:],
            first = False
        else:
            print >> out, line,
    
