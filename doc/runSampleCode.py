#!/usr/bin/env python
#
# This script get sys.argv[1], open it and run as if in a 
# interactive session. Then, the output is separated accorindg
# to the instructions within.
#
# 
import code, sys, os

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
output = open(sys.argv[1]+'.out', 'w')

# set stdin, stderr, stdout
sys.stdin = wrapper(sys.stdin)
#sys.stderr = sys.stdout
sys.stderr = output
sys.stdout = output


b = runScriptInteractively(filename = sys.argv[1])
b.interact(None)

output.close()

# separate this file using a perl script
perlscript = open(sys.argv[1]+'.pl', 'w')
perlscript.write('''\
#!/usr/bin/perl
$output = 0;
$first = 0;
$file = "";
while(<>)
{
  if ( /^(>>>|\.\.\.) #file/ )
  {
    ($file) = $_ =~ /.*file\s*(.*)\s*/;
    open(FILE, ">$file");
    $first = 1;
    $output = 1;
    next;
  }
  elsif ( /^>>> #end/ )
  {
    $first = 0;
    $output = 0;
    print FILE ">>>  \n";
    close(FILE);
    next;
  }
  elsif (/^(>>>|\.\.\.)\s*#PS / )
  {
    ($pre, $cmd) = $_ =~ /^(>>>|\.\.\.)\s*#PS (.*)$/;
    system("$cmd");
  }
  if( $output )
  {
    if( $first )
    {
      print FILE ">>> ", substr($_, 4);
      $first = 0;
    }
    else
    {
      print FILE "$_";
    }
  }
}
close(FILE);
''')

perlscript.close()

# run this perl script
os.system('perl ' + sys.argv[1]+'.pl' + ' ' + sys.argv[1]+'.out')

os.remove(sys.argv[1]+'.pl')
os.remove(sys.argv[1]+'.out')
