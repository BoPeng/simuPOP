"""
simuPOP utilities.
This module provides some commonly used operators
and some utility class for plotting.
"""

from scipy import gplt
from simuPOP import pyEval
from simuUtil import Aggregator

class VarPlotter:
  """
  plot given variable. 
  
  """
  skip = 0
  title = ""
  xtitle = ""
  ytitle = ""
  legend = []
  
  def __init__(self, win=0, update=1, title="plotting variable",
               xtitle="generation", ytitle="", legend=[], saveAs=""):
    """
    win: window size. I.e., maximum generations of data to plot
    update: update figure every 'update' calls.
    title: plot title
    xtitle, ytitle, x,y titles
    legend, an array of strings. If size is one, 1,2 etc
      will be appended to other names.
    """
    self.data = Aggregator(win)
    self.update = update
    self.skip = 0
    self.title = title
    self.xtitle = xtitle
    self.ytitle = ytitle
    self.legend = legend
    self.saveAs = saveAs
  
  def __del__(self) :
    gplt.close()

  def save(self, filename, format="png"):
    gplt.output(filename, format)
    
  def plot(self, gen, data ):
    # Does passed data match legend length?
    if len(self.legend) == 0:
      for i in range( len(data) ):
        self.legend.append( 'var' + str(i) )
    elif len(self.legend) == 1:
      for i in range(1, len(data)):
        self.legend.append( self.legend[0] + str(i) )
    
    if len(self.legend) != len(data):
      raise ValueError, "Data and legend size should equal"

    self.data.push(gen, data)
    
    # plot data
    if len( self.data.gen) > 1:
      self.skip += 1
      if self.update == 0 or self.skip == self.update:
        self.skip = 0
    
        # generate command. The most difficulty lies
        # in the proper setting of legends.
        cmd = 'gplt.plot(self.data.gen, self.data.data[0],\'t "' \
              + self.legend[0] + '" with lines\''
        for i in range( 1, len( self.legend ) ):
          cmd += ', self.data.gen, self.data.data[' + str(i)  \
            + '], \'t "' + self.legend[i] + '" with lines\''
        cmd += ')'
        # print cmd
        eval(cmd)

        gplt.title(self.title)
        gplt.xtitle(self.xtitle)
        gplt.ytitle(self.ytitle)
  
        if self.saveAs != "":
          self.save(self.saveAs)
          
# wrapper of VarPlotter class
def varPlotter(expr, win=0, update=1, 
               title="variable", xtitle="generation", ytitle="", 
               legend=[], saveAs="", *args, **kwargs):
  # deal with additional arguments
  parm = ''
  for (k,v) in kwargs.items():
    parm += str(k) + '=' + str(v) + ', '
  #
  # pyExec( preStmts = 'rpv=SimuPlotter(options....)',
  #   stmts = 'rpv.plot(gen,rep,expr', ... other options... )
  cmd = ('''pyEval(exposePop=1, %s \
      preStmts= \"\"\"spv=VarPlotter(win=%d, update=%d, 
      title=\'\'\'%s\'\'\', xtitle=\'\'\'%s\'\'\',  
      ytitle=\'\'\'%s\'\'\', legend=''' + str(legend) + ''', 
      saveAs=\'\'\'%s\'\'\')\"\"\",
      stmts=r\"\"\"spv.plot(gen, eval(\'\'\'%s\'\'\'))\"\"\")''') \
      % ( parm, win, update, title, xtitle, ytitle, saveAs, expr) 
  #print cmd
  return eval(cmd)
