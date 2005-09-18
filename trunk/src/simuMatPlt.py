"""
simuPOP plotting with matplotlib


"""

import Numeric
import matplotlib
matplotlib.use('TkAgg')
import pylab


pylab.plot([1,2,3,4])
pylab.show

class varPlotter:
  """
  plot given variable. 
  
  """
  y_data = Numeric.array([])
  x_gen = Numeric.array([])
  skip = 0
  title = ""
  xtitle = ""
  ytitle = ""
  legend = []
  
  def __init__(self, win=0, update=1, title="plotting variable",
               xtitle="generation", ytitle="", legend=[]):
    """
    win: window size. I.e., maximum generations of data to plot
    update: update figure every 'update' calls.
    title: plot title
    xtitle, ytitle, x,y titles
    legend, an array of strings. If size is one, 1,2 etc
      will be appended to other names.
    """
    self.win = win
    self.update = update
    self.skip = 0
    self.title = title
    self.xtitle = xtitle
    self.ytitle = ytitle
    self.legend = legend
  
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

    # first add data to allData
    if len(self.y_data) == 0:
      self.y_data = Numeric.array([ data ])
    else:
      self.y_data = Numeric.concatenate( (self.y_data, [tuple(data)]))
    if len(self.x_gen) == 0:
      self.x_gen = Numeric.array([ gen ])
    else:
      self.x_gen = Numeric.concatenate( (self.x_gen, [ gen ] ))

    # trim data if necessary
    if self.win > 0 :
      if len(self.x_gen) > self.win:
        self.x_gen = self.x_gen[1:]
        self.y_data = self.y_data[1:]

    # plot data
    if len( self.y_data) > 1:
      self.skip += 1
      if self.update == 0 or self.skip == self.update :
        self.skip = 0
    
        # generate command. The most difficulty lies
        # in the proper setting of legends.
        cmd = 'gplt.plot(self.x_gen, self.y_data[:,0],\'t "' \
              + self.legend[0] + '" with lines\''
        for i in range( 1, len( self.legend ) ):
          cmd += ', self.x_gen, self.y_data[:,' + str(i)  \
            + '], \'t "' + self.legend[i] + '" with lines\''
        cmd += ')'
        eval(cmd)

        gplt.title(self.title)
        gplt.xtitle(self.xtitle)
        gplt.ytitle(self.ytitle)
  
