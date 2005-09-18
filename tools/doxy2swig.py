#!/usr/bin/env python
"""Doxygen XML to SWIG docstring converter.

Converts Doxygen generated XML files into a file containing docstrings
that can be used by SWIG-1.3.x.  Note that you need to get SWIG
version > 1.3.23 or use Robin Dunn's docstring patch to be able to use
the resulting output.

Usage:

 doxy2swig.py input.xml output.i

input.xml is your doxygen generated XML file and output.i is where the
output will be written (the file will be clobbered).

"""

# This code is implemented using Mark Pilgrim's code as a guideline:
#   http://www.faqs.org/docs/diveintopython/kgp_divein.html
#
# Author: Prabhu Ramachandran
# Modified by: Bo Peng
# License: BSD style

from xml.dom import minidom
import re
import textwrap
import sys
import types
import os.path

def my_open_read(source):
   if hasattr(source, "read"):
       return source
   else:
       return open(source)

def my_open_write(dest):
   if hasattr(dest, "write"):
       return dest
   else:
       return open(dest, 'w')

class Doxy2SWIG:
   """Converts Doxygen generated XML files into a file containing
   docstrings that can be used by SWIG-1.3.x that have support for
   feature("docstring").  Once the data is parsed it is stored in
   self.pieces.

   """

   def __init__(self, src):
       """Initialize the instance given a source object (file or
       filename).

       """
       f = my_open_read(src)
       self.my_dir = os.path.dirname(f.name)
       self.xmldoc = minidom.parse(f).documentElement
       f.close()

       self.pieces = []
       self.pieces.append('\n// File: %s\n'%\
                          os.path.basename(f.name))

       self.space_re = re.compile(r'\s+')
       self.lead_spc = re.compile(r'^(%feature\S+\s+\S+\s*?)"\s+(\S)')
       self.multi = 0
       self.indent = 0
       self.maxChar = 70
       self.curCol = 0
       self.ignores = ('inheritancegraph', 'param', 'listofallmembers',
                       'innerclass', 'name', 'declname', 'incdepgraph',
                       'invincdepgraph', 'programlisting', 'type',
                       'references', 'referencedby', 'location',
                       'collaborationgraph', 'reimplements',
                       'reimplementedby', 'derivedcompoundref',
                       'basecompoundref', 'header', 'includes')
       #self.generics = []

   def generate(self):
       """Parses the file set in the initialization.  The resulting
       data is stored in `self.pieces`.

       """
       self.parse(self.xmldoc)

   def parse(self, node):
       """Parse a given node.  This function in turn calls the
       `parse_<nodeType>` functions which handle the respective
       nodes.

       """
       # print node.__class__.__name__
       pm = getattr(self, "parse_%s" %node.__class__.__name__)
       pm(node)

   def parse_Document(self, node):
       self.parse(node.documentElement)

   def parse_Text(self, node):
       txt = node.data
       txt = txt.replace('\\', r'\\\\')
       txt = txt.replace('"', r'\"')
       # ignore pure whitespace
       m = self.space_re.match(txt)
       if m and len(m.group()) == len(txt):
           pass
       else:
           self.add_text( self.wrap_text(txt, self.curCol - self.indent))

   def parse_Element(self, node):
       """Parse an `ELEMENT_NODE`.  This calls specific
       `do_<tagName>` handers for different elements.  If no handler
       is available the `generic_parse` method is called.  All
       tagNames specified in `self.ignores` are simply ignored.

       """
       name = node.tagName
       ignores = self.ignores
       if name in ignores:
           return
       attr = "do_%s" % name
       if hasattr(self, attr):
           handlerMethod = getattr(self, attr)
           handlerMethod(node)
       else:
           self.generic_parse(node)
           #if name not in self.generics: self.generics.append(name)

   def add_text(self, value):
       """Adds text corresponding to `value` into `self.pieces`."""
       if type(value) in (types.ListType, types.TupleType):
           self.pieces.extend(value)
       else:
           self.pieces.append(value)

   def get_specific_nodes(self, node, names):
       """Given a node and a sequence of strings in `names`, return a
       dictionary containing the names as keys and child
       `ELEMENT_NODEs`, that have a `tagName` equal to the name.

       """
       nodes = [(x.tagName, x) for x in node.childNodes \
                if x.nodeType == x.ELEMENT_NODE and \
                x.tagName in names]
       return dict(nodes)

   def generic_parse(self, node, pad=0):
       """A Generic parser for arbitrary tags in a node.

       Parameters:

        - node:  A node in the DOM.
        - pad: `int` (default: 0)

          If 0 the node data is not padded with newlines.  If 1 it
          appends a newline after parsing the childNodes.  If 2 it
          pads before and after the nodes are processed.  Defaults to
          0.

       """
       npiece = 0
       if pad:
           npiece = len(self.pieces)
           if pad == 2:
               self.add_text('\n'+' '*self.indent)
       for n in node.childNodes:
           self.parse(n)
       if pad:
           if len(self.pieces) > npiece:
               self.add_text('\n'+' '*self.indent)

   def space_parse(self, node):
       self.add_text(' ')
       self.generic_parse(node)

   do_ref = space_parse
   do_emphasis = space_parse
   do_bold = space_parse
   do_computeroutput = space_parse
   do_formula = space_parse

   def do_compoundname(self, node):
       self.add_text('\n\n')
       data = node.firstChild.data
       self.add_text('%%feature("docstring") %s "\n'%data)

   def do_compounddef(self, node):
       kind = node.attributes['kind'].value
       if kind in ('class', 'struct'):
           prot = node.attributes['prot'].value
           if prot <> 'public':
               return
           names = ('compoundname', 'briefdescription',
                    'detaileddescription', 'includes')
           first = self.get_specific_nodes(node, names)
           for n in names:
               if first.has_key(n):
                   self.parse(first[n])
           self.add_text(['";','\n'])
           for n in node.childNodes:
               if n not in first.values():
                   self.parse(n)
       elif kind in ('file', 'namespace'):
           nodes = node.getElementsByTagName('sectiondef')
           for n in nodes:
               self.parse(n)

  # def do_includes(self, node):
  #     # self.add_text('C++ includes: ')
  #     self.generic_parse(node, pad=1)

   def do_parameterlist(self, node):
       # self.add_text(['\n', '\n', 'Parameters:', '\n'])
       if( node.hasChildNodes):
          self.add_text(['\nArguments:\n'])
          self.indent = 0
          self.generic_parse(node, pad=1)
          self.curCol = 2
          self.indent = 2          

   def do_para(self, node):
       # self.add_text('\n')
       # txt = ''
       #for n in node.parentNode.childNodes:
       #  txt = txt + n.data
       self.generic_parse(node, pad=1)
       self.curCol = self.indent
       # self.add_text( self.wrap_text(txt, self.curCol - self.indent))

   def do_parametername(self, node):
       self.indent = 2
       self.curCol = 2
       self.add_text('\n  ')
       self.add_text(self.wrap_text("%s:" % node.firstChild.data.strip(), 0))
       self.add_text('  ')
     
   def do_parameterdescription(self, node):
       self.indent = 6
       self.generic_parse(node, pad=0)
       self.indent = 2

   def do_detaileddescription(self, node):
       self.add_text('\nDetails:\n  ')
       self.indent = 2
       self.curCol = 2
       self.generic_parse(node, pad=1)
       self.add_text('\n\n')

   def do_definition(self, node):
       # self.indent = 2
       # self.curCol = 2
       # self.add_text('\nDescription:\n  ')
       # self.generic_parse(node, pad=1)
       # self.add_text('\n\n')
       #
       data = node.firstChild.data
       data = data.replace('"',r'\"')
       data = data.replace('\\',r'\\\\"')
       self.add_text('%s "\n%s'%(data, data))

   def do_briefdescription(self, node):
       self.indent = 2
       self.curCol = 2
       self.add_text('\nDescription:\n  ')
       self.generic_parse(node, pad=1)
       self.add_text('\n\n')

   def do_memberdef(self, node):
       prot = node.attributes['prot'].value
       id = node.attributes['id'].value
       kind = node.attributes['kind'].value
       tmp = node.parentNode.parentNode.parentNode
       compdef = tmp.getElementsByTagName('compounddef')[0]
       cdef_kind = compdef.attributes['kind'].value

       if prot == 'public':
           first = self.get_specific_nodes(node, ('definition', 'name'))
           name = first['name'].firstChild.data
           if name[:8] == 'operator': # Don't handle operators yet.
               return

           # defn = first['definition'].firstChild.data
           defn = ''
           for n in first['definition'].childNodes:
             defn = defn + n.data

           self.add_text('\n')
           # a new function
           self.add_text('%feature("docstring") ')

           # its type determined by parents
           anc = node.parentNode.parentNode
           if cdef_kind in ('file', 'namespace'):
               ns_node = anc.getElementsByTagName('innernamespace')

               if not ns_node and cdef_kind == 'namespace':
                   ns_node = anc.getElementsByTagName('compoundname')
               if ns_node:
                   ns = ns_node[0].firstChild.data
                   self.add_text(' %s::%s " \nUsage:\n  ' %(ns, name) )
                   self.indent = 2
                   self.curCol = 2
                   self.add_text( self.wrap_text( defn.split(' ')[-1], 0) )
               #   self.add_text(' %s::%s "\n%s'%(ns, name, defn.split(' ')[-1]))
               else:
                   self.add_text(' %s " \nUsage:\n  ' % name )
                   self.indent = 2
                   self.curCol = 2
                   self.add_text( self.wrap_text( defn.split(' ')[-1], 0) )
               #   self.add_text(' %s "\n%s'%(name, defn))
           elif cdef_kind in ('class', 'struct'):
               # Get the full function name.
               anc_node = anc.getElementsByTagName('compoundname')
               cname = anc_node[0].firstChild.data
               self.add_text(' %s::%s " \nUsage:\n  ' %(cname, name) )
               self.indent = 2
               self.curCol = 2
               defName = defn.split(' ')[-1]
               # force the first character lower case
               if( defName == cname.split(':')[-1] ): # constructor
                 if (len(defName) > 1 ):
                   defName = defName[0].lower() + defName[1:]
                 self.add_text( self.wrap_text( defName, 0) )
               else: 
                 self.add_text( 'x.' )
                 self.add_text( self.wrap_text( defName, 0) )

           for n in node.childNodes:
               if n not in first.values():
                   self.parse(n)
           self.add_text(['";', '\n'])

   def do_sectiondef(self, node):
       kind = node.attributes['kind'].value
       if kind in ('public-func', 'func', 'user-defined'):
           self.generic_parse(node)

   def do_xrefsect(self, node):
       # first child
       # print node.firstChild.tagName
       if(node.firstChild.firstChild.data == 'Test'):
           self.add_text(['\nExamples:\n'])
           # get file name
           filename = node.firstChild.nextSibling.firstChild.firstChild.data
           # get content as string
           try:
              try:
                file = open(filename.strip())
              except:
                file = open('../doc/' + filename.strip() )
              cont = file.read()
              file.close()
              cont = "  " + ("\n  ".join( cont.split('\n') ))
              self.add_text(cont)
           except:
              print "File " + filename + " does not exist\n"
              sys.exit(1)
              self.add_text("  " + filename + "does not exist\n")
       self.add_text("\nDetails:\n  ");
       self.indent = 2;
       self.curCol = 2;

   def do_simplesect(self, node):
       kind = node.attributes['kind'].value
       if kind in ('date', 'rcs', 'version'):
           pass
       elif kind == 'warning':
           self.add_text('\nWARNING:\n  ')
           self.indent = 2
           self.curCol = 2
           self.generic_parse(node)
           self.add_text(['\n  \n'])
           self.indent = 2
           self.curCol = 2
       elif kind == 'see':
           self.add_text('\nSee Also:\n  ')
           self.indent = 2
           self.curCol = 2
           self.generic_parse(node)
           self.indent = 2
           self.curCol = 2
       elif kind == 'note':
           self.add_text('\nNote:\n  ')
           self.indent = 2
           self.curCol = 2
           self.generic_parse(node)
           self.indent = 2
           self.curCol = 2
       elif kind == 'return':
           self.add_text('\nValue:\n  ')
           self.indent = 2
           self.curCol = 2
           self.generic_parse(node)
           self.indent = 2
           self.curCol = 2
       else:
           self.generic_parse(node)
       self.add_text("\nDetails:\n  ");
       self.indent = 2;
       self.curCol = 2;

   def do_argsstring(self, node):
       txt = ''
       for n in node.childNodes:
         txt = txt + n.data
       
       # print  'Ori: ', txt, '\n'
       # replace the trailing const
       txt = txt.replace(')  const',')')
       txt = txt.replace(') const',')')
       txt = txt.replace(')const',')')
       args = txt.split(',')
       out=[]
       for s in args:
         piece = s.split('=')
         var = piece[0].split(' ')[-1].split(')')[0].split('(')[-1]
         if( len( piece ) == 2 ): 
           defVal = piece[1].split('(')[0].split(')')[0].split(')')[0]
           defVal = defVal.replace('vectorlu','[]')
           defVal = defVal.replace('vectoru','[]')
           defVal = defVal.replace('vectorl','[]')
           defVal = defVal.replace('vectorf','[]')
           defVal = defVal.replace('vectora','[]')
           defVal = defVal.replace('vectorop','[]')
           defVal = defVal.replace('vectorstr','[]')
           defVal = defVal.replace('dictionary','{}')
           defVal = defVal.replace('matrix','[]')
           defVal = defVal.replace('true','True')
           defVal = defVal.replace('false','False')
           defVal = defVal.replace('\\', r'\\\\')
           defVal = defVal.replace('"', r'\"')
           out.append( var + '=' + defVal  )
         else:
           out.append( var )

       self.indent = 4
       # print  '(' + (', '.join(out)) + ')\n'
       self.add_text( self.wrap_text( '(' + (', '.join(out)) + ')\n', self.curCol - self.indent ) + '\n')
       self.add_text( 'OriArgString: ' + txt + '\n')
       
       # remove type and default value.
       # self.add_text(txt)
       # self.generic_parse(node, pad=1)

   def do_member(self, node):
       kind = node.attributes['kind'].value
       refid = node.attributes['refid'].value
       if kind == 'function' and refid[:9] == 'namespace':
           self.generic_parse(node)

   def do_doxygenindex(self, node):
       self.multi = 1
       comps = node.getElementsByTagName('compound')
       for c in comps:
           refid = c.attributes['refid'].value
           fname = refid + '.xml'
           if not os.path.exists(fname):
               fname = os.path.join(self.my_dir,  fname)
           print "parsing file: %s"%fname
           p = Doxy2SWIG(fname)
           p.generate()
           # self.pieces.extend(self.clean_pieces(p.pieces))
           self.pieces.extend(p.pieces)

   def write(self, fname):
       o = my_open_write(fname)
       if self.multi:
           o.write("".join(self.pieces))
       else:
           # o.write("".join(self.clean_pieces(self.pieces)))
           o.write("".join(self.pieces))
       o.close()

   def wrap_text(self, text, start_pos ):
       """ wrap text given current indent """ 
       strs = textwrap.wrap(text.lstrip('\n '), width=self.maxChar, 
         initial_indent=' '*(start_pos+self.indent),
         subsequent_indent = ' '*self.indent)
       if( len(strs) == 1 ):
         self.curCol = len(strs[-1]) + start_pos + self.indent
       elif( len(strs) > 1 ):
         self.curCol = len(strs[-1])
       else:
         self.curCol = 0
       if( self.curCol < self.indent):
         self.curCol = self.indent
       # print 'wrap \'' + text + '\' to \'' + ('\n'.join(strs)).lstrip() + '\''
       return  ('\n'.join(strs)).lstrip()

def main(input, output):
   p = Doxy2SWIG(input)
   p.generate()
   p.write(output)

if __name__ == '__main__':
   if len(sys.argv) != 3:
       print __doc__
       sys.exit(1)
   main(sys.argv[1], sys.argv[2])
