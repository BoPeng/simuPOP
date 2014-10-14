#!/usr/bin/env python
"""Doxygen XML to SWIG docstring converter.

Converts Doxygen generated XML files into a file containing docstrings
that can be used by SWIG-1.3.x. It also generate a latex file with
reference manauls for all classes and functions. These definitions
can be included into the simuPOP reference manual easily.

Usage:

 doxy2swig.py input.xml output.i output.tex

input.xml is your doxygen generated XML file; output.i and output.tex
are where the output interface and latex files will be written.

"""

# This code is implemented using Mark Pilgrim's code as a guideline:
#     http://www.faqs.org/docs/diveintopython/kgp_divein.html
#
# Author: Prabhu Ramachandran
# Modified by: Bo Peng
# Last Modified, Apr, 2007
# License: BSD style

# Reference links (may be obsolete soon):
#
#    XML DOM in general: http://www.w3schools.com/dom/default.asp
#        Note that Element and Text class names. This file does
#        not handle Attribute class.
#    Python minidom: http://docs.python.org/lib/module-xml.dom.minidom.html
#


from xml.dom import minidom
import re, textwrap, sys, types, os.path, sets, inspect
from pydoc import *

from docutils import core
from docutils import writers
from docutils.writers.latex2e import Writer
from docutils.writers.latex2e import LaTeXTranslator

overrides = {'input_encoding': 'ascii',
             'output_encoding': 'latin-1'}

import pprint
try:
    import simuPOP as sim
except:
    # this script needs simuPOP but I will have to let it run without simuPOP from time to time
    pass

# customized latex output
class myLaTeXTranslator(LaTeXTranslator):
    def visit_definition_list(self, node):
        self.body.append( '\n{\leftskip 0.3in \parindent=-0.3in ' )

    def depart_definition_list(self, node):
        self.body.append( '\\par}\n' )

    def visit_term(self, node):
        self.body.append('\\emph{')

    def depart_term(self, node):
        # definition list term.
        # \leavevmode results in a line break if the term is followed by a item list.
        self.body.append(': } ')

    #def visit_document(self, node):
    #    self.body_prefix.append('\\begin{document}\n')
    #    # REMOVE THIS FROM THE INITIAL IMPLEMENTAION
    #    #self.body.append('\n\\setlength{\\locallinewidth}{\\linewidth}\n')

    def visit_enumerated_list(self, node):
        # STOP USING SELF-DEFINED ENUMERATION LIST
        self.body.append('\\begin{enumerate}\n')

    def depart_enumerated_list(self, node):
        self.body.append('\\end{enumerate}\n')

    def literal_block_env(self, begin_or_end):
        # THIS WILL BE WRONG BECAUSE THIS CAN NOT
        # EXTEND OVER SEVERAL LINES.
        if begin_or_end == 'begin':
            return '\\lstinline!'
        return '!\n'

    def visit_doctest_block(self, node):
        self.literal_block = 1
        self.insert_none_breaking_blanks = 1
        self.body.append('{\\ttfamily \\raggedright \\noindent\n')
        #self.body.append( '\\begin{verbatim}' )
        #self.verbatim = 1

    def depart_doctest_block(self, node):
        #self.body.append( '\\end{verbatim}\n' )
        #self.verbatim = 0
        self.body.append('\n}')
        self.insert_none_breaking_blanks = 0
        self.literal_block = 0


class myWriter(Writer):
    def __init__(self):
        Writer.__init__(self)
        self.translator_class = myLaTeXTranslator



class Doxy2SWIG:
    """Converts Doxygen generated XML files into a data struture
    (self.content) that can be written in doctring and latex format
    Docstrings that can be used by SWIG-1.3.x that have support for
    feature("docstring").
    """

    def __init__(self, src):
        """Initialize the instance given a source filename
        """
        f = open(src)
        self.my_dir = os.path.dirname(src)
        # get everything into this XMLDOC
        self.xmldoc = minidom.parse(f).documentElement
        f.close()

        # all the information will be in content,
        # interface and tex files will be generated from content
        self.content = []

        # ignore these tags
        self.ignores = ('inheritancegraph', 'param', 'listofallmembers',
                        'innerclass', 'name', 'declname', 'incdepgraph',
                        'invincdepgraph', 'programlisting', 'type',
                        'references', 'referencedby', 'location',
                        'collaborationgraph', 'reimplements',
                        'reimplementedby', 'derivedcompoundref',
                        'basecompoundref', 'header', 'includes')

        # match one or more space/tab etc
        self.space_re = re.compile(r'\s+')

        self.lead_spc = re.compile(r'^(%feature\S+\s+\S+\s*?)"\s+(\S)')
        self.maxChar = 70
        # current field in self.content
        self.curField = ''
        self.uniqueName = []


    def generate(self):
        """Parses the file set in the initialization. The resulting
        data is stored in `self.content`.
        """
        self.parse(self.xmldoc)


    def parse(self, node):
        """Parse a given node. This function in turn calls the
        `parse_<nodeType>` functions which handle the respective
        nodes.
        """
        # every node in a XML tree has __class__.__name__ e.g. Text, Document, etc.
        # This function get this information and call the relevant function to parse this node.
        # This is a dispatch function that is by nature recursive.
        #
        # in fact, we will only see Element and Text classes.
        handlerMethod = getattr(self, "parse_%s" % node.__class__.__name__)
        handlerMethod(node)


    def parse_Element(self, node):
        """Parse an `ELEMENT_NODE`. This calls specific
        `do_<tagName>` handers for different elements. If no handler
        is available the `parse_childnodes` method is called. All
        tagNames specified in `self.ignores` are simply ignored.
        """
        name = node.tagName
        if name in self.ignores:
            return
        # a bunch of 'tagName's will be handled by do_XXX functions
        # We defined do_ref, do_emphasis, do_bold, do_computeroutput
        # do_formula (all using space_parse), do_compindname, do_compounddef
        # do_parameterlist, do_para, do_parametername, do_parameterdescription
        # etc.....
        #
        attr = "do_%s" % name
        if hasattr(self, attr):
            handlerMethod = getattr(self, attr)
            handlerMethod(node)
        else:
            self.parse_childnodes(node)

    def parse_Text(self, node):
        txt = node.data
        # ignore pure whitespace
        m = self.space_re.match(txt)
        # for example '     ' yields 5 blank matches
        if m is not None and len(m.group()) == len(txt):
            pass
        else:
            self.add_text(txt)

    def add_text(self, value):
        """Adds text corresponding to `value` into `self.content`."""
        # each argument is a list
        if self.curField == 'Arguments':
            if type(value) in (types.ListType, types.TupleType):
                self.content[-1][self.curField][-1]['Description'] += ' '.join(value)
            else:
                self.content[-1][self.curField][-1]['Description'] += value
        else:
            if type(value) in (types.ListType, types.TupleType):
                self.content[-1][self.curField] += ' '.join(value)
            else:
                self.content[-1][self.curField] += value

    def get_specific_nodes(self, node, names):
        """Given a node and a sequence of strings in `names`, return a
        dictionary containing the names as keys and child
        `ELEMENT_NODEs`, that have a `tagName` equal to the name.

        """
        nodes = [(x.tagName, x) for x in node.childNodes \
                         if x.nodeType == x.ELEMENT_NODE and \
                         x.tagName in names]
        return dict(nodes)


    def parse_childnodes(self, node):
        """parse all child nodes.
        """
        #self.parse(node)
        for n in node.childNodes:
            self.parse(n)


    def space_parse(self, node):
        self.add_text(' ')
        self.parse_childnodes(node)

    do_ref = space_parse
    do_emphasis = space_parse
    do_bold = space_parse
    do_computeroutput = space_parse
    do_formula = space_parse


    def do_compoundname(self, node):
        data = node.firstChild.data
        self.content.append({'Name': data, 'type': 'class'})


    def do_emphasis(self, node):
        #if node.firstChild is not None:
        #    self.add_text(r'<em>%s</em>' % node.firstChild.data)
        self.add_text('<em>')
        self.parse_childnodes(node)
        self.add_text('</em>')


    def do_computeroutput(self, node):
        self.add_text('<tt>')
        self.parse_childnodes(node)
        self.add_text('</tt>')


    def do_bold(self, node):
        self.add_text('<bf>')
        self.parse_childnodes(node)
        self.add_text('</bf>')

    def do_linebreak(self, node):
        #newline tag indicates the linebreak we do want
        self.add_text(r'</newline>')


    def do_itemizedlist(self, node):
        if( node.hasChildNodes ):
            self.add_text(r'<itemize>')
            self.parse_childnodes(node)
            self.add_text(r'</itemize>')


    def do_listitem(self, node):
        self.add_text('<item>')
        self.parse_childnodes(node)
        self.add_text('</item>')


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
            for n in node.childNodes:
                if n not in first.values():
                    self.parse(n)
        elif kind in ('file', 'namespace'):
            nodes = node.getElementsByTagName('sectiondef')
            for n in nodes:
                self.parse(n)

    # def do_includes(self, node):
    #         # self.add_text('C++ includes: ')
    #         self.parse_childnodes(node, pad=1)

    def do_parameterlist(self, node):
        if( node.hasChildNodes):
             self.curField = 'Arguments'
             if not self.content[-1].has_key('Arguments'):
                 self.content[-1]['Arguments'] = []
             self.parse_childnodes(node)
        self.curField = 'Details'


    def do_para(self, node):
        if self.curField == 'Details' and self.content[-1]['Details'] != '':
            self.content[-1]['Details'] += '\n\n'
        self.parse_childnodes(node)


    def do_parametername(self, node):
        if node.firstChild.nodeName == 'ref':
            parameter_name = node.firstChild.firstChild.data.strip()
        else:
            parameter_name = node.firstChild.data.strip()
        assert self.curField == 'Arguments'
        self.content[-1][self.curField].append({'Name': parameter_name, 'Description': ''})


    def do_detaileddescription(self, node):
        self.curField = 'Details'
        if not self.content[-1].has_key('Details'):
            self.content[-1]['Details'] = ''
        self.parse_childnodes(node)


    def do_briefdescription(self, node):
        self.curField = 'Description'
        self.content[-1]['Description'] = ''
        self.parse_childnodes(node)


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

            # its type determined by parents
            anc = node.parentNode.parentNode
            if cdef_kind in ('file', 'namespace'):
                ns_node = anc.getElementsByTagName('innernamespace')

                if not ns_node and cdef_kind == 'namespace':
                    ns_node = anc.getElementsByTagName('compoundname')
                if ns_node:
                    ns = ns_node[0].firstChild.data
                    #
                    # FIXME: when the sandbox namespace is introduced to sandbox.h/cpp, it somehow
                    # manifests itself to all namespaces in this script (doxygen output is correct).
                    # I do not know how to thoroughly fix this problem, and have to manually remove
                    # sandbox from the namespace. Hopefully a full solution could be found later.
                    #
                    if ':' in ns:
                        ns = ns.split(':')[0]
                    func_name = '%s::%s' %(ns, name)
                    self.content.append({'Name': func_name, 'type': 'global_function'})
                    self.content[-1]['Usage'] = func_name.split(':')[-1]
                    self.curField = 'Usage'
                else:
                    self.content.append({'Name': name})
                    print "Content type unknown??? ", name
                    self.content[-1]['Usage'] = ''
                    self.curField = 'Usage'
                    self.add_text( defn.split(' ')[-1] )
            elif cdef_kind in ('class', 'struct'):
                # Get the full function name.
                anc_node = anc.getElementsByTagName('compoundname')
                cname = anc_node[0].firstChild.data
                classname = cname.split('::')[-1]
                if classname == name:
                    self.content.append({'Name': '%s::%s' %(cname, name),
                                         'type': 'constructorofclass_%s' % cname})
                else:
                    self.content.append({'Name': '%s::%s' %(cname, name),
                                         'type': 'memberofclass_%s' % cname})
                self.content[-1]['Usage'] = ''
                self.curField = 'Usage'
                defName = defn.split(' ')[-1]
                if( defName == cname.split(':')[-1] ): # constructor
                    self.add_text( self.swig_text( defName, 0, 0 ) )
                else:
                    self.add_text( 'x.' )
                    self.add_text( self.swig_text( defName, 0, 0 ) )

            for n in node.childNodes:
                if n not in first.values():
                    self.parse(n)


    def do_sectiondef(self, node):
        kind = node.attributes['kind'].value
        if kind in ('public-func', 'func', 'user-defined'):
                self.parse_childnodes(node)


    def do_xrefsect(self, node):
        # first child
        # print node.firstChild.tagName
        if(node.firstChild.firstChild.data == 'Test'):
            self.curField = 'Examples'
            self.content[-1]['Examples'] = ''
            self.parse_childnodes(node)


    def do_simplesect(self, node):
        kind = node.attributes['kind'].value
        if kind in ('date', 'rcs', 'version'):
            pass
        if kind in ('warning', 'see', 'note', 'return'):
            self.curField = kind
            self.content[-1][kind] = ''
            self.parse_childnodes(node)
        else:
            self.parse_childnodes(node)
        self.curField = 'Details'


    def do_argsstring(self, node):
        txt = ''
        for n in node.childNodes:
            txt = txt + n.data
        if txt.split('=')[-1].strip() == '0':
             # remove =0 from function definition (WHY?)
             txt = '='.join(txt.split('=')[:-1])
        self.content[-1]['cppArgs'] = txt
        # @ is used instead of , to avoid separation of replaced text, it will be replaced back to ,
        txt = re.sub(r'uintList\(NULL\)', r'UNSPECIFIED', txt)
        txt = re.sub(r'lociList\(NULL\)', r'UNSPECIFIED', txt)
        txt = re.sub(r'stringList\((["\w]+),\s*(["\w]+)\)', r'[\1@ \2]', txt)
        txt = re.sub(r'vectorstr\(1, ([^)]+)\)', r'[\1]', txt)
        txt = re.sub(r'(.*)vector(str|u|i)\(1,\s*([\w"\d]+)\)(.*)', r'\1\3\4', txt)
        txt = re.sub(r'\)\s*const\s*$', ')', txt)
        #txt = txt.replace(')    const',')')
        args = txt.split(',')
        out=[]
        for s in args:
            #use @ in order to split the arguments correctly
            s = s.replace('@', r',')
            piece = s.split('=')
            var = piece[0].split(' ')[-1].split(')')[0].split('(')[-1]
            #delete &,*
            var = var.replace('&', '')
            var = var.replace('*', '')
            if( len( piece ) == 2 ):
                defVal = piece[1].split('(')[0].split(')')[0].split(')')[0]
                # re function used to repalce the following sentances
                vect = re.compile('vector(info|u|i|f|a|op|splitter|str|vsp)')
                defVal = vect.sub('[]', defVal)
                defVal = defVal.replace('strDict','{}')
                defVal = defVal.replace('intDict','{}')
                defVal = defVal.replace('intMatrix','[]')
                defVal = defVal.replace('floatMatrix','[]')
                defVal = defVal.replace('vspID','[]')
                defVal = defVal.replace('true','True')
                defVal = defVal.replace('false','False')
                defVal = defVal.replace('NULL','None')
                defVal = defVal.replace('uintListFunc','[]')
                defVal = defVal.replace('uintList','ALL_AVAIL')
                defVal = defVal.replace('lociList','ALL_AVAIL')
                defVal = defVal.replace('intList','ALL_AVAIL')
                defVal = defVal.replace('floatListFunc','[]')
                defVal = defVal.replace('floatList','[]')
                defVal = defVal.replace('stringList','ALL_AVAIL')
                defVal = defVal.replace('stringMatrix','[]')
                defVal = defVal.replace('subPopList', 'ALL_AVAIL')
                defVal = defVal.replace('opList','[]')
                defVal = defVal.replace('string','""')
                defVal = defVal.replace('UnnamedSubPop','""')
                #defVal = defVal.replace('""', "''")
                out.append(var + '=' + defVal)
            else:
                out.append(var)
        self.add_text('(' + ', '.join(out) + ')')


    def do_member(self, node):
        kind = node.attributes['kind'].value
        refid = node.attributes['refid'].value
        if kind == 'function' and refid[:9] == 'namespace':
            self.parse_childnodes(node)


    def do_doxygenindex(self, node):
        '''Parse files included as compound, member and refid, like

            <doxygenindex xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="index.xsd" version="1.3.9.1">
              <compound refid="a00153" kind="class"><name>simuPOP::affectedSibpairSample</name>
              <member refid="a00153_1a0" kind="function"><name>affectedSibpairSample</name></member>
              </compound>
            </doxygenindex>

            In this case, we process each compound using a Doxy2SWIG class and
            obtain its results.
        '''
        comps = node.getElementsByTagName('compound')
        for c in comps:
            refid = c.attributes['refid'].value
            fname = refid + '.xml'
            if not os.path.exists(fname):
                fname = os.path.join(self.my_dir, fname)
            print "parsing file: %s" % fname
            #try:
            p = Doxy2SWIG(fname)
            p.generate()
            self.content.extend(p.content)
            #except Exception, e:
            #    print "This file can not be parsed, something wrong with file format"
            #    print e
            #    pass

    def fix_entries(self):
        def myhash(entry):
            'encode an entry to a string for easy comparison'
            ret = entry['Name']
            if entry.has_key('Description'):
                ret += entry['Description']
            if entry.has_key('Details'):
                ret += entry['Details']
            return ret
        seen = []
        self.content = [x for x in self.content if not (myhash(x) in seen or seen.append(myhash(x)))]
        print "Unique entries: ", len(self.content)
        #
        for entry in self.content:
            # latex names cannot have _ or number (for all_ind, getInfo2 etc)
            entry['Name'] = entry['Name'].replace('_', '').replace('2','two').replace('3', 'three').replace('0', 'zero').replace('9', 'nine')
            if (entry.has_key('Details') and '<group>' in entry['Details']):
                piece1 = entry['Details'].split('<group>')
                try:
                    piece2 = piece1[1].split('</group>')
                    entry['Details'] = piece1[0] + piece2[1]
                except:
                    print 'WRONG GROUP INFORMATION: ', entry['Name'], entry['Details']
                entry['group'] = piece2[0].strip()
            elif (entry.has_key('Description') and '<group>' in entry['Description']):
                piece1 = entry['Description'].split('<group>')
                piece2 = piece1[1].split('</group>')
                entry['Description'] = piece1[0] + piece2[1]
                entry['group'] = piece2[0].strip()
            else:
                entry['group'] = ''
            #
            if not entry.has_key('Doc'):
                entry['Doc'] = ''
                if entry.has_key('Description') and entry['Description'] != '':
                    entry['Doc'] += entry['Description']
                if entry.has_key('Details') and entry['Details'] != '':
                    entry['Doc'] += entry['Details']
            #
            if entry['Doc'] == '':
                print('Warning: No documentation for ', entry['Name'])
                entry['Doc'] = 'FIXME: No document'
            #
            entry['ignore'] = 'CPPONLY' in entry['Doc']
            entry['hidden'] = 'HIDDEN' in entry['Doc']
            #
            if entry['ignore'] and '~' in entry['Name']:
                print "Desctructor of %s has CPPONLY. Please correct it." % entry['Name']
                sys.exit(1)
        # handle Examples
        for entry in self.content:
            if entry.has_key('Examples') and entry['Examples'].strip() != '':
                # get file name
                filename = entry['Examples'][4:]
                title = ''
                # if the filename has an optional title, separate them
                if ' ' in filename:
                    pieces = filename.split(' ')
                    filename = pieces[0]
                    title = ' '.join(pieces[1:])
                # get content as string
                try:
                    try:
                        # usual ../doc/log directory
                        file = open(os.path.join('..', 'doc', 'log', filename))
                        entry['ExampleFile'] = os.path.join('..', 'doc', 'log', filename)
                    except:
                        # local file
                        file = open(filename)
                        entry['ExampleFile'] = filename
                    cont = file.read()
                    file.close()
                    # add file content to 'Examples'
                    entry['Examples'] = cont
                    entry['ExampleTitle'] = title
                except:
                    entry['ExampleFile'] = None
                    entry['ExampleTitle'] = title
                    print "File " + filename + " does not exist\n"

    def post_process(self):
        # remove duplicate entry
        # They might be introduced if a function is list both under 'file' and under 'namespace'
        # Add additional entries manually
        self.content.extend([
            {'Name': u'simuPOP::Population::dvars',
             'type': u'memberofclass_simuPOP::Population',
             'Description': '',
             'Details': ur'<group>9-var1</group>' \
                'Return a wrapper of Python dictionary returned by <tt>vars(subPop)</tt> ' \
                'so that dictionary keys can be accessed as attributes.',
             'cppArgs': u'(vspID subPop=[])',
             'Usage': u'x.dvars(subPop=[])',
             },
            {'Name': u'simuPOP::Simulator::dvars',
             'type': u'memberofclass_simuPOP::Simulator',
             'Description': '',
             'Details': ur'<group>9-var1</group>' \
                'Return a wrapper of Python dictionary returned by <tt>vars(rep, subPop)</tt> ' \
                'so that dictionary keys can be accessed as attributes.',
             'cppArgs': u'(int rep, vspID subPop=[])',
             'Usage': u'x.dvars(rep, subPop=[])',
            }
        ])
        print "Number of entries: ", len(self.content)
        self.fix_entries()
        #
        self.content.sort(lambda x, y: x['Name'] > y['Name'])

    def write_swig(self, out):
        for entry in self.content:
            if entry['ignore']:
                if entry.has_key('cppArgs'):
                    print >> out, '%%ignore %s%s;\n' % (entry['Name'], entry['cppArgs'].strip())
                else:
                    print >> out, '%%ignore %s;\n' % entry['Name'].strip()
                continue
            if entry['hidden']:
                print >> out, '%%feature("docstring") %s "Obsolete or undocumented function."\n' % entry['Name']
                continue
            print >> out, '%%feature("docstring") %s "\n' % entry['Name']
            if entry.has_key('Description') and entry['Description'] != '':
                print >> out, 'Description:'
                print >> out, '\n    %s\n' % self.swig_text(entry['Description'], 0, 4)
            if entry.has_key('Usage') and entry['Usage'] != '':
                print >> out, 'Usage:'
                print >> out, '\n    %s\n' % self.swig_text(entry['Usage'], 0, 6)
            if entry.has_key('Details') and entry['Details'] != '':
                print >> out, 'Details:'
                print >> out, '\n    %s\n' % self.swig_text(entry['Details'], 0, 4)
            if entry.has_key('Arguments') and entry['Arguments'] != '':
                print >> out, 'Arguments:\n'
                for arg in entry['Arguments']:
                    print >> out, '    %-16s%s' % (arg['Name']+':', self.swig_text(arg['Description'], 0, 20))
                print >> out
            if entry.has_key('note') and entry['note'] != '':
                print >> out, 'Note:'
                print >> out, '\n    %s\n' % self.swig_text(entry['note'], 0, 4)
            if entry.has_key('Examples') and entry['Examples'] != '':
                print >> out, 'Example:'
                print >> out, '\n%s\n' % entry['Examples'].replace('\\', r'\\\\').replace('"', r'\"')
            #if len(str) > 2048:
            #    #print 'Entry %s is too long (%d)' % (entry['Name'], len(str))
            #    str = str[0:1970] + '...' + '\n\nPlease refer to the reference manual for more details.\n\n'
            print >> out, '\"; \n'


    def replacePythonConstants(self, text):
        '''replace True by ALL_AVAIL etc'''
        try:
            text = text.replace('subPops=%s' % sim.ALL_AVAIL, 'subPops=ALL_AVAIL')
            for sexMode in ['NO_SEX', 'RANDOM_SEX', 'PROB_OF_MALES', 'NUM_OF_MALES', 'NUM_OF_FEMALES']:
                text = text.replace('sexMode=%s' % eval('sim.' + sexMode), 'sexMode=%s' % sexMode)
            # for polySex and alphaSex
            for sex in ['MALE', 'FEMALE']:
                text = text.replace('Sex=%s' % eval('sim.' + sex), 'Sex=%s' % sex)
            for transmitter in [ 'HaplodiploidGenoTransmitter', 'MendelianGenoTransmitter', 'MitochondrialGenoTransmitter',
                'SelfingGenoTransmitter', 'CloneGenoTransmitter', 'Recombinator']:
                text = text.replace('ops=%s' % eval('sim.' + transmitter+'()'), 'ops=%s' % (transmitter+'()'))
            for ms in ['MatingScheme']:
                text = text.replace('matingScheme=%s' % eval('sim.' + ms+'()'),
                    'matingScheme=%s' % (ms+'()'))
        except Exception, e:
            print e
            pass
        return text

    def scan_module(self, module):
        ''' scan python module file and retrieve function definitions '''
        # add module entry
        # load module
        try:
            # resolve is defined in pydoc
            object, name = resolve(module, True)
            print 'Scanning module %s' % module
            #exec('import ' + module)
        except Exception, e:
            print "Module %s failed to load. It is description is not documented." % module
            print "Please compile simuPOP and rerun this program"
            print
            print e
            return
        #object = eval(module)
        synop, desc = splitdoc(getdoc(object))
        self.content.append({'type': 'docofmodule_' + module, 
            'Name': module, 'ignore': False, 'hidden': False, 'module': module})
        self.content[-1]['Description'] = desc
        self.content[-1]['Doc'] = desc
        for key, value in inspect.getmembers(object, inspect.isclass):
            if (inspect.getmodule(value) or object) is object:
                if not visiblename(key):
                    continue
                self.content.append({'type': 'module_class', 'module': module})
                self.content[-1]['Name'] = key
                self.content[-1]['Doc'] = getdoc(value)
                self.content[-1]['InitDoc'] = getdoc(value.__init__)
                self.content[-1]['ignore'] = 'CPPONLY' in self.content[-1]['Doc']
                self.content[-1]['hidden'] = 'HIDDEN' in self.content[-1]['Doc']
                # these is no details...
                try:
                    args, varargs, varkw, defaults = inspect.getargspec(value.__init__)
                except:
                    print('Failed to get args for {}'.format(value.__init__))
                    self.content.pop()
                    continue
                self.content[-1]['Usage'] = key + inspect.formatargspec(
                    args, varargs, varkw, defaults)
                self.content[-1]['Usage'] = self.content[-1]['Usage'].replace('self, ', '').replace('self)', ')')
                self.content[-1]['Usage'] = self.replacePythonConstants(self.content[-1]['Usage'])
                members = []
                for key1, value1 in inspect.getmembers(value, inspect.isroutine):
                    if not visiblename(key1) or key1.startswith('_'):
                        continue
                    try:
                        args, varargs, varkw, defaults = inspect.getargspec(value1)
                    except:
                        continue
                    member = {'Name': key1}
                    member['Usage'] = key1 + inspect.formatargspec(
                        args, varargs, varkw, defaults)
                    member['Usage'] = member['Usage'].replace('self, ', '').replace('self)', ')')
                    member['Usage'] = re.sub(r'(<simuPOP.(\w+)>)', r'\2()', member['Usage'])
                    member['Doc'] = getdoc(value1)
                    members.append(member)
                self.content[-1]['Members'] = members
        for key, value in inspect.getmembers(object, inspect.isroutine):
            if inspect.isbuiltin(value) or inspect.getmodule(value) is object:
                if not visiblename(key) or not inspect.isfunction(value):
                    continue
                self.content.append({'type': 'module_func', 'module': module})
                self.content[-1]['Name'] = key
                args, varargs, varkw, defaults = inspect.getargspec(value)
                self.content[-1]['Usage'] = key + inspect.formatargspec(
                    args, varargs, varkw, defaults)
                self.content[-1]['Usage'] = self.replacePythonConstants(self.content[-1]['Usage'])
                self.content[-1]['Description'] = getdoc(value)
                # latex can not yet handle numbers in function name
                self.content[-1]['ignore'] = 'CPPONLY' in self.content[-1]['Description'] or \
                    not self.content[-1]['Name'].isalpha()
                self.content[-1]['hidden'] = 'HIDDEN' in self.content[-1]['Description'] or \
                    not self.content[-1]['Name'].isalpha()
                # these is no details...
                self.content[-1]['Doc'] = self.content[-1]['Description']


    def latexName(self, name):
        # function name can overload
        uname = None
        for suffix in [''] + [chr(x) for x in range(ord('a'), ord('z'))]:
            if name + suffix not in self.uniqueName:
                uname = name + suffix
                self.uniqueName.append(uname)
                break
        if uname is None:
            print name, ' has too many overload names!'
            print self.uniqueName
            sys.exit(0)
        return uname.replace(':', '').replace('~', 'tld').replace('_', 'us').replace('0', 'o').replace('1', 'l').replace('2', 'z')


    def latex_text(self, text):
        """ wrap text given current indent """
        # handle the in-text and displayed math formulas
        newtext = ''
        inmath = False
        indispmath = False
        dispstart = False
        for i in range(len(text)):
            if dispstart:
                # skip [ or ] after \[, \] for display math mode
                dispstart = False
                continue
            if text[i] == '$':
                inmath = not inmath
                newtext += '$'
            elif i+1 < len(text) and text[i:i+2] == r'\[':
                indispmath = True
                newtext += '$$'
                dispstart = True
            elif i+1 < len(text) and text[i:i+2] == r'\]':
                indispmath = False
                newtext += '$$'
                dispstart = True
            else:
                if not inmath and not indispmath and text[i] in ['\\', '&', '$', '~', '%', '#', '_', '{', '}', '^']:
                    newtext += '\\' + text[i]
                else:
                    newtext += text[i]
        text = newtext
        #text = re.compile(r'%s' % text)
        text = text.replace('<em>', r'{\em ')
        text = text.replace('</em>', '}')
        text = text.replace('<tt>', r'{\tt ')
        text = text.replace('</tt>', '}')
        text = text.replace('<bf>', r'{\bf ')
        text = text.replace('</bf>', '}')
        text = text.replace('</newline>', '\n\n')
        text = text.replace('<itemize>', r'\begin{itemize} ')
        text = text.replace('</itemize>', r'\end{itemize}')
        text = text.replace('<item>', r'\item ')
        text = text.replace('</item>', ' ')
        text = text.replace('>','>{}')
        #text = text.replace('<par>', r'\par ')
        #text = text.replace('</par>', ' ')
        return text


    def latex_formatted_text(self, text):
        """format text according to some simple rules"""
        txt = core.publish_parts(source=text, writer = myWriter())['body']
        # some versions of textutil adds this automatically, which is a pain.
        txt = txt.replace(r'\setlength{\locallinewidth}{\linewidth}', '')
        if 'unexpected unindent' in txt:
            print text
            return ''
        else:
            return txt.lstrip()
                
    def swig_text(self, text, start_pos, indent):
        """ wrap text given current indent """
        #delete format sign used in latex file
        text = text.replace('<em>', '')
        text = text.replace('</em>', '')
        text = text.replace('<tt>', '')
        text = text.replace('</tt>', '')
        text = text.replace('<bf>', '')
        text = text.replace('</bf>', '')
        text = text.replace('<itemize>', '')
        text = text.replace('</itemize>', '')
        #itemize items in swig file
        text = text.replace('<item>', r'</newline> * ')
        text = text.replace('</item>', '')
        #displayed formulas in swig file
        text = text.replace('\[', r'</newline> $')
        text = text.replace('\]', '$')
        for i in range(10):
            text = text.replace('  ', ' ')
        #used in swig file output
        text = text.split(r'</newline>')
        strs = []
        for txt in text:
            strs.extend(textwrap.wrap(txt.lstrip('\n '), width=self.maxChar,
                initial_indent=' '*(start_pos+indent),
                subsequent_indent = ' '*indent))
        return  ('\n'.join(strs)).lstrip().replace('\\', '\\\\').replace('"', r'\"')


    def write_latex(self, out):
        # first handle glocal functions
        for entry in [x for x in self.content if x['type'] == 'global_function' and not x['ignore'] and not x['hidden'] \
                and 'test' not in x['Name']]:
            funcname = self.latexName(entry['Name'].replace('simuPOP::', '', 1))
            print >> out, '\\newcommand{\\%sRef}{' % funcname
            print >> out, '\n\\subsection{Function \\texttt{%s}\index{%s}}' % (funcname, funcname)
            if entry.has_key('Usage') and entry['Usage'] != '':
                func_name = entry['Usage'].split('(')[0]
                func_body = entry['Usage'][len(func_name):].lstrip('(').rstrip(')')
                print >> out, '\\par\n\\begin{funcdesc}{%s}{%s}\n\\par' % \
                    (self.latex_text(func_name), self.latex_text(func_body))
            else:
                print >> out, '\\par\n\\begin{funcdesc}{%s}{}\n\\par' % funcname
            if entry['Doc'] == '':
                print('Warning: no documentation for ', entry['Name'])
                print >> out, r'FIXME: No document.\par'
            else:
                print >> out, r'\MakeUppercase %s\par' % self.latex_text(entry['Doc'])
            if entry.has_key('Arguments') and entry['Arguments'] != '':
                for arg in entry['Arguments']:
                    print >> out, r'{\leftskip 0.3in \parindent=-0.3in \emph{%s: }\MakeUppercase %s\par}' \
                        % (self.latex_text(arg['Name']), self.latex_text(arg['Doc']))
            if entry.has_key('note') and entry['note'] != '':
                print >> out, '\\par\n\\strong{Note: }'
                print >> out, '    %s\n' % self.latex_text(entry['note'])
            if entry.has_key('ExampleFile') and entry['ExampleFile'] is not None:
                label = os.path.split(cons['ExampleFile'])[-1].split('_')[-1]
                title = self.latex_text(cons['ExampleTitle'])
                print >> out, '\\lstinputlisting[caption={%s},label={%s}]{%s}' % \
                    (title, label, entry['ExampleFile'].replace('\\', '/'))
            print >> out, '\\end{funcdesc}\n}\n'
        # then python modules
        modules = set(
            [x['module'] for x in self.content if (x['type'].startswith('module') or x['type'].startswith('docofmodule')) \
                and not x['ignore'] and not x['hidden']])
        for module in modules:
            # module functions
            #mod = [x for x in self.content if x['type'] == 'docofmodule_' + module][0]

            print >> out, '\\newcommand{\\%sRef}{\\index{%s}' % (module.replace('.', ''), module)
            doc = [x['Doc'] for x in self.content if x['type'] == 'docofmodule_' + module][0]
            print >> out, self.latex_formatted_text(doc)
            print >> out, '}\n'

            # module functions
            funcs = [x for x in self.content if x['type'] == 'module_func' and x['module'] == module and not x['ignore'] and not x['hidden']]
            # sort it
            funcs.sort(lambda x, y: cmp(x['Name'], y['Name']))
            # print all functions
            #print >> out, '\\begin{description}'
            for mem in funcs:
                print >> out, '\\newcommand{\\%s%sRef}{\\index{%s!%s}' % (module.replace('.',''), mem['Name'], module, mem['Name'])
                print >> out, '\n\\subsection{Function \\texttt{%s}\index{%s}}' % (mem['Name'], mem['Name'])
                if mem.has_key('Usage') and mem['Usage'] != '':
                    func_name = mem['Usage'].split('(')[0]
                    func_body = mem['Usage'][len(func_name):].lstrip('(').rstrip(')')
                    print >> out, '\\begin{funcdesc}{%s}{%s}\n' % (self.latex_text(func_name),
                        self.latex_text(func_body))
                else:
                    print >> out, '\\begin{funcdesc}{%s}{}\n' % mem['Name']
                print >> out, r' %s' % self.latex_formatted_text(mem['Doc'])
                if mem.has_key('note') and mem['note'] != '':
                    print >> out, '\\par\n\\strong{Note: }%s\\par' % self.latex_text(mem['note'])
                print >> out, '\\end{funcdesc}\n}\n'
            # module classes
            classes = [x for x in self.content if x['type'] == 'module_class' and x['module'] == module and not x['ignore'] and not x['hidden']]
            # sort it
            classes.sort(lambda x, y: cmp(x['Name'], y['Name']))
            # print all functions
            for cls in classes:
                print >> out, '\\newcommand{\\%s%sRef}{' % (module.replace('.', ''), cls['Name'])
                print >> out, '\n\\subsection{Class \\texttt{%s}\index{module!%s}' % (cls['Name'], cls['Name'])
                print >> out, '}\n'
                print >> out, r' %s' % self.latex_formatted_text(cls['Doc'])
                if cls.has_key('Usage') and cls['Usage'] != '':
                    func_name = cls['Usage'].split('(')[0]
                    func_body = cls['Usage'][len(func_name):].lstrip('(').rstrip(')')
                    print >> out, '\\begin{classdesc}{%s}{%s}\n' % (self.latex_text(func_name),
                        self.latex_text(func_body))
                else:
                    print >> out, '\\begin{classdesc}{%s}{}\n' % cls['Name']
                print >> out, self.latex_formatted_text(cls['InitDoc'])
                if not cls.has_key('Members') or len(cls['Members']) == 0:
                    print >> out, '\\end{classdesc}\n}\n'
                    continue
                else:
                    print >> out, '\\par\n\\vspace{3pt}\\par\n'
                for mem in cls['Members']:
                    if mem['Doc'] == '':
                        continue
                    usage = self.latex_text(mem['Usage'])
                    usageName = usage.split('(')[0]
                    usageParam = usage[len(usageName):].lstrip('(').rstrip(')')
                    print >> out, r'\begin{methoddesc}{%s}{%s}' % (usageName, usageParam)
                    if mem['Doc'] != '':
                        print >> out, '\n\\vspace{3pt} \\MakeUppercase %s\\par\n' % self.latex_formatted_text(mem['Doc'])
                    else:
                        print >> out, '\\hspace{0pt}\\par\n'
                    print >> out, '\\end{methoddesc}\n'
                print >> out, '\\end{classdesc}\n}\n'
        # then classes
        for entry in [x for x in self.content if x['type'] == 'class' and not x['ignore'] and not x['hidden']]:
            print >> out, '\\newcommand{\\%sRef}{' % self.latexName(entry['Name'].replace('simuPOP::', '', 1).replace('.', ''))
            classname = self.latex_text(entry['Name'].replace('simuPOP::', '', 1))
            print >> out, '\n\\subsection{Class \\texttt{%s}\\index{class!%s}}\n' % (classname, classname)
            print >> out, '\\par \\MakeUppercase %s' % self.latex_text(entry['Doc'])
            if entry.has_key('note') and entry['note'] != '':
                print >> out, '\\par\n\\strong{Note: }\n\\par'
                print >> out, '%s' % self.latex_text(entry['note'])
            # only use the first constructor
            constructor = [x for x in self.content if x['type'] == 'constructorofclass_' + entry['Name'] and not x['ignore'] and not x['hidden']]
            if len(constructor) == 0:
                print >> out, '}\n'
                continue
            elif len(constructor) > 1:
                print "Warning: multiple constructors: %s" % entry['Name']
            #print >> out, '\\vspace{6pt}\n'
            cons = constructor[0]
            #
            #print >> out, '\\par\n\\strong{Initialization}\n\\par'
            if cons.has_key('Usage') and cons['Usage'] != '':
                usage = self.latex_text(cons['Usage'])
                usageName = usage.split('(')[0]
                usageParam = usage[len(usageName):].lstrip('(').rstrip(')')
                print >> out, r'\begin{classdesc}{%s}{%s}' % (usageName, usageParam)
            else:
                print >> out, r'\begin{classdesc}{%s}{}' % entry['Name']
            if cons['Doc'] != '':
                print >> out, '\n\\vspace{3pt} \\MakeUppercase %s\\par\n' % self.latex_text(cons['Doc'])
            else:
                print >> out, '\\hspace{0pt}\\par\n'
            if cons.has_key('Arguments') and len(cons['Arguments']) > 0:
                print >> out, '\\par\n\n'
                # cons['Arguments'].sort(lambda x, y: cmp(x['Name'], y['Name']))
                for arg in cons['Arguments']:
                    #print >> out, r'{\emph{%s: }\MakeUppercase %s\par}' \
                    print >> out, r'{\leftskip 0.3in \parindent=-0.3in \emph{%s: }\MakeUppercase %s\par}' \
                        % (self.latex_text(arg['Name']), self.latex_text(arg['Description']))
            if cons.has_key('note') and cons['note'] != '':
                print >> out, '\\par\n\\strong{Note} '
                print >> out, '%s' % self.latex_text(cons['note'])
            members = [x for x in self.content if x['type'] == 'memberofclass_' + entry['Name'] and \
                       not x['ignore'] and not x['hidden'] and not '~' in x['Name'] and not '__' in x['Name']]
            def sort_member(x, y):
                if not x.has_key('group'):
                    print x
                if not y.has_key('group'):
                    print y
                res = cmp(x['group'], y['group'])
                if res == 0:
                    return cmp(x['Name'], y['Name'])
                else:
                    return res
            members.sort(sort_member)
            if len(members) == 0:
                print >> out, '\\end{classdesc}\n}\n'
                continue
            group = ''
            for mem in members:
                #print "MEMBER %s, GROUP '%s'" % (mem['Name'], mem['group'])
                if group != mem['group']:
                    if group != '':
                        print >> out, '\\vspace{10pt}\n'
                    group = mem['group']
                if mem.has_key('Usage') and mem['Usage'] != '':
                    usage = self.latex_text(mem['Usage'])
                    assert usage.startswith('x.')
                    usageName = usage.split('(')[0][2:]
                    usageParam = usage[len(usageName)+2:].lstrip('(').rstrip(')')
                    print >> out, r'\begin{methoddesc}{%s}{%s}' % (usageName, usageParam)
                else:
                    print >> out, r'\begin{methoddesc}{%s}{}' % mem['Name'].split(':')[-1]
                print >> out, '\\MakeUppercase %s\n' % self.latex_text(mem['Doc'])
                if mem.has_key('Arguments') and mem['Arguments'] != '':
                    # mem['Arguments'].sort(lambda x, y: cmp(x['Name'], y['Name']))
                    for arg in mem['Arguments']:
                        print >> out, r'{\leftskip 0.3in \parindent=-0.3in \emph{%s: }\MakeUppercase %s\par}' % \
                            (self.latex_text(arg['Name']), self.latex_text(arg['Description']))
                if mem.has_key('note') and mem['note'] != '':
                    print >> out, '\\par\n\\strong{Note:} %s\\par' % self.latex_text(mem['note'])
                print >> out, r'\end{methoddesc}'
            if cons.has_key('ExampleFile') and cons['ExampleFile'] is not None:
                print >> out, '\\strong{Example}\n'
                # ../log/ref_xxx.log => xxx.log
                label = os.path.split(cons['ExampleFile'])[-1].split('_')[-1]
                title = self.latex_text(cons['ExampleTitle'])
                print >> out, '\\lstinputlisting[caption={%s},label={%s}]{%s}' % \
                    (title, label, cons['ExampleFile'].replace('\\', '/'))
            print >> out, '\\end{classdesc}\n}\n'

    def wrap_reST(self, input, initial_indent = '   ', subsequent_indent = '   '):
        text = input
        for tag in ['<tt>', '<bf>', '<em>', '<item>']:
            text = re.sub('%s\s*' % tag, tag, text)
        # highlight stuff
        for keyword in self.auto_keywords.keys():
            for tag in ['tt', 'bf']:
                for cls in self.auto_keywords[keyword]:
                    text = text.replace('<%s>%s</%s>' % (tag, cls, tag),
                        ':%s:`%s`' % (keyword, cls))
        text = text.replace('<em>', '*')
        text = text.replace('</em>', '*')
        text = text.replace('<tt>', '``')
        text = text.replace('</tt>', '``')
        text = text.replace('<bf>', '**')
        text = text.replace('</bf>', '**')
        text = text.replace('</newline>', '\n\n')
        text = text.replace('<itemize>', '\n')
        text = text.replace('</itemize>', '')
        text = text.replace('<item>', '\n+ ')
        text = text.replace('</item>', '\n')
        text = text.split('\n')
        txt = []
        for t in text:
            if t.strip() == '':
                continue
            txt.append('')
            # itemize
            if t.startswith('+'):
                txt.extend(textwrap.wrap(t.strip(), initial_indent = initial_indent,
                    subsequent_indent = '  ' + subsequent_indent))
            else:
                txt.extend(textwrap.wrap(t.strip(), initial_indent = initial_indent,
                    subsequent_indent=subsequent_indent))
        return ('\n'.join(txt)).strip('\n')

    def shiftText(self, txt, shift='   '):
        # highlight stuff
        text = txt
        for keyword in self.auto_keywords.keys():
            for tag in ['``', '**']:
                for cls in self.auto_keywords[keyword]:
                    text = text.replace('%s%s%s' % (tag, cls, tag),
                        ':%s:`%s`' % (keyword, cls))
        return '\n'.join([shift + x for x in text.split('\n')])

    def write_reST(self, dir):
        # first handle glocal functions
        for entry in [x for x in self.content if x['type'] == 'global_function' and not x['ignore'] and not x['hidden'] \
                and 'test' not in x['Name']]:
            refName = '%sRef.ref' % entry['Name'].replace('simuPOP::', '', 1).replace('.', '')
            print 'Writing reference for global function ', refName
            out = open(os.path.join(dir, refName), 'w')
            funcname = entry['Name'].replace('simuPOP::', '', 1)
            print >> out, '\nFunction %s' % funcname
            print >> out, '-' * (9 + len(funcname))
            print >> out
            print >> out, '\n.. function::',
            if entry.has_key('Usage') and entry['Usage'] != '':
                print >> out, entry['Usage']
            else:
                print >> out, '%s()' % funcname
            #
            print >> out
            if entry['Doc'] != '':
                print >> out, self.wrap_reST(entry['Doc'])

            if entry.has_key('note') and entry['note'] != '':
                print >> out, '**Note**:'
                print >> out, self.wrap_reST(entry['note'])
            out.close()
        # then python modules
        modules = set(
            [x['module'] for x in self.content if (x['type'].startswith('module') or x['type'].startswith('docofmodule')) \
                and not x['ignore'] and not x['hidden']])
        for module in modules:
            # module functions
            # MODULE DOC
            refName = '%sRef.ref' % module.replace('.', '')
            out = open(os.path.join(dir, refName), 'w')
            #print >> out, ':mod:`%s`' % module
            #print >> out, '='*(len(module)+7)
            if module != 'simuPOP':
                print >> out
                print >> out, '.. module:: %s\n' % module
            doc = [x['Doc'] for x in self.content if x['type'] == 'docofmodule_' + module][0]
            print >> out, doc
            out.close()
            #
            # module functions
            funcs = [x for x in self.content if x['type'] == 'module_func' \
                and x['module'] == module and not x['ignore'] and not x['hidden']]
            # sort it
            funcs.sort(lambda x, y: cmp(x['Name'], y['Name']))
            # print all functions
            #print >> out, '\\begin{description}'
            for mem in funcs:
                refName = '%s%sRef.ref' % (module.replace('.', ''), mem['Name'])
                print 'Writing reference for function ', refName
                out = open(os.path.join(dir, refName), 'w')
                funcname = mem['Name']
                print >> out, '\nFunction %s' % funcname
                print >> out, '-' * (9 + len(funcname))
                print >> out
                print >> out, '\n.. function::',
                if mem.has_key('Usage') and mem['Usage'] != '':
                    print >> out, mem['Usage']
                else:
                    print >> out, '%s()' % funcname
                print >> out
                print >> out, self.shiftText(mem['Doc'])
                if mem.has_key('note') and mem['note'] != '':
                    print >> out, '**Note**'
                    print >> out, self.shiftText(mem['note'])
                out.close()
            # module classes
            classes = [x for x in self.content if x['type'] == 'module_class' and x['module'] == module and not x['ignore'] and not x['hidden']]
            # sort it
            classes.sort(lambda x, y: cmp(x['Name'], y['Name']))
            # print all functions
            for cls in classes:
                refName = '%s%sRef.ref' % (module.replace('.', ''), cls['Name'])
                print 'Writing reference for class ', refName
                out = open(os.path.join(dir, refName), 'w')
                classname = cls['Name']
                print >> out, '\nclass %s' % classname
                print >> out, '-' * (6 + len(classname))
                print >> out
                print >> out, '.. class:: %s\n' % classname
                print >> out, self.shiftText(cls['Doc'])
                print >> out
                if cls.has_key('Usage') and cls['Usage'] != '':
                    print >> out, '   .. method:: %s.%s' % (classname, cls['Usage'])
                else:
                    print >> out, '   .. method:: %s.%s()' % (classname, classname)
                print >> out
                print >> out, self.shiftText(cls['InitDoc'], ' '*6)
                print >> out
                #
                for mem in cls['Members']:
                    usage = self.latex_text(mem['Usage'])
                    if mem.has_key('Usage') and mem['Usage'] != '':
                        print >> out, '   .. method:: %s.%s' % (classname, mem['Usage'])
                    else:
                        print >> out, '   .. method:: %s.%s()' % (classname, mem['Name'])
                    print >> out
                    if mem['Doc'] != '':
                        print >> out, self.shiftText(mem['Doc'], ' '*6)
                    print >> out
                out.close()
        # then classes
        for entry in [x for x in self.content if x['type'] == 'class' and not x['ignore'] and not x['hidden']]:
            refName = '%sRef.ref' % entry['Name'].replace('simuPOP::', '', 1).replace('.', '')
            print 'Writing reference for class ', refName
            out = open(os.path.join(dir, refName), 'w')
            classname = self.latex_text(entry['Name'].replace('simuPOP::', '', 1))
            print >> out, '\nclass %s' % classname
            print >> out, '-' * (6 + len(classname))
            print >> out
            print >> out, '.. class::', classname
            print >> out
            print >> out, self.wrap_reST(entry['Doc'])
            if entry.has_key('note') and entry['note'] != '':
                print >> out, '\n   .. note::\n'
                print >> out, '%s' % self.wrap_reST(entry['note'].strip(), ' '*6, ' '*6)
            print >> out
            # only use the first constructor
            constructor = [x for x in self.content if x['type'] == 'constructorofclass_' + entry['Name'] and not x['ignore'] and not x['hidden']]
            if len(constructor) == 0:
                continue
            elif len(constructor) > 1:
                print "Warning: multiple constructors: %s" % entry['Name']
            cons = constructor[0]
            #
            #print >> out, '\\par\n\\strong{Initialization}\n\\par'
            if cons.has_key('Usage') and cons['Usage'] != '':
                usage = self.latex_text(cons['Usage'])
                print >> out, '\n   .. method:: %s\n' % cons['Usage']
            else:
                print >> out, '\n   .. method:: %s()\n' % entry['Name']
            if cons['Doc'] != '':
                print >> out, '\n%s\n' % self.wrap_reST(cons['Doc'], ' '*6, ' '*6)
            if cons.has_key('Arguments') and len(cons['Arguments']) > 0:
                # cons['Arguments'].sort(lambda x, y: cmp(x['Name'], y['Name']))
                for arg in cons['Arguments']:
                    #print >> out, r'{\emph{%s: }\MakeUppercase %s\par}' \
                    print >> out, '      %s\n         %s\n' % (arg['Name'], self.wrap_reST(arg['Description'], '', ' '*9))
            if cons.has_key('note') and cons['note'] != '':
                print >> out, '      .. note::\n'
                print >> out, '%s' % self.wrap_reST(cons['note'], ' '*9, ' '*9)
            #
            members = [x for x in self.content if x['type'] == 'memberofclass_' + entry['Name'] and \
                       not x['ignore'] and not x['hidden'] and not '~' in x['Name'] and not '__' in x['Name']]
            def sort_member(x, y):
                if not x.has_key('group'):
                    print x
                if not y.has_key('group'):
                    print y
                res = cmp(x['group'], y['group'])
                if res == 0:
                    return cmp(x['Name'], y['Name'])
                else:
                    return res
            members.sort(sort_member)
            if len(members) == 0:
                continue
            group = ''
            for mem in members:
                if group != mem['group']:
                    group = mem['group']
                if mem.has_key('Usage') and mem['Usage'] != '':
                    print >> out, '\n   .. method:: %s.%s\n' % (classname, mem['Usage'].lstrip('x.'))
                else:
                    print >> out, '\n   .. method:: %s.%s()\n' % (classname, mem['Name'].split(':')[-1])
                print >> out, self.wrap_reST(mem['Doc'], ' '*6, ' '*6)
            out.close()

    def write_latex_testfile(self, out, ref_file):
        print >> out, r'''\documentclass[oneside,english]{manual}
\usepackage[T1]{fontenc}
\usepackage[latin9]{inputenc}
\usepackage{listings}
\usepackage{url}
\lstset{basicstyle={\ttfamily},
language=Python,
showspaces=false,
showstringspaces=false,
showtabs=false,
xleftmargin=15pt}
\setcounter{secnumdepth}{3}
\setcounter{tocdepth}{3}
\usepackage{verbatim}
\usepackage{amsmath}
\usepackage{makeidx}
\makeindex
\usepackage{listings}

\newlength{\admonitionwidth} 
\setlength{\admonitionwidth}{0.7\textwidth}


\makeatletter
\usepackage{babel}
\makeatother

\begin{document}

\title{simuPOP Reference Manual}

\maketitle

\lstlistoflistings
\include{%s}''' % os.path.basename(os.path.splitext(ref_file)[0])
        global_funcs = [x for x in self.content if x['type'] == 'global_function' and not x['ignore'] and not x['hidden'] \
                and 'test' not in x['Name']]
        global_funcs.sort(lambda x, y: cmp(x['Name'], y['Name']))
        for entry in global_funcs:
             print >> out, r'\%sRef' % self.latexName(entry['Name'].replace('simuPOP::', '', 1))
             print >> out, r'\vspace{.5in}\par\rule[.5ex]{\linewidth}{1pt}\par\vspace{0.3in}'
        # modules
        modules = sets.Set(
            [x['module'] for x in self.content if x['type'].startswith('module') and not x['ignore'] and not x['hidden']])
        for module in modules:
            funcs = [x for x in self.content if x['type'] == 'module_func' and x['module'] == module and not x['ignore'] and not x['hidden']]
            # sort it
            funcs.sort(lambda x, y: cmp(x['Name'], y['Name']))
            # print all functions
            for mem in funcs:
                print >> out, r'\%s%sRef' % (module.replace('.',''), mem['Name'])
                print >> out, r'\vspace{.5in}\par\rule[.5ex]{\linewidth}{1pt}\par\vspace{0.3in}'
            classes = [x for x in self.content if x['type'] == 'module_class' and x['module'] == module and not x['ignore'] and not x['hidden']]
            # print all functions
            for cls in classes:
                print >> out, r'\%s%sRef' % (module.replace('.', ''), cls['Name'])
                print >> out, r'\vspace{.5in}\par\rule[.5ex]{\linewidth}{1pt}\par\vspace{0.3in}'
        for entry in [x for x in self.content if x['type'] in ['class'] and not x['ignore'] and not x['hidden']]:
             print >> out, r'\%sRef' % self.latexName(entry['Name'].replace('simuPOP::', '', 1).replace('.',''))
             print >> out, r'\vspace{.1in}\par\rule[.3ex]{\linewidth}{1pt}\par\vspace{0.1in}'
        print >> out, r'\end{document}'


    def write(self, output, type, ref_file=''):
        # write in binary (b) mode so that \n is written as \n, not \r\n under windows.
        if type == 'reST':
            self.write_reST(output)
            return
        fout = open(output, 'wb')
        if type == 'swig':
            self.write_swig(fout)
        elif type == 'latex_single':
            self.write_latex(fout)
        elif type == 'latex_all':
            self.write_latex_testfile(fout, ref_file)
        fout.close()

       
if __name__ == '__main__':
    if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
        print __doc__
        sys.exit(1)
    # to make my life a bit easier, provide some default parameters
    if os.path.isfile('doxy2swig.py'):
        src_path = '..'
    else:
        src_path = '.'
    # doxygen generated xml file
    if len(sys.argv) < 2:
        xml_file = os.path.join(src_path, 'doxygen_doc', 'xml', 'index.xml')
    else:
        xml_file = sys.argv[1]
    # output interface file
    if len(sys.argv) < 3:
        interface_file = os.path.join(src_path, 'src', 'simuPOP_doc.i')
    else:
        interface_file = sys.argv[2]
    # output ref_single.tex file
    if len(sys.argv) < 4:
        latex_file = os.path.join(src_path, 'doc', 'simuPOP_ref.tex')
    else:
        latex_file = sys.argv[3]
    # output ref_all.tex file
    if len(sys.argv) < 5:
        latex_testfile = os.path.join(src_path, 'doc', 'simuPOP_ref_test.tex')
    else:
        latex_testfile = sys.argv[4]
    if not os.path.isdir(os.path.join(src_path, 'doc', 'build')):
        os.mkdir(os.path.join(src_path, 'doc', 'build'))
    refFile = os.path.join(src_path, 'doc', 'build', 'reflist.py')
    # read the XML file (actually a index.xml file that contains all others)
    p = Doxy2SWIG(xml_file)
    # generate interface file.
    p.generate()
    # clean up, and process CPPONLY etc
    p.post_process()
    # write interface file to output interface file.
    print 'Writing SWIG interface file to', interface_file
    p.write(interface_file, type='swig')
    sys.path = [os.path.join(src_path, 'src')] + sys.path
    p.scan_module('simuPOP')
    p.scan_module('simuOpt')
    p.scan_module('simuPOP.utils')
    p.scan_module('simuPOP.demography')
    p.scan_module('simuPOP.plotter')
    p.scan_module('simuPOP.sampling')
    p.scan_module('simuPOP.sandbox')
    p.scan_module('simuPOP.gsl')
    #
    for func,realClass,realFunc,group in [
        ('evolve_pop', 'Population', 'evolve', '7-evolve'),
        ('all_individuals', 'Population', 'allIndividuals', '4-ind'),
        ('as_pedigree', 'Population', 'asPedigree', '1-popConvert'),
        ('as_population', 'Pedigree', 'asPopulation', '1-pedConvert')]:
        entry = [x for x in p.content if func in x['Name']][0]
        entry['ignore'] = True
        entry['hidden'] = True
        p.content.append({'Name': u'simuPOP::%s::%s' % (realClass, realFunc),
             'type': u'memberofclass_simuPOP::%s' % realClass,
             'Description': '',
             'Details': ur'<group>%s</group>' % group + entry['Description'].replace('\n', ' '),
             'Usage': 'x.' + entry['Usage'].replace(func, realFunc).replace('self, ', '').replace('self)', ')'),
        })
    p.fix_entries()
    # hange a few usages:
    #
    print 'Writing latex reference file to', latex_file
    p.write(latex_file, type='latex_single')
    p.uniqueName = []
    #
    # list of classes and functions....
    #
    module_classes = set(['%s.%s' % (x['module'], x['Name']) for x in p.content if x['type'] == 'module_class' \
        and not x['ignore'] and not x['hidden']])
    module_funcs = set([('%s.%s') % (x['module'], x['Name']) for x in p.content if x['type'] == 'module_func' \
        and not x['ignore'] and not x['hidden']])
    simuPOP_classes = set([x.lstrip('simuPOP.') for x in module_classes if x.startswith('simuPOP.')])
    simuPOP_funcs = set([x.lstrip('simuPOP.') for x in module_funcs if x.startswith('simuPOP.')])
    modules = set([x['Name'] for x in p.content if x['type'].startswith('docofmodule_')])
    global_funcs = set([x['Name'].replace('simuPOP::', '') for x in p.content if x['type'] == 'global_function' and not x['ignore'] and not x['hidden'] \
            and 'test' not in x['Name']])
    classes = set([x['Name'].replace('simuPOP::', '') for x in p.content if x['type'] == 'class' and not x['ignore'] and not x['hidden']])
    members = set([x['Name'].replace('simuPOP::', '').replace('::', '.') \
            for x in p.content if x['type'].startswith('memberofclass_') and \
            not x['ignore'] and not x['hidden'] and not '~' in x['Name'] and not '__' in x['Name']])
    # class members;
    module_class_members = []
    for mod in modules:
        for cls in [x for x in p.content if x['type'] == 'module_class' and x['module'] == mod and not x['ignore'] and not x['hidden']]:
            if cls.has_key('Members'):
                module_class_members.extend(['%s.%s.%s' % (mod, cls['Name'], x['Name']) for x in cls['Members']])
                #print 'Extending', ['%s.%s' % (cls['Name'], x['Name']) for x in cls['Members']]
    p.auto_keywords =  {'func': list(global_funcs | module_funcs | simuPOP_funcs),
        'class': list(classes | module_classes | simuPOP_classes),
        'meth': list(members) + module_class_members,
        'mod': list(modules)}
    try:
        lst = open(refFile, 'w')
        print >> lst, 'auto_keywords = \\'
        pprint.pprint(p.auto_keywords, stream = lst)
        lst.close()
    except:
        print 'Failed to write a list file for cross referencing purposes'
    #
    builddir = os.path.join(src_path, 'doc', 'build')
    if not os.path.isdir(builddir):
        os.mkdir(builddir)
    p.write(builddir, type='reST')
    # clear unique name
    p.uniqueName = []
    print 'Writing latex test file to', latex_testfile
    p.write(latex_testfile, type='latex_all', ref_file=latex_file)
    # generating sample document
    os.chdir(os.path.dirname(latex_testfile))
    os.system('pdflatex %s' % os.path.split(latex_testfile)[1])
    # ending statement
    print 'Done.'
