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
import re, textwrap, sys, types, os.path, sets

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
        if node.firstChild is not None:
            self.add_text(r'<em>%s</em>' % node.firstChild.data)
            

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
                    self.add_text( self.format_text( defName, 0, 0 ) )
                else:
                    self.add_text( 'x.' )
                    self.add_text( self.format_text( defName, 0, 0 ) )

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
            # get file name
            filename = node.firstChild.nextSibling.firstChild.firstChild.data
            # get content as string
            try:
                try:
                    # usual ../doc/log directory
                    file = open(os.path.join('..', 'doc', 'log', filename.strip()))
                    self.content[-1]['ExampleFile'] = os.path.join('..', 'doc', 'log', filename.strip())
                except:
                    # local file
                    file = open(filename.strip() )
                    self.content[-1]['ExampleFile'] = filename.strip()
                cont = file.read()
                file.close()
                self.add_text(cont)
            except:
                self.content[-1]['ExampleFile'] = None
                print "File " + filename + " does not exist\n"
                self.add_text(filename + "does not exist\n")


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
        self.content[-1]['cppArgs'] = txt
        # replace the trailing const
        # @ is used instead of , to avoid separation of replaced text, it will be replaced back to ,
        txt = txt.replace('vectorstr(TAG_InheritFields, TAG_InheritFields+2)',
            '["paternal_tag"@ "maternal_tag"]')
        txt = txt.replace('vectorstr(TAG_ParentsFields, TAG_ParentsFields+2)',
            '["father_idx"@ "mother_idx"]')
        txt = txt.replace('vectorstr(ASC_AS_Fields, ASC_AS_Fields+2)',
            '["father_idx"@ "mother_idx"]')
        # re function used to replace the following sentances
        vec1 = re.compile('(.*)vectorstr\(1,\s*([\w"]+)\)(.*)')
        txt = vec1.sub(r'\1[\2]\3', txt)
        #txt = txt.replace('vectorstr(1, "qtrait")', '["qtrait"]')
        con1 = re.compile('\)\s*const\s*$')
        txt = con1.sub(')', txt)
        #txt = txt.replace(')    const',')')
        args = txt.split(',')
        out=[]
        for s in args:
            #use @ in order to split the arguments correctly
            s = s.replace('@', r',')
            piece = s.split('=')
            var = piece[0].split(' ')[-1].split(')')[0].split('(')[-1]
            #delete &
            var = var.replace('&', '')
            if( len( piece ) == 2 ):
                defVal = piece[1].split('(')[0].split(')')[0].split(')')[0]
                # re function used to repalce the following sentances
                vect = re.compile('vector(lu|u|l|i|f|a|op|str)')
                defVal = vect.sub('[]', defVal)
                #defVal = defVal.replace('vectorlu','[]')
                #defVal = defVal.replace('vectorstr','[]')
                defVal = defVal.replace('strDict','{}')
                defVal = defVal.replace('intDict','{}')
                defVal = defVal.replace('matrix','[]')
                defVal = defVal.replace('intMatrix','[]')
                defVal = defVal.replace('true','True')
                defVal = defVal.replace('false','False')
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
            print "parsing file: %s"%fname
            try:
                p = Doxy2SWIG(fname)
                p.generate()
                self.content.extend(p.content)
            except:
                print "This file can not be parsed, something wrong with file format"
                pass


    def post_process(self):
        # remove duplicate entry
        # They might be introduced if a function is list both under 'file' and under 'namespace'
        print "Number of entries: ", len(self.content)
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
        # remove all empty entries 
        def notEmpty(entry):
            'judge if an entry is empty'
            return (entry.has_key('Description') and entry['Description'] != '') or \
                (entry.has_key('Details') and entry['Details'] != '')
        for entry in self.content:
            if not (entry.has_key('Description') and entry['Description'] != '') and \
                not (entry.has_key('Details') and entry['Details'] != ''):
                entry['Description'] = entry['Name']
        # mark all entries with 'CPPONLY' in description or details as ignore
        for entry in self.content:
            if (entry.has_key('Description') and 'CPPONLY' in entry['Description']) or \
               (entry.has_key('Details') and 'CPPONLY' in entry['Details']):
                entry['ignore'] = True
            else:
                entry['ignore'] = False
        # add funcForm key to content and delete function form from Details 
        for entry in self.content:
            if (entry.has_key('Details') and '<funcForm>' in entry['Details']):
                piece1 = entry['Details'].split('<funcForm>')
                piece2 = piece1[1].split('</funcForm>')
                entry['Details'] = piece1[0] + piece2[1]
                entry['funcForm'] = '<tt>' + piece2[0] + '</tt>'
        # destructor can not be CPPONLY, this will cause memory leak. Let us check this
        for entry in self.content:
            if entry['ignore'] and '~' in entry['Name']:
                print "Desctructor of %s has CPPONLY. Please correct it." % entry['Name']
                sys.exit(1)
        #for entry in self.content:
        #    if entry.has_key('description') and 'PLOIDY:' in entry['description'];
        #        if 'PLOIDY:ALL' in entry['description']:
        #            self.content['ploidy'] = 'all'
        #            self.content['description'].remove('PLOIDY:ALL')
        # sort the entries
        self.content.sort(lambda x, y: x['Name'] > y['Name'])
    

    def write_swig(self, out):
        for entry in self.content:
            if entry['ignore']:
                if entry.has_key('cppArgs'):
                    print >> out, '%%ignore %s%s;\n' % (entry['Name'], entry['cppArgs'])
                else:
                    print >> out, '%%ignore %s;\n' % entry['Name']
                continue          
            print >> out, '%%feature("docstring") %s "\n' % entry['Name']
            if entry.has_key('funcForm'):
                print >> out, 'Function form:'
                print >> out, '\n    %s\n' % self.format_text(entry['funcForm'], 0, 4)
            if entry.has_key('Description') and entry['Description'] != '':
                print >> out, 'Description:'
                print >> out, '\n    %s\n' % self.format_text(entry['Description'], 0, 4)            
            if entry.has_key('Usage') and entry['Usage'] != '':
                print >> out, 'Usage:'
                print >> out, '\n    %s\n' % self.format_text(entry['Usage'], 0, 6)
            if entry.has_key('Details') and entry['Details'] != '':
                print >> out, 'Details:'
                print >> out, '\n    %s\n' % self.format_text(entry['Details'], 0, 4)
            if entry.has_key('Arguments') and entry['Arguments'] != '':
                print >> out, 'Arguments:\n'
                for arg in entry['Arguments']:
                    print >> out, '    %-16s%s' % (arg['Name']+':', self.format_text(arg['Description'], 0, 20))
                print >> out
            if entry.has_key('note') and entry['note'] != '':
                print >> out, 'Note:'
                print >> out, '\n    %s\n' % self.format_text(entry['note'], 0, 4)
            if entry.has_key('Examples') and entry['Examples'] != '':
                print >> out, 'Examples:'
                print >> out, '\n%s\n' % entry['Examples'].replace('\\', r'\\\\').replace('"', r'\"')
            #if len(str) > 2048:
            #    #print 'Entry %s is too long (%d)' % (entry['Name'], len(str))
            #    str = str[0:1970] + '...' + '\n\nPlease refer to the reference manual for more details.\n\n'
            print >> out, '\"; \n'


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


    def write_latex(self, out):
        # first handle glocal functions
        for entry in [x for x in self.content if x['type'] == 'global_function' and not x['ignore'] \
                and 'test' not in x['Name']]:
            print >> out, '\\newcommand{\\%sRef}{' % self.latexName(entry['Name'].replace('simuPOP::', '', 1))
            if entry.has_key('Description') and entry['Description'] != '':
                print >> out, '\\par\n\\strong{Function \\texttt{%s}}\n\\par' % self.latex_text(entry['Name'].replace('simuPOP::', '', 1))
                print >> out, '%s\par' % self.latex_text(entry['Description'])
            if entry.has_key('Usage') and entry['Usage'] != '':
                print >> out, '\\begin{quote}\\function{%s}\\end{quote}' % self.latex_text(entry['Usage'])
            if entry.has_key('Details') and entry['Details'] != '':
                print >> out, '\\par\n\\strong{Details}\n\\par'
                print >> out, '    %s\n' % self.latex_text(entry['Details'])
            if entry.has_key('Arguments') and entry['Arguments'] != '':
                print >> out, '\\par\n\\strong{Arguments}'
                print >> out, '\\begin{description}'
                for arg in entry['Arguments']:
                    print >> out, '\\item [{%s}]%s' % (arg['Name'], self.latex_text(arg['Description']))
                print >> out, '\\end{description}'
            if entry.has_key('note') and entry['note'] != '':
                print >> out, '\\par\n\\strong{Note}\n\\par'
                print >> out, '    %s\n' % self.latex_text(entry['note'])
            if entry.has_key('Examples') and entry['Examples'] != '':
                print >> out, '\\strong{Examples}\n\\begin{algorithm}[h]'
                print >> out, '\\caption{\\label{alg:%s}Example for function %s}' % (entry['Name'], entry['Name'])
                if entry['ExampleFile'] is not None:
                    print >> out, '\\verbatiminput{%s}' % entry['ExampleFile']
                print >> out, '\\end{algorithm}'
            print >> out, '}\n'
        # then classes
        for entry in [x for x in self.content if x['type'] == 'class' and not x['ignore']]:
            print >> out, '\\newcommand{\\%sRef}{' % self.latexName(entry['Name'].replace('simuPOP::', '', 1))
            if entry.has_key('Description') and entry['Description'] != '':
                print >> out, '\\par\n\\strong{Class \\texttt{%s}}' % self.latex_text(entry['Name'].replace('simuPOP::', '', 1))
                if entry.has_key('funcForm'):
                    print >> out, '  (Function form: %s)\n\\par' % self.latex_text(entry['funcForm'])
                print >> out, '%s' % self.latex_text(entry['Description'])
            if entry.has_key('Details') and entry['Details'] != '':
                print >> out, '\\par\n\\strong{Details}\n\\par'
                print >> out, '%s' % self.latex_text(entry['Details'])
            if entry.has_key('note') and entry['note'] != '':
                print >> out, '\\par\n\\strong{Note}\n\\par'
                print >> out, '%s' % self.latex_text(entry['note'])
            # only use the first constructor
            constructor = [x for x in self.content if x['type'] == 'constructorofclass_' + entry['Name'] and not x['ignore']]          
            if len(constructor) == 0:
                print >> out, '}\n'
                continue
            elif len(constructor) > 1:
                print "Warning: multiple constructors: %s" % entry['Name']
            cons = constructor[0]
            print >> out, '\\par\n\\strong{Initialization}\n\\par'
            if cons.has_key('Description') and cons['Description'] != '':
                print >> out, '%s\par' % self.latex_text(cons['Description'])
            if cons.has_key('Usage') and cons['Usage'] != '':
                print >> out, '\\begin{quote}\\function{%s}\\end{quote}' % self.latex_text(cons['Usage'])
            if cons.has_key('Details') and cons['Details'] != '':
                print >> out, '%s\par' % self.latex_text(cons['Details'])
            #if cons.has_key('Details') and cons['Details'] != '':
            #    print >> out, '\\par\n\\strong{Details}\n\\par'
            #    print >> out, '%s' % self.latex_text(cons['Details'])
            if cons.has_key('Arguments') and cons['Arguments'] != '':
                print >> out, '\\par\n'
                print >> out, '\\begin{description}'
                cons['Arguments'].sort(lambda x, y: cmp(x['Name'], y['Name']))
                for arg in cons['Arguments']:
                    print >> out, '\\item [{%s}]{\\leftskip 0.5cm %s\par}' % (self.latex_text(arg['Name']), self.latex_text(arg['Description']))
                print >> out, '\\end{description}'
            if cons.has_key('note') and cons['note'] != '':
                print >> out, '\\par\n\\strong{Note}\n\\par'
                print >> out, '%s' % self.latex_text(cons['note'])
            members = [x for x in self.content if x['type'] == 'memberofclass_' + entry['Name'] and \
                       not x['ignore'] and not '~' in x['Name'] and not '__' in x['Name']]
            members.sort(lambda x, y: cmp(x['Name'], y['Name']))
            if len(members) == 0:
                print >> out, '}\n'
                continue
            print >> out, '\\par\n\\strong{Member Functions}\n\\par'
            print >> out, '\\begin{description}'
            for mem in members:
                if mem.has_key('Usage') and mem['Usage'] != '':
                    print >> out, '\\item [\\function{%s}] ' % self.latex_text(mem['Usage'])
                if mem.has_key('Description') and mem['Description'] != '':
                    print >> out, '%s' % self.latex_text(mem['Description'])
                if mem.has_key('Details') and mem['Details'] != '':
                    # if we have short description, use a separate paragraph for details.
                    if mem.has_key('Description') and mem['Description'] != '':
                        print >> out, '\\par\n%s\\par' % self.latex_text(mem['Details'])
                    # otherwise, use details as description
                    else:
                        print >> out, '%s' % self.latex_text(mem['Details'])
                if mem.has_key('Arguments') and mem['Arguments'] != '':
                    print >> out, '\\begin{description}'
                    mem['Arguments'].sort(lambda x, y: cmp(x['Name'], y['Name']))
                    for arg in mem['Arguments']:
                        print >> out, '\\item [{%s}]%s' % (self.latex_text(arg['Name']), self.latex_text(arg['Description']))
                    print >> out, '\\end{description}'
                if mem.has_key('note') and mem['note'] != '':
                    print >> out, '\\par\n\\strong{Note:} %s\\par' % self.latex_text(mem['note'])
            print >> out, '\\end{description}'
            if cons.has_key('Examples') and cons['Examples'] != '':
                print >> out, '\\strong{Examples}\n'
                #print >> out, '\\caption{\\label{alg:%s}Example for function %s}' % (cons['Name'], cons['Name'])
                if cons['ExampleFile'] is not None:
                    print >> out, '\\lstinputlisting{%s}' % cons['ExampleFile'].replace('\\', '/')
            print >> out, '}\n'


    def write_latex_testfile(self, out, ref_file):
        print >> out, r'''\documentclass[oneside,english]{manual}
\usepackage[T1]{fontenc}
\usepackage[latin9]{inputenc}
\setcounter{secnumdepth}{3}
\setcounter{tocdepth}{3}
\usepackage{verbatim}
\usepackage{float}
\usepackage{listings}
\makeatletter

\providecommand{\tabularnewline}{\\}
\floatstyle{ruled}
\newfloat{algorithm}{tbp}{loa}
\floatname{algorithm}{Algorithm}

\newenvironment{lyxcode}
{\begin{list}{}{
\setlength{\rightmargin}{\leftmargin}
\setlength{\listparindent}{0pt}
\raggngth{\itemsep}{0pt}
\setlength{\parsep}{0pt}
\normalfont\ttfamily}
 \item[]}
{\end{list}}
\floatname{algorithm}{Example}

\lstset{language=Python,
   basicstyle=\ttfamily,
   showstringspaces=false,
   xleftmargin=30pt,
   xrightmargin=10pt,
   resetmargins=false,
   showspaces=false,
   showtabs=false,
   frame=none,
   columns=fixed,
   keepspaces=true,
}

\usepackage{babel}
\makeatother
\begin{document}
\include{%s}''' % os.path.basename(os.path.splitext(ref_file)[0])
        for entry in [x for x in self.content if x['type'] in ['global_function'] and not x['ignore'] \
            and 'test' not in x['Name']]:
             print >> out, r'\%sRef' % self.latexName(entry['Name'].replace('simuPOP::', '', 1))
             print >> out, r'\vspace{.5in}\par\rule[.5ex]{\linewidth}{1pt}\par\vspace{0.3in}'
        for entry in [x for x in self.content if x['type'] in ['class'] and not x['ignore']]:
             print >> out, r'\%sRef' % self.latexName(entry['Name'].replace('simuPOP::', '', 1))
             print >> out, r'\vspace{.1in}\par\rule[.3ex]{\linewidth}{1pt}\par\vspace{0.1in}'
        print >> out, r'\end{document}'


    def write(self, output, type, ref_file=''):
        fout = open(output, 'w')
        if type == 'swig':
            self.write_swig(fout)
        elif type == 'latex_single':
            self.write_latex(fout)
        elif type == 'latex_all':
            self.write_latex_testfile(fout, ref_file)
        fout.close()


    def format_text(self, text, start_pos, indent):
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
        #used in swig file output
        text = text.split(r'</newline>')
        strs = []
        for txt in text:
            strs.extend(textwrap.wrap(txt.lstrip('\n '), width=self.maxChar,
                initial_indent=' '*(start_pos+indent),
                subsequent_indent = ' '*indent))
        return  ('\n'.join(strs)).lstrip().replace('\\', r'\\\\').replace('"', r'\"')              


if __name__ == '__main__':
    if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
        print __doc__
        sys.exit(1)
    # to make my life a bit easier, provide some default parameters
    # doxygen generated xml file
    if len(sys.argv) < 2:
        xml_file = os.path.join('..', 'doxygen_doc', 'xml', 'index.xml')
    else:
        xml_file = sys.argv[1]
    # output interface file
    if len(sys.argv) < 3:
        interface_file = os.path.join('..', 'src', 'simuPOP_doc.i')
    else:
        interface_file = sys.argv[2]
    # output ref_single.tex file
    if len(sys.argv) < 4:
        latex_file = os.path.join('..', 'doc', 'simuPOP_ref.tex')
    else:
        latex_file = sys.argv[3]
    # output ref_all.tex file
    if len(sys.argv) < 5:
        latex_testfile = os.path.join('..', 'doc', 'simuPOP_ref_test.tex')
    else:
        latex_testfile = sys.argv[4]
    # read the XML file (actually a index.xml file that contains all others)
    p = Doxy2SWIG(xml_file)
    # generate interface file.
    p.generate()
    # clean up, and process CPPONLY etc
    p.post_process()
    # write interface file to output interface file.
    print 'Writing SWIG interface file to', interface_file
    p.write(interface_file, type='swig')
    print 'Writing latex reference file to', latex_file
    p.write(latex_file, type='latex_single')
    # clear unique name
    p.uniqueName = []
    print 'Writing latex test file to', latex_testfile
    p.write(latex_testfile, type='latex_all', ref_file=latex_file)
    # generating sample document
    os.chdir(os.path.dirname(latex_testfile))
    os.system('pdflatex %s' % os.path.split(latex_testfile)[1])
    # ending statement
    print 'Done.'
