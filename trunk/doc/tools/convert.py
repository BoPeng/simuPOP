# -*- coding: utf-8 -*-
"""
    Convert the simuPOP documentation to Sphinx
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 2007-2008 by Georg Brandl.
    :license: BSD.
"""

import sys
import os, re


import sys
import os
import glob
import shutil
import codecs
from os import path

from converter.tokenizer import Tokenizer
from converter.latexparser import DocParser
from converter.restwriter import RestWriter
from converter.filenamemap import (fn_mapping, copyfiles_mapping, newfiles_mapping,
                          rename_mapping, dirs_to_make, toctree_mapping,
                          amendments_mapping)
from converter.console import red, green

from converter.docnodes import CommentNode, RootNode, NodeList, ParaSepNode, \
     TextNode, EmptyNode, NbspNode, SimpleCmdNode, BreakNode, CommandNode, \
     DescLineCommandNode, InlineNode, IndexNode, SectioningNode, \
     EnvironmentNode, DescEnvironmentNode, TableNode, VerbatimNode, \
     ListNode, ItemizeNode, EnumerateNode, DescriptionNode, \
     DefinitionsNode, ProductionListNode

from converter.util import umlaut, empty, text
from converter.latexparser import ParserError

class MyDocParser(DocParser):

    def __init__(self, *args, **kwargs):
        DocParser.__init__(self, *args, **kwargs)

  
    def mk_metadata_handler(self, name, mdname=None, arg='M'):
        if mdname is None:
            mdname = name
        def handler(self):
            data = self.parse_args('\\'+name, arg)
            self.rootnode.params[mdname] = data[0]
            return EmptyNode()
        return handler

    handle_color = mk_metadata_handler(None, 'color', None, 'M')
    handle_lstset = mk_metadata_handler(None, 'lstset')
    handle_setcounter = mk_metadata_handler(None, 'setcounter', None, 'MM')
    handle_hypersetup = mk_metadata_handler(None, 'hypersetup')
    handle_definecolor = mk_metadata_handler(None, 'definecolor', None, 'MMM')
    handle_sectionfont = mk_metadata_handler(None, 'sectionfont', None, 'O')
    handle_subsectionfont = mk_metadata_handler(None, 'subsectionfont', None, 'O')
    handle_subsubsectionfont = mk_metadata_handler(None, 'subsubsectionfont', None, 'O')
    handle_makeatother = mk_metadata_handler(None, 'makeatother', None, 'O')
    handle_totalheight = mk_metadata_handler(None, 'totalheight', None, 'O')
    handle_columnwidth = mk_metadata_handler(None, 'columnwidth', None, 'O')
    handle_vspace = mk_metadata_handler(None, 'vspace', None, 'M')
    handle_hspace = mk_metadata_handler(None, 'hspace', None, 'M')
    handle_hrule = mk_metadata_handler(None, 'hrule', None, 'O')
    handle_lstlistoflistings = mk_metadata_handler(None, 'lstlistoflistings', None, 'O')
    handle_centering = mk_metadata_handler(None, 'centering', None, 'M')
    handle_textwidth = mk_metadata_handler(None, 'textwidth', None, 'O')
    handle_end = mk_metadata_handler(None, 'end', None, 'O')
    handle_textendash = mk_metadata_handler(None, 'textendash', None, 'O')

    #handle_item = mk_metadata_handler(None, 'item', None, 'O')
    handle_textmd = mk_metadata_handler(None, 'textmd', None, 'O')
    handle_normalsize = mk_metadata_handler(None, 'normalsize', None, 'O')
    handle_textcompwordmark = mk_metadata_handler(None, 'textcompwordmark', None, 'O')
    handle_citep = mk_metadata_handler(None, 'citep', None, 'O')
    handle_citet = mk_metadata_handler(None, 'citet', None, 'O')
    handle_citeyearpar = mk_metadata_handler(None, 'citeyearpar', None, 'O')
    handle_bibliographystyle = mk_metadata_handler(None, 'bibliographystyle', None, 'O')
    handle_bibliography = mk_metadata_handler(None, 'bibliography', None, 'O')
    handle_printindex = mk_metadata_handler(None, 'printindex', None, 'O')

    def handle_minipage_env(self):
        # Ignore the minipage part.
        txt = ''
        while not txt.endswith('end{minipage}'):
            nextl, nextt, nextv, nextr = self.tokens.pop()
            txt += nextv
        return EmptyNode()

    def handle_include(self):
        data = self.parse_args('\\include', 'M')[0].text
        return EmptyNode()

    def handle_newenvironment(self, numOpt=3):
        txt = ''
        opt = 0
        depth = 0
        while True:
            nextl, nextt, nextv, nextr = self.tokens.pop()
            if nextr == '{' or nextr == '[':
                depth += 1
            elif nextr == '}' or nextr == ']':
                depth -= 1
            if nextr == '}' and depth == 0:
                opt += 1
            if opt == numOpt:
                break
        return EmptyNode()

    def handle_newcommand(self):
        return self.handle_newenvironment(2)

    def handle_unrecognized(self, name, line):
        def handler():
            #self.unrecognized.add(name)
            return InlineNode('include', name)
        return handler

        
    handle_small = mk_metadata_handler(None, '\\small', None, 'O')
    handle_ttfamily = mk_metadata_handler(None, '\\small', None, 'O')
    handle_textsf = mk_metadata_handler(None, '\\textsf', None, 'O')
    handle_slshape = mk_metadata_handler(None, '\\small', None, 'O')
    handle_bf = mk_metadata_handler(None, '\\small', None, 'O')
    handle_makeatletter = mk_metadata_handler(None, 'makeatletter', None, 'O')
    handle_lyxline = mk_metadata_handler(None, 'lyxline', None, 'O')
    handle_par = mk_metadata_handler(None, 'par', None, 'O')
    handle_rule = mk_metadata_handler(None, 'rule', None, 'O')
    handle_hfill = mk_metadata_handler(None, 'hfill', None, 'O')
    handle_sloppy = mk_metadata_handler(None, 'sloppy', None, 'O')
    handle_lstlistingname = mk_metadata_handler(None, 'lstlistingname', None, 'O')
    handle_lstlistlistingname = mk_metadata_handler(None, 'lstlistlistingname', None, 'O')


class MyRestWriter(RestWriter):
    def __init__(self, dir = '.', auto_keywords = {}, *args, **kwargs):
        RestWriter.__init__(self, *args, **kwargs)
        self.dirname = dir
        if self.dirname == '':
            self.dirname = '.'
        self.auto_keywords = auto_keywords

    def visit_InlineNode(self, node):
        cmdname = node.cmdname
        if not node.args:
            self.curpar.append(self.simplecmd_mapping[cmdname])
            return
        content = node.args[0]
        if cmdname == 'include':
            file = node.args
            for dir in ['.', self.dirname, 'build']:
                for suffix in ['', '.ref', '.rst', '.txt']:
                    filename = os.path.join(dir, file + suffix)
                    if os.path.isfile(filename):
                        txt = open(filename).read()
                        self.write(txt)
                        return
            print 'Warning: Failed to find included file for filename "%s".' % file
            return
        elif cmdname == 'ref':
            ref_text = text(node.args[0]).lower().replace(':', '-')
            if ref_text.startswith('cha-'):
                self.curpar.append(':ref:`%s%s`' % (self.labelprefix, ref_text))
            else:
                self.curpar.append('`%s%s`_' % (self.labelprefix, ref_text))
            return
        sym = ''
        if cmdname in ('code', 'bfcode', 'samp', 'texttt', 'regexp'):
            sym = '``'
        elif cmdname in ('strong', 'textbf'):
            sym = '**'
        if sym != '':
            cnt = self.get_textonly_node(content, 'code', warn=1)
            done = False
            if isinstance(cnt, TextNode):
                for keyword in self.auto_keywords.keys():
                    txt = text(cnt).split('(')[0]
                    leftover = text(cnt)[len(txt):]
                    if txt in self.auto_keywords[keyword]:
                        self.curpar.append(':%s:`%s`' % (keyword, txt))
                        if leftover != '':
                            if leftover.startswith('()'):
                                # sphinx ignores them!!
                                self.curpar.append("()")
                                if len(leftover) > 2:
                                    self.curpar.append("%s%s%s" % (sym, leftover[2:], sym))
                            elif leftover.startswith('(') and leftover.endswith(')'):
                                self.curpar.append("(%s%s%s)" % (sym, leftover[1:-1], sym))
                            else:
                                self.curpar.append("%s%s%s" % (sym, leftover, sym))
                        return
        # the regular emph/strong case and other stuff
        RestWriter.visit_InlineNode(self, node)

    def visit_CommandNode(self, node):
        cmdname = node.cmdname
        if cmdname == 'figure':
            self.write('.. figure:: %s\n   :width: 670\n' % text(node.args[0]))
            return
        else:
            RestWriter.visit_CommandNode(self, node)

    def visit_CommentNode(self, node):
        # no inline comments -> they are all output at the start of a new paragraph
        pass #self.comments.append(node.comment.strip())

    
def convert_file(infile, outfile, doraise=True, splitchap=True,
                 toctree=None, deflang=None, labelprefix=''):
    inf = codecs.open(infile, 'r', 'latin1')
    p = MyDocParser(Tokenizer(inf.read()).tokenize(), infile)
    if not splitchap:
        outf = codecs.open(outfile, 'w', 'utf-8')
    else:
        outf = None
    refFile = os.path.join(os.path.dirname(infile), 'reflist.py')
    if os.path.isfile(refFile):
        execfile(refFile, globals(), globals())
    r = MyRestWriter(os.path.dirname(infile), auto_keywords, outf, splitchap, toctree, deflang, labelprefix)
    try:
        r.write_document(p.parse())
        if splitchap:
            outf = codecs.open(outfile, 'w', 'utf-8')
            outf.write('.. toctree::\n\n')
            for i, chapter in enumerate(r.chapters[1:]):
                dir = path.dirname(outfile)
                if dir == '':
                    dir = '.'
                filename = '%s/%d_%s' % (dir, i+1, path.basename(outfile))
                outf.write('   %s\n' % filename.lstrip('%s/' % dir))
                coutf = codecs.open(filename, 'w', 'utf-8')
                coutf.write(chapter.getvalue())
                coutf.close()
            outf.close()
        else:
            outf.close()
        p.finish()
        return 1, r.warnings
    except Exception, err:
        if doraise:
            raise
        return 0, str(err)


if __name__ == '__main__':
    convert_file(*sys.argv[1:])
