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

from converter.util import umlaut, empty

class MyDocParser(DocParser):
    def __init__(self, *args, **kwargs):
        DocParser.__init__(self, *args, **kwargs)
    #
    
    def parse_args(self, cmdname, argspec):
        """ Helper to parse arguments of a command. """
        # argspec: M = mandatory, T = mandatory, check text-only,
        #          O = optional, Q = optional, check text-only
        args = []
        def optional_end(type, value, bracelevel):
            return type == 'eoptional' and bracelevel == 0

        for i, c in enumerate(argspec):
            assert c in 'OMTQ'
            nextl, nextt, nextv, nextr = self.tokens.pop()
            while nextt == 'comment' or (nextt == 'text' and nextv.isspace()):
                nextl, nextt, nextv, nextr = self.tokens.pop()

            if c in 'OQ':
                if nextt == 'boptional':
                    arg = self.parse_until(optional_end)
                    if c == 'Q' and not isinstance(arg, TextNode):
                        raise ParserError('%s: argument %d must be text only' %
                                          (cmdname, i), nextl)
                    args.append(arg)
                else:
                    # not given
                    args.append(EmptyNode())
                    self.tokens.push((nextl, nextt, nextv, nextr))
                continue

            if nextt == 'bgroup':
                arg = self.parse_until(None, endatbrace=True)
                if c == 'T' and not isinstance(arg, TextNode):
                    raise ParserError('%s: argument %d must be text only' %
                                      (cmdname, i), nextl)
                args.append(arg)
            else:
                if nextt == 'command':
                    nextl, nextt, nextv, nextr = self.tokens.pop()
                elif nextt != 'text':
                    raise ParserError('%s: non-grouped non-text arguments not '
                                      'supported' % cmdname, nextl)
                args.append(TextNode(nextv[0]))
                self.tokens.push((nextl, nextt, nextv[1:], nextr[1:]))
        return args

    def mk_metadata_handler(self, name, mdname=None, arg='M'):
        if mdname is None:
            mdname = name
        def handler(self):
            data = self.parse_args('\\'+name, arg)
            self.rootnode.params[mdname] = data[0]
            return EmptyNode()
        return handler

    handle_lstset = mk_metadata_handler(None, 'lstset')
    handle_color = mk_metadata_handler(None, 'color', None, 'M')
    handle_setcounter = mk_metadata_handler(None, 'setcounter')
    handle_hypersetup = mk_metadata_handler(None, 'hypersetup')
    handle_newcommand = mk_metadata_handler(None, 'newcommand')
    handle_renewcommand = mk_metadata_handler(None, 'renewcommand', None, 'QMM')
    handle_definecolor = mk_metadata_handler(None, 'definecolor', None, 'MMM')
    handle_newenvironment = mk_metadata_handler(None, 'newenvironment', None, 'QTT')
    handle_sectionfont = mk_metadata_handler(None, 'sectionfont', None, 'O')
    handle_subsectionfont = mk_metadata_handler(None, 'subsectionfont', None, 'O')
    handle_subsubsectionfont = mk_metadata_handler(None, 'subsubsectionfont', None, 'O')
    handle_makeatother = mk_metadata_handler(None, 'makeatother', None, 'O')
    handle_totalheight = mk_metadata_handler(None, 'totalheight', None, 'O')
    handle_columnwidth = mk_metadata_handler(None, 'columnwidth', None, 'O')
    handle_vspace = mk_metadata_handler(None, 'vspace', None, 'M')
    handle_hrule = mk_metadata_handler(None, 'hrule', None, 'O')
    handle_lstlistoflistings = mk_metadata_handler(None, 'lstlistoflistings', None, 'O')
    handle_caption = mk_metadata_handler(None, 'caption', None, 'O')
    handle_centering = mk_metadata_handler(None, 'centering', None, 'O')
    handle_includegraphics = mk_metadata_handler(None, 'includegraphics', None, 'O')
    handle_textwidth = mk_metadata_handler(None, 'textwidth', None, 'O')
    handle_end = mk_metadata_handler(None, 'end', None, 'O')
    handle_textendash = mk_metadata_handler(None, 'textendash', None, 'O')
    handle_textquoteleft = mk_metadata_handler(None, 'textquoteleft', None, 'O')
    handle_item = mk_metadata_handler(None, 'item', None, 'O')
    handle_textmd = mk_metadata_handler(None, 'textmd', None, 'O')
    handle_normalsize = mk_metadata_handler(None, 'normalsize', None, 'O')
    handle_textcompwordmark = mk_metadata_handler(None, 'textcompwordmark', None, 'O')
    handle_citep = mk_metadata_handler(None, 'citep', None, 'O')
    handle_citet = mk_metadata_handler(None, 'citet', None, 'O')
    handle_citeyearpar = mk_metadata_handler(None, 'citeyearpar', None, 'O')
    handle_bibliographystyle = mk_metadata_handler(None, 'bibliographystyle', None, 'O')
    handle_bibliography = mk_metadata_handler(None, 'bibliography', None, 'O')
    handle_printindex = mk_metadata_handler(None, 'printindex', None, 'O')

    handle_ne = mk_metadata_handler(None, 'ne', None, 'O')
    handle_sum = mk_metadata_handler(None, 'sum', None, 'O')
    handle_mu = mk_metadata_handler(None, 'mu', None, 'O')
    handle_infty = mk_metadata_handler(None, 'infty', None, 'O')
    handle_gg = mk_metadata_handler(None, 'gg', None, 'O')
    handle_left = mk_metadata_handler(None, 'left', None, 'O')
    handle_right = mk_metadata_handler(None, 'right', None, 'O')
    handle_times = mk_metadata_handler(None, 'times', None, 'O')
    handle_rightarrow = mk_metadata_handler(None, 'rightarrow', None, 'O')

    def handle_textquoteright(self):
        data = self.parse_args('\\textquoteright', 'M')
        return TextNode("'")

    def handle_minipage_env(self):
        # Ignore the minipage part.
        text = ''
        while not text.endswith('end{minipage}'):
            nextl, nextt, nextv, nextr = self.tokens.pop()
            text += nextv
        return EmptyNode()

    def handle_lyxcode_env(self):
        return VerbatimNode(TextNode('a'))

    def handle_lstinputlisting(self):
        text = ''
        while not text.endswith('\n'):
            nextl, nextt, nextv, nextr = self.tokens.pop()
            text += nextr
        try:
            label = re.search('label={([^}]*)}', text).groups()[0]
        except:
            label = ''
        try:
            caption = re.search('caption={([^}]*)}', text).groups()[0]
            caption = caption.replace(label, '')
        except:
            caption = ''
        try:
            listing = re.search('\]{([^{]*)}', text).groups()[0]
        except:
            listing = ''
        # do not use this one yet
        #return InlineNode('listing', '%s\n%s\n%s' % (listing, label, caption))
        try:
            src = open(listing)
        except:
            print("File ", listing, " can not be opened")
            return EmptyNode()
        text = src.read()
        if label != '':
            return NodeList([CommandNode('label', [TextNode(label)]),
                VerbatimNode(TextNode(text))])
        return VerbatimNode(TextNode(text))


    def handle_figure_env(self):
        text = ''
        while not text.endswith('end{figure}'):
            nextl, nextt, nextv, nextr = self.tokens.pop()
            text += nextr
        print re.match('label', text)
        try:
            label = re.search('\\label\{([^}]*)}', text).groups()[0]
        except:
            label = ''
        try:
            caption = re.search('\\caption{(.*)}', text).groups()[0]
            caption = caption.replace('\\label{%s}' % label, '')
        except:
            caption = ''
        try:
            figure = re.search('\\includegraphics(\[.*\]){([^}]*)}', text).groups()[1]
        except:
            figure = ''
        try:
            legend = re.search('\\includegraphics[^\n]*([^\z]*)', text, re.M).groups()[0]
            legendlist = legend.split()
            legend = ''
            for line in legendlist:
                if not line.startswith('end{'):
                    legend += line + ' '
        except:
            legend = ''
        if label != '':
            return NodeList([CommandNode('label', [TextNode(label)]),
                InlineNode('figure', '%s\n%s\n%s\n%s' % (figure, label, caption, legend))])
        else:
            return InlineNode('figure', '%s\n%s\n%s\n%s' % (figure, label, caption, legend))

    def handle_lstlisting_env(self):
        text = ''
        opt = False
        while not text.endswith('end{lstlisting}'):
            nextl, nextt, nextv, nextr = self.tokens.pop()
            if nextr == '[':
                opt = True
            elif opt and nextr == ']':
                opt = False
                continue
            if not opt:
                text += nextr
        textlist = text.split('\n')
        text = ''
        for t in textlist:
            if t.startswith('[cap') or '\\end{lstlisting}' in t:
                continue
            text += t + '\n'
        return VerbatimNode(TextNode(text))


    def handle_centering_env(self):
        return VerbatimNode(TextNode('a'))

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
    def __init__(self, *args, **kwargs):
        RestWriter.__init__(self, *args, **kwargs)

    def visit_InlineNode(self, node):
        cmdname = node.cmdname
        if cmdname == 'figure':
            content = node.args[0]
            #self.flush_par()
            text = node.args.split('\n')
            figure = text[0]
            label = text[1]
            caption = text[2]
            legend = text[3]
            self.write('.. figure:: %s' % figure)
            self.write('   %s' % caption)
            self.write('')
            self.write(legend)
        elif cmdname == 'listing':
            content = node.args[0]
            #self.flush_par()
            text = node.args.split('\n')
            file = text[0]
            label = text[1]
            caption = text[2]
            self.write('.. literalinclude:: %s' % file)
            self.write('   :language: python')
            self.write('   %s\n' % caption)
        else:
            RestWriter.visit_InlineNode(self, node)


def convert_file(infile, outfile, doraise=True, splitchap=False,
                 toctree=None, deflang=None, labelprefix=''):
    inf = codecs.open(infile, 'r', 'latin1')
    p = MyDocParser(Tokenizer(inf.read()).tokenize(), infile)
    if not splitchap:
        outf = codecs.open(outfile, 'w', 'utf-8')
    else:
        outf = None
    r = MyRestWriter(outf, splitchap, toctree, deflang, labelprefix)
    try:
        r.write_document(p.parse())
        if splitchap:
            for i, chapter in enumerate(r.chapters[1:]):
                coutf = codecs.open('%s/%d_%s' % (
                    path.dirname(outfile), i+1, path.basename(outfile)),
                                    'w', 'utf-8')
                coutf.write(chapter.getvalue())
                coutf.close()
        else:
            outf.close()
        return 1, r.warnings
    except Exception, err:
        if doraise:
            raise
        return 0, str(err)


if __name__ == '__main__':
    convert_file('userGuide.tex', 'userGuide.rst')
