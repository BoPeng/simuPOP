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
    #def handle_ifhtml1(self):
    #    txt = ''
    #    while not txt.endswith('\\fi'):
    #        nextl, nextt, nextv, nextr = self.tokens.pop()
    #        txt += nextr
    #    return EmptyNode()

    def __init__(self, *args, **kwargs):
        DocParser.__init__(self, *args, **kwargs)
    #    self.handle_ifhtml = self.handle_ifhtml1
    #
    def handle_unrecognized(self, name):
        def handler():
            print 'Unrecognized name: %s, use include directive to try to include an external file' % name
            return InlineNode('include', name)
        return handler

    def parse_until(self, condition=None, endatbrace=False):
        nodelist = NodeList()
        bracelevel = 0
        for l, t, v, r in self.tokens:
            if condition and condition(t, v, bracelevel):
                return nodelist.flatten()
            if t == 'command':
                if len(v) == 1 and not v.isalpha():
                    nodelist.append(self.handle_special_command(v))
                    continue
                handler = getattr(self, 'handle_' + v, None)
                if not handler:
                    handler = self.handle_unrecognized(v)
                nodelist.append(handler())
            elif t == 'bgroup':
                bracelevel += 1
            elif t == 'egroup':
                if bracelevel == 0 and endatbrace:
                    return nodelist.flatten()
                bracelevel -= 1
            elif t == 'comment':
                nodelist.append(CommentNode(v))
            elif t == 'tilde':
                nodelist.append(NbspNode())
            elif t == 'mathmode':
                pass # ignore math mode
            elif t == 'parasep':
                nodelist.append(ParaSepNode())
            else:
                # includes 'boptional' and 'eoptional' which don't have a
                # special meaning in text
                nodelist.append(TextNode(v))
        return nodelist.flatten()
 
    def parse_args(self, cmdname, argspec):
        """ Helper to parse arguments of a command. """
        # argspec: M = mandatory, T = mandatory, check text-only,
        #          O = optional, Q = optional, check text-only
        args = []
        # ignore \[ \] for now
        
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


    def handle_special_command(self, cmdname):
        if cmdname == '[':
            print 'MATH START'
            while True:
                nextl, nextt, nextv, nextr = self.tokens.next()
                print 'MATH %s %s' % (nextr, nextt)
                if nextr == '\\]':
                    break
            print 'MATH END'
            return EmptyNode()
        return DocParser.handle_special_command(self, cmdname)

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
    handle_centering = mk_metadata_handler(None, 'centering', None, 'O')
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

    handle_ne = mk_metadata_handler(None, 'ne', None, 'O')
    handle_sum = mk_metadata_handler(None, 'sum', None, 'O')
    handle_mu = mk_metadata_handler(None, 'mu', None, 'O')
    handle_py = mk_metadata_handler(None, 'py', None, 'O')
    handle_infty = mk_metadata_handler(None, 'infty', None, 'O')
    handle_gg = mk_metadata_handler(None, 'gg', None, 'O')
    handle_left = mk_metadata_handler(None, 'left', None, 'O')
    handle_right = mk_metadata_handler(None, 'right', None, 'O')
    handle_times = mk_metadata_handler(None, 'times', None, 'O')
    handle_rightarrow = mk_metadata_handler(None, 'rightarrow', None, 'O')

    def handle_textquoteright(self):
        data = self.parse_args('\\textquoteright', 'M')
        return TextNode("'")

    def handle_textquoteleft(self):
        data = self.parse_args('\\textquoteleft', 'M')
        return TextNode("'")

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
    
    #def handle_renewcommand(self):
    #    return self.handle_newenvironment(2)

    #def handle_title(self):
    #    data = self.parse_args('\\title', 'M')
    #    return SectioningNode('chapter', data)

    def handle_lstinputlisting(self):
        txt = ''
        while not txt.endswith('\n'):
            nextl, nextt, nextv, nextr = self.tokens.pop()
            txt += nextr
        try:
            label = re.search('label={([^}]*)}', txt).groups()[0]
            label = label.split(':')[-1]
        except:
            label = ''
        try:
            caption = re.search('caption={([^}]*)}', txt).groups()[0]
            caption = caption.replace(label, '')
        except:
            caption = ''
        try:
            listing = re.search('\]{([^{]*)}', txt).groups()[0]
        except:
            listing = ''
        # do not use this one yet
        #return InlineNode('listing', '%s\n%s\n%s' % (listing, label, caption))
        try:
            src = open(listing)
        except:
            print("File ", listing, " can not be opened")
            return EmptyNode()
        txt = src.read()
        return NodeList([CommandNode('label', [TextNode(label)]),
                InlineNode('strong', [TextNode('Example')]),
                TextNode(': ' + caption),
                VerbatimNode(TextNode(txt))])


    def handle_figure_env(self):
        txt = ''
        while not txt.endswith('end{figure}'):
            nextl, nextt, nextv, nextr = self.tokens.pop()
            txt += nextr
        print re.match('label', txt)
        try:
            label = re.search('\\label{([^}]*)}', txt).groups()[0]
            label = label.split(':')[-1]
        except:
            label = ''
        try:
            txt = txt.replace('\\label{fig:%s}' % label, '')
            print 'LABEL', label
            print txt
            caption = re.search('\\caption{([^}]*)}', txt).groups()[0]
        except:
            caption = ''
        try:
            figure = re.search('\\includegraphics(\[.*\]){([^}]*)}', txt).groups()[1]
        except:
            figure = ''
        try:
            legend = re.search('\\includegraphics[^\n]*([^\z]*)', txt, re.M).groups()[0]
            legendlist = legend.split()
            legend = ''
            for line in legendlist:
                if not line.startswith('\\end{'):
                    legend += line + ' '
        except:
            legend = ''
        return NodeList([CommandNode('label', [TextNode(label)]),
                InlineNode('strong', [TextNode('Figure')]),
                TextNode(': '),
                InlineNode('emph', [TextNode(caption.strip())]),
                #BreakNode(), 
                InlineNode('figure', figure),
                TextNode('\n\n' + legend)])

    def handle_lstlisting_env(self):
        txt = ''
        opt = False
        optText = ''
        while not txt.endswith('end{lstlisting}'):
            nextl, nextt, nextv, nextr = self.tokens.pop()
            if nextr == '[':
                opt = True
            elif opt and nextr == ']':
                opt = False
                continue
            if opt:
                optText += nextr
            else:
                txt += nextr
        # find label
        try:
            label = re.search('label=\{([^}]*)}', optText).groups()[0]
            label = label.split(':')[-1]
        except:
            label = ''
        try:
            caption = re.search('caption=\{([^}]*)}', optText).groups()[0]
            caption = caption.split(':')[-1]
        except:
            caption = ''
        txtlist = txt.split('\n')
        txt = ''
        for t in txtlist:
            if t.startswith('[cap') or '\\end{lstlisting}' in t:
                continue
            txt += t + '\n'
        if label != '':
            return NodeList([CommandNode('label', [TextNode(label)]),
                InlineNode('strong', [TextNode('Example')]),
                TextNode(': ' + caption),
                VerbatimNode(TextNode(txt))])
        return VerbatimNode(TextNode(txt))


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
            figure = node.args
            for suffix in ['', '.pdf', '.png', '.jpg', '.eps']:
                file = os.path.join('..', figure + suffix)
                if os.path.isfile(file):
                    figure = file
                    break
            self.write('.. image:: %s\n' % figure )
        elif cmdname == 'listing':
            content = node.args[0]
            #self.flush_par()
            txt = node.args.split('\n')
            file = txt[0]
            label = txt[1]
            caption = txt[2]
            self.write('.. literalinclude:: %s' % file)
            self.write('   :language: python')
            self.write('   %s\n' % caption)
        elif cmdname == 'include':
            file = node.args
            for suffix in ['', '.ref', '.rst', '.txt']:
                if os.path.isfile(file + suffix):
                    file += suffix
                    txt = open(file).read()
                    self.write(txt)
                    break
        elif cmdname == 'ref':
            self.curpar.append('`%s%s`_' % (self.labelprefix,
                                                text(node.args[0]).lower().split(':')[-1]))
        else:
            RestWriter.visit_InlineNode(self, node)

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
    r = MyRestWriter(outf, splitchap, toctree, deflang, labelprefix)
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
        return 1, r.warnings
    except Exception, err:
        if doraise:
            raise
        return 0, str(err)


if __name__ == '__main__':
    convert_file(*sys.argv[1:])
