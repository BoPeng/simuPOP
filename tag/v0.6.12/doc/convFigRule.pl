#!/usr/bin/perl


# with this file, I do not have to worry about verbatim float algorithm problem.
# just use in the usual way and later on, chance or do not change using
# this file.
#
undef $/;

$whole = <>;

$whole =~ s/\\begin_inset Include \\verbatiminput{(.*?)\.log}\s*preview false\s*\\end_inset/\\begin_inset ERT\n status Collapsed \n\n\\layout Standard\n\\backslash\nvspace{.5cm}\n\\backslash\nhrule\n\\layout Standard\n\n\\end_inset\n\n \\series bold\nExample: \n\\series default\n\1\n\\layout Standard\n\n\n\\begin_inset ERT\n status Collapsed\n\n\\layout Standard\n\\layout Standard\n\\backslash \nvspace{.1cm} \n\\backslash \nhrule \n\\layout Standard\n\n\\end_inset\n\\begin_inset Include \\verbatiminput{\1\.log}\n\npreview false\n\n\\end_inset \n\n\\begin_inset ERT\n status Collapsed\n \n\\layout Standard\n\n\\layout Standard\n\n\\backslash \nhrule \n\\backslash \nvspace{.5cm} \n\\layout Standard\n\n\\end_inset\n /gsi; 

#$whole =~ s/\\IfFileExists{url.sty}{\\usepackage{url}}\n\s*{\\newcommand{\\url}{\\texttt}}//si;

#$whole =~ s/\\begin_preamble/\\begin_preamble\n\\usepackage{html}\n\\usepackage{url}\n\\latex{\\renewcommand{\\htmladdnormallink}[2]{#1 (\\url{#2})}}/gsi;
#$whole =~ s/\\begin_inset LatexCommand \\url\[(.*?)\]\{(.*?)\}/\\begin_inset ERT\nstatus Collapsed\n\n\\layout Standard\n\n\n\\backslash\nhtmladdnormallink\{\1\}\{\2\}\n/gsi;
#$whole =~ s/\\begin_inset LatexCommand \\htmlurl\[(.*?)\]\{(.*?)\}/\\begin_inset ERT\nstatus Collapsed\n\n\\layout Standard\n\n\n\\backslash\nhtmladdnormallink\{\1\}\{\2\}\n/gsi;
#$whole =~ s/\.eps/.png/gsi;
#$whole =~ s/type=eps/type=png/gsi;


print $whole;