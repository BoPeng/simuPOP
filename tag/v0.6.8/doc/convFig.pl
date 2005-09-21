#!/usr/bin/perl


# with this file, I do not have to worry about verbatim float algorithm problem.
# just use in the usual way and later on, chance or do not change using
# this file.
#
undef $/;

$whole = <>;

$whole =~ s/\\begin_inset Include \\verbatiminput{(.*?)\.log}\s*preview false\s*\\end_inset/\\begin_inset Float algorithm\nwide false\ncollapsed false\n\\layout Caption\n\n\1\n\\layout Standard\n\n\n\\begin_inset Include \\verbatiminput{\1\.log}\n\npreview false\n\n\\end_inset \n\n\\end_inset/gsi;

$whole =~ s/\\begin_preamble/\\begin_preamble\n\\usepackage{html}\n\\usepackage{url}\n\\latex{\\renewcommand{\\htmladdnormallink}[2]{#1 (\\url{#2})}}/gsi;
$whole =~ s/\\begin_inset LatexCommand \\url\[(.*?)\]\{(.*?)\}/\\begin_inset ERT\nstatus Collapsed\n\n\\layout Standard\n\n\n\\backslash\nhtmladdnormallink\{\1\}\{\2\}\n/gsi;
$whole =~ s/\\begin_inset LatexCommand \\htmlurl\[(.*?)\]\{(.*?)\}/\\begin_inset ERT\nstatus Collapsed\n\n\\layout Standard\n\n\n\\backslash\nhtmladdnormallink\{\1\}\{\2\}\n/gsi;
print $whole;
