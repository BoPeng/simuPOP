#!/usr/bin/perl


# with this file, I do not have to worry about verbatim float algorithm problem.
# just use in the usual way and later on, chance or do not change using
# this file.
#
undef $/;

$whole = <>;

$whole =~ s/\\begin_inset Include \\verbatiminput{(.*?)\.log}\s*preview false\s*\\end_inset/\\begin_inset ERT\n status Collapsed \n\n\\layout Standard\n\\backslash\nvspace{.5cm}\n\\backslash\nhrule\n\\layout Standard\n\n\\end_inset\n\n \\series bold\nExample: \n\\series default\n\1\n\\layout Standard\n\n\n\\begin_inset ERT\n status Collapsed\n\n\\layout Standard\n\\layout Standard\n\\backslash \nvspace{.1cm} \n\\backslash \nhrule \n\\layout Standard\n\n\\end_inset\n\\begin_inset Include \\verbatiminput{\1\.log}\n\npreview false\n\n\\end_inset \n\n\\begin_inset ERT\n status Collapsed\n \n\\layout Standard\n\n\\layout Standard\n\n\\backslash \nhrule \n\\backslash \nvspace{.5cm} \n\\layout Standard\n\n\\end_inset\n /gsi; 
print $whole;
