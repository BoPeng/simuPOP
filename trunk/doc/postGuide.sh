#!/bin/sh

# user's guide
#LYX=/usr/local/bin/lyx14
LYX=/usr/bin/lyx
./convFigRule.pl userGuide.lyx >  userGuidePost.lyx
$LYX --export latex userGuidePost.lyx
mv -f userGuidePost.tex userGuide.tex
perl -pi.bak -e 's/\\IfFileExists{url.sty}{\\usepackage{url}}//' userGuide.tex
perl -pi.bak -e 's/{\\newcommand{\\url}{\\texttt}}//' userGuide.tex

tools/mkhowto --pdf userGuide.tex
/bin/cp -f userGuide.pdf /var/www/html/simuPOP_doc
/bin/rm -f userGuidePost.lyx

mkdir -p /var/www/html/simuPOP_doc/userGuide

latex2html  -show_section_numbers -auto_navigation \
  -local_icons -dir /var/www/html/simuPOP_doc/userGuide  userGuide.tex
/bin/cp -f userGuide.css /var/www/html/simuPOP_doc/userGuide
/bin/cp -f userGuide.tex userGuide_single.tex
latex2html  -split 0 -show_section_numbers -no_navigation \
   userGuide_single.tex
/bin/cp userGuide_single/userGuide_single.html /var/www/html/simuPOP_doc/userGuide_single.html
/bin/cp -f *.png /var/www/html/simuPOP_doc/
perl -pi.bak -e 's,\/home\/bpeng\/research\/simupop\/doc\/\/,,g' /var/www/html/simuPOP_doc/*.html

# reference manual

./convFigRule.pl refManual.lyx >  refManualPost.lyx
$LYX --export latex refManualPost.lyx
/bin/mv -f refManualPost.tex refManual.tex
perl -pi.bak -e 's/\\IfFileExists{url.sty}{\\usepackage{url}}//' refManual.tex
perl -pi.bak -e 's/{\\newcommand{\\url}{\\texttt}}//' refManual.tex

tools/mkhowto --pdf refManual.tex
/bin/cp -f refManual.pdf /var/www/html/simuPOP_doc
rm -f refManualPost.lyx


mkdir -p /var/www/html/simuPOP_doc/refManual
latex2html  -show_section_numbers -auto_navigation \
   -local_icons -dir /var/www/html/simuPOP_doc/refManual  refManual.tex
/bin/cp -f refManual.css /var/www/html/simuPOP_doc/refManual
/bin/cp -f refManual.tex refManual_single.tex
latex2html  -split 0 -show_section_numbers -no_navigation \
   refManual_single.tex
/bin/cp refManual_single/refManual_single.html /var/www/html/simuPOP_doc/refManual_single.html
/bin/cp -f *.png /var/www/html/simuPOP_doc/
perl -pi.bak -e 's,\/home\/bpeng\/research\/simupop\/doc\/\/,,g' /var/www/html/simuPOP_doc/*.html


/bin/cp -f  ../INSTALL /var/www/html/simuPOP_doc/INSTALL.txt

/bin/cp -f ../ChangeLog /var/www/html/simuPOP_doc/ChangeLog
/bin/cp -f ../RoadMap /var/www/html/simuPOP_doc/RoadMap
