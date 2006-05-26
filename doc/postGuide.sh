#!/bin/sh

# user's guide
#LYX=/usr/local/bin/lyx14
LYX=/usr/local/bin/lyx
$LYX --export latex userGuide.lyx
tools/mkhowto --pdf userGuide.tex
/bin/cp -f userGuide.pdf /var/www/html/simuPOP_doc

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

$LYX --export latex refManual.lyx
tools/mkhowto --pdf refManual.tex
/bin/cp -f refManual.pdf /var/www/html/simuPOP_doc

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
/bin/cp -f ../TODO /var/www/html/simuPOP_doc/RoadMap
