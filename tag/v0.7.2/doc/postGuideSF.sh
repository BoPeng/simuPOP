#!/bin/sh

scp  /var/www/html/simuPOP_doc/* simupop@shell.sf.net:/home/groups/s/si/simupop/htdocs/simuPOP_doc
scp -r  /var/www/html/simuPOP_doc/userGuide simupop@shell.sf.net:/home/groups/s/si/simupop/htdocs/simuPOP_doc
scp -r /var/www/html/simuPOP_doc/refManual simupop@shell.sf.net:/home/groups/s/si/simupop/htdocs/simuPOP_doc/regManual
