# backup from webpage
get:
	rsync -v -rlt -z simupop,simupop@web.sourceforge.net:htdocs ${HOME}/simuPOP/web --exclude='src_doc' --exclude='manual'
	rsync -v -rlt -z simupop,simupop@web.sourceforge.net:/home/persistent/s/si/simupop/web ${HOME}/simuPOP/web/persistent --exclude='sessions'
	rsync -v -rlt -z simupop,simupop@web.sourceforge.net:/home/persistent/s/si/simupop/cookbook ${HOME}/simuPOP/web/persistent --exclude='sessions'

# figure out what has been added or changed
stat:
	@svn stat | grep -v -P '\.htaccess|\.flock|\.lastmod|\.page|\.stat|RecentChanges|config.php|skins|\.notify';\
	if [ "$$?" = "1" ];\
	then\
		echo "No difference";\
	fi

