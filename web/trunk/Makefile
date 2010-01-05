# backup from webpage
get:
	rsync -v -rlt -z simupop,simupop@web.sourceforge.net:htdocs ${HOME}/simuPOP/web --exclude='src_doc' --exclude='manual_release' --exclude='manual_svn' --exclude='.svn'
	rsync -v -rlt -z simupop,simupop@web.sourceforge.net:/home/persistent/s/si/simupop/web ${HOME}/simuPOP/web/persistent --exclude='sessions' --exclude='.svn'
	rsync -v -rlt -z simupop,simupop@web.sourceforge.net:/home/persistent/s/si/simupop/cookbook ${HOME}/simuPOP/web/persistent --exclude='sessions' --exclude='.svn'

put:
	rsync -v -rlt -z ${HOME}/simuPOP/web/persistent/cookbook/uploads simupop,simupop@web.sourceforge.net:/home/persistent/s/si/simupop/cookbook/ \
	--exclude='.svn' --exclude='.pyc' --exclude='.pyo'

clean:
	rm -f *.png sample* *.dat *.txt *.log  persistent/cookbook/uploads/Cookbook/SimuRecHotSpot/simuRecHotSpot.cfg
	find . -name '*.bak' -exec rm -f {} \;
	rm -rf cdcv simu
	
# figure out what has been added or changed
diff:
	@svn stat | grep -v -E '\.htaccess|\.flock|\.lastmod|\.page|\.stat|RecentChanges|config.php|skins|\.notify|SiteAdmin|\.pyc|\.pyo';\
	if [ "$$?" = "1" ];\
	then\
		echo "No difference";\
	fi

exec:
	for dir in `ls -d persistent/cookbook/uploads/Cookbook/*`;\
	do \
		for file in `ls $$dir/*.py`; \
		do \
			echo;\
			echo;\
			echo;\
			SIMUGUI=False; SIMUDEBUG=DBG_BATCHTESTING; export SIMUGUI SIMUDEBUG; \
			grep BATCHTESTING $$file; \
			if [ "$$?" -ne "0" ]; \
			then \
				echo "Running $$file"; \
				python $$file; \
			else \
				OPTION="`grep BATCHTESTING $$file`"; \
				OPTION=`expr substr "$$OPTION" 15 1000`; \
				if [ "$$OPTION" == "SKIP" ]; \
				then \
					echo; \
					echo; \
				    echo "Skipping testing $$file"; \
					echo; \
				else \
					echo "Running $$file $$OPTION"; \
				    python $$file $$OPTION; \
				fi; \
			fi; \
			if [ "$$?" -ne "0" ]; \
			then \
				echo "Execution of $$file failed. Please modify it."; \
				exit 0; \
			fi; \
		done; \
	done

