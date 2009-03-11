<?php
/* 

 A bit more documentation of google search boxes can be found here:
 http://www.google.com/support/customsearch/bin/answer.py?hl=en&answer=71640

 Markup:

	(:googlesearchbox:)

	(:googlesearchresults:)

Idea and HTML by Rex C. Eastbourne, implementation by Christian R.

*/

Markup('googlesearchbox', 'directives',
       '/\\(:googlesearchbox(.*?)(.*?)\\s*:\\)/ei',
       "GoogleSearchBox(\$pagename, PSS('$1'), PSS('$2'))");

Markup('googlesearchresults', 'directives',
       '/\\(:googlesearchresults(.*?)(.*?)\\s*:\\)/ei',
       "GoogleSearchResults(\$pagename, PSS('$1'), PSS('$2'))");


function GoogleSearchBox($pagename, $path, $opt) {
  $s = '
<!-- Google CSE Search Box Begins  -->
<form action="" id="cse-search-box">
 <input type="hidden" name="cx" value="008972033992887248866:oarvwprmgvo" />
 <input type="hidden" name="cof" value="FORID:11" />
 <input type="text" name="q" size="25" />
 <input type="submit" name="sa" value="Search" />
</form>
<script type="text/javascript" src="http://www.google.com/coop/cse/brand?form=cse-search-box&lang=en"></script>
<!-- Google CSE Search Box Ends -->
';
  return Keep($s);
}

function GoogleSearchResults($pagename, $path, $opt) {
  $s = '
<!-- Google Search Result Snippet Begins -->
<div id="cse-search-results"></div>
<script type="text/javascript">
 var googleSearchIframeName = "cse-search-results";
 var googleSearchFormName = "cse-search-box";
 var googleSearchFrameWidth = 600;
 var googleSearchDomain = "www.google.com";
 var googleSearchPath = "/cse";
</script>
<script type="text/javascript" src="http://www.google.com/afsonline/show_afs_search.js"></script>

<!-- Google Search Result Snippet Ends -->
';
  return Keep($s);
}