<?php
/*  
Copyright by Ben Woodruff 2006.  
Modifications by Patrick R. Michaud, 2006.

You may use this code as you want, and change it as much as you want.  
The jsMath package offers many benefits over MimeTeX, such as 
automatic resizing of mathematics.  It also uses standard TeX commands, 
and offers more.  If you install the appropriate fonts, then jsMath
doesn't require downloading pictures at all.

The jsMath package is available from http://www.math.union.edu/~dpvc/jsMath/ .

In order to use this recipe:

1.  Download and install the jsMath package into your PmWiki's 
    pub/ directory, as pub/jsMath/ .  Or, you can install jsMath
    wherever you wish, and set $JsMathUrl to the url of the jsMath
    directory.

2.  Add the following line to a local customization file:

    include_once('cookbook/jsmath.php');

That's it!  The script adds {$$ ... $$} and {$...$} markups that
display LaTeX-style math equations, the first form centers the 
equation, while the second generates the equation "inline".

- You will find your experience to be better if you actually download 
  and install the jsMath fonts, but the site can be viewed just the 
  same without.

- The math graphic for the GUI toolbar is available at
  http://www.pmichaud.com/pmwiki/pub/guiedit/math.gif .

*/

SDV($RecipeInfo['Cookbook.JsMath']['Version'], '20061228');

//  $JSMathUrl contains the url to the jsmath directory on the server.
//  Defaults to pub/jsMath/ .
SDV($JSMathUrl, "$PubDirUrl/jsMath");

// This line gives you LaTeX $$ $$ display equations in the center
Markup('{$$', '<{$',
  '/\\{\\$\\$(.*?)\\$\$\\}/e',
  "Keep('<div class=\"math\">'.PSS('$1').'</div>')");

//  This line gives you $ $ equations in line.  You can then use 
//  \displaystyle as normal to get pretty print equations inline.
Markup('{$', 'directives',
  '/\\{\\$(.*?)\\$\\}/e',
  "Keep('<span class=\"math\">'.PSS('$1').'</span>')");


$HTMLHeaderFmt['jsMath'] = '
  <script> jsMath = {Controls: {cookie: {scale: 120}}} </script>
  <script src="$JSMathUrl/plugins/autoload.js"></script>
';
$HTMLFooterFmt['jsMath'] = '
  <script>
    jsMath.Autoload.Check();
    jsMath.Process(document);
  </script>
';


//  The graphic is available from 
//  http://www.pmichaud.com/pmwiki/pub/guiedit/math.gif .
SDV($GUIButtons['math'],array(1000, '{$ ', ' $}', '\\\\sqrt{n}', 
  '$GUIButtonDirUrlFmt/math.gif"$[Math formula (LaTeX/MimeTeX)]"'));


