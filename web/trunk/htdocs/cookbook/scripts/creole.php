<?php if (!defined('PmWiki')) exit();
/*  Copyright 2007 Patrick R. Michaud (pmichaud@pobox.com)
    This file is part of PmWiki; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  See pmwiki.php for full details.

    This script adds Creole v0.4 markup (http://www.wikicreole.org/)
    to PmWiki.  To activate this script, simply add the following into
    a local customization file:

        include_once('scripts/creole.php');

*/

## **strong**
Markup('**', 'inline', 
  '/^\\*\\*(?>(.+?)\\*\\*)(?!\\S)|(?<!^)\\*\\*(.+?)\\*\\*/',
  '<strong>$1$2</strong>');

## //emphasized//
Markup('//', 'inline', 
  '/(?<!http:|ftp:)\\/\\/(.*?)\\/\\//',
  '<em>$1</em>');

## == Headings ==
Markup('^=', 'block',
  '/^(={1,6})\\s?(.*?)(\\s*=*\\s*)$/e',
  "'<:block,1><h'.strlen('$1').PSS('>$2</h').strlen('$1').'>'");

## Line breaks
Markup('\\\\', 'inline', '/\\\\\\\\/', '<br />');

## Preformatted
Markup('^{{{', '[=',
  "/^\\{\\{\\{\n(.*?\n)\\}\\}\\}[^\\S\n]*\n/sme",
  "Keep(PSS('<pre class=\"escaped\">$1</pre>'))");
Markup('{{{', '>{{{',
  '/\\{\\{\\{(.*?)\\}\\}\\}/se',
  "Keep(PSS('<code class=\"escaped\">$1</code>'))");

## Tables
Markup('|-table', '>^||',
  '/^\\|(.*)$/e',
  "FormatTableRow(PSS('$0'), '\\|')");

## Images
Markup('{{', 'inline',
  '/\\{\\{(?>(\\L))([^|\\]]*)(?:\\|\\s*(.*?)\\s*)?\\}\\}/e',
  "Keep(\$GLOBALS['LinkFunctions']['$1'](\$pagename, '$1', '$2', '$3',
     '$1$2', \$GLOBALS['ImgTagFmt']),'L')");


## GUIButtons
SDVA($GUIButtons, array(
  'em'       => array(100, "//", "//", '$[Emphasized]',
                  '$GUIButtonDirUrlFmt/em.gif"$[Emphasized (italic)]"',
                  '$[ak_em]'),
  'strong'   => array(110, "**", "**", '$[Strong]',
                  '$GUIButtonDirUrlFmt/strong.gif"$[Strong (bold)]"',
                  '$[ak_strong]'),
  'h2'       => array(400, '\\n== ', ' ==\\n', '$[Heading]',
                  '$GUIButtonDirUrlFmt/h.gif"$[Heading]"'),

  ));

