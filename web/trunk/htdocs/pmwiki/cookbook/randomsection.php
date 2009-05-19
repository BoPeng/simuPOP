<?php if (!defined('PmWiki')) exit();
/*  Copyright 2006 Patrick R. Michaud (pmichaud@pobox.com)
    This file is randomsection.php; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  

    To install this script, simply add the line below to a local
    customization file:

        include_once('cookbook/randomsection.php');

    This script adds a {$RandomSection} page variable markup to PmWiki.
    Its primary use is to implement RandomQuotes -- i.e., to grab
    random quotes from other pages.  The typical usage is inside
    an (:include:) directive, as in:

        (:include {Main.Quotes$RandomSection}:)

    This grabs a random section from Main.Quotes and includes
    it within the current text.  Sections are denoted by anchors --
    i.e.: the Main.Quotes page could look like:

        [[#first]]
        This is a quote.

        [[#second]]
        Another quote

        [[#third]]
        A third quote.
 
*/

SDV($RecipeInfo['RandomQuote']['Version'], '2007-02-15');

$FmtPV['$RandomSection'] = 'RandomSection($pn)';

function RandomSection($pagename) {
  static $anchors;
  $npat = '[[:alpha:]][-\\w]*';
  $page = RetrieveAuthPage($pagename, 'read', false, READPAGE_CURRENT);
  $a = '';
  if ($page 
      && preg_match_all("/\\[\\[(#$npat)\\]\\]/", $page['text'], $alist)) {
    while ($alist[1]) {
      list($a) = array_splice($alist[1], rand(0, count($alist[1])-1), 1);
      if (@$anchors[$a]++ == 0) break;
    }
  }
  return "$pagename$a";
}

