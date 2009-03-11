<?php if (!defined('PmWiki')) exit();
/*
 * Beautifier - a source code syntax highlighter for PmWiki
 * Copyright 2003-2004 by Jason Perkins (jason@sim8.com)
 * Beautifier module written by and (C) Mike Jewell
 * Rewritten 2004-12-02 for PmWiki 2.0 by Patrick R. Michaud 
 *    (pmichaud@pobox.com)
 *
 * See http://www.pmichaud.com/wiki/Cookbook/Beautifier for info.
 */

SDV($BEAUT_PATH, dirname(__FILE__)."/php");
SDV($BeautifierOutputType, 'HTML');

function register_beautifier($language)
{
  global $KeepToken;
  Markup("^=beaut_$language", "block", 
    "/^=$language\\s+$KeepToken(\\d.*?)$KeepToken/e",
    "'<:block,1>'.do_beautify('$language', '$1')");
}

function do_beautify($language, $keepindex)
{
  global $BEAUT_PATH, $KPV, $BeautifierOutputType, $HTMLStylesFmt;
  include_once("$BEAUT_PATH/Beautifier/Init.php");
  include_once("$BEAUT_PATH/Output/Output_HTML.php");
  include_once("$BEAUT_PATH/Output/Output_css.php");
  include_once("$BEAUT_PATH/HFile/HFile_$language.php");
  $code = str_replace(array('<:vspace>', '&lt;', '&gt;', '&amp;'), 
                      array('', '<', '>', '&'), $KPV[$keepindex]);
  if ($BeautifierOutputType != 'css') $BeautifierOutputType = 'HTML';
  eval("\$beautifier = new HFile_$language();");
  eval("\$h1 = new Core(\$beautifier, new Output_$BeautifierOutputType());");
  $results = "<pre class='sourcecode'>" 
      . $h1->highlight_text($code) . '</pre>';
  return Keep($results);
}
   
?>
