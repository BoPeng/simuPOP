<?php if (!defined('PmWiki')) exit();
/*  Copyright 2004 Patrick R. Michaud (pmichaud@pobox.com)
    This file is distributed under the terms of the GNU General Public 
    License as published by the Free Software Foundation; either 
    version 2 of the License, or (at your option) any later version.  
/*---------------------------------------------------------------

  * Copyright *
  This PmWiki addon was written on March 20, 2004 by Steven Leite
  (steven_leite@kitimat.net).  The same terms and conditions that
  apply to PmWiki also apply to this script.

  * Special Thanks *
  Thanks to  Pm (Patrick Michaud, www.pmichaud.com), for creating
  PmWiki, and for sharing his knowledge  and insightfulness which
  is what made this script possible.

  * Description *
  This add-on allows a PmWiki user to include the contents of any
  external webpage into their WikiPage.

  * Features *
  The script won't overwrite any  existing configuration  date in
  your PmWiki install.

  * Warning *
  Since there's  no filtering of  the imported content (yet), you
  should use caution  if you use this  script.  Users can include
  ANY external page in your Wiki.

  * Installation Instructions *
  1.  Create a new folder: /local/scripts.
  2.  Save this file (x-include.php) in the new directory.
  3.  Add this one line to your config.php file:
      include_once("local/scripts/x-include.php");

  * Usage Instructions *
  Using this script is easy  and intuitive, it uses  the existing
  PmWiki [[Square Bracket]] notation to include an external page.

  Example:  (:includeSite http://www.yahoo.com border=3 scroll=no:)
  or
  Example:  (:includeSite http://www.yahoo.com border=3 scroll=yes:)


  Supported (optional) fields are:

    height = default/pixels/%
    width = default/pixels/%
    border = default/pixels
    scroll = default/yes/no
    align = default/left/right/center

  If the fields  are not specified, the script  defaults  will be
  used.  You can set these defaults in the script (see below).

  * History *
  March 08, 2005 - movement to pmwiki-2
  March 30, 2004 - Slight enhancement to the way the script
                   preserves the url in the <iframe> tag so that


    To use this module, simply place this file in the cookbook/ directory
    and add the following line into config.php:

        include_once('cookbook/includeSite.php');

*/

$DefaultWidth  = '600';        // pixels or %
$DefaultHeight = '400';        // pixels or %
$DefaultAlign  = 'default';    // default, left, right, center
$DefaultScroll = 'default';    // default, yes, no
$DefaultBorder = 'default';    // pixels

function includeSite($str) {
  global $DefaultWidth, $DefaultHeight, $DefaultAlign;
  global $DefaultScroll, $DefaultBorder;

  // parse the [[Double Bracket]] syntax
  $pieces = explode(" ", $str);

  // seperate the pieces in to $key and $value  pairs
  foreach($pieces as $piece) {
    list($key,$value) = explode("=", $piece);
    $vars["$key"] = $value;
   }

  // initialize the <iframe> parameters
  $url         = $vars["url"];
  $width       = $vars["width"];
  $height      = $vars["height"];
  $align       = $vars["align"];
  $scroll      = $vars["scroll"];
  $border      = $vars["border"];


  SDV($url, $pieces[0]);
  SDV($width, $DefaultWidth);
  SDV($height, $DefaultHeight);
  SDV($align, $DefaultAlign);
  SDV($scroll, $DefaultScroll);
  SDV($border, $DefaultBorder);

  $Output .= "\n\n<!-- X-include -->\n\n";
  $Output .= "<iframe width=$width height=$height align=$align frameborder=$border scroll=$scroll src=" . Keep($str) . "></iframe>";
  $Output .= "\n\n<!--/ X-include -->\n\n";
  return $Output;
}

Markup('includeSite', 'directives', "/\\(:includeSite\\s+(http:[^$UrlExcludeChars]*?)\\s*:\\)/e", "includeSite('$1')");


?>
