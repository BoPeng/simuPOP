<?php
/*  Copyright 2004, 2007 Patrick R. Michaud (pmichaud@pobox.com)

    This script enables the ?skin= and ?setskin= parameters for
    wiki pages.  The ?skin= parameter causes the current page to be 
    displayed with some alternate skin (defined by the WikiAdministrator),
    the ?setskin= parameter sets a cookie on the user's browser to
    always display pages with the requested skin.  

    The easy way to use this script is to do following in 
    local/config.php: 

        $EnableAutoSkinList = 1;
        include_once('cookbook/skinchange.php');

    This will allow browser to choose any skin that exists
    in the pub/skins/ or $FarmD/pub/skins/ directory.

    To only allow a specific set of skins, create an array
    called $PageSkinList that maps names from ?skin= or
    ?setskin= into the corresponding skin to be used for that
    name.  For example:
    
        $PageSkinList = array(
           'pmwiki' => 'pmwiki',
           'myskin' => 'myskin',
           'classic' => 'myclassicskin');
        include_once('cookbook/skinchange.php');

    If $EnableAutoSkinList is also set, then entries in
    $PageSkinList take precedence.  If a browser requests
    a skin that isn't available, then PmWiki defaults to
    the skin already defined by $Skin.

    By default, the setskin cookie that is created will expire after
    one year.  You can set it to expire at the end of the browser
    session by setting $SkinCookieExpires=0;
*/

SDV($RecipeInfo['SkinChange']['Version'], '2007-02-19');
SDV($SkinCookie, $CookiePrefix.'setskin');
SDV($SkinCookieExpires, $Now+60*60*24*365);

if (isset($_COOKIE[$SkinCookie])) $sk = $_COOKIE[$SkinCookie];
if (isset($_GET['setskin'])) {
  $sk = $_GET['setskin'];
  setcookie($SkinCookie, $sk, $SkinCookieExpires, '/');
}
if (isset($_GET['skin'])) $sk = $_GET['skin'];

##  If $EnableAutoSkinList is set, then we accept any skin that
##  exists in pub/skins/ or $FarmD/pub/skins/ .
if (IsEnabled($EnableAutoSkinList, 0) 
    && preg_match('/^[-\\w]+$/', $sk)
    && (is_dir("pub/skins/$sk") || is_dir("$FarmD/pub/skins/$sk")))
  $Skin = $sk;

##  If there's a specific mapping in $PageSkinList, we use it no
##  matter what.
if (@$PageSkinList[$sk]) $Skin = $PageSkinList[$sk];


