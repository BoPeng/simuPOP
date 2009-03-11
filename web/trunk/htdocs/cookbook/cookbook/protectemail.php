<?php if (!defined('PmWiki')) exit();
                                                                    /*
  * Copyright *
      Copyright 2005, by Stefan Schimanski (sts@1stein.org).
      Copyright 2004, by Steven Leite (steven_leite@kitimat.net).
      Changes 2004, by Karl Loncarek (Klonk)

  * License *
  Same as PmWiki (GNU GPL)

  * Special Thanks *
  Thanks to Steve Leite and Karl Loncarek for the e-protect code this
  script was based on.

  Thanks to Pm (Patrick Michaud, www.pmichaud.com), for creating
  PmWiki, and for sharing his knowledge and insightfulness with the
  PmWiki community.

  * Description *
  protectemail is an email obfuscation add-on for PmWiki. It rewrites mailto: 
  links during display into some javescript code which hides the address

  While e-protect encodes mail addresses even in the Wiki code of a
  page, this add-on does not. Only while displaying mails are
  hidden. Be aware that this might require some protection of the
  actions edit, source and so on if you want to protect your addresses
  from email grabbers.

  * Installation Instructions *
  1. Copy this script (protectemail.php) to cookbook/
  2. In your config.php file, add the following line:
     include_once('cookbook/protectemail.php');
  3. That's it!

  * History *
  Oct 2 2005 - * written, based on e-protect

  * Configuration *
  There aren't (yet) any configuration variables for this script.    */

//----------------------------------------------------------------------

$RecipeInfo['ProtectEmail']['Version'] = '$Rev$';

## [[mailto:target]]
Markup('protectedmailto','<links',
  "/\\[\\[mailto:([^\\s$UrlExcludeChars]*)\\s*\\]\\]($SuffixPattern)/e",
    "protectEmail('$1','',1)");

## [[nolinkmailto:target]]
Markup('protectednolinkmailto','<links',
  "/\\[\\[nolinkmailto:([^\\s$UrlExcludeChars]*)\\s*\\]\\]($SuffixPattern)/e",
    "protectEmail('$1','',0)");

## [[mailto:target | text]]
Markup('protectedmailto|','<links',
  "/\\[\\[mailto:([^\\s$UrlExcludeChars]*)\\s*\\|\\s*(.*?)\\s*\\]\\]($SuffixPattern)/e",
    "protectEmail('$1','$2',1)");

## [[ text -> mailto:target]]
Markup('-protectedmailto','<links',
  "/\\[\\[(.*?)\\s*-+&gt;\\s*mailto:([^\\s$UrlExcludeChars]*)\\s*\\]\\]($SuffixPattern)/e",
    "protectEmail('$2','$1',1)");

## Add decoding script to Header
$HTMLHeaderFmt['eProtect']= "\n<script type='text/JavaScript'>\n<!--\nNix={map:null,convert:function(a){Nix.init();var s='';for(i=0;i<a.length;i++){var b=a.charAt(i);s+=((b>='A'&&b<='Z')||(b>='a'&&b<='z')?Nix.map[b]:b);}return s;},init:function(){if(Nix.map!=null)return;var map=new Array();var s='abcdefghijklmnopqrstuvwxyz';for(i=0;i<s.length;i++)map[s.charAt(i)]=s.charAt((i+13)%26);for(i=0;i<s.length;i++)map[s.charAt(i).toUpperCase()]=s.charAt((i+13)%26).toUpperCase();Nix.map=map;},decode:function(a){document.write(Nix.convert(a));}}\n//-->\n</script>\n";

//----------------------------------------------------------------------
function encodeStr($s) {
/* str_rot13, extended to recode digits and @#. */
//----------------------------------------------------------------------
  return strtr ($s,
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'nopqrstuvwxyzabcdefghijklmNOPQRSTUVWXYZABCDEFGHIJKLM');
}

//----------------------------------------------------------------------
function protectEmail ($email,$AlternateText,$makeLink) {
//----------------------------------------------------------------------
    $email = encodeStr( $email );
    $html = "";
    if ($makeLink==1)
	if ($AlternateText=='')
	    $html .= "<script type='text/JavaScript'>\n<!--\nNix.decode" .
		"(\"<n pynff='heyyvax' uers='znvygb:$email'>$email</n>\");" . "\n//-->\n</script>";
	else
	    $html .= "<script type='text/JavaScript'>\n<!--\nNix.decode" .
		"(\"<n pynff='heyyvax' uers='znvygb:$email'>\");" . "\n//-->\n</script>" . $AlternateText . 
		"<script type='text/JavaScript'><!--\nNix.decode" .
		"(\"</n>\");" . "\n//-->\n</script>";
    else
	$html .= "<span class=\"nolinkmailto\">" .
	    "<script type='text/JavaScript'>\n<!--\nNix.decode" .
	    "(\"$email\");" . "\n//-->\n</script></span>";
	
    return Keep($html);
}

?>
