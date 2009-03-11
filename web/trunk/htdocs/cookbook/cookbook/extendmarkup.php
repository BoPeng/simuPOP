<?php if (!defined('PmWiki')) exit();
/*
# This adds a range of markup extensions to PmWiki 2. It combines into
# a single script a number of separate items from PmWiki 1 and adds more

    Version 2.0.59a

    Copyright 2004-2006 John Rankin (john.rankin@affinity.co.nz)
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    
    lazy web links Copyright 2004 Patrick R. Michaud (pmichaud@pobox.com)  
*/

SDV($MarkupCss,false);
if ($MarkupCss) $HTMLHeaderFmt[] = 
 "<link rel='stylesheet' href='\$FarmPubDirUrl/css/markup.css' type='text/css' />\n";
else $HTMLStylesFmt['extend'] = "
a.createlink { color: red; }
#wikitext { line-height: 1.4em; }
#wikitext sub, #wikitext sup { line-height: 0; }
div.footnote { 
    width: 160px; 
    border-bottom: 1px solid blue;
	margin-bottom: 0.5em;
}
p.footnote {
	text-indent: -1em;
	margin-right: 3em;
	margin-left: 3em;
	margin-top: 0px;
	margin-bottom: 0.5em;
	font-size: smaller;
}
p.qanda:first-letter {
    float: left;
    font-family: Old English, Georgia, serif;
    color: #777777;
    font-size: 250%;
    line-height: 0.85em;
    margin-right: 0.3em;
    margin-bottom:-0.25em;
}
p.drop:first-letter {
    float: left;
    font-family: Old English, Georgia, serif;
    font-size: 290%;
    line-height: 0.85em;
    margin-right: 0.1em;
    margin-bottom:-0.25em;
}
del { color: red; }
ins { background-color: yellow; }
div.inote {
    font-size: 10px;
    line-height: 1.2em;
    float: right;
    padding: 2px;
    margin-left: 10px;
    margin-bottom: 10px;
    width: 200px;
    border-top: 1px dotted gray;
    border-bottom: 1px dotted gray;
    background-color: #ffffa1;
}
div.inote h1 {
    background-color: #ffe53e;
    font-size: 10px;
    font-weight: normal;
    margin-top: 0px;
    padding-bottom: 1px;
    margin-bottom: 3px;
}
div.inote h1 span.inote {
    float: right;
}
div.inote ul, div.inote ol { margin-left: -1.5em; }
div.inote p.vspace { margin-top:0.5em; }
span.stickynote {
    font-size: smaller;
    float: right;
    padding: 8px;
    margin-left: 10px;
    margin-bottom: 10px;
    width: 175px;
    border-top: 2px solid gray;
    border-bottom: 2px solid gray;
    text-align: center;
    color: navy;
    background-color: #dcdcdc;
}
span.smallcaps { font-variant: small-caps; }
dfn  { font-style: normal; cursor: help; }
abbr { font-style: italic; cursor: help; }
abbr, dfn.definition { border-bottom: 1px dotted; }
h5.runin { display: run-in; font-size: 100%; border: none; }
div.figure {
    border: thin silver solid;
    padding: 0.3em;
}
div.figure p {
    text-align: center;
    font-style: italic;
    font-size: smaller;
    padding-top: 0.2em;
    margin: 0;
}
dd, li p { margin-bottom: 0.5em }
b.selflink { border-bottom: 1px dotted; }
@media screen{ #wikitext b.selflink { color: #e66e31; } }
";

SDV($MarkupExtensionsFmt,
    array("inote abbr `A `. `- `s `: `f -d ... aquo mac '/ '@ '; [^ copy",
    "q&a A; {|} =| {= revisions ^!! fig :: para lazyweb spaced squo links"));
foreach(explode(' ',implode(' ',$MarkupExtensionsFmt)) as $me)
    SDV($MarkupExtensions[$me], true);

## (:group:) -- removed as not compatible with PmWiki 2.1
/*
if ($MarkupExtensions['group']) {
    Markup('group','directives','/\\(:group\\s(.*?):\\)/ei',
      "PZZ(\$GLOBALS['PCache'][\$pagename]['group']=PSS('$1'))");
    $FmtP['/\\$Groupspaced/e'] = 
  '(@$PCache[$pagename]["group"]) ? $PCache[$pagename]["group"] : $AsSpacedFunction(@$match[1])';
}
*/

if ($MarkupExtensions['inote']) {
    SDV($InoteTextFmt, "[[\$FullName|\$Title]]" . 
        " <span class='inote'>([[\$FullName?action=edit|edit]])</span>");
    SDV($InoteExpiredFmt, "Expired \$LastModified");
    if ($action=="print" || $action=="publish")
        Markup('inote','>if',"/\\(:inote\\s+.*?:\\)/",'');
    else 
        Markup('inote','>if',
    "/\(:inote\s+(?:days=(\d+)\s+)?((?:$GroupPattern(?:[\/.]))?$NamePattern)(.*?):\)/e",
    "PRR().IncludeNoteText(\$pagename,array('days' => '$1', 'fmt' => \$GLOBALS['InoteExpiredFmt']),'$2','$3',\$GLOBALS['InoteTextFmt'])");
}

function IncludeNoteText($pagename,$t,$page,$opts,$fmt) {
  global $Now,$PCache;
  $age = $Now - 86400 * (($t['days']) ? $t['days'] : 365);
  $p = MakePageName($pagename,$page);
  PCache($p,RetrieveAuthPage($p, 'read', false, READPAGE_CURRENT));
  $i = ($PCache[$p]['time'] >= $age) ? preg_replace('/\(:title\s+.*?:\)/', '',
        IncludeText($pagename,"include $p$opts")) : FmtPageName($t['fmt'],$p);
  return "<div class='inote'>\n!".FmtPageName($fmt,$p)."\n$i\n<:block></div>";
}

/*
## require 2 lower and upper case characters for a WikiWord
$WikiWordPattern = '[[:upper:]][[:alnum:]]*(?:[[:upper:]][[:lower:]0-9][[:lower:]0-9]|[[:lower:]0-9][[:lower:]0-9][[:upper:]]|[[:lower:]0-9][[:upper:]]+[[:lower:]0-9])[[:alnum:]]*';
*/
## prevent wikiwords with only one lower case letter
SDV($AbbreviationPattern,
  "[[:upper:]]+(?:[[:upper:]][[:lower:]0-9]|[[:lower:]0-9][[:upper:]])[[:upper:]]*");
if ($MarkupExtensions['abbr']) {
    $AbbreviationEnabled = true;
    Markup("abbr",'<`wikiword',"/`?\\b($AbbreviationPattern)\\b/e",
      "Keep(((PageExists(MakePageName(\$pagename,'$1'))) ? MakeLink(\$pagename, '$1', '$1') : '$1'), 'L')");
#      "Keep('$1')");
    Markup("mc",'<`wikiword',"/`?\\b(Ma?c[[:upper:]][[:lower:]]+)\\b/e",
      "Keep('$1')");
}

#### escape character (backtick) ####
## prevent WikiWords with Wiki`Word  and `WikiWord markups
if ($MarkupExtensions['`A'])
    Markup("`A",'>links','/([[:alnum:].\/])?`([[:upper:]])/','$1$2');

## '`.' (invisible stop)
if ($MarkupExtensions['`.'])
    Markup("`.",'>links',"/`\./",'');

## '`-' (en dash) 
if ($MarkupExtensions['`-'])
    Markup("`-",'inline',"/`-/",'&ndash;');

## '` ' (nonbreaking space)
if ($MarkupExtensions['`s'])
    Markup("`s",'inline',"/`\s/",'&nbsp;');

## '`:' (middot)
if ($MarkupExtensions['`:'])
    Markup("`:",'inline',"/`:/",'&middot;');

## simple fractions (quarter, half, three quarters)
if ($MarkupExtensions['`f']) {
    Markup("1/4",'inline',"/`1\/?4/",'&#188;');
    Markup("1/2",'inline',"/`1\/?2/",'&#189;');
    Markup("3/4",'inline',"/`3\/?4/",'&#190;');
}

## em dash, en dash, plus or minus, and minus
if ($MarkupExtensions['-d']) {
    Markup("--",'>[+',"/(^|[^!-;])--([^-&>]|$)/",'$1&mdash;$2');
    Markup("d-d",'>links',"/(\\d)-(\\d)/",'$1&ndash;$2');
    Markup("dxd",'>links',"/(\\d)x(\\d)/",'$1&times;$2');
    Markup("+-",'<-d',"/\+\/?-/",'&plusmn;');
    Markup("-d",'>d-d',"/([^'\"[:alpha:]])-(\\d)/",'$1&minus;$2');
}

## ellipsis ...
if ($MarkupExtensions['...'])
    Markup("...",'inline',"/\.\.\./",'&hellip;');

if ($MarkupExtensions['aquo']) {
## left and right arrows
    Markup("<->",'<<-',"/&lt;--?&gt;/",'&harr;');
    Markup("<-",'<lsa',"/&lt;--?/",'&larr;');
    Markup("->",'>^->',"/--?&gt;/",'&rarr;');

## angle brackets
    Markup("aquo",'>links',"/&lt;&lt;(.*?)&gt;&gt;/",'&laquo;$1&raquo;');
    Markup("lsa",'>aquo',"/&lt;(.*?)\|/",'&lsaquo;$1|');
    Markup("rsa",'>->',"/\|(.*?)&gt;/",'|$1&rsaquo;');
}

## long vowels (macrons)
if ($MarkupExtensions['mac']) {
    $LongVowels = array (
    'A' => '&#256;',
    'a' => '&#257;',
    'E' => '&#274;',
    'e' => '&#275;',
    'I' => '&#298;',
    'i' => '&#299;',
    'O' => '&#332;',
    'o' => '&#333;',
    'U' => '&#362;',
    'u' => '&#363;');
    Markup("mac",'>&',"/(?:&|{)([AaEeIiOoUu])(?:m;|})/e","Macron('$1')");
    Markup('title','>&','/\\(:title\\s(.*?):\\)/ei',
    "PZZ(PCache(\$pagename, 
    array('title' => SetProperty(\$pagename,'title',Macronise(PSS('$1'))))))");
/*  (only required for pmwiki 2.0, not for pmwiki 2.1)
    $FmtP['/\\$Title/e'] = '(@$PCache[$pagename]["title"]) ? Macronise($PCache[$pagename]["title"]) : (($GLOBALS["SpaceWikiWords"]) ? \'$Namespaced\' : \'$Name\')';
*/
    $FmtPV['$Titlespaced'] = 
    '@$page["title"] ? Macronise($page["title"]) : $AsSpacedFunction($name)';
    $FmtPV['$Title'] = 
    '@$page["title"] ? Macronise($page["title"]) : ($GLOBALS["SpaceWikiWords"]
       ? $AsSpacedFunction($name) : $name)';
}

function Macron($vowel) {
  global $LongVowels;
  return $LongVowels[$vowel];
}

function Macronise($text) {
  return preg_replace("/(?:&|{)([AaEeIiOoUu])(?:m;|})/e","Macron('$1')",$text);
}

#### inline markups ####
## '/cite/'
if ($MarkupExtensions["'/"])
    Markup("'/","<'''''","/'\/(.*?)\/'/",'<cite>$1</cite>');

## '@keyboard@'
if ($MarkupExtensions["'@"])
    Markup("'@","<'''''","/'@(.*?)@'/",'<kbd>$1</kbd>');

## ';small caps;'
if ($MarkupExtensions["';"])
    Markup("';","<'''''","/';(.*?);'/",'<span class=\'smallcaps\'>$1</span>');

## [^footnote text^] and [^#^] to list footnotes
## includes a style to tidy line spacing
if ($MarkupExtensions['[^']) {
    Markup("[^",'>links','/\[\^(.*?)\^\]/e',"Footnote('$1')");
    Markup("^[^",'<[^','/^\[\^#\^\]$/e',"'<:block>'.Footnote('#')");
}

function Footnote($foottext) {
  static $fngroup, $fncount, $fntext;
  if ($foottext == "#") {
     $fncount = 0; $fngroup++;
     $r = "<div class='footnote'>&nbsp;</div>$fntext";
     $fntext = '';
  } else {
     $fncount++; $fnid = $fngroup+1 . '_' . $fncount;
     $r = "<a name='fnr$fnid' id='fnr$fnid'></a>".
        "<sup><a href='#fn$fnid'>$fncount</a></sup>";
     $foottext = stripslashes($foottext);
     $fntext .= "<p class='footnote'><a name='fn$fnid' id='fn$fnid'></a>".
        "<sup>$fncount</sup> $foottext <a href='#fnr$fnid'>(&uarr;)</a></p>";
  }
  return $r;
}

## copyright and related entities
if ($MarkupExtensions['copy']) {
    Markup("copy",'inline',"/\([cC]\)/",'&copy;');
    Markup("trade",'inline',"/\((?:tm|TM)\)/",'&trade;');
    Markup("reg",'inline',"/\([rR]\)/",'&reg;');
}

## Q: and A: markup
if ($MarkupExtensions['q&a']) {
    $HTMLStylesFmt['q&a'] = "
p.question { margin-top: 2.0em; }
p.question:first-letter {
    float: left;
    font-family: Old English, Georgia, serif;
    color: #777777;
    font-size: 200%;
    line-height: 1.0em;
    margin-right: 0.2em;
}";
    Markup('^Q:', 'block', '/^Q:(.*)$/', "<:block><p class='question'>Q$1</p>");
/*
    Markup('q&a','<block','/^([QA]):(.*)$/',
    '<:block><p class=\'qanda\'>$1$2</p>');
*/
}
## Z; dropcaps markup
if ($MarkupExtensions['A;'])
    Markup('A;','block','/^([[:upper:]]);(([^;&]*(&[^;]+;)*)*);(.*)$/',
    '<:block><p class=\'drop\'>$1<span class=\'smallcaps\'>$2</span>$5</p>');

## {abbr|abbreviations}, {:term:definitions}, =< left & =>right aligned text
if ($action=="print" || $action=="publish") {
    if ($MarkupExtensions['{|}']) {
        Markup("{|}",'>links',"/\{(.*?)\|(.*?\}?)\}/",'$1 ($2)');
        Markup("{:}",'>&',
        "/\{:((?:\[\[[^\]]+\]\])?\{?[^:\}]*?\}?):(.*?\}?)\}/e",
        "PSS('$1').' ('.Keep('$2').')'");
    }
    if ($MarkupExtensions['=|'])
        Markup('^=>','block','/^=&[gl]t;(.*)$/','<:block>');
    $hide = 2;
} else {
    if ($MarkupExtensions['{|}']) {
        Markup("{|}",'>links',
        "/\{(.*?)\|(.*?\}?)\}/",'<abbr title=\'$2\'>$1</abbr>');
        Markup("{:}",'>&',
        "/\{:((?:\[\[[^\]]+\]\])?\{?[^:\}]*?\}?):(.*?\}?)\}/e",
        "'<dfn title='.Keep(PSS(DfnTitle('$2','$1'))).'>'.PSS('$1').'</dfn>'");
    }
    if ($MarkupExtensions['=|']) {
        Markup('^=>','block','/^=&gt;(.*)$/',
        '<:block><p style=\'text-align: right\'>$1</p>');
        Markup('^=<','block','/^=&lt;(.*)$/',
        '<:block><p style=\'text-align: left\'>$1</p>');
    }
}

function DfnTitle($title,$text) {
    $title = str_replace('"','&quot;',$title);
    $title = (strstr($title,"'")) ? '"'.$title.'"' : "'$title'";
    return (preg_match('/^[[:alnum:]].*$/',$text)) ? 
        "$title class='definition'" : $title;
}

## =| centred text
if ($MarkupExtensions['=|'])
    Markup('^=|','block','/^=\|(.*)$/',
    '<:block><p style=\'text-align: center\'>$1</p>');

## {+insertions+}, {-deletions-}, (:revisions:) and {=sticky notes=}
SDV($hide, isset($_GET['hide']) ? $_GET['hide'] : 0);
$pgnum = isset($_GET['p']) ? "?p=".$_GET['p'] : '';
if ($hide) {
    if ($MarkupExtensions['{='])
        Markup("{=",'inline',"/{=(.*?)=}/",'');
    if ($MarkupExtensions['revisions']) {
        Markup("{+",'inline',"/{\+(.*?)\+}/",'$1');
        Markup("{-",'inline',"/{\-(.*?)\-}/",'');
        if ($hide==1) Markup('revisions','<${fmt}','/\(:revisions:\)/',
            '[[{$Name}?hide=0'.$pgnum.' | Show revisions]]');
        else Markup('revisions','directives','/\(:revisions:\)/','');
    }
} else { 
    if ($MarkupExtensions['{='])
        Markup("{=",'inline',"/{=(.*?)(?:\|\s*(.*?))?=}/e",
        "'<span class=\'stickynote\''.NoteStyle('$2').PSS('>$1</span>')");
    if ($MarkupExtensions['revisions'])
        Markup('revisions','<${fmt}','/\\(:revisions:\)/',
        '[[{$Name}?hide=1'.$pgnum.' | Hide revisions]]');
}

function NoteStyle($color) {
    $colors = array(
            'yellow' => array('ffffa1','ffe53e'),
            'green'  => array('b2ffa1','95ff95'),
            'blue'   => array('71ffff','3ee5ff'),
            'purple' => array('b2c7ff','91b8ff'),
            'pink'   => array('ffc7c7','ffb2b2'),
            'grey'   => array('eeeeee','d4d4d4')
            );
    return ($colors[$color][0]) ? 
        " style='background-color:#" . $colors[$color][0] . 
        "; border-top:2px solid #" . $colors[$color][1] . 
        "; border-bottom:2px solid #" . $colors[$color][1] .";'"
        : '';
}

## !run-in heads!and text
if ($MarkupExtensions['^!!']) {
## add one extra <:vspace> after !headings
    Markup('!!vspace', '<!vspace', "/^(!(?>[^!\n]+![^\n]+)\n)/m", '$1<:vspace>');

    Markup('^!!','<^!','/^!([^!]+)!(.*?)$/',
    '<:block><h5 class=\'runin\'>$1.</h5><p> $2</p>');
}
## figure captions
if ($MarkupExtensions['fig'])
    Markup('fig','<links',
    "/^=figure\s+((?:\[\[)?.*?$ImgExtPattern\"([^\"]*)\"(?:\]\])?)\s*(.*?)$/e",
    "'<:block><div class=\'figure\'><p>'.PSS('$1').'</p><p>'.
    (('$3'=='') ? PSS('$2') : PSS('$3')).'</p></div>'");

if ($MarkupExtensions['::']) {
## tidy :: used merely to indent
    Markup('^::2: :','<^: :2->','/^(:+)(:[^:]+)$/','$1 $2');
    Markup('^: :2->','<^::','/^(:+)\\s+:/e',
    "str_replace(':','-','$1').'&gt;'");

## :: or :+ for multiple <dd> per <dt> and multipar item lists
    Markup('::$','<\\$',"/:[:+]\n/",':+');
    Markup(':+','<block','/^(:+.*?:)((?:.*?:\+.*?)+)$/e',
    "PSS('$1').str_replace(':+','</dd><dd>',PSS('$2'))");
    Markup(':+*','<block','/^([#*]+)((?:.*?:\+.*?)+)$/e',
    "PSS('$1').'<p>'.str_replace(':+','</p><p>',PSS('$2')).'</p>'");
    Markup(':+P','>:+*','/:\+/','<br />&nbsp;&nbsp;&nbsp;');
}

## teaser markups T[:*#] Name#id and (:para Name#id:)
if ($MarkupExtensions['para']) {
    Markup('para','directives',
    "/\(:para\s+(.+?)(?:#([^:\s]+))?(?:\s+(more|edit))?:\)/e",
    "TeaseParagraph(\$pagename,'$1','$2','$3')");
    Markup('tfl','directives',"/^T([:*#]+)\s*(\[\[.+?\]\])/e",
    "TeaserFL(\$pagename,'$1','$2')");
    Markup('tww','directives',
    "/^T([:*#]+)\s*((?:$GroupPattern([\/.]))?$WikiWordPattern)/e",
    "Teaser(\$pagename,'$1','$2')");
    SDV($ParaBadAnchorFmt,"'''\$Anchor''' \$[not found in] \$FullName\n");
    SDV($DefaultTeaserAnchor,'teaser');
    SDV($TeaserMoreFmt,' ([[$FullName | more]])');
    SDV($TeaserEditFmt,' ([[$FullName?action=edit | edit]])');
    SDV($DefaultTeaserTextFmt,'Page [[$Group/$Namespaced]] is undefined.');
}

function TeaserFL($pagename,$markup,$linkword) {
  global $UrlExcludeChars,$DefaultTeaserAnchor;
  if (preg_match('/#wikipublisher\\.[^\\|]+\\|([^\\]]+)/',$linkword,$match))
      $link = $match[1];
  else
      $link = FLRef($linkword);
  if (preg_match("/^\\[\\[(.+?)#([^\\s$UrlExcludeChars]*)/",$linkword,$m)) {
      $link = str_replace('#'.$m[2],'',$link);
      $linkword = str_replace($m[1].'#'.$m[2],$m[1],$linkword);
      $anch = ($m[2]=='') ? $DefaultTeaserAnchor : $m[2];
  } else $anch = '';
  return "$markup$linkword: " . TeaseParagraph($pagename,$link,$anch,'');
}

function FLRef($linkword) {
  $l = preg_replace('/\\s*\\|[^\\]]+/','',$linkword);
  $l = preg_replace('/[^\\]]+-+&gt;\\s*/','',$l);
  $l = preg_replace('/[()]/','',$l);
  return preg_replace('/[#?][^\\s]+/','',$l);
}

function Teaser($pagename,$markup,$linkword) {
  return "$markup$linkword: " . TeaseParagraph($pagename,$linkword,'','');
}

function TeaseParagraph($pagename,$teasername,$teaseranch,$act=NULL) {
  global $ParaBadAnchorFmt,$TeaserMoreFmt,$TeaserEditFmt,$DefaultTeaserAnchor,
    $DefaultTeaserTextFmt;
  $tname = MakePageName($pagename,$teasername);
  if ($tname==$pagename) return "''self reference omitted''";
  if ($act=='edit') $taction = str_replace('$FullName',$tname,$TeaserEditFmt);
  else $taction = '';
  $tpage=RetrieveAuthPage($tname,'read',false,'');
  if (isset($tpage['text'])) $ttext = $tpage['text'];
  else return FmtPageName($DefaultTeaserTextFmt,$teasername);
  $tgroup = FmtPageName('$Group',$tname);
  if ($teaseranch=='') {
      $tpara = CleanParagraph($pagename,$tgroup,
                    substr($ttext,0,strpos($ttext."\n","\n")));
      if ($act=='more') $taction=str_replace('$FullName',$tname,$TeaserMoreFmt);
  } elseif (preg_match("/\\[\\[#+$teaseranch\\]\\]\\n?([^\\n]+)/",$ttext,$m)) {
      $tpara = CleanParagraph($pagename,$tgroup,$m[1]);
      if ($act=='more') 
        $taction = str_replace('$FullName',"$tname#$teaseranch",$TeaserMoreFmt);
  } elseif ($teaseranch==$DefaultTeaserAnchor)
      $tpara = CleanParagraph($pagename,$tgroup,
                    substr($ttext,0,strpos($ttext."\n","\n")));
  else
      $tpara = str_replace('$Anchor',$teaseranch,
                    FmtPageName($ParaBadAnchorFmt,$tname));
  return htmlspecialchars($tpara,ENT_NOQUOTES).$taction;
}

function CleanParagraph($pagename,$group,$para) {
  global $GroupPattern,$WikiWordPattern;
  if (preg_match('/^\\|\\|/',$para)) return "''tabular material omitted''";
  $pgroup = FmtPageName('$Group',$pagename);
  $p = preg_replace("/^[#*!]+\s*/","",$para);
  $p = preg_replace("/^:.*?:/","",$p);
  $p = preg_replace("/^([[:upper:]]);(.*?);/","$1$2",$p);
  $p = preg_replace("/`\\..*?$/","...",$p);
  $p = preg_replace("/\\[@(.*?)@\\]/","@@[=$1=]@@",$p);
  $p = preg_replace("/\\[=(.*?)=\\]/e",'Keep(PSS("$1"))',$p);
  $p = preg_replace("/\\(:title.*?:\\)/","",$p);
  $p = preg_replace("/\\[\\[#.*?\\]\\]/","",$p);
  $p = preg_replace("/([`:\/]?)\\b(($GroupPattern([\\/.]))?$WikiWordPattern)/e",
          'QualifyWLink($pgroup,$group,"$1","$2")',$p);
  $p = preg_replace("/\\[\\[(.*?)\\]\\]/e",
          "'[['.QualifyFLink('$pgroup','$group','$1').']]'",$p);
  $p = str_replace('::','',$p);
  return FmtPageName(preg_replace("/{(\\$.*?)}/",'$1',$p),$pagename);
}

function QualifyWLink($pgroup,$group,$esc,$link) {
  global $WikiWordCount,$WikiWordCountMax,$AbbreviationEnabled,
    $AbbreviationPattern;;
  if ($esc) return "$esc$link";
  if ($pgroup==$group) return $link;
  $wwcount = (isset($WikiWordCount[$link])) ? $WikiWordCount[$link] : 
    $WikiWordCountMax;
  if ($wwcount==0) return $link;
  if ($AbbreviationEnabled && preg_match("/^$AbbreviationPattern$/",$link))
    return $link;
  return (preg_match("/[.\\/]/",$link)) ? $link : QualifiedLink($group,$link);
}

function QualifyFLink($pgroup,$group,$link) {
  if ($pgroup==$group) return $link;
  $l = FLRef($link);
  return (preg_match("/[~!:.\\/]/",$l)) ? $link : 
            str_replace("$l",QualifiedLink($group,$l),$link);
}

function QualifiedLink($grp,$ref) {
  return ($grp.'1'==FmtPageName('$Group',MakePageName($grp.'1.'.$grp,$ref))) ?
         "$grp/$ref" : $ref;
}

## lazy web links (an alternative to the one from Pm)
if ($MarkupExtensions['lazyweb'])
    Markup('lazyweb','<wikilink',
    "/\\bwww\\.[^\\s$UrlExcludeChars]*[^\\s.,?!$UrlExcludeChars]/e",
    "Keep(MakeLink(\$pagename,'http://$0','$0'),'L')");

## enhanced AsSpaced function
if ($MarkupExtensions['spaced']) {
    $SpaceWikiWords = 1;
    $AsSpacedFunction = 'SpaceWikiWords';
    $SpaceWikiWordsFunction = 'SpaceWikiWords';
    $RecentChangesFmt['$SiteGroup.AllRecentChanges'] =
        '* [[$FullName | $Group.$Title]]  . . . $CurrentTime $[by] $AuthorLink'.
        ': [=$ChangeSummary=]';
    $RecentChangesFmt['$Group.RecentChanges'] =
        '* [[$FullName | $Title]]  . . . $CurrentTime $[by] $AuthorLink'.
        ': [=$ChangeSummary=]';
    $DefaultPageTextFmt = 'Describe [[$Group/$Title]] here.';
 #   $FPLByGroupIFmt = "<dd><a href='\$PageUrl'>\$Title</a></dd>";
    $StopList = array(
		'A',
		'An',
		'And',
		'But',
		'By',
		'For',
		'From',
		'In',
		'Is',
		'It',
		'Of',
		'On',
		'Or',
		'The',
		'To',
		'With',
            );
    $UnspacedList = array(
        'Mac ',
        'Mc ',
        'Pm Wiki',
        'Side Bar'
            );
}

function SpaceWikiWords($text) {
  global $StopList,$UnspacedList;
  $text = AsSpaced($text);
#  $text = preg_replace('/([[:lower:]])([[:upper:]\\d])/','$1 $2',$text);
#  $text = preg_replace('/([[:upper:]\\d])([[:upper:]][[:lower:]\\d])/',
#    '$1 $2',$text);
  foreach((array)$StopList as $s)
    $text = preg_replace("/(\\s$s\\s)/e","strtolower('$1')",$text);
  foreach((array)$UnspacedList as $u)
    $text = str_replace($u,str_replace(' ','',$u),$text);
  return $text;
}

## automatic smart quotes
if ($MarkupExtensions['squo']) {
    Markup('nl>','<<nl',"/\s?\n\s*([^<]+?>)/",' $1');
    Markup('<nl','<squo',"/(<[^>]+?)\s*\n\s?/",'$1 ');
    Markup('squo','>style',"/(<.*?>['\"]*)|(.?['\"]+)/e",
    "BypassHTML(PSS('$1'),PSS('$2'))");
    Markup('sq|','>inline',"/(\\[\\[[^|\\]]+\\|)(.*?)(\\]\\])/e",
    "'$1'.SmartenLinkText(PSS('$2')).'$3'");
    Markup('sq->','>inline',"/(\\[\\[)([^\\]]+?)(-+&gt;.*?\\]\\])/e",
    "'$1'.SmartenLinkText(PSS('$2')).'$3'");
}

function SmartenLinkText($txt) {
  global $LinkPattern,$UrlExcludeChars,$ImgExtPattern;
  if (!preg_match("/($LinkPattern)([^$UrlExcludeChars]+$ImgExtPattern)/",$txt)) 
        $txt = preg_replace("/(<.*?>['\"]*)|(.?['\"]+)/e",
                "BypassHTML(PSS('$1'),PSS('$2'))",$txt);
  return $txt;
}

function BypassHTML($hstring,$qstring) {
  if ($qstring=='') {
     $qstring = preg_replace("/.*>/",'',$hstring);
     $hstr = preg_replace("/>.*/",'>',$hstring);
     if ($qstring=='') return $hstr;
     else { if (strstr($hstr,"</")) $qstring = "`" . $qstring;
            return $hstr . SmartenQuotes($qstring); }
  }
  else return SmartenQuotes($qstring);
}

function SmartenQuotes($chars) {
  $s = 0;  $r = '';
  if ($chars[0] =="'" || $chars[0] == '"') {
      $quotes = $chars;
      $char = '';
  } else {
      $quotes = substr($chars,1);
      $char = $chars[0];
      if (strstr("0123456789",$char)) {
         $p = ($quotes[0]=="'") ? "p" : "P";
         $r = "$char&$p" . "rime;";
         $s = 1;
         $char = "`";
      }
  }
  $hands = array('l','r');
  if ($char=="" || strstr(" =-[(",$char)) $hi = 0;
  else $hi = 1;
  if ($char=="`") $char = "";
  $r .= $char;
  $prevq = "x";
  for ($i=$s;$i<strlen($quotes);$i++) {
      $q = ($quotes[$i]=="'") ? "s" : "d";
      if ($q==$prevq) $hi = 1 - $hi;
      $r .= "&$hands[$hi]$q" . "quo;";
      $prevq = $q;
  }
  return $r;
}

## page self-reference format and tool-tip format
SDVA($LinkCleanser, array(
    '/`\..*?$/' => '...',
    "/\\{(\\$.*?)\\}/" => '$1',
    "/\\[\\[([^|\\]]+)\\|\\s*(.*?)\\]\\]($SuffixPattern)/e" =>
            "MakeLink(\$pagename,PSS('$1'),PSS('$2'),'$3','\$LinkText')",
    "/\\[\\[([^\\]]+?)\\s*-+&gt;\\s*(.*?)\\]\\]($SuffixPattern)/e" =>
            "MakeLink(\$pagename,PSS('$2'),PSS('$1'),'$3','\$LinkText')",
    '/\\[\\[#([A-Za-z][-.:\\w]*)\\]\\]/' => "",
    "/\\[\\[(.*?)\\]\\]($SuffixPattern)/e" =>
            "MakeLink(\$pagename,PSS('$1'),NULL,'$2','\$LinkText')",
    '/[\\[\\{](.*?)\\|(.*?)[\\]\\}]/' => '$1',
    "/`(($GroupPattern([\\/.]))?($WikiWordPattern))/" => '$1',
    "/$GroupPattern\\/($WikiWordPattern)/" => '$1',
            ));
if ($MarkupExtensions['links']) {
    SDV($WikiStylePattern,'%%|%[A-Za-z][-,=:#\\w\\s\'"().]*%');
    $oLinkPageFunction = $LinkFunctions['<:page>'];
    $LinkFunctions['<:page>'] = 'LinkPageTitle';
    if ($action=='browse') {
        $LinkPageSelfFmt = "<b class='selflink'>\$LinkText</b>";
        $HTMLStylesFmt['selfref'] = "li.browse b { font-weight: normal; }";
    } elseif ($action=='edit') 
        $HTMLStylesFmt['selfref'] = "li.edit a { border-bottom: 1px dotted; }";
    elseif ($action=='diff')
        $HTMLStylesFmt['selfref'] = "li.diff a { border-bottom: 1px dotted; }";
    elseif ($action=='upload')
        $HTMLStylesFmt['selfref'] = "li.upload a { border-bottom: 1px dotted; }";
    elseif ($action=='searchinsitu')
        $HTMLStylesFmt['selfref'] = "form a.search { border-bottom: 1px dotted; }";
    $LinkPageExistsTitleFmt = 
    "<a class='wikilink' href='\$LinkUrl' title=\"\$ToolTip\">\$LinkText</a>";
    $LinkPageCreateFmt = 
    "<a class='createlinktext' href='\$PageUrl?action=edit' title='Create page'>\$LinkText</a><a 
  class='createlink' href='\$PageUrl?action=edit'>?</a>";
#    $oLinkUrlFunction = $LinkFunctions['http:'];
#    $LinkFunctions['http:'] = 'LinkUrlImg';
#    $UrlLinkImgFmt = 
#    "<a class='urllinkimg' href='\$LinkUrl' rel='nofollow'>\$LinkText</a>";
}

function LinkPageTitle($pagename,$imap,$path,$title,$txt,$fmt=NULL) {
    global  $oLinkPageFunction,$LinkPageExistsTitleFmt,$UrlExcludeChars;
    if ($fmt!='') 
        return $oLinkPageFunction($pagename,$imap,$path,$title,$txt,$fmt);
    if (preg_match("/^([^#?]+)(?:#([^\\s$UrlExcludeChars]*))?$/",$path,$match)) {
        $tgtname = MakePageName($pagename,$match[1]); $anch=@$match[2];
        if (PageExists($tgtname) && $tgtname!=$pagename) {
            $title = TitleParagraph($tgtname,$anch);
            if ($title) 
                $fmt = str_replace('$ToolTip',$title,$LinkPageExistsTitleFmt);
        }
    }
    return $oLinkPageFunction($pagename,$imap,$path,$title,$txt,$fmt);
}

function TitleParagraph($pagename,$anch) {
    global $LinkCleanser, $WikiStylePattern, $ParaBadAnchorFmt;
    $refpage = ReadPage($pagename); $para = '';
    $title = ($anch=='') ? 
        preg_match("/^(?:!+|:.*?:)\\s*(?:\\[\\[#.*?\\]\\])?([^\\n]+)/",
                $refpage['text'],$match) :
        preg_match("/\\[\\[#+$anch\\]\\]\\n?([^\\n]+)/",
                $refpage['text'],$match);
    if ($title) {
        $para = preg_replace("/!.*?$/",'',$match[1]);
        $para = preg_replace("/(''+|@@)(.*?)\\1/",'$2',$para);
        $para = preg_replace("/'([-_^;+\\/])(.*?)\\1'/",'$2',$para);
        $para = preg_replace("/\\[([@=]|[-+]+)(.*?)\\1\\]/",'$2',$para);
        $para = preg_replace("/$WikiStylePattern/",'',$para);
        foreach ($LinkCleanser as $p => $c) $para = preg_replace($p,$c,$para);
        $para = 
            htmlentities(str_replace('"','&quot;',str_replace('`','',$para)));
    } elseif ($anch!='') 
#    $para = str_replace("'''","'",
#    str_replace('$Anchor',$anch,FmtPageName($ParaBadAnchorFmt,$pagename)));
        $para = $anch;
    return $para;
}
/*
function LinkUrlImg($pagename,$imap,$path,$title,$txt,$fmt=NULL) {
  global $UrlLinkImgFmt, $oLinkUrlFunction;
  if (isset($fmt))
     return  $oLinkUrlFunction($pagename,$imap,$path,$title,$txt,$fmt);
  if (preg_match("/^<img/",$txt)) 
     return $oLinkUrlFunction($pagename,$imap,$path,$title,$txt,$UrlLinkImgFmt);
  return $oLinkUrlFunction($pagename,$imap,$path,$title,$txt,$fmt);
}
*/
## ||table attributes
Markup('^||','>^||||','/^\\|\\|(.*)$/e',
  "PZZ(\$GLOBALS['BlockMarkups']['table'][0] = PSS('<table '. QuoteAttrs('$1') . '>'))");

function QuoteAttrs($attr) {
  return preg_replace('/([a-zA-Z])\\s*=\\s*([^\'"]\\S*)/',"\$1='\$2'",$attr);
}

?>
