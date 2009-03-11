<?php if (!defined('PmWiki')) exit();
/*
    The pagetoc script adds support for automatically generating
    a table of contents for a wiki page.

    Version 2.0.31 (development version; works with PmWiki 2.1.0 or above)

    Copyright 2004, 2005 John Rankin (john.rankin@affinity.co.nz)
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
*/
SDV($TocSize,'smaller');
SDV($TocFloat,false);
$HTMLStylesFmt['toc'] = "
span.anchor {
	float: left;
	font-size: 10px;
	margin-left: -10px;
	width: 10px;
    position:relative; top:-0.1em;
	text-align: center;
}
span.anchor a { text-decoration: none; }
span.anchor a:hover { text-decoration: underline; }
ol.toc { text-indent:-20px; list-style: none; }
ol.toc ol.toc { text-indent:-40px; }";
$HTMLStylesFmt['tocf'] = "
div.tocfloat { font-size: $TocSize; margin-bottom: 10px;
    border-top: 1px dotted #555555; border-bottom: 1px dotted #555555;
    padding-top: 5px; padding-bottom: 5px; 
    width: 38%; float: right; margin-left: 10px; clear: right;
    margin-right:-13px; padding-right: 13px; padding-left: 13px;
    background-color: #eeeeee; }
div.toc { font-size: $TocSize; 
    padding: 5px; border: 1px dotted #cccccc;
    background: #f7f7f7;
    margin-bottom: 10px; }
div.toc p { background-color: #f9f6d6;
    margin-top:-5px;   padding-top: 5px;
    margin-left:-5px;  padding-left: 5px;
    margin-right:-5px; padding-right: 5px;
    padding-bottom: 3px;
    border-bottom:  1px dotted #cccccc; }"; 
SDV($ToggleText, array('hide', 'show'));
$HTMLHeaderFmt['toggle'] = "<script type=\"text/javascript\">
function toggle(obj) {
    var elstyle = document.getElementById(obj).style;
    var text    = document.getElementById(obj + \"tog\");
    if (elstyle.display == 'none') {
        elstyle.display = 'block';
        text.innerHTML = \"{$ToggleText[0]}\";
    } else {
        elstyle.display = 'none';
        text.innerHTML = \"{$ToggleText[1]}\";
    }
}
</script>";

## in-page cross-references
Markup('[[#|#','>nl1','/\[\[#([A-Za-z][-.:\w]*)\s*\|\s*#\]\]/e',
    "'[[#$1 | '.CrossReference(\$pagename,\$x,'$1').']]'");
Markup('[[#|*','<[[|','/\[\[#([A-Za-z][-.:\w]*)\s*\|\s*\*\]\]/',
    '[[#$1 | $1]]');
Markup('[[#|+','<[[|','/\[\[#([A-Za-z][-.:\w]*)\s*\|\s*\+\]\]/',
    '[[#$1 | Back to $1]]');
Markup("[^#",'<[[#|#','/\[\^#([A-Za-z][-.:\w]*)\^\]/e',"Shortcut(\$x,'$1')");
SDVA($LinkCleanser, array(
    '/`\..*?$/' => '...',
    "/\\[\\[([^|\\]]+)\\|\\s*(.*?)\\]\\]($SuffixPattern)/e" =>
            "MakeLink(\$pagename,PSS('$1'),PSS('$2'),'$3','\$LinkText')",
    "/\\[\\[([^\\]]+?)\\s*-+&gt;\\s*(.*?)\\]\\]($SuffixPattern)/e" =>
            "MakeLink(\$pagename,PSS('$2'),PSS('$1'),'$3','\$LinkText')",
    '/\\[\\[#([A-Za-z][-.:\\w]*)\\]\\]/' => "",
    "/\\[\\[(.*?)\\]\\]($SuffixPattern)/e" =>
            "MakeLink(\$pagename,PSS('$1'),NULL,'$2','\$LinkText')",
    '/[\\[\\{](.*?)\\|(.*?)[\\]\\}]/' => '$1',
    "/`(($GroupPattern([\\/.]))?($WikiWordPattern))/" => '$1',
    "/$GroupPattern\\/($WikiWordPattern)/" => '$1'
            ));

function CrossReference($pagename,$text,$anchor) {
  global $LinkCleanser;
  $r = Shortcut($text,$anchor);
  foreach ($LinkCleanser as $p => $c) $r = preg_replace($p,$c,$r);
  return $r;
}

function Shortcut($text,$anchor) {
  if (preg_match("/\\[\\[#+$anchor\\]\\]\\n?([^\n]+)/",$text,$match)) {
    return preg_replace("/^[#*!:]+\s*/","",
                preg_replace("/([^!]+)!.+/","$1",$match[1]));
  } else {
    return "<em>$anchor</em> not found";
  }
}

## [[##visibleanchor]]
SDV($VisibleAnchor,'&sect;');
SDV($VisibleAnchorLinks,false);
SDV($DefaultTocAnchor,'toc');
$RefOrTitle = ($VisibleAnchorLinks) ? 'href' : 'title';

## autonumber anchors
Markup('^!#','<links','/^(!+|Q?:)#(#?)/e',"'$1'.TocAnchor('$2')");
   
function TocAnchor($visible) {
  global $DefaultTocAnchor;
  static $toccounter;
  return "[[$visible#$DefaultTocAnchor" . ++$toccounter . "]]";
}

## (:markup:) that excludes heading markup examples
function HMarkupMarkup($pagename, $lead, $texta, $textb) {
  return "$lead<:block>" .
    Keep("<table class='markup' align='center'><tr><td class='markup1'><pre>" .
      wordwrap($texta, 70) .  "</pre></td></tr><tr><td class='markup2'>") .
    "\n$textb\n(:divend:)</td></tr></table>\n";
}

Markup('`markup','<markup',
  "/(^|\\(:nl:\\))\\(:markup:\\)[^\\S\n]*\\[([=@])((?:\n`\\.!+.*?)+)\\2\\]/seim",
  "HMarkupMarkup(\$pagename, '$1', PSS(str_replace('`.','','$3')), PSS('$3'))");
#Markup("`.",'>links',"/`\./",''); ## included in extendmarkup.php

SDV($format, ''); # Hack - /Christian
## page table of contents
$IdPattern = "[A-Za-z][-.:\w]*";
if ($format=='pdf') {
    SDV($DefaultTocTitle,'Contents');
    SDV($TocHeaderFmt,
        '[[#toc]]<tbook:visual markup="bf">$TocTitle</tbook:visual>');
    SDV($RemoteTocFmt,
        '<tbook:visual markup="bf">Contents of [[$Toc(#toc)]]</tbook:visual>');
} else {
    SDV($DefaultTocTitle,'On this page...');
    SDV($TocHeaderFmt,'[[#toc]]<b>$TocTitle</b>');
    SDV($RemoteTocFmt,'<b>On page [[$Toc(#toc)]]...</b>');
}
SDV($NumberToc,true);
SDV($L1TocChar, '.');
SDV($OmitQMarkup,false);

if ($action=="print" || $action=="publish") {
    Markup('[[##','<[[#','/\[\[##([A-Za-z][-.:\w]*)\]\]/','[[#$1]]');
    if ($action=='publish') Markup('toc','>include',
        '/\(:([#\*])?toc(?:-(float|hide))?(?:\s+anchors=(v)isible)?(?:\s+(.*?))?:\)/', '');
    Markup('tocback','directives','/\(:toc-back(?:\s+(.*?))?:\)/','');
} else {
    Markup('[[##','<[[#','/\[\[##([A-Za-z][-.:\w]*)\]\]/e',
    "Keep(\"<span class='anchor'><a name='$1' id='$1' $RefOrTitle='#$1'>$VisibleAnchor</a></span>\",
    'L')");
}
Markup('toc','>nl1',
    '/\(:([#\*])?toc(?:-(float|hide))?(?:\s+anchors=(v)isible)?(?:\s+(.*?))?(?:\s+(Q))?:\)(.*)$/se', 
    "TableOfContents(\$pagename,'$1','$2',PSS('$4'),'$5',PSS('$6')).
    TocEntryAnchors('$3',PSS('$6'))");
SDV($TocBackFmt,'&uarr; Contents');
Markup('tocback','directives','/\(:toc-back(?:\s+(.*?))?:\)/e',
    "'[[#toc | '.TocLinkText(PSS('$1')).']]'");
Markup('tocpage','directives','/\(:toc-page\s+(.*?)(?:\s+self=([01]))?:\)/e',
    "RemoteTableOfContents(\$pagename,'$1','$2')");

function RemoteTableOfContents($pagename,$ref,$self=0) {
    global $TocHeaderFmt,$RemoteTocFmt;
    $oTocHeader = $TocHeaderFmt;
    $TocHeaderFmt = str_replace('$Toc',$ref,$RemoteTocFmt);
    $tocname = MakePageName($pagename,$ref);
    if ($tocname==$pagename && $self==0) return '';
    $tocpage=RetrieveAuthPage($tocname,'read',false);
    $toctext=@$tocpage['text'];
    if (preg_match('/\(:([#\*])?toc(?:-(float|hide))?(?:\s+anchors=(v)isible)?(?:\s+(.*?))?(?:\s+(Q))?:\)(.*)$/se',$toctext,$m))
        $toc = str_replace('[[#',"[[$ref#",
            TableOfContents($tocname,$m[1],'page','',$m[5],PSS($m[6])));
    $TocHeaderFmt = $oTocHeader;
    return $toc;
}

function TocLinkText($text) {
    global $TocBackFmt;
    if ($text) $TocBackFmt = $text;
    return $TocBackFmt;
}

function TocEntryAnchors($visible,$text) {
    global $IdPattern;
    return preg_replace("/\n(!+|Q:)((\[\[#+$IdPattern\]\])|##?)?/e",
                '"\n$1".InsertAnchor($visible,"$1","$2")',$text);
}

function InsertAnchor($visible,$h,$mark) {
  global $OmitQMarkup, $NumberToc, $L1TocChar;
  static $l1,$l2,$toc1,$toc2;
  if ($h=='Q:' && $OmitQMarkup) return $mark;
  if ($mark=='') $visibility = ($visible=='') ? '#' : '##';
  else $visibility = $mark;
  if ($h=='Q:') return $visibility;
  $r = '';
  $len = strlen($h);
  if ($l1==0) { $l1 = $len; }
  else if ($len!=$l1 && $l2==0) { $l2 = $len; }
#  if ($l1==$len || $l2==$len) $r = $visibility;
  if ($l1==$len) { 
    $toc1++; $toc2 = 0; $r = $visibility;
    if ($NumberToc) $r .= "$toc1$L1TocChar&ensp; ";
  } elseif ($l2==$len) { 
    $toc2++; $r = $visibility;
    if ($NumberToc) $r .= "$toc1.$toc2&ensp; "; 
  }
  return $r;
}

function TableOfContents($pagename,$number,$float,$title,$includeq,$text) {
    global $DefaultTocTitle,$TocHeaderFmt,$IdPattern,$NumberToc,$OmitQMarkup,
        $format,$L1TocChar,$DefaultTocAnchor,$TocFloat,$HTMLHeaderFmt,
        $ToggleText;
    if ($includeq)    $OmitQMarkup = (!$OmitQMarkup);
    if ($float=='float') $TocFloat = (!$TocFloat);
    $l1 = 0; $l2 = 0; $l3 = 0;
    $q = 0; $prelen = 1; $counter = 0;
    $r = ''; $toc1 = 0;
    if (!$title) $title = $DefaultTocTitle;
    $toc = str_replace('$TocTitle',$title,$TocHeaderFmt);
    if ($number=='*') $NumberToc = false;
    elseif ($number=='#') $NumberToc = true;
    $closel = 0;
    if ($format=='pdf') {
        $l = 'tbook:item'; 
        $s = ($NumberToc) ? 'tbook:enumerate' : 'tbook:itemize'; $sc = $s;
        $toc = "<tbook:group class='toc'><tbook:p>$toc</tbook:p>".
            "<$sc><$l>\$List</$l></$s></tbook:group>";
    } elseif ($float=='hide') { return '';
    } else {
        $tocid = ($float=='page') ? 'ptocid' : 'tocid';  // remote toc?
        $toggle = " (<a id=\"{$tocid}tog\" href=\"javascript:toggle('$tocid');\">{$ToggleText[0]}</a>)";
        $l = 'li'; $s = ($NumberToc) ? 'ol' : 'ul'; 
        $sc = "$s class='toc'";
        $f = ($TocFloat) ? 'float' : '';
        $toc = "<div class='toc$f'><p>$toc$toggle</p>" . 
               "<$sc id='$tocid'><$l>\$List</$l></$s></div>";
    }
    preg_match_all("/\n(!+|Q?:)\s*(\[\[#+$IdPattern\]\]|#*)([^\n]*)/",$text, $match);
    for ($i=0;$i<count($match[0]);$i++) {
         if ($match[1][$i]==':' || ($match[1][$i]=='Q:' && $OmitQMarkup)) {
            if ($match[2][$i] && $match[2][$i][0]=='#')  $counter++; 
            continue; }
         $len = ($match[1][$i]=='Q:') ? 9 : strlen($match[1][$i]);
         if ($len==9) $q = 1;
         if ($l1==0) { $l1 = $len; $prelen = $l1; }
         if ($len!=$l1 && $l2==0) { $l2 = $len; }
         if ($len!=$l1 && $len!=$l2 && $q==1 && $l3==0) { $l3 = $len; }
         if ($len==$l2 && $l1==9) { $len = $l1; }
         if ($len==$l3) { $len = $l2; }
         if ($len!=$prelen) {
            if ($len==$l1) { $r .= "</$l></$s>"; }
            else if ($len==$l2) { $r .= "<$sc><$l>"; 
                $toc2 = 0;
                $closel = 0;
            }
         }
         if ($len==$l1 || $len==$l2) {
            $prelen = $len;
            if ($len==$l1) {
                $toc1++;
                $tocout = ($NumberToc) ? "$toc1$L1TocChar&ensp;" : '';
            } else {
                $toc2++;
                $tocout = ($NumberToc) ? "$toc1.$toc2&ensp;" : '';
            }
            if ($format=='pdf') $tocout = '';
            $m = preg_replace("/^(\\[\\[#)#/","$1",$match[2][$i]);
            $m = preg_replace("/^#+/",'',$m);
            $t = preg_replace('/%(center|right)%/','',$match[3][$i]);
            if ($closel==1) $r .= "</$l><$l>";
            $closel = 1;
            if (strpos($m,'[#')==1) 
                $r .= $tocout.str_replace(']]',' | '.
                    CrossReference($pagename,"$m$t",
                    preg_replace("/\[\[#(.*?)\]\]/","$1",$m)).']]',$m);
            else {
                $counter++;
                $r .= $tocout."[[#$DefaultTocAnchor$counter | ".
                    CrossReference($pagename,"[[#]]$t","").']]';
            }
         }
    }
    if ($prelen==$l2) $r .= "</$l></$s>";
    if ($r!='') $r = str_replace('$List',$r,$toc);
    return $r;
}

?>