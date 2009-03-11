<?php

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

function QualifyWLink($pgroup,$group,$esc,$link) {
  global $WikiWordCount,$WikiWordCountMax;
  if ($esc) return "$esc$link";
  if ($pgroup==$group) return $link;
  $wwcount = (isset($WikiWordCount[$link])) ? $WikiWordCount[$link] : 
    $WikiWordCountMax;
  if ($wwcount==0) return $link;
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

function FLRef($linkword) {
  $l = preg_replace('/\\s*\\|[^\\]]+/','',$linkword);
  $l = preg_replace('/[^\\]]+-+&gt;\\s*/','',$l);
  return preg_replace('/#[^\\s]+/','',$l);
}

?>