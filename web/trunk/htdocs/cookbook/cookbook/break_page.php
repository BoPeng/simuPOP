<?php if (!defined('PmWiki')) exit();
/*
    original (:breakpage:) script
    Copyright 2004 Patrick R. Michaud (pmichaud@pobox.com)
    
    modifications to use ____ and edit page chunks 
    version 2.0.7
    Copyright 2004 John Rankin (john.rankin@affinity.co.nz)
    
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.    
*/
$HTMLStylesFmt['breakpage'] = "
div.breaklist { text-align: right; }
div.breaklist strong { background-color: yellow; }
div.breakpage { text-align: right; }";
SDVA($LinkCleanser, array(
    '/`\..*?$/' => '...',
    "/\\[\\[([^|\\]]+)\\|(.*?)\\s*\\]\\]($SuffixPattern)/e" =>
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

function PageBreak($pagename,$text) {
  $p = explode('____',$text."\n\n");
  $n = @$_REQUEST['p'];
  if ($n<1) $n=1;  if ($n>count($p)) $n=count($p);
  $out[] = "<div class='breaklist'>Page";
  for($i=1;$i<=count($p);$i++) {
    $tooltip = PageBreakTip($p[$i-1]);
    if ($i==$n) $out[] = " <strong>$n</strong>";
    else $out[] = FmtPageName(" <a href='\$PageUrl?p=$i' title=$tooltip>$i</a>",$pagename);
  }
  $out[] = '</div>';
  $ed = FmtPageName("<a href='\$PageUrl?p=$n&action=editpage'>edit</a>", $pagename);
  return '<:block>'.Keep("<div class='breakpage'>Page $n ($ed)</div>").
    "(:nl:)".$p[$n-1].'<:block>'.Keep(implode('',$out));
}

function PageBreakTip($text) {
  global $LinkCleanser;
  $title = preg_match("/^\\n*(?:\\(:.*?:\\))?(?:!+|:.*?:)?\\s*(?:\\[\\[#.*?\\]\\])?([^\\n]+)/",
                $text,$match);
  if ($title) {
    $para = preg_replace("/!.*?$/",'',$match[1]);
    $para = preg_replace("/(''+|@@)(.*?)\\1/",'$2',$para);
    $para = preg_replace("/'([-_^;+\\/])(.*?)\\1'/",'$2',$para);
    $para = preg_replace("/\\[([@=]|[-+]+)(.*?)\\1\\]/",'$2',$para);
    foreach ($LinkCleanser as $p => $c) $para = preg_replace($p,$c,$para);
    $para = str_replace('"','&quot;',str_replace('`','',$para));
  } else $para = 'Go to page';
  return (strstr($para,"'")) ? '"'.$para.'"' : "'$para'";
}
  
function HandleEditPage($pagename) {
  global $IsPagePosted,$EditFields,$ChangeSummary,$EditFunctions,$FmtV,$Now,
    $HandleEditFmt,$PageStartFmt,$PageEditFmt,$PagePreviewFmt,$PageEndFmt,
    $GroupHeaderFmt,$GroupFooterFmt,$EnablePost;
  if ($_REQUEST['cancel']) { Redirect($pagename); return; }
  $IsPagePosted = false;
  $PageEditFmt = "<div id='wikiedit'>
  <a id='top' name='top'></a>
  <h1 class='wikiaction'>$[Editing <a href='\$PageUrl' accesskey='b'>\$Group \$Name</a>]</h1>
  <form method='post' action='\$PageUrl?action=editpage'>
  <input type='hidden' name='action' value='editpage' />
  <input type='hidden' name='n' value='\$FullName' />
  <input type='hidden' name='basetime' value='\$EditBaseTime' />
  <input type='hidden' name='prechunk' value=\"\$PreChunk\" />
  <input type='hidden' name='p' value='\$PNum' />
  <input type='hidden' name='postchunk' value=\"\$PostChunk\" />
  \$EditMessageFmt
  <textarea name='text' rows='25' cols='64'
    onkeydown='if (event.keyCode==27) event.returnValue=false;'
    >\$EditText</textarea><br />
  $[Summary]: <input type='text' name='csum' value =\"\$CSum\" size='56' /><br />
  $[Author]: <input type='text' name='author' value='\$Author' />
  <input type='checkbox' name='diffclass' value='minor' \$DiffClassMinor />
    $[This is a minor edit]<br />
  <input type='submit' name='post' value=' $[Save] ' />
  <input type='submit' name='preview' value=' $[Preview] ' />
  <input type='submit' name='cancel' value=' $[Cancel] ' /></form></div>";
  if (@$_REQUEST['preview']) 
    $PagePreviewFmt = array(
  "<a id='preview' name='preview'></a>",
  "<p id='wikisubtitle'>$[Page is unsaved]</p>
   <h1 id='headbar'>$[Preview] \$Group \$Name</h1><div id='wikipreview'>",
  "\$PreviewText",
  "</div><div id='footbar'>",
  "<div class='footnavright'><a href='#top'>$[Top]</a></div>",
  "<div class='footnavleft'><b>$[End of preview: remember to save]</b></div></div>"
    );
  $PageChunks = array('prechunk','postchunk');
  Lock(2);
  $page = RetrieveAuthPage($pagename,'edit');
  if (!$page) Abort("?cannot edit $pagename"); 
  PCache($pagename,$page);
  $new = $page;
  $p = explode('____',$new['text']);
  $n = @$_REQUEST['p'];
  if ($n<1) $n=1;  if ($n>count($p)) $n=count($p);
  $new['text'] = $p[$n-1];
  foreach((array)$EditFields as $k) 
    if (isset($_POST[$k])) $new[$k]=str_replace("\r",'',stripmagic($_POST[$k]));
  foreach((array)$PageChunks as $c) {
    $$c = '';
    if (@$_POST[$c]) $$c  = str_replace("\r",'',stripmagic($_REQUEST[$c]));
  }
  if ($ChangeSummary) $new["csum:$Now"] = $ChangeSummary;
  $EnablePost &= (@$_POST['post'] || @$_POST['postedit']);
  if (@$_POST['post'])
    $new['text'] = $prechunk."\n".$new['text']."\n".$postchunk;
  elseif (@$_POST['preview']) {
    if ($n>1) $GroupHeaderFmt = '';
    $GroupHeaderFmt .= '=&gt;(Page '.@$_REQUEST['p'].' of '.count($p).')(:nl:)';
    if ($n<count($p)) $GroupFooterFmt = '';
  } else
    for ($i=1;$i<=count($p);$i++)
      if ($i<$n) $prechunk .= $p[$i-1] . '____';
      elseif ($i>$n) $postchunk .= '____' . $p[$i-1];
  foreach((array)$EditFunctions as $fn) $fn($pagename,$page,$new);
  if ($IsPagePosted) { Redirect($pagename,"\$PageUrl?p=$n"); return; }
  $FmtV['$PreChunk'] = str_replace('"','&quot;',
        str_replace('$','&#036;',htmlspecialchars($prechunk, ENT_NOQUOTES)));
  $FmtV['$PostChunk'] = str_replace('"','&quot;',
        str_replace('$','&#036;',htmlspecialchars($postchunk, ENT_NOQUOTES)));
  $FmtV['$CSum'] = str_replace('"','&quot;',
    str_replace('$','&#036;',htmlspecialchars(@$_POST['csum'], ENT_NOQUOTES)));
  $FmtV['$PNum'] = $n;
  $FmtV['$DiffClassMinor'] = 
    (@$_POST['diffclass']=='minor') ?  "checked='checked'" : '';
  $FmtV['$EditText'] = 
    str_replace('$','&#036;',htmlspecialchars(@$new['text'],ENT_NOQUOTES));
  $FmtV['$EditBaseTime'] = $Now;
  SDV($HandleEditFmt,array(&$PageStartFmt,
    &$PageEditFmt,&$PagePreviewFmt,&$PageEndFmt));
  PrintFmt($pagename,$HandleEditFmt);
}

if ($action=='browse') 
  Markup('breakpage','>include','/^.*____.*$/se',
    "PageBreak(\$pagename,PSS('$0'))");
elseif ($action=='edit' || $action=='editpage') {
  Markup('breakpage','>include','/____/',
    "\n<div class='breakpage'>&mdash; <em>Page Break</em> &mdash;</div>\n");
  if ($action=='editpage') $HandleActions['editpage'] = 'HandleEditPage';
} elseif ($action=='print')
  Markup('breakpage','>include','/____/','(:div class="section":)');
else
  Markup('breakpage','directives','/____/','');

?>