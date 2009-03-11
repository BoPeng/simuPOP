<?php if (!defined('PmWiki')) exit();
/*
    page rename script, version 2.0.16
    Copyright 2004, 2005 John Rankin (john.rankin@affinity.co.nz)
    This file adds page rename capability to pmwiki 2, inter or intragroup,
    via action=rename.
    As a by-product it qualifies unqualified links in cross-group includes,
    with an (:includeg Group.PageName:) directive.
    It also provides a list of links on a page with the (:linkslist page:)
    directive or via action=links.
    
    This file extends PmWiki; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  See pmwiki.php for full details.
*/

SDV($RedirectToRenameFmt,'(:redirect [[$RenameText]]:)');
SDV($HandleActions['links'],'HandlePageLinks');
SDV($PageRenameFmt,
    "<h1 class='wikiaction'>$[Rename] <a href='\$PageUrl'>\$FullName</a></h1>
    <form action='\$PageUrl' method='post'>
    <input type='hidden' name='n' value='\$FullName' />
    <input type='hidden' name='action' value='postrename' />
    \$RenameGroup
    <input type='text' name='renametext' value='\$Name' size='25' /><br/>
    $[Add:] <input type='radio' name='addgroup' value='' checked='checked' /><em>$[nothing]</em> 
    <input type='radio' name='addgroup' value='\$Group' /><em>\$Group</em> 
    <input type='radio' name='addgroup' value='new' /><em>$[new group]</em> $[to any Unqualified Links]
    <input type='submit' value='$[Rename]' /></form>
    <p class='vspace'></p><p>$[Unqualified Links]:</p>
    \$UnqualifiedLinks");
SDV($HandleActions['rename'],'HandleRename');
SDV($HandleActions['postrename'],'HandlePostRename');
$FPLFunctions['pickgroup'] = 'FPLPickGroup';
$FPLFunctions['pglist'] = 'FPLPageLinks';
$FPLFunctions['pgbygroup'] = 'FPLLinksByGroup';

Markup('linkslist', 'directives','/\\(:linkslist\\s+(?:(publish)=)?(.*?)\\s*:\\)/ei',
  "'<:block>'.Keep(FmtLinksList(MakePageName(\$pagename,'$2'), array('o' => 'fmt=pgbygroup','list'=>'all','action'=>'$1')))");

Markup('includeg','>if',
    "/\\(:includeg\\s+($GroupPattern(?:[\\/.])$NamePattern)(.*?):\\)/e",
    "PRR().IncludeGText(\$pagename,'$1','$2')");

function IncludeGText($pagename,$page,$opts) {
  global $GroupPattern,$WikiWordPattern;
  $ogroup = FmtPageName('$Group',$page);
  $ngroup = FmtPageName('$Group',$pagename);
  return 
    preg_replace("/(\\[[=@].*?[=@]\\])|(\\[\\[.*?\\]\\])|([`:\/\$])?\\b(($GroupPattern([\\/.]))?$WikiWordPattern)/e",
    'QualifyUnqualifiedLinks($ngroup,$ogroup,"$0")',
    IncludeText($pagename,"include $page$opts"));
}

function HandlePageLinks($pagename) {
  global $HandlePageLinksFmt,$PageStartFmt,$PageLinksFmt,$PageEndFmt,$HTMLVSpace;
  SDV($PageLinksFmt,array(
    "<h1 class='wikiaction'>$[Links on] <a href='\$PageUrl'>\$Group: \$Title</a></h1>
    $HTMLVSpace",
    'markup:(:linkslist publish=$FullName:)'));
  SDV($HandlePageLinksFmt,array(&$PageStartFmt,&$PageLinksFmt,&$PageEndFmt));
  PrintFmt($pagename,$HandlePageLinksFmt);
}

function HandleRename($pagename) {
  global $HandleRenameFmt,$PageStartFmt,$PageRenameFmt,$PageEndFmt;
  SDV($HandleRenameFmt,array(&$PageStartFmt,&$PageRenameFmt,&$PageEndFmt));
  $PageRenameFmt = str_replace('$UnqualifiedLinks',
    FmtLinksList($pagename,array('o'=>'fmt=pglist','list'=>'unqualified')),
        str_replace('$RenameGroup',
    FmtGroupList($pagename,array('o'=>'fmt=pickgroup')),$PageRenameFmt));
  PrintFmt($pagename,$HandleRenameFmt);
}

function HandlePostRename($pagename) {
  global $RedirectToRenameFmt,$GroupPattern,$WikiWordPattern;
  $newpagename = MakePageName($pagename,
    stripmagic($_POST['group'].'.'.$_POST['renametext']));
  if (PageExists($newpagename)) { Abort("'$newpagename' already exists"); }
  else { 
    Lock(2);
    $page = RetrieveAuthPage($pagename,"edit");
    if ($page) $ntext = $page['text']; else Abort("cannot get '$pagename'");
    $ogroup = FmtPageName('$Group',$pagename);
    $ngroup = FmtPageName('$Group',$newpagename);
#    Abort('stop for testing');
    if ($_POST['addgroup']) {
        if ($_POST['addgroup']=='new') 
            {$h = $ogroup; $ogroup = $ngroup; $ngroup = $h;}
        if ($ogroup==$ngroup) $ngroup = $ngroup.'1';
        $ntext= 
    preg_replace("/(\\[[=@].*?[=@]\\])|(\\[\\[[^#].*?\\]\\])|([`:\/])?\\b(($GroupPattern([\\/.]))?$WikiWordPattern)/e",
            'QualifyUnqualifiedLinks($ngroup,$ogroup,"$0")',$ntext);
#        Abort(str_replace("\n",'<br/>',str_replace("\n\n",'</p><p>',$ntext)));
#        $page = RetrieveAuthPage($newpagename,"edit");
    }        
    $page['text'] = $ntext;
    WritePage($newpagename,$page);
    $page['text'] = 
                str_replace('$RenameText',$newpagename,$RedirectToRenameFmt);
    WritePage($pagename,$page);
  }
  Redirect($pagename); 
}

function QualifyUnqualifiedLinks($ngroup,$ogroup,$link) {
  global $GroupPattern,$WikiWordPattern;
  $link = htmlentities(stripmagic($link));
  if (preg_match("/^\\[[=@]/",$link)) return $link;
  preg_match("/^(\\[\\[(.*?)\\]\\])|([`:\/\$])?\\b(($GroupPattern([\\/.]))?$WikiWordPattern)$/",
    $link,$m);
  return ($m[0][0]=='[') ? '[['.QualifyFLink($ngroup,$ogroup,$m[2]).']]' :
                            QualifyWLink($ngroup,$ogroup,$m[3],$m[4]);
}

function FmtGroupList($pagename,$opt) {
  global $SearchPatterns,$FPLFunctions;
  $pagelist = ListPages((array)@$SearchPatterns['normal']);
  sort($pagelist);
  $matches = array();
  foreach ($pagelist as $pagefile) $matches[] = array('pagename' => $pagefile);
  if (preg_match('/^([^=]*)=(.*?)$/',$opt['o'],$mat)) $f[$mat[1]] = $mat[2];
  $fmtfn = @$FPLFunctions[$f['fmt']];
  if (!function_exists($fmtfn)) $fmtfn='FPLPickGroup';
  return $fmtfn($pagename,$matches,$opt);
}

function FPLPickGroup($pagename,&$pagelist,$opt) {
  global $GroupPickListFmt;
  SDV($GroupPickListFmt,'<option$Select>$Group</option>');
  $currentgroup = FmtPageName($GroupPickListFmt,$pagename);
  $out = array();
  foreach($pagelist as $item) {
    $pgroup = FmtPageName($GroupPickListFmt,$item['pagename']);
    if (@!$seen[$pgroup]++) {
        $s = ($pgroup==$currentgroup) ? " selected='selected'" : '';
        $out[] = str_replace('$Select',$s,$pgroup);
    }
  }
  return "<select name='group'>" . implode('',$out) . "</select> . ";
}

function FmtLinksList($pagename,$opt) {
  global $FPLFunctions,$FPLSimpleIFmt,$FPLSimpleSFmt,$PublishSearchChecked;
  global $PDFCheckboxFmt,$PDFTypesetFmt,$PDFOptionsFmt,$HTMLVSpace;
  global $FPLByGroupStartFmt,$FPLByGroupEndFmt,$FPLByGroupIFmt,$FPLByGroupCFmt,
    $FPLByGroupSFmt;
  $FPLSimpleIFmt = "<li><a href='\$PageUrl'>\$Title</a></li>\n";
  $FPLSimpleSFmt = "<li><b class='selflink'>\$Title</b></li>\n";
  $FPLByGroupSFmt = "<dd><b class='selflink'>\$Title</b></dd>\n";
  if ($opt['action']=='publish') {
      $checked = ($PublishSearchChecked) ? "checked='checked'" : '';
      $FPLByGroupStartFmt =
    "<form class='publish' action='\$ScriptUrl' method='get'>
    <input type='hidden' name='n' value='\$FullName' />
    <dl class='fplbygroup'>";
      $FPLByGroupEndFmt = "</dl>$HTMLVSpace
    <input type='hidden' name='action' value='publish' />
    <input type='hidden' name='ptype' value='search' />
    <input type='hidden' name='page' value='\$Group/\$Name' />
    <input type='hidden' name='include' value='include' />";
      $FPLByGroupEndFmt .= "$PDFCheckboxFmt$PDFTypesetFmt$PDFOptionsFmt</form>";
      $FPLByGroupIFmt = 
    "<dd><input type='checkbox' name='pagearray[]' value='\$Group/\$Name' $checked/>
     <a href='\$PageUrl'>\$Title</a></dd>\n";
/*
      $FPLByGroupCFmt =
    "<dd><input type='checkbox' name='pagearray[]' value='\$Group/\$Name' $checked/>
     <a class='createlinktext' href='\$PageUrl?action=edit'>\$Title</a><a 
     class='createlink' href='\$PageUrl?action=edit'>?</a></dd>\n";
*/
  } else {
      SDV($FPLByGroupStartFmt, "<dl class='fplbygroup'>");
      SDV($FPLByGroupEndFmt, '</dl>');
      SDV($FPLByGroupGFmt, "<dt><a href='\$ScriptUrl/\$Group'>\$Group</a></dt>\n");
      SDV($FPLByGroupIFmt, "<dd><a href='\$PageUrl'>\$Title</a></dd>\n");
      SDV($FPLByGroupCFmt,
   "<dd><a class='createlinktext' href='\$PageUrl?action=edit'>\$Title</a><a 
   class='createlink' href='\$PageUrl?action=edit'>?</a></dd>\n");
      SDV($FPLByGroupSFmt, "<dd><b class='selflink'>\$Title</b></dd>\n");
  }
  $pagelist = ListPageLinks($pagename,$opt['list']);
  sort($pagelist);
  $matches = array();
  foreach ($pagelist as $pagefile) if ($pagefile!=@$lpage) {
    $matches[] = array('pagename' => $pagefile);
    $lpage = $pagefile;
  }
  if (count($matches)==0) return MarkupToHTML($pagename,"-&gt;''$[No links found].''");
  if (preg_match('/^([^=]*)=(.*?)$/',$opt['o'],$mat)) $f[$mat[1]] = $mat[2];
  $fmtfn = @$FPLFunctions[$f['fmt']];
  if (!function_exists($fmtfn)) $fmtfn='FPLPageLinks';
  return $fmtfn($pagename,$matches,$opt);
}

function FPLLinksByGroup($pagename, &$matches, $opt) {
  global $FPLByGroupStartFmt, $FPLByGroupEndFmt, $FPLByGroupGFmt,
    $FPLByGroupIFmt, $FPLByGroupCFmt, $FPLByGroupSFmt;
  SDV($FPLByGroupStartFmt,"<dl class='fplbygroup'>");
  SDV($FPLByGroupEndFmt,'</dl>');
  SDV($FPLByGroupGFmt,"<dt><a href='\$ScriptUrl/\$Group'>\$Group</a></dt>\n");
  SDV($FPLByGroupIFmt,"<dd><a href='\$PageUrl'>\$Title</a></dd>\n");
  SDV($FPLByGroupCFmt,
   "<dd><a class='createlinktext' href='\$PageUrl?action=edit'>\$Title</a><a 
   class='createlink' href='\$PageUrl?action=edit'>?</a></dd>\n");
  SDV($FPLByGroupSFmt,"<dd><b class='selflink'>\$Title</b></dd>\n");
  $out = array();
  foreach($matches as $pc) {
    $pgroup = FmtPageName($FPLByGroupGFmt, $pc['pagename']);
    if ($pgroup != @$lgroup) { $out[] = $pgroup; $lgroup = $pgroup; }
    $out[] = ($pc['pagename']==$pagename) ? 
        FmtPageName($FPLByGroupSFmt, $pc['pagename']) :
        ((PageExists($pc['pagename'])) ?
        FmtPageName($FPLByGroupIFmt, $pc['pagename']) :
        FmtPageName($FPLByGroupCFmt, $pc['pagename']));
  }
  return FmtPageName($FPLByGroupStartFmt, $pagename) . implode('', $out) .
             FmtPageName($FPLByGroupEndFmt, $pagename);
}

function FPLPageLinks($pagename, &$matches, $opt) {
  global $FPLSimpleStartFmt, $FPLSimpleIFmt, $FPLSimpleSFmt, $FPLSimpleEndFmt;
  SDV($FPLSimpleStartFmt, "<ul class='fplsimple'>");
  SDV($FPLSimpleEndFmt, "</ul>");
  SDV($FPLSimpleIFmt, "<li><a href='\$PageUrl'>\$FullName</a></li>\n");
  SDV($FPLSimpleSFmt,"<li><b class='selflink'>\$FullName</b></li>\n");
  $out = array();
  foreach($matches as $pc) 
    $out[] = ($pc['pagename']==$pagename) ?
        FmtPageName($FPLSimpleSFmt, $pc['pagename']) :
        FmtPageName($FPLSimpleIFmt, $pc['pagename']);
  return FmtPageName($FPLSimpleStartFmt, $pagename) . implode('', $out) .
             FmtPageName($FPLSimpleEndFmt, $pagename);
}
   
function ListPageLinks($pagename,$list) {
  global $GroupPattern,$WikiWordPattern;
  $g = FmtPageName('$Group',$pagename);
  $dummy = 'AB101BA';
  $r = array();
  $page = RetrieveAuthPage($pagename,(($list=='all') ? 'read' : 'edit'));
  if ($page) $text = $page['text']; else Abort("cannot rename '$pagename'");
  PCache($pagename,$page);
  preg_match_all("/(\\[[=@].*?[=@]\\])|(\\[\\[([^#].*?)\\]\\])|([`:\/\$])?\\b(($GroupPattern([\\/.]))?$WikiWordPattern)/",
    $text,$m);
  for ($i=0;$i<count($m[0]);$i++) {
    $p = QualifyUnqualifiedLinks($g,$dummy,preg_replace("/\\{(\\$.*?)\\}/e",
        "FmtPageName('$1',\$pagename)",$m[0][$i]));
    if (strstr($p,$dummy)) {
        $p = str_replace("$dummy/","$g.",$p);
        $r[] = ($m[0][$i][0]=='[') ? MakePageName($pagename,FLRef($p)) : $p;
    } elseif ($list=='all')
     
        if (preg_match("/^($GroupPattern([\\/.])$WikiWordPattern)/",$p,$w))
            $r[] = str_replace('/','.',$w[1]);
        elseif (preg_match("/^\\[\\[(.*?)\\]\\]/",$p,$f)) {
            $l = FLRef($f[1]);
            if (preg_match("/^[^~!:.\\/]+([.\\/][^.\\/]+)?$/",$l,$fl))
                $r[] = MakePageName($pagename,$l);
        }
    
  }
  return $r;
}

?>