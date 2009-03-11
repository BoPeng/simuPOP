<?php if (!defined('PmWiki')) exit();
/*
	Copyright 2004, 2005 John Rankin john.rankin@affinity.co.nz
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    
    Adds (:commentbox:) and (:commentboxchrono:) markups.
    Version 2.0.10 (arbitrary start point)
    Version 2.0.11 (adds block styles to entries)
    Version 2.0.12 (adds security check on action=comment)
    Version 2.0.13 (adds an addslashes check and tidies code)
    Version 2.0.14 (adds a 'post to' option to commentbox)
    Version 2.0.15 (adds max external links per post and access code)
    Version 2.0.16 (adds classes to box fields and $PostCount)
    Version 2.0.17 (adds basic encryption to access code)
    Version 2.0.18 (restricts ability to post to any page)
    Version 2.0.19 (adds cross-group posting control)
*/

$HTMLStylesFmt['forms'] = "
th.prompt { text-align: right; vertical-align: top; 
    font-size: smaller;
}";
$HTMLStylesFmt['commentbox'] = "
h4.journalhead, h4.messagehead { background:#ffb; border-top:1px dotted #999; }
.message, .journal, .diary {

    padding:3px;

    border:1px dotted #999;
    background:#ffc;

}
.message h4 { margin-top: 0px; }
em.access { font-style: normal; color: #FF2222; }";
SDV($commentboxaccesscode,false);
SDV($DiaryBoxFmt,"<div id='diary'><form action='\$PageUrl' method='post'>
	<input type='hidden' name='n' value='\$FullName' />
	<input type='hidden' name='action' value='comment' />
    <input type='hidden' name='accesscode' value='\$CryptCode' />
	<table width='95%'><tr>
	<th class='prompt'>$[New entry]</th>
	<td><textarea name='text' rows='6' cols='50'></textarea>".
	($commentboxaccesscode ? "</td></tr>
	<tr><th class='prompt'>$[Enter code] <em class='access'>\$AccessCode</em></th>
	<td><input type='text' size='4' maxlength='3' name='access' value='' /> " 
        : "<input type='hidden' name='access' value='\$AccessCode' /><br />").
	"<input type='submit' name='post' value=' $[Post] ' />
	<input type='reset' value='$[Reset]' /></td></tr></table></form></div>");
SDV($CommentBoxFmt,
    "<div class='commentbox'><form action='\$PageUrl' method='post'>
	<input type='hidden' name='n' value='\$FullName' />
	<input type='hidden' name='action' value='comment' />
	<input type='hidden' name='order' value='\$Chrono' />
	<input type='hidden' name='postto' value='\$PostTo' />
    <input type='hidden' name='accesscode' value='\$CryptCode' />
	<table width='95%'><tr>
	<th class='prompt'>$[Comment]</th>
	<td><textarea name='text' rows='6' cols='50'></textarea>
	</td></tr><tr><th class='prompt'>$[Author]</th>
	<td><input type='text' name='author' value='\$Author' size='32' class='inputbox' />".
	($commentboxaccesscode ? "</td></tr>
	<tr><th class='prompt'>$[Enter code] <em class='access'>\$AccessCode</em></th>
	<td><input type='text' size='4' maxlength='3' name='access' value='' /> " 
        : "<input type='hidden' name='access' value='\$AccessCode' />").

	"<input type='submit' name='post' value=' $[Post] ' class='inputbutton' />
	<input type='reset' value='$[Reset]' class='inputbutton' /></td></tr></table></form></div>");
SDV($JournalDateFmt,'%d %B %Y');
SDV($JournalTimeFmt,'%H:%M');
SDV($JournalPattern,'/Journal$/');
SDV($DiaryPattern,'/Diary$/');

SDV($FmtPV['$PostCount'], 'PageStringCount($pn,">>message<<")');

function PageStringCount($pagename,$find) {

   $page = ReadPage($pagename, READPAGE_CURRENT);

   return substr_count($page['text'], $find);

}


if ($action == 'comment') {
    if (auditJP($MaxLinkCount)) 
         SDV($HandleActions['comment'], 'HandleCommentPost');
    else Redirect($pagename);
} else if ($action=='print' || $action=='publish') 
    Markup('cbox','<links','/\(:commentbox(chrono)?.*?:\)/','');
else {  
    Markup('cbox','<links','/\(:commentbox(chrono)?(?:\\s+(\\S.*?))?:\)/e',

    "'<:block>'.Keep(str_replace(array('\$Chrono','\$PostTo'), array('$1','$2'), 
        RandomAccess(FmtPageName(\$GLOBALS['CommentBoxFmt'],\$pagename))))");
    Markup('dbox','<block','/\(:diarybox:\)/e',
        "'<:block>'.RandomAccess(FmtPageName(\$GLOBALS['DiaryBoxFmt'],
        \$pagename))");
    if (preg_match($JournalPattern,$pagename) ||
        preg_match($DiaryPattern,$pagename)) {
            $GroupHeaderFmt .= '(:diarybox:)(:nl:)';
            if (!PageExists($pagename)) $DefaultPageTextFmt = '';
    }
}

function CryptCode($code, $salt=NULL) {
  $pw = substr('MNBVCXZLKJHGFDSAPOIUYTREWQ', $code%26, 1) . "$code";
  return $salt ? crypt($pw, $salt) : crypt($pw);
}

function RandomAccess($fmt) {
  $code  = rand(100,999);
  return str_replace(array('$AccessCode','$CryptCode'),
    array($code,CryptCode($code)),$fmt);
}

function auditJP($MaxLinkCount) {
  SDV($MaxLinkCount, 1);
  if (!(@$_POST['post'] && @$_POST['access'] && 
    (CryptCode($_POST['access'],$_POST['accesscode'])==$_POST['accesscode'])))
     return false;
  preg_match_all('/https?:/',$_POST['text'],$match);
  return (count($match[0])>$MaxLinkCount) ? false : true;
}

function HandleCommentPost($pagename) {
  global $_POST,$JournalPattern,$DiaryPattern,$Author;
  global $AuthFunction,$oAuthFunction,$EditRedirectFmt,$EnablePostToAnyPage,
    $EnableRedirectToCommentPage, $AllowedPostToGroups;
  if (@$_POST['author']=='') $Author = 'anon'; 
  if (@$_POST['postto']) { 
    SDV($EnableRedirectToCommentPage, false);
    if (!$EnableRedirectToCommentPage) SDV($EditRedirectFmt, $pagename);
    SDV($EnablePostToAnyPage, false);
    SDVA($AllowedPostToGroups, array());
    if ($EnablePostToAnyPage) {
        $pg = MakePageName($pagename, $_POST['postto']);
        $gp = FmtPageName('$Group', $pg);
        if ($gp == FmtPageName('$Group', $pagename) || in_array($gp, $AllowedPostToGroups)) 
            $pagename = $pg;
    }
  } 

  SDV($AuthFunction,'PmWikiAuth');
  $oAuthFunction = $AuthFunction;
  $AuthFunction = 'BypassAuth';
  $page = RetrieveAuthPage($pagename, "read");
  if(get_magic_quotes_gpc()==1) $page['text'] = addslashes($page['text']);

  $HandleCommentFunction = (preg_match($JournalPattern,$pagename)) ? 'Journal' : 
    ((preg_match($DiaryPattern,$pagename)) ? 'Diary'   : 'Message');
  $HandleCommentFunction = 'Handle' . $HandleCommentFunction . 'Post';
  $HandleCommentFunction($pagename, $page['text']);
  HandleEdit($pagename);
  exit;
}

function BypassAuth($pagename,$level,$authprompt=true) {
    global $AuthFunction,$oAuthFunction;
    if ($level=='edit') $AuthFunction = $oAuthFunction;
    return $oAuthFunction($pagename,"read",$authprompt);
}

function FormatDateHeading($txt,$datefmt,$fmt) {
  return str_replace($txt,strftime($datefmt,time()),$fmt);
}

## This function handles the layout and posting of a Journal entry
function HandleJournalPost($pagename, $pagetext) {
   global $_POST, $JournalDateFmt, $JournalTimeFmt, $JPDateFmt, $JPTimeFmt,
     $JPItemEndFmt;
   SDV($JPDateFmt,'!!!!%block class=journalhead%$Date');
   SDV($JPTimeFmt,"\n\n>>journal<<\n!\$Time!");
   SDV($JPItemEndFmt,"\n>><<");
   $date = FormatDateHeading('$Date',$JournalDateFmt,$JPDateFmt);
   $entry = $date . FormatDateHeading('$Time',$JournalTimeFmt,$JPTimeFmt) . 
        $_POST['text'] . $JPItemEndFmt;
   $_POST['text'] = (strstr($pagetext, $date)) ?
        str_replace($date, $entry, $pagetext) : 
        "$entry\n\n" . $pagetext;
}

## This function handles the layout and posting of a Diary entry
function HandleDiaryPost($pagename, $pagetext) {
   global $_POST, $JournalDateFmt, $DPDateFmt, $DPItemFmt, $DPItemEndFmt;
   SDV($DPDateFmt,">>diary<<\n!$Date!"); SDV($DPItemFmt,"\n\n: :");
   SDV($DPItemEndFmt,"\n>><<");
   $date = FormatDateHeading('$Date',$JournalDateFmt,$DPDateFmt);
   $entry = $date . $_POST['text'];
   $_POST['text'] = (strstr($pagetext, $date)) ?
        str_replace($date, $entry.$DPItemFmt, $pagetext) :
        "$entry$DPItemEndFmt\n\n" . $pagetext;
}

## This function handles the layout and posting of a Comment entry
function HandleMessagePost($pagename, $pagetext) {
   global $_POST,$JournalDateFmt,$JournalTimeFmt,$MPDateFmt,$MPTimeFmt,
        $MPItemFmt,$MPItemEndFmt,$MPDateTimeFmt,$MultipleItemsPerDay;
   SDV($MPDateFmt,'!!!!%block class=messagehead%$Date');
   SDV($MPTimeFmt,"\n\n>>message<<\n!\$Time by '''\$Author'''!");
   SDV($MPItemFmt,'');
   SDV($MPItemEndFmt,"\n>><<");
   SDV($MPDateTimeFmt,
        ">>message<<\n!!!!%block class=messagehead%\$Date at \$Time by '''\$Author'''\n\n");
   SDV($MultipleItemsPerDay,true);
   $name = @$_POST['author'];
   if (@$_POST['author']=='') $_POST['author'] = 'anon';
   $name = ($name=='') ? 'anonymous' : '[[~' . $name . ']]';
   if ($MultipleItemsPerDay) {
        $date = FormatDateHeading('$Date',$JournalDateFmt,$MPDateFmt);
        $entry = str_replace('$Author',$name,
            FormatDateHeading('$Time',$JournalTimeFmt,$MPTimeFmt));
   } else {
        $date = '';
        $entry = FormatDateHeading('$Date',$JournalDateFmt,
            str_replace('$Author',$name,
            FormatDateHeading('$Time',$JournalTimeFmt,$MPDateTimeFmt)));
   }
   $entry.= $_POST['text'].str_replace('$Author',$name,
            FormatDateHeading('$Time',$JournalTimeFmt,$MPItemFmt)).
            $MPItemEndFmt;
   $order= @$_POST['order'];
   if ($order=='') {
      if (strstr($pagetext,'(:commentbox:)')) {
         $pos = strpos($pagetext,'(:commentbox:)');

         $len = strlen('(:commentbox:)');

         $before = substr($pagetext,0,$pos+$len);

         $after  = substr($pagetext,$pos+$len);
      } else {
         $before = '';
         $after  = $pagetext;
      }

      $entry = "$date$entry";
      $after = ($MultipleItemsPerDay && strstr($after, $date)) ? 
            str_replace($date, $entry, $after) : "\n$entry\n$after";
   } else {
      $entry .= "\n";
      if (strstr($pagetext,'(:commentboxchrono:)')) {
         $pos = strpos($pagetext,'(:commentboxchrono:)');

         $before = substr($pagetext,0,$pos);

         $after  = substr($pagetext,$pos);
      } else {
         $before = $pagetext;
         if ($before[strlen($before)-1]!='\n') $before .="\n";
         $after  = '';
      }

      $before .= ($MultipleItemsPerDay && strstr($before, $date)) ? 
            substr($entry,1) : "\n$date$entry";
   }
   $_POST['text'] = "$before\n$after";
}

# This function not used at present
function IndentParagraphs($text) {
   $text = preg_replace("/\n([[:alnum:]])/","\n->$1",$text);
   $text = preg_replace("/\n([#*])([^#*])/","\n$1$1$2",$text);
   return preg_replace("/\n:([^:]+):/","\n::$1:",$text);
}

?>