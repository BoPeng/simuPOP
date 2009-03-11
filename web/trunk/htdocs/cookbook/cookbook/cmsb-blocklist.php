<?php if (!defined('PmWiki')) exit();
/*
Copyright 2004 Patrick R. Michaud (pmichaud@pobox.com)
Copyright 2006 Hagan Fox
Many significant contributions by Crisses.

This file is cmsb-blocklist.php for PmWiki, you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.  

This script adds blocklisting capabilities to PmWiki.  By default,
it searches through $SiteGroup.Blocklist, Main.Blocklist, and
$SiteGroup.FarmBlocklist to find blocklist items.  You can override
the default Blocklist pages by specifying your own list of one or
more pages in a configuration file.  Here's an example:

   $BlocklistPages[] = 'MySpecialGroup.Blocklist';
   $BlocklistPages[] = 'SomeOtherGroup.Blocklist';

The script blocks by source IP address or by matching string patterns
in the wikitext of a page.

BLOCKING BY IP ADDRESS

For blocking by IP address, the script scans the Blocklist page for IP
addresses of the form

   block:a.b.c.d
   block:a.b.c.*
   block:a.b.*.*
   
and blocks any postings coming from a host matching the IP address or
within the address range.  Someone from a blocked IP address will not
be able to see page histories (?action=diff is ignored).

BLOCKING BY STRING PATTERN

For blocking by string pattern, the script looks for strings to be
excluded from a posting in the form of "block:something".  The string
can contain a space for multiple words (leading and trailing spaces
are ignored), and will match even for partial words unless it's in
double-quotes. For example, if a Blocklist page contains

   block:192.168.123.*
   block:spam.com
   block:buy drugs here
   block:"urge"
   block:aspirin

then any attempt to save a page with text that contains "spam.com",
"buy drugs here", or the whole word "urge" because it's in double-quotes
(so "burger" wouldn't match) will be disallowed.  Anyone trying to edit
he page from an IP address of 192.168.123.0 to 192.168.123.255 will not
be able to save the page no matter what text it contains.

UNBLOCKLIST PAGES

The script also reads Unblocklist pages looking for strings to be prevented
from blocking in the form of "unblock:something".  The default Unblocklist
pages are $SiteGroup.Unblocklist and Main.Unblocklist but you can override
them the same way, with e.g.

   $UnblocklistPages[] = 'MySpecialGroup.Unblocklist';

Any unblock item will be ignored in blocklist pages.  This way a wiki in
a farm can override something that's included in a farm-wide blocklist page.
If an Unblocklist page contains

   unblock:aspirin

Then block:aspirin in a Blocklist page will have no effect.

INSTALLATION AND USAGE

To use the script, simply place it in the cookbook/ directory and add the
following lines to local/config.php:

  # $EnableWhyBlocked = 1;  # Place a why-blocked message on the edit page.
  # $EnableSimpleBlock = 1; # Simple mode; just quietly block.
  if ($action=='edit' || $action=='comment' || $action=='history') {
  include_once("$FarmD/cookbook/blocklist.php"); }

Here are variables you can set prior to including the script:

  $EnableSimpleBlock is a switch to enable "simple blocking", where
   $action=edit or $action=diff is silently switched to to $action=browse.
  $EnableWhyBlocked is a switch for whether or not to report to the
    author the reasons why their post has been blocked.
  $BlocklistPage is the wiki page where the block list will be stored.
  $BlocklistPages is an array of page names containing block lists.
  $UnblocklistPages is an array of page names containing unblock lists.
  $BlocklistMessageFmt is the format for the "you've been blocked"
    statement on the edit screen.
  $BlockIPMessageFmt is the message format preceding the blocked IP
    address in mails and on the edit screen.
  $BlockTextMessageFmt is the message format preceding the blocked text
    content in mails and on the edit screen.
  $EnableCommentMessages is a switch that causes blocked comments to
    send the author to an edit page so the author will see messages.
  $EnableLogBlocked is a switch to turn on logging of blocked posts
    to a wiki page.
  $BlockedLogPageName is name of the log page.  The default page is
    $SiteGroup.Blocklog.  You will almost certainly  want to add a read
    password to the log page once it's created.
  $BlockedLogLinesMax is the maximum number of entries on the log page.

The script sets the some variables you can use to take actions other than
simply blocking the post, such as sending an email message or logging:

  $Blocklisted is the number of "offenses" in the attempted post. 
  $WhyBlockedFmt is an array containing string messages for the reasons
    the post was blocked.  This makes for easy emails to an admin or data
    transfer to a security script.
*/

// Avoid conflict with another blocklist recipe.
if (defined('BLOCKLIST_VERSION')) return;

// Define the blocklist version.
define('BLOCKLIST_VERSION', '2.5.0-cmsb'); // (cmsb stands for "CMS Bundle".)

// Settings you can override in a configuration file.
SDV($EnableSimpleBlock, 0);     // Use a "silent drop".
SDV($EnableWhyBlocked, 0);      // Default is not to report the reason.
SDV($EnableCommentMessages, 0); // Let commenters see a block messages.
SDV($EnableLogBlocked, 0);      // Write a log of blocked posts.
SDV($BlocklistMessageFmt,
  "<h3 class='wikimessage blockmessage'>"
  ."$[This post has been blocked by the administrator]</h3>");
SDV($BlockIPMessageFmt, "$[Remote host matches blocked IP address: ]");
SDV($BlockTextMessageFmt, "Content matches blocked pattern: ");
SDV($BlockedLogPageName, "$SiteGroup.Blocklog");
SDV($BlockedLogLinesMax, 100); // Number of entries to keep.

// Should we even bother?
if (!($action == 'edit' || $action == 'comment' || $action == 'diff')) { return; }

// From a banned IP?  No diff for you!
if ($action == 'diff') { $EnableSimpleBlock = 1; }

// For backward compatibility
if (isset($BlocklistPage)) { $BlocklistPages[] = $BlocklistPage; }

// Set default whitelist pages to $SiteGroup.Unblocklist, Main.Unlocklist,
// if no Unblocklist pages have been specified.
if (!isset($UnblocklistPages)) {
  $UnblocklistPages = array("$SiteGroup.Unblocklist", 'Main.Unblocklist');
}

// Set default Blocklist pages to $SiteGroup.Blocklist, Main.Blocklist,
// and $SiteGroup.FarmBlocklist if no Blocklist pages have been specified.
if (!isset($BlocklistPages)) {
  $BlocklistPages =
    array("$SiteGroup.Blocklist", 'Main.Blocklist', "$SiteGroup.FarmBlocklist");
}

$Blocklisted = 0;
$Unblock_ar = array();
$page_name = ResolvePageName($pagename);

foreach ($UnblocklistPages as $UnblocklistPage) {
  // Make sure the Blocklist page exists.
  if (! PageExists($UnblocklistPage)) { continue; }

  // Get the page contents.
  $pn = FmtPageName($UnblocklistPage, $page_name);
  $page = ReadPage($pn);

  // Build list of unblocked string patterns.
  $text_ar = explode("\n", $page['text']);
  foreach($text_ar as $k => $line) {
    if (strpos($line, 'unblock:') === FALSE) continue;
    $line = trim($line);
    $line_ar = explode('unblock:', $line);
    foreach($line_ar as $term) { 
      $term = trim($term);
      if (!empty($term)) { $Unblock_ar[$term] = $term; }
    }
  }
}

foreach ($BlocklistPages as $BlocklistPage) {
  // Make sure the Blocklist page exists and this page isn't it.
  if (! PageExists($BlocklistPage)
    || $page_name == $BlocklistPage
    || $page_name == $BlockedLogPageName
    || in_array($page_name, $UnblocklistPages))
  {
    continue;
  }

  // Get the page contents.
  $pn = FmtPageName($BlocklistPage, $page_name);
  $page = ReadPage($pn);

  // Check for blocked IP address.
  list($ipa, $ipb, $ipc, $ipd) = explode('.', $_SERVER['REMOTE_ADDR']);
  if (preg_match("/(\\D$ipa\\.$ipb\\.(\\D|$ipc)\\.($ipd)?\\D)/", $page['text'],
    $matchedip))
  {
    $Blocklisted++;
    $WhyBlockedFmt .= $BlockIPMessageFmt.$matchedip[1]."\n";
  }

  // Build list of blocked string patterns.
  $text_ar = explode("\n", $page['text']);
  foreach($text_ar as $k => $line) {
    if (strpos($line, 'block:') === FALSE) continue;
    $line = trim($line);
    $line_ar = explode('block:', $line);
    foreach($line_ar as $term) {
      $term = trim($term);
      if (!empty($term) && !$Unblock_ar[$term]) {
        $Block_ar[$term] = $term;
      }
    }
  }
}

// Look for matches.
if ($Block_ar) { array_walk($Block_ar, 'BlocklistFindBlocked'); }

// Disallow a blocked post.
if (@$Blocklisted) {
  if ($EnableSimpleBlock == 1) { $action = 'browse'; return; }
  if ($action == 'comment' && $EnableCommentMessages == 1) {
    $page = ReadPage($page_name);
      $_POST['text'] = $page['text'];
    $action='edit';
    $CommentBlocked = 1;
  }
  $EnablePost = 0;
  unset($_POST['post']);
  unset($_POST['postattr']);
  unset($_POST['postedit']);
  // Let the poster know they've been blocked.
  $MessagesFmt[] = $BlocklistMessageFmt;
  // Include the reason(s) for the blockage, if enabled.
  if ($EnableWhyBlocked == 1)  {
    $MessagesFmt[] = "<pre class='blocklistmessage'>$WhyBlockedFmt</pre>";
    if (@$CommentBlocked) {
      $MessagesFmt[] = "<tt class='blocklistmessage'>Please go <a
      href='#' onClick='history.go(-1)'>back</a>
      and submit your comment again, but without the blocked content.</tt>";
    }
  }
}

// Write an entry to the log once the page has been sent to the browser.
if ($Blocklisted && $EnableLogBlocked == 1) {
  register_shutdown_function(LogBlocked, getcwd());
}

/**
* Check the posted page content against blocklist.
*
* Matches only on first, not repeated, offenses-per-criteria.
*/
function BlocklistFindBlocked(&$blockitem) {
  global $Blocklisted, $WhyBlockedFmt, $BlockTextMessageFmt;
    if (preg_match('!^".+"$!', $blockitem)) {
      $trimmed = trim($blockitem, '"');
      if (!preg_match("!(^|[[:space:]])($trimmed)([[:punct:]]|[[:space:]]|$)!i",
        @$_POST['text'])) { return; }
    } elseif (stristr(@$_POST['text'], $blockitem) === FALSE) {
      return;
    }
    $Blocklisted++;
    $WhyBlockedFmt .="$BlockTextMessageFmt$blockitem\n";
}

/**
* Write information about the blocked post to the log page.
*
* Adapted from the ActionLog recipe Copyright 2005 by D.Faure
* (dfaure@cpan.org) See http://www.pmwiki.org/wiki/Cookbook/ActionLog
*/
function LogBlocked($dir) {
  global $pagename, $action, $Author, $BlockedLogPageName,
    $BlockedLogSelfExclude, $BlockedLogLineFmt, $BlockedLogLinesMax,
  $WhyBlockedFmt;
  chdir($dir);
  $BlockedAuthor = (!empty($Author)) ? '\\\\'."\n($Author)" : '';
  $WhyBlockedFmt = substr($WhyBlockedFmt, 0, -1);
  $WhyBlockedFmt =
    preg_replace('!'.PHP_EOL.'!', '\\\\\\\\'."\n", $WhyBlockedFmt);
  $pagename = FmtPageName('$Group.$Name', $pagename);
  $logpagename = FmtPageName('$Group.$Name', $BlockedLogPageName);
  $page = ReadPage($logpagename); if(!$page) return;
  Lock(2);
  if(substr($page['text'], -1, 1) != "\n") $page['text'] .= "\n";
  $blocklogtime = time();
  $blocklogtime = strftime("%m-%d\\\\\n%H:%M:%S", $blocklogtime);
  $page['text'] = "|| [-$blocklogtime-] ||[-[[$pagename]]-]"
    ."|| [-{$_SERVER['REMOTE_ADDR']}-][-$BlockedAuthor-] ||"
    ."[-$WhyBlockedFmt-]||\n".$page['text'];
  if (@$BlockedLogLinesMax > 0)
    $page['text'] = implode("\n", array_slice(explode("\n",
      $page['text'], $BlockedLogLinesMax + 1), 0, $BlockedLogLinesMax));
  WritePage($logpagename, $page);
  Lock(0);
}

/* vim: set expandtab tabstop=2 shiftwidth=2: */
