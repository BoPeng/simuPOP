<?php if (!defined('PmWiki')) exit();
/*
cmsmode.php - Adds some CMS Mode features for PmWiki
Copyright 2006 by Hagan Fox
This program is distributed under the terms of the
GNU General Public License, Version 2

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License, Version 2 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

This recipe script adds the following to a PmWiki site:

* Certain changes to the site's behavior when the visitor isn't permitted
  to edit pages.
  - Only essential site-related actions are permitted.  This way, for
    example, ?action=diff and ?action=source have no effect.
  - Only your site's pages will appear in search results.
  - Most PmWiki-related pages will be inaccessible.
  - The rel="nofollow" attribute is removed from external links, which
    will signal search engines to index them.
  - Page modification times are simplified to just a date.

* (All)RecentChanges pages are modified so they're more appropriate to
  be used as public web feeds.
  - RecentChanges pages other than $SiteGroup.AllRecentChanges will
    have a web feed-friendly format, with PmWiki-related pages excluded.
  - A Main.AllRecentChanges page is created that is appropriate to be used
    as a public site-wide web feed page.

Copy the script to the cookbook/ directory and add the following to
your configuration file:

   ## Put the wiki in "CMS Mode".
   $pagename = ResolvePageName($pagename);
   include_once("$FarmD/cookbook/cmsb-cmsmode.php");

Version: 2006-06-05 - Updated method of checking for edit permissions and added
                      'comment' and 'approvesites' as defaault allowed actions.
Version: 2006-06-05 - Added new $CMSAllowedActions setting.
Version: 2006-06-01 - Added documentation.
Version: 2006-03-28 - Initial version.
*/

## Some default settings
SDV($FeedAllRecentChanges, 'Main.AllRecentChanges'); // Public A.R.C. feed
SDV($EnableCMSMode, TRUE);
SDV($CMSAllowedActions, // Actions allowed for non-authors
  array('browse', 'print', 'search', 'edit', 'login', 'rss', 'atom', 'addlink',
    'skin', 'setskin', 'comment', 'approvesites', 'download'));

## Allow the recipe to be disabled in config.php or another recipe.
if ($EnableCMSMode == FALSE) return;

## Create some special restrictions for read-only pages, which is nearly
## all pages if there's a site-wide edit password.
$page = RetrieveAuthPage($pagename, 'read', false, READPAGE_CURRENT); // DEPRECATED
if (! CondAuth($pagename, 'edit')) {
  ##  Allow only essential site-related actions.
  // global $View; SDV($View, 'site');
  if (! in_array($action, $CMSAllowedActions)) { $action='browse'; }
  ## Exclude Certain pages / groups from search results.
  global $SearchPatterns;
  $SearchPatterns['default'][] = '!\\.(All)?Recent(Changes|Uploads)$!';
  $SearchPatterns['default'][] = '!\\.(Group|Print)(Header|Footer)$!';
  $SearchPatterns['default'][] = '!\\.(GroupAttributes|WikiSandbox)$!';
  $SearchPatterns['default'][] = '!^(Test|Site|PmWiki|Profiles)\\.!';
  ## Deny access to certain groups / pages.
  if (preg_match('!^(PmWiki|Test|Profiles)\\.!', $pagename)
    || (preg_match("!^$SiteGroup".'\\.!', $pagename)
      && !preg_match('!\\.(Search|ActionLog)$!', $pagename))
    || preg_match('!\\.(GroupHeader|GroupFooter|GroupAttributes|WikiSandbox|AuthUser)$!', $pagename)) {
    global $DefaultPage; Redirect($DefaultPage);
  }
  ## Lose the rel='nofollow' attribute for external links.
  global $UrlLinkFmt; $UrlLinkFmt =
    "<a class='urllink' href='\$LinkUrl'>\$LinkText</a>";
  ## Format the time string for last-updated notices.
  global $TimeFmt; $TimeFmt = '%B %d, %Y';
}

## Create custom (All)RecentChanges pages to use for web feeds.
if (preg_match('!\\.(GroupHeader|GroupFooter|GroupAttributes|WikiSandbox)$!', $pagename)
    || preg_match('!^(PmWiki|Test|Site|Profiles)\\.!', $pagename)) {
  $RecentChangesFmt['$Group.RecentChanges'] = '';
  $RecentChangesFmt[$FeedAllRecentChanges] = '';
} else {
  $RecentChangesFmt[$FeedAllRecentChanges] =
    '* [[$Group.$Name]]  $[was modified] $CurrentTime. [=$ChangeSummary=]';
  $RecentChangesFmt['Main.AllRecentChanges'] =
    '* [[{$Group}.{$Name}]]  $[was modified] $CurrentTime. [=$ChangeSummary=]';
}

