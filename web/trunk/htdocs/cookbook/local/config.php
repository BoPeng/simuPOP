<?php if (!defined('PmWiki')) exit();
##  This is a sample config.php file.  To use this file, copy it to
##  local/config.php, then edit it for whatever customizations you want.
##  Also, be sure to take a look at http://www.pmichaud.com/wiki/Cookbook
##  for more details on the types of customizations that can be added
##  to PmWiki.  

##  $WikiTitle is the name that appears in the browser's title bar.
$WikiTitle = 'the simuPOP cookbook';

##  $ScriptUrl is your preferred URL for accessing wiki pages
##  $PubDirUrl is the URL for the pub directory.
$ScriptUrl = 'http://simupop.sourceforge.net/cookbook/pmwiki.php';
$PubDirUrl = 'http://simupop.sourceforge.net/cookbook/pub';

# This follows suggestions from
#
# http://www.pmwiki.org/wiki/Cookbook/SourceForgeServers
#
session_save_path('/home/groups/s/si/simupop/persistent/cookbook/sessions');
$WorkDir = "/home/groups/s/si/simupop/persistent/cookbook/wiki.d";
$WikiDir = new PageStore('/home/groups/s/si/simupop/persistent/cookbook/wiki.d/$FullName');

$UploadDir = '/home/groups/s/si/simupop/persistent/cookbook/uploads';
$EnableDirectDownload = 0;

$Skin = 'simupop';


##  If you want to use URLs of the form .../pmwiki.php/Group/PageName
##  instead of .../pmwiki.php?p=Group.PageName, try setting
##  $EnablePathInfo below.  Note that this doesn't work in all environments,
##  it depends on your webserver and PHP configuration.  You might also 
##  want to check http://www.pmwiki.org/wiki/Cookbook/CleanUrls more
##  details about this setting and other ways to create nicer-looking urls.
$EnablePathInfo = 1;

## $PageLogoUrl is the URL for a logo image -- you can change this
## to your own logo if you wish.
# $PageLogoUrl = "$PubDirUrl/skins/pmwiki/pmwiki-32.gif";

## You'll probably want to set an administrative password that you
## can use to get into password-protected pages.  Also, by default 
## the "attr" passwords for the PmWiki and Main groups are locked, so
## an admin password is a good way to unlock those.  See PmWiki.Passwords
## and PmWiki.PasswordsAdmin.
# $DefaultPasswords['admin'] = crypt('secret');

include_once("cookbook/pagetoc.php");
$DefaultPasswords['admin'] = '$1$a96bFAxN$D0x6FmtDUb3QGRTUNv1Kv0';
$DefaultPasswords['edit'] = crypt('cookbook');
$DefaultPasswords['upload'] = crypt('cookbook');


## If you're running a publicly available site and allow anyone to
## edit without requiring a password, you probably want to put some
## blocklists in place to avoid wikispam.  See PmWiki.Blocklist.
$EnableBlocklist = 1;                    # enable manual blocklists
$EnableBlocklist = 10;                   # enable automatic blocklists

##  PmWiki comes with graphical user interface buttons for editing;
##  to enable these buttons, set $EnableGUIButtons to 1.  
# $EnableGUIButtons = 1;

##  To enable markup syntax from the Creole common wiki markup language
##  (http://www.wikicreole.org/), include it here:
# include_once('scripts/creole.php');

##  Some sites may want leading spaces on markup lines to indicate
##  "preformatted text blocks", set $EnableWSPre=1 if you want to do
##  this.  Setting it to a higher number increases the number of
##  space characters required on a line to count as "preformatted text".
# $EnableWSPre = 0;                        # PmWiki 2.2.0 default (disabled)
# $EnableWSPre = 1;                        # lines beginning with space are preformatted
# $EnableWSPre = 4;                        # lines with 4 spaces are preformatted

##  If you want uploads enabled on your system, set $EnableUpload=1.
##  You'll also need to set a default upload password, or else set
##  passwords on individual groups and pages.  For more information
##  see PmWiki.UploadsAdmin.
$EnableUpload = 1;                       
# control how upload works
$EnableUploadOverwrite = 0;
# Organize uploads by page
$UploadPrefixFmt = '/$Group/$Name';
# Keep the old versions! Good for me.
$EnableUploadVersions=1;
# limit total uploads to 10000K
$UploadDirQuota = 10240000;
# allow upload of Python file
$UploadExts['py'] = 'text/x-py';

# $DefaultPasswords['upload'] = crypt('secret');

##  Setting $EnableDiag turns on the ?action=diag and ?action=phpinfo
##  actions, which often helps others to remotely troubleshoot 
##  various configuration and execution problems.
# $EnableDiag = 1;                         # enable remote diagnostics

##  By default, PmWiki doesn't allow browsers to cache pages.  Setting
##  $EnableIMSCaching=1; will re-enable browser caches in a somewhat
##  smart manner.  Note that you may want to have caching disabled while
##  adjusting configuration files or layout templates.
# $EnableIMSCaching = 1;                   # allow browser caching

##  Set $SpaceWikiWords if you want WikiWords to automatically 
##  have spaces before each sequence of capital letters.
# $SpaceWikiWords = 1;                     # turn on WikiWord spacing

##  Set $EnableWikiWords if you want to allow WikiWord links.
##  For more options with WikiWords, see scripts/wikiwords.php .
$EnableWikiWords = 0;                      # enable WikiWord links

##  $DiffKeepDays specifies the minimum number of days to keep a page's
##  revision history.  The default is 3650 (approximately 10 years).
# $DiffKeepDays=30;                        # keep page history at least 30 days

## By default, viewers are prevented from seeing the existence
## of read-protected pages in search results and page listings,
## but this can be slow as PmWiki has to check the permissions
## of each page.  Setting $EnablePageListProtect to zero will
## speed things up considerably, but it will also mean that
## viewers may learn of the existence of read-protected pages.
## (It does not enable them to access the contents of the
## pages.)
# $EnablePageListProtect = 0;

##  The refcount.php script enables ?action=refcount, which helps to
##  find missing and orphaned pages.  See PmWiki.RefCount.
# if ($action == 'refcount') include_once('scripts/refcount.php');

##  The feeds.php script enables ?action=rss, ?action=atom, ?action=rdf,
##  and ?action=dc, for generation of syndication feeds in various formats.
# if ($action == 'rss') include_once('scripts/feeds.php');   # RSS 2.0
# if ($action == 'atom') include_once('scripts/feeds.php');  # Atom 1.0
# if ($action == 'dc') include_once('scripts/feeds.php');    # Dublin Core
# if ($action == 'rdf') include_once('scripts/feeds.php');   # RSS 1.0

##  In the 2.2.0-beta series, {$var} page variables are absolute by
##  default, but a future version will make them relative.  This setting
##  sets them out as relative to begin with.  (If you're starting a new
##  site, it's probably best to leave this setting alone.)
$EnableRelativePageVars = 1;

##  By default, pages in the Category group are manually created.
##  Uncomment the following line to have blank category pages
##  automatically created whenever a link to a non-existent
##  category page is saved.  (The page is created only if
##  the author has edit permissions to the Category group.)
# $AutoCreate['/^Category\\./'] = array('ctime' => $Now);

##  PmWiki allows a great deal of flexibility for creating custom markup.
##  To add support for '*bold*' and '~italic~' markup (the single quotes
##  are part of the markup), uncomment the following lines. 
##  (See PmWiki.CustomMarkup and the Cookbook for details and examples.)
# Markup("'~", "inline", "/'~(.*?)~'/", "<i>$1</i>");        # '~italic~'
# Markup("'*", "inline", "/'\\*(.*?)\\*'/", "<b>$1</b>");    # '*bold*'

##  If you want to have to approve links to external sites before they
##  are turned into links, uncomment the line below.  See PmWiki.UrlApprovals.
##  Also, setting $UnapprovedLinkCountMax limits the number of unapproved
##  links that are allowed in a page (useful to control wikispam).
# include_once('scripts/urlapprove.php');
# $UnapprovedLinkCountMax = 10;

##  The following lines make additional editing buttons appear in the
##  edit page for subheadings, lists, tables, etc.
# $GUIButtons['h2'] = array(400, '\\n!! ', '\\n', '$[Heading]',
#                     '$GUIButtonDirUrlFmt/h2.gif"$[Heading]"');
# $GUIButtons['h3'] = array(402, '\\n!!! ', '\\n', '$[Subheading]',
#                     '$GUIButtonDirUrlFmt/h3.gif"$[Subheading]"');
# $GUIButtons['indent'] = array(500, '\\n->', '\\n', '$[Indented text]',
#                     '$GUIButtonDirUrlFmt/indent.gif"$[Indented text]"');
# $GUIButtons['outdent'] = array(510, '\\n-<', '\\n', '$[Hanging indent]',
#                     '$GUIButtonDirUrlFmt/outdent.gif"$[Hanging indent]"');
# $GUIButtons['ol'] = array(520, '\\n# ', '\\n', '$[Ordered list]',
#                     '$GUIButtonDirUrlFmt/ol.gif"$[Ordered (numbered) list]"');
# $GUIButtons['ul'] = array(530, '\\n* ', '\\n', '$[Unordered list]',
#                     '$GUIButtonDirUrlFmt/ul.gif"$[Unordered (bullet) list]"');
# $GUIButtons['hr'] = array(540, '\\n----\\n', '', '',
#                     '$GUIButtonDirUrlFmt/hr.gif"$[Horizontal rule]"');
# $GUIButtons['table'] = array(600,
#                       '||border=1 width=80%\\n||!Hdr ||!Hdr ||!Hdr ||\\n||     ||     ||     ||\\n||     ||     ||     ||\\n', '', '', 
#                     '$GUIButtonDirUrlFmt/table.gif"$[Table]"');


if ("$pagename"=="")
   $pagename = "Main.HomePage";

## Activate the RenamePage recipe.
if ($action == 'rename' || $action == 'postrename' || $action == 'links' ) {
  #include_once("cookbook/renamehelper.php");
  include_once("cookbook/rename.php");
}

## $pagename = ResolvePageName($pagename);
## include_once('cookbook/cmsb-cmsmode.php');
## 
## ## Enable the IP- and content-banning recipe.
## if ($action == 'edit' || $action == 'comment' || $action == 'diff') {
##   include_once("cookbook/cmsb-blocklist.php"); }
## 
## ## Apply a CMS Look for non-authors (ver. 2.1.10 or newer PmWiki default skin).
## if (! @$page['=auth']['edit']) {
##   $HTMLStylesFmt[CMSLook] = '
##  .headnav, #wikicmds, .pagegroup, .footnav { display:none }
##  .pagetitle { margin-top:8px; } div.lastmod { color:#cccccc; } ';
## }
## 

include_once('cookbook/beautifier/beautifier.php');
register_beautifier('cpp');
# register_beautifier('csharp');
# register_beautifier('javascript');
# register_beautifier('lua');
# register_beautifier('php3');
register_beautifier('python');
# register_beautifier('vbdotnet');
# register_beautifier('xml');


include_once("cookbook/totalcounter.php");

# time out session in 30 min
$InactivityTimeout = 1800;
include_once('cookbook/InactiveTimeout.php');

# notify me of changes
#
# NOTIFY IS NOT SUPPORTED BY SOURCEFORGE.NET
#
#$EnableNotify = 1;
#$NotifySquelch = 8*3600;
#$NotifyDelay = 0;
#$NotifyBodyFmt = "All recent posts: "
#    . "http://simupop.sourceforge.net/cookbook/pmwiki.php/\$SiteGroup/AllRecentChanges"
#    . "\n\n\$NotifyItems\n";
#$NotifyItemFmt =
#    '* http://simupop.sourceforge.net/cookbook/{$FullName} . . . $PostTime by {$LastModifiedBy}';
#$NotifyList[] = 'notify=ben.bob@gmail.com';
#$NotifyList[] = 'notify=bpeng@mdanderson.org';
#$NotifyFile = '/home/groups/s/si/simupop/persistent/cookbook/wiki.d/.notifylist';

# RSS feed to get updates
include_once('scripts/feeds.php');

$MarkupCss = true;
SDV($MarkupExtensionsFmt, array("inote abbr `A `. `- `s `: `f -d ... aquo mac '/ '@ '; [^", "q&a A; {|} =| {= revisions ^!! fig :: para lazyweb spaced squo")); 
include_once("cookbook/extendmarkup.php");
include_once("cookbook/break_page.php"); 

$commentboxaccesscode = true;
include_once('cookbook/commentbox.php');
