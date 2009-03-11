<?php if (!defined('PmWiki')) exit ();

/**
 * This script breaks a page in several editable sections.
 *
 * @package sectionedit
 * @author Patrick R. Michaud <pmichaud@pobox.com>
 * @author John Rankin <john.rankin@affinity.co.nz>
 * @author Sebastian Siedentopf <schlaefer@macnews.de>
 * @author Karl Loncarek <dh2mll@web.de>
 * @author Aidin Abedi <fooguru@msn.com>
 * @version 2.1.7 (2008-02-03)
 * @link http://www.pmwiki.org/wiki/Cookbook/SectionEdit
 * @copyright by the authors 2005
 * @license http://www.gnu.org/copyleft/gpl.html GNU General Public License
 */

/***************************************************************************
For the Admins
==============
Requirements
------------
* Requires PmWiki 2.1.6 and higher
* Last tested on PmWiki 2.2.0-beta45

Install
-------
1. put this script into your cookbook folder
2. include it into your farmconfig.php (when using a Farm-setup:
          include_once("$FarmD/cookbook/sectionedit.php");
   or when using default setup add the following to your config.php:
          include_once("cookbook/sectionedit.php");

Customization
-------------
### Layout
The position and layout of the section links are controlled by the div.sectionedit class in
the script.

### XLPage Strings
* $[(Edit Section &#x2193;)]
* $[Section] x $[of] y

###  Horizontal-line Sectioning
This means that every horizontal line ---- is used for starting a new section.
To enable this option set the variable
    $SectionEditHorzLines = true;
By default the value is false;

### Autosectioning
This means that every heading is used for starting a new section. Which
headings are used depends on the variable $SectionEditAutoDepth.
To only use ==== and no headers for starting sections set the variable
    $SectionEditWithoutHeaders = true;
By default the value is false;


### Depth of Autosectioning
The following variable defines which range of headings is used for starting a
new section. Default is '6'. Here some examples:
    $SectionEditAutoDepth = 1;
means that only the heading "!" is used for a new section.
    $SectionEditAutoDepth = 3;
means that only the headings "!", "!!", and "!!!" are used for a new section.

### Autosectioning as it works at MediaWiki
    $SectionEditMediaWikiStyle = true; //default value
With this setting the automatic sectioning works as in MediaWiki, e.g.
everything between two "!" headings is one section. When set to FALSE the
sections are delimited by headings or the ==== markup.
Examples:
    value:TRUE  (sections)     FALSE (sections)   meanings:
          Text                 Text               Start of section:  \
          !Head1          \    !Head1    \        Middle of section: |
          Text            |    Text      /        End of section:    /
          !!Head2       \ |    !!Head2   \
          Text          | |    Text      /
          !!!Head3    \ | |    !!!Head3  \
          Text        | | |    Text      /
          ====      \ | | |    ====      \
          Text      / | | |    Text      /
          !!!!Head4 \ | | |    !!!!Head4 \
          Text      / / | |    Text      /
          !!!Head5  \   | |    !!!Head5  \
          Text      /   / /    Text      /
          !Head6    \          !Head6    \
          Text      /          Text      /

Usage
=====
To break a page use "====" on your wikipage *before* the wanted section.
Alternatively set $SectionEditAuto or use (:autosections:)

For Developers
==============
What's going on here?
---------------------
1. On browsing the text is searched and all occurences of ==== and headings
   are found and preceded by a (:sectionedit <type> <filename> <reference>:)
   where <type> will be e.g. '====' or e.g.'!!'. <filename> is reserved for
   usage with includes (the file where it is contained) <reference> points to
   the main document. Also headers get a style so the margin before the headers
   that are used for autosectioning are removed.
2. include Markup is changed to first point to our own function which works
   the same way as 1.
3. Markup (:nosections:) (:autosections:) (:noautosections:) and
   (:nosectionsinincludes:) as directives setting some variables
4. Replacement, numbering, check and display of editlinks with Markup
   (:sectionedit <sectionno> <reference>:). editlinks will be
   <filename>?action=edit&s=<sectionno>&auto=<1|>&from=<reference>
5. edit-action is changed to an own function when edit contains an 's'
   parameter. Calculate correct section

ToDo
----

Author/Contributors
-------------------
original (:breakpage:) script
Copyright 2004 Patrick R. Michaud (pmichaud@pobox.com)

modifications to use ____ and edit page chunks
Copyright 2004 John Rankin (john.rankin@affinity.co.nz)

modifications to use ==== to seperate a page in several
editable sections on one single page
Copyright 2005 Sebastian Siedentopf (schlaefer@macnews.de)

some enhancements and bugfixes
Copyright 2005-2006 Karl Loncarek (dh2mll@web.de) <Klonk>

(optional) modifications to also use ---- to seperate a page in several
editable sections on one single page
Copyright 2008 Aidin Abedi (fooguru@msn.com)

Version History
---------------
* 2.1.8 - 2008-03-03 - Klonk
** fixed bug in parameter handling
* 2.1.7 - 2007-05-24 - Aidin
** added $SectionEditHorzLines, (:horzsections:), (:nohorzsections:)
* 2.1.6 - 2007-05-24 - Klonk
** speed optimization
* 2.1.5 - 2007-05-23 - Klonk
** fixed bug in optimization introduced with 2.1.4 (when editing was password protected)
* 2.1.4 - 2007-05-09 - Petko
** optimization of the conditional display
* 2.1.3 - 2006-09-22 - Klonk
** fixed some wrong regex, and some other minor stuff
* 2.1.2 - 2006-09-22 - Klonk
** fixed silly bug... forgot ?= in an regex
** [bugfix] (:...:) markup now still remains is not removed anymore
* 2.1.1 - 2006-09-22 - Klonk
** [bugfix] (:if ...:)====(:ifend:) lead to an empty editwindow
** now anything:)==== works also
* 2.1.0 - 2006-09-21 - Klonk
** section edit links are now only shown, when edit rigths are available
** [bugfix] Part before section could appear several times, when no author is provided
* 2.0.5 - 2006-09-20 - Klonk
** fixed bug in creation of editlinks for includes
** [bugfix] "Save and Edit" replacing the original file
* 2.0.4 - 2006-09-18 - Klonk
** added support for "action=view"
** fixed some stuff for use with PHP5
** should work now also in Auth stuff
* 2.0.3 - 2005-10-15 - Schlaefer
** [bugfix]: include parameter for lines and between anchors now working
** [bugfix]: media wiki style broke after 1.4.1
* 2.0.2 - 2005-10-09 - Schlaefer
** Change: changed & to &amp; in URLs
* 2.0.1 - 2005-09-30 - Klonk
** Change: actual pagename is now determined during the first check for sections
* 2.0 - 2005-09-28 - Klonk
** complete rewrite of browsing part, without programming around given tools
	in PmWiki. This should give greatly enhanced compatibility
** Feature: added (:noautosections:) for disabling automatic sections generation
** Change: provided directives now work also within GroupHeader and GroupFooter
** Change: The core splitting code staid the same, but I did a lot of modularizations for
	easier maintainance
** Feature: With the rewrite e.g. (:nosections:) now also works within conditional
	markup
** Change: Creation of editlinks now happens when links markup is processed
** Change: Now using 'auto' parameter in the URL for detection of section editing
	thus using a different edit handler
** Change: Adopted Edithandler to new parameters
***************************************************************************/

/*** defines the layout of the section edit links */
SDV($HTMLStylesFmt['sectionedit'], "
div.sectionedit { text-align:right;font-size:smaller;clear:both;}
");

/*** setting default values of global variables */
SDV($SectionEditWithoutHeaders, FALSE);
SDV($SectionEditAutoDepth, '6');
SDV($SectionEditMediaWikiStyle, TRUE);
SDV($SectionEditInIncludes, TRUE);
SDV($SectionEditHorzLines, FALSE);
/*** Internal gloabl variable, should not be changed */
$SectionEditDisable = FALSE;
$SectionEditCounter = 0;
$SectionEditSectionCounter = Array();
$SectionEditIncludes = array();
$SectionEditActualPage = '';

/*** markup for setting variables*/
Markup('nosections', 'directives', '/\\(:nosections:\\)/ei', "PZZ(\$GLOBALS['SectionEditDisable'] = true)");
Markup('autosections', '>nosections', '/\\(:autosections\s*(\d)*:\\)/ei', "PZZ(\$GLOBALS['SectionEditWithoutHeaders'] = false)");
Markup('noautosections', '>autosections', '/\\(:noautosections\s*(\d)*:\\)/ei', "PZZ(\$GLOBALS['SectionEditWithoutHeaders'] = true)");
Markup('nosectionsinincludes','directives','/\\(:nosectionsinincludes:\\)/ei', "PZZ(\$GLOBALS['SectionEditInIncludes'] = false)");
Markup('horzsections', '>nosections', '/\\(:horzsections\s*(\d)*:\\)/ei', "PZZ(\$GLOBALS['SectionEditHorzLines'] = true)");
Markup('nohorzsections', '>horzsections', '/\\(:nohorzsections\s*(\d)*:\\)/ei', "PZZ(\$GLOBALS['SectionEditHorzLines'] = false)");

/*** initial markup handling*/
Markup('removesectionmarker', 'block', '/====+/', '');

/*** Replacement for the supplied (:include:) markup*/
Markup('include', '>if',
  '/\\(:include\\s+(\\S.*?):\\)/ei',
  "PRR(SectionEditIncludeText(\$pagename, '$1'))");

/*** Conversion if temporary Markup (:sectionedit <section> <reference>:) to link*/
Markup('sectionedit','links','/\\(:sectionedit\\s+(\\S.*?)\\s+(\\S.*?)\\s+(\\S.*?):\\)/ei',"SectionEditCreateLink('$1','$2','$3')");

/*** Convert temporary markup (:sectionedit...:) into clickable editlinks*/
function SectionEditCreateLink($pagename, $type, $from='') {
	global $SectionEditWithoutHeaders, $SectionEditDisable, $SectionEditSectionCounter,$ScriptUrl;
	/*exit, when disabled or not allowed*/
	if ($SectionEditDisable || ($SectionEditWithoutHeaders && ($type == "!")))
		return '';
	/*set counter to 0 when $pagename and $from differ*/
	if (($pagename != $from) && !isset($SectionEditSectionCounter[$pagename]))
		$SectionEditSectionCounter[$pagename] = 0;
	/*reset counterarray, when  working on included file is done*/
	/*create editlink with anchor*/
	$editlink = FmtPageName("<a name='s$pagename".'_'."$SectionEditSectionCounter[$pagename]'></a><a href='\$PageUrl?action=edit&amp;s=".
		$SectionEditSectionCounter[$pagename]."&amp;auto=".
		($SectionEditWithoutHeaders ? "n" : "y")."&amp;from=$from'>$[(Edit Section &#x2193;)]</a>", $pagename);
	$out = Keep("<div class='sectionedit'>$editlink</div>");
	$SectionEditSectionCounter[$pagename]++;
	return $out;
}

/***Take text, add temporary markup before sections, then give back wikitext (helper function)*/
function SectionEditMarkup($mainpage,$pagename, $text) {
	global $SectionEditWithoutHeaders, $SectionEditAutoDepth, $SectionEditActualPage, $SectionEditHorzLines;
	$out = '';
	/*editing is not allowed anyway thus exit immediately*/
	if (!CondAuth($pagename, 'edit'))
		return ($text);
	/*check whether $pagename contains '#' and remove parts afterwards*/
	$pagenameparts = explode('#',$pagename);
	$pagename = $pagenameparts[0];
	/*make $pagename a valid wikipage name including group*/
	$pagename = MakePageName($mainpage,$pagename);
	/*check whether exisiting sections contain headers or ==== and split them*/
#	$p = preg_split('/((?<=header:\))|(?m:^)|(?:.*:\)))((?=!{1,'.$SectionEditAutoDepth.'}[^!])|(?=====))/', $text);
	$horzline_match = '';
	if ($SectionEditHorzLines) $horzline_match = '|(?=----)';
	$p = preg_split('/((?m:^)|(?<=:\)))((?=!{1,'.$SectionEditAutoDepth.'}[^!])|(?=====)'.$horzline_match.')/', $text);
	/*check for special PHP version and fix strange problem with preg_split*/
	if (phpversion()=='4.1.2') {
		if (substr($p[(count($p)-1)],0,1)=='=')
			$p[(count($p)-1)] = '='.$p[(count($p)-1)];
		else
			$p[(count($p)-1)] = '!'.$p[(count($p)-1)];
	}
	/*Set style to remove the space before the headings*/
	if (!$SectionEditWithoutHeaders)
		$p = preg_replace('/((?m:^)|(?<=:\)))((!{1,'.$SectionEditAutoDepth.'})([^!]))/','$1$3%block margin-top=0px%$4',$p);
	/*creation of temporary markup*/
	for ($i = 0; $i < count($p); $i ++) {
		$horzline = (substr($p[$i],0,1) == '-');
		
		if ($horzline) {
			$out .= '----';
		}
	
		if (($i != 0) || ($pagename != $SectionEditActualPage) || $horzline) {
			/*add temporary link before text, for included files*/
			$out .= "<:block>(:sectionedit $pagename ".substr($p[$i],0,1)." $SectionEditActualPage:)\n";
		}
		
		if ($horzline) {
			$p[$i] = substr($p[$i],4);
		}
		
		$out .= $p[$i];
		
	}
	
	return ($out);
}

/*** Creation of temporary markup in front of the sections for the actual page*/
function SectionEditFirstTime($pagename, $text) {
	global $SectionEditCounter,$SectionEditIncludes,$InclCount;
	global $SectionEditActualPage,$SectionEditSectionCounter;
	/*set name of actual page in global variable*/
	if ($SectionEditActualPage == ''){
		$SectionEditActualPage = $pagename;
		$SectionEditSectionCounter[$pagename] = 1;
	}
	/*allow adding of markup only for the first function call*/
	if ((++$SectionEditCounter > 1) || ($InclCount>0))
		return $text;
	/*modify text, add temporary markup*/
	$out = SectionEditMarkup($pagename,$pagename,$text);
	/*check for includes and add filenames to global variable*/
	SectionEditCheckForIncludes('',$out);
	return $out;
}

/***Check text for includes within and collect pagenames (helper function)*/

function SectionEditCheckForIncludes($pagename,$text) {
	global $SectionEditIncludes;
	/*check whether checked page is in list, then add included pages*/
	if (in_array($pagename,$SectionEditIncludes) || ($pagename =='')) {
		/*get all occurences of (:include...:)*/
		preg_match_all('/\\(:include\\s+(\\S.*?):\\)/',$text,$args);
		/*create list of pagenames*/
		$args = ParseArgs(implode(" ",$args[1]));
		/*create array only containing pagenames*/
		$SectionEditIncludes = array_merge((array)$SectionEditIncludes,(array)$args['']);
		/*remove double page entries*/
		$SectionEditIncludes = array_unique((array)$SectionEditIncludes);
	}
	return;
}

/***Add temporary markup to included text*/
function SectionEditIncludeText($pagename,$inclspec) {
	global $SectionEditInIncludes,$SectionEditIncludes;
	/*exit when Includes should not have sections*/
	if (!$SectionEditInIncludes)
		return IncludeText($pagename,$inclspec);
	/*get all names that should be included*/
	$args = ParseArgs($inclspec);
	/* keeps all parameters but pagenames*/
	$argsparam = '';
	foreach($args as $k => $v) $argsparam .= ($k != '#' && !is_array($v) ) ? " ".$k."=".$v : "";
	/*keep only array with pagenames*/
	$args = $args[''];
	for ($i = 0; $i < count($args); $i ++) {
		/*get text from includes*/
		$text = IncludeText($pagename,$args[$i].$argsparam);
		/*work on text only if present and valid include*/
		if (in_array($args[$i],$SectionEditIncludes) && ($text != '')) {
			/*check for recursive includes, get pagenames and save them*/
			SectionEditCheckForIncludes($args[$i],$text);
			/*add temporary markup to text from includes*/
			$out[] = SectionEditMarkup($pagename,$args[$i], $text);
		}
		else
			$out[] = $text;
	}
	return implode('',$out);
}

/* This function handles the edit, preview and saving of sections.
 * It derived from the standard HandleEdit() function defined in pmwiki.php. */
function HandleEditSection($pagename, $auth = 'edit') {
	global $IsPagePosted, $EditFields, $ChangeSummary, $EditFunctions, $FmtV, $Now, $HandleEditFmt;
	global $PageStartFmt, $PageEditFmt, $PagePreviewFmt, $PageEndFmt, $GroupHeaderFmt, $GroupFooterFmt;
	global $PageEditForm, $EnablePost, $InputTags, $SectionEditWithoutHeaders, $SectionEditHorzLines;
	global $SectionEditMediaWikiStyle, $SectionEditAutoDepth, $MessageFmt;

	// we need some additional values in the edit form for section editing.
	// To respect Site.EditForm we replace the standard PmWiki edit form
	// e_form defined in /scripts/form.php with this
	$InputTags['e_form'] = array (":html" => "<form method='post' action='\$PageUrl?action=edit&amp;s=\$PNum&amp;auto=\$AutoS&amp;from=\$FromP'>
											<input type='hidden' name='action' value='edit' />
											<input type='hidden' name='n' value='\$FullName' />
											<input type='hidden' name='basetime' value='\$EditBaseTime' />
											<input type='hidden' name='prechunk' value=\"\$PreChunk\" />
											<input type='hidden' name='s' value='\$PNum' />
											<input type='hidden' name='auto' value='\$AutoS' />
											<input type='hidden' name='from' value='\$FromP' />
											<input type='hidden' name='postchunk' value=\"\$PostChunk\" />");

	/* standard code from HandleEdit()*/
	if ($_REQUEST['cancel']) {
		Redirect($pagename);
		return;
	}
	Lock(2);
	$IsPagePosted = false;
	$page = RetrieveAuthPage($pagename, $auth, true);
	if (!$page)
		Abort("?cannot edit $pagename");
	PCache($pagename, $page);
	$new = $page;
	/*if a preview previously took place the currently not edited text sections are obtained*/
	$PageChunks = array ('prechunk', 'postchunk');
	foreach ((array) $PageChunks as $c) {
		$$c = '';
		if (@ $_POST[$c]) {
			$$c = str_replace("\r", '', stripmagic($_REQUEST[$c]));
			//remove additonal character again that was added before
			$$c = substr($$c,0,strlen($$c)-1);
		}
	}
	/*combine chunks again, when previewing for correct*/
	if (@ $_POST['text']){
		$new['text'] = str_replace("\r", '', stripmagic($_REQUEST['text']));
		if ($postchunk != '')
			$new['text'] .= "\n";
		$new['text'] = $prechunk.$new['text'].$postchunk;
	}

	/*set variable for autosectioning depending on parameter given, when calling function*/
	if (@ $_GET['auto'] == 'y')
		$SectionEditWithoutHeaders = false;
	else
		$SectionEditWithoutHeaders = true;
	$originalpage = @ $_GET['from'];

	/*disable sectioning when simultaneous edits*/
	if (@!$_POST['basetime'] || !PageExists($pagename) || $_POST['basetime']>=$page['time']) {
		/* splits the page text and sets the currently edited section */
		$horzline_match = '';
		if ($SectionEditHorzLines) $horzline_match = '|(?=----)';
		
		if ($SectionEditWithoutHeaders)
			$p = preg_split('/((?m:^)|(?<=:\)))(?=====)'.$horzline_match.'/', $new['text']);
		else {
			//now check whether exisiting sections contain headers and split them again
			$p = preg_split('/((?m:^)|(?<=:\)))((?=!{1,'.$SectionEditAutoDepth.'}[^!])|(?=====)'.$horzline_match.')/', $new['text']);
		}
		//check for PHP version and fix strange problem
		if (phpversion()=='4.1.2') {
			if (substr($p[(count($p)-1)],0,1)=='=')
				$p[(count($p)-1)] = '='.$p[(count($p)-1)];
			else
				$p[(count($p)-1)] = '!'.$p[(count($p)-1)];
		}

	}

	/*$n holds the number of the section to be edited*/
	$n = @ $_REQUEST['s'];
	if ($n < 0)
		$n = 0;
	if ($n > count($p))
		$n = count($p);

      /* here the merging of all sub-headings is done to achieve autosectioning a la MediaWiki */
	if ($SectionEditMediaWikiStyle && !$SectionEditWithoutHeaders) {
		//check what heading level started this section
            if (preg_match('/(!{1,'.$SectionEditAutoDepth.'})[^!]/',$p[$n],$parts)) {
			for ($i = $n+1; $i < count($p); $i ++) {
				// search in following sections
				if (preg_match('/((?m:^)|(?<=:\)))!{1,'.strlen($parts[1]).'}[^!]/',$p[$i]))
					//if higher or same level exists then exit
					break;
				else {
					//add currently checked section to active section
					$p[$n] .= $p[$i];
				}
			}
			// combine parts of old array into a new array
			$p = array_merge_recursive(array_slice($p,0,$n+1),array_slice($p,$i));
		}
	}

	/* here the section for editing is selected*/
	$new['text'] = $p[$n];
	/*remove last newline except for last section*/
	if ($n<(count($p)-1))
		$new['text'] = substr($p[$n],0,strlen($p[$n])-1);

	/*get fields when previewing */
	foreach ((array) $EditFields as $k)
		if (isset ($_POST[$k]))
			$new[$k] = str_replace("\r", '', stripmagic($_POST[$k]));
	if ($ChangeSummary)
		$new["csum:$Now"] = $ChangeSummary;

	if (@ $_POST['post']) {//the currently not edited sections are added
		//combine string, add newline that was removed before
		$newedit['text'] = $prechunk.$new['text'];
		//add newline to section except for last section
		if ($n<(count($p)-1))
			$newedit['text'] .= "\n";
		$newedit['text'] .= $postchunk;
	}
	elseif (@ $_POST['preview']) { //page header contains info which section is edited
		if ($n > 0)
			$GroupHeaderFmt = '';
		$GroupHeaderFmt .= '<:block>=&gt; ($[Section] '.(@ $_REQUEST['s']+1).' $[of] '.count($p).')(:nl:)';
		if ($n < count($p))
			$GroupFooterFmt = '';
	}
	elseif (@ !$_POST['postedit']) {
		/* if the section is edited, the not edited sections go into $prechunk and
		 * $postchunk and retained here while editing/previewing until saving the whole page
		 */
		for ($i = 0; $i < count($p); $i ++) {
			if ($i < $n) {
				$prechunk .= $p[$i];
			}
			elseif ($i > $n) {
				$postchunk .= $p[$i];
			}
		}

	}

	/* modified code from HandleEdit() */
	$EnablePost &= (@ $_POST['post'] || @ $_POST['postedit']);
	foreach ((array) $EditFunctions as $fn) {
		$newedit = (array)$new;
		if (@ $_POST['postedit'] || @ $_POST['post']) {
			//combine string, add newline that was removed before
			$newedit['text'] = $prechunk.$new['text'];
			//add newline to section except for last section
			if ($n<(count($p)-1))
				$newedit['text'] .= "\n";
			$newedit['text'] .= $postchunk;
		}
		$fn ($pagename, $page, $newedit);
	}

	//add addtional char at end for keeping last newline when saved through form
	$postchunk .= "-";
	$prechunk .= "-";

	Lock(0);
	if ($IsPagePosted && !@ $_POST['postedit']) {
		//jump directly to section that was last edited
		Redirect($originalpage, "\$PageUrl#s".$pagename."_".$n);
		return;
	}
	$FmtV['$DiffClassMinor'] = (@ $_POST['diffclass'] == 'minor') ? "checked='checked'" : '';
	$FmtV['$EditText'] = str_replace('$', '&#036;', htmlspecialchars(@ $new['text'], ENT_NOQUOTES));
	$FmtV['$EditBaseTime'] = $Now;

	/*additional FmtV for this script, stuff that has to be saved for "save and edit" */
	$FmtV['$PreChunk'] = str_replace('"', '&quot;', str_replace('$', '&#036;', htmlspecialchars($prechunk, ENT_NOQUOTES)));
	$FmtV['$PostChunk'] = str_replace('"', '&quot;', str_replace('$', '&#036;', htmlspecialchars($postchunk, ENT_NOQUOTES)));
	$FmtV['$PNum'] = $n;
	$FmtV['$AutoS'] = @ $_GET['auto'];
	$FmtV['$FromP'] = $originalpage;

	/* standard code from HandleEdit() */
	if ($PageEditForm) {
		$form = ReadPage(FmtPageName($PageEditForm, $pagename), READPAGE_CURRENT);
		$FmtV['$EditForm'] = MarkupToHTML($pagename, $form['text']);
	}
	SDV($HandleEditFmt, array (& $PageStartFmt, & $PageEditFmt, & $PageEndFmt));
	PrintFmt($pagename, $HandleEditFmt);
}

/***special Action handling*/
if (($action == 'browse') | ($action == 'view')) {
	$horzline_match = '';
	if ($SectionEditHorzLines) $horzline_match = '|(----+)';

	/* only show Links for section editing when browsing */
	Markup('editsecmarkgen', '<if', '/^.*(\(:nl:\)|header:\\)|\n)((!{1,6})|(====+)'. $horzline_match .').*$/se', "SectionEditFirstTime(\$pagename,PSS('$0'))");
}
elseif (($action == 'edit') && (@ $_REQUEST['auto'])) {
	/* change Edithandler only if a parameter "s" is supplied (s=section number) */
	$HandleActions['edit'] = 'HandleEditSection';
}
