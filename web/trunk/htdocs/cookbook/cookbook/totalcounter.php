<?php
if (!defined('PmWiki'))
	exit ();

$RecipeInfo['TotalCounter']['Version']='2008-01-22 - v1.9.1';
/*
	TotalCounter 1.9
	statistic counter for PmWiki
	copyright (c) 2005/2006 Yuri Giuntoli (www.giuntoli.com)

	This PHP script is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation.

	This PHP script is not part of the standard PmWiki distribution.

	0.1 - 23.06.2005
		First version, counts page views and total views.
	0.2 - 20.11.2005
		Added action=totalcounter which displays a page with statistics summary.
	0.3 - 24.11.2005
		Added logging of users, browsers, operating systems, referers and locations.
	0.4 - 28.11.2005
		Optimization of the detection routines.
		Improved detection of the user.
		Added logging of web bots.
	0.5 - 02.12.2005
		Added possibility to blacklist specific items from being logged.
		Modified regex for better referer and location detection.
		Added extended description of location in statistic summary.
	0.6 - 14.12.2005
		Added possibility to DNS lookup the location in case the server doesn't do it automatically.
		Added detection of location when user is sitting behind a proxy server.
		Added possibility to blacklist with regexes for pages, users, referers and locations.
		Listed pages now are link to the actual page.
		Added possibility to assign a password authorization level (edit, admin, etc).
	1.0 - 21.12.2005
		Corrected a bug when the page is the default page.
		Corrected a bug which assigned a browser when pages were crawled by a web bot.
		Optimization of array routines.
		Public release.
	1.1 - 03.01.2006
		Fixed a bug when no bots are present yet.
		Now users work with both UserAuth and AuthUser.
		Added recognition for other popular web bots.
		Added configuration of bars color in the statistics page.
		Added numbers on items (configurable) in the statistics page.
	1.1b - 05.01.2006
		Fixed a bug with empty blacklist array.
		Fixed an alignment problem in the statistics page.
		Fixed a problem which treated Group/Page different from Group.Page.
		Added version display in the statistics page.
	1.1c - 17.01.2006
		Fixed a problem with the markup to work with 2.1.beta20.
	1.2 - 24.01.2006
		Added links to profile pages for the users.
		Reduced locking loop to 5 seconds.
	1.3 - 30.01.2006
		Suppressed the modification to $pagename, now uses internal variable.
		Fixed a bug when remote location is in upper case.
		Changed creation of lock directory to lock file, to prevent problems with some providers.
	1.4 - 31.01.2006
		Optimized the detection of the current page (using ResolvePageName).
		Added statistic count of languages (when used with the MultiLanguage recipe).
	1.4b - 20.02.2006
		Added blacklist support for languages.
		Some fixes about arrays.
	1.5 - 07.03.2006
		Added {$PageViews} page variable.
		Fixed a problem when ResolvePageName function does not exist (earlier versions of PmWiki).
		Fixed a problem with PHP version <4.3.
	1.6 - 11.06.2006
		Florian Xaver:
		 Added os: "DOS"
		 Added browser: "Arachne GPL"
		 Added browser: "Blazer"
		 Changed 'palmos' to 'palm'
		Schlaefer: a daily page counter, a short input field to set the $TotalCounterMaxItems. Changes he mades have a ## comment.
	1.7 - 26.07.2006
		Florian Xaver:
		 Fixed bug, which resets counter. Now there should be no problems
                 with slow servers anymore.


        IMPORTANT: If you get errors on your server, please change creating and deleting
                   of the directory $lock with creating and deleting of a file. This code
                   is commented.
	1.8 - 2007-01-01 - Dave Carver
		Added ($TotalCounterGEOIP) variable.
		Added ($TotalCounterEnableGEOIP) - Set to 1 to use MaxMind's GEOIP Database
		   for country identification. Make sure to turn off Lookup (set to 0).
		Added code to get Location by looking up GEOIP
	        Added code to hopefully fix resets of the file.
	        Added ignore_user_abort(true) to keep file from reseting.
	        Defaults to 'admin' level for viewing of stats.
	        Minor code refactoring to only open the file in write mode when action=browse

        1.8a - 2007-01-21 - Florian Xaver
                Improved/Fixed handling of userlanguage plug-in: (uses $userlang2 instead of $userlang)
                Fixed handling of "File Downloads" (no "." at the filename)
	1.9 - 2007-10-01 - Mateusz Czaplinski
		Added time statistics (last day, last month,...).
		Chmods can be disabled via configuration option.
	1.9.1 - 2008-01-22 - Mateusz Czaplinski
		A fix which tries to ensure that the site won't get locked up by TC's lockfile.
		Added $TotalCounterFile & $TotalCounterLockfile configuration variables.
*/

define(TOTALCOUNTER, '1.9.1');

SDV($TotalCounterAction, 'totalcounter');
SDV($TotalCounterAuthLevel, 'admin');
SDV($TotalCounterMaxItems, 30);
SDV($TotalCounterEnableLookup, 0);
SDV($TotalCounterBarColor, '#5af');
SDV($TotalCounterShowNumbers, 1);
SDV($TotalCounterEnableGEOIP, 1);
SDV($TotalCounterGeoIPData, "GeoIP.dat");
SDV($TotalCounterEnableDownload, 1);
SDV($TotalCounterDownloadManager, ".download.manager");
SDV($TotalCounterEnableChmods, 1);
SDV($TotalCounterFile, "$WorkDir/totalcounter.stat");
SDV($TotalCounterLockfile, "$WorkDir/totalcounter.lock");

SDV($HTMLStylesFmt['TotalCounter'],
	".TCprogress {background-color:$TotalCounterBarColor;height:13px;width:13px;color:#fff}\n".
	"table.TotalCounter td {font-size:x-small;text-align:center}");
	
SDVA($TotalCounterMonthsShort,
	array('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'));

SDV($TotalCounterBlacklist['Pages'], array ());
SDV($TotalCounterBlacklist['Users'], array ());
SDV($TotalCounterBlacklist['Browsers'], array ());
SDV($TotalCounterBlacklist['OSes'], array ());
SDV($TotalCounterBlacklist['Referers'], array ());
SDV($TotalCounterBlacklist['Locations'], array ());
SDV($TotalCounterBlacklist['Bots'], array ());
SDV($TotalCounterBlacklist['Languages'], array ());

## by MateuszCzaplinski
## last day, last week, ... - data & display descriptions
SDVA($TotalCounterTimeBins, array(
	'LastDay' => array(             # LastDay = 24 hours; 1 hour = 60*60sec
		'max'=>24, 'atom'=>60*60,
		'fmt'=> 'date("G",$now-$atom*($maxnr-1-$nr))' ),
	'LastWeek' => array(            # LastWeek = 7 days
		'max'=>7,  'atom'=>24*60*60,
		'fmt'=> 'date("D",$now-$atom*($maxnr-1-$nr))' ),
	'LastMonth' => array(
		'max'=>30, 'atom'=>24*60*60,
		'fmt'=> 'date("j",$now-$atom*($maxnr-1-$nr))' ),
	'LastYear' => array(            # date('n') is the month of the year
		'max'=>12, 'atom'=>'n',
		'fmt'=> '$TotalCounterMonthsShort[(12+(int)date("n",$now)-$maxnr+$nr)%12]' ),
	'LastYears' => array(
		'max'=>30, 'atom'=>'Y',
		'fmt'=> '(int)date("Y",$now)-($maxnr-1-$nr)' )
));

SDVA($HandleActions, array (
	$TotalCounterAction => 'HandleTotalCounter'
));
SDVA($HandleAuth, array (
	$TotalCounterAction => $TotalCounterAuthLevel
));

global $TotalCounter;
if ($TotalCounterMaxItems <= 0)
	$TotalCounterMaxItems = 1;

$file = $TotalCounterFile;
$lock = $TotalCounterLockfile;
$geoIpFile = "$WorkDir/$TotalCounterGeoIPData";
clearstatcache();
ignore_user_abort(true);

#	while (file_exists($lock)) {
#		$st = stat($lock);
#		if ((time()-$st['mtime']) > 15) {
#                        Abort("Too many visitors at the moment...please try later!");
#			break;
#		}
#	}

//------------------------------------------------------------------------------------

if (!function_exists("file_get_contents")) {
	function file_get_contents($filename) {
		if (($contents = file($filename))) {
			$contents = implode('', $contents);
			return $contents;
		} else
			return '';
	}
}

if (function_exists('ResolvePageName')) {
	$tc_pagename = ResolvePageName($pagename);
} else {
	$tc_pagename = str_replace('/', '.', $pagename); /* line changed by Chris Morison 9/3/06 */
}

if ($tc_pagename == '')
	$tc_pagename = "$DefaultGroup.$DefaultName";

if ($action == 'browse') {

	//find users
	if (isset ($AuthId)) {
		$tc_user = $AuthId;
	} else {
		if (isset ($Author)) {
			$tc_user = $Author;
		} else {
			@ session_start();
			if (isset ($_SESSION['authid'])) {
				$tc_user = $_SESSION['authid'][0];
			} else {
				$tc_user = 'Guest (not authenticated)';
			}
		}
	}

	//find web bot
	if (eregi('ia_archiver', $_SERVER['HTTP_USER_AGENT']))
		$tc_bot = 'Alexa';
	elseif (eregi('ask jeeves', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Ask Jeeves';
	elseif (eregi('baiduspider', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Baidu';
	elseif (eregi('libcurl', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'cURL';
	elseif (eregi('gigabot', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Gigablast';
	elseif (eregi('googlebot', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Google';
	elseif (eregi('grub-client', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Grub';
	elseif (eregi('slurp@inktomi.com', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Inktomi';
	elseif (eregi('msnbot', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'MSN';
	elseif (eregi('scooter', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Altavista';
	elseif (eregi('wget', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'wget';
	elseif (eregi('yahoo! slurp', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Yahoo!';
	elseif (eregi('becomebot', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Become';
	elseif (eregi('fast', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'Fast/Alltheweb';
	elseif (eregi('zyborg', $_SERVER['HTTP_USER_AGENT']) || eregi('zealbot', $_SERVER['HTTP_USER_AGENT'])) $tc_bot = 'WiseNut!';

	//not a bot, so find the browser
	elseif (eregi('arachne', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Arachne GPL';
	elseif (eregi('blazer', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Blazer';
	elseif (eregi('opera', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Opera';
	elseif (eregi('webtv', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'WebTV';
	elseif (eregi('camino', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Camino';
	elseif (eregi('netpositive', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'NetPositive';
	elseif (eregi('internet explorer', $_SERVER['HTTP_USER_AGENT']) || eregi('msie', $_SERVER['HTTP_USER_AGENT']) || eregi('mspie', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'MS Internet Explorer';
	elseif (eregi('avant browser', $_SERVER['HTTP_USER_AGENT']) || eregi('advanced browser', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Avant Browser';
	elseif (eregi('galeon', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Galeon';
	elseif (eregi('konqueror', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Konqueror';
	elseif (eregi('icab', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'iCab';
	elseif (eregi('omniweb', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'OmniWeb';
	elseif (eregi('phoenix', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Phoenix';
	elseif (eregi('firebird', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Firebird';
	elseif (eregi('firefox', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Firefox';
	elseif (eregi('minimo', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Minimo';
	elseif (eregi("mozilla", $_SERVER['HTTP_USER_AGENT']) && eregi("rv:[0-9].[0-9][a-b]", $_SERVER['HTTP_USER_AGENT']) && !eregi("netscape", $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Mozilla';
	elseif (eregi("mozilla", $_SERVER['HTTP_USER_AGENT']) && eregi("rv:[0-9].[0-9]", $_SERVER['HTTP_USER_AGENT']) && !eregi("netscape", $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Mozilla';
	elseif (eregi("libwww", $_SERVER['HTTP_USER_AGENT'])) {
		if (eregi("amaya", $_SERVER['HTTP_USER_AGENT'])) {
			$tc_browser = 'Amaya';
		} else {
			$tc_browser = 'Text browser';
		}
	}
	elseif (eregi('safari', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Safari';
	elseif (eregi('elinks', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'ELinks';
	elseif (eregi('offbyone', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Off By One';
	elseif (eregi('playstation portable', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'PlayStation Portable';
	elseif (eregi('netscape', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Netscape';
	elseif (eregi('mozilla', $_SERVER['HTTP_USER_AGENT']) && !eregi("rv:[0-9]\.[0-9]\.[0-9]", $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Firefox';
	elseif (eregi('links', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Links';
	elseif (eregi('ibrowse', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'iBrowse';
	elseif (eregi('w3m', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'w3m';
	elseif (eregi('aweb', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'AWeb';
	elseif (eregi('voyager', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Voyager';
	elseif (eregi('oregano', $_SERVER['HTTP_USER_AGENT'])) $tc_browser = 'Oregano';
	else
		$tc_browser = 'Unknown';

	//find operating system
	if (eregi('linux', $_SERVER['HTTP_USER_AGENT']))
		$tc_os = 'Linux';
	elseif (eregi('irix', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'IRIX';
	elseif (eregi('hp-ux', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'HP-Unix';
	elseif (eregi('os/2', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'OS/2';
	elseif (eregi('beos', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'BeOS';
	elseif (eregi('sunos', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'SunOS';
	elseif (eregi('palm', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'PalmOS';
	elseif (eregi('cygwin', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'Cygwin';
	elseif (eregi('amiga', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'Amiga';
	elseif (eregi('unix', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'Unix';
	elseif (eregi('qnx', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'QNX';
	elseif (eregi('win', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'Windows';
	elseif (eregi('mac', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'Mac';
	elseif (eregi('risc', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'RISC';
	elseif (eregi('dreamcast', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'Dreamcast';
	elseif (eregi('freebsd', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'FreeBSD';
	elseif (eregi('dos', $_SERVER['HTTP_USER_AGENT'])) $tc_os = 'dos';
	else
		$tc_os = 'Unknown';

	//find referrer domain
	preg_match("/^(http:\/\/)?([^\/:]+)/i", $_SERVER['HTTP_REFERER'], $matches);
	$host = $matches[2];
	//if (preg_match("/[^\.\/]+\.*[^\.\/]+$/", $host, $matches) != 0) $host = $matches[0];
	if ($matches[2] != '') {
		$tc_referer = $matches[2];
	} else {
		$tc_referer = 'Unknown';
	}

	//find location
	if ($_SERVER['HTTP_X_FORWARDED_FOR'] != '') {
		if (strstr($_SERVER['HTTP_X_FORWARDED_FOR'], ', ')) {
			$ips = explode(', ', $_SERVER['HTTP_X_FORWARDED_FOR']);
			$thehost = $ips[0];
		} else {
			$thehost = $_SERVER['HTTP_X_FORWARDED_FOR'];
		}
	} else {
		if (strstr($_SERVER['REMOTE_HOST'], ', ')) {
			$ips = explode(', ', $_SERVER['REMOTE_HOST']);
			$thehost = $ips[0];
		} else {
			$thehost = $_SERVER['REMOTE_HOST'];
		}
	}
	if (preg_match("/[^\.0-9]+$/", $thehost, $matches) != 0)
		$loc = $matches[0];
	if ($loc != '') {
		$tc_location = $loc;
	} else {
		if ($TotalCounterEnableLookup == 1) {
			$thehost = @ gethostbyaddr($_SERVER['REMOTE_ADDR']);
			if (preg_match("/[^\.0-9]+$/", $thehost, $matches) != 0)
				$loc = $matches[0];
			if ($loc != '') {
				$tc_location = $loc;
			} else {
				$tc_location = 'Unknown';
			}
		} else
			if ($TotalCounterEnableGeoIP == 1) {
                                include ('geoip/geoip.inc');
				$gi = geoip_open($geoIpFile, GEOIP_STANDARD);
				$tc_location = strtolower(geoip_country_code_by_addr($gi, $_SERVER['REMOTE_ADDR']));
				geoip_close($gi);
			} else {
				$tc_location = 'Unknown';
			}
	}
	if ($tc_location != 'Unknown')
		$tc_location = strtolower($tc_location);
}

//------------------------------------------------------------------------------------

$oldumask = umask(0);
//mkdir($lock, 0777);
#	touch($lock);
#	fixperms($lock);

$downloadfile = "$WorkDir/$TotalCounterDownloadManager";
if (file_exists($downloadfile)) {
	$TotalCounterDownloads = unserialize(file_get_contents($downloadfile));
}

if (file_exists($file)) {
	$TotalCounter = unserialize(file_get_contents($file));
} else {
        touch($file);
	$TotalCounter['Total'] = 0;
	$TotalCounter['Pages'][$tc_pagename] = 0;
}

if (($action == 'browse') && ($tc_pagename != '')) {
	if( dblock($file) ) {
		$TotalCounter = unserialize(file_get_contents($file));
		$TotalCount = ++ $TotalCounter['Total'];

		if (!@ in_array($tc_pagename, $TotalCounterBlacklist['Pages'])) {
			$blacklisted = false;
			if (is_array($TotalCounterBlacklist['Pages']))
				foreach ($TotalCounterBlacklist['Pages'] as $value)
					if (substr($value, 0, 1) == '/')
						if (preg_match($value, $tc_pagename) > 0)
							$blacklisted = true;

			if (!$blacklisted) {
				$PageCount = ++ $TotalCounter['Pages'][$tc_pagename];
				## handles the daily counter
				if ($TotalCounter['PagesTodayDay'][$tc_pagename] == date("%y%m%d"))
					$PageCountToday = ++ $TotalCounter['PagesTodayCounter'][$tc_pagename];
				else {
					$TotalCounter['PagesTodayDay'][$tc_pagename] = date("%y%m%d");
					$TotalCounter['PagesTodayCounter'][$tc_pagename] = 1;
				}
			} else {
				$PageCount = 0;
			}
		}

		if (!@ in_array($tc_user, $TotalCounterBlacklist['Users'])) {
			$blacklisted = false;
			if (is_array($TotalCounterBlacklist['Users']))
				foreach ($TotalCounterBlacklist['Users'] as $value)
					if (substr($value, 0, 1) == '/')
						if (preg_match($value, $tc_user) > 0)
							$blacklisted = true;

			if (!$blacklisted)
				$TotalCounter['Users'][$tc_user]++;
		}

		if (defined('MULTILANGUAGE'))
			if (isset ($userlang2))
				$TotalCounter['Languages'][$userlang2]++;

		if (isset ($tc_browser) && !@ in_array($tc_browser, $TotalCounterBlacklist['Browsers'])) {
			$TotalCounter['Browsers'][$tc_browser]++;
		}
		if (isset ($tc_bot) && !@ in_array($tc_bot, $TotalCounterBlacklist['Bots'])) {
			$TotalCounter['Bots'][$tc_bot]++;
		}
		if (!@ in_array($tc_os, $TotalCounterBlacklist['OSes'])) {
			$TotalCounter['OSes'][$tc_os]++;
		}

		if (!@ in_array($tc_referer, $TotalCounterBlacklist['Referers'])) {
			$blacklisted = false;
			if (is_array($TotalCounterBlacklist['Referers']))
				foreach ($TotalCounterBlacklist['Referers'] as $value)
					if (substr($value, 0, 1) == '/')
						if (preg_match($value, $tc_referer) > 0)
							$blacklisted = true;

			if (!$blacklisted)
				$TotalCounter['Referers'][$tc_referer]++;
		}

		if (!@ in_array($tc_location, $TotalCounterBlacklist['Locations'])) {
			$TotalCounter['Locations'][$tc_location]++;
		}

		if (defined('MULTILANGUAGE'))
			if (!@ in_array($tc_location, $TotalCounterBlacklist['Languages']))
				$TotalCounter['Languages'][$userlang2]++;
				
		## by MateuszCzaplinski
		## last day, last week, ... - collect data
		if (!isset ($tc_bot)) {
			$TCnow = time();
			foreach ($TotalCounterTimeBins as $n=>$a)
				TCbins($n, $a['max'], $a['atom']);
			$TotalCounter['LastTimestamp'] = $TCnow;
		}

		dbexport_unlock($file, serialize($TotalCounter), 'w');
	} else { // could not acquire a lockfile
		// check if the lockfile isn't a stale one, try to delete it if so
		dblock_remove_stale($file);
	}

} else {
	$TotalCount = $TotalCounter['Total'];
	$PageCount = $TotalCounter['Pages'][$tc_pagename];
	## by Schlaefer
	$TotalCounter['PagesTodayDay'][$tc_pagename] == date("%y%m%d") ? $PageCountToday = $TotalCounter['PagesTodayCounter'][$tc_pagename] : $PageCountToday = 1;
}

//rmdir($lock);
#	unlink($lock);
#	umask($oldumask);

//add the {$PageCount} and {$TotalCount} markup
Markup('{$PageCount}', '<{$var}', '/\\{\\$PageCount\\}/e', $PageCount);
Markup('{$TotalCount}', '<{$var}', '/\\{\\$TotalCount\\}/e', $TotalCount);

## by Schlaefer
## adds vars for the input form
Markup('{$PopularPagesItems}', '<{$var}', '/{\\$TotalCounterMaxItems}/', $_REQUEST['TotalCounterMaxItems'] ? $_REQUEST['TotalCounterMaxItems'] : $TotalCounterMaxItems);

//add the {$PageViews} page variable
$FmtPV['$PageViews'] = '$GLOBALS["TotalCounter"]["Pages"][$pagename]';

## by Schlaefer
## add the {$PagesTodayCounter} page variable
$FmtPV['$PageCountToday'] = '$GLOBALS["TotalCounter"]["PagesTodayCounter"][$pagename]';

function HandleTotalCounter($pagename, $auth = 'read') {
	global $Action, $TotalCounter, $TotalCounterMaxItems, $TotalCounterBarColor, $TotalCounterShowNumbers, $TotalCount, $TotalCounterDownloads, $TotalCounterTimeBins, $TotalCounterBinsFmt;
	global $PageStartFmt, $PageEndFmt;

	//$page = RetrieveAuthPage($pagename, $auth, true, READPAGE_CURRENT);
	$page = RetrieveAuthPage($pagename, $auth);
	if (!$page)
		Abort("?you are not permited to perform this action");

	$alllocations = array (
		'localhost' => 'localhost',
		'Unknown' => 'Unknown',

		'com' => 'Commercial',
		'net' => 'Networks',
		'org' => 'Organizations',
		'aero' => 'Aviation',
		'biz' => 'Business organizations',
		'coop' => 'Co-operative organizations',
		'edu' => 'Educational',
		'gov' => 'US Government',
		'info' => 'Info',
		'int' => 'International organizations',
		'mil' => 'US Dept of Defense',
		'museum' => 'Museums',
		'name' => 'Personal',
		'travel' => 'Travelling',

		'ac' => 'Ascension Island',
		'ad' => 'Andorra',
		'ae' => 'United Arab Emirates',
		'af' => 'Afghanistan',
		'ag' => 'Antigua & Barbuda',
		'ai' => 'Anguilla',
		'al' => 'Albania',
		'am' => 'Armenia',
		'an' => 'Netherlands Antilles',
		'ao' => 'Angola',
		'aq' => 'Antarctica',
		'ar' => 'Argentina',
		'as' => 'American Samoa',
		'at' => 'Austria',
		'au' => 'Australia',
		'aw' => 'Aruba',
		'az' => 'Azerbaijan',

		'ba' => 'Bosnia & Herzegovina',
		'bb' => 'Barbados',
		'bd' => 'Bangladesh',
		'be' => 'Belgium',
		'bf' => 'Burkina Faso',
		'bg' => 'Bulgaria',
		'bh' => 'Bahrain',
		'bi' => 'Burundi',
		'bj' => 'Benin',
		'bm' => 'Bermuda',
		'bn' => 'Brunei Darussalam',
		'bo' => 'Bolivia',
		'br' => 'Brazil',
		'bs' => 'Bahamas',
		'bt' => 'Bhutan',
		'bv' => 'Bouvet Island',
		'bw' => 'Botswana',
		'by' => 'Belarus',
		'bz' => 'Belize',

		'ca' => 'Canada',
		'cc' => 'Cocos (Keeling) Islands',
		'cd' => 'Democratic republic of Congo',
		'cf' => 'Central African Republic',
		'cg' => 'Congo',
		'ch' => 'Switzerland',
		'ci' => 'Ivory Coast',
		'ck' => 'Cook Islands',
		'cl' => 'Chile',
		'cm' => 'Cameroon',
		'cn' => 'China',
		'co' => 'Colombia',
		'cr' => 'Costa Rica',
		'cs' => 'Czechoslovakia',
		'cu' => 'Cuba',
		'cv' => 'Cape Verde',
		'cx' => 'Christmas Island',
		'cy' => 'Cyprus',
		'cz' => 'Czech Republic',

		'de' => 'Germany',
		'dj' => 'Djibouti',
		'dk' => 'Denmark',
		'dm' => 'Dominica',
		'do' => 'Dominican Republic',
		'dz' => 'Algeria',

		'ec' => 'Ecuador',
		'ee' => 'Estonia',
		'eg' => 'Egypt',
		'eh' => 'Western Sahara',
		'er' => 'Eritrea',
		'es' => 'Spain',
		'et' => 'Ethiopia',
		'eu' => 'European Union',

		'fi' => 'Finland',
		'fj' => 'Fiji',
		'fk' => 'Falkland Islands',
		'fm' => 'Micronesia',
		'fo' => 'Faroe Islands',
		'fr' => 'France',

		'ga' => 'Gabon',
		'gb' => 'United Kingdom',
		'gd' => 'Grenada',
		'ge' => 'Georgia',
		'gf' => 'French Guiana',
		'gg' => 'Guernsey',
		'gh' => 'Ghana',
		'gi' => 'Gibraltar',
		'gl' => 'Greenland',
		'gm' => 'Gambia',
		'gn' => 'Guinea',
		'gp' => 'Guadeloupe',
		'gq' => 'Equatorial Guinea',
		'gr' => 'Greece',
		'gs' => 'South Georgia & South Sandwich Islands',
		'gt' => 'Guatemala',
		'gu' => 'Guam',
		'gw' => 'Guinea-Bissau',
		'gy' => 'Guyana',

		'hk' => 'Hong Kong',
		'hm' => 'Heard & McDonald Islands',
		'hn' => 'Honduras',
		'hr' => 'Croatia',
		'ht' => 'Haiti',
		'hu' => 'Hungary',

		'id' => 'Indonesia',
		'ie' => 'Ireland',
		'il' => 'Israel',
		'im' => 'Isle of Man',
		'in' => 'India',
		'io' => 'British Indian Ocean Territory',
		'iq' => 'Iraq',
		'ir' => 'Iran',
		'is' => 'Iceland',
		'it' => 'Italy',

		'je' => 'Jersey',
		'jm' => 'Jamaica',
		'jo' => 'Jordan',
		'jp' => 'Japan',

		'ke' => 'Kenya',
		'kg' => 'Kyrgyzstan',
		'kh' => 'Cambodia',
		'ki' => 'Kiribati',
		'km' => 'Comoros',
		'kn' => 'Saint Kitts & Nevis',
		'kp',
		'North Korea',
		'kr' => 'South Korea',
		'kw' => 'Kuwait',
		'ky' => 'Cayman Islands',
		'kz' => 'Kazakhstan',

		'la' => 'Laos',
		'lb' => 'Lebanon',
		'lc' => 'Saint Lucia',
		'li' => 'Liechtenstein',
		'lk' => 'Sri Lanka',
		'lr' => 'Liberia',
		'ls' => 'Lesotho',
		'lt' => 'Lithuania',
		'lu' => 'Luxembourg',
		'lv' => 'Latvia',
		'ly' => 'Libyan Arab Jamahiriya',

		'ma' => 'Morocco',
		'mc' => 'Monaco',
		'md' => 'Moldova',
		'mg' => 'Madagascar',
		'mh' => 'Marshall Islands',
		'mk' => 'Macedonia',
		'ml' => 'Mali',
		'mm' => 'Myanmar',
		'mn' => 'Mongolia',
		'mo' => 'Macau',
		'mp' => 'Northern Mariana Islands',
		'mq' => 'Martinique',
		'mr' => 'Mauritania',
		'ms' => 'Montserrat',
		'mt' => 'Malta',
		'mu' => 'Mauritius',
		'mv' => 'Maldives',
		'mw' => 'Malawi',
		'mx' => 'Mexico',
		'my' => 'Malaysia',
		'mz' => 'Mozambique',

		'na' => 'Namibia',
		'nc' => 'New Caledonia',
		'ne' => 'Niger',
		'nf' => 'Norfolk Island',
		'ng' => 'Nigeria',
		'ni' => 'Nicaragua',
		'nl' => 'The Netherlands',
		'no' => 'Norway',
		'np' => 'Nepal',
		'nr' => 'Nauru',
		'nu' => 'Niue',
		'nz' => 'New Zealand',

		'om' => 'Oman',

		'pa' => 'Panama',
		'pe' => 'Peru',
		'pf' => 'French Polynesia',
		'pg' => 'Papua New Guinea',
		'ph' => 'Philippines',
		'pk' => 'Pakistan',
		'pl' => 'Poland',
		'pm' => 'St. Pierre & Miquelon',
		'pn' => 'Pitcairn',
		'pr' => 'Puerto Rico',
		'ps' => 'Palestine',
		'pt' => 'Portugal',
		'pw' => 'Palau',
		'py' => 'Paraguay',

		'qa' => 'Qatar',

		're' => 'Reunion',
		'ro' => 'Romania',
		'ru' => 'Russia',
		'rw' => 'Rwanda',

		'sa' => 'Saudi Arabia',
		'sb' => 'Solomon Islands',
		'sc' => 'Seychelles',
		'sd' => 'Sudan',
		'se' => 'Sweden',
		'sg' => 'Singapore',
		'sh' => 'St. Helena',
		'si' => 'Slovenia',
		'sj' => 'Svalbard & Jan Mayen Islands',
		'sk' => 'Slovakia',
		'sl' => 'Sierra Leone',
		'sm' => 'San Marino',
		'sn' => 'Senegal',
		'so' => 'Somalia',
		'sr' => 'Surinam',
		'st' => 'Sao Tome & Principe',
		'su' => 'USSR',
		'sv' => 'El Salvador',
		'sy' => 'Syrian Arab Republic',
		'sz' => 'Swaziland',

		'tc' => 'The Turks & Caicos Islands',
		'td' => 'Chad',
		'tf' => 'French Southern Territories',
		'tg' => 'Togo',
		'th' => 'Thailand',
		'tj' => 'Tajikistan',
		'tk' => 'Tokelau',
		'tm' => 'Turkmenistan',
		'tn' => 'Tunisia',
		'to' => 'Tonga',
		'tp' => 'East Timor',
		'tr' => 'Turkey',
		'tt' => 'Trinidad & Tobago',
		'tv' => 'Tuvalu',
		'tw' => 'Taiwan',
		'tz' => 'Tanzania',
		'ua' => 'Ukraine',
		'ug' => 'Uganda',
		'uk' => 'United Kingdom',
		'um' => 'United States Minor Outlying Islands',
		'us' => 'United States',
		'uy' => 'Uruguay',
		'uz' => 'Uzbekistan',

		'va' => 'Vatican City',
		'vc' => 'Saint Vincent & the Grenadines',
		've' => 'Venezuela',
		'vg' => 'British Virgin Islands',
		'vi' => 'US Virgin Islands',
		'vn' => 'Vietnam',
		'vu' => 'Vanuatu',

		'wf' => 'Wallis & Futuna Islands',
		'ws' => 'Samoa',

		'ye' => 'Yemen',
		'yt' => 'Mayotte',
		'yu' => 'Yugoslavia',

		'za' => 'South Africa',
		'zm' => 'Zambia',
		'zr' => 'Zaire',
		'zw' => 'Zimbabwe',
		
	);

	$Action = 'TotalCounter statistics';

	## by Schlaefer
	## sets the max items if provided by the form
	if ($_REQUEST['TotalCounterMaxItems'])
		$TotalCounterMaxItems = $_REQUEST['TotalCounterMaxItems'];

	//------------------------------------------------------------------------------------------------------------
	// PAGES

	$html = '<h1>Total Counter $[statistics]</h1>' .
	'<br /><hr />' .
	'<h2>$[Page views]</h2>' .
	'<table border=\'0\'>' .
	'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Pages]&nbsp;</b></td><td colspan=\'2\'><b>$[Percent]</b></td><td align=\'right\'><b>$[Count]</b></td></tr>';

	@ arsort($TotalCounter['Pages']);
	$tar = @ array_slice($TotalCounter['Pages'], 0, $TotalCounterMaxItems);
	$tot = $TotalCount;
	$max = @ current($tar);

	$i = 0;
	if (is_array($tar) && $tot) // by Florian Xaver
		foreach ($tar as $pn => $cnt) {
			$html .= '<tr>' .
			 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
			"<td><a href='\$ScriptUrl/$pn'>$pn</a>&nbsp;</td><td>" . Round(100 * $cnt / $tot) . "%</td><td><div style='background-color:$TotalCounterBarColor;height:13px;width:" . Round(200 * $cnt / $max) . "px;color:#fff'></div></td><td align='right'>&nbsp;$cnt</td>" .
			'</tr>';
		}

	## by Schlaefer

	//------------------------------------------------------------------------------------------------------------
	## PAGES daily

	$html .= '</table>' .
	'<br /><hr />' .
	'<h2>$[Page views] $[today]</h2>' .
	'<table border=\'0\'>' .
	'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Pages]&nbsp;</b></td><td colspan=\'2\'><b>$[Percent]</b></td><td align=\'right\'><b>$[Count]</b></td></tr>';
	$tar = array ();
	foreach ($TotalCounter['PagesTodayCounter'] as $pn => $cnt) {
		if ($TotalCounter['PagesTodayDay'][$pn] === date("%y%m%d"))
			$tar[$pn] = $cnt;
	}
	@ arsort($tar);
	$tot = @ array_sum($tar);
	$tar = @ array_slice($tar, 0, $TotalCounterMaxItems);
	$max = @ current($tar);

	$i = 0;
	if (is_array($tar))
		foreach ($tar as $pn => $cnt) {
			$html .= '<tr>' .
			 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
			"<td><a href='\$ScriptUrl/$pn'>$pn</a>&nbsp;</td><td>" . Round(100 * $cnt / $tot) . "%</td><td><div style='background-color:$TotalCounterBarColor;height:13px;width:" . Round(200 * $cnt / $max) . "px;color:#fff'></div></td><td align='right'>&nbsp;$cnt</td>" .
			'</tr>';
			if ($i == $TotalCounterMaxItems)
				break;
		}

	//------------------------------------------------------------------------------------------------------------
	// USERS

	$html .= '</table>' .
	'<br /><hr />' .
	'<h2>$[Users]</h2>' .
	'<table border=\'0\'>' .
	'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Users]&nbsp;</b></td><td colspan=\'2\'><b>$[Percent]</b></td><td align=\'right\'><b>$[Count]</b></td></tr>';

	@ arsort($TotalCounter['Users']);
	$tar = @ array_slice($TotalCounter['Users'], 0, $TotalCounterMaxItems);
	$max = @ current($tar);
	$tot = @ array_sum($tar);

	$i = 0;
	if (is_array($tar))
		foreach ($tar as $pn => $cnt) {
			$html .= '<tr>' .
			 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
			'<td>' .
			 ($pn != 'Guest (not authenticated)' ? "<a href='\$ScriptUrl/\$AuthorGroup/$pn'>$pn</a>" : $pn) .
			"&nbsp;</td><td>" . Round(100 * $cnt / $tot) . "%</td><td><div style='background-color:$TotalCounterBarColor;height:13px;width:" . Round(200 * $cnt / $max) . "px;color:#fff'></div></td><td align='right'>&nbsp;$cnt</td>" .
			'</tr>';
		}

	//------------------------------------------------------------------------------------------------------------
	// LANGUAGES

	if (defined('MULTILANGUAGE')) {
		$html .= '</table>' .
		'<br /><hr />' .
		'<h2>$[Languages]</h2>' .
		'<table border=\'0\'>' .
		'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Languages]&nbsp;</b></td><td colspan=\'2\'><b>$[Percent]</b></td><td align=\'right\'><b>$[Count]</b></td></tr>';

		@ arsort($TotalCounter['Languages']);
		$tar = @ array_slice($TotalCounter['Languages'], 0, $TotalCounterMaxItems);
		$max = @ current($tar);
		$tot = @ array_sum($tar);

		$i = 0;
		if (is_array($tar))
			foreach ($tar as $pn => $cnt) {
				$html .= '<tr>' .
				 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
				"<td>$pn&nbsp;</td><td>" . Round(100 * $cnt / $tot) . "%</td><td><div style='background-color:$TotalCounterBarColor;height:13px;width:" . Round(200 * $cnt / $max) . "px;color:#fff'></div></td><td align='right'>&nbsp;$cnt</td>" .
				'</tr>';
			}
	}

	//------------------------------------------------------------------------------------------------------------
	// BROWSERS

	$html .= '</table>' .
	'<br /><hr />' .
	'<h2>$[Browsers]</h2>' .
	'<table border=\'0\'>' .
	'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Browsers]&nbsp;</b></td><td colspan=\'2\'><b>$[Percent]</b></td><td align=\'right\'><b>$[Count]</b></td></tr>';

	@ arsort($TotalCounter['Browsers']);
	$tar = @ array_slice($TotalCounter['Browsers'], 0, $TotalCounterMaxItems);
	$max = @ current($tar);
	$tot = @ array_sum($tar);

	$i = 0;
	if (is_array($tar))
		foreach ($tar as $pn => $cnt) {
			$html .= '<tr>' .
			 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
			"<td>$pn&nbsp;</td><td>" . Round(100 * $cnt / $tot) . "%</td><td><div style='background-color:$TotalCounterBarColor;height:13px;width:" . Round(200 * $cnt / $max) . "px;color:#fff'></div></td><td align='right'>&nbsp;$cnt</td>" .
			'</tr>';
		}

	//------------------------------------------------------------------------------------------------------------
	// OPERATING SYSTEMS

	$html .= '</table>' .
	'<br /><hr />' .
	'<h2>$[Operating systems]</h2>' .
	'<table border=\'0\'>' .
	'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Operating systems]&nbsp;</b></td><td colspan=\'2\'><b>$[Percent]</b></td><td align=\'right\'><b>$[Count]</b></td></tr>';

	@ arsort($TotalCounter['OSes']);
	$tar = @ array_slice($TotalCounter['OSes'], 0, $TotalCounterMaxItems);
	$max = @ current($tar);
	$tot = @ array_sum($tar);

	$i = 0;
	if (is_array($tar))
		foreach ($tar as $pn => $cnt) {
			$html .= '<tr>' .
			 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
			"<td>$pn&nbsp;</td><td>" . Round(100 * $cnt / $tot) . "%</td><td><div style='background-color:$TotalCounterBarColor;height:13px;width:" . Round(200 * $cnt / $max) . "px;color:#fff'></div></td><td align='right'>&nbsp;$cnt</td>" .
			'</tr>';
		}

	//------------------------------------------------------------------------------------------------------------
	// REFERERS

	$html .= '</table>' .
	'<br /><hr />' .
	'<h2>$[Referers]</h2>' .
	'<table border=\'0\'>' .
	'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Referers]&nbsp;</b></td><td colspan=\'2\'><b>$[Percent]</b></td><td align=\'right\'><b>$[Count]</b></td></tr>';

	@ arsort($TotalCounter['Referers']);
	$tar = @ array_slice($TotalCounter['Referers'], 0, $TotalCounterMaxItems);
	$max = @ current($tar);
	$tot = @ array_sum($tar);

	$i = 0;
	if (is_array($tar))
		foreach ($tar as $pn => $cnt) {
			$html .= '<tr>' .
			 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
			"<td>$pn&nbsp;</td><td>" . Round(100 * $cnt / $tot) . "%</td><td><div style='background-color:$TotalCounterBarColor;height:13px;width:" . Round(200 * $cnt / $max) . "px;color:#fff'></div></td><td align='right'>&nbsp;$cnt</td>" .
			'</tr>';
		}

	//------------------------------------------------------------------------------------------------------------
	// LOCATIONS

	$html .= '</table>' .
	'<br /><hr />' .
	'<h2>$[Locations]</h2>' .
	'<table border=\'0\'>' .
	'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Locations]&nbsp;</b></td><td colspan=\'2\'><b>$[Percent]</b></td><td align=\'right\'><b>$[Count]</b></td></tr>';

	@ arsort($TotalCounter['Locations']);
	$tar = @ array_slice($TotalCounter['Locations'], 0, $TotalCounterMaxItems);
	$max = @ current($tar);
	$tot = @ array_sum($tar);

	$i = 0;
	if (is_array($tar))
		foreach ($tar as $pn => $cnt) {
			$html .= '<tr>' .
			 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
			'<td>' . ($alllocations[$pn] == '' ? 'Unknown' : $alllocations[$pn]) . ' ' .
			 ($pn == 'Unknown' || $pn == 'localhost' ? '' : "(.$pn)") . '&nbsp;</td>' .
			'<td>' . Round(100 * $cnt / $tot) . '%</td>' .
			'<td><div style=\'background-color:$TotalCounterBarColor;height:13px;width:' . Round(200 * $cnt / $max) . "px;color:#fff'></div></td><td align='right'>&nbsp;$cnt</td>" .
			'</tr>';
		}

	//------------------------------------------------------------------------------------------------------------
	// WEB BOTS

	$html .= '</table>' .
	'<br /><hr />' .
	'<h2>$[Web bots]</h2>' .
	'<table border=\'0\'>' .
	'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Web bots]&nbsp;</b></td><td colspan=\'2\'><b>$[Percent]</b></td><td align=\'right\'><b>$[Count]</b></td></tr>';

	@ arsort($TotalCounter['Bots']);
	$tar = @ array_slice($TotalCounter['Bots'], 0, $TotalCounterMaxItems);
	$max = @ current($tar);
	$tot = @ array_sum($tar);

	$i = 0;
	if (is_array($tar))
		foreach ($tar as $pn => $cnt) {
			$html .= '<tr>' .
			 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
			"<td>$pn&nbsp;</td><td>" . Round(100 * $cnt / $tot) . "%</td><td><div style='background-color:$TotalCounterBarColor;height:13px;width:" . Round(200 * $cnt / $max) . "px;color:#fff'></div></td><td align='right'>&nbsp;$cnt</td>" .
			'</tr>';
		}


//------------------------------------------------------------------------------------------------------------
	// Downloads

	$html .= '</table>' .
	'<br /><hr />' .
	'<h2>$[File Downloads]</h2>' .
	'<table border=\'0\'>' .
	'<tr><td' . ($TotalCounterShowNumbers ? ' colspan=\'2\'' : '') . '><b>$[Downloads]&nbsp;</b></td><td colspan=\'2\'></td><td align=\'right\'><b>$[Count]</b></td></tr>';

	@ arsort($TotalCounterDownloads);
	$max = count($TotalCounterDownloads);
	$tot = @ array_sum($TotalCounterDownloads);

	$i = 0;
	if (is_array($TotalCounterDownloads)) {

		for ($row = 0; $row < $max; $row++) {
			$tablerow = each($TotalCounterDownloads);
			$value = $tablerow['value'];
			$html .= '<tr>' .
			 ($TotalCounterShowNumbers ? '<td align=\'right\' valign=\'bottom\'><small>' . ++ $i . '.</small></td>' : '') .
			'<td>' . $tablerow['key'] . '</td>' .
			'<td></td>' .
			'<td></td><td align="right">&nbsp;' . $value . '</td>' .
			'</tr>';
		}
	}

//------------------------------------------------------------------------------------------------------------
	// Time statistics
	## by MateuszCzaplinski
	
	foreach( $TotalCounterTimeBins as $n=>$a ) {
		$html .= '</table>' .
		'<br /><hr />' .
		"<h2>$[$n]</h2><table border='0' class='TC-$n TotalCounter'>";

		SDVA($TotalCounterBinsFmt,array(
			'"$count"',
			'"<div class=\"TCprogress\" style=\"$direction:".Round(1+200*$count/$maxcount)."px;\"></div>"',
			'date("G",$now-$atom*($maxnr-1-$nr))' ));
		$fmt = $a['fmt'];
		if( is_string($fmt) ) {
			$tmp = $fmt;
			$fmt = $TotalCounterBinsFmt;
			$fmt[2] = $tmp;
		}
		if( is_array($fmt) ) {
			$rows = array();
			## Variables used in 'fmt'
			$maxcount = @max( $TotalCounter[$n] );
			$direction = 'height';
			$maxnr = $a['max'];
			$atom = $a['atom'];
			$now = time();
			for( $nr=0; $nr<$a['max']; $nr++ ) {
				for( $j=0; $j<count($fmt); $j++ ) {
					$count = $TotalCounter[$n][$nr];
					$rows[$j] = (string)@$rows[$j] .
						"<td valign='bottom' class='seq$j'>".
						(string)eval("global \$TotalCounterMonthsShort; return ({$fmt[$j]});").
						"</td>\n";
				}
			}
		}
		$html .= '<tr>'.implode('</tr><tr>',$rows).'</tr>';
	}
	
		
	$html .= '</table><hr /><p align=\'right\'>TotalCounter v' . TOTALCOUNTER . '</p>';

	PrintFmt($pagename, array (
		& $PageStartFmt,
		$html,
		& $PageEndFmt
	));
}

## by MateuszCzaplinski
## Manages an array of counters, each for a specified time interval.
## In $TotalCounter[$name] array there are $max counters. Each counter
## is for time interval of $atom length.
## Note: if $atom is a number, it is a length of interval measured
## in seconds. If $atom is a string, it means date($atom) is executed 
## and the result is the index of an interval.
## NOTE: See TODO below
function TCbins($name,$max,$atom) {
	global $TotalCounter, $TCnow;
	$last = $TotalCounter['LastTimestamp'];
	if( $TCnow < $last ) return; // some error?
	if( !$last ) $TotalCounter[$name] = array_fill(0,$max,0);
	if( is_string($atom) ) {
		$diff = (int)date($atom,$TCnow) - (int)date($atom,$last);
		if( $diff < 0 ) $diff += $max;
		# TODO: handle time delta > $max
		# Until fixed, if the site has no visitor for about a
		# year, statistics will get falsified (empty years will compress)
	}
	else
		$diff = (int)($TCnow/$atom) - (int)($last/$atom);
		
	if( $diff < 0 ) return; 
	if( $diff > 0 ) {
		$a = @array_slice($TotalCounter[$name], $diff, max(0,$max-$diff));
		if(!$a) $a = array();
		$a = array_pad($a, $max, 0);
		$TotalCounter[$name] = $a;
	}
	$TotalCounter[$name][$max-1]++; // put our visit in last bin
}

//Windows 9x FLOCK Alternative - Chozo4
//Passes multiple instance stress testing
//http://mechresource.myvnc.com/board
//Modified (breaks and returns 0 on failure,
// or returns 1 on success) by Mateusz Czaplinski, 22.01.2008

function aquirelock($wp) {
	//Check if lock doesn't exist or our target is unwritable
	if(file_exists("$wp.l") || !is_writable($wp))
		return 0;

	//create the lock - hide warnings and pass empty if already created from racing
	return @ fopen("$wp.l", 'x');
}

function dblock($wp) {
	global $TotalCounterEnableChmods;
	//Check for lockfile handle - if empty , another process raced the lock so report a failure
	$ftw = aquirelock($wp);
	if( !$ftw )
		return 0;

	if($TotalCounterEnableChmods) chmod($wp, 0444); //set the target file to read-only
	fwrite($ftw, 'lock'); //write the lockfile with 4bytes
	if($TotalCounterEnableChmods) chmod("$wp.l", 0444); //set the lockfile to read only (OPTIONAL)
	fclose($ftw); //close our lockfile
	clearstatcache(); //Clear the stat cache
	return 1;
}

// Note: don't call it if 'dblock()' returned 0 !
function dbexport_unlock($wp, $data, $meth) {
	global $TotalCounterEnableChmods;
	if($TotalCounterEnableChmods) chmod($wp, 0666); //Set the target file to read+write

	//Write the passed string to the target file then close
	fwrite($ftw = fopen($wp, $meth), $data);
	fclose($ftw);

	//Validate the written data ujsing a string comparison
	$check = file_get_contents($wp);
	if ($check != $data)
		echo "Data Mismatch - Locking FAILED!<br>";

	chmod("$wp.l", 0666); //Set the lockfile to read+write (OPTIONAL)
	unlink("$wp.l"); //Release the lockfile by removing it
}

function dblock_remove_stale($wp) {
	$t=filemtime("$wp.l");
	// 75 minutes - to make absolutely sure we're not tricked by Daylight
	// Savings on Windows - see http://www.php.net/manual/en/function.stat.php#58404
	if( $t+(75*60) < time())
		@unlink("$wp.l");
}
?>
