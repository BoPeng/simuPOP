<?php if (!defined('PmWiki')) exit();
/*  Copyright 2006 Patrick R. Michaud (pmichaud@pobox.com)
    This file is part of PmWiki; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  See pmwiki.php for full details.
*/

SDV($DraftSuffix, '-Draft');
if ($DraftSuffix) 
  SDV($SearchPatterns['normal']['draft'], "!$DraftSuffix\$!");

##  set up a 'publish' authorization level, defaulting to 'edit' authorization
SDV($DefaultPasswords['publish'], '');
SDV($AuthCascade['publish'], 'edit');
SDV($FmtPV['$PasswdPublish'], 'PasswdVar($pn, "publish")');
if ($AuthCascade['attr'] == 'edit') $AuthCascade['attr'] = 'publish';

## Add a 'publish' page attribute if desired
if (IsEnabled($EnablePublishAttr, 0))
  SDV($PageAttributes['passwdpublish'], '$[Set new publish password:]');

##  with drafts enabled, the 'post' operation requires 'publish' permissions
if ($action == 'edit' && $_POST['post'] && $HandleAuth['edit'] == 'edit')
  $HandleAuth['edit'] = 'publish';

$basename = preg_replace("/$DraftSuffix\$/", '', $pagename);
##  if no -Draft page, switch to $basename
if (!PageExists($pagename) && PageExists($basename)) $pagename = $basename;

##  set edit form button labels to reflect draft prompts
SDVA($InputTags['e_savebutton'], array('value' => ' '.XL('Publish').' '));
SDVA($InputTags['e_saveeditbutton'], array('value' => ' '.XL('Save draft and edit').' '));
SDVA($InputTags['e_savedraftbutton'], array(
    ':html' => "<input type='submit' \$InputFormArgs />",
    'name' => 'postdraft', 'value' => ' '.XL('Save draft').' ',
    'accesskey' => XL('ak_savedraft')));

##  disable the 'publish' button if not authorized to publish
if (!CondAuth($basename, 'publish')) 
  SDVA($InputTags['e_savebutton'], array('disabled' => 'disabled'));

##  add the draft handler into $EditFunctions
if ($action == 'edit') array_unshift($EditFunctions, 'EditDraft');
function EditDraft(&$pagename, &$page, &$new) {
  global $WikiDir, $DraftSuffix, $DeleteKeyPattern, 
    $DraftRecentChangesFmt, $RecentChangesFmt;
  SDV($DeleteKeyPattern, "^\\s*delete\\s*$");
  $basename = preg_replace("/$DraftSuffix\$/", '', $pagename);
  $draftname = $basename . $DraftSuffix;
  if ($_POST['postdraft'] || $_POST['postedit']) $pagename = $draftname; 
  else if ($_POST['post'] && !preg_match("/$DeleteKeyPattern/", $new['text'])) { 
    $pagename = $basename; 
    $page = ReadPage($basename);
    $WikiDir->delete($draftname);
  }
  else if (PageExists($draftname) && $pagename != $draftname)
    { Redirect($draftname, '$PageUrl?action=edit'); exit(); }
  if ($pagename == $draftname && isset($DraftRecentChangesFmt))
    $RecentChangesFmt = $DraftRecentChangesFmt;
}


