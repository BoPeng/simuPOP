<?php if (!defined('PmWiki')) exit();

/* Extension to PmWiki 2. Copyright Christophe David 2007.
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
 */

$RecipeInfo['InactivityTimeout']['Version'] = '2007-07-28';

SDV($InactivityTimeout, 600);

if (isset($_SESSION['LastActivity']) && (time() - $_SESSION['LastActivity'] > $InactivityTimeout))
   {
   HandleLogoutA($pagename);
   }
else
   {
   $_SESSION['LastActivity'] = time();
   }
?>
