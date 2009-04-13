<?php if (!defined('PmWiki')) exit ();

/*  copyright 2007-8 Jon Haupt.
    Build on code from swf.php copyright 2004 Patrick R. Michaud
    and from quicktime.php copyright 2006 Sebastian Siedentopf.
    This file is distributed under the terms of the GNU General Public 
    License as published by the Free Software Foundation; either 
    version 2 of the License, or (at your option) any later version.  

    This module enables embedding of Google Video, Vimeo, and YouTube 
    movies into wiki pages. simply use:
    (:googlevideo whatever:)
    (:youtube whatever:)
    (:vimeo whatever:)
    where 'whatever' is the unique number given to your movie.  so if
    the URL to play the video is:
    http://video.google.com/videoplay?docid=--348928491823948192,
    you would do (:googlevideo -348928491823948192:).
    
    Flickr requires additional parameters.  You have to include the
    width, height, id, and secret.  Width and height have defaults 400
    by 300.  If you are viewing the photo embed page for the video, 
    you should be able to extract this information from what is 
    displayed there.

*/
# Version date
$RecipeInfo['SWFSites']['Version'] = '2008-04-09';

Markup('googlevideo', '<img', "/\\(:googlevideo (.*?)\\s*:\\)/e", "ShowGoogleVideo('$1')");
Markup('vimeo', '<img', "/\\(:vimeo (.*?)\\s*:\\)/e", "ShowVimeoVideo('$1')");
Markup('youtube', '<img', "/\\(:youtube (.*?)\\s*:\\)/e", "ShowYouTube('$1')");
Markup('flickrvid', '<img', "/\\(:flickrvid (.*?)\\s*:\\)/e", "ShowFlickrVid('$1')");

function ShowGoogleVideo($url) {

      $out .= "\n<object type='application/x-shockwave-flash' ";
      $out .= "data='http://video.google.com/googleplayer.swf?docId=$url' width='400' height='326' class='VideoPlayback'>";
      $out .= "\n  <param name='movie' value='http://video.google.com/googleplayer.swf?docId=$url'/>";
      $out .= "\n  <param name='allowScriptAccess' value='sameDomain' />";
      $out .= "\n  <param name='quality' value='best' />";
      $out .= "\n  <param name='bgcolor' value='#ffffff' />";
      $out .= "\n  <param name='scale' value='noScale' />";
      $out .= "\n  <param name='salign' value='TL' />";
      $out .= "\n  <param name='wmode' value='transparent' />";
      $out .= "\n  <param name='FlashVars' value='playerMode=embedded' />";
      $out .= "\n</object>";
      return Keep($out);
  }

function ShowVimeoVideo($id) {  

  $out = "\n<object type='application/x-shockwave-flash' ";
  $out .= "width='480' height='360' ";
  $out .= "data='http://vimeo.com/moogaloop.swf?clip_id=$id";
  $out .= "&amp;server=vimeo.com&amp;fullscreen=1&amp;show_title=1";
  $out .= "&amp;show_byline=1&amp;show_portrait=1&amp;color=00ADEF'>";
  $out .= "\n  <param name='quality' value='best' />";
  $out .= "\n  <param name='allowfullscreen' value='true' />";
  $out .= "\n  <param name='scale' value='showAll' />";
  $out .= "\n  <param name='movie' value='http://vimeo.com/moogaloop.swf?clip_id=$id&amp;server=vimeo.com&amp;fullscreen=1&amp;show_title=1&amp;show_byline=1&amp;show_portrait=1&amp;color=00ADEF' />";
  $out .= "\n</object>";
  return Keep($out);

}

function ShowYouTube($url) {

  $out = "\n<object type='application/x-shockwave-flash' style='width:425px; height:350px;' ";
  $out .= "data='http://www.youtube.com/v/$url'>";
  ## transparent so you can see through it if necessary
  $out .= "\n  <param name='wmode' value='transparent' />";
  $out .= "\n  <param name='movie' value='http://www.youtube.com/v/$url' />";
  $out .= "\n</object>";
  return Keep($out);
  }
  
function ShowFlickrVid($input) {

  $defaults = array(
    'width' => '400',
    'height' => '300');
  $args = array_merge($defaults, ParseArgs($input));

  $out = "\n<object type='application/x-shockwave-flash' width='".$args['width']."' height='".$args['height']."' data='http://www.flickr.com/apps/video/stewart.swf?v=1.167' classid='clsid:D27CDB6E-AE6D-11cf-96B8-444553540000'>";
  $out .= "\n<param name='flashvars' value='intl_lang=en-us&amp;photo_secret=".$args['secret']."&amp;photo_id=".$args['id']."'></param>";
  $out .= "\n<param name='movie' value='http://www.flickr.com/apps/video/stewart.swf?v=1.167'></param>";
  $out .= "\n<param name='bgcolor' value='#000000'></param>";
  $out .= "\n<param name='allowFullScreen' value='true'></param>";
  $out .= "\n<embed type='application/x-shockwave-flash' src='http://www.flickr.com/apps/video/stewart.swf?v=1.167' bgcolor='#000000' allowfullscreen='true' flashvars='intl_lang=en-us&amp;photo_secret=".$args['secret']."&amp;photo_id=".$args['id']."' height='".$args['height']."' width='".$args['width']."'></embed>";
  $out .= "\n</object>";
  return Keep($out);
}