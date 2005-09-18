#!/usr/bin/perl
$output = 0;
$first = 0;
$file = "";
while(<>)
{
  if ( /^(>>>|\.\.\.) #file/ )
  {
    ($file) = $_ =~ /.*file\s*(.*)\s*/;
    open(FILE, ">$file");
    $first = 1;
    $output = 1;
    next;
  }
  elsif ( /^>>> #end/ )
  {
    $first = 0;
    $output = 0;
    print FILE ">>>  \n";
    close(FILE);
    next;
  }
  elsif (/^(>>>|\.\.\.)\s*#PS / )
  {
    ($pre, $cmd) = $_ =~ /^(>>>|\.\.\.)\s*#PS (.*)$/;
    system("$cmd");
  }
  if( $output )
  {
    if( $first )
    {
      print FILE ">>> ", substr($_, 4);
      $first = 0;
    }
    else
    {
      print FILE "$_";
    }
  }
}
close(FILE);
