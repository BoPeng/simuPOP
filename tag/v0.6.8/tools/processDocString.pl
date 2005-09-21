#!/usr/bin/perl

$text = "";
$valid = 1;
$ignore = 0;

$alltext = "";

@cont = ("","","","","","","","","","" );
$head = 0;
$brief = 1;
$usage = 2;
$argument = 3;
$detail = 4;
$return = 5;
$note = 6;
$see = 7;
$example = 8;
$end = 9;

$ArgString = "";

@Pattern = (
  '^.*feature\("docstring"\)', 
  "^Description:",
  "^Usage:",
  "^Arguments:",
  "^Details:",
  "^Value:",
  "^Note:",
  "^See Also:",
  "^Examples:",
  '^";'
 );


$cur = 0;

while(<>)
{
  # %feature("docstring") something "text"
  # 
  next if /^\s*\/\/.*/;

  if( /^  / || /^\s*$/ )
  { 
    $valid = 0 if /NODOCSTRING/ ;
    $ignore = 1 if /CPPONLY/ ;

    die "out of range : $cur" if $cur < 0 || $cur > 9;

    $cont[ $cur] .= $_;
    next;
  }
  elsif( /^OriArgString/ )
  {
    ($ArgString) = $_ =~ /^OriArgString:\s*(.*)$/;
  }
  else
  {
    $cur = -1;
    for($i=0; $i<10; $i++)
    { 
      $cur = $i if /$Pattern[$i]/;
    }

    die "can not find any pattern?\n" . $_ unless $cur != -1;

    if( $cur == 0 )
    {
      for($i=0; $i<9; $i++)
      {  
        $empty =  $Pattern[$i] . '\s*$';
        $cont[$i] = "" if( $cont[$i] =~ /$empty/s );
      }

      next if $cont[0] eq "";

      $cont[$brief] =~ s/^  //mg;
      
      $text = join("\n", @cont);

      $alltext .=  $text if($valid == 1 && $ignore == 0);

      if( $ignore == 1 )
      {
        $alltext .= '%ignore ';
        ($func) = $cont[$head] =~ /.*docstring"\)\s*(.*?)\s*"/;
        $alltext .= $func . $ArgString . ";\n\n";
      }


      @cont = ("","","","","","","","","","" );
      $valid = 1;
      $ignore = 0;
      $ArgString = "";
    }
    # allow addition of new content
    # but we do not want repeative directives!
    $cont[$cur] .= $_ unless $cur == $brief || 
       ( $_ =~ $Pattern[$cur] && $cont[$cur] =~ $Pattern[$cur] );
  }
}


$alltext .= $text if( $valid == 1 && $ignore == 0);
$alltext =~ s/DEVONLY{.*?}//sg;
$alltext =~ s/\n\s*\n\s*\n/\n\n/sg;
$alltext =~ s/\\\s+"/\\"/sg;
print $alltext;
