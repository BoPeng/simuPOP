#!/usr/bin/perl

undef $/;

$whole = <>;

while( $whole =~ /for\([^;=]*=[^;=]*;[^;\(\)]*\([^;\(\)]*\)\s*;/si )
{
($Pre, $init, $cond1, $op, $cond2, $After) = 
  $whole =~
 /^(.*)for\(([^;=]*=[^;=]*);\s*(\w+)\s*(==|<=|!=|<|>)([^;\(\)]*\([^;\(\)]*\))\s*;(.*)/si;
$whole = "$Pre" . "for($init, $cond1" . "End =$cond2; $cond1 $op $cond1" . "End; " . $After;
#print "$init  $cond1 $op $cond2 \n";
#print "$init, $cond1", "End =$cond2; $cond1 $op $cond1", "End\n";

}

print $whole;
