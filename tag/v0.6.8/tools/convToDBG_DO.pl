#!/usr/bin/perl

undef $/;

$whole = <>;

while( $whole =~ /#ifdef\s+DEBUG\s+if\(\s*debug[^;]*;\s*#endif/si )
{
($Pre, $code, $Expr, $After) = 
  $whole =~
 /^(.*)#ifdef\s+DEBUG\s+if\(\s*debug\(\s*(\w+)\s*\)\s*\)\s*([^;]*?);\s*#endif(.*)/si;
$whole = "$Pre\nDBG_DO($code, $Expr);\n$After";
}
while( $whole =~ /#ifdef\s+DEBUG\s+if\([^;]*?\)\s*throw[^;]*;\s*#endif/si )
{
($Pre, $cond, $excep, $Expr, $After) = 
  $whole =~ /^(.*)#ifdef\s+DEBUG\s+if\(([^;]*?)\)\s*throw\s+(\w+)\(([^;]*?)\);\s*#endif(.*)/si;
$whole = "$Pre\nDBG_FAILIF($cond, $excep, $Expr);\n$After";
}
print $whole;
