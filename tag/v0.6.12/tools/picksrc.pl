#!/usr/bin/perl

# get objs
@objs = split(/[ \n\t\r]/, `find boost -type f -name '*' -print`);

foreach $file (@objs){
  next unless -f $file;

  # compile a file and see if it passes
  $cmd = "cd src;g++ -DHAVE_CONFIG_H -I. -I. -I..  " . 
    " -I/usr/site/python/include/python2.3 " . 
    " -I/usr/site/matlab/share/extern/include " . 
    "-I/home/bpeng/temp/simupop-0.5.2 -I..   " .
    " -E -MT simuPOP_wrap.o -MD -MP -MF '.deps/simuPOP_wrap.Tpo' " .
    " -o simuPOP_wrap.o simuPOP_wrap.cpp";
  print "dealing with $file \n";
#  print "/bin/mv -f " . $file . " ~/temp \n";
  die "can not move away file " if system("/bin/mv -f /home/bpeng/temp/simupop-0.5.2/" . $file . " /home/bpeng/temp") != 0;
  if( system($cmd) == 0 ){ # success
    print "$file removed \n";
  }else{ 
   ($filename) = $file =~ /^.*\/([^\/]*)$/;
#    print "/bin/mv -f /home/bpeng/temp/$filename  /home/bpeng/temp/$file\n";
    die "can not move file back " if system("/bin/mv -f /home/bpeng/temp/$filename  /home/bpeng/temp/simupop-0.5.2/$file") != 0;
 }
}



