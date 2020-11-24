#;-*- Perl -*-
#
#
##################  Alberto Roldan 9-2013  ###############

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;

#----------------------------------------------------------------------------------------------------
   $filename=$ARGV[0];
   $numargs=$#ARGV+1;
  if ($numargs lt 1 ) { if (-e "CONTCAR") { $filename="CONTCAR"; }; };

  print "\n ---  Loaded file is $filename  ---\n ";
#----------------------------------------------------------------------------------------------------
 system("cartesian.pl $filename");
 open IN, "$filename.cart";
   while (<IN>) { $file.= $_; };
 close (IN);
    @file=split(/\n/,$file);
#---------------------------------------------------------------------------------------------------- Z
    @line=split(/\s+/,@file[4]);
      foreach $l (@line) { if ($l) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
    @Z[0]=@line[0]; @Z[1]=@line[1]; @Z[2]=@line[2]+11.1671985392;
#---------------------------------------------------------------------------------------------------- Newpos
    $n=1;
      for ($i=9; $i<=56+8; $i++;) {
      foreach $l (@line) { if ($l) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
     @pos[$n]=(@line[0], @line[1], @line[2]+11.1671985392, "T T T");
      $n++;
     };
#---------------------------------------------------------------------------------------------------- print
 open OUT, ">>tempZ";
     print OUT "Fe  S Ni\n";
    for ($i=1; $i<=3; $i++) { print OUT "@file[$i]\n"; };
     print OUT "   @Z\n   ";
     print OUT " Fe  S Ni\n";
     print OUT " 47 64  1\n";
 close OUT;
 system("cat tempZ o-Td >> tempZ-1");
 open OUT, ">>tempZ-1";
     print OUT "@file[7]\n";
     print OUT "@file[8]\n";
    for ($i=9; $i<=9+7; $i++) { print OUT "@pos[$i]\n"; };
 close OUT;
 system("cat tempZ-1 o-Oh >> tempZ-2");
 open OUT, ">>tempZ-2";
    for ($i=9+8; $i<=9+8+16; $i++) { print OUT "@pos[$i]\n"; };
 close OUT;
 system("cat tempZ-2 o-S >> tempZ-3");
 open OUT, ">>tempZ-3";
    for ($i=9+8+17; $i<=9+8+34; $i++) { print OUT "@pos[$i]\n"; };
 close OUT;
#----------------------------------------------------------------------------------------------------
 system("mv tempZ-3 $filename.z2; rm tempZ tempZ-1 tempZ-2");













