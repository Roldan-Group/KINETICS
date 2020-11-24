#;-*- Perl -*-
#
#
##################  Alberto Roldan 9-2013  ###############

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;


   $filename=$ARGV[0];
   $numargs=$#ARGV+1;
  if ($numargs lt 1 ) {
    print "\n ---  Provide a valid file  --- \n\n "; exit 0; 
   }elsif (-e "profile.dat") { $filename="profile.dat"; };
  print "\n ---  Loaded file is $filename  --- \n\n "; 
#----------------------------------------------------------------------------------------------------
 open IN, $filename;
   while (<IN>) { $file.= $_; };
 close (IN);
    @file=split(/\n/,$file);
      foreach $f (@file) { if ($f) { push(@Nfile,$f); }elsif ($f =~ 0) { push(@Nfile,"0.00"); };}; @file=@Nfile; @Nfile=();
       $plot=0; @data=(); $mMax=0;
    &graceHead_sub($filename);
    &gracePlot_sub($filename,$plot);
      $i=0;  $n=0.0; $m=0.25;
   while (@file[$i]) {
     @line=split(/\s+/,@file[$i]);
      foreach $l (@line){ if (($l) or ($l =~ 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
    if ((@line[0] =~ 0) or (@line[0] =~ /[0-9]$/)) {
       push(@data,"$n $line[0]");
       push(@data,"$m $line[0]");
       $n=$n+0.5; $m=$m+0.5;
      if ($m > $mMax) { $mMax=$m; };
     }else{ 
        if (@data) {
          open OUT, ">>data.tmp.agr";
            print OUT '@target G0.S'."$plot\n".'@type xy'."\n";
             foreach $d (@data) { @v=split(/\s+/,$d); printf OUT " %.3f %.3f\n",@v[0],@v[1]; };
            print OUT "&\n";
             $plot=$plot+1;
            print OUT '@target G0.S'."$plot\n".'@type xy'."\n";
             foreach $d (@data) { @v=split(/\s+/,$d); printf OUT " %.3f %.3f\n",@v[0],@v[1]; };
            print OUT "&\n";
           close OUT;
           $n=0; $m=0.25; @data=();
             $plot=$plot+1;
            &gracePlot_sub($filename,$plot);
     }; }; #if
    $i++;
   }; #while
        if (@data) {
          open OUT, ">>data.tmp.agr";
            print OUT '@target G0.S'."$plot\n".'@type xy'."\n";
             foreach $d (@data) { @v=split(/\s+/,$d); printf OUT " %.3f %.3f\n",@v[0],@v[1]; };
            print OUT "&\n";
             $plot=$plot+1;
            print OUT '@target G0.S'."$plot\n".'@type xy'."\n";
             foreach $d (@data) { @v=split(/\s+/,$d); printf OUT " %.3f %.3f\n",@v[0],@v[1]; };
            print OUT "&\n";
           close OUT;
             $n=0; $m=0.25; @data=();
        }; #if data
      &graceLine_sub($filename,$mMax);
#----------------------------------------------------------------------------------------------------
system ("cat $filename.tmp.agr data.tmp.agr >> $filename.agr; rm data.tmp.agr $filename.tmp.agr");
   sub graceLine_sub {
           ($name,$max)=@_; $max=$max-0.25;
 open OUT, ">>$name.tmp.agr";
   print OUT '@with line'."\n".'@    line on'."\n".'@    line loctype world'."\n".'@    line g0'."\n";
   print OUT '@    line -0.25, 0, '."$max".', 0'."\n".'@    line linewidth 0.5'."\n".'@    line linestyle 3'."\n".'@    line color 1'."\n";
   print OUT '@  line arrow 0'."\n".'@  line arrow type 0'."\n".'@  line arrow length 1.0'."\n";
   print OUT '@  line arrow layout 1.0, 1.0'."\n".'@line def'."\n";
 close OUT;
    return ();
 }; #--> sub graceLine
#==============================================================================================================================
   sub gracePlot_sub {
           ($name,$nameplot)=@_;
            $nameplot2=$nameplot+1;
 open OUT, ">>$name.tmp.agr";
   print OUT '@  s'."$nameplot".' type xy'."\n".'@  s'."$nameplot".' line type 4'."\n".'@    s0 line linestyle 1'."\n".'@    s0 line linewidth 1.5'."\n";
   print OUT '@  s'."$nameplot".' comment "setdata='."$nameplot".'"'."\n".'@    s'."$nameplot".' legend  "'."$nameplot".'"'."\n";
   print OUT '@  s'."$nameplot2".' type xy'."\n".'@  s'."$nameplot2".' line type 1'."\n".'@    s0 line linestyle 2'."\n".'@    s0 line linewidth 1.0'."\n";
   print OUT '@  s'."$nameplot2".' comment "setdata='."$nameplot".'"'."\n";
 close OUT;
    return ();
 }; #--> sub gracePlot
#==============================================================================================================================
   sub graceHead_sub {
           ($name)=@_;
 open OUT, ">>$name.tmp.agr";
   print OUT "\n# --- Alberto Roldan --- \n";
   print OUT "\n\n".'@page size 950, 600'."\n".'@page scroll 5%'."\n".'@page inout 5%'."\n";
   print OUT '@map font 29 to "StandardSymbolsL-Regular", "StandardSymbolsL-Regular"'."\n";
   print OUT '@map font 12 to "Symbol", "Symbol"'."\n".'@map font 31 to "Symbol-Regular", "Symbol-Regular"'."\n";
   print OUT '@default sformat "%.8g"'."\n";
# ------------------------------- AXIS
   print OUT '@g0 on'."\n".'@g0 type XY'."\n".'@with g0'."\n".'@    world 0, -2, 5, 2'."\n";
   print OUT '@    stack world 0, 0, 0, 0'."\n".'@    znorm 1'."\n".'@    view 0.168750, 0.150, 1.293750, 0.850'."\n";
   print OUT '@    xaxes scale Normal'."\n".'@    yaxes scale Normal'."\n".'@    xaxis  off'."\n";
   print OUT '@    yaxis  on'."\n".'@    yaxis  label "\xD\f{}E / eV"'."\n".'@    yaxis  label char size 1.750'."\n";
   print OUT '@    yaxis  tick on'."\n".'@    yaxis  tick major 0.2'."\n".'@    yaxis  tick minor ticks 1'."\n";
   print OUT '@    yaxis  tick in'."\n".'@    yaxis  tick major size 0.50'."\n".'@    yaxis  tick minor size 0.250'."\n";
 close OUT;
    return ();
 }; #--> sub graceHead






       
