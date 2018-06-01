#!/usr/bin/perl
print "Uzorak,Duljina uzorka,Algoritam,Broj ponavljanja,Vrijeme,Udaljenost\n";
while (defined($line = <>))
{
    chomp $line;
    if( $line =~ /Query: (.*)/ ){
        $prefix = "$1";
    }
    if( $line =~ /Query length: ([^]*)/ ){
        $queryLength = "$1";
    }
    if( $line =~ /=>\s([^,]+),\s(\d+)[^=]+=([^s]*)s,[^=]+=(\d+)/ ){
        print "$prefix,$queryLength,$1,$2,$3,$4\n";
    }

}