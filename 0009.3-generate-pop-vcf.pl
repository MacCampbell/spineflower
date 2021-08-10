#! /usr/bin/perl -w

#08062021
#Generating some vcf files.
`mkdir vcf`;
#First making a poplist
#`cut -f 4,5 -d ',' ../meta/182.csv | perl -pe 's/Location,DateCollected\n//g' | sort | uniq > vcf/poplist.txt`;
#Adding in "Newhall" and "UpperLasVirgenesCanyon" to poplist for measurements of larger area
#Read in poplist

my $infile = "vcf/poplist.txt";

open (INFILE, "<", $infile);

while(<INFILE>) {
  chomp;
  my $pop=$_;
    print "Working on $pop \n";
  my $cleanpop = $pop;
  $cleanpop =~ s/,|\s|\"//g;
  #Lazily not making this general
  `grep $pop ../meta/182.csv | cut -f 2 -d ',' > vcf/$cleanpop.txt`;
  #Generating subsetted vcf with filtering
  `vcf-subset -c vcf/$cleanpop.txt recode.reheadered.pruned.vcf > vcf/$cleanpop.vcf`;
  #Filter for missing data and minmaf
  `bcftools view -i 'MAF > 0.1 && F_MISSING<0.05' vcf/$cleanpop.vcf -Ov -o vcf/$cleanpop.maf.vcf`;
}

close(INFILE)


