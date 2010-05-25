#!/usr/bin/env perl

use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);

#Specific Files and Default Parameters
# VersionTag  1.00


#initialisation
srand (time());
%code = qw (TTT F TTC F TTA L TTG L TCT S TCC S TCA S TCG S TAT Y TAC Y TAA X TAG X TGT C TGC C TGA X TGG W CTT L CTC L CTA L CTG L CCT P CCC P CCA P CCG P CAT H CAC H CAA E CAG E CGT R CGC R CGA R CGG R ATT I ATC I ATA I ATG M ACT T ACC T ACA T ACG T AAT N AAC N AAA K AAG K AGT S AGC S AGA R AGG R GTT V GTC V GTA V GTG V GCT A GCC A GCA A GCG A GAT D GAC D GAA Q GAG Q GGT G GGC G GGA G GGG G);
@codons=keys %code;
@amino_acids=qw (A C D E F G H I K L M N P Q R S T V W X Y);
foreach $c (@codons)
  {
    $n=$number_codons{$code{$c}}++;
    $codon_list{$code{$c}}{$n}=$c;
    
  }
#Done


while (<>)
  {
    if ( /^>(\S+)/)
      {print $_;}
    else 
      {
	@sequence=/(\S)/g;
	
	foreach $amino_acid (@sequence)
	  {
	    $upper_case_aa= uc $amino_acid;

	    $n_codons=$number_codons{$upper_case_aa};
	    $chosen=int rand ($n_codons);
	    print "$codon_list{$upper_case_aa}{$chosen}",
	  }
	print "\n";
      }
  }

      
