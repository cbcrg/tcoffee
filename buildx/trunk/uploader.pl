#!/usr/bin/perl
use File::stat;

use strict;

my $source = $ARGV[0];
my $target = $ARGV[1];

my $SCPARGS=$ENV{'SCPARGS'};

if( -d $source && check($source) ) {
	
	my $scp = "scp -r $SCPARGS $source $target"; 
	print "$scp \n";
	my $r = `$scp`;
	`rm -rf $source`;
}


#
# Check that the file has not change in the elapsed sleep time
#
sub check() {
	my ($file) = (@_);

	my $s1 = fsize($file);
	# wait for 4 second and get again the size to check if it is changed
	sleep(4);
	my $s2 = fsize($file);
	
	my $result = ($s1 == $s2);
	
	if( !$result ) { print "Upload skipped because the folder content is changing .. \n"; }
	return $result;
}


sub fsize() {
	my ($path) = (@_);
	
	# check if it is a single file 
	if( -f $path ) {
		# just upload the file 
		my $size = -s $path;
		return $size;
	}
		
	# .. or a directory, in this case copy all the content 
	elsif( -d $path ) {
		opendir(IMD, "$path");
		my @thefiles= readdir(IMD);
		closedir(IMD);
			
		my $tot=0;
		for my $_file (@thefiles) {
			if( $_file ne "." && $_file ne "..") {
				$tot = $tot + fsize("$path/$_file");
			}
			
		}
		return $tot;
	}
}



