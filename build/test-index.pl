#!/usr/bin/perl

# 
# Generate html index page for T-Coffee integration tests.
#
# It generate a "index.html" file listing the tests result. 
# The script have to used in the same folder containing the tests result
# and it will list all files having the suffix ".passed.html" or ".failed.html" 
# If a file named "result.failed" exists it means that at least a test has failed.
#
# Author: Paolo Di Tommaso
# Date: 16 June 2010
#

use strict;

open HTML, ">index.html" or die $!;

#
# write out the html header 
#
print HTML "<html>\n";
print HTML "<body>\n";
print HTML "<head>\n";
print HTML "<style>\@import url(test-index.css);</style> \n";
print HTML "</head>\n";

#
# check if exists the "result.failed" marked file
#
my $passed = ( -e "result.failed" ) ? 0 : 1;

#
# create the page hader 
#
my $headclass = $passed ? "passed" : "failed";
print HTML "<div id=\"header\" class=\"block $headclass\">\n";
print HTML "<h1>Test Results</h1>\n";           
print HTML "<p>Check the result of the latest T-Coffee integration build</p>\n";
print HTML "</div>\n";

#
# open the current directory 
#
opendir(IMD, ".");
my @thefiles= readdir(IMD);
closedir(IMD);


#
# FAILED tests
#
print HTML "<div id=\"tests\" class=\"block\"> \n"; 
print HTML "<h2><span>Failed tests</span></h2> \n";
print HTML "<ul>\n"; 
my $failedCount=0;
for my $_file (@thefiles) {
	next if !($_file =~ m/.+failed\.html/);

	print HTML "<li>\n"; 
	print HTML "<div class=\"test failed\" > \n"; 
	print HTML "<a href=\"$_file\">$_file</a> \n"; 
	print HTML "<div class=\"testResult\"></div> \n"; 
	print HTML "</div> \n"; 
	print HTML "</li> \n"; 
	
	$failedCount++;
}

if( $failedCount eq 0 ) { print HTML "<b>(none)</b>"; } 


print HTML "</ul></div>\n"; 

#
# PASSED tests
#
print HTML "<div id=\"tests\" class=\"block\"> \n"; 
print HTML "<h2><span>Passed tests</span></h2> \n";
print HTML "<ul>\n"; 
my $passedCount=0;

for my $_file (@thefiles) {
	next if !($_file =~ m/.+passed\.html/);

	print HTML "<li>\n"; 
	print HTML "<div class=\"test passed\" > \n"; 
	print HTML "<a href=\"$_file\">$_file</a> \n"; 
	print HTML "<div class=\"testResult\"></div> \n"; 
	print HTML "</div> \n"; 
	print HTML "</li> \n"; 
	
	$passedCount++;
}

if( $passedCount eq 0 ) { print HTML "<b>(none)</b>"; } 

print HTML "</ul></div>\n"; 


print HTML "</body>";

close HTML;
