#!/usr/bin/env perl
my $WSDL = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSWUBlast.wsdl';

use strict;
use SOAP::Lite;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;

my $checkInterval = 15;

my $numOpts = scalar(@ARGV);
my ($outfile, $outformat, $help, $async, $polljob, $status, $ids, $jobid, $trace, $sequence);
my %params= ( # Defaults
	      'async' => 1, # Force into async mode
	      'exp' => 10.0, # E-value threshold
	      'numal' => 50, # Maximum number of alignments
	      'scores' => 100, # Maximum number of scores
            );
GetOptions( # Map the options into variables
    "program|p=s"     => \$params{program}, # BLAST program
    "database|D=s"    => \$params{database}, # Search database
    "matrix|m=s"      => \$params{matrix}, # Scoring matrix
    "exp|E=f"         => \$params{exp}, # E-value threshold
    "echofilter|e"    => \$params{echofilter}, # Display filtered sequence
    "filter|f=s"      => \$params{filter}, # Low complexity filter name
    "alignments|b=i"  => \$params{numal}, # Number of alignments
    "scores|s=i"      => \$params{scores}, # Number of scores
    "sensitivity|S=s" => \$params{sensitivity}, # Search sensitivity
    "sort|t=s"	      => \$params{sort}, # Sort hits by...
    "stats|T=s"       => \$params{stats}, # Scoring statistic to use
    "strand|d=s"      => \$params{strand}, # Strand to use in DNA vs. DNA search
    "topcombon|c=i"   => \$params{topcombon}, # Consistent sets of HSPs
    "outfile=s"       => \$outfile, # Output file
    "outformat|o=s"   => \$outformat, # Output format
    "help|h"	      => \$help, # Usage info
    "async|a"	      => \$async, # Asynchronous mode
    "polljob"	      => \$polljob, # Get results
    "status"	      => \$status, # Get job status
    "ids"             => \$ids, # Get ids from result
    "jobid|j=s"       => \$jobid, # JobId
    "email=s"         => \$params{email}, # E-mail address
    "trace"           => \$trace, # SOAP trace
    "sequence=s"      => \$sequence, # Query sequence
    );

my $scriptName = basename($0, ());
if($help || $numOpts == 0) {
    &usage();
    exit(0);
}

if($trace){
    print STDERR "Tracing active\n";
    SOAP::Lite->import(+trace => 'debug');
}

my $soap = SOAP::Lite
    ->service($WSDL)
    ->proxy('http://localhost/',
    #proxy => ['http' => 'http://your.proxy.server/'], # HTTP proxy
    timeout => 600, # HTTP connection timeout
    )
    ->on_fault(sub { # SOAP fault handler
        my $soap = shift;
        my $res = shift;
        # Throw an exception for all faults
        if(ref($res) eq '') {
            die($res);
        } else {
            die($res->faultstring);
        }
        return new SOAP::SOM;
    }
               );

if( !($polljob || $status || $ids) &&
    !( defined($ARGV[0]) || defined($sequence) )
    ) {
    print STDERR 'Error: bad option combination', "\n";
    &usage();
    exit(1);
}
elsif($polljob && defined($jobid)) {
    print "Getting results for job $jobid\n";
    getResults($jobid);
}
elsif($status && defined($jobid)) {
    print STDERR "Getting status for job $jobid\n";
    my $result = $soap->checkStatus($jobid);
    print STDOUT "$result\n";
    if($result eq 'DONE') {
	print STDERR "To get results: $scriptName --polljob --jobid $jobid\n";
    }
}  
elsif($ids && defined($jobid)) {
    print STDERR "Getting ids from job $jobid\n";
    getIds($jobid);
}
else {
    # Prepare input data
    my $content;
    my (@contents) = ();
    if(-f $ARGV[0] || $ARGV[0] eq '-') {	
	$content={type=>'sequence',content=>read_file($ARGV[0])};	
    }
    if($sequence) {	
	if(-f $sequence || $sequence eq '-') {	
	    $content={type=>'sequence',content=>read_file($ARGV[0])};	
	} else {
	    $content={type=>'sequence',content=>$sequence};
	}
    }
    push @contents, $content;

    # Submit the job
    my $paramsData = SOAP::Data->name('params')->type(map=>\%params);
    my $contentData = SOAP::Data->name('content')->value(\@contents);
    # For SOAP::Lite 0.60 and earlier parameters are passed directly
    if($SOAP::Lite::VERSION eq '0.60' || $SOAP::Lite::VERSION =~ /0\.[1-5]/) {
        $jobid = $soap->runWUBlast($paramsData, $contentData);
    }
    # For SOAP::Lite 0.69 and later parameter handling is different, so pass
    # undef's for templated params, and then pass the formatted args.
    else {
        $jobid = $soap->runWUBlast(undef, undef,
				   $paramsData, $contentData);
    }

    # Asynchronous mode: output jobid and exit.
    if (defined($async)) {
	print STDOUT $jobid, "\n";
        print STDERR "To check status: $scriptName --status --jobid $jobid\n";
    }
    # Synchronous mode: try to get results
    else {
        print STDERR "JobId: $jobid\n";
        sleep 1;
        getResults($jobid);
    }
}

sub getIds($) {
    my $jobid = shift;
    my $results = $soap->getIds($jobid);
    for my $result (@$results){
	print "$result\n";
    }
}

sub clientPoll($) {
    my $jobid = shift;
    my $result = 'PENDING';
    # Check status and wait if not finished
    while($result eq 'RUNNING' || $result eq 'PENDING') {
        $result = $soap->checkStatus($jobid);
        print STDERR "$result\n";
        if($result eq 'RUNNING' || $result eq 'PENDING') {
            # Wait before polling again.
            sleep $checkInterval;
        }
    }
}

sub getResults($) {
    my $jobid = shift;
    my $res;
    # Check status, and wait if not finished
    clientPoll($jobid);
    # Use JobId if output file name is not defined
    unless(defined($outfile)) {
        $outfile=$jobid;
    }
    # Get list of data types
    my $resultTypes = $soap->getResults($jobid);
    # Get the data and write it to a file
    if(defined($outformat)) { # Specified data type
	if($outformat eq 'xml') {$outformat = 'toolxml';}
	if($outformat eq 'txt') {$outformat = 'tooloutput';}
        my $selResultType;
        foreach my $resultType (@$resultTypes) {
            if($resultType->{type} eq $outformat) {
                $selResultType = $resultType;
            }
        }
        $res=$soap->poll($jobid, $selResultType->{type});
	if($outfile eq '-') {
	     write_file($outfile, $res);
	} else {
	    write_file($outfile.'.'.$selResultType->{ext}, $res);
	}
    } else { # Data types available
        # Write a file for each output type
        for my $resultType (@$resultTypes){
            #print STDERR "Getting $resultType->{type}\n";
            $res=$soap->poll($jobid, $resultType->{type});
	    if($outfile eq '-') {
		write_file($outfile, $res);
	    } else {
		write_file($outfile.'.'.$resultType->{ext}, $res);
	    }
        }
    }
}

sub read_file($) {
    my $filename = shift;
    my ($content, $buffer);
    if($filename eq '-') {
	while(sysread(STDIN, $buffer, 1024)) {
	    $content .= $buffer;
	}
    }
    else { # File
	open(FILE, $filename) or die "Error: unable to open input file";
	while(sysread(FILE, $buffer, 1024)) {
	    $content .= $buffer;
	}
	close(FILE);
    }
    return $content;
}

sub write_file($$) {
    my ($filename, $data) = @_;
    print STDERR 'Creating result file: ' . $filename . "\n";
    if($filename eq '-') {
	print STDOUT $data;
    }
    else {
	open(FILE, ">$filename") or die "Error: unable to open output file";
	syswrite(FILE, $data);
	close(FILE);
    }
}

sub usage {
    print STDERR <<EOF
WU-BLAST
========

Rapid sequence database search programs utilizing the BLAST algorithm.
   
[Required]

      --email       : str  : user email address 
  -p, --program	    : str  : BLAST program to use: blastn, blastp, blastx, 
                             tblastn or tblastx
  -D, --database    : str  : database to search
  seqFile           : file : query sequence data file ("-" for STDIN)

[Optional]

  -m, --matrix	    : str  : scoring matrix
  -E, --exp	    : real : 0<E<= 1000. Statistical significance threshold
                             for reporting database sequence matches.
  -e, --echofilter  :      : display the filtered query sequence in the output
  -f, --filter	    : str  : activates filtering of the query sequence
  -b, --alignments  : int  : number of alignments to be reported
  -s, --scores	    : int  : number of scores to be reported
  -S, --sensitivity : str  :
  -t, --sort	    : str  :
  -T, --stats       : str  :
  -d, --strand      : str  : DNA strand to search with in DNA vs. DNA searches 
  -c, --topcombon   :      :

[General]	

  -h, --help       :      : prints this help text
  -a, --async      :      : forces to make an asynchronous query
      --status     :      : poll for the status of a job
      --polljob    :      : poll for the results of a job
  -j, --jobid      : str  : jobid that was returned when an asynchronous job 
                            was submitted.
  -O, --outfile    : str  : name of the file results should be written to 
                            (default is based on the jobid; "-" for STDOUT)
  -o, --outformat  : str  : txt or xml output (no file is written)
      --trace	   :      : show SOAP messages being interchanged 

Synchronous job:

  The results/errors are returned as soon as the job is finished.
  Usage: $scriptName --email <your\@email> [options...] seqFile
  Returns: saves the results to disk

Asynchronous job:

  Use this if you want to retrieve the results at a later time. The results 
  are stored for up to 24 hours. 
  The asynchronous submission mode is recommended when users are submitting 
  batch jobs or large database searches	
  Usage: $scriptName --async --email <your\@email> [options...] seqFile
  Returns : jobid

  Use the jobid to query for the status of the job. 
  Usage: $scriptName --status --jobid <jobId>
  Returns : string indicating the status of the job:
    DONE - job has finished
    RUNNING - job is running
    NOT_FOUND - job cannot be found
    ERROR - the jobs has encountered an error

  When done, use the jobid to retrieve the status of the job. 
  Usage: $scriptName --polljob --jobid <jobId> [--outfile string]
  Returns: saves the results to disk

[Help]

For more detailed help information refer to 
http://www.ebi.ac.uk/blast2/WU-Blast2_Help_frame.html
 
EOF
;
}

