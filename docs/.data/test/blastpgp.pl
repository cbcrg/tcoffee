#!/usr/bin/env perl

my $WSDL = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSBlastpgp.wsdl';

use SOAP::Lite;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;

my $checkInterval = 15;

my %params=(
	    'async' => '1', # Use async mode and simulate sync mode in client
	    );
GetOptions(
    "mode=s"           => \$params{mode}, # Search mode: PSI-Blast or PHI-Blast
    "database|d=s"     => \$params{database}, # Database to search
    "matrix|M=s"       => \$params{matrix},# Scoring maxtrix
    "exp|e=f"          => \$params{exp}, # E-value
    "expmulti|h=f"     => \$params{expmulti}, # E-value
    "filter|F=s"       => \$params{filter}, # Low complexity filter
    "dropoff|X=i"      => \$params{dropoff}, # Dropoff score
    "finaldropoff|Z=i" => \$params{finaldropoff}, # Final dropoff score
    "scores|v=i"       => \$params{scores}, # Max number of scores
    "align=i"          => \$params{align}, # Alignment view
    "startregion|S=i"  => \$params{startregion}, # Start of region in query
    "endregion|H=i"    => \$params{endregion}, # End of region in query
    "maxpasses|j=i"    => \$params{maxpasses}, # Number of PSI iterations
    "opengap|G=i"      => \$params{opengap}, # Gap open penalty
    "extendgap|E=i"    => \$params{extendgap}, # Gap extension penalty
    "pattern=s"        => \$params{pattern}, # PHI-BLAST pattern
    "usagemode|p=s"    => \$params{usagemode}, # PHI-BLAST program
    "appxml=s"         => \$params{appxml}, # Application XML
    "sequence=s"       => \$sequence, # Query sequence
    "help"	       => \$help, # Usage info
    "polljob"	       => \$polljob, # Get results
    "status"	       => \$status, # Get status
    "ids"      	       => \$ids, # Get ids from result
    "jobid=s"          => \$jobid, # JobId
    "outfile=s"        => \$outfile, # Output filename
    "outformat|o=s"    => \$outformat, # Output file format
    "async|a"	       => \$async, # Async submission
    "email=s"          => \$params{email}, # User e-mail address
    "trace"            => \$trace, # Show SOAP messages
    );

my $scriptName = basename($0, ());
if($help) {
    &usage();
    exit(0);
}

if ($trace){
    print "Tracing active\n";
    SOAP::Lite->import(+trace => 'debug');
}

my $soap = SOAP::Lite
    ->service($WSDL)
    ->on_fault(sub {
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
    !( (defined($ARGV[0]) && -f $ARGV[0]) || defined($sequence) )
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
    print STDOUT $result, "\n";
    if($result eq 'DONE') {
	print STDERR "To get results: $scriptName --polljob --jobid $jobid\n";
    }
}  
elsif($ids && defined($jobid)) {
    print STDERR "Getting ids from job $jobid\n";
    getIds($jobid);
}
else {
    if(-f $ARGV[0]) {	
	$content={type=>'sequence', content=>read_file($ARGV[0])};	
    }
    if($sequence) {	
	if(-f $sequence) {
	    $content={type=>'sequence', content=>read_file($sequence)};	
	} else {
	    $content={type=>'sequence', content=>$sequence};
	}
    }
    push @content, $content;

    my $jobid;
    my $paramsData = SOAP::Data->name('params')->type(map=>\%params);
    my $contentData = SOAP::Data->name('content')->value(\@content);
    # For SOAP::Lite 0.60 and earlier parameters are passed directly
    if($SOAP::Lite::VERSION eq '0.60' || $SOAP::Lite::VERSION =~ /0\.[1-5]/) {
        $jobid = $soap->runBlastpgp($paramsData, $contentData);
    }
    # For SOAP::Lite 0.69 and later parameter handling is different, so pass
    # undef's for templated params, and then pass the formatted args.
    else {
        $jobid = $soap->runBlastpgp(undef, undef,
				    $paramsData, $contentData);
    }

    if (defined($async)) {
	print STDOUT $jobid, "\n";
        print STDERR "To check status: $scriptName --status --jobid $jobid\n";
    } else { # Synchronous mode
        print STDERR "JobId: $jobid\n";
        sleep 1;
        getResults($jobid);
    }
}

sub getIds($) {
    $jobid = shift;
    my $results = $soap->getIds($jobid);
    for $result (@$results){
	print "$result\n";
    }
}

sub clientPoll($) {
    my $jobid = shift;
    my $result = 'PENDING';
    # Check status and wait if not finished
    #print STDERR "Checking status: $jobid\n";
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
    $jobid = shift;
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
        my $selResultType;
        foreach my $resultType (@$resultTypes) {
            if($resultType->{type} eq $outformat) {
                $selResultType = $resultType;
            }
        }
        $res=$soap->poll($jobid, $selResultType->{type});
        write_file($outfile.'.'.$selResultType->{ext}, $res);
    } else { # Data types available
        # Write a file for each output type
        for my $resultType (@$resultTypes){
            #print "Getting $resultType->{type}\n";
            $res=$soap->poll($jobid, $resultType->{type});
            write_file($outfile.'.'.$resultType->{ext}, $res);
        }
    }
}

sub read_file($) {
    my $filename = shift;
    open(FILE, $filename);
    my $content;
    my $buffer;
    while(sysread(FILE, $buffer, 1024)) {
	$content.= $buffer;
    }
    close(FILE);  
    return $content;
}

sub write_file($$) {
    my ($tmp,$entity) = @_;
    print STDERR "Creating result file: ".$tmp."\n";
    unless(open (FILE, ">$tmp")) {
	return 0;
    }
    syswrite(FILE, $entity);
    close (FILE);
    return 1;
}

sub usage {
    print STDERR <<EOF
Blastpgp
========
   
The blastpgp program implements the PSI-BLAST and PHI-BLAST variations
of NCBI BLAST.

For more detailed help information refer to
http://www.ebi.ac.uk/blastpgp/blastpsi_help_frame.html
 
Blastpgp specific options:

[Required]

      --mode            : str  : search mode to use: PSI-Blast or PHI-Blast
  -d, --database        : str  : protein database to search
  seqFile               : file : query sequence

[Optional]

  -M, --matrix          : str  : scoring matrix
  -e, --exp             : real : Expectation value
  -h, --expmulti        : real : threshold (multipass model)
  -F, --filter          : str  : filter query sequence with SEG [T,F]
  -m, --align           : int  : alignment view option:
                                 0 - pairwise, 1 - M/S identities,
                                 2 - M/S non-identities, 3 - Flat identities,
                                 4 - Flat non-identities
  -G, --opengap         : int  : cost to open a gap
  -E, --extendgap       : int  : cost to extend a gap
  -g, --gapalign        : str  : Gapped [T,F]
  -v, --scores          : int  : number of scores to be reported
  -j, --maxpasses       : int  : number of iterations
  -X, --dropoff         : int  : Dropoff score
  -Z, --finaldropoff    : int  : Dropoff for final alignment
  -S, --startregion     : int  : Start of required region in query
  -H, --endregion       : int  : End of required region in query
  -k, --pattern         : str  : Hit File (PHI-BLAST only)
  -p, --usagemode       : str  : Program option (PHI-BLAST only):
                                 blastpgp, patseedp, seedp

[General]

      --help            :      : prints this help text
  -a, --async           :      : forces to make an asynchronous query
      --status          :      : poll for the status of a job
      --polljob         :      : poll for the results of a job
      --jobid           : str  : jobid of an asynchronous job
      --ids             :      : get hit identifiers for result 
  -O, --outfile         : str  : name of the file results should be written to
                                 (default is based on the jobid)
  -o, --outformat       : str  : txt or xml output (no file is written)
      --trace           :      : show SOAP messages being interchanged

Synchronous job:

  The results/errors are returned as soon as the job is finished.
  Usage: blastpgp.pl --email <your@email> [options...] seqfile
  Returns: saves the results to disk

Asynchronous job:

  Use this if you want to retrieve the results at a later time. The results
  are stored for up to 24 hours.
  The asynchronous submission mode is recommended when users are submitting
  batch jobs or large database searches
  Usage: blastpgp.pl --email <your@email> --async [options...] seqFile
  Returns: jobid

  Use the jobid to query for the status of the job.
  Usage: blastpgp.pl --status --jobid <jobId>
  Returns: string indicating the status of the job
    DONE - job has finished
    RUNNING - job is running
    NOT_FOUND - job cannot be found
    ERROR - the jobs has encountered an error

  When done, use the jobid to retrieve the results of the job.
  Usage: blastpgp.pl --polljob --jobid <jobId> [--outfile <fileName>]
  Returns: saves the results to disk
EOF
;
}

