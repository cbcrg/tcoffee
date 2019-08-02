#!/usr/bin/env perl
my $WSDL = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSDaliLite.wsdl';

use SOAP::Lite;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;

my $checkInterval = 5;

my %params=(
	    'async' => '1', # Use async mode and simulate sync mode in client
	    );
GetOptions(
    'pdb1=s'     => \$params{'sequence1'},
    'chainid1=s' => \$params{'chainid1'},
    'pdb2=s'     => \$params{'sequence2'},
    'chainid2=s' => \$params{'chainid2'},
    "help|h"	 => \$help, # Usage info
    "async|a"	 => \$async, # Asynchronous submission
    "polljob"	 => \$polljob, # Get results
    "status"	 => \$status, # Get status
    "jobid|j=s"  => \$jobid, # JobId
    "email|S=s"  => \$params{email}, # E-mail address
    "trace"      => \$trace, # SOAP messages
    "sequence=s" => \$sequence, # Input PDB
    );

my $scriptName = basename($0, ());
if($help) {
    &usage();
    exit(0);
}

if($trace) {
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

if( !($polljob || $status) &&
    !( defined($params{'sequence1'}) && defined($params{'sequence2'}) )
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
    print STDOUT "$result", "\n";
    if($result eq 'DONE') {
	print STDERR "To get results: $scriptName --polljob --jobid $jobid\n";
    }
}
else {
    if(-f $params{'sequence1'}) {
	$params{'sequence1'} = read_file($params{'sequence1'});
    }
    if(-f $params{'sequence2'}) {
	$params{'sequence2'} = read_file($params{'sequence2'});
    }

    my $jobid;
    my $paramsData = SOAP::Data->name('params')->type(map=>\%params);
    # For SOAP::Lite 0.60 and earlier parameters are passed directly
    if($SOAP::Lite::VERSION eq '0.60' || $SOAP::Lite::VERSION =~ /0\.[1-5]/) {
        $jobid = $soap->runDaliLite($paramsData);
    }
    # For SOAP::Lite 0.69 and later parameter handling is different, so pass
    # undef's for templated params, and then pass the formatted args.
    else {
        $jobid = $soap->runDaliLite(undef,
				     $paramsData);
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
DaliLite
========

Pairwise comparison of protein structures

[Required]

  --pdb1                : str  : PDB ID for structure 1
  --pdb2                : str  : PDB ID for structure 2

[Optional]

  --chain1              : str  : Chain identifer in structure 1
  --chain2              : str  : Chain identifer in structure 2

[General]

  -h, --help            :      : prints this help text
  -S, --email           : str  : user email address
  -a, --async           :      : asynchronous submission
      --status          :      : poll for the status of a job
      --polljob         :      : poll for the results of a job
  -j, --jobid           : str  : jobid for an asynchronous job
  -O, --outfile         : str  : file name for results (default is jobid)
      --trace	        :      : show SOAP messages being interchanged 

Synchronous job:

  The results/errors are returned as soon as the job is finished.
  Usage: $scriptName --email <your\@email> [options] pdbFile [--outfile string]
  Returns: saves the results to disk

Asynchronous job:

  Use this if you want to retrieve the results at a later time. The results 
  are stored for up to 24 hours. 
  The asynchronous submission mode is recommended when users are submitting 
  batch jobs or large database searches	
  Usage: $scriptName --email <your\@email> --async [options] pdbFile
  Returns: jobid

  Use the jobid to query for the status of the job. 
  Usage: $scriptName --status --jobid <jobId>
  Returns: string indicating the status of the job:
    DONE - job has finished
    RUNNING - job is running
    NOT_FOUND - job cannot be found
    ERROR - the jobs has encountered an error

  When done, use the jobid to retrieve the status of the job. 
  Usage: $scriptName --polljob --jobid <jobId> [--outfile string]

[Help]

  For more detailed help information refer to
  http://www.ebi.ac.uk/DaliLite/
EOF
;
}

