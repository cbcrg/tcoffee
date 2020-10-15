#!/usr/bin/env perl

=head1 NAME

ncbiblast.pl

=head1 DESCRIPTION

NCBI Blast (REST) web service Perl client using L<LWP>.

Tested with:

=over

=item *
L<LWP> 6.35, L<XML::Simple> 2.25 and Perl 5.22.0 (MacOS 10.13.6)

=back

For further information see:

=over

=item *
L<https://www.ebi.ac.uk/Tools/webservices/>

=back

=head1 LICENSE

Copyright 2012-2018 EMBL - European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Perl Client Automatically generated with:
https://github.com/ebi-wp/webservice-clients-generator

=cut

# ======================================================================
# Enable Perl warnings
use strict;
use warnings;

# Load libraries
use English;
use LWP;
use XML::Simple;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use Data::Dumper;
use Time::HiRes qw(usleep);

# Base URL for service
my $baseUrl = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast';
my $version = '2020-02-28 14:43';

# Set interval for checking status
my $checkInterval = 3;

# Set maximum number of 'ERROR' status calls to call job failed.
my $maxErrorStatusCount = 3;

# Output level
my $outputLevel = 1;

# Process command-line options
my $numOpts = scalar(@ARGV);
my %params = (
    'debugLevel' => 0,
    'maxJobs'    => 1
);
my @database;
# Default parameter values (should get these from the service)
GetOptions(
    # Tool specific options
    'program=s'       => \$params{'program'},        # The BLAST program to be used for the Sequence Similarity Search.
    'task=s'          => \$params{'task'},           # Task option (only selectable for blastn)
    'matrix=s'        => \$params{'matrix'},         # (Protein searches) The substitution matrix used for scoring alignments when searching the database.
    'alignments=i'    => \$params{'alignments'},     # Maximum number of match alignments reported in the result output.
    'scores=i'        => \$params{'scores'},         # Maximum number of match score summaries reported in the result output.
    'exp=s'           => \$params{'exp'},            # Limits the number of scores and alignments reported based on the expectation value. This is the maximum number of times the match is expected to occur by chance.
    'dropoff=i'       => \$params{'dropoff'},        # The amount a score can drop before gapped extension of word hits is halted
    'match_scores=s'  => \$params{'match_scores'},   # (Nucleotide searches) The match score is the bonus to the alignment score when matching the same base. The mismatch is the penalty when failing to match.
    'gapopen=i'       => \$params{'gapopen'},        # Penalty taken away from the score when a gap is created in sequence. Increasing the gap openning penalty will decrease the number of gaps in the final alignment.
    'gapext=i'        => \$params{'gapext'},         # Penalty taken away from the score for each base or residue in the gap. Increasing the gap extension penalty favors short gaps in the final alignment, conversly decreasing the gap extension penalty favors long gaps in the final alignment.
    'filter=s'        => \$params{'filter'},         # Filter regions of low sequence complexity. This can avoid issues with low complexity sequences where matches are found due to composition rather than meaningful sequence similarity. However in some cases filtering also masks regions of interest and so should be used with caution.
    'seqrange=s'      => \$params{'seqrange'},       # Specify a range or section of the input sequence to use in the search. Example: Specifying '34-89' in an input sequence of total length 100, will tell BLAST to only use residues 34 to 89, inclusive.
    'gapalign'        => \$params{'gapalign'},       # This is a true/false setting that tells the program the perform optimised alignments within regions involving gaps. If set to true, the program will perform an alignment using gaps. Otherwise, if it is set to false, it will report only individual HSP where two sequence match each other, and thus will not produce alignments with gaps.
    'wordsize=i'      => \$params{'wordsize'},       # Word size for wordfinder algorithm
    'taxids=s'        => \$params{'taxids'},         # Specify one or more TaxIDs so that the BLAST search becomes taxonomically aware.
    'compstats=s'     => \$params{'compstats'},      # Use composition-based statistics.
    'align=i'         => \$params{'align'},          # Formating for the alignments
    'transltable=i'   => \$params{'transltable'},    # Query Genetic code to use in translation
    'stype=s'         => \$params{'stype'},          # Indicates if the sequence is protein or DNA/RNA.
    'sequence=s'      => \$params{'sequence'},       # The query sequence can be entered directly into this form. The sequence can be in GCG, FASTA, EMBL (Nucleotide only), GenBank, PIR, NBRF, PHYLIP or UniProtKB/Swiss-Prot (Protein only) format. A partially formatted sequence is not accepted. Adding a return to the end of the sequence may help certain applications understand the input. Note that directly using data from word processors may yield unpredictable results as hidden/control characters may be present.
    'database=s'      => \@database,                 # Database
    # Generic options
    'email=s'         => \$params{'email'},          # User e-mail address
    'title=s'         => \$params{'title'},          # Job title
    'outfile=s'       => \$params{'outfile'},        # Output file name
    'outformat=s'     => \$params{'outformat'},      # Output file type
    'jobid=s'         => \$params{'jobid'},          # JobId
    'help|h'          => \$params{'help'},           # Usage help
    'asyncjob'        => \$params{'asyncjob'},       # Asynchronous submission
    'polljob'         => \$params{'polljob'},        # Get results
    'pollFreq=f'      => \$params{'pollFreq'},       # Poll Frequency
    'resultTypes'     => \$params{'resultTypes'},    # Get result types
    'status'          => \$params{'status'},         # Get status
    'params'          => \$params{'params'},         # List input parameters
    'paramDetail=s'   => \$params{'paramDetail'},    # Get details for parameter
    'multifasta'      => \$params{'multifasta'},     # Multiple fasta input
    'useSeqId'        => \$params{'useSeqId'},       # Seq Id file name
    'maxJobs=i'       => \$params{'maxJobs'},        # Max. parallel jobs

    'verbose'         => \$params{'verbose'},        # Increase output level
    'version'         => \$params{'version'},        # Prints out the version of the Client and exit.
    'quiet'           => \$params{'quiet'},          # Decrease output level
    'debugLevel=i'    => \$params{'debugLevel'},     # Debugging level
    'baseUrl=s'       => \$baseUrl,                  # Base URL for service.
);
if ($params{'verbose'}) {$outputLevel++}
if ($params{'quiet'}) {$outputLevel--}
if ($params{'pollFreq'}) {$checkInterval = $params{'pollFreq'} * 1000 * 1000}
if ($params{'baseUrl'}) {$baseUrl = $params{'baseUrl'}}
$params{"database"} = [@database];

# Debug mode: LWP version
&print_debug_message('MAIN', 'LWP::VERSION: ' . $LWP::VERSION,
    1);

# Debug mode: print the input parameters
&print_debug_message('MAIN', "params:\n" . Dumper(\%params), 11);

# LWP UserAgent for making HTTP calls (initialised when required).
my $ua;

# Get the script filename for use in usage messages
my $scriptName = basename($0, ());

# Print usage and exit if requested
if ($params{'help'} || $numOpts == 0) {
    &usage();
    exit(0);
}

# Debug mode: show the base URL
&print_debug_message('MAIN', 'baseUrl: ' . $baseUrl, 1);
if (
    !(
        $params{'polljob'}
            || $params{'resultTypes'}
            || $params{'status'}
            || $params{'params'}
            || $params{'paramDetail'}
            || $params{'version'}
    )
        && !(defined($ARGV[0]) || defined($params{'sequence'}))
) {

    # Bad argument combination, so print error message and usage
    print STDERR 'Error: bad option combination', "\n";
    &usage();
    exit(1);
}
# Get parameters list
elsif ($params{'params'}) {
    &print_tool_params();
}

# Get parameter details
elsif ($params{'paramDetail'}) {
    &print_param_details($params{'paramDetail'});
}

# Print Client version
elsif ($params{'version'}) {
  print STDOUT 'Revision: ' . $version, "\n";
  exit(1);
}

# Job status
elsif ($params{'status'} && defined($params{'jobid'})) {
    &print_job_status($params{'jobid'});
}

# Result types
elsif ($params{'resultTypes'} && defined($params{'jobid'})) {
    &print_result_types($params{'jobid'});
}

# Poll job and get results
elsif ($params{'polljob'} && defined($params{'jobid'})) {
    &get_results($params{'jobid'});
}

# Submit a job
else {
    # Multiple input sequence mode, assume fasta format.
    if (defined($params{'multifasta'}) && $params{'multifasta'}) {
        &multi_submit_job();
    }

    # Entry identifier list file.
    elsif ((defined($params{'sequence'}) && $params{'sequence'} =~ m/^\@/)
        || (defined($ARGV[0]) && $ARGV[0] =~ m/^\@/)) {
        my $list_filename = $params{'sequence'} || $ARGV[0];
        $list_filename =~ s/^\@//;
        &list_file_submit_job($list_filename);
    }
    # Default: single sequence/identifier.
    else {
        # Warn for invalid batch only option use.
        if (defined($params{'useSeqId'}) && $params{'useSeqId'}) {
            print STDERR "Warning: --useSeqId option ignored.\n";
            delete $params{'useSeqId'};
        }
        if (defined($params{'maxJobs'}) && $params{'maxJobs'} > 1) {
            print STDERR "Warning: --maxJobs option ignored.\n";
            $params{'maxJobs'} = 1;
        }
        # Load the sequence data and submit.
        &submit_job(&load_data());
    }
}



=head1 FUNCTIONS

=cut

### Wrappers for REST resources ###

=head2 rest_user_agent()

Get a LWP UserAgent to use to perform REST requests.

  my $ua = &rest_user_agent();

=cut

sub rest_user_agent() {
    print_debug_message('rest_user_agent', 'Begin', 21);
    # Create an LWP UserAgent for making HTTP calls.
    my $ua = LWP::UserAgent->new();
    # Set 'User-Agent' HTTP header to identifiy the client.
    my $revisionNumber = 0;
    $revisionNumber = "Revision: " . $version;
    $ua->agent("EBI-Sample-Client/$revisionNumber ($scriptName; $OSNAME) " . $ua->agent());
    # Configure HTTP proxy support from environment.
    $ua->env_proxy;
    print_debug_message('rest_user_agent', 'End', 21);
    return $ua;
}

=head2 rest_error()

Check a REST response for an error condition. An error is mapped to a die.

  &rest_error($response, $content_data);

=cut

sub rest_error() {
    print_debug_message('rest_error', 'Begin', 21);
    my $response = shift;
    my $contentdata;
    if (scalar(@_) > 0) {
        $contentdata = shift;
    }
    if (!defined($contentdata) || $contentdata eq '') {
        $contentdata = $response->content();
    }
    # Check for HTTP error codes
    if ($response->is_error) {
        my $error_message = '';
        # HTML response.
        if ($contentdata =~ m/<h1>([^<]+)<\/h1>/) {
            $error_message = $1;
        }
        #  XML response.
        elsif ($contentdata =~ m/<description>([^<]+)<\/description>/) {
            $error_message = $1;
        }
        # die 'http status: ' . $response->code . ' ' . $response->message . '  ' . $error_message;
    }
    print_debug_message('rest_error', 'End', 21);
}

=head2 rest_request()

Perform a REST request (HTTP GET).

  my $response_str = &rest_request($url);

=cut

sub rest_request {
    print_debug_message('rest_request', 'Begin', 11);
    my $requestUrl = shift;
    print_debug_message('rest_request', 'URL: ' . $requestUrl, 11);

    # Get an LWP UserAgent.
    $ua = &rest_user_agent() unless defined($ua);
    # Available HTTP compression methods.
    my $can_accept;
    eval {
        $can_accept = HTTP::Message::decodable();
    };
    $can_accept = '' unless defined($can_accept);
    # Perform the request
    my $response = $ua->get($requestUrl,
        'Accept-Encoding' => $can_accept, # HTTP compression.
    );
    print_debug_message('rest_request', 'HTTP status: ' . $response->code,
        11);
    print_debug_message('rest_request',
        'response length: ' . length($response->content()), 11);
    print_debug_message('rest_request',
        'request:' . "\n" . $response->request()->as_string(), 32);
    print_debug_message('rest_request',
        'response: ' . "\n" . $response->as_string(), 32);
    # Unpack possibly compressed response.
    my $retVal;
    if (defined($can_accept) && $can_accept ne '') {
        $retVal = $response->decoded_content();
    }
    # If unable to decode use orginal content.
    $retVal = $response->content() unless defined($retVal);
    # Check for an error.
    &rest_error($response, $retVal);
    print_debug_message('rest_request', 'retVal: ' . $retVal, 12);
    print_debug_message('rest_request', 'End', 11);

    # Return the response data
    return $retVal;
}

=head2 rest_get_parameters()

Get list of tool parameter names.

  my (@param_list) = &rest_get_parameters();

=cut

sub rest_get_parameters {
    print_debug_message('rest_get_parameters', 'Begin', 1);
    my $url = $baseUrl . '/parameters/';
    my $param_list_xml_str = rest_request($url);
    my $param_list_xml = XMLin($param_list_xml_str);
    my (@param_list) = @{$param_list_xml->{'id'}};
    print_debug_message('rest_get_parameters', 'End', 1);
    return(@param_list);
}

=head2 rest_get_parameter_details()

Get details of a tool parameter.

  my $paramDetail = &rest_get_parameter_details($param_name);

=cut

sub rest_get_parameter_details {
    print_debug_message('rest_get_parameter_details', 'Begin', 1);
    my $parameterId = shift;
    print_debug_message('rest_get_parameter_details',
        'parameterId: ' . $parameterId, 1);
    my $url = $baseUrl . '/parameterdetails/' . $parameterId;
    my $param_detail_xml_str = rest_request($url);
    my $param_detail_xml = XMLin($param_detail_xml_str);
    print_debug_message('rest_get_parameter_details', 'End', 1);
    return($param_detail_xml);
}

=head2 rest_run()

Submit a job.

  my $job_id = &rest_run($email, $title, \%params );

=cut

sub rest_run {
    print_debug_message('rest_run', 'Begin', 1);
    my $email = shift;
    my $title = shift;
    my $params = shift;
    $email = '' if (!$email);
    print_debug_message('rest_run', 'email: ' . $email, 1);
    if (defined($title)) {
        print_debug_message('rest_run', 'title: ' . $title, 1);
    }
    print_debug_message('rest_run', 'params: ' . Dumper($params), 1);

    # Get an LWP UserAgent.
    $ua = &rest_user_agent() unless defined($ua);

    # Clean up parameters
    my (%tmp_params) = %{$params};
    $tmp_params{'email'} = $email;
    $tmp_params{'title'} = $title;
    foreach my $param_name (keys(%tmp_params)) {
        if (!defined($tmp_params{$param_name})) {
            delete $tmp_params{$param_name};
        }
    }

    # Submit the job as a POST
    my $url = $baseUrl . '/run';
    my $response = $ua->post($url, \%tmp_params);
    print_debug_message('rest_run', 'HTTP status: ' . $response->code, 11);
    print_debug_message('rest_run',
        'request:' . "\n" . $response->request()->as_string(), 11);
    print_debug_message('rest_run',
        'response: ' . length($response->as_string()) . "\n" . $response->as_string(), 11);

    # Check for an error.
    &rest_error($response);

    # The job id is returned
    my $job_id = $response->content();
    print_debug_message('rest_run', 'End', 1);
    return $job_id;
}

=head2 rest_get_status()

Check the status of a job.

  my $status = &rest_get_status($job_id);

=cut

sub rest_get_status {
    print_debug_message('rest_get_status', 'Begin', 1);
    my $job_id = shift;
    print_debug_message('rest_get_status', 'jobid: ' . $job_id, 2);
    my $status_str = 'UNKNOWN';
    my $url = $baseUrl . '/status/' . $job_id;
    $status_str = &rest_request($url);
    print_debug_message('rest_get_status', 'status_str: ' . $status_str, 2);
    print_debug_message('rest_get_status', 'End', 1);
    return $status_str;
}

=head2 rest_get_result_types()

Get list of result types for finished job.

  my (@result_types) = &rest_get_result_types($job_id);

=cut

sub rest_get_result_types {
    print_debug_message('rest_get_result_types', 'Begin', 1);
    my $job_id = shift;
    print_debug_message('rest_get_result_types', 'jobid: ' . $job_id, 2);
    my (@resultTypes);
    my $url = $baseUrl . '/resulttypes/' . $job_id;
    my $result_type_list_xml_str = &rest_request($url);
    my $result_type_list_xml = XMLin($result_type_list_xml_str);
    (@resultTypes) = @{$result_type_list_xml->{'type'}};
    print_debug_message('rest_get_result_types',
        scalar(@resultTypes) . ' result types', 2);
    print_debug_message('rest_get_result_types', 'End', 1);
    return(@resultTypes);
}

=head2 rest_get_result()

Get result data of a specified type for a finished job.

  my $result = rest_get_result($job_id, $result_type);

=cut

sub rest_get_result {
    print_debug_message('rest_get_result', 'Begin', 1);
    my $job_id = shift;
    my $type = shift;
    print_debug_message('rest_get_result', 'jobid: ' . $job_id, 1);
    print_debug_message('rest_get_result', 'type: ' . $type, 1);
    my $url = $baseUrl . '/result/' . $job_id . '/' . $type;
    my $result = &rest_request($url);
    print_debug_message('rest_get_result', length($result) . ' characters',
        1);
    print_debug_message('rest_get_result', 'End', 1);
    return $result;
}

### Service actions and utility functions ###

=head2 print_debug_message()

Print debug message at specified debug level.

  &print_debug_message($method_name, $message, $level);

=cut

sub print_debug_message {
    my $function_name = shift;
    my $message = shift;
    my $level = shift;
    if ($level <= $params{'debugLevel'}) {
        print STDERR '[', $function_name, '()] ', $message, "\n";
    }
}

=head2 print_tool_params()

Print list of tool parameters.

  &print_tool_params();

=cut

sub print_tool_params {
    print_debug_message('print_tool_params', 'Begin', 1);
    my (@param_list) = &rest_get_parameters();
    foreach my $param (sort (@param_list)) {
        print $param, "\n";
    }
    print_debug_message('print_tool_params', 'End', 1);
}

=head2 print_param_details()

Print details of a tool parameter.

  &print_param_details($param_name);

=cut

sub print_param_details {
    print_debug_message('print_param_details', 'Begin', 1);
    my $paramName = shift;
    print_debug_message('print_param_details', 'paramName: ' . $paramName, 2);
    my $paramDetail = &rest_get_parameter_details($paramName);
    print $paramDetail->{'name'}, "\t", $paramDetail->{'type'}, "\n";
    print $paramDetail->{'description'}, "\n";
    if (defined($paramDetail->{'values'}->{'value'})) {
        if (ref($paramDetail->{'values'}->{'value'}) eq 'ARRAY') {
            foreach my $value (@{$paramDetail->{'values'}->{'value'}}) {
                &print_param_value($value);
            }
        }
        else {
            &print_param_value($paramDetail->{'values'}->{'value'});
        }
    }
    print_debug_message('print_param_details', 'End', 1);
}

=head2 print_param_value()

Print details of a tool parameter value.

  &print_param_details($param_value);

Used by print_param_details() to handle both singluar and array values.

=cut

sub print_param_value {
    my $value = shift;
    print $value->{'value'};
    if ($value->{'defaultValue'} eq 'true') {
        print "\t", 'default';
    }
    print "\n";
    print "\t", $value->{'label'}, "\n";
    if (defined($value->{'properties'})) {
        foreach
        my $key (sort ( keys(%{$value->{'properties'}{'property'}}) )) {
            if (ref($value->{'properties'}{'property'}{$key}) eq 'HASH'
                && defined($value->{'properties'}{'property'}{$key}{'value'})
            ) {
                print "\t", $key, "\t",
                    $value->{'properties'}{'property'}{$key}{'value'}, "\n";
            }
            else {
                print "\t", $value->{'properties'}{'property'}{'key'},
                    "\t", $value->{'properties'}{'property'}{'value'}, "\n";
                last;
            }
        }
    }
}

=head2 print_job_status()

Print status of a job.

  &print_job_status($job_id);

=cut

sub print_job_status {
    print_debug_message('print_job_status', 'Begin', 1);
    my $jobid = shift;
    print_debug_message('print_job_status', 'jobid: ' . $jobid, 1);
    if ($outputLevel > 0) {
        print STDERR 'Getting status for job ', $jobid, "\n";
    }
    my $result = &rest_get_status($jobid);
    print "$result\n";
    if ($result eq 'FINISHED' && $outputLevel > 0) {
        print STDERR "To get results: perl $scriptName --polljob --jobid " . $jobid
            . "\n";
    }
    print_debug_message('print_job_status', 'End', 1);
}

=head2 print_result_types()

Print available result types for a job.

  &print_result_types($job_id);

=cut

sub print_result_types {
    print_debug_message('result_types', 'Begin', 1);
    my $jobid = shift;
    print_debug_message('result_types', 'jobid: ' . $jobid, 1);
    if ($outputLevel > 0) {
        print STDERR 'Getting result types for job ', $jobid, "\n";
    }
    my $status = &rest_get_status($jobid);
    if ($status eq 'PENDING' || $status eq 'RUNNING') {
        print STDERR 'Error: Job status is ', $status,
            '. To get result types the job must be finished.', "\n";
    }
    else {
        my (@resultTypes) = &rest_get_result_types($jobid);
        if ($outputLevel > 0) {
            print STDOUT 'Available result types:', "\n";
        }
        foreach my $resultType (@resultTypes) {
            print STDOUT $resultType->{'identifier'}, "\n";
            if (defined($resultType->{'label'})) {
                print STDOUT "\t", $resultType->{'label'}, "\n";
            }
            if (defined($resultType->{'description'})) {
                print STDOUT "\t", $resultType->{'description'}, "\n";
            }
            if (defined($resultType->{'mediaType'})) {
                print STDOUT "\t", $resultType->{'mediaType'}, "\n";
            }
            if (defined($resultType->{'fileSuffix'})) {
                print STDOUT "\t", $resultType->{'fileSuffix'}, "\n";
            }
        }
        if ($status eq 'FINISHED' && $outputLevel > 0) {
            print STDERR "\n", 'To get results:', "\n",
                "  perl $scriptName --polljob --jobid " . $params{'jobid'} . "\n",
                "  perl $scriptName --polljob --outformat <type> --jobid "
                    . $params{'jobid'} . "\n";
        }
    }
    print_debug_message('result_types', 'End', 1);
}

=head2 submit_job()

Submit a job to the service.

  &submit_job($seq);

=cut

sub submit_job {
    print_debug_message('submit_job', 'Begin', 1);

    # Set input sequence
    $params{'sequence'} = shift;
    my $seq_id = shift;

    # Load parameters
    &load_params();

    # Submit the job
    my $jobid = &rest_run($params{'email'}, $params{'title'}, \%params);

    # Asynchronous submission.
    if (defined($params{'asyncjob'})) {
        print STDOUT $jobid, "\n";
        if ($outputLevel > 0) {
            print STDERR
                "To check status: perl $scriptName --status --jobid $jobid\n";
        }
    }

    # Simulate synchronous submission serial mode.
    else {
        if ($outputLevel > 0) {
            print STDERR "JobId: $jobid\n";
        } else {
            print STDERR "$jobid\n";
        }
        usleep($checkInterval);
        # Get results.
        &get_results($jobid, $seq_id);

    }
    print_debug_message('submit_job', 'End', 1);
    return $jobid;
}
=head2 multi_submit_job()

Submit multiple jobs assuming input is a collection of fasta formatted sequences.

  &multi_submit_job();

=cut

sub multi_submit_job {
    print_debug_message('multi_submit_job', 'Begin', 1);
    my (@filename_list) = ();

    # Query sequence
    if (defined($ARGV[0])) {                  # Bare option
        if (-f $ARGV[0] || $ARGV[0] eq '-') { # File
            push(@filename_list, $ARGV[0]);
        }
        else {
            warn 'Warning: Input file "' . $ARGV[0] . '" does not exist';
        }
    }
    if ($params{'sequence'}) {                                      # Via --sequence
        if (-f $params{'sequence'} || $params{'sequence'} eq '-') { # File
            push(@filename_list, $params{'sequence'});
        }
        else {
            warn 'Warning: Input file "'
                . $params{'sequence'}
                . '" does not exist';
        }
    }

    # Job identifier tracking for parallel execution.
    my @jobid_list = ();
    my $job_number = 0;
    $/ = '>';
    foreach my $filename (@filename_list) {
        my $INFILE;
        if ($filename eq '-') { # STDIN.
            open($INFILE, '<-')
                or die 'Error: unable to STDIN (' . $! . ')';
        }
        else { # File.
            open($INFILE, '<', $filename)
                or die 'Error: unable to open file '
                . $filename . ' ('
                . $! . ')';
        }
        while (<$INFILE>) {
            my $seq = $_;
            $seq =~ s/>$//;
            if ($seq =~ m/(\S+)/) {
                my $seq_id = $1;
                print STDERR "Submitting job for: $seq_id\n"
                    if ($outputLevel > 0);
                $seq = '>' . $seq;
                &print_debug_message('multi_submit_job', $seq, 11);
                $job_number++;
                my $job_id = &submit_job($seq, $seq_id);

                my $job_info_str = sprintf('%s %d %d', $job_id, 0, $job_number);

                push(@jobid_list, $job_info_str);
            }

            # Parallel mode, wait for job(s) to finish to free slots.
            while ($params{'maxJobs'} > 1
                && scalar(@jobid_list) >= $params{'maxJobs'}) {
                &_job_list_poll(\@jobid_list);
                print_debug_message('multi_submit_job',
                    'Remaining jobs: ' . scalar(@jobid_list), 1);
            }
        }
        close $INFILE;
    }

    # Parallel mode, wait for remaining jobs to finish.
    while ($params{'maxJobs'} > 1 && scalar(@jobid_list) > 0) {
        &_job_list_poll(\@jobid_list);
        print_debug_message('multi_submit_job',
            'Remaining jobs: ' . scalar(@jobid_list), 1);
    }
    print_debug_message('multi_submit_job', 'End', 1);
}


=head2 _job_list_poll()

Poll the status of a list of jobs and fetch results for finished jobs.

  while(scalar(@jobid_list) > 0) {
    &_job_list_poll(\@jobid_list);
  }

=cut

sub _job_list_poll {
    print_debug_message('_job_list_poll', 'Begin', 1);
    my $jobid_list = shift;
    print_debug_message('_job_list_poll', 'Num jobs: ' . scalar(@$jobid_list),
        11);

    # Loop though job Id list polling job status.
    for (my $jobNum = (scalar(@$jobid_list) - 1); $jobNum > -1; $jobNum--) {
        my ($jobid, $seq_id, $error_count, $job_number) =
            split(/\s+/, $jobid_list->[$jobNum]);
        print_debug_message('_job_list_poll', 'jobNum: ' . $jobNum, 12);
        print_debug_message('_job_list_poll',
            'Job info: ' . $jobid_list->[$jobNum], 12);

        # Get job status.
        my $job_status = &rest_get_status($jobid);
        print_debug_message('_job_list_poll', 'Status: ' . $job_status, 12);

        # Fetch results and remove finished/failed jobs from list.
        if (
            !(
                $job_status eq 'RUNNING'
                    || $job_status eq 'PENDING'
                    || ($job_status eq 'ERROR'
                    && $error_count < $maxErrorStatusCount)
            )
        ) {
            if ($job_status eq 'ERROR' || $job_status eq 'FAILED') {
                print STDERR
                    "Warning: job $jobid failed for sequence $job_number: $seq_id\n";
            }
            # Duplicated getting results.
            #&get_results($jobid, $seq_id);
            splice(@$jobid_list, $jobNum, 1);
        }
        else {

            # Update error count, increment for new error or clear old errors.
            if ($job_status eq 'ERROR') {
                $error_count++;
            }
            elsif ($error_count > 0) {
                $error_count--;
            }

            # Update job tracking info.
            my $job_info_str = sprintf('%s %s %d %d',
                $jobid, $seq_id, $error_count, $job_number);
            $jobid_list->[$jobNum] = $job_info_str;
        }
    }
    print_debug_message('_job_list_poll', 'Num jobs: ' . scalar(@$jobid_list),
        11);
    print_debug_message('_job_list_poll', 'End', 1);
}

=head2 list_file_submit_job()

Submit multiple jobs using a file containing a list of entry identifiers as
input.

  &list_file_submit_job($list_filename)

=cut

sub list_file_submit_job {
    print_debug_message('list_file_submit_job', 'Begin', 1);
    my $filename = shift;

    # Open the file of identifiers.
    my $LISTFILE;
    if ($filename eq '-') { # STDIN.
        open($LISTFILE, '<-')
            or die 'Error: unable to STDIN (' . $! . ')';
    }
    else { # File.
        open($LISTFILE, '<', $filename)
            or die 'Error: unable to open file ' . $filename . ' (' . $! . ')';
    }

    # Job identifier tracking for parallel execution.
    my @jobid_list = ();
    my $job_number = 0;

    # Iterate over identifiers, submitting each job
    while (<$LISTFILE>) {
        my $line = $_;
        chomp($line);
        if ($line ne '') {
            &print_debug_message('list_file_submit_job', 'line: ' . $line, 2);
            if ($line =~ m/\w:\w/) {
                # Check this is an identifier
                my $seq_id = $line;
                print STDERR "Submitting job for: $seq_id\n"
                    if ($outputLevel > 0);
                $job_number++;
                my $job_id = &submit_job($seq_id, $seq_id);
                my $job_info_str =
                    sprintf('%s %s %d %d', $job_id, $seq_id, 0, $job_number);
                push(@jobid_list, $job_info_str);
            }
            else {
                print STDERR
                    "Warning: line \"$line\" is not recognised as an identifier\n";
            }

            # Parallel mode, wait for job(s) to finish to free slots.
            while ($params{'maxJobs'} > 1
                && scalar(@jobid_list) >= $params{'maxJobs'}) {
                &_job_list_poll(\@jobid_list);
                print_debug_message('list_file_submit_job',
                    'Remaining jobs: ' . scalar(@jobid_list), 1);
            }
        }
    }
    close $LISTFILE;

    # Parallel mode, wait for remaining jobs to finish.
    while ($params{'maxJobs'} > 1 && scalar(@jobid_list) > 0) {
        &_job_list_poll(\@jobid_list);
        print_debug_message('list_file_submit_job',
            'Remaining jobs: ' . scalar(@jobid_list), 1);
    }
    print_debug_message('list_file_submit_job', 'End', 1);
}


=head2 load_data()

Load sequence data from file or option specified on the command-line.

  &load_data();

=cut

sub load_data {
    print_debug_message('load_data', 'Begin', 1);
    my $retSeq;

    # Query sequence
    if (defined($ARGV[0])) {                  # Bare option
        if (-f $ARGV[0] || $ARGV[0] eq '-') { # File
            $retSeq = &read_file($ARGV[0]);
        }
        else { # DB:ID or sequence
            $retSeq = $ARGV[0];
        }
    }
    if ($params{'sequence'}) {                                      # Via --sequence
        if (-f $params{'sequence'} || $params{'sequence'} eq '-') { # File
            $retSeq = &read_file($params{'sequence'});
        }
        else { # DB:ID or sequence
            $retSeq = $params{'sequence'};
        }
    }
    print_debug_message('load_data', 'End', 1);
    return $retSeq;
}

=head2 load_params()

Load job parameters from command-line options.

  &load_params();

=cut

sub load_params {
    print_debug_message('load_params', 'Begin', 1);

    # Pass default values and fix bools (without default value)
    if ($params{'stype'} eq 'protein') {
        if (!$params{'task'}) {
            $params{'task'} = 'blastp'
        }
    }
    if ($params{'stype'} eq 'nucleotide') {
        if (!$params{'task'}) {
            $params{'task'} = 'blastn'
        }
    }
    if ($params{'stype'} eq 'vector') {
        if (!$params{'task'}) {
            $params{'task'} = 'blastn'
        }
    }

    if ($params{'stype'} eq 'protein') {
        if (!$params{'matrix'}) {
            $params{'matrix'} = 'BLOSUM62'
        }
    }
    if ($params{'stype'} eq 'nucleotide') {
        if (!$params{'matrix'}) {
            $params{'matrix'} = 'NONE'
        }
    }
    if ($params{'stype'} eq 'vector') {
        if (!$params{'matrix'}) {
            $params{'matrix'} = 'NONE'
        }
    }

    if (!$params{'alignments'}) {
        $params{'alignments'} = 50
    }

    if (!$params{'scores'}) {
        $params{'scores'} = 50
    }

    if (!$params{'exp'}) {
        $params{'exp'} = '10'
    }

    if (!$params{'dropoff'}) {
        $params{'dropoff'} = 0
    }

    if ($params{'stype'} eq 'nucleotide') {
        if (!$params{'match_scores'}) {
            $params{'match_scores'} = '1,-3'
        }
    }
    if ($params{'stype'} eq 'vector') {
        if (!$params{'match_scores'}) {
            $params{'match_scores'} = '1,-3'
        }
    }

    if (!$params{'gapopen'}) {
        $params{'gapopen'} = -1
    }

    if (!$params{'gapext'}) {
        $params{'gapext'} = -1
    }

    if ($params{'stype'} eq 'protein') {
        if (!$params{'filter'}) {
            $params{'filter'} = 'F'
        }
    }
    if ($params{'stype'} eq 'nucleotide') {
        if (!$params{'filter'}) {
            $params{'filter'} = 'T'
        }
    }
    if ($params{'stype'} eq 'vector') {
        if (!$params{'filter'}) {
            $params{'filter'} = 'T'
        }
    }

    if (!$params{'gapalign'}) {
        $params{'gapalign'} = 'true'
    }

    if (!$params{'compstats'}) {
        $params{'compstats'} = 'F'
    }

    if (!$params{'align'}) {
        $params{'align'} = 0
    }

    if (!$params{'transltable'}) {
        $params{'transltable'} = 1
    }

    print_debug_message('load_params', 'End', 1);
}

=head2 client_poll()

Client-side job polling.

  &client_poll($job_id);

=cut

sub client_poll {
    print_debug_message('client_poll', 'Begin', 1);
    my $jobid = shift;
    my $status = 'PENDING';

    # Check status and wait if not finished. Terminate if three attempts get "ERROR".
    my $errorCount = 0;
    while ($status eq 'RUNNING'
        || $status eq 'PENDING'
        || ($status eq 'ERROR' && $errorCount < 2)) {
        $status = rest_get_status($jobid);
        print STDERR "$status\n" if ($outputLevel > 0);
        if ($status eq 'ERROR') {
            $errorCount++;
        }
        elsif ($errorCount > 0) {
            $errorCount--;
        }
        if ($status eq 'RUNNING'
            || $status eq 'PENDING'
            || $status eq 'ERROR') {

            # Wait before polling again.
            usleep($checkInterval);
        }
    }
    print_debug_message('client_poll', 'End', 1);
    return $status;
}

=head2 get_results()

Get the results for a job identifier.

  &get_results($job_id);

=cut

sub get_results {
    print_debug_message('get_results', 'Begin', 1);
    my $jobid = shift;
    print_debug_message('get_results', 'jobid: ' . $jobid, 1);
    my $seq_id = shift;
    print_debug_message('get_results', 'seq_id: ' . $seq_id, 1) if ($seq_id);

    my $output_basename = $jobid;

    # Verbose
    if ($outputLevel > 1) {
        print 'Getting results for job ', $jobid, "\n";
    }

    # Check status, and wait if not finished
    client_poll($jobid);

    # Default output file names use JobId, however the name can be specified...
    if (defined($params{'outfile'})) {
        $output_basename = $params{'outfile'};
    }
    # Or use sequence identifer.
    elsif (defined($params{'useSeqId'} && defined($seq_id) && $seq_id ne '')) {
        $output_basename = $seq_id;

        # Make safe to use as a file name.
        $output_basename =~ s/\W/_/g;
    }

    # Use JobId if output file name is not defined
    else {
        unless (defined($params{'outfile'})) {
            #$params{'outfile'} = $jobid;
            $output_basename = $jobid;
        }
    }

    # Get list of data types
    my (@resultTypes) = rest_get_result_types($jobid);


    # Get the data and write it to a file
    if (defined($params{'outformat'})) {
        # Specified data type
        # check to see if there are multiple formats (comma separated)
        my $sep = ",";
        my (@multResultTypes);
        if ($params{'outformat'} =~ /$sep/) {
            @multResultTypes = split(',', $params{'outformat'});
        }
        else {
            $multResultTypes[0] = $params{'outformat'};
        }
        # check if the provided formats are recognised
        foreach my $inputType (@multResultTypes) {
            my $expectation = 0;
            foreach my $resultType (@resultTypes) {
                if ($resultType->{'identifier'} eq $inputType && $expectation eq 0) {
                    $expectation = 1;
                }
            }
            if ($expectation ne 1) {
                die 'Error: unknown result format "' . $inputType . '"';
            }
        }
        # if so get the files
        my $selResultType;
        foreach my $resultType (@resultTypes) {
            if (grep {$_ eq $resultType->{'identifier'}} @multResultTypes) {
                $selResultType = $resultType;
                my $result = rest_get_result($jobid, $selResultType->{'identifier'});
                if (defined($params{'outfile'}) && $params{'outfile'} eq '-') {
                    write_file($params{'outfile'}, $result);
                }
                else {
                    write_file(
                        $output_basename . '.'
                            . $selResultType->{'identifier'} . '.'
                            . $selResultType->{'fileSuffix'},
                        $result
                    );
                }
            }
        }
    }
    else { # Data types available
        # Write a file for each output type
        for my $resultType (@resultTypes) {
            if ($outputLevel > 1) {
                print STDERR 'Getting ', $resultType->{'identifier'}, "\n";
            }
            my $result = rest_get_result($jobid, $resultType->{'identifier'});
            if (defined($params{'outfile'}) && $params{'outfile'} eq '-') {
                write_file($params{'outfile'}, $result);
            }
            else {
                write_file(
                    $output_basename . '.'
                        . $resultType->{'identifier'} . '.'
                        . $resultType->{'fileSuffix'},
                    $result
                );
            }
        }
    }
    print_debug_message('get_results', 'End', 1);
}

=head2 read_file()

Read a file into a scalar. The special filename '-' can be used to read from
standard input (STDIN).

  my $data = &read_file($filename);

=cut

sub read_file {
    print_debug_message('read_file', 'Begin', 1);
    my $filename = shift;
    print_debug_message('read_file', 'filename: ' . $filename, 2);
    my ($content, $buffer);
    if ($filename eq '-') {
        while (sysread(STDIN, $buffer, 1024)) {
            $content .= $buffer;
        }
    }
    else {
        # File
        open(my $FILE, '<', $filename)
            or die "Error: unable to open input file $filename ($!)";
        while (sysread($FILE, $buffer, 1024)) {
            $content .= $buffer;
        }
        close($FILE);
    }
    print_debug_message('read_file', 'End', 1);
    return $content;
}

=head2 write_file()

Write data to a file. The special filename '-' can be used to write to
standard output (STDOUT).

  &write_file($filename, $data);

=cut

sub write_file {
    print_debug_message('write_file', 'Begin', 1);
    my ($filename, $data) = @_;
    print_debug_message('write_file', 'filename: ' . $filename, 2);
    if ($outputLevel > 0) {
        print STDERR 'Creating result file: ' . $filename . "\n";
    }
    if ($filename eq '-') {
        print STDOUT $data;
    }
    else {
        open(my $FILE, '>', $filename)
            or die "Error: unable to open output file $filename ($!)";
        syswrite($FILE, $data);
        close($FILE);
    }
    print_debug_message('write_file', 'End', 1);
}

=head2 usage()

Print program usage message.

  &usage();

=cut

sub usage {
    print STDERR <<EOF
EMBL-EBI NCBI Blast Perl Client:

Sequence similarity search with NCBI Blast.

[Required (for job submission)]
  --email               E-mail address.
  --program             The BLAST program to be used for the Sequence Similarity
                        Search.
  --stype               Indicates if the sequence is protein or DNA/RNA.
  --sequence            The query sequence can be entered directly into this form.
                        The sequence can be in GCG, FASTA, EMBL (Nucleotide only),
                        GenBank, PIR, NBRF, PHYLIP or UniProtKB/Swiss-Prot (Protein
                        only) format. A partially formatted sequence is not
                        accepted. Adding a return to the end of the sequence may
                        help certain applications understand the input. Note that
                        directly using data from word processors may yield
                        unpredictable results as hidden/control characters may be
                        present.
  --database            Database.

[Optional]
  --task                Task option (only selectable for blastn).
  --matrix              (Protein searches) The substitution matrix used for scoring
                        alignments when searching the database.
  --alignments          Maximum number of match alignments reported in the result
                        output.
  --scores              Maximum number of match score summaries reported in the
                        result output.
  --exp                 Limits the number of scores and alignments reported based on
                        the expectation value. This is the maximum number of times
                        the match is expected to occur by chance.
  --dropoff             The amount a score can drop before gapped extension of word
                        hits is halted.
  --match_scores        (Nucleotide searches) The match score is the bonus to the
                        alignment score when matching the same base. The mismatch is
                        the penalty when failing to match.
  --gapopen             Penalty taken away from the score when a gap is created in
                        sequence. Increasing the gap openning penalty will decrease
                        the number of gaps in the final alignment.
  --gapext              Penalty taken away from the score for each base or residue
                        in the gap. Increasing the gap extension penalty favors
                        short gaps in the final alignment, conversly decreasing the
                        gap extension penalty favors long gaps in the final
                        alignment.
  --filter              Filter regions of low sequence complexity. This can avoid
                        issues with low complexity sequences where matches are found
                        due to composition rather than meaningful sequence
                        similarity. However in some cases filtering also masks
                        regions of interest and so should be used with caution.
  --seqrange            Specify a range or section of the input sequence to use in
                        the search. Example: Specifying '34-89' in an input sequence
                        of total length 100, will tell BLAST to only use residues 34
                        to 89, inclusive.
  --gapalign            This is a true/false setting that tells the program the
                        perform optimised alignments within regions involving gaps.
                        If set to true, the program will perform an alignment using
                        gaps. Otherwise, if it is set to false, it will report only
                        individual HSP where two sequence match each other, and thus
                        will not produce alignments with gaps.
  --wordsize            Word size for wordfinder algorithm.
  --taxids              Specify one or more TaxIDs so that the BLAST search becomes
                        taxonomically aware.
  --compstats           Use composition-based statistics.
  --align               Formating for the alignments.
  --transltable         Query Genetic code to use in translation.

[General]
  -h, --help            Show this help message and exit.
  --asyncjob            Forces to make an asynchronous query.
  --title               Title for job.
  --status              Get job status.
  --resultTypes         Get available result types for job.
  --polljob             Poll for the status of a job.
  --pollFreq            Poll frequency in seconds (default 3s).
  --jobid               JobId that was returned when an asynchronous job was submitted.
  --outfile             File name for results (default is JobId; for STDOUT).
  --multifasta          Treat input as a set of fasta formatted sequences.
  --useSeqId            Use sequence identifiers for output filenames.
                        Only available in multi-fasta and multi-identifier modes.
  --maxJobs             Maximum number of concurrent jobs. Only
                        available in multifasta or list file modes.
  --outformat           Result format(s) to retrieve. It accepts comma-separated values.
  --params              List input parameters.
  --paramDetail         Display details for input parameter.
  --quiet               Decrease output.
  --verbose             Increase output.
  --version             Prints out the version of the Client and exit.
  --baseUrl             Base URL. Defaults to:
                        https://www.ebi.ac.uk/Tools/services/rest/ncbiblast

Synchronous job:
  The results/errors are returned as soon as the job is finished.
  Usage: perl $scriptName --email <your\@email.com> [options...] <SeqFile|SeqID(s)>
  Returns: results as an attachment

Asynchronous job:
  Use this if you want to retrieve the results at a later time. The results
  are stored for up to 24 hours.
  Usage: perl $scriptName --asyncjob --email <your\@email.com> [options...] <SeqFile|SeqID(s)>
  Returns: jobid

Check status of Asynchronous job:
  Usage: perl $scriptName --status --jobid <jobId>

Retrieve job data:
  Use the jobid to query for the status of the job. If the job is finished,
  it also returns the results/errors.
  Usage: perl $scriptName --polljob --jobid <jobId> [--outfile string]
  Returns: string indicating the status of the job and if applicable, results
  as an attachment.

Further information:
  https://www.ebi.ac.uk/Tools/webservices and
    https://github.com/ebi-wp/webservice-clients

Support/Feedback:
  https://www.ebi.ac.uk/support/
EOF
}

=head1 FEEDBACK/SUPPORT

Please contact us at L<https://www.ebi.ac.uk/support/> if you have any
feedback, suggestions or issues with the service or this client.

=cut
