#!/usr/bin/env perl

######################################################################
#2006-06-26                                                          #
#Olivier POIROT, poirot [AT] igs.cnrs-mrs.fr, Labo IGS, UPR2589      #
#Sebastien MORETTI, moretti.sebastien [AT] gmail.com, SIB            #
#Web interface for Tcoffee and other programs                        #
#Need feedback.cgi, configuration_file.txt                           #
######################################################################

use CGI;
use Time::localtime;
use Time::Local;
use Sys::Hostname;

#############################################################################################################################################
use env_param;                                                                                                                              #
my $location = &hostname();                                                                                                                   #
&env_param::variables($location);                                                                                                           #
#############################################################################################################################################

my $web_base_link_tmp = $env_param::web_base_link_tmp;
my $web_base_images   = $env_param::web_base_images;
my $web_base          = $env_param::web_base;
my $dir_cgi           = $env_param::dir_cgi;
my $pg_source         = $env_param::pg_source;
my $scratch_area      = $env_param::scratch_area;
my $tmpOldFiles       = $env_param::tmpOldFiles;


$|           = 1;
$query       = new CGI;
$child       = $query->param('child');
$result_file = $query->param('result_file');
$pg_source   = $query->param('pg_source');
$email       = $query->param('email');
$level       = $query->param('level');
$mode        = $query->param('mode');
$time_init   = $query->param('time_init');

my $appli          = substr($mode,0,3).substr($level,0,1);
$scratch_area      = "$scratch_area/$appli";
$web_base_link_tmp = "$web_base_link_tmp/$appli";

# see if child process is done
$status = `ps up $child |grep -v PID`;
if ( $status =~ /$pg_source/ && length($pg_source) > 0 ){
        # child process is still running, set up this output page to refresh so we can check things again
        $me = $query->self_url;
        $me =~ s/alia\.igs\.univ/www\.igs\.cnrs/;
        print $query->header (-Refresh=>'3; URL=',"$me");
        $done = 0;
}
else {
    if ( -e "$scratch_area/$result_file" && -s "$scratch_area/$result_file" ){
       print $query->header (-Refresh=>"2; URL=$web_base_link_tmp/$result_file");
       $done = 1;
       $me = $query->self_url;
       $me =~ s/alia\.igs\.univ/www\.igs\.cnrs/;
    }
    else{
       &ifResultFileIsOlderThenRemoved();
       exit 1;
    }
}
# start the output page
#print $query->start_html(-title=>"$mode :: $level", -meta=>{'Prev'=>"./$pg_source?stage1=1&amp;daction=${mode}::${level}"}, -BGCOLOR=>'#F6F6FF');
print "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">
<html>
<head>
    <title>$mode :: $level</title>
    <link rel=\"prev\" href=\"./$pg_source?stage1=1&amp;daction=${mode}::${level}\">
    <link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"$web_base_images/tcoffee.ico\">
</head>

<body bgcolor=\"#F6F6FF\">\n";
if ( $done ) {
    print 'Your job is finished';
    #test IP submision
    #system("rm /usr/local/httpd/htdocs/Tcoffee/Tmp/$who");
    #print "<br/>child=$child<br/>status=$status<br/>me=$me<br/>te=$te<br/>result_file=$result_file<br/><h1>done=$done</h1><br/>pgsource=$pg_source" ;
}
else {
    print "\n<DIV>Processing, please wait...<BR/><BR/>\n<I>You can bookmark this page <BR/>to retrieve your results when done</I>\n";
    print "\n<CENTER><img src=\"$web_base_images/l5.gif\" alt=\"Be patient\" /></CENTER>\n";
    $tm = localtime;
    ( $SEC, $MIN, $HOUR, $DAY, $MONTH, $YEAR ) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, $tm->mon, $tm->year);
    $time_now = timelocal($SEC,$MIN,$HOUR,$DAY,$MONTH,$YEAR);
    $diff     =  $time_now - $time_init;
    printf("time: %s seconds<br/></DIV>", $diff);
    if( &ValidEmailAddr($email) ){
        print "\n<DIV>";
        print $query->start_form(-action=>"$web_base/$dir_cgi/$pg_source");
        print 'You can close this window, you will receive your results by Email';
        print "\n<br/>Return to the server: ";
        print $query->hidden('stage1','1');
        print $query->submit(-name=>'daction',-value=>$mode.'::'.$level);
        print $query->endform;
        print "\n</DIV>\n";
    }
}
print $query->end_html;
exit 0;

sub ValidEmailAddr {

    my $mail = shift;
    return 0 if ( $mail !~ /^[\w\.\-\_]+\@[\w\.\-\_]+\.[A-Za-z][A-Za-z][A-Za-z]?$/ );
    return 1;
}

sub ifResultFileIsOlderThenRemoved{
    print "Content-type: text/html\n
<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">
<HTML>
<HEAD>
  <TITLE>Result file too old !</TITLE>
</HEAD>

<BODY>
  <DIV align=center style=\"font-size:x-large\;\"><br/><br/><br/><br/>
The result file you are requesting is too old (older than $tmpOldFiles days).<br/>
So, it has been removed from our server!<br/><br/>
---&gt; <A href=\"$web_base/$dir_cgi/$pg_source\"> Homepage </A>
  </DIV>
</BODY>
</HTML>
";

    return;
}

