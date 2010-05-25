#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $q = new CGI;

my $oldmessages = '';
if ( -e 'Info.txt' ){
    if ( -z 'Info.txt' ){
        unlink('Info.txt');
    }
    else{
        open( my $MES, '<', 'Info.txt' );
        $oldmessages=join("\n", map( m/^([\#\@].*)$/, <$MES>));
        close $MES;
    }
}

my @servers;
if ( -e 'configuration_file.txt' ){
    open( my $CONF, '<', 'configuration_file.txt' );
    @servers=map( m/^server::([\w\(\)]+::\w+)::/,<$CONF> );
    close $CONF;
}

print "Content-Type\: text/html

<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">
<HTML>
<HEAD>
  <TITLE></TITLE>
</HEAD>\n
<BODY>
<FORM name=\"messenger\" method=\"post\" action=\"applyMessages.cgi\">
<DIV align=\"center\">

",$q->textarea(-name=>'messages',-columns=>130,-rows=>15,-value=>$oldmessages),"

<br/>
<br/>

Password: ",$q->password_field(-name=>'secret',-maxlength=>20),"
<br/>
<br/>

",$q->hidden(-name=>'oldmes',-value=>$oldmessages),"
<input type=\"submit\" value=\"Apply these messages\" />
<input type=\"reset\" value=\"Re-Use the current messages file\" />

</DIV>
</FORM>
</BODY>
</HTML>\n";

exit 0;

