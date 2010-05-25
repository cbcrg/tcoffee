#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $q = new CGI;


my $passw = $q->param('secret') || '';
if ( $passw ne 'GoodPassW4fun!' ){
    &sortir();
    exit 0;
}


my $oldmessages = $q->param('oldmes')   || '';
my $messages    = $q->param('messages') || '';
if ( $oldmessages eq $messages ){
    &sortir();
    exit 0;
}
if ( $messages eq '' ){
    unlink('Info.txt');
}
else{
    open( my $MES, '>', 'Info.txt');
    print {$MES} $messages;
    close $MES;
}
print "Content-Type\: text/html

<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">
<HTML>
<HEAD>
  <TITLE></TITLE>
  <meta http-equiv=\"refresh\" content=\"8; URL=/Tcoffee/tcoffee_cgi/index.cgi\" />
</HEAD>\n
<BODY>
Former message:<br/>[ $oldmessages ]<br/><br/><br/>New message:<br/>[ $messages ]
</BODY>
</HTML>\n";





sub sortir{
    print "Content-Type\: text/html

<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">
<HTML>
<HEAD>
  <TITLE></TITLE>
  <meta http-equiv=\"refresh\" content=\"0; URL=/Tcoffee/tcoffee_cgi/index.cgi\" />
</HEAD>\n

<BODY>
</BODY>
</HTML>\n";
    return;
}

exit 0;

