#!/usr/bin/env perl

######################################################################
#2006-06-26                                                          #
#Olivier POIROT, poirot [AT] igs.cnrs-mrs.fr, Labo IGS, UPR2589      #
#Sebastien MORETTI, moretti.sebastien [AT] gmail.com, SIB            #
#Web interface for Tcoffee and other programs                        #
#Need feedback.cgi, configuration_file.txt                           #
######################################################################

use File::stat;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Time::Local;
use GD::Graph::bars;
use Sys::Hostname;


#############################################################################################################################################
use env_param;                                                                                                                              #
my $location = &hostname();                                                                                                                   #
&env_param::variables($location);                                                                                                           #
#############################################################################################################################################

$| = 1;

my $tcoffee_homepage  = $env_param::tcoffee_homepage;
my $home_cgi          = $env_param::home_cgi;
my $web_cgi_base      = $env_param::web_cgi_base;
my $web_base_doc      = $env_param::web_base_doc;
my $web_base_link_tmp = $env_param::web_base_link_tmp;
my $web_base_images   = $env_param::web_base_images;
my $logo              = $env_param::logo;
my $webmaster         = $env_param::webmaster;
my $scratch_area      = $env_param::scratch_area;

my $pg_source         = $env_param::pg_source;

my $q = new CGI;

print $q->header;
print "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">
<HTML>
<HEAD>
  <TITLE>Tcoffee</TITLE>
  <META name=\"keywords\" content=\"Multiple sequence alignment,T-COFFEE,TCoffee,t_coffee,MSA,alignment,protogene,proto-gene,3dcoffee,3d-coffee,expresso,APDB,iRMSD,IGS,Vital-IT\" />
  <META http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\" />
</HEAD>

<body background=\"$web_base_images/pattern.gif\" style=\"margin: 0\">
";
#print $q->start_html(-background=>"$web_base_images/pattern.gif",-leftmargin=>'0',-topmargin=>'0',marginwidth=>'0',-marginheight=>'0',-title=>"Tcoffee",-meta=>{'keywords'=>'Multiple sequence,alignment T-COFFEE'});
print &print_mise_en_page_deb();

$la_date = $q->param('date');
$ip      = $q->param('ip');
print "<a href='date.cgi?date=$la_date'>Retour</a>
          </td>
        </tr>
      </table>
    </td>
  </tr>
  <tr>
    <td align=center>
      <table width='730' align='center' border='0' cellpadding='0' cellspacing='0' summary=\"job(s) for ...\">
        <tr>
          <td align='center'>\n";
&analyse_pid_date($pid,$la_date);

#print &print_mise_en_page_fin();
print $q->end_html;

sub analyse_pid_date($ip,$la_date) {
    $pays = `$home_cgi/ip.pl $ip`;
    $/    = "<user-call>";
    open ( my $EVOL, '<', "$home_cgi/userlog.xml" );
    LIST_IP:
    while (<$EVOL>) {
        if( /$la_date/ && /$ip/ ){
            $nb_tot++;
            #print $_,"<hr/>";
            @rec = split(/\n/, $_);
            $ip  = $rec[2];
            $ip  =~ /<remote-host>(.*)<\//;
            $ip  = $1;
            $d_date = $rec[6];
#           $mail   = $rec[3];
            $file   = $rec[1];
            $file   =~ /<file>(.*?)<\//;
            $file   = $1;
            $tab_ass_ip_date{$d_date} = $file;
        }
    }

    close $EVOL;
    $/ = "\n";
    chomp($pays);
    print "            <hr/>$la_date, $ip, $pays, $nb_tot hits<hr/>
          </td>
        </tr>
      </table>
    </td>
  </tr>
  <tr>
    <td align=center>
      <table width='730' align='center' border='0' cellpadding='0' cellspacing='0' summary=\"hits per day per IP\">
";
    my @tlist     = keys(%tab_ass_ip_date);
    my @time_list = reverse( sort(@tlist) );
    LINKS:
    for (my $y=0; $y <=$#time_list; $y++){
        my $crash   = '<td>&nbsp; </td>';
        my ($appli) = $tab_ass_ip_date{$time_list[$y]} =~ /tcf(....).*/;
        if ( -e "$scratch_area/$appli/$tab_ass_ip_date{$time_list[$y]}.fatal" ){
            $crash = '<td><b>Crash</b></td>';
        }
        my $datum = $time_list[$y];
        $datum    =~ s/ *<\/?time-stamp> *//g;
        print "        <tr style='font-size:small'>
          <td align='right' width='50%'>$datum &nbsp;&nbsp;</td>
          <td align='left' width='25%'><a href=\"$web_base_link_tmp/$appli/$tab_ass_ip_date{$time_list[$y]}.file_result.html\" style=\"text-decoration:none\">$tab_ass_ip_date{$time_list[$y]}</a></td>
          $crash
        </tr>\n";
    }
    print "      </table><br/>
    </td>
  </tr>
</table>\n";
  #print "<table border=1><tr><td colspan=15 align=center><img src=$web_base_link_tmp/picture.png alt=\"picture\" /></td></tr><tr></table>\n";
    return;
}

######################################################################
#Look html ###########################################################
######################################################################
sub print_mise_en_page_deb() {

    $chaine = '';
    $chaine.= "\n<table summary=\"menu\" align=\"center\" bgcolor=\"#ffffff\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"792\">
  <tr valign=\"top\">
    <td width=\"750\" align=\"center\">
      <table summary=\"logo\" align=\"center\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"792\">
        <tr>
          <td colspan=\"5\" align=\"center\" nowrap=\"nowrap\">$logo
          </td>
        </tr>
        <tr>
           <td align=\"left\" bgcolor=\"#000000\" width=\"1\">
             <a href=\"$web_cgi_base/evol2.cgi?log_file=userlog.xml\" style=\"text-decoration: none; color: #000000;\">.</a>
           </td>
           <td align=\"right\" bgcolor=\"#000000\" height=\"6\" valign=\"top\" style=\"color: #ffffff;\" colspan=\"2\" width=\"791\">
             <a href=\"$web_cgi_base/$pg_source\" style=\"color: #ffffff;\">home</a>&nbsp;|&nbsp;
             <a href=\"$tcoffee_homepage\" style=\"color: #ffffff;\">references</a>&nbsp;|&nbsp;
             <a href=\"$web_base_doc/doc3.html\" target=\"_blank\" style=\"color: #ffffff;\">help</a>&nbsp;|&nbsp;\n";
    $chaine .= '             '.&env_param::transpose_eMails_delayed($webmaster, "$web_base_images/mail1.gif")."&nbsp;&nbsp;&nbsp;
             <a href=\"$web_cgi_base/${pg_source}?load_cfg_file=1\" style=\"text-decoration: none; color: #000000;\">.</a>
           </td>
        </tr>
      </table>
    </td>
  </tr>
  <tr>
    <td align=\"center\">
      <table width=\"730\" align=\"center\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\" summary=\"back\">
        <tr>
          <td align=\"center\"><br/>";
    return $chaine;
}

sub print_mise_en_page_fin() {

    $chaine = '';
    $chaine.= '</td>';
    $chaine.= '</tr>';
    $chaine.= "</table>\n";
    $chaine.= '</td>';
    $chaine.= "<td><img src=\"$web_base_images/px_trans.gif\" width=\"1\" height=\"1\" alt=\"background\" /></td>\n";
    $chaine.= "<td width=\"1\" bgcolor=\"#000000\"><img src=\"$web_base_images/px_trans.gif\" width=1 height=1 alt=\"background\" /></td>\n";
    $chaine.= "<td width=\"20\" background=\"$web_base_images/patterndroite.gif\">&nbsp;</td>\n";
    $chaine.= '</tr>';
    $chaine.= "</table>\n";
    return $chaine;
}

#sub sort_numer {  $a <=> $b;}

