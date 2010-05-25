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
my $pg_source         = $env_param::pg_source;

my $q        = new CGI;

my $log_file = $q->param('log_file');

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

print "<A HREF=\"evol2.cgi?log_file=$log_file\">Retour</A></td>
        </tr>
      </table>
    </td>
  </tr>
  <tr>
    <td align=center>
      <table width=730 align=center border=0 cellpadding=0 cellspacing=0 summary=\"hits per day\">
        <tr>
";
my $la_date = $q->param('date');
my $serv    = '';
if ( defined($q->param('serv')) ){
    $serv   = $q->param('serv');
}
&analyse_date($la_date,$serv);

#print &print_mise_en_page_fin();
print $q->end_html;

sub analyse_date{

    my ( $la_date, $serv ) = @_;

    $/ = "<user-call>";
    open ( my $EVOL, '<', "$log_file" );
    LIST_IP_COUNTRY:
    while (<$EVOL>){
        if(/$la_date/){
            if( $serv eq 'nul' ){
                $nb_tot++;
                @rec  = split(/\n/,$_);
                $ip   = $rec[2];
                $ip   =~ /<remote-host>(.*)<\//;
                $ip   = $1;
                $file = $rec[1];
                $tab_ass_ip_occ{$ip}++;
                if( ! $tab_ass_ip_pays{$ip} ){
                    $com  = "$home_cgi/ip.pl $ip";
                    $pays =`$com`;
                    $tab_ass_ip_pays{$ip} = $pays;
                }
            }
            else {
                if(/$serv/){
                    $nb_tot++;
                    @rec  = split(/\n/,$_);
                    $ip   = $rec[2];
                    $ip   =~ /<remote-host>(.*)<\//;
                    $ip   = $1;
                    $file = $rec[1];
                    $tab_ass_ip_occ{$ip}++;
                    if( ! $tab_ass_ip_pays{$ip} ){
                        $com = "$home_cgi/ip.pl $ip";
                        $pays=`$com`;
                        $tab_ass_ip_pays{$ip} = $pays;
                    }
                }
            }
        }
    }

    close $EVOL;
    $/ = "\n";
    print "          <td align=center colspan=3><hr/>$la_date, $nb_tot hits<hr/></td>\n        </tr>
      </table>
    </td>
  </tr>
  <tr>
    <td align=center>
      <table width=730 align=center border=0 cellpadding=0 cellspacing=0 summary=\"back\">\n";
    LINKS:
    foreach $k(sort {$tab_ass_ip_occ{$b} <=> $tab_ass_ip_occ{$a}}keys %tab_ass_ip_occ) {
        chomp($tab_ass_ip_pays{$k});
        $tab_ass_ip_occ{$k} = '&nbsp;'x(4-length($tab_ass_ip_occ{$k})).$tab_ass_ip_occ{$k} if ( length($tab_ass_ip_occ{$k}) < 4 );
        my $noSpaceInIP     = $k;
        $noSpaceInIP        =~ s/ /\%20/g;
        print "        <tr style=\"font-size:small\">
          <td align=right width=\"33%\"><A HREF=\"ip.cgi?ip=$noSpaceInIP&amp;date=$la_date\" style=\"text-decoration:none\">$k</A></td>
          <td align=center width=\"33%\">$tab_ass_ip_pays{$k}</td>
          <td align=left width=\"33%\" style=\"font-family:courier\">$tab_ass_ip_occ{$k}  hits</td>
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
    <td align=center>
      <table width=730 align=center border=0 cellpadding=0 cellspacing=0 summary=\"back\">
        <tr>
          <td align=center><br/>";
    return $chaine;
}

sub print_mise_en_page_fin() {

    $chaine = '';
    $chaine.= '</td>';
    $chaine.= '</tr>';
    $chaine.= "</table>\n";
    $chaine.= '</td>';
    $chaine.= "<td><img src=\"$web_base_images/px_trans.gif\" width=1 height=1 alt=\"background\" /></td>\n";
    $chaine.= "<td width=1 bgcolor=\"#000000\"><img src=\"$web_base_images/px_trans.gif\" width=1 height=1 alt=\"background\" /></td>\n";
    $chaine.= "<td width=20 background=\"$web_base_images/patterndroite.gif\">&nbsp;</td>\n";
    $chaine.= '</tr>';
    $chaine.= "</table>\n";
    return $chaine;
}

sub sort_numer {  $a <=> $b;}

