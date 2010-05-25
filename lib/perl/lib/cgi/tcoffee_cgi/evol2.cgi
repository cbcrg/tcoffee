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
use GD::Graph::pie;
use Date::Calc qw(Delta_Days);
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
my $scratch_area      = $env_param::scratch_area;
my $logo              = $env_param::logo;
my $webmaster         = $env_param::webmaster;
my $pg_source         = $env_param::pg_source;

my @t_mois = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec');


#my $life = 72;

my $q = new CGI;

my $log_file = $q->param('log_file') || 'userlog.xml';
my $monf     = $q->param('monf');
my $service  = $q->param('serv') || '';
   $service  = "&amp;serv=$service" if ($service ne '');

my $crash = $log_file eq 'userlog.xml' ? "<a href='evol2.cgi?log_file=crashlog.xml$service'>crashed jobs</a>" :
                                         "<a href='evol2.cgi?log_file=userlog.xml$service'>full jobs</a>";

print $q->header;
print "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">
<html lang='en'>
<head>
    <title>Tcoffee</title>
    <meta http-equiv='Content-Type' content='text/html; charset=utf-8'>
</head>

<body background='$web_base_images/pattern.gif' style='margin: 0'>
";
#print $q->start_html(-background=>"$web_base_images/pattern.gif",-leftmargin=>'0',-topmargin=>'0',marginwidth=>'0',-marginheight=>'0',-title=>"Tcoffee",-meta=>{'keywords'=>'Multiple sequence alignment, T-COFFEE'});
print &print_mise_en_page_deb();
print "<hr>\n";
&evol($monf);
print "        <tr>
          <td align='center'><br>
            $crash<hr>
          </td>
        </tr>
      </table>
    </td>
  </tr>
  <tr>
    <td align='center'>
      <table width='730' align='center' border='0' cellpadding='0' cellspacing='0' summary='jobs per server'>
        <tr>
          <td align='center'>\n";
LINKS:
while ( ($k,$v) = each(%tab_serveur_oc) ){
    $ks = substr($k,3);
    print "            &nbsp;&nbsp;<a href='evol2.cgi?serv=$ks&amp;log_file=userlog.xml' style='text-decoration:none;'>$ks</a>&nbsp;&nbsp;\n";
}
print "          <hr></td>
        </tr>
      </table>
    </td>
  </tr>
  <tr>
    <td align='center'>
      <table width='730' align='center' border='0' cellpadding='0' cellspacing='0' summary='jobs per day'>\n";
#print "<table width=80%><tr><td>";
@tab_list_date = reverse @tab_list_date;

my $link_serv  = '&amp;serv=nul';
if ( defined($q->param('serv')) ){
    $link_serv = '&amp;serv='.$q->param('serv');
}
LIST_DATE:
foreach $d(@tab_list_date){
    $tab_ass{$d}='&nbsp;'x(4-length($tab_ass{$d})).$tab_ass{$d} if ( length($tab_ass{$d}) < 4 );
    if( $monf == 13 ) {
        print "        <tr style='font-size:small;'>
          <td align='right' width='50%'>
            mois $d
          </td>
          <td align='left' width='50%'><code>&nbsp;$tab_ass{$d}</code> hits
          </td>
        </tr>\n";
    }
    else {
        print "        <tr style='font-size:small;'>
          <td align='right' width='50%'>
            <a href='date.cgi?date=$d&amp;log_file=$log_file"."$link_serv' style='text-decoration:none;'>$d</a>
          </td>
          <td align='left' width='50%'><code>&nbsp;$tab_ass{$d}</code> hits
          </td>
        </tr>\n";
    }
}
print "      </table><br>
    </td>
  </tr>
</table>\n";


#print &print_mise_en_page_fin();
print $q->end_html;

sub evol($monf){

    if( $monf ) {
        ($mm,$yy) = split(/_/,$monf);
        if( $monf == 0 ) {
            print 'Year ',$yy;
        }
        else {
            print $t_mois[$mm-1],' ', $yy;
        }
    }
    else {
        print '31 last days';
    }
    my (@tab_fic,@tab_nb_fic);
    open ( my $EVOL, '<', "$log_file");
    while (<$EVOL>){
        if(/<file>/){
            $tag = 0;
            $_ = ~ /<file>(.*)<\//;
            $file = $1;
            $file =~ /^(......).*/;
            $serveur = $1;
            if (defined($q->param('serv'))){
                $serv = $q->param('serv');
                #print "<hr> $serv";
                if ($serveur =~ /$serv/){
                    $tab_serveur_oc{$serveur}++;
                    $tag = 1;
                }
            }
            else {
                $tab_serveur_oc{$serveur}++;
                $tag = 1;
                #print " $tag";
            }
        }
        if (/time-stamp/ && $tag == 1){
            $_ =~ /.* (.*)<\//;
            $date = $1;
            $date =~ /(.*)\/(.*)\/(.*)/;
            $day1 = $1;
            $month1 = $2;
            $year1 = $3;
            ($day2,$month2,$year2) = (localtime)[3,4,5];
            $month2 += 1;
            $year2 += 1900;
            $month1 =~ s/^0//;
            $delta_days = &Delta_Days( $year1, $month1, $day1, $year2, $month2, $day2);
            #$temps = timelocal(1,1,1,$day1,$month1,$year1);
            #($se,$mi,$ho,$jm,$m,$a,$js,$ja,$he) = localtime($temps);
            $tab_ass_month{$month1.'_'.$year1}++;
            #print  $month1,"--",$monf,"<br>";
            if( $monf ) {
                ($mm,$yy) = split(/_/,$monf);
                if($mm == 0) {
                    if($yy == $year1) {
                        if (! $tab_ass{$month1.'_'.$year1}){
                            $tab_list_date[$cpt] = $month1.'_'.$year1;
                            $cpt++;
                        }
                        $tab_ass{$month1.'_'.$year1}++;
                        $tab_cpt_date[$cpt]++;
                        $tab_ass_cpt_date{$cpt} = $month1.'_'.$year1;
                    }
                }
                elsif($month1 == $mm && $year1 == $yy) {
                    if (! $tab_ass{$date}){
                        $tab_list_date[$cpt] = $date;
                        $cpt++;
                    }
                    $tab_ass{$date}++;
                    $tab_cpt_date[$cpt]++;
                    $tab_ass_cpt_date{$cpt} = $date;
                }
            }
            elsif($delta_days <=31) {
                if (! $tab_ass{$date}){
                    $tab_list_date[$cpt] = $date;
                    $cpt++;
                }
                $tab_ass{$date}++;
                $tab_cpt_date[$cpt]++;
                $tab_ass_cpt_date{$cpt} = $date;
            }
        }
    }
    close $EVOL;

    my (@tab_serv,@tab_nb_serv);
    SERVERS:
    foreach $k(keys %tab_serveur_oc){
        push (@tab_serv, $k);#print "$k nb= $tab_cpt{$k}<br>";
        push (@tab_nb_serv, $tab_serveur_oc{$k});
    }

    DATES:
    foreach $nb_date(@tab_cpt_date){
        push (@tab_fic, $tab_ass_cpt_date{$z});
        push (@tab_nb_fic, $nb_date);
        #print "$z  $nb_date <br>\n";
        $z++;
    }
    shift @tab_nb_fic;
    shift @tab_fic;
    my @data_fic = ([@tab_fic],[@tab_nb_fic]);
    my $mygraph = GD::Graph::bars->new(500, 200);
    $mygraph->set(
                  x_labels_vertical => '1',
                  x_label           => 'Days/Month',
                  y_label           => 'nb_hits',
                  long_ticks        => '1',
                  title             => '# Hits',
                  y_tick_number     => 5,
                  ) or warn $mygraph->error;

    my $myimage = $mygraph->plot(\@data_fic) or die $mygraph->error;

    open( my $PICTURE, '>', "$scratch_area\/picture.png") or die("Cannot open file for writing");
    binmode $PICTURE;
    print {$PICTURE} $myimage->png;
    close $PICTURE;

    if(!$monf) {
        my @data_serv = ([@tab_serv],[@tab_nb_serv]);
        my $mygraph2  = GD::Graph::pie->new(200, 200);
        $mygraph2->set(
               'title'       => 'Serveurs',
               '3d'          => 1,
        ) or warn $mygraph->error;

        my $myimage2 = $mygraph2->plot(\@data_serv) or die $mygraph2->error;

        open( my $PICTURE2, '>', "$scratch_area\/picture2.png") or die("Cannot open file for writing");
        binmode $PICTURE2;
        print {$PICTURE2} $myimage2->png;
        close $PICTURE2;
    }
    MONTHS:
    foreach $k(sort keys %tab_ass_month) {
        ($mm,$yy) = split(/_/,$k);
        if(length($mm) == 1) {
            $mm = '0'.$mm;
        }
        $t_ass_month_sorted{$yy."_".$mm} = $k;
    }
    SORTED_MONTHS:
    foreach $k(sort keys %t_ass_month_sorted) {
        ($mm,$yy) = split(/_/,$t_ass_month_sorted{$k});
        $t_ass_deja_vu_y{$yy}++;
        if( $t_ass_deja_vu_y{$yy} == 1 ) {
            $liste_month = sprintf("%s <a href='evol2.cgi?log_file=userlog.xml&amp;monf=0_%s'><strong>%s</strong></a>", $liste_month, $yy, $yy);
        }
        $liste_month = sprintf("%s <a href='evol2.cgi?log_file=userlog.xml&amp;monf=%s'>%s</a>", $liste_month, $t_ass_month_sorted{$k}, $mm);
    }
    #$liste_month = sprintf("%s <a href='evol2.cgi?log_file=userlog.xml&amp;monf=13'>2006</a>", $liste_month);
    print "            <table border='1' summary='statistics graphics'>
              <tr>
                <td colspan='15' align='center'><img src='$web_base_link_tmp/picture.png' alt='by date'></td>";
    if( !$monf ) {
        print "<td><img src='$web_base_link_tmp/picture2.png' alt='by server'></td>";
    }
    print "</tr>
            </table>
             $liste_month
            </center><hr>
          </td>
        </tr>
      </table>
    </td>
  </tr>
  <tr>
    <td align='center'>
      <table width='730' align='center' border='0' cellpadding='0' cellspacing='0' summary='crashed jobs'>\n";
    return;
}


######################################################################
#Look html ###########################################################
######################################################################
sub print_mise_en_page_deb() {

    $chaine = '';
    $chaine.= "\n<table summary='menu' align='center' bgcolor='#ffffff' border='0' cellpadding='0' cellspacing='0' width='792'>
  <tr valign='top'>
    <td width='750' align='center'>
      <table summary='logo' align='center' border='0' cellpadding='0' cellspacing='0' width='792'>
        <tr>
          <td colspan='5' align='center' nowrap='nowrap'>$logo
          </td>
        </tr>
        <tr>
           <td align='left' bgcolor='#000000' width='1'>
             <a href='$web_cgi_base/evol2.cgi?log_file=userlog.xml' style='text-decoration: none; color: #000000;'>.</a>
           </td>
           <td align='right' bgcolor='#000000' height='6' valign='top' style='color: #ffffff;' colspan='2' width='791'>
             <a href='$web_cgi_base/$pg_source' style='color: #ffffff;'>home</a>&nbsp;|&nbsp;
             <a href='$tcoffee_homepage' style='color: #ffffff;'>references</a>&nbsp;|&nbsp;
             <a href='$web_base_doc/doc3.html' target='_blank' style='color: #ffffff;'>help</a>&nbsp;|&nbsp;\n";
    $chaine .= '             '.&env_param::transpose_eMails_delayed($webmaster, "$web_base_images/mail1.gif")."&nbsp;&nbsp;&nbsp;
             <a href='$web_cgi_base/${pg_source}?load_cfg_file=1' style='text-decoration: none; color: #000000;'>.</a>
           </td>
        </tr>
      </table>
    </td>
  </tr>
  <tr>
    <td align='center'>
      <table width='730' align='center' border='0' cellpadding='0' cellspacing='0' summary='statistics'>
        <tr>
          <td>
            <center><br>";
    return $chaine;
}

sub print_mise_en_page_fin() {

    $chaine = '';
    $chaine.= '</td>';
    $chaine.= '</tr>';
    $chaine.= "</table>\n";
    $chaine.= '</td>';
    $chaine.= "<td><img src='$web_base_images/px_trans.gif' width='1' height='1' alt='background'></td>\n";
    $chaine.= "<td width='1' bgcolor='#000000'><img src='$web_base_images/px_trans.gif' width='1' height='1' alt='background'></td>\n";
    $chaine.= "<td width='20' background='$web_base_images/patterndroite.gif'>&nbsp;</td>\n";
    $chaine.= '</tr>';
    $chaine.= "</table>\n";
    return $chaine;
}

sub sort_numer {  $a <=> $b;}

