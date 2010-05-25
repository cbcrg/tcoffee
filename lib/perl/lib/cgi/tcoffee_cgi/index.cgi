#!/usr/bin/env perl

######################################################################
#2006-06-26                                                          #
#Olivier POIROT, poirot [AT] igs.cnrs-mrs.fr, Labo IGS, UPR2589, CNRS#
#Sebastien MORETTI, moretti.sebastien [AT] gmail.com, SIB            #
#Web interface for Tcoffee and other programs                        #
#Need feedback.cgi, configuration_file.txt                           #
######################################################################

use CGI ':standard';
use CGI qw/:any/ ;
use CGI::Carp qw(fatalsToBrowser);
use locale;
use Time::localtime;
use Time::Local;
use GD;
use Mail::Send;
use Sys::Hostname;


#############################################################################################################################################
use env_param;                                                                                                                              #
my $location = &hostname();                                                                                                                   #
&env_param::variables($location);                                                                                                           #
#############################################################################################################################################

$CGI::POST_MAX        = 1024 * 10000;  # max 100K posts
$CGI::DISABLE_UPLOADS = 0;  # set to a non-zero value, this will disable file uploads completely.

$| = 1;

my $tcoffee_homepage  = $env_param::tcoffee_homepage;
##############  Configuration  ##############
my $config_infile     = $env_param::config_infile;
my $pg_source         = $env_param::pg_source;
my $tmpOldFiles       = $env_param::tmpOldFiles;
my $dir_pdb           = $env_param::dir_pdb;
my $webmaster         = $env_param::webmaster;
my $fromEMail         = $env_param::fromEMail;
my $ppath             = $env_param::ppath;

###############  Directories  ###############
my $dir_tcoffee       = $env_param::dir_tcoffee;
my $dir_tmp           = $env_param::dir_tmp;
my $dir_exe           = $env_param::dir_exe;
my $dir_images        = $env_param::dir_images;
my $programme         = $env_param::programme;
my $specific_lib      = $env_param::specific_lib;

############  Server directories ############
my $home_web          = $env_param::home_web;
my $home_html         = $env_param::home_html;
my $home_cgi          = $env_param::home_cgi;
my $scratch_area      = $env_param::scratch_area;
my $scratch_area2     = $env_param::scratch_area2;

#############  Web directories #############
my $web_base          = $env_param::web_base;
my $web_cgi_base      = $env_param::web_cgi_base;
my $web_base_images   = $env_param::web_base_images;
my $web_base_doc      = $env_param::web_base_doc;
my $web_base_link_tmp = $env_param::web_base_link_tmp;

##################  Look  ##################
my $colorTitle        = $env_param::colorTitle;
my $colorAppli        = $env_param::colorAppli;
my $Rgb               = $env_param::Rgb;
my $rGb               = $env_param::rGb;
my $rgB               = $env_param::rgB;
my $logo              = $env_param::logo;


my $nbr = $colorAppli;
$nbr    =~ s/^.//;
my $goodColor = '#';
SET_COLORS:
for (my $u=0; $u<length($nbr); $u+=2){
    my $valueColor = substr($nbr,$u,2);
    my $hexa       = hex($valueColor)-20;
    $goodColor    .= sprintf "%x", $hexa;
}


my @tab_list_input_files;

my $q = new CGI;

######################################################################
#Choice form #########################################################
######################################################################
#upload/paste fichier de config
if ( defined($q->param('load_cfg_file')) ){
    print $q->header;
    print &std_HTML_header();
    print &print_mise_en_page_deb();
    &load_cfg_file($q);
    print $q->end_html;
}
elsif ( defined($q->param('stage2')) ){
#process form server
    &process_form($q);
}
elsif ( defined($q->param('stage1')) ){
#print selected server
    print $q->header('text/html', -cache-control=>'NO-CACHE', -expires=>'-1');
    my $scriptjs="<SCRIPT language=\"JavaScript\" type=\"text/javascript\">
    var submitcount=0;

    function verifForm(formulaire) {
        if (submitcount == 0){
            submitcount++;
            location.replace();
            return true;
        }
        else {
            alert('already submitted');
            return false;
        }

        if(formulaire.email.value == \"\" && formulaire.mand_email.value == \"mandatory\") { /* on detecte si email est vide */
            alert('Email not valid !!'); /* dans ce cas on lance un message d'alerte */
            return false;
         }
         else {
            formulaire.submit(); /* sinon on envoie le formulaire */
            return false;
        }
    }
    </SCRIPT>";
    $mo = $q->param('daction');
    print &std_HTML_header($scriptjs);
    print &print_mise_en_page_deb();
    if( defined($q->param('config_infile')) ){
        @x = $q->param('config_infile');
        $config_infile = join(', ', @x);
    }
    &print_form_2($q);
    print $q->end_html;
}
else {
#print server list
    &reorder_results($scratch_area2);
    if ( defined $q->param('personal_server_file') ){
        $file = "$scratch_area/".'config_perso';
        if ( $q->param('pasted_cfg_file') ){
            $txt = $q->param('pasted_cfg_file');
            open( my $TXT_FILE, '>', "$file");
            print {$TXT_FILE} $txt;
            close $TXT_FILE;
        }
        else {
            $fic = $q->param('uploaded_cfg_file');
            $fh  = $q->upload($fic);
            if ( !$fh && $q->cgi_error ){
                print $q->header(-status=>$q->cgi_error);
            }
            open ( my $OUTFILE, '>', "$file");
            while ( $bytesread = read($fic,$buffer,1024) ){
                $buffer =~ s///gs;
                $buffer =~ s/\r/\n/gs;
                if ( $buffer =~ /^\W+$/ ){
                     $buffer =~ s/\*/X/gs;
                }
                print {$OUTFILE} $buffer;
            }
            close $OUTFILE;
        }
        $config_infile = $file;
    }
    else {
        $config_infile = "$home_cgi/$config_infile";
    }
    #Reads list of servers
    open ( my $CONFIG_FILE, '<', "$config_infile") or die "unable to find config_file\n";
    while (<$CONFIG_FILE>) {
        if ( $_ =~ /^server::(.*)/ || $_ =~ /^(section::.*)/ || $_ =~ /^(mirrors::.*)/ ){
            $str = $1;
            chomp $str;
            push @servers, $str;
        }
    }
    close $CONFIG_FILE;
    print $q->header(-type=>'text/html');
    print &std_HTML_header();
    print &print_mise_en_page_deb();
    &print_form_1($q,$config_infile);
    print &print_sponsor();
    print $q->end_html;
}

######################################################################
#Load Config_file#####################################################
######################################################################
sub load_cfg_file() {

    my ($q) = @_;
    print '              '.$q->start_multipart_form(-action=>"$web_cgi_base/$pg_source");
    print "              <center>
                <H3>Config-server file</H3><HR noshade size=1 width='95%' align='center'>
                <table bgcolor=\"#FFFFFF\" border=\"0\" summary=\"uploadcfg\" width=\"100%\">
                  <tr>
                    <td width=\"30%\">Upload your config-server file</td>
                    <td align=\"left\" width=\"70%\">
                      <input type=\"file\" name=\"uploaded_cfg_file\" value=\"Input your config_server file\" size=\"50\" maxlength=\"80\">
                    </td>
                  </tr>
                </table>
                <table bgcolor=\"#FFFFFF\" border=\"0\" summary=\"pastecfg\" width=\"100%\">
                  <tr>
                    <td width=\"30%\">or paste it</td>
                    <td align=\"left\" width=\"70%\">
                      <textarea name=\"pasted_cfg_file\" rows=\"5\" cols=\"60\"></textarea>
                    </td>
                  </tr>
                  <tr>
                    <td>
                      <input type=\"submit\" name=\"Submit\" value=\"Submit\">
                      <input type=\"reset\" name=\".reset\">
                    </td>
                  </tr>
                </table>
                <input type=\"hidden\" name=\"pid\" value=\"$$\">
                <input type=\"hidden\" name=\"personal_server_file\" value=\"yes\">
              </center>
              </form>

            <hr noshade size='1' width='95%' align='center'>
            <div align=center>You can see <a href='$web_base_doc/$config_infile'>here</a> the default config-server file.</div>
          </td>
        </tr>
      </table>
    </td>
  </tr>
</table>\n";
    return;
}

######################################################################
#print_mise_en_page_deb###############################################
######################################################################
sub print_mise_en_page_deb() {

    my $chaine = "\n<table summary=\"menu\" align=\"center\" bgcolor=\"#ffffff\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"792\" style=\"border-style:solid;border-width:1px;border-color:#000000;border-top-style:none;\">
  <tr valign=\"top\">
    <td width=\"750\" align=\"center\">
      <table summary=\"logo\" align=\"center\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"792\">
        <tr>
          <td colspan='2' align='center' nowrap='nowrap'>$logo</td>
        </tr>
        <tr bgcolor='#000000'>
           <td onclick=\"window.open('evol2.cgi', '_self')\">.</td>
           <td align='right' height='6' valign='top' style='color: #ffffff;' width='791'>
             <a href=\"$web_cgi_base/$pg_source\" style=\"color: #ffffff;\">HOME</a>&nbsp;|&nbsp;
             <a href=\"$tcoffee_homepage\" style=\"color: #ffffff;\">references</a>&nbsp;|&nbsp;
             <a href=\"$web_base_doc/doc3.html\" target=\"_blank\" style=\"color: #ffffff;\">help</a>&nbsp;|&nbsp;\n";

    $chaine .= '             '.&env_param::transpose_eMails_delayed($webmaster, "$web_base_images/mail1.gif")."&nbsp;&nbsp;&nbsp;
             <a href=\"$web_cgi_base/${pg_source}?load_cfg_file=1\" style=\"text-decoration: none; color: #000000;\">.</a>
           </td>
        </tr>
      </table>
      <table width=730 align=center border=0 cellpadding=0 cellspacing=0 summary=\"appli\">
        <tr>
          <td>
            \n";

  return $chaine;
}

######################################################################
#print sponsor########################################################
######################################################################
sub print_sponsor() {

    my $chaine = "            <table cellpadding=0 cellspacing=0 border=0 bgcolor=\"#C0C0C0\" summary=\"sponsor\">
              <tr>
                <td>
                  <table cellpadding=3 cellspacing=1 border=0 summary=\"book\">
                    <tr>
                      <td bgcolor=\"$colorTitle\" align=\"center\" style=\"font-family:arial;color:#333333;font-size:small\">
                        Featured in
                      </td>
                    </tr>
                    <tr>
                      <td bgcolor=\"$colorAppli\">
                        <a HREF=\"http://www.amazon.com/gp/product/0470089857?ie=UTF8&amp;tag=wwwtcoffeeorg-20&amp;linkCode=as2&amp;camp=1789&amp;creative=9325&amp;creativeASIN=0470089857\" target=\"_blank\">
                        <img border=0 src=\"$web_base_images/cover2.jpg\" alt=\"Bioinformatics for dummies\"></a>
                      </td>
                    </tr>
                  </table>
                </td>
              </tr>
            </table>
          </td>
          <td width=\"5%\">&nbsp;&nbsp;</td>
          <td align=\"left\" width=\"40%\">
            <table bgcolor=\"#ffffff\" border=0 summary=\"homepage\">
              <tr>
                <td style=\"font-family:arial;font-size:x-small\" align=\"left\">
                  All these packages are opensource freeware. Download them and look for extra documentation on
                  <A HREF=\"$tcoffee_homepage\">The Tcoffee Homepage</A>
                  and on my personal <A HREF=\"http://www.tcoffee.org/homepage.html\">Homepage</A>.
                  <br><br>
                  The Book <B><A href=\"http://www.amazon.com/gp/product/0470089857?ie=UTF8&amp;tag=wwwtcoffeeorg-20&amp;linkCode=as2&amp;camp=1789&amp;creative=9325&amp;creativeASIN=0470089857\" target=\"_blank\">Bioinformatics For Dummies</A></B> contains a T-Coffee tutorial and many other <A href=\"http://www.dummies.com/WileyCDA/DummiesTitle/productCd-0470089857,page-1.html\" target=\"_blank\"><B>online resources</B></A>.
                  <br><br>
                  The following <A href=\"http://www.tcoffee.org/Books/\">Text-Books</A> can also help you interpret your results.
                </td>
              </tr>
            </table>
          </td>
          <td width=\"15%\"></td>
        </tr>
      </table><br>
    </td>
  </tr>
</table><br>\n";
    return $chaine;
}

######################################################################
#form init ###########################################################
######################################################################
sub print_form_1 {

    my ( $q, $config_infile ) = @_;
    #Liste des serveurs dispo#############################################

    print '            '.$q->start_form(-action=>"$web_cgi_base/$pg_source");
    print "            <table bgcolor='$colorTitle' border='1' align='center' rules='all' cellpadding='3' summary='form'>
              <tr>
                <td align='center'>
                  <a href='$tcoffee_homepage'>
                  <img src='$web_base_images/t-coffee4_big.png' border='0' alt='T-Coffee'></a><br>
                  <b><font color='#333333' face='arial' size='4'>A collection of tools for Computing, Evaluating and Manipulating Multiple Alignments of DNA, RNA, Protein Sequences and Structures</font></b>
                </td>
              </tr>
            </table>\n"; #CDCEFF#CBA574#5dbdba=vert1
    print '            ',$q->hidden('stage1','1'),"
            <table align='center' bgcolor='$colorAppli' border='2' cellspacing='0' cellpadding='0' width='70%' summary='appli form' style='border-style:none;'>\n"; #CDCEFF#0080C0#CEBCA4#92d7d5=vert2
    #print "<center><table bgcolor=#e7c7a8 border=1 cellspacing=1 cellpadding=1 rules=all width=60% >";#CDCEFF#0080C0#CEBCA4#92d7d5=vert2
    my $compteur = -1;
    SERVERS:
    foreach $s (@servers){
        $compteur++;
        if ( $s =~ /^mirrors::(.+)$/ ){
            my $mirrors_country = $1;
               $mirrors_country =~ s{\r}{}g;
            my @mirrors_lab     = split(/::/, $mirrors_country);
            my $mirrors_place   = '';
            MIRRORS:
            for (my $i=0; $i<=$#mirrors_lab; $i++){
                if ( $mirrors_lab[$i] =~ /(IGS)/ ){
                    $mirrors_place .= "&nbsp;<a href='".&env_param::goodFlag("$1")."' title='IGS lab., CNRS'><img src='${web_base_images}/cnrs.jpg' border='0' align='middle' alt='IGS lab., CNRS'></a>&nbsp;";
                }
                elsif ( $mirrors_lab[$i] =~ /(HOME)/i ){
                    $mirrors_place .= "&nbsp;<a href='".&env_param::goodFlag("$1")."' title='www.tcoffee.org'><img src='${web_base_images}/t-coffee4.png' align='middle' border='0' alt='www.tcoffee.org'></a>&nbsp;";
                }
                elsif ( $mirrors_lab[$i] =~ /(Vital-IT)/i ){
                    $mirrors_place .= "&nbsp;<a href='".&env_param::goodFlag("$1")."' title='Vital-IT, Swiss Institute of Bioinformatics (SIB)'><img src='${web_base_images}/logo_sib.jpg' align='middle' border='0' alt='Vital-IT, Swiss Institute of Bioinformatics (SIB)'></a>&nbsp;";
                }
                elsif ( $mirrors_lab[$i] =~ /(EBI)/ ){
                    $mirrors_place .= "&nbsp;<a href='".&env_param::goodFlag("$1")."' title='EBI, European Bioinformatics Institute'><img src='${web_base_images}/EBI.gif' align='middle' border='0' alt='EBI, European Bioinformatics Institute'></a>&nbsp;";
                }
                else{
                    #USE: ::mirror_name#mirror_address#mirror_logo
                    #mirror logo can either be a we address (http://...) or a file in ${web_base_image}
                    #BETTER if images are 36x25 pixels to avoid reformating by browsers and image pixelisation

                    my ( $mirror_name, $mirror_address, $mirror_logo ) = split (/\#/, $mirrors_lab[$i]);
                    $mirror_logo    = "${web_base_images}/$mirror_logo" if ( $mirror_logo !~ /^(ht|f)tp/i );

                    $mirrors_place .= "&nbsp;<a href='$mirror_address' title='$mirror_name'><img src='$mirror_logo' align='middle' border='0' alt='$mirror_name'></a>&nbsp;";
                }
            }
            print "\n              <tr>
                <th align='center' colspan='5' bgcolor='#ffffff' height='30' style='border-style:none;vertical-align: middle;'>
                  <em style='vertical-align: middle;'>Mirror sites:</em> $mirrors_place
                </th>
              </tr>\n";
            next SERVERS;
        }
        if ( $s =~ /^section::(.+)$/ ){
            print "\n              <tr>
                <th align=\"center\" colspan=5 bgcolor=\"$goodColor\" height=30 valign=\"middle\" style=\"border-style:none\">
                  <B>$1</B>
                </th>
              </tr>\n";
            next SERVERS;
        }
        $s     = trim($s);
        ( $mode, $level, $desc, $PMID ) = split(/::/, $s);
        $level = trim($level);
        $mode  = trim($mode);
        $publi = '';
        $publi = "<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=PubMed&amp;list_uids=$PMID\" target=\"_blank\">cite</a>" if($PMID);
        chomp($desc);
        my $jsdesc = $desc;
        $jsdesc =~ s/ /%20/g;
        if ( $tab_ass_mode{$mode}++ == 0 ){
            print "\n              <tr align=center>
                <td width=\"40%\" height=45 align=\"center\">\n";
            $im         = new GD::Image(220,35);
            $vertpale   = $im->colorAllocate(146,215,213);
            $rougepale  = $im->colorAllocate(251,101,111);
            $marron     = $im->colorAllocate(231,199,168);
            $palef      = $im->colorAllocate(206,188,164);
            $white      = $im->colorAllocate(255,255,255);
            $black      = $im->colorAllocate(0,0,0);
            $red        = $im->colorAllocate(255,0,0);
            $blue       = $im->colorAllocate(0,0,255);
            my $myColor = $im->colorAllocate($Rgb,$rGb,$rgB);
            $im->filledRectangle(0,0,220,35,$myColor);
            #my $courier = GD::Font->load('./courierR12.fnt') or die "Can't load font";
            $im->string(gdGiantFont,8,10,"$mode",$blue);
            binmode STDOUT;
            $look2 = $im->png;
            open ( my $OUT, '>', "${home_html}/$dir_tcoffee/$dir_images/$n.lookbis.png") or die "ERROR writing to file: $!\n";
            print {$OUT} "$look2";
            close $OUT;
            print "                  <A href=\"javascript:alert('$jsdesc')\" title=\"$desc\">
                  <img src=\"${web_base_images}/$n.lookbis.png\" border=0 alt=\"program\"></A>
                </td>\n";
            $n++;
            print "                <td align=\"center\">\n";

        }
        else {
            print "                <td align=center>\n";
        }

        if (`diff $config_infile $home_cgi/configuration_file.txt`){
            $strr = $web_cgi_base."/".$pg_source."?stage1=1&daction=".$mode."::".$level."&config_infile=".$config_infile;
        }
        else{
            $strr = $web_cgi_base."/".$pg_source."?stage1=1&daction=".$mode."::".$level;
        }
        my $referer = $ENV{HTTP_REFERER};
        if ( $referer =~ /myhits/i ){
            $referer = 'myhits';
        }
        elsif ( $referer =~ /ch.embnet/i ){
            $referer = 'embnet';
        }
        else {
            $referer = '';
        }
        $strr .= '&referer0='.$referer if ( $referer ne '' );
        print '                  '.button(-name   => 'daction',
                                          -value  => $level,
                                          -onClick=> "parent.location=\"$strr\"",
        );
        print "\n                </td>\n";
        print "                <td align=\"center\">
                  $publi
                </td>
                <TD align=\"center\">
                  &nbsp;<A href=\"javascript:alert('$jsdesc')\" title=\"$desc\">?</A>&nbsp;
                </TD>
              </tr>\n" if ($tab_ass_mode{$mode} == 2);
        print "                <TD align=\"center\">&nbsp;
                </TD>
                <td align=\"center\">
                  $publi
                </td>
                <TD align=\"center\">
                  &nbsp;<A href=\"javascript:alert('$jsdesc')\" title=\"$desc\">?</A>&nbsp;
                </TD>
              </tr>\n" if ($tab_ass_mode{$mode} == 1 and $servers[$compteur+1] !~ /Advanced/);
    }

    print "            </table>
            </form>

          </td>
        </tr>
      </table>
      <table summary=\"foot page\" align=\"center\" width='730' border='0' cellpadding='0' cellspacing='0'>
        <tr>
          <td align=\"right\">\n";
    return;
}

######################################################################
#Form server #########################################################
######################################################################
sub print_form_2 {

    my ($q) = @_;
    if ( defined $q->param('daction') ){
        $daction = $q->param('daction');
    }
    else {
        $daction = 'Perso_server';
    }
    print '            '.$q->start_multipart_form(-action=>"$web_cgi_base/$pg_source",-name=>"my_form",-onSubmit=>"return verifForm(my_form)");
    my $referer = $q->param('referer0');
    print "              <center>
              <input type=\"hidden\" name=\"stage2\" value=\"2\">
              <input type=\"hidden\" name=\"daction\" value=\"$daction\">
              <input type=\"hidden\" name=\"pid\" value=\"$$\">
              <input type='hidden' name='referer' value='", $referer, "'>
";
    my ( $code_htmll, $sendto, $maxnseq, $maxlseq ) = &print_html_form($q, $daction);
    $code_htmll =~ s/^ +//;
    $code_htmll =~ s{\r}{}g;
    print "\n              ".$q->hidden('code_ht', "$code_htmll");
    print "\n              <input type=\"hidden\" name=\"send_to\" value=\"$sendto\">
              </center>
            </form>

          </td>
        </tr>
      </table>
    </td>
  </tr>
</table><br>\n";
    return;
}

######################################################################
#Process Form ########################################################
######################################################################
sub process_form {

    my ($q) = @_;
    $my_pid        = $q->param('pid');
    $daction       = $q->param('daction');
    $code_htm      = $q->param('code_ht');
    my $sendto     = $q->param('send_to');
    ($mode,$level) = split(/::/,$daction);
    $rand_number   = int(rand 100000)+1;
    $run_name      = 'tcf'.substr($mode,0,3).substr($level,0,1).$rand_number.'_'.$my_pid;
    my $appli          = substr($mode,0,3).substr($level,0,1);
    $scratch_area      = "$scratch_area/$appli";
    $web_base_link_tmp = "$web_base_link_tmp/$appli";
    #@list_visible_outfile = split /,/,$q->param('list_visible_outfile');
    $list_visible_outfile = $q->param('list_visible_outfile');
    $cptfile = 0;
    PARAM:
    foreach $key ($q->param){
        @values = $q->param($key);
        if ( $key eq '-infile' || $key eq '-in' || $key eq '-aln' || $key eq '--msa' || $key eq '-template_file' || $key eq '--template' ){
            PARAM_VAL:
            foreach $txt (@values){
                if ($txt){
                    if ( $txt !~ /^M.*_pair/ && $txt !~ /M.*_msa/ && $txt !~ /SCRIPT_/ && $txt !~ /SELF_/ && $txt !~/MODE_/ ){
                        $list_param .= " $key ";
                        $fil         = $run_name.'.in'.$cptfile;
                        $file        = "$scratch_area/".$fil;
                        open( my $TXT_FILE, '>', "$file");
                        #sous Unix, in bash/tcsh, press Ctrl-V then Ctrl-M , ca fait \r \n
                        $txt =~ s/\r\n/\n/gs;
                        $txt =~ s/\r/\n/gs;
                        #p.206 perl en action
                        if ( $txt =~ /^>/ ){
                            #fasta file with * make tcoffee crash!
                            $txt =~ s/\*/X/gs;
                            #$txt =~ s/\//_/gs;
                            #$txt =~ s/\|/_/gs;A verifier
                        }
                        print {$TXT_FILE} $txt."\n";
                        close $TXT_FILE;
                        if ( -B "$file" ){
                            print $q->header();
                            print $q->start_html(-title=>"$mode :: $level", -BGCOLOR=>'#F6F6FF');
                            print 'You have submited a binary file, certainly a Microsoft Office Word Document. You must only submit an ascii - flat file format - file';
                            print $q->end_html;
                            exit(1);
                        }
                        $cptfile++;
                        $list_param .= " $file ";
                        push @tab_list_input_files, $fil;
                    }
                    else {
                        $list_param .= " $key $txt ";
                    }
                }
            }
        }
        elsif ( $key =~ /executable/ ){
            $programme = " $dir_exe".'/'.$values[0];
        }
        elsif ( $key =~ /-uploaded_file/ ){
            $key2 = $key;
            $key2 =~ s/-uploaded_file//g;
            $key2 =~ s/_pdb/-pdb/g;
            foreach $fic (@values){
                if ($fic){
                    $list_param .= " $key2 ";
                    $fh = $q->upload($fic);
                    $fil = $run_name.'.in'.$cptfile;
                    $file = "$scratch_area/".$fil;
                    if ( !$fh && $q->cgi_error ){
                        print $q->header(-status=>$q->cgi_error);
                    }
                    open ( my $OUTFILE, '>>', "$file");
                    while ( $bytesread = read($fic,$buffer,1024) ){
                        #$buffer =~ s//\n/gs;
                        #sous Unix, in bash/tcsh, press Ctrl-V then Ctrl-M , ca fait \r \n
                        $buffer =~ s/\r\n/\n/gs;
                        $buffer =~ s/\r/\n/gs;
                        if ( $buffer =~ /^\W+$/ ){
                            $buffer =~ s/\*/X/gs;
                        }
                        if ( $buffer =~ /^>/ ){
                            $buffer =~ s/\*/X/gs;
                            #$buffer =~ s/\//_/gs;
                            #$buffer =~ s/\|/_/gs;
                        }
                        print {$OUTFILE} $buffer;
                    }
                    close $OUTFILE;
                    if ( -B "$file" ){
                        print $q->header();
                        print $q->start_html(-title=>"$mode :: $level", -BGCOLOR=>'#F6F6FF');
                        print 'You have submited a binary file, certainly a Microsoft Office Word Document. You must only submit an ascii - flat file format - file';
                        print $q->end_html;
                        exit(1);
                    }
                    push @tab_list_input_files, $fil;
                    if ( $key =~ /pdb/ ){
                        #$list_param .= ' P'.$run_name.'.in'.$cptfile.' ';
                        $list_param .= ' '."$scratch_area/".$run_name.'.in'.$cptfile.' ';
                    }
                    else {
                        $list_param .= " $file ";
                    }
                    $cptfile++;
                }
            }
        }
        else {
            if ( $key =~ /^-/ ){
                if ( $values[0] ){
                    $list_param .= " $key ".join(' ',@values);
                }
            }
        }
    }
    if ( $cptfile > 0 ){
        $email = $q->param('email');
        if ($email){
            if ( ! &ValidEmailAddr($email) ){
                print $q->header();
                print $q->start_html(-title=>"$mode :: $level", -BGCOLOR=>'#FFF0F5');
                print $q->h4('Error: Email not Valid...');
                print $q->end_html;
                exit(1);
            }
        }

        #/home/igs/Tools/T-COFFEE_web/Struct contient le fugue_client qui fonctionne en local
        #Appelle fugue_align.sh qui utilise magicPDB.pl pour ne garder que la chaine A des fichiers pdb
        #export PDB_DIR=/usr/local/httpd/htdocs/Tcoffee/Tmp;  A rajouter au debut de $env_var si recup pdbFile via net
        #"si /Sequences out"
        #HOME necessaire pour .tmp
        $env_var = "export NO_ERROR_REPORT_4_TCOFFEE=1;export HOME=$scratch_area;export PDB_DIR=$dir_pdb;export PATH=$dir_exe:$ppath:$dir_exe;\n";

        # si pdb cherche sur le net
        #$env_var = "export NO_REMOTE_PDB_DIR=0;export PATH=/home/igs/public_html/Tcoffee/Exe:/usr/local/bin/:\$PATH;\n";
        #$env_var = "export NO_REMOTE_PDB_DIR=0;export PATH=/home/vital-it/bnyffele/projects/tcoffee/bin:\$PATH;\n";

        $run_name_param = " -run_name=$scratch_area/$run_name";
        $opt_check_pdb  = '';
        if ( $mode =~ /3D/ ){
            $opt_check_pdb = ' -check_pdb_status ';
        }
        $default_option = "  $opt_check_pdb -cache=no -remove_template_file=1 -quiet=stdout >$scratch_area/$run_name.tc_log 2>$scratch_area/$run_name.tc_log2 ";
        if ( $programme =~ /ProtoGene/ ){
            $default_option = " >$scratch_area/$run_name.tc_log 2>$scratch_area/$run_name.tc_log2 ";
        }
        if ( $mode =~ /apdb/i ){
            $default_option = "  $opt_check_pdb  -quiet=stdout >$scratch_area/$run_name.tc_log 2>$scratch_area/$run_name.tc_log2 ";
        }

        my $who = $ENV{HTTP_X_FORWARDED_FOR}; #for through a proxy
        $who    = $ENV{REMOTE_ADDR} if ( $who eq '' );
        open( my $WHO_FILE, '>', "$scratch_area/$run_name.tc_log_who");
        print {$WHO_FILE} $who;
        close $WHO_FILE;

        my $referer = $q->param('referer');
        $referer    = $ENV{HTTP_REFERER} if ( $referer eq '' );
        if ( $referer ne '' && $referer =~ /tcoffee/i ){
            $referer = '';
        }

        #test IP submision
        #if( -e "$scratch_area/$who") {
        if ( 1==2 ){
            #IP already exists, please submit later
            print $q->header();
            print $q->start_html(-title=>"$mode :: $level", -BGCOLOR=>"#F6F6FF");
            print "You have already submited a job, please wait untill its completion to submit another one...<br> $scratch_area/$who<br>Thank you for your comprehension";
            print $q->end_html;
            exit(1);
        } #IP doesn't exists, ok to run
        else {
            @t_who = split(/ /,$who);
            $f_who = $t_who[0];
            open( my $WHO_FILE, '>', "$scratch_area/$f_who");
            print {$WHO_FILE} $who;
            close $WHO_FILE;
        }

        $command = "cd $scratch_area;";
        $command .= "$env_var /usr/bin/env time -p -o $scratch_area/$run_name.titime $programme  $list_param $run_name_param $default_option";

        if ( $daction =~ /Combine/ ){
            $command =~ s/infile/in/g;
        }
        #$pfam_command = " && \nat now <<EOF\nexec /home/igs/Tools/ESPript_associated_tools/ESPript_display2.sh $scratch_area $run_name pfam 2>/dev/null\nEOF\n;";
        #$command .= $pfam_command;
        $command =~ s/\n/ /gs;
        $command =~ s/\r/ /gs;

        my $jobname_validity = 0;
        CHECK_JOBNAME:
        for my $srv (&env_param::server_list($config_infile)){
            $jobname_validity++ if ( $run_name =~ /tcf$srv\d+_\d+/ );
            last CHECK_JOBNAME  if ( $jobname_validity==1 );
        }
        die "Invalid server and/or job name\n" if ( $jobname_validity==0 );
        undef $jobname_validity;
        unless ( $f = fork ){
            # child process - close input and output  # pipes so parent isn't left hanging
            close (STDERR);
            close (STDOUT);
            # run the background job

            &tmpManagement($scratch_area);

            if ( $location !~ /^frt/ ){ #CNRS IGS Marseille & other
                system($command);
            }
            elsif ( $location =~ /^frt/ ){ #SIB Vital-IT
                my $bsub="\#\!\/bin\/bash\n\n\#BSUB -L \/bin\/bash\n\#BSUB -o $scratch_area/$run_name.tc_log\n\#BSUB -e $scratch_area/$run_name.tc_log2\n";
                #$bsub.="\#BSUB -m Xeon\n" if ($run_name =~ /^tcfMCO/);
                $bsub   .= "\#BSUB -J $run_name\n";

                my $cmdLine = $command;
                $command    =~ s/cd $scratch_area;//;
                $command    =~ s/\/usr\/bin\/env time -p -o $scratch_area\/$run_name.titime //; #No need to get duration because it is in bsub log
                $command    =~ s/ >$scratch_area\/$run_name.tc_log//;
                $command    =~ s/ 2>$scratch_area\/$run_name.tc_log2//;
                $command    =~ s/ *\; */\n/g;
                $command    = $bsub."\n".$command;

                my $ScratchTmp = $scratch_area2.'/'.$appli;
                $command =~ s/$scratch_area/$ScratchTmp/g; #Because frontal see $scratch_area but nodes can only see $ScratchTmp

                `echo \"$command\" >$scratch_area/$run_name.wwww`;
                my $bjob = `. /mnt/common/lsf/conf/profile.lsf; echo \"$command\" | bsub -q normal 2>>$scratch_area/$run_name.wwww`; # On nodes of a computers cluster
                chomp($bjob);
                $bjob =~ s/^Job \<(\d+)\> is submitted to queue.*$/$1/;
                #`echo -e \"\n[Job $bjob is submitted to queue ...]\n\" >>$scratch_area/$run_name.wwww`;

                my $Status=1;
                CHECK_JOB_STATUS:
                while ( $Status =! 0 ){
                    my $ProcessStatus = `. /mnt/common/lsf/conf/profile.lsf; bjobs -p $bjob |grep -v 'JOBID'`;
                    chomp($ProcessStatus);
                    if ( $ProcessStatus =~ /DONE/ || $ProcessStatus =~ /EXIT/ ){ #Job has terminated successful: DONE   --   Job has failed: EXIT
                        $Status=0;
                        last CHECK_JOB_STATUS;
                    }
                    if ( $ProcessStatus =~ /ZOMBI/ ){
                        system(". /mnt/common/lsf/conf/profile.lsf; bkill $bjob");
                        $Status=0;
                    }
                    sleep 3;
                }
                $command = $cmdLine;
                $command =~ s/^.+\.titime //;
            }

            unlink("$scratch_area/$f_who");
            &tc_log_analyse($run_name,$command,$daction,$code_htm,$sendto);
            #system("mv $scratch_area/${run_name}* ${home_html}/${dir_tmp}/"); #Add this move if there is a html download quota between apache and SFS volume

            # email the user here once the job is done
            if ( defined ($email) ){
                my $msg = new Mail::Send Subject=>"$mode :: $level job done", To=>"$email";
                $msg->add('From', $fromEMail);
                $msg->add('Reply-To', $webmaster);
                $msg->add('X-Mailer', "Mail::Send ($Mail::Send::VERSION) Perl module");

                my $fh = $msg->open;
                print $fh "$mode :: $level\nYour job is done, you can consult your results at this URL : ${web_base_link_tmp}/${run_name}.file_result.html\nYour data will remain available on this server over the next $tmpOldFiles days. It will then be deleted.";
                $fh->close;
            }

            #system ("echo \"$mode :: $level\nYour job is done, you can consult your results at this URL : $web_base_link_tmp/$run_name.file_result.html\nYour data will remain available on this server over the next $tmpOldFiles days. It will then be deleted.\" 2>&1 | /usr/bin/env Mail  $email -s \"$mode :: $level job done\" ") if (defined ($email)); # -R $webmaster option doesn't exist on frt Mail program
        }
        else {
            # parent process -- show output page and exit

            $tm = localtime;
            ( $SEC, $MIN, $HOUR, $DAY, $MONTH, $YEAR ) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, $tm->mon, $tm->year);
            $time_init = timelocal($SEC,$MIN,$HOUR,$DAY,$MONTH,$YEAR);
            $SEC       = '0'.$SEC   if ( length($SEC) < 2 );
            $MIN       = '0'.$MIN   if ( length($MIN) < 2 );
            $HOUR      = '0'.$HOUR  if ( length($HOUR) < 2 );
            $DAY       = '0'.$DAY   if ( length($DAY) < 2 );
            $MONTH     = $MONTH+1;
            $MONTH     = '0'.$MONTH if ( length($MONTH) < 2 );
            #$uurl = '';
            $level = trim($level);
            print $q->header(-Refresh=>"2; URL=feedback.cgi?level=$level&child=$f&result_file=$run_name.file_result.html&pg_source=$0&time_init=$time_init&email=$email&mode=$mode");
            print $q->start_html(-title=>"$mode :: $level", -BGCOLOR=>"#F6F6FF");
            print '<div>Processing, please wait...</div>';
            print "\n<div>E-mail will be sent to $email when the job completes.</div>\n" if ($email);
            print $q->end_html;


            open( my $FH, '>>', "$home_cgi/userlog.xml");
            print {$FH} "<user-call>\n <file>",$run_name,"</file>\n <remote-host>",$who,"</remote-host>\n <from>", $referer, "</from>\n <user-agent>",$q->user_agent(),"</user-agent>\n <mail>",$email,"</mail>\n <time-stamp>",$HOUR,":",$MIN,":",$SEC," ",$DAY,"/",$MONTH,"/",$YEAR+1900,"</time-stamp>\n</user-call>\n";
            close $FH;


            exit (0);
        }
    }
    else {
        print $q->header();
        print $q->start_html(-title=>"$mode :: $level", -BGCOLOR=>"#FFF0F5");
        print $q->h4("Error: You didn't provide any data...");
        print $q->end_html;
        exit(1);
    }
    return;
}

######################################################################
#print tag/form ######################################################
######################################################################
sub print_elem_form {

    my ( $q, $param, $comment, $tag, $list_values, $default_values, $docu ) = @_;
#    $ref = $q->self_url."#$param";
    @val     = split(/,/, $list_values);
    @default = split(/,/, $default_values);
    $docu    =~ /doc_[a-z]*:(.*)$/;
    $docu    = $1;
    $docu    =~ s{\r}{}g;
    chomp($docu);
    my $jsdocu = $docu;
    $jsdocu    =~ s/ /%20/g;
    $doc = "<a href=\"javascript:alert('$jsdocu')\" title='$docu'>$comment</a>";
    #$doc = "<font title=\"$docu\" color=blue>$comment</font>";

    SWITCH: {
    ###################################################
    #checkbox_group ###################################
    ###################################################
        if ( $tag eq 'checkbox_group' ){
            my $rows = 3;
            $rows = $#val+1 if ( ($#val+1) < 3 );
            print "                    <table bgcolor='#FFFFFF' border='0' summary='Checkbox' width='100%'>
                      <tr>
                        <td width='30%'>
                          $doc
                        </td>
                        <td align='left' width='70%'>
                          ", $q->checkbox_group(-name=>$param, -values=>[@val], -default=>[@default], -rows=>$rows, -columns=>4),"
                        </td>
                      </tr>
                    </table>\n";
            last SWITCH;
        }
    ###################################################
    #scrolling_list ###################################
    ###################################################
        if ( $tag eq 'scrolling_list' ){
            #vestige de t_coffee -p, utilise par aucun parametre pour l'instant
            if (1 == 1){
                $mul = 'true';
            }
            else {
                $mul = 'false';
            }
            print $q->table({-border=>0,-width=>'100%',-bgcolor=>'#FFFFFF',-summary=>'Scroll'},$q->Tr($q->td({-width=>'20%'},$doc),$q->td({-width=>'70%',-align=>'left'},$q->scrolling_list(-name=>$param,-values=>[@val],-default=>[@default],-multiple=>$mul))));
            last SWITCH;
        }
    ###################################################
    #upload ###########################################
    ###################################################
        if ( $tag eq 'upload' ){
            my $example = '';
            $example    = "<a href='$web_base_doc/$list_values' target='_blank'>Example</a>" if ( $list_values ne 'empty' );
            $param2 = '-uploaded_file'.$param;
            print "                    <table bgcolor='#FFFFFF' border='0' summary='UploadSeq' width='100%'>
                      <tr>
                        <td width='30%' style='text-align:justify'>
                          $doc &nbsp; $example
                        </td>
                        <td align='left' width='70%'>
                          ",$q->filefield(-name=>$param2,-default=>'Input your set of sequences',-size=>'50',-maxlength=>'80'),"
                        </td>
                      </tr>
                    </table>\n";
            $docc = "<font title='$docu' color=blue>or paste data</font>";
            #$docc = "<A href=\"javascript: alert('$docu')\">or paste data</A>";
            print "                    <table bgcolor='#FFFFFF' border='0' summary='PasteSeq' width='100%'>
                      <tr>
                        <td width='30%'>
                          $docc
                        </td>
                        <td align='left' width='70%'>
                          ",$q->textarea(-name=>$param,-rows=>5,-columns=>60),"
                        </td>
                      </tr>
                    </table>\n";
            last SWITCH;
        }
    ###################################################
    #upload pdb########################################idem upload mais sans text_area
    ###################################################
        if ( $tag eq 'upload_pdb' ){
            print "                    <table bgcolor=\"#FFFFFF\" border=\"0\" summary=\"UploadPDB\" width=\"100%\">
                      <tr>
                        <td width=\"30%\">
                          $doc
                        </td>
                        <td align=\"left\" width=\"70%\">
                          ",$q->filefield(-name=>'-uploaded_file_pdb',-default=>'Input your set of sequences',-size=>50,-maxlength=>80),"
                        </td>
                      </tr>
                    </table>\n";
            last SWITCH;
        }
    ###################################################
    #popup_menu ########################################1FEFD5
    ###################################################
        if ( $tag eq 'popup_menu' ){
            print "                    <table bgcolor='#FFFFFF' border='0' summary='Popup' width='100%'>
                      <tr>
                        <td width='30%'>
                          $doc
                        </td>
                        <td align='left' width='70%'>
                          ", $q->popup_menu(-name=>$param, -values=>[@val], -default=>[@default]), "
                        </td>
                      </tr>
                    </table>\n";
            last SWITCH;
        }
        if ( $tag eq 'hidden' ){
            print '                    '.$q->hidden(-name=>$param, -default=>[@default]);
        }
        last SWITCH;

        #$nothing = 1;
    }
    return $str;
}

######################################################################
#print_paragraph######################################################
######################################################################
sub print_paragraph {

    my ( $q, $where, $str1, $str2, $maxNseq, $maxLseq ) = @_;

    if ($str2){
        if ( $str1 eq 'Description' ){
            $col  = '#F4F4FF';
            $si   = '2';
            $str2 =~ s/_maxNseq_/$maxNseq/;
            $str2 =~ s/_maxLseq_/$maxLseq/;
        }
        else {
            $col = '#FFFFCC';
            $si  = '1';
        }
        $mem_str2    = $str2;
        $flag_anchor = 0;
        while ( $str2 =~ /<doc_>(.*?)<\/doc_>/gs ){
            ( $fic_doc, $anchor ) = split(/ /, $1);
            $l           = "<a href='$web_base_doc/$fic_doc#$anchor' target='_blank'>$anchor</a>";
            $mem_str2    =~ /<doc_>(.*?)<\/doc_>/;
            $fi          = $';
            $str3       .= $`.$l;
            $mem_str2    = $';
            $flag_anchor = 1;
        }
        if ($flag_anchor){
            $str2 = $str3.$fi;
        }
        if( $where =~ /in/ ){
            print "                    <table bgcolor='$col' summary='Desc' width='100%'>
                      <tr align='left'>
                        <td style='font-family: arial; color:#333366;'>
                          <font size='$si'>$str2</font>
                        </td>
                      </tr>
                    </table>";
        }
        else {
            $code .="              <table width='100%' bgcolor='$col' summary='Desc'>
                <tr align=left>
                  <td style='font-family: arial; color:#333366; font-size:10pt'>$str2</td>
                </tr>
              </table>\n";
            return $code;
        }
    }
    else {
        if ( $where =~ /in/ ){
            my $les_arg     = "stage1=1&amp\;daction=$mode".'::Regular';
            my $other_arg   = "stage1=1&amp\;daction=$mode".'::Advanced';
            my $linkA       = "$web_cgi_base/$pg_source?$les_arg";
            my $linkB       = "$web_cgi_base/$pg_source?$other_arg";
            my $cssA        = $level eq 'Regular' ? ' style="display:none"'
                            :                       '';
            my $cssB        = $level eq 'Regular' ? ''
                            :                       ' style="display:none"';

            print "                    <br>
                    <table bgcolor=\"#99FFFF\" summary=\"Title\" width=\"100%\">
                      <tr>
                        <td width='30%' style='font-family:arial; font-size:8pt' align='left'><span$cssA>switch to <a href='$linkA'>$mode :: Regular</a>&nbsp;</span></td>
                        <td align=\"center\" style=\"font-family: arial; color:#817f82; font-size:14pt\">$str1</td>
                        <td width='30%' style='font-family:arial; font-size:8pt' align='right'><span$cssB>switch to <a href='$linkB'>$mode :: Advanced</a></span></td>
                      </tr>
                    </table>";
        }
        else {
            $code .="              <table width=\"100%\" bgcolor=\"#99FFFF\" summary=\"Title\">
                <tr align=center>
                  <td style=\"font-family: arial; color:#817f82; font-size:14pt\">$str1</td>
                </tr>
              </table>\n";
            return $code;
        }
    }
    return;
}


######################################################################
#print form/server ###################################################
######################################################################
sub print_html_form {

    my ($q) = @_;
    my $daction = $_[1];
    my $sendto  = '';
    my ( $maxNseq, $maxLseq ) = ( '', '' );

    ( $mode, $level ) = split(/::/, $daction);
    print "              <table summary='server form'>
                <tr>
                  <td align='center'>
                    <h2 style='font-family: arial;color:#3AC6C9'>$daction</h2>";

    my @messages;
    @messages = &alertMessages if ( -e 'Info.txt' );
    if ( exists($messages[0]) ){
        my $safe_mode = $mode;
        $safe_mode =~ s/\(.*$//;
        my @serverMessages = grep(/\G^$safe_mode\E/i, @messages);
        if ( exists($serverMessages[0]) ){
            my @warnings = map( m/\G^$safe_mode\:\:(.*)\E/i, @serverMessages);
            print "\n                    <h4 style=\"color:#FF0000\">",join('<br>', @warnings),"</h4>" if ( exists($warnings[0]) );
        }
    }

    $/ = 'server::';
    open ( my $CONFIG_FILE, '<', "$config_infile") or die "Unable to find config_file\n";
    <$CONFIG_FILE>;
    CONFIG:
    while (<$CONFIG_FILE>){
        next CONFIG if (/^#/);
        @conf_server = split("\n", $_);
        chop $conf_server[0];
        $serv        = shift @conf_server;
        ( $s1, $s2, $s3 ) = split(/::/, $serv);
        $serv        = $s1.'::'.$s2;
        $serv        = trim($serv);
        $daction     = trim($daction);
#       $code_html   = '';

        if ( $daction eq $serv || $daction eq 'Perso_server' ){
            foreach $spec(@conf_server){
                if ( $spec =~ /^parameter::(.*)/ ){
                    if ( $spec =~ /-maxnseq::.*::(\d+)::empty$/ || $spec =~ /--lim::.*::(\d+)::empty$/ ){
                        $maxNseq=$1;
                    }
                    if ( $spec =~ /-maxlen::.*::(\d+)::empty$/ ){
                        $maxLseq=$1;
                    }
                }
            }
        }

        if ( $daction eq $serv || $daction eq 'Perso_server' ){
            foreach $spec(@conf_server){
                print "\n";
                if ( $spec =~ /^paragraph_in::(.*)/ ){
                    @tab = split /::/, $spec;
                    &print_paragraph($q,'in',$tab[1],$tab[2],$maxNseq,$maxLseq);
                }
                elsif ( $spec =~ /^paragraph_out::(.*)/ ){
                    @tab = split /::/, $spec;
                    $code_h = &print_paragraph($q,'out',$tab[1],$tab[2]);
                }
                elsif ( $spec =~ /^sendto::(.*)$/ ){
                    $sendto = $1;
                }
                elsif ( $spec =~ /^parameter::(.*)/ ){
                    @tab = split /::/, $spec;
                    &print_elem_form($q,$tab[1],$tab[2],$tab[3],$tab[4],$tab[5],$tab[6]);
                }
                elsif ( $spec =~ /^outfile::(.*)$/ ){
                    $spec =~ s{\r}{}g;
                    @tab = split /::/, $spec;
                    $list_visible_outfile = $tab[1];
                    print '                    '.$q->hidden('list_visible_outfile', $list_visible_outfile);
                }
                elsif ( $spec =~ /^config::/ ){
                    $spec = trim($spec);
                    @tab  = split /::/, $spec;
                    if ( $tab[1] eq 'executable' ){
                        $executable = $tab[2];
                        print '                    '.$q->hidden('executable',$executable);
                    }
                    elsif ( $tab[1] eq 'email' ){
                        $mand_email = $tab[2];
                        print '                    '.$q->hidden('mand_email',$mand_email);
                    }
                }
                else {}
            }
        }
    }
    close $CONFIG_FILE;
    $/ = "\n";
    print "                    <hr noshade size=1 width=\"95%\" align=\"center\">
                    <table bgcolor=\"#E1F5FF\" summary=\"PasteEMail\" width=\"100%\">
                      <tr>
                        <td style=\"font-family: arial; color:#333366; font-size:10pt\">You may paste your e-mail address:
                          ",$q->textfield(-name=>'email',-size=>30,-maxlength=>80),"
                        </td>
                      </tr>
                      <tr>
                        <td align=\"right\">
                          <input type=\"submit\" name=\"Submit\" value=\"Submit\">
                          <input type=\"reset\"  name=\"Reset\" value=\"Reset\">
                        </td>
                      </tr>
                    </table>
                  </TD>
                </TR>
              </TABLE>\n";
    return ( $code_h, $sendto );
}

######################################################################
#Error####### ########################################################
######################################################################
sub output_error_page {

    my ( $myerror, $back, $mmode ) = @_;

    print "Content-type: text/html\n\n";
    print "<html>\n";
    print "<head>\n";
    print "<h1><pre>ERROR</pre></h1>\n";
    print "<h3><pre>$myerror\n\n</pre></h3>\n";
    $back =~ s/\/\/\//\//g;
    print "<a href=$back>BACK</a></pre>\n";
    print "</head>\n";
    print "</html>\n";
    exit 1;
    return;
}

######################################################################
#tc_log_analyse ######################################################
######################################################################
sub tc_log_analyse {

    my ( $run_name, $command, $daction, $code_htmle, $sendto ) = @_;

    $file_in    = $run_name.'.tc_log';
    $file_error = $run_name.'.tc_log2';
    $file_out   = $run_name.'.file_result.html';
    if ( -e "$scratch_area/$file_error" ){
        open( my $ERROR_FILE, '<', "$scratch_area/$file_error");
        my $log2beginning='';
        while (<$ERROR_FILE>){
            if ( ( $_ =~/FATAL/ && $_ !~ /blastall/ ) || $_ =~ /TERMINATION STATUS: FAILURE/ || $_ =~/Abnormal Program Termination/ || $flag_error == 1 || $_ =~/Job NOT Completed/ ){
                if ( !-e "$scratch_area/$run_name.fatal" ){
                    `touch $scratch_area/$run_name.fatal`;
                    $tm = localtime;
                    ( $SEC, $MIN, $HOUR, $DAY, $MONTH, $YEAR ) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, $tm->mon, $tm->year);
                    $time_init = timelocal($SEC,$MIN,$HOUR,$DAY,$MONTH,$YEAR);
                    $SEC       = '0'.$SEC   if ( length($SEC) < 2 );
                    $MIN       = '0'.$MIN   if ( length($MIN) < 2 );
                    $HOUR      = '0'.$HOUR  if ( length($HOUR) < 2 );
                    $DAY       = '0'.$DAY   if ( length($DAY) < 2 );
                    $MONTH     = $MONTH+1;
                    $MONTH     = '0'.$MONTH if ( length($MONTH) < 2 );
                    my $who    = `cat $scratch_area/$run_name.tc_log_who`;
                    chomp($who);
                    open( my $FH, '>>', "$home_cgi/crashlog.xml");
                    print {$FH} "<user-call>\n <file>",$run_name,"</file>\n <remote-host>",$who,"</remote-host>\n <user-agent>",$q->user_agent(),"</user-agent>\n <mail>",$email,"</mail>\n <time-stamp>",$HOUR,":",$MIN,":",$SEC," ",$DAY,"/",$MONTH,"/",$YEAR+1900,"</time-stamp>\n</user-call>\n";
                    close $FH;
                }
                if ( $cccpt++ == 0 ){
                    $fatal_error .="<pre>\n$log2beginning";
                }
                my $cleanError = $_;
                $cleanError    =~ s/$home_web\/$dir_tcoffee\///g;
                $cleanError    =~ s/$home_html\///g;
                $cleanError    =~ s/$scratch_area2\/....\///g;
                $cleanError    =~ s/$dir_exe\///g;
                $cleanError    =~ s/Exe\///g;
                $cleanError    =~ s/Tmp\///g;
                $fatal_error  .= $cleanError.'<br>';
                $flag_error    = 1;
            }
            $log2beginning .= $_ if ( !$fatal_error);
        }
        $fatal_error .='</pre>' if ($fatal_error);
        close $ERROR_FILE;
    }

    open( my $TC_COMMAND, '>', "$scratch_area/$file_out".'_command_line');
    $pub_command2 = $command;
#   $pub_command  =~ /^.+;(.+$programme.+).*;*.*$/;
#   $pub_command2 = $command if ($location =~ /^alia/);
#   $pub_command2 = $command if ($location =~ /^frt/);
    $pub_command2 =~ s/$scratch_area\///g;
    $pub_command2 =~ s/$home_web\/$dir_tcoffee\///g;
    $pub_command2 =~ s/$home_html\///g;
    $pub_command2 =~ s/$dir_exe\///g;
    $pub_command2 =~ s/Exe\///g;
    $pub_command2 =~ s/Tmp\///g;
    $pub_command2 =~ s/^.*\.titime//;
    $pub_command2 =~ s/2>.*//;
    $pub_command2 =~ s/tc_log /tc_LOG /;
    print {$TC_COMMAND} $pub_command2;
    close $TC_COMMAND;
    if ( !-s "$scratch_area/$file_in" ){
        $log_empty = 1;
    }

    #test de bonne terminaison de Tcoffee:le fichier de log doit se terminer par T-COFFEE CPU Usage
    open( my $TC_LOG, '<', "$scratch_area/$file_in");
    @tab_file = <$TC_LOG>;
    if ( $location =~ /^frt/ ){
        my $flaggy = 0;
        my @tabFile;
        LOG_TAB:
        for(my $dbut=0; $dbut <= $#tab_file; $dbut++){
            last LOG_TAB if ( $flaggy == 1 && $tab_file[$dbut] =~ /^PS:/ );
            @tabFile = (@tabFile,$tab_file[$dbut]) if ( $flaggy == 1 );
            $flaggy  = 1 if ( $flaggy == 0 && $tab_file[$dbut] =~ /^The output \(if any\) follows:/ );
        }
        @tab_file = @tabFile;
    }
    $y = join('', @tab_file);
    $x = $tab_file[$#tab_file];
    close $TC_LOG;
    #$log_unfinished = 1;
    if ( $x =~ /T-COFFEE CPU Usage/ ){
        #$log_unfinished = 0;
    }
    $y =~s/$home_web\/$dir_tcoffee\///g;
    $y =~s/$home_html\///g;
    $y =~s/$scratch_area2\/....\///g;
    $y =~ s/[^ ]+\/Tmp\///g;
    $y =~ s/$dir_exe\///g;
    $y =~ s/Exe\///g;
    $y =~ s/Tmp\///g;
    open( my $TC_LOG2, '>', "$scratch_area/${run_name}.tc_LOG");
    print {$TC_LOG2} $y;
    close $TC_LOG2;

    my ( $Ok4files, $Ok4receivers ) = ( 0, 0 );
    open( my $TC_LOG, '<', "$scratch_area/$file_in" );
    if ( open(  my $RES, '>', "$scratch_area/$file_out" ) ){
        if ($fatal_error){
            print {$RES} "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
            print {$RES} "<HTML><HEAD><title>Error in $mode :: $level job</title></HEAD>\n\n<BODY BGCOLOR=#FFF0F5>\n";
            print {$RES} "<DIV align=left style=\"font-size: 14pt\">$email</DIV><br>\n" if ($email);
            print {$RES} $fatal_error;
        }
        elsif ( $log_empty || $log_unfinished ){
            print {$RES} "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
            print {$RES} "<HTML><HEAD><title>Crash in $mode :: $level job</title></HEAD>\n\n<BODY BGCOLOR=#FFF0F5>\n<P>";
            print {$RES} "<DIV align=left style=\"font-size: 14pt\">$email</DIV><br>\n" if ($email);
            print {$RES} "Sorry, your job crashed.<br>Either your input was <A HREF=$web_base_doc/doc3.html>not conform</A> or the server is overloaded.<br>Please, check your input, or retry later.</P>\n";
        }
        else {
            $javascriptAVT = "<script language=\"JavaScript\">
function openWin( u ) {
  atv_window = open(\"\", \"atv_window\",
    \"width=300,height=150,status=no,toolbar=no,menubar=no,resizable=yes\");

  atv_window.document.open();

  atv_window.document.write( \"<HTML><HEAD><TITLE>ATV\" );
  atv_window.document.write( \"</TITLE></HEAD><BODY>\" );
  atv_window.document.write( \"<BODY TEXT =\\\"#FFFFFF\\\" BGCOLOR =\\\"#000000\\\">\" );
  atv_window.document.write( \"<FONT FACE = \\\"HELVETICA, ARIAL\\\">\" );
  atv_window.document.write( \"<CENTER><B>\" );
  atv_window.document.write( \"Please do not close this window<br>as long as you want to use ATV.\" );
  atv_window.document.write( \"<APPLET ARCHIVE = \\\"ATVapplet.jar\\\"\" );
  atv_window.document.write( \" CODE = \\\"forester.atv_awt.ATVapplet.class\\\" CODEBASE=\\\"$web_base\\\"\" );
  atv_window.document.write( \" WIDTH = 200 HEIGHT = 50>\" );
  atv_window.document.write( \"<PARAM NAME = url_of_tree_to_load\" );
  atv_window.document.write( \" VALUE = \" );
  atv_window.document.write( \" http://\" + u + \">\" );
  atv_window.document.write( \"</APPLET>\" );
  atv_window.document.write( \"</BODY></HTML>\" );

  atv_window.document.close();
}

</script>
";


            print {$RES} "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";

            my ( $send2section, $form4send2 ) = ( '', '' );
            if ($sendto){
                $form4send2 .="\n<!-- Send to form  -->\n<form method=\"post\" action=\"\">\n";

                $send2section .="\n        <!-- Send results to  -->
        <tr align=\"center\"><td bgcolor=\"#FFFFFF\">&nbsp;</td></tr>
        <tr align=\"center\"><th bgcolor=\"#c0c0c0\" style=\"font-size:15pt\">SEND RESULTS</th></tr>
        <tr align=\"center\"><td>
           <table style=\"empty-cells: hide;\" summary=\"send2\" bgcolor=\"#c0c0c0\" border=\"0\" cellpadding=\"3\" cellspacing=\"1\" width=\"100%\">
";

                my @sendtos = split(/::/, $sendto);
                for (my $u=0;$u <= $#sendtos; $u++){
                    my $my_sent_format = $sendtos[$u];
                    $my_sent_format    =~ s/^.*@//;
                    my $receiver       = $sendtos[$u];
                    $receiver          =~ s/@.*$//;
                    if ( $receiver =~ /^(ProtoGene)/i && -e "$scratch_area/$run_name.$my_sent_format" ){
                        my $title      = 'PROTOGENE: turning amino acid alignments into bona fide CDS nucleotide alignments';
                        $send2section .="             <TR align=center>
               <th align=\"center\" bgcolor=\"#e4eeff\" nowrap width=\"180\">
                 <A href=\"$web_cgi_base/$pg_source\" title=\"$title\"><img src=\"$web_base_images/ProtoGene.gif\" alt=\"$title\" border=0></A>
               </th>
               <td bgcolor=\"#ffffe0\" nowrap valign=\"middle\" width=\"200\">
                 <input name=\"to_ProtoGene\" value=\" to ProtoGene \" type=\"button\" onclick=\"document.forms[0].action='$web_cgi_base/$pg_source'; submit();\">
               </td>
               <td align=left style=\"font-size:10pt\">$title</td>
             </tr>\n";
                        $Ok4receivers = 1;

                        if ( -s "$scratch_area/$run_name.$my_sent_format" ){
                            $form4send2 .= "<INPUT type=hidden name=\"stage1\" value=\"1\">\n";
                            $form4send2 .= "<INPUT type=hidden name=\"daction\" value=\"PROTOGENE::Regular\">\n";
                            $form4send2 .= "<INPUT type=hidden name=\"email\" value=\"$email\">\n" if ($email);
                            open( my $DATA, '<', "$scratch_area/$run_name.$my_sent_format");
                            my $try = join '',<$DATA>;
                            close $DATA;
                            $form4send2 .= hidden(-name=>'--msa',-value=>$try)."\n";
                            $Ok4files = 1;
                        }
                    }
                    elsif ( $receiver =~ /^(MyHits)/i && ( -e "$scratch_area/$run_name.$my_sent_format" || -e "$scratch_area/$run_name.clustalw_aln" ) ){
                        my $title      = 'MyHits: a new interactive resource for protein annotation and domain identification';
                        $send2section .="             <TR align=center>
               <th align=\"center\" bgcolor=\"#e4eeff\" nowrap width=\"180\">
                 <A href=\"http://myhits.isb-sib.ch/\" title=\"$title\"><img src=\"$web_base_images/MyHits-small.gif\" alt=\"$title\" border=0></A>
               </th>
               <td bgcolor=\"#ffffe0\" nowrap valign=\"middle\" width=\"200\">
                 <input name=\"to_MyHits\" value=\" to MSA hub \" type=\"button\" onclick=\"document.forms[0].action='http://myhits.isb-sib.ch/cgi-bin/msa_hub'; submit();\">
               </td>
               <td align=left style=\"font-size:10pt\">$title</td>
             </tr>\n";
                        $Ok4receivers = 1;

                        if ( -s "$scratch_area/$run_name.$my_sent_format" ){
                            open( my $DATA, '<', "$scratch_area/$run_name.$my_sent_format");
                            my $try = join '',<$DATA>;
                            close $DATA;
                            $form4send2 .= hidden(-name=>'text',-value=>$try)."\n";
                            $Ok4files = 1;
                        }
                        elsif ( -s "$scratch_area/$run_name.clustalw_aln" ){
                            open( my $DATA, '<', "$scratch_area/$run_name.clustalw_aln");
                            my $try = join '',<$DATA>;
                            close $DATA;
                            $form4send2 .= hidden(-name=>'text',-value=>$try)."\n";
                            $Ok4files = 1;
                        }
                        if ( -s "$scratch_area/$run_name.dnd" ){
#                           $form4send2 .= "<INPUT type=hidden name=\"dnd_file\"       value=\"";
#                           open( my $DATA, '<', "$scratch_area/$run_name.dnd");
#                           $form4send2 .= join '',<$DATA>, "\">\n";
#                           close $DATA;
                        }
                        if ( -s "$scratch_area/$run_name.score_ascii" ){
#                           $form4send2 .= "<INPUT type=hidden name=\"ascii_file\"       value=\"";
#                           open( my $DATA, '<', "$scratch_area/$run_name.score_ascii");
#                           $form4send2 .= join '',<$DATA>, "\">\n";
#                           close $DATA;
                        }
                    }
                }


                $send2section .= "           </table>
        </td></tr>\n";
            }

            print {$RES} "<html>
<head>
    <title>... $mode :: $level ...</title>
    <link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"$web_base_images/tcoffee.ico\">
    <meta http-equiv=\"Content-Style-Type\" content=\"text/css\">
    <meta http-equiv=\"Content-Script-Type\" content=\"text/javascript\">
</head>

<body background=\"$web_base_images/pattern.gif\" link=\"#014294\" vlink=\"#014294\" alink=\"#014294\" style=\"margin: 0\">\n";

            print {$RES} $form4send2 if ($Ok4files==1 and $Ok4receivers==1);

            my $e_mail = '<br>';
            $e_mail    = "<div align=\"center\" style=\"font-size: 14pt\"><br>$email</div>" if ($email);

            my $noTemplate = '';
            #$noTemplate="<DIV style=\"color:red\">WARNING: NO STRUCTURAL TEMPLATE FOUND WITH A BLAST AGAINST PDB !<br>THE SERVER RAN THE REGULAR T-COFFEE</DIV><br>" if ( $mode =~ /^EXPRESSO/ && (! -e "${scratch_area}/${run_name}_1.template_list" || -z "${scratch_area}/${run_name}_1.template_list"));
            $noTemplate    = "<div style=\"color:red\">WARNING: NO STRUCTURAL TEMPLATE FOUND WITH A BLAST AGAINST PDB !<br>THE SERVER RAN THE REGULAR T-COFFEE</div><br>" if ( $mode =~ /^EXPRESSO/ && !&log_contains_template_file("$scratch_area/$file_in") );
            print {$RES} "
<table summary=\"menu\" align=\"center\" bgcolor=\"#ffffff\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"792\" style=\"border-style:double;border-width:5px;border-color:#c0c0c0;border-top-style:none;\">
  <tr valign=\"top\" align=\"center\">
    <td width=\"750\" align=\"center\">
      <table summary=\"logo\" align=\"center\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"792\">
        <tr>
          <td align='center' nowrap>$logo</td>
        </tr>
        <tr align='center'>
          <td align='right' bgcolor='#000000' height='6' valign='top' style='color: #ffffff;'>
            <a href=\"$web_cgi_base/$pg_source\" style=\"color: #ffffff;\">HOME</a>&nbsp;|&nbsp;
            <a href=\"$tcoffee_homepage\" style=\"color: #ffffff;\">references</a>&nbsp;|&nbsp;
            <a href=\"$web_base_doc/doc3.html\" target=\"_blank\" style=\"color: #ffffff;\">help</a>&nbsp;|&nbsp;\n";

            print {$RES} '            ', &env_param::transpose_eMails_delayed($webmaster, "$web_base_images/mail1.gif"), "&nbsp;&nbsp;&nbsp;
            <a href=\"$web_cgi_base/${pg_source}?load_cfg_file=1\" style=\"text-decoration: none; color: #000000;\">.</a>
          </td>
        </tr>
      </table>

    </td></tr><tr valign=\"top\" align=\"center\">
    <td width=\"750\" align=\"center\">

      <div align=center style=\"font-size: 20pt\"> $e_mail <b>$mode :: $level</b> </div>\n";

            my @messages;
            @messages = &alertMessages if ( -e 'Info.txt' );
            if ( exists($messages[0]) ){
                my $safe_mode = $mode;
                $safe_mode =~ s/\(.*$//;
                my @serverMessages = grep(/\G^$safe_mode\E/i, @messages);
                if ( exists($serverMessages[0]) ){
                    my @warnings = map( m/\G^$safe_mode\:\:(.*)\E/i, @serverMessages);
                    print {$RES} "      <h4 style=\"color:#FF0000\">",join('<br>', @warnings),"</h4>\n" if ( exists($warnings[0]) );
                }
                else { print {$RES} "      <br>\n"; }
            }
            else { print {$RES} "      <br>\n"; }

            print {$RES} "      <div align=center style=\"font-size:small\">
      Your data will remain available on this server over the next $tmpOldFiles days. It will then be deleted.<br>
      Do not forget to bookmark this <a href=\"$file_out\">URL</a> or save it for further reference.
      <br><br></div>
    </td></tr><tr valign=\"top\" align=\"center\">
    <td width=\"750\" align=\"center\">

      $noTemplate
      <table summary=\"results\" bgcolor=\"#c0c0c0\" border=\"1\" cellpadding=\"0\" cellspacing=\"0\" width=\"100%\" align=\"center\">
        <tr align=\"center\"><th style=\"font-size:15pt\">RESULTS</th></tr>
        <tr align=\"center\"><td>

        <!-- Result tab -->
           <table style=\"empty-cells: hide;\" summary=\"resultats\" bgcolor=\"#c0c0c0\" border=\"0\" cellpadding=\"3\" cellspacing=\"1\" width=\"100%\">
";

            $flag = 0;
            $dssp = '';
            while (<$TC_LOG>){

                if (/\s*Input File (.*) Seq .* type PROTEIN Struct Yes PDBID (.*)/){
                    $pdb_file = $1;
                    $pdb_id   = $2;
                    $dssp    .= $pdb_file.'__'.$pdb_id.'::';
                }
                #elsif(/OUTPUT RESULTS/) {
                    #$flag++;
                #}
                else {
                    if (/File Type/){
                        $line   = $_;
                        $line   =~ /.*Type=\s*(.*)\s*Format=\s*(.*)\s*Name=\s*([^\s]*)/;
                        $type   = $1;
                        $format = $2;
                        $nom    = $3;
                        $format =~ s/ //;
                        if ( $nom =~ /.*\/(.*)$/ ){
                            $nom = $1;
                        }

                        $res_template = '';
                        if ( $type =~ /Template_List/ ){
                            $res_template = "<TD bgcolor=\"#FFFFE0\" valign=top nowrap><A HREF=\"$web_base_link_tmp/"."$nom\">Template_file</A></TD>\n";
                        }
                        if ( $type =~ /MSA/ || $type=~ /APDB_RESULT/ ){
                            $res_table_msa .= "               <TD bgcolor=\"#FFFFE0\" valign=top nowrap><A HREF=\"$web_base_link_tmp/"."$nom\">$format</A></TD>\n";
                        }
                        elsif ( $type =~ /TREE/ && ($line !~ /NOT PRODUCED/ && $line !~ /Name= no/) ){
                            if ( $type =~ /GUIDE_TREE/ ){
                                $res_table_tree .= "<TD bgcolor=\"#FFFFE0\" valign=top nowrap><A HREF=\"$web_base_link_tmp/"."$nom\">dnd file</A></TD>\n";
                            }
                            else {
                                $res_table_tree .= "               <TD bgcolor=\"#FFFFE0\" valign=top nowrap><A HREF=\"$web_base_link_tmp/"."$nom\">ph file</A></TD>\n";
                            }
                        }
                    }
                }
            }
            if ($res_table_msa){
                print {$RES} "             <TR align=center>
               <TH bgcolor=\"#E4EEFF\" align=\"left\" nowrap width=\"180\" style=\"font-family: arial; color:#817f82; font-size:14pt\">Multiple Alignment</TH>
$res_table_msa             </TR>\n";
            }
            if ( ($daction =~ /Make/ || $daction =~ /3D/) && -e "${scratch_area}/${run_name}.dnd" && -s "${scratch_area}/${run_name}.dnd" ){
            }
            if ($res_table_tree){
                print {$RES} "             <TR align=center>
               <TH bgcolor=\"#E4EEFF\" align=left nowrap style=\"font-family: arial; color:#817f82; font-size:14pt\">Tree</TH>
               $res_table_tree             </TR>\n";
            }
            print {$RES} "             <TR align=center>
               <TH bgcolor=\"#E4EEFF\" align=left nowrap style=\"font-family: arial; color:#817f82; font-size:14pt\">System files</TH>
               <TD bgcolor=\"#FFFFE0\" valign=top nowrap><A HREF=\"${run_name}.tc_LOG\">LOG</A></TD>
               <TD bgcolor=\"#FFFFE0\" valign=top nowrap><A HREF=\"$file_out"."_command_line\">Command line</A></TD>
               $res_template             </TR>\n";

            foreach $item_input (@tab_list_input_files){
                $list_link_items .= "               <TD bgcolor=\"#FFFFE0\" width=\"200\"><A HREF=\"$item_input\">$item_input</A></TD>\n"
            }
            print {$RES} "             <TR align=center>
               <TH bgcolor=\"#E4EEFF\" align=left nowrap style=\"font-family: arial; color:#817f82; font-size:14pt\">Inputs</TH>
$list_link_items             </TR>\n";

            print {$RES} "           </table>
        </td></tr>";

            print {$RES} $send2section if ($Ok4files==1 and $Ok4receivers==1);

            print {$RES} "      </table>
      <table summary='paragraph_out' align='center' bgcolor='white' width='100%'>
        <tr><td align='center'><br>\n";
            if ($code_htmle){
                print {$RES} "$code_htmle
            </td></tr>
        <tr><td align='center'>\n";
            }


            $les_arg = "stage1=1&amp\;daction=$mode".'::'.$level;
            my $other_level = $level eq 'Regular' ? 'Advanced'
                            :                       'Regular';
            my $other_arg   = "stage1=1&amp\;daction=$mode".'::'.$other_level;
            $lien    = "            <br>Home Server: \n            <a href='$web_cgi_base/$pg_source?$les_arg'>$mode :: $level</a> | <a href='$web_cgi_base/$pg_source?$other_arg'>$mode :: $other_level</a>";
            print {$RES} "          <center style='font-size:14pt'>\n$lien
          </center></td>
        </tr>
      </table>
    </td>
  </tr>
</table><br>\n";

            #test postprocessing
            if ( !$fatal_error ){
                if (1==2){
                    print {$RES} "<CENTER><HR noshade size=1 width='95%' align='center'><H2><font face=helvetica>Post processing</font></H2>";
                    print {$RES} "<table cellpadding=0 cellspacing=0 border=0 bgcolor=\"#C0C0C0\" summary=\"results\"><tr><td><table cellpadding=3 cellspacing=1 border=0 summary=\"resultats\">";
                    if ($dssp){
                        $dssp =~ s/\s//g;
                        if ( $daction =~ /3D Coffee::Regular/ ){
                            $arg_dssp = "&dssp=$dssp&mode=regular";
                        }
                        else {
                            $arg_dssp = "&dssp=$dssp";
                        }
                    }
                }
            }
        }
        print {$RES} "</form>\n" if (($Ok4files==1 and $Ok4receivers==1) and !$fatal_error and !$log_empty and !$log_unfinished);
        print {$RES} "\n</BODY>\n</HTML>\n";
        close $RES;
    }
    else {
        die "Unable to write file: $!\n";
        exit(1);
    }
    close $TC_LOG;
    return;
}

######################################################################
#Verif Mail###########################################################
######################################################################
sub ValidEmailAddr {

    my $mail = shift;

    #characters allowed on name: 0-9a-Z-._ on host: 0-9a-Z-. on between: @
    #return 0 if ( $mail !~ /^[0-9a-zA-Z\.\-\_]+\@[0-9a-zA-Z\.\-]+$/ );
    return 0 if ( $mail !~ /^[\w\.\-\_]+\@[\w\.\-\_]+\.[a-zA-Z][a-zA-Z][a-zA-Z]?$/ );
    return 1;
}

######################################################################
#trim#################################################################
######################################################################
sub trim {

    my @out = @_;
    for (@out) {
        s/^\s+//;
        s/\s+$//;
    }
    return wantarray ? @out : $out[0];
}

######################################################################
#Erase Old Files######################################################
######################################################################

sub tmpManagement{

    my @temp_files = glob(${scratch_area}.'/tcf*.in0');
    my @tempFiles  = grep { -M $_ > $tmpOldFiles } @temp_files; #Sort list only with $tmpOldFiles days old files

    if ( exists ($tempFiles[0]) ){
        while(<@tempFiles>){
            my $appli = $_;
            $appli    =~ s/\.in0$//;
            unlink( glob($appli.'*') );
        }
    }
    return;
}

sub reorder_results {

    my ($tmp) = @_;

    my @applications = &env_param::server_list($config_infile); #Add automaticaly servers from config file

    my $status = 0;
    SUBDIR_PER_APPLI:
    for my $appli (@applications){
        unlink "$tmp/$appli" if ( -e "$tmp/$appli" && -f "$tmp/$appli" );
        mkdir "$tmp/$appli" if ( !-e "$tmp/$appli" );
        $status += system("ln -s \"$tmp/index.html\" \"$tmp/$appli/\"") if ( !-l "$tmp/$appli/index.html" );
        $status += system("ln -s \"$tmp/.t_coffee\" \"$tmp/$appli/\"") if ( !-l "$tmp/$appli/.t_coffee" );
    }
    print {*STDERR} "Links for subdirectories failed\n" if ( $status );

    return;
}

######################################################################
#Print alert messages if any##########################################
######################################################################

sub alertMessages {

    if ( -z 'Info.txt' ){
        unlink('Info.txt');
        return;
    }
    else{
        open( my $ALERT, '<', 'Info.txt');
        my @messages;
        ALERT_LINES;
        while(<$ALERT>){
            if ( $_ =~ /^@(.+::.*)$/ ){
                @messages=(@messages,$1);
            }
        }
        close $ALERT;
        return @messages;
    }
}


######################################################################
# HTML standard header, from doctype to body #########################
######################################################################

sub std_HTML_header {

    my ($script) = @_;
    my $HTML_header = '';

    $HTML_header .= "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">
<html>
<head>
    <title>Tcoffee</title>
    <link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"$web_base_images/tcoffee.ico\">
    <meta name=\"keywords\" content=\"Multiple sequence alignment,T-COFFEE,TCoffee,t_coffee,MSA,alignment,protogene,proto-gene,3dcoffee,3d-coffee,expresso,APDB,iRMSD,RCoffee,R-Coffee\">
    <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">
    <meta http-equiv=\"Content-Script-Type\" content=\"text/javascript\">
    <meta http-equiv=\"Content-Style-Type\" content=\"text/css\">
    $script
</head>

<body background=\"$web_base_images/pattern.gif\" style=\"margin: 0\">
";

    return $HTML_header;
}


######################################################################
# Tree in PDF ########################################################
######################################################################

#sub treeInPDF{
#    my ($file,$tp_dir,$mime_type)=@_;
#
#    if ($file ne '' and $tp_dir ne ''){ # and $mime_type ne ''){
#       unlink("$tp_dir/plotfile","$tp_dir/outtree","$tp_dir/intree");
#       if ( -e "$tp_dir/$file"){
#          system("/usr/bin/env perl -i -pe 's/ //g;' $tp_dir/$file");
#          open (CF2,">$tp_dir/$file.conf2");
#          print CF2 "$tp_dir/$file\n$tp_dir/../fontfile\nV\nN\nS\nC\nH\nY\n" if ($location =~ /^frt/);
#      print CF2 "$tp_dir/$file\n$tp_dir/../fontfile\nL\nN\n1\n4\n90\nY\n" if ($location =~ /^alia/);
#          close CF2;
#
#          system("export LD_LIBRARY_PATH=$specific_lib:\$LD_LIBRARY_PATH; $dir_exe/drawgram < $tp_dir/$file.conf2 >/dev/null; mv plotfile $tp_dir/$file.ps");
#          system("/usr/bin/env ps2pdf $tp_dir/$file.ps $tp_dir/$file.pdf");
#          unlink("$tp_dir/plotfile","$tp_dir/outtree","$tp_dir/$file.conf2","$tp_dir/$file.ps","$tp_dir/intree");
#
#          return("<TD bgcolor=\"#FFFFE0\" valign=top nowrap><A HREF=$web_base_link_tmp/"."$file.pdf>pdf file</A></TD>");
#       }
#       else{
#          return('');
#       }
#    }
#    else{
#       return('');
#    }
#}

######################################################################
# Log Analysis
######################################################################
sub log_contains_template_file{

    my ($log) = @_;

    open ( my $LOG_TEST, '<', "$log");
    while (<$LOG_TEST>){
        if ( /Template_List/ ){
            close $LOG_TEST;
            return 1;
        }
    }
    close $LOG_TEST;
    return 0;
}

