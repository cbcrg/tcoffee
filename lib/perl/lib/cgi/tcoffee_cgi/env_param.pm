#File env_param.pm  for 'IGS' and 'Vital-It'

######################################################################
#2006-06-26                                                          #
#Sebastien MORETTI, moretti.sebastien [AT] gmail.com, SIB            #
#Olivier POIROT, poirot [AT] igs.cnrs-mrs.fr, Labo IGS, UPR2589      #
#Web interface for Tcoffee and other programs                        #
#Need feedback.cgi, configuration_file.txt                           #
######################################################################

package env_param;

our $tcoffee_homepage = 'http://www.tcoffee.org/Projects_home_page/t_coffee_home_page.html';
##############  Configuration  ##############
our ( $config_infile, $pg_source,    $tmpOldFiles,     $dir_pdb,      $webmaster,        $fromEMail, $ppath                    ) = ( '', '', '', '', '', '', '' );
###############  Directories  ###############
our ( $dir_tcoffee,   $dir_cgi,      $dir_tmp,         $dir_images,   $dir_doc,          $dir_exe,   $programme, $specific_lib ) = ( '', '', '', '', '', '', '', '' );
############  Server directories ############
our ( $home_web,      $home_html,    $home_dir,        $home_cgi,     $scratch_area,     $scratch_area2                        ) = ( '', '', '', '', '', '' );
#############  Web directories #############
our ( $web_base,      $web_cgi_base, $web_base_images, $web_base_doc, $web_base_link_tmp                                       ) = ( '', '', '', '', '' );
##################  Look  ##################
our ( $colorTitle,    $colorAppli,   $Rgb,             $rGb,          $rgB,              $logo                                 ) = ( '', '', '', '', '', '' );

sub variables{
    my ($hostname) = @_;

    if ( $hostname =~ /^alia/ ){ #for IGS lab computer, alia host
        #############################################
        ##############  Configuration  ##############
        #############################################
        $config_infile = 'configuration_file.txt';
        $pg_source     = 'index.cgi';
        $tmpOldFiles   = 9; #Keep only files more recent than $tmpOldFiles days old
        $dir_pdb       = '/banks/_Structures/PDB/all/pdb';
        $webmaster     = 't_coffee_webmaster@igs.cnrs-mrs.fr,cedric.notredame@europe.com,moretti.sebastien@gmail.com';
        $fromEMail     = $webmaster;
        $ppath         = "\$PATH";

        #############################################
        ###############  Directories  ###############
        #############################################
        $dir_tcoffee  = 'Tcoffee';
        $dir_cgi      = 'tcoffee_cgi';
        $dir_tmp      = 'Tmp';
        $dir_images   = 'Images';
        $dir_doc      = 'Doc';
        $dir_exe      = '/home/igs/public_html/Tcoffee/Exe';
        $programme    = $dir_exe.'/t_coffee';
        $specific_lib = '';

        #############################################
        ############  Server directories ############
        #############################################
        $home_web      = '/home/igs/public_html';
        $home_html     = '/home/igs/public_html'; #Add for Vital-It
        $home_dir      = "$home_web/$dir_tcoffee";
        $home_cgi      = "$home_dir/$dir_cgi";
        $scratch_area  = "${home_dir}/${dir_tmp}";
        $scratch_area2 = $scratch_area;

        ############################################
        #############  Web directories #############
        ############################################
        $web_base          = "http://www.igs.cnrs-mrs.fr/$dir_tcoffee";
        $web_cgi_base      = "$web_base/$dir_cgi";
        $web_base_images   = "$web_base/$dir_images";
        $web_base_doc      = "$web_base/$dir_doc";
        $web_base_link_tmp = "$web_base/$dir_tmp";

        ############################################
        ##################  Look  ##################
        ############################################
        $colorTitle = '#cc9966';
        $colorAppli = '#e7c7a8';
        $Rgb        = 231;
        $rGb        = 199;
        $rgB        = 168;
        $logo = '<a href="http://www.cnrs.fr/"><img src="'.$web_base_images.'/logo_cnrs.gif" width=160 alt="Menu " border="0"></a><a href="http://www.igs.cnrs-mrs.fr/"><img src="'.$web_base_images.'/logo_igs.jpg" width=590 alt=" Menu" border="0"></a>';
    }
    elsif ( $hostname =~ /^frt/ ){ #for Vital-IT computer, on frt host
        #############################################
        ##############  Configuration  ##############
        #############################################
        $config_infile = 'configuration_file.txt';
        $pg_source     = 'index.cgi';
        $tmpOldFiles   = 9;    #Keep only files more recent than $tmpOldFiles days old
        $dir_pdb       = '/db/pdb';
        $webmaster     = 'tcoffee@vital-it.ch,cedric.notredame@europe.com,moretti.sebastien@gmail.com';
        $fromEMail     = 'tcoffee@vital-it.ch';
        $ppath         = "\\\$PATH";

        #############################################
        ###############  Directories  ###############
        #############################################
        $dir_tcoffee  = 'Tcoffee';
        $dir_cgi      = 'tcoffee_cgi';
        $dir_tmp      = 'Tmp';
        $dir_images   = 'Images';
        $dir_doc      = 'Doc';
        $dir_exe      = '/mnt/local/bin';
        $programme    = $dir_exe.'/t_coffee';
        $specific_lib = '/mnt/local/lib';

        #############################################
        ############  Server directories ############
        #############################################
        $home_web      = '/var/www/cgi-bin/tcoffee';
        $home_html     = '/var/www/html/tcoffee'; #Add for Vital-It
        $home_dir      = "$home_web/$dir_tcoffee";
        $home_cgi      = "$home_dir/$dir_cgi";
        $scratch_area  = "${home_html}/${dir_tmp}";        #Use different scratch places because Vital-IT frontal can see both
        $scratch_area2 = '/scratch/frt/tcoffee/'.$dir_tmp; # (more one than the other) but nodes can only see one of them (the other)

        ############################################
        #############  Web directories #############
        ############################################
#       $web_base          = "http://www.vital-it.ch/prd/smoretti/cgi-bin/$dir_tcoffee";
        $web_base          = "http://tcoffee.vital-it.ch/cgi-bin/$dir_tcoffee";
        $web_cgi_base      = "$web_base/$dir_cgi";
#       $web_base_images   = "http://www.vital-it.ch/prdpub/smoretti/$dir_images";
        $web_base_images   = "http://tcoffee.vital-it.ch/$dir_images";
#       $web_base_doc      = "http://www.vital-it.ch/prdpub/smoretti/$dir_doc";
        $web_base_doc      = "http://tcoffee.vital-it.ch/$dir_doc";
#       $web_base_link_tmp = "http://www.vital-it.ch/prdpub/smoretti/$dir_tmp";
        $web_base_link_tmp = "http://tcoffee.vital-it.ch/$dir_tmp";

        ############################################
        ##################  Look  ##################
        ############################################
        $colorTitle = '#84c2e5';
        $colorAppli = '#a5daf8';
        $Rgb        = 165;
        $rGb        = 218;
        $rgB        = 248;
        $logo = '
            <table width="720" border="0" cellspacing="0" cellpadding="0" bgcolor="#FFFFFF" align="center" summary="bandeauA">
              <tr>
                <td height="7"><img src="'.$web_base_images.'/top-couleurs.gif" alt="couleurs" width="720" height="7"></td>
              </tr>
              <tr>
                <td height="26"><img src="'.$web_base_images.'/top-titre-sib.gif" alt="Titre" width="720" height="26"></td>
              </tr>
            </table>
            <table width="720" border="0" cellspacing="0" cellpadding="0" align="center" bgcolor="#FFFFFF" summary="bandeauB" style="height:81">
              <tr>
                <td width="20"><img src="'.$web_base_images.'/topdessin-gauche.jpg" alt="lf" width="20" height="81"></td>
                <td width="160" valign="middle" align="center" style="background-image:url('.$web_base_images.'/topdessin-fondlogo.jpg);background-repeat:no-repeat;">
                  <a href="http://www.isb-sib.ch/"><img src="'.$web_base_images.'/top-logoSIB10.gif" alt="SIB_" width="112" height="70" border=0></a></td>
                <td width="540"><a href="http://www.vital-it.ch/"><img src="'.$web_base_images.'/topdessin-droit.jpg" border=0 alt="_Vital-IT" width="540" height="81"></a></td>
              </tr>
            </table>';
    }
    # With your_host is the first part of the 'hostname' command result on your computer
    # e.g.: myhost.mydomain -> your_host = myhost
    #       myhost          -> your_host = myhost
    elsif ( $hostname =~ /^your_host/ ){
        #############################################
        ##############  Configuration  ##############
        #############################################
        $config_infile = 'configuration_file.txt';
        $pg_source     = 'index.cgi';
        $tmpOldFiles   = 9;    #Keep only files more recent than $tmpOldFiles days old
        $dir_pdb       = '/banks/_Structures/PDB/all/pdb';
        $webmaster     = 't_coffee_webmaster@igs.cnrs-mrs.fr';
        $fromEMail     = $webmaster;
        $ppath         = "\$PATH";

        #############################################
        ###############  Directories  ###############
        #############################################
        $dir_tcoffee  = 'Tcoffee';
        $dir_cgi      = 'tcoffee_cgi';
        $dir_tmp      = 'Tmp';
        $dir_images   = 'Images';
        $dir_doc      = 'Doc';
        $dir_exe      = '/home/igs/public_html/Tcoffee/Exe';
        $programme    = $dir_exe.'/t_coffee';
        $specific_lib = '';

        #############################################
        ############  Server directories ############
        #############################################
        $home_web     = '/home/igs/public_html';
        $home_html    = '/home/igs/public_html'; #Add for Vital-It
        $home_dir     = "$home_web/$dir_tcoffee";
        $home_cgi     = "$home_dir/$dir_cgi";
        $scratch_area = "${home_dir}/${dir_tmp}";

        ############################################
        #############  Web directories #############
        ############################################
        $web_base          = "http://www.igs.cnrs-mrs.fr/$dir_tcoffee";
        $web_cgi_base      = "$web_base/$dir_cgi";
        $web_base_images   = "$web_base/$dir_images";
        $web_base_doc      = "$web_base/$dir_doc";
        $web_base_link_tmp = "$web_base/$dir_tmp";

        ############################################
        ##################  Look  ##################
        ############################################
        $colorTitle = '#7f31e0';
        $colorAppli = '#b47df8';
        $Rgb        = 180;
        $rGb        = 125;
        $rgB        = 248;
        $logo = '<a href="http://www.cnrs.fr/"><img src="'.$web_base_images.'/logo_cnrs.gif" width=160 alt="Menu " border="0"></a><a href="http://www.igs.cnrs-mrs.fr/"><img src="'.$web_base_images.'/logo_igs.jpg" width=590 alt=" Menu" border="0"></a>';
    }
    return;
}


sub goodFlag{

    my ($hostname) = @_;

    if ( $hostname =~ /IGS/i ){
        return('http://www.igs.cnrs-mrs.fr/Tcoffee/');
    }
    elsif ( $hostname =~ /HOME/i ){
        return('http://www.tcoffee.org/');
    }
    elsif ( $hostname =~ /Vital-IT/i || $hostname eq 'Suisse' ){
        return('http://tcoffee.vital-it.ch/');
    }
    elsif ( $hostname =~ /EBI/i ){
        return('http://www.ebi.ac.uk/t-coffee/');
    }
    return;
}

######################################################################
# List servers #######################################################
######################################################################

sub server_list {
    my ($config_file) = @_;

    my @list;
    if ( -r "$config_file" ){
        open(my $CONF, '<', "$config_file") or die "Cannot open the configuration file\n";
        @list = grep( s{^server::([^:][^:][^:])[^:]*::([^:])[^:]*::.*$}{$1$2}, <$CONF>);
        close $CONF;
    }
    @list = (@list, 'PerR', 'PerA');
    chomp(@list);
    return @list;
}

######################################################################
# Hide e-mail addresses ##############################################
######################################################################

#cf. http://www.pgregg.com/projects/encode/htmlemail2.php & http://www.pgregg.com/projects/encode/htmlemail.php
sub transpose_eMails_delayed {
    my ($emails, $imgPath) = @_;

    my $href = "mailto:$emails"; #For mailto href here only
    my $text = '<img src="'.$imgPath.'" alt="e-mail" border="0">';
    my $code = sprintf("function seb_transpose2(h) {var s='%s';var r='';for(var i=0;i<s.length;i++,i++){r=r+s.substring(i+1,i+2)+s.substring(i,i+1)}h.href=r;}document.write('<a href=\"#\" onMouseOver=\"javascript:seb_transpose2(this)\" onFocus=\"javascript:seb_transpose2(this)\">%s</a>');", &transpose($href), $text);
    my $userCode = sprintf("%s%s%s", "<script type=\"text/javascript\">eval(unescape('", &escapeencode($code), "'))</script><noscript>&nbsp;</noscript>");

    return $userCode;
}

sub escapeencode {
    my ($strg) = @_;

    my $ret = '';
    my @arr = unpack('C*', $strg);
    for my $char (@arr){
        $ret .= sprintf("%%%X", $char);
    }

    return $ret;
}

sub transpose {
    my ($strg) = @_;

    # function takes the string and swaps the order of each group of 2 chars
    my $len = length($strg);
    my $ret = '';
    for (my $i=0; $i<$len; $i=$i+2){
        if ($i+1 == $len){
            $ret .= substr($strg, $i, 1);
        }
        else{
            $ret .= sprintf("%s%s", substr($strg, $i+1, 1), substr($strg, $i, 1));
        }
    }

    return $ret;
}

1;


