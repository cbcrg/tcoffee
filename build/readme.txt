==========================================================
                 T-Coffee binary installer
==========================================================

The installer has finished the installation process. 


READ THIS BEFORE START
----------------------

T-Coffee is a command line tool and you will use it from the shell terminal typing the command 't_coffee'. 

To test that T-Coffee has been installed successfully make sure to open a new shell terminal (so that changes made in the installation procedure to the system environment will be effective) and type the following command followed by the enter key: 

t_coffee -version 

The application must respond with the version number of the T-Coffee release you have just installed on your system. 


To run your first basic T-Coffee alignment, you may try the following example: 

t_coffee <your-sequences-file>


Check the home page to get more information on how to use T-Coffee http://www.tcoffee.org/Projects_home_page/t_coffee_home_page.html


WHAT HAS BEEN CHANGED ON MY SYSTEM 
----------------------------------

The installer copies all required T-Coffee binaries and third party tools to the target folder that you have chosen during the installation procedure. 

The following environment configuration files are updated (depending on your system): 

[Linux]
  ~/.bashrc	 		
	
[Mac OSX]	
  ~/.profile  
  ~/.MacOSX/environment.plist	


The following environment variables are added or updated:    

  DIR_4_COFFEE: the T-Coffee base installation directory;
  CACHE_4_TCOFFEE: points to cache folder used by T-Coffee;
  EMAIL_4_TCOFFEE: contains the email address you have entered, required by the remote BLAST server; 
  PERL5LIB: is updated to include the directory containing Perl modules required by T-Coffee;
  MAFFT_BINARIES: points to the directory containing MAFFT tool binaries;
  DYLD_LIBRARY_PATH: gfortran runtime libraries (only on Mac OS X);
  PATH: is prepended with the path to the newly installed T-coffee version.



WHAT IS INSTALLED 
-----------------

T-Coffee is installed to the directory $HOME/tcoffee/<version-number> unless you have specified a different one. 

The installation directory contains the following folders: 

	/bin		: contains the main T-Coffee application file;
	/methods	: (empty, only for compatibility with previous versions)
	/mcoffee	: data files required by M-Coffee special mode;
	/perl		: perl modules required by T-Coffee;
	/plugins	: contains third party tools used by T-Coffee special modes;
	/gfortran	: gfortran runtime library (only on Mac OS X);
	/tmp		: where T-Coffee temporary files will be placed.
	
	
	
HOW TO UNINSTALL T-COFFEE 
----------------------

The installer backs up your environment configuration file, making a copy of it adding to the original name the suffix .bak<n>. 

In order to un-install T-Coffee rollback your configuration restoring your previous environment configuration file.

Moreover you may want also to remove all installer binaries deleting the T-Coffee installation directory.
