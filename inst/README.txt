
README file for the maanova package
----------------------------------------------------------------------

This is a quick document to explain the installation and use
of the "maanova" package for R


Installation (Windows)

  1. Unzip the "maanova.zip" file into the directory $RHOME\library
     ($RHOME is something like c:\Program Files\R\rw1030)
     Note that this should create a directory $RHOME\library\maanova
     containing the R source code and the compiled dll

  2. Start Rgui 

  3. Type "link.html.help()" to get the help files for the maanova 
     package added to the help indices

  4. Note that the source code is in the file "maanova_*.tar.gz"


Installation (Unix)

  1. Go into the directory containing "maanova_*.tar.gz"

  2. Type "R INSTALL maanova". Note that this command will install
     the package in $RHOME directory (normally /usr/lib/R/library).
     To install the package in user's own directory, type 
     "R INSTALL -l ~/Rlib/ maanova". Then create a file .Renviron 
     in user's home directory with the following line:
         R_LIBS=/home/<user>/Rlib 


----------------------------------------------------------------------

end of README.txt
