## this file tests the gtdownload command...
## /work/01839/dakotaz/cghub/bin/gtdownload is the executable
## -v is the verbose option so that I know what the program is doing
## -c is the credential to download protected data, good until May 1, 2016
## -p is the path to where we want to store the file being downloaded
## -d is the file to download, tested here with a one file manifest

manifest=$1
path=$2

time /work/01839/dakotaz/cghub/bin/gtdownload -v -c /work/01839/dakotaz/cghub.key -p $path -d $manifest