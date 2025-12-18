#HTC condor job

#Executable  = /bin/zsh
#arguments = one.sh

Executable  = ./one

request_cpus = 16

Log         = one.log
Output      = one.out
Error       = one.err

getenv = True

# as long as jobs run at DESY Zeuthen batch, please disable file transfer feature
should_transfer_files = no
# request 2GB RAM
request_memory = 16GB

Queue 1
