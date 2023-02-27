#set path=(. $HOME/anaconda3/bin $HOME/bin/ $HOME/scripts /programs/bin/x86_64 /programs/bin/noarch /programs $HOME/ ~/Programs ~/Programs/ncbi-blast-2.2.29+/bin $path)
#setenv PATH $PATH\:$HOME/bin/\:$HOME/.local/bin/\:$HOME/scripts\:/programs/bin/x86_64\:/programs/bin/noarch\:/programs\:$HOME/\:~/Programs\:~/Programs/ncbi-blast-2.2.29+/bin\:/data02/gloiseau/bin\:.
#set path=(. $HOME/bin/ $HOME/scripts /programs/bin/x86_64 /programs/bin/noarch /programs $HOME/ ~/Programs ~/Programs/ncbi-blast-2.2.29+/bin $path)
#set verbose
#alias
#eval `dircolors --csh ~/.dircolors`
#alias ls 'ls --color=tty'
#alias grep 'grep --color=auto'
#alias condor 'ssh bmueller@submit.chtc.wisc.edu'
alias mv 'mv -i'
alias cp 'cp -i'

set ignorecase
set smartcase

#TODO: CHANGE LINES 18, 25, AND 26 IN FILE IF MSLIB IS EVER MOVED TO A NEW DIRECTORY
# MSL ENVIRONMENTAL VARIABLES FOR INCLUSION OF EXTERNAL LIBRARIES
setenv MSL_DIR "/mnt/d/github/mslib/trunk_AS"
setenv MSL_DIR "/mnt/d/github/mslib/trunk_AS"

setenv MSL_GSL T
setenv MSL_GLPK F
setenv MSL_BOOST F
setenv MSL_DEBUG F
setenv MSL_R F
setenv MSL_EXTERNAL_LIB_DIR /mnt/d/github/mslib/ext_libs
setenv MSL_EXTERNAL_INCLUDE_DIR /mnt/d/github/mslib/ext_includes

#setenv MSL_CHARMM_TOP "/exports/home/scondon/mslib/trunk/toppar/charmm22.top"
#setenv MSL_CHARMM_PAR "/exports/home/scondon/mslib/trunk/toppar/charmm22.par"
#setenv MSL_CHARMM_SOLV "/exports/home/scondon/mslib/trunk/library/charmmTopPar/solvpar22.inp"
#setenv MSL_HBOND_PAR "/exports/home/scondon/mslib/trunk/toppar/scwrl4hb/par_hbond_2.txt"
#setenv MSL_HBOND_CA_PAR "/exports/home/scondon/mslib/trunk/toppar/scwrl4hb/par_hbond_CA_2.txt"
#setenv MSL_ROTLIB "/exports/home/scondon/mslib/trunk/library/EBL_11-2011_CHARMM22.txt"
###########################################################
#
#       MAKE PROMPT "hostname path>"

# VERSION 1 simple
#alias setprompt 'set prompt="`hostname -s` `pwd`> "'
#if ($?prompt) then
#       setprompt
#       alias cd 'cd \!*; setprompt'
#endif

# VERSION 2 colored 
#alias setprompt 'set prompt="%{^[[34m%}`hostname -s` %{^[[33m%}`pwd`%{^[[31m%}>%{^[[0m%} "'
#if ($?prompt) then
#       setprompt
#       alias cd 'cd \!*; setprompt'
#endif

# VERSION 3 colored and truncated to 50 characters
set dollarone = '$1'
set dotdotdot = '"..."'
#alias trimmedpath "pwd|gawk '{if (length($dollarone) <= 50) print $dollarone; else print substr($dollarone,1,7) $dotdotdot substr($dollarone,length($dollarone) -39,40)}'"//gawk no longer works on 20.04 to show pwd on the terminal, but removing it does
alias trimmedpath "pwd '{if (length($dollarone) <= 50) print $dollarone; else print substr($dollarone,1,7) $dotdotdot substr($dollarone,length($dollarone) -39,40)}'"
alias setprompt 'set prompt="%{^[[34m%}`hostname -s` %{^[[33m%}`trimmedpath`%{^[[31m%}>%{^[[0m%} "'
if ($?prompt) then
        setprompt
        alias cd 'cd \!*; setprompt'
endif

###########################################################
#	Idiot-proofing the terminal
###########################################################
set noclobber
set rmstar

###########################################################
#	History settings
###########################################################
#set history=10000
#set savehist=10000 merge
#set histfile=~/.csh_history
#
## Permanently store commands in logfiles
#set FULL_CMD_LOG=$HOME/.cmdlog/.tcshlog/tcsh-history-`date +%Y-%m`.log
#alias precmd 'eval "if  `id -u` != 0 echo "[`date +%Y-%m-%d` `date +%T`] $USER@`hostname`:`pwd` `history 1`" >>! ${FULL_CMD_LOG}"'

#if ( "$(id -u)" -ne 0 ) alias (
#        FULL_CMD_LOG="$HOME/.cmdlog/.tcshlog/tcsh-history-$(date -u "+%Y-%m").log"
#	echo "$USER" >> ${FULL_CMD_LOG}
#	#echo "$USER@`hostname`:`pwd` [$(date -u)] `\history 1`" >> ${FULL_CMD_LOG}
#    endif
#}
#
###########################################################

#setenv PATH /exports/home/sabs/mirror/mpich2/bin:$PATH
#setenv LD_LIBRARY_PATH /exports/home/sabs/mirror/mpich2/lib

###########################################################
##
####       Gurobi
##
#############################################################
#
#setenv GUROBI_HOME /exports/home/scondon/gurobi702/linux64
#setenv PATH ${PATH}:${GUROBI_HOME}/bin
#setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:${GUROBI_HOME}/lib
#
#
#
###########################################################
##
####       Mosek Optimization Software
##
#############################################################
#
#setenv MSKHOME /exports/home/scondon/mosek
#setenv PATH ${PATH}:${MSKHOME}/mosek/8/tools/platform/linux64x86/bin/
#
#
#
