# .cshrc file with environmental variables and aliases for MSL, as well as setting up the terminal layout
alias mv 'mv -i'
alias cp 'cp -i'

set ignorecase
set smartcase

#TODO: CHANGE LINES 10, 17, AND 18 IN FILE IF MSLIB IS EVER MOVED TO A NEW DIRECTORY
# MSL ENVIRONMENTAL VARIABLES FOR INCLUSION OF EXTERNAL LIBRARIES
setenv MSL_DIR "/home/loiseau@ad.wisc.edu/github/mslib/trunk_AS"

setenv MSL_GSL T
setenv MSL_GLPK F
setenv MSL_BOOST F
setenv MSL_DEBUG F
setenv MSL_R F
setenv MSL_EXTERNAL_LIB_DIR /home/loiseau@ad.wisc.edu/github/mslib/ext_libs
setenv MSL_EXTERNAL_INCLUDE_DIR /home/loiseau@ad.wisc.edu/github/mslib/ext_includes

###########################################################
#
#       MAKE PROMPT "hostname path>"

# colored and truncated to 50 characters
set dollarone = '$1'
set dotdotdot = '"..."'
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