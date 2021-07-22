# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
export PATH=$PATH:/home/jbueno/flex/bin/
export PETSC_DIR=/home/jbueno/petsc/
#export PETSC_ARCH=optimO3native
export PETSC_ARCH=debug
export PETIGA_DIR=/home/jbueno/PetIGA/



alias cdpet='cd ~/petsc/'
alias cdiga='cd ~/PetIGA/demo'
alias cdigakit='cd ~/igakit/demo'
alias e='emacs'