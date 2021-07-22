ls
emacs oa.sh.o5231
sh copy 
cd post/
ls
emacs step.dat
emacs database_3D.f
sh launch.sh 
cd ..
ls
cd cil_unity1/
ls
emacs oa.sh.o5230
sh copy 
cd post/
ls
emacs step.dat
emacs database_3D.f
sh launch.sh 
cd ..
cd fsi_bd/
ls
cd bd_unity0
ls
emacs oa.sh.o5216 
sh copy 
ls
cd post/
emacs step.dat
emacs database_3D.f
sh compile
./a.out 
cd ..
ls
cd ..
ls
cd bd_unity1/
ls
emacs oa.sh.o5217
sh copy 
cd post/
emacs step.dat
emacs database_3D.f
sh compile
./a.out 
ls
cd ..
ls
cd bd_unity2/
ls
emcas oa.sh.o5222 
cd ..
ls
ls -lsa
cd .bash
emacs .bash
cd bin/
ls
ls -lsa
cd ..
ls
emacs .bashrc 
cdiga
e 
cd petsc/
ls
cd ..
emacs .bashrc
cdiga
cd
ls
exit
cdiga
emacs .bashrc
exit
cdiga
ls
e Implosion
ls
cd
ls
cd fsi_bd/
ls
cd bd_unity2/
ls
e oa.sh.o5222 
sh copy 
cd post/
ls
cd ..
e oa.sh.o5222 
cd post/
e step.dat
ls
e database_3D.f
sh compile
./a.out 
exit
ls
cd fsi_bd/
ls
cd ..
ls
qstat
qdel 5394
ls
cd cube/
ls
emacs driver.f
sh compileVENUS
cd cube/
ls
sh copy 
rm restart.*
ls
ls -lsa output.out 
qstat
sh send_mpi.sh
qstat
ls
emacs oa.sh.o5407
ls
cd post/
sh del
cd ..
sh copy 
ls
cd post/
e step.dat
sh launch.sh 
ls
cd fsi_par
ls
cd 3thin0
ls
emacs oa.sh.o5005
emacs oa.sh
ls
cd fsi_par
cd fsi_cil
w
cd ..
cd fsi_cil
ls
emacs driver.f
ls
cd fsi_cil/
ls
cd cil_unity0/
ls
emacs oa.sh.o5231
ls
cd fsi_cil
cd cil_unity0/
ls
qstat
cd ..
cd cil_unity1/
ls
emacs oa.sh.o5230
qstat
emacs oa.sh.o5230
ls
qstat
cd fsi_par
cd real0/
emacs meshu.1.dat 
ls
cd fsi_cil/
ls
cd cil_unity1/
ls
qstat
emacs oa.sh.o5230
ls -lsa 
ls -lsa restart.*
qstat
emacs oa.sh.o5230
sh copy 
cd post/
e step.dat
sh launch.sh 
cd ..
cd cil_unity0/
ls
emacs oa.sh.o5231
sh copy 
cd post/
emacs step.dat
sh launch.sh 
cd ..
ls
cd ..
qstat
cd fsi_par/
ls
cd 3thin2/
ls
emacs oa.sh.o5160 
sh copy 
cd post/
e step.dat
sh launch.sh 
cd ..
qstat
qdel 5160
cd ..
cd fsi_bd/
ls
cd bd_unity0/
ls
emacs oa.sh.o5216 
cd post/
emacs step.dat
qdel 5216
qstat
cd ..
cd bd_unity1
cd ..
cd bd_unity1/
emacs oa.sh.o5217
cd post/
emacs step.dat
qstat
qdel 5217
ls
cd ..
ls
cd ..
ls
cd bd_unity2/
emacs oa.sh.o5222 
cd post/
emacs step.dat
cd ..
sh copy 
cd post/
emacs step.dat
ls
sh compile
./a.out 
ls
qstat
qstat -u "*"
ls
dk -f
df -k
ls
qstat -u "*"
ls
cd fsi_cil
ls
cd cil_unity1/
ls
ls -lsa restart.*
emacs oa.sh.o5230
qstat
emacs oa.sh.o5230
ls -lsa restart.4000.*
ls -lsa restart.3900.*
ls
cd ..
ls
cd cil_unity0/
ls
qstat
cd ..
cd --
cd fsi_bd
ls
cd bd_unity1/
ls
qstat
cd ..
cd bd_unity0
cd ..
ls
cd bd_unity2/
ls
cd ..
lks
ls
cd ..
ls
cd fsi_par
ls
cd ..
ls
cd cube/
ls
cd cube/
ls
qstat
qdel 5407
ls
rm restart.*
ls
rm restart.2*
rm restart.1*
rm restart.3*
rm restart.4*
rm restart.*
ls
cd ..
ls
dk -f
df -k
cd ..
ls
qstat
w
cd ..
ls
cd jbueno/
ls
w
exit
ls
cd fsi_par
ls
emacs e3Rhs_3D.f
ls
exit
ls
conex -f
ls -lsa
emacs .Xauthority 
emacs .bashrc
ls
rm .hg/
cd .hg/
ls
cd ..
rm -r .hg/
ls
ls -lsa 
cd .pki/
ls
cd nssdb/
ls
ls -lsa
cd ..
ls
ls -lsa
cd ..
ls
emacs .bash_history 
emacs .bash_profile 
ls
cd PetIGA
ls
make all
cd demo/
ls
sh del.sh 
ls
rm out*.txt
ls
rm *~
ls
mkdir implosion
emacs Implosion.c 
emacs send_mpi 
emacs oa.sh 
qstat
make Implosion
cd implosion/
ls
cd ..
sh send_mpi
qstat
ls
emacs oa.sh.o5483 
emacs send_mpi
emacs oa.sh
sh send_mpi
qstat
emacs oa.sh.o5484
ls
emacs out.txt 
emacs oa.sh.o5484
qstat
qdel 5484
ls
emacs oa.sh
qstat
sh send_mpi
qstat
cd implosion/
ls
emacs out.txt 
cd ..
emacs oa.sh.o5485
cd implosion/
emacs out.txt 
ls
emacs out.txt 
ls
emacs out.txt 
ls
ls -lsa
cd ..
ls
cd implosion/
ls
ls -lsa 
ls
qstat
ls -lsa
cd ./ssh
cd .ssh/
ls
emacs authorized_keys 
ls
emacs id_rsa
emacs known_hosts 
conx -f
conex -f
ls
mkdir copy
ls
cp -r .ssh/ ./copy/
cd copy/
ls
ls -lsa 
cd ./ssh
cd .ssh
ls
cd ..
ls
cd ..
ls
ls -lsa
cd .ssh/
ls
cd ..
ls
cd PetIGA
ls
conex -f
ls
 cd bin/
ls
emacs conex 
ls
emacs jb\@ubuntu 
ls -lsa
cd ..
ls
ls -lsa 
rm .ssh/
rm -r .ssh/
ls
conex -f
conex -v
conex -f
ls
ls -lsa 
conex -v
conex -f
ls
ls -lsa 
cd .ssh/
ls
emacs authoized_keys 
emacs authorized_keys
emacs id_rsa
emacs id_rsa.pub 
ls
emacs known_hosts 
emacs authorized_keys
cd ..
ls
wget http://prdownloads.sourceforge.net/flex/flex-2.5.33.tar.gz?download
ls
tar -xvzf flex-2.5.33.tar.gz
ls
cd flex-2.5.33
ls
cd ..
ls
rm -r flex-2.5.33
rm -r flex-2.5.33.tar.gz 
ls
cd PetIGA
ls
cd ..
lsls
ls
cd PetIGA_acel/
ls
cd demo/
ls
cp oa.sh send_mpi ../../PetIGA/demo/
ls
exit
ls
ls -lsa
cd .ssh/
ls
emacs id_rsa
emacs id_rsa.pub
emacs authorized_keys
emacs authoized_keys 
rm authoized_keys 
ls
exit
ls
cd copy/
ls
ls -lsa
cd ..
ls
rm -r .ssh/
ls
ls -lsa
ls
ls -lsa
mkdir .ssh
cp copy/.ssh/* ./.ssh/
cd .ssh/
ls
exit
ls
ls -lsa
emacs .Xauthority 
emacs .kshrc 
cd .gnome2/
ls
cd keyrings/
ls
emacs login.keyring 
cd ..
cd .mozilla/
ls
emacs extensions/
cd extensions/
ls
ls -lsa
cd ..
cd plugins/
ls
ls -lsa 
ls .
ls ./
ls
cd ..
ls
ls -lsa
ssh 
ssh -F
ls
cd .ssh/
emacs authorized_keys 
exit
ls
ssh-keygen -t dsa
ls
exit
ls
qstat
w
cd ..
ls
qstat -u "*"
ls
cd 
ls
cd PetIGA_acel/
ls
cd demo/
ls
cd ..
ls
cd ..
cd PetIGA
ls
cd demo/
ls
cd ..
ls
cd PetIGA_
cd PetIGA_acel/
ls
cd demo/
ls
cd ..
ls
rm -r PetIGA
ls
cd
ls
mkdir PetIGA
mv -r ./arch-linux2-c-debug ./PetIGA
mv ./arch-linux2-c-debug ./PetIGA
ls
cd PetIGA
ls
cd ..
mv ./conf/ ./PetIGA
mv ./demo/ ./PetIGA
mv ./flex/ ./PetIGA
mv ./include/ ./PetIGA
mv makefile ./PetIGA
rm out*
mv README.rst ./PetIGA
mv ./test/ ./PetIGA
mv CMakeLists.txt  ./PetIGA
mv ./docs/ ./PetIGA
mv LICENSE.rst ./PetIGA
mv ./src/ ./PetIGA
ls
cd PetIGA
ls
cd demo/
ls
cd ..
ls
cd ..
ls
rm -r PetIGA
ls
mkdir PetIGA
ls
cd .ssh/
ls
e id_rsa
/etc/init.d/ssh restart
cd
ls
chmod 700 ~/.ssh 
ls
exit
ls
cd .ssh/
ls
cd
emacs /etc/ssh/sshd_config 
/etc/ssh/sshd_config 
ls
cd .ssh/
ls
emacs authorized_keys
exit
wget http://prdownloads.sourceforge.net/flex/flex-2.5.33.tar.gz?download
ls
tar -xvzf flex-2.5.33.tar.gz
cd flex-2.5.33
ls
./configure --prefi/flex
cd ..
ls
./configure --prefix=~/flex
./configure --prefix=/home/jbueno/flex
ls
cd flex-2.5.33
./configure --prefix=/home/jbueno/flex
make
make install-sh 
make install
ls
cd ..
ls
rm flex-2.5.33.tar.gz 
ls
exit
conex -f
exit
ls
emacs .bashrc
emacs /etc/ssh/ssh_config 
ls
ls -lsa
emacs /etc/ssh/ssh_known_hosts 
emacs /etc/ssh/ssh_host_rsa_key
emacs /etc/ssh/ssh_host_key.pub 
emacs /etc/ssh/ssh_host_key
emacs /etc/ssh/ssh_host_rsa_key
emacs /etc/ssh/ssh_host_dsa_key
emacs ssh_config
cd ..
ls
cd ..
ls
cd etc/
ls
cd ..
ls
cd
ls
w
ls
rm -r .ssh/
mkdir .ssh
cp copy/.ssh/* .ssh/
ls
cd .ssh/
ls
ls -lsa 
chmod 777 *
ls -lsa 
cd ..
chmod 777 .ssh/
ls -lsa 
exit
conex -f
ls
chmod 700 .ssh/
conex -f
cd .ssh/
ls
chmod 700 *
conex -f
cp * ../
ls
cd ..
ls
exit
cd .ssh/
ls
cd ..
cp copy/.ssh/* ./.ssh/
cd .ssh/
ls
e config
ls -lsa
chmod 700 *
ls
ls -lsa
ssh-agent .ssh/id_rsa
ssh-agent id_rsa
ssh-agent ~/.ssh/id_rsa
ls
cd ..
ls
cd .ssh/
ls
rm config 
ls
conex -f
ls
rm ./.ssh/
rm -r ./.ssh/
ls
ls -lsa
ssh jbuenofenix.udc.es
ssh jbueno@fenix.udc.es
ls -lsa
cd .ssh/
ls
emacs known_hosts 
conex -f
ls
emacs authorized_keys 
rm autho*
rm id_rsa*
rm known_hosts 
ls
cd .ssh/
ls
rm authoized_keys 
ls
emacs authorized_keys 
emcas id_rsa
e id_rsa
e id_rsa.pub 
e known_hosts 
conex -f
pwd
ls
e /var/log/auth.log
tail -f /var/log/auth.log
ls -lsa
ls -lsa .*
ls -lsa
ls -lsa
cd PetIGA
ls
cd demo/
ls
cd implosion/
ls
ssh -i linceiberico jbueno@fenix.udc.es


cp posp_vel.py ./implosion/
ls
mkdir implosion1
cp oa.sh send_mpi ./implosion1/
emacs oa.sh
emacs send_mpi
cd ..
cd src/
e petigats2.c 
e petigaelem.c 
cd ..
make all
cdiga
make Implosion
cp implosion ./implosion1/
cd implosion1/
ls
cd ..
ls -lsa
cp Implosion ./implosion1
cd implosion1/
ls
emacs oa.sh 
cd ..
emacs oa.sh
sh send_mpi
qstat
cd implosion1/
ls
cd ..
ls
cd implosion1/
ls
emacs out.txt 
cd ..
e oa.sh.o5486
emacs oa.sh
qstat
sh send_mpi
cd implosion1/
ls
emacs out.txt 
ls
emacs out.txt 
cd ..
ls
emacs oa.sh.o5487
emacs Implosion.c
make Implosion
cp Implosion ./implosion1/
cd implosion1/
ls -lsa out.txt 
ls -lsa Implosion 
ls
rm oa.sh
rm send_mpi 
cd ..
ls
sh send_mpi
qstat
cd implosion1/
ls
emacs out.txt 
e out.txt 
cd ..
ls
emacs out.txt 
ls -lsa out.txt 
emacs oa.sh.o5488
cd implosion1/
ls
emacs out.txt 
ls -lsa out.txt 
ls
cd ..
ls
mkdir Geo
cp *.dat ./Geo/
cd Geo/
ls
cd ..
emacs send_mpi
emacs oa.sh
ls
ls -lsa oa.sh.o5488
ls
emacs makefile 
cd implosion1/
ls
cd ..
make implosion2
mkdir implosion2
ls
emacs makefile
make Implosion
emacs makefile
make Implosion
emacs makefile
make Implosion
e makefile
make Implosion
cd implosion2/
ls
cd ..
ls
ls -lsa Implosion
emacs makefile
make Implosion
cd implosion2/
ls
cd ..
ls -lsa Implosion
ls
emacs makefile
make Implosion
e makefile
e compile
sh compile 
cd implosion2/
ls
cd ..
ls
emacs oa.sh
cd implosion1/
ls
cd ..
ls
psswd
passwd
ls
cd PetIGA
ls
cd src/
ls
e petigaelem.c
cd ..
make all
cdiga
ls
e oa.sh
e compile
make Implosion
cd implosion2/
ls -lsa ou
ls -lsa Implosion 
ls
rm Implosion 
cd ..
sh compile 
cd implosion2/
ls
ls -lsa ou
ls -lsa Implosion 
cd ..
sh compile 
cd implosion2/
ls -lsa Implosion 
cd ..
ls
emacs send_mpi
w
sh send_mpi
qstat
ls
cd implosion2/
ls
emacs out.txt 
e out.txt 
qstat -u*
qstat -u *
qstat -u "*"
w
exit
w
ls
exit
ls
cd PetIGA
ls
cd demo
ls
cd implosion2/
ls
qstat
emacs out.txt 
ls
cd ..
ls
emacs oa.sh
e Implosion.c
cd implosion2/
e out.txt 
cd ..
e Implosion.c
sh compile 
cd implosion2/
cd ..
ls
cd implosion2
ls -lsa Implosion 
cd ..
cd implosion2/
ls
rm disp*.dat
ls
rm iga*.dat
ls
rm vel*.dat
ls
cd ..
sh send_mpi
qstat
cd implosion2/
ls
e out.txt 
