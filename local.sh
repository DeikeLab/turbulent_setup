#The executable name
EXE=windwave_turb
#Parameter value
RE_tau=720 #Default 40000
BO=200
MAXLEVEL=8
MINLEVEL=3
c_=0.5
ak=0.2
TIME=38.5001
uemaxRATIO=0.3
alterMU=16
end_sim=1000
do_eta_loc=1 
do_profile=1
do_fields=1
do_tagging=0
tout_eta=16.0
tout_slc=16.0
tout_pro=16.0
tout_res=2.0
tout_mov=4.0
prt_res=9
from_pr=0

mkdir ./profiles
mkdir ./eta
mkdir ./eta/eta_loc
mkdir ./field
mkdir ./images
mkdir ./restart_bk
mkdir ./tagging_f1
mkdir ./tagging_f2

./windwave_turb/$EXE $RE_tau $BO $MAXLEVEL $MINLEVEL $c_ $ak $TIME $uemaxRATIO $alterMU $end_sim $do_eta_loc $do_profile $do_fields $do_tagging $tout_eta $tout_slc $tout_pro $tout_res $tout_mov $prt_res $from_pr
