#!/bin/sh
set -x
date

HOMEhafs=${HOMEhafs:-/lfs/h2/emc/hur/noscrub/${USER}/save/hafs_202502}
source ${HOMEhafs}/ush/hafs_pre_job.sh.inc

cd ${HOMEhafs}/rocoto
EXPT=$(basename ${HOMEhafs})
#opts="-t -s sites/${WHERE_AM_I:-wcoss2}.ent -f"
opts="-t -f"
#===============================================================================

 # HAFSv2.1B phase3 configuration
 confopts="config.EXPT=${EXPT} config.SUBEXPT=${EXPT}_v2p1b_phase3 ../parm/hafsv2p1b_phase3.conf \
     dir.COMrtofs=/lfs/h2/emc/hur/noscrub/hafs-input/COMRTOFSv2.5_new" # wcoss2
#    dir.COMrtofs=/work/noaa/hwrf/noscrub/hafs-input/COMRTOFSv2.5_new" # orion/hercules
#    dir.COMrtofs=/lfs5/HFIP/hwrf-data/hafs-input/COMRTOFSv2.5_new" # jet
#    dir.COMrtofs=/scratch1/NCEPDEV/hwrf/noscrub/hafs-input/COMRTOFSv2.5_new" # hera
#    dir.COMrtofs=/gpfs/f6/drsa-hurr1/world-shared/noscrub/hafs-input/COMRTOFSv2.5_new" # gaeaC6

 # Technical testing
#./run_hafs.py ${opts} 2024092506-2024092512 09L HISTORY ${confopts} # Helene

#===============================================================================
 # 2024 NATL Storms
#./run_hafs.py ${opts} 2024061712-2024062012 01L HISTORY ${confopts} # Alberto
#./run_hafs.py ${opts} 2024062800-2024070818 02L HISTORY ${confopts} # Beryl
#./run_hafs.py ${opts} 2024062812-2024070100 03L HISTORY ${confopts} # Chris
#./run_hafs.py ${opts} 2024080218-2024080812 04L HISTORY ${confopts} # Debby
#./run_hafs.py ${opts} 2024081118-2024082006 05L HISTORY ${confopts} # Ernesto
#./run_hafs.py ${opts} 2024090806-2024091200 06L HISTORY ${confopts} # Francine
#./run_hafs.py ${opts} 2024091106-2024091806 07L HISTORY ${confopts} # Gordon
#./run_hafs.py ${opts} 2024091500-2024091612 08L HISTORY ${confopts} # Eight
#./run_hafs.py ${opts} 2024092312-2024092712 09L HISTORY ${confopts} # Helene
#./run_hafs.py ${opts} 2024092506-2024093006 10L HISTORY ${confopts} # Isaac
#./run_hafs.py ${opts} 2024092612-2024100100 11L HISTORY ${confopts} # Joyce
#./run_hafs.py ${opts} 2024092906-2024100706 12L HISTORY ${confopts} # Kirk
#./run_hafs.py ${opts} 2024100112-2024101212 13L HISTORY ${confopts} # Leslie
#./run_hafs.py ${opts} 2024100500-2024101012 14L HISTORY ${confopts} # Milton
#./run_hafs.py ${opts} 2024101812-2024102006 15L HISTORY ${confopts} # Nadine
#./run_hafs.py ${opts} 2024101900-2024102212 16L HISTORY ${confopts} # Oscar
#./run_hafs.py ${opts} 2024110118-2024110406 17L HISTORY ${confopts} # Patty
#./run_hafs.py ${opts} 2024110312-2024111012 18L HISTORY ${confopts} # Rafael
#./run_hafs.py ${opts} 2024111318-2024111800 19L HISTORY ${confopts} # Sara

 # 2023 NATL Storms
#./run_hafs.py ${opts} 2023060106-2023060312 02L HISTORY ${confopts} # Arlene
#./run_hafs.py ${opts} 2023061900-2023062412 03L HISTORY ${confopts} # Bret
#./run_hafs.py ${opts} 2023062200-2023062600 04L HISTORY ${confopts} # Cindy
#./run_hafs.py ${opts} 2023071212-2023071606 05L HISTORY ${confopts} # Don part 1
#./run_hafs.py ${opts} 2023071618-2023072406 05L HISTORY ${confopts} # Don part 2
#./run_hafs.py ${opts} 2023081818-2023082212 06L HISTORY ${confopts} # Gert part 1
#./run_hafs.py ${opts} 2023083100-2023090412 06L HISTORY ${confopts} # Gert part 2
#./run_hafs.py ${opts} 2023081818-2023082112 07L HISTORY ${confopts} # Emily part 1
#./run_hafs.py ${opts} 2023082212-2023082512 07L HISTORY ${confopts} # Emily part 2
#./run_hafs.py ${opts} 2023082000-2023090112 08L HISTORY ${confopts} # Franklin
#./run_hafs.py ${opts} 2023082106-2023082212 09L HISTORY ${confopts} # Harold
#./run_hafs.py ${opts} 2023082618-2023090218 10L HISTORY ${confopts} # Idalia
#./run_hafs.py ${opts} 2023082900-2023090118 11L HISTORY ${confopts} # Jose
#./run_hafs.py ${opts} 2023083112-2023090412 12L HISTORY ${confopts} # Katia
#./run_hafs.py ${opts} 2023090500-2023091606 13L HISTORY ${confopts} # Lee part 1
#./run_hafs.py ${opts} 2023091618-2023091712 13L HISTORY ${confopts} # Lee part 2
#./run_hafs.py ${opts} 2023090700-2023091712 14L HISTORY ${confopts} # Margot
#./run_hafs.py ${opts} 2023091418-2023092206 15L HISTORY ${confopts} # Nigel
#./run_hafs.py ${opts} 2023092112-2023092318 16L HISTORY ${confopts} # Ophelia
#./run_hafs.py ${opts} 2023092300-2023100612 17L HISTORY ${confopts} # Philippe
#./run_hafs.py ${opts} 2023092606-2023100200 18L HISTORY ${confopts} # Rina
#./run_hafs.py ${opts} 2023101012-2023101518 19L HISTORY ${confopts} # Sean
#./run_hafs.py ${opts} 2023101712-2023102900 20L HISTORY ${confopts} # Tammy
#./run_hafs.py ${opts} 2023102318-2023102406 21L HISTORY ${confopts} # Twenty-on
#./run_hafs.py ${opts} 2023111618-2023111718 22L HISTORY ${confopts} # Twenty-tw

 # 2022 NATL Storms
#./run_hafs.py ${opts} 2022060118-2022060618 01L HISTORY ${confopts} # Alex
#./run_hafs.py ${opts} 2022062612-2022070206 02L HISTORY ${confopts} # Bonnie
#./run_hafs.py ${opts} 2022070118-2022070300 03L HISTORY ${confopts} # Colin
#./run_hafs.py ${opts} 2022082000-2022082006 04L HISTORY ${confopts} # Four
#./run_hafs.py ${opts} 2022090100-2022090812 05L HISTORY ${confopts} # Danielle
#./run_hafs.py ${opts} 2022083012-2022091018 06L HISTORY ${confopts} # Earl
#./run_hafs.py ${opts} 2022091406-2022092418 07L HISTORY ${confopts} # Fiona
#./run_hafs.py ${opts} 2022092006-2022092518 08L HISTORY ${confopts} # Gaston
#./run_hafs.py ${opts} 2022092012-2022100106 09L HISTORY ${confopts} # Ian
#./run_hafs.py ${opts} 2022092306-2022092500 10L HISTORY ${confopts} # Hermine
#./run_hafs.py ${opts} 2022092718-2022092912 11L HISTORY ${confopts} # Eleven
#./run_hafs.py ${opts} 2022100412-2022100618 12L HISTORY ${confopts} # Twelve
#./run_hafs.py ${opts} 2022100412-2022100912 13L HISTORY ${confopts} # Julia
#./run_hafs.py ${opts} 2022101106-2022101506 14L HISTORY ${confopts} # Karl
#./run_hafs.py ${opts} 2022103006-2022110506 15L HISTORY ${confopts} # Lisa
#./run_hafs.py ${opts} 2022103118-2022110312 16L HISTORY ${confopts} # Martin
#./run_hafs.py ${opts} 2022110612-2022111100 17L HISTORY ${confopts} # Nicole

#===============================================================================
 # 2024 EPAC storms
#./run_hafs.py ${opts} 2024070412-2024070518 01E HISTORY ${confopts} # Aletta
#./run_hafs.py ${opts} 2024072418-2024072618 02E HISTORY ${confopts} # Bud
#./run_hafs.py ${opts} 2024073112-2024080600 03E HISTORY ${confopts} # Carlotta
#./run_hafs.py ${opts} 2024080312-2024080518 04E HISTORY ${confopts} # Daniel
#./run_hafs.py ${opts} 2024080412-2024080818 05E HISTORY ${confopts} # Emilia
#./run_hafs.py ${opts} 2024080518-2024080712 06E HISTORY ${confopts} # Fabio
#./run_hafs.py ${opts} 2024081812-2024083000 07E HISTORY ${confopts} # Gilma
#./run_hafs.py ${opts} 2024082518-2024082906 08E HISTORY ${confopts} # Hector
#./run_hafs.py ${opts} 2024091212-2024091506 09E HISTORY ${confopts} # Ileana
#./run_hafs.py ${opts} 2024092218-2024092412 10E HISTORY ${confopts} # John part 1
#./run_hafs.py ${opts} 2024092500-2024092718 10E HISTORY ${confopts} # John part 2
#./run_hafs.py ${opts} 2024100118-2024100312 11E HISTORY ${confopts} # Eleven
#./run_hafs.py ${opts} 2024102118-2024102706 12E HISTORY ${confopts} # Kristy
#./run_hafs.py ${opts} 2024110118-2024110306 13E HISTORY ${confopts} # Lane
#./run_hafs.py ${opts} 2024110612-2024110712 14E HISTORY ${confopts} # Fourteen

 # 2023 EPAC storms
#./run_hafs.py ${opts} 2023062718-2023070212 01E HISTORY ${confopts} # Adrian
#./run_hafs.py ${opts} 2023062900-2023070112 02E HISTORY ${confopts} # Beatriz
#./run_hafs.py ${opts} 2023071118-2023071912 03E HISTORY ${confopts} # Calvin
#./run_hafs.py ${opts} 2023072106-2023072206 04E HISTORY ${confopts} # Four
#./run_hafs.py ${opts} 2023073118-2023081300 05E HISTORY ${confopts} # Dora
#./run_hafs.py ${opts} 2023080512-2023080712 06E HISTORY ${confopts} # Eugene
#./run_hafs.py ${opts} 2023081218-2023081706 07E HISTORY ${confopts} # Fernanda
#./run_hafs.py ${opts} 2023081400-2023081718 08E HISTORY ${confopts} # Greg
#./run_hafs.py ${opts} 2023081612-2023082018 09E HISTORY ${confopts} # Hilary
#./run_hafs.py ${opts} 2023082700-2023082912 10E HISTORY ${confopts} # Irwin
#./run_hafs.py ${opts} 2023090418-2023091018 11E HISTORY ${confopts} # Jova
#./run_hafs.py ${opts} 2023091518-2023091618 12E HISTORY ${confopts} # Twelve
#./run_hafs.py ${opts} 2023091912-2023092212 13E HISTORY ${confopts} # Kenneth
#./run_hafs.py ${opts} 2023092318-2023092418 14E HISTORY ${confopts} # Fourteen
#./run_hafs.py ${opts} 2023100306-2023101100 15E HISTORY ${confopts} # Lidia
#./run_hafs.py ${opts} 2023100800-2023101000 16E HISTORY ${confopts} # Max
#./run_hafs.py ${opts} 2023101718-2023102306 17E HISTORY ${confopts} # Norma
#./run_hafs.py ${opts} 2023102212-2023102512 18E HISTORY ${confopts} # Otis
#./run_hafs.py ${opts} 2023102818-2023110600 19E HISTORY ${confopts} # Pilar
#./run_hafs.py ${opts} 2023112312-2023112612 20E HISTORY ${confopts} # Ramon

 # 2022 EPAC storms
#./run_hafs.py ${opts} 2022052800-2022053106 01E HISTORY ${confopts} # Agatha
#./run_hafs.py ${opts} 2022061412-2022062012 02E HISTORY ${confopts} # Blas
#./run_hafs.py ${opts} 2022061618-2022062812 03E HISTORY ${confopts} # Celia
#./run_hafs.py ${opts} 2022071212-2022071912 04E HISTORY ${confopts} # Bonnie
#./run_hafs.py ${opts} 2022070918-2022071618 05E HISTORY ${confopts} # Darby
#./run_hafs.py ${opts} 2022071512-2022072112 06E HISTORY ${confopts} # Estelle
#./run_hafs.py ${opts} 2022072606-2022080212 07E HISTORY ${confopts} # Frank
#./run_hafs.py ${opts} 2022072712-2022080318 08E HISTORY ${confopts} # Georgette
#./run_hafs.py ${opts} 2022080612-2022081018 09E HISTORY ${confopts} # Howard
#./run_hafs.py ${opts} 2022081318-2022081618 10E HISTORY ${confopts} # Ivette p1
#./run_hafs.py ${opts} 2022082012-2022082100 10E HISTORY ${confopts} # Ivette p2
#./run_hafs.py ${opts} 2022082112-2022082118 10E HISTORY ${confopts} # Ivette p3
#./run_hafs.py ${opts} 2022090118-2022090400 11E HISTORY ${confopts} # Javier
#./run_hafs.py ${opts} 2022090412-2022090918 12E HISTORY ${confopts} # Kay
#./run_hafs.py ${opts} 2022091518-2022091712 13E HISTORY ${confopts} # Lester
#./run_hafs.py ${opts} 2022091718-2022092012 14E HISTORY ${confopts} # Madeline
#./run_hafs.py ${opts} 2022092118-2022092518 15E HISTORY ${confopts} # Newton p1
#./run_hafs.py ${opts} 2022092718-2022092818 15E HISTORY ${confopts} # Newton p2
#./run_hafs.py ${opts} 2022092900-2022100318 16E HISTORY ${confopts} # Orlene
#./run_hafs.py ${opts} 2022100318-2022100512 17E HISTORY ${confopts} # Paine
#./run_hafs.py ${opts} 2022100918-2022101012 18E HISTORY ${confopts} # Julia
#./run_hafs.py ${opts} 2022102000-2022102318 19E HISTORY ${confopts} # Roslyn

#===============================================================================
 # 2023 WPAC Storms
#./run_hafs.py ${opts} 2023042012-2023042206 01W HISTORY ${confopts} # Sanvu
#./run_hafs.py ${opts} 2023052012-2023060306 02W HISTORY ${confopts} # Mawar
#./run_hafs.py ${opts} 2023060600-2023061206 03W HISTORY ${confopts} # Guchol
#./run_hafs.py ${opts} 2023071500-2023071806 04W HISTORY ${confopts} # Talim
#./run_hafs.py ${opts} 2023072112-2023072806 05W HISTORY ${confopts} # Doksuri
#./run_hafs.py ${opts} 2023072706-2023081018 06W HISTORY ${confopts} # Khanun
#./run_hafs.py ${opts} 2023080718-2023081712 07W HISTORY ${confopts} # Lan
#./run_hafs.py ${opts} 2023082306-2023082806 08W HISTORY ${confopts} # Damrey
#./run_hafs.py ${opts} 2023082318-2023090312 09W HISTORY ${confopts} # Saola
#./run_hafs.py ${opts} 2023082812-2023090418 10W HISTORY ${confopts} # Haikui
#./run_hafs.py ${opts} 2023083000-2023090400 11W HISTORY ${confopts} # Kirogi
#./run_hafs.py ${opts} 2023090518-2023090818 12W HISTORY ${confopts} # Yun-Yeung
#./run_hafs.py ${opts} 2023092506-2023092518 13W HISTORY ${confopts} # Thirteen
#./run_hafs.py ${opts} 2023092912-2023101000 14W HISTORY ${confopts} # Koinu
#./run_hafs.py ${opts} 2023100700-2023101412 15W HISTORY ${confopts} # Bolaven
#./run_hafs.py ${opts} 2023101800-2023102100 16W HISTORY ${confopts} # Sanba
#./run_hafs.py ${opts} 2023111212-2023111512 17W HISTORY ${confopts} # Seventeen
#./run_hafs.py ${opts} 2023121718-2023121812 18W HISTORY ${confopts} # Jelawat

 # 2022 WPAC Storms
#./run_hafs.py ${opts} 2022033018-2022033100 01W HISTORY ${confopts} # One
#./run_hafs.py ${opts} 2022040700-2022041506 02W HISTORY ${confopts} # Malakas
#./run_hafs.py ${opts} 2022040906-2022041218 03W HISTORY ${confopts} # Megi
#./run_hafs.py ${opts} 2022062918-2022070212 04W HISTORY ${confopts} # Chaba
#./run_hafs.py ${opts} 2022063018-2022070900 05W HISTORY ${confopts} # Aere
#./run_hafs.py ${opts} 2022072900-2022080100 06W HISTORY ${confopts} # Songda
#./run_hafs.py ${opts} 2022080100-2022080106 07W HISTORY ${confopts} # Trases
#./run_hafs.py ${opts} 2022080400-2022080406 08W HISTORY ${confopts} # Eight
#./run_hafs.py ${opts} 2022081106-2022081412 09W HISTORY ${confopts} # Meari
#./run_hafs.py ${opts} 2022082206-2022082512 10W HISTORY ${confopts} # Ma-on
#./run_hafs.py ${opts} 2022082206-2022082512 11W HISTORY ${confopts} # Tokage
#./run_hafs.py ${opts} 2022082806-2022090606 12W HISTORY ${confopts} # Hinnamnor
#./run_hafs.py ${opts} 2022083018-2022090100 13W HISTORY ${confopts} # Thirteen
#./run_hafs.py ${opts} 2022090618-2022091600 14W HISTORY ${confopts} # Muifa
#./run_hafs.py ${opts} 2022091018-2022091500 15W HISTORY ${confopts} # Merbok
#./run_hafs.py ${opts} 2022091218-2022091918 16W HISTORY ${confopts} # Nanmado
#./run_hafs.py ${opts} 2022092118-2022092400 17W HISTORY ${confopts} # Talas
#./run_hafs.py ${opts} 2022092200-2022092806 18W HISTORY ${confopts} # Noru
#./run_hafs.py ${opts} 2022092518-2022092900 19W HISTORY ${confopts} # Kulap
#./run_hafs.py ${opts} 2022092806-2022100500 20W HISTORY ${confopts} # Roke
#./run_hafs.py ${opts} 2022101212-2022101506 21W HISTORY ${confopts} # Twenty-on
#./run_hafs.py ${opts} 2022101400-2022101500 22W HISTORY ${confopts} # Sonca
#./run_hafs.py ${opts} 2022101418-2022102012 23W HISTORY ${confopts} # Nesat
#./run_hafs.py ${opts} 2022101800-2022101906 24W HISTORY ${confopts} # Haitang
#./run_hafs.py ${opts} 2022102000-2022102300 25W HISTORY ${confopts} # Twenty-fi
#./run_hafs.py ${opts} 2022102700-2022110312 26W HISTORY ${confopts} # Nalgae
#./run_hafs.py ${opts} 2022103018-2022103118 27W HISTORY ${confopts} # Banyan
#./run_hafs.py ${opts} 2022111212-2022111418 28W HISTORY ${confopts} # Yamaneko
#./run_hafs.py ${opts} 2022121018-2022121218 29W HISTORY ${confopts} # Pakhar

#===============================================================================
 # 2023 NIO storms
#./run_hafs.py ${opts} 2023051018-2023051418 01B HISTORY ${confopts} # Mocha
#./run_hafs.py ${opts} 2023060600-2023061612 02A HISTORY ${confopts} # Biparjoy
#./run_hafs.py ${opts} 2023060912-2023061000 03B HISTORY ${confopts} # Three
#./run_hafs.py ${opts} 2023073112-2023080112 04B HISTORY ${confopts} # Four
#./run_hafs.py ${opts} 2023102006-2023102406 05A HISTORY ${confopts} # Tej
#./run_hafs.py ${opts} 2023102306-2023102500 06B HISTORY ${confopts} # Hamoon
#./run_hafs.py ${opts} 2023111618-2023111712 07B HISTORY ${confopts} # Midhili
#./run_hafs.py ${opts} 2023120300-2023120418 08B HISTORY ${confopts} # Michaung

 # 2022 NIO storms
#./run_hafs.py ${opts} 2022030418-2022030612 01B HISTORY ${confopts} # One
#./run_hafs.py ${opts} 2022050706-2022051112 02B HISTORY ${confopts} # Asani
#./run_hafs.py ${opts} 2022081206-2022081312 03A HISTORY ${confopts} # Three
#./run_hafs.py ${opts} 2022081818-2022082000 04B HISTORY ${confopts} # Four
#./run_hafs.py ${opts} 2022102312-2022102418 05B HISTORY ${confopts} # Sitrang
#./run_hafs.py ${opts} 2022120800-2022121000 06B HISTORY ${confopts} # Mandous
#./run_hafs.py ${opts} 2022121500-2022121718 07A HISTORY ${confopts} # Seven

#===============================================================================

date

echo 'cronjob done'
