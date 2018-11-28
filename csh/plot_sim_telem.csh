#!/bin/csh

# plotting the SORCE SIM telemetry 
# if no arguments are provided, plot telemetry from yesterday (24 hours worth)
# if one argument is provided, plot the last $arg[1] hours of data (ending now)
# if two arguments are provided, assume they are the start and stop times in SORCE days
#
# SBeland 2013/12/19

# IDL Setup
#setenv IDL_PATH "<IDL_DEFAULT>"
#setenv IDL_PATH ${IDL_PATH}:.:+/sorce/tools:+/data_systems/tools/astron:+/data_systems/tools/dfanning
#
# Tools Knapp Setup
#setenv bgkroot /data_systems/tools/knapp  
#setenv IDL_PATH ${IDL_PATH}:+${bgkroot}/idl 

# Set up the new query_database.
#set querydb_dir = /data_systems/tools/IdlJavaDbBridge
#setenv IDL_PATH ${IDL_PATH}:${querydb_dir} 

#set campaign_web_page_home = /sorce_data/internal_web/status/monitor/campaign_watcher
set campaign_web_page_home = /Users/sbeland/SORCE/data
echo "Executing: IDL SIM PLOT_TELEM "

#echo $IDL_PATH

if ($#argv == 0) then 
    set d1=`date +"%d"`
    @ d0 = $d1 - 1
    set t=`date +"-%b-%Y"`
    set t0="'$d0$t'"
    set t1="'$d1$t'"
    set time_format='/vms'
    set fname=`date +"%Y%m%d_%H%M"`
else if ($#argv == 1) then 
    set h1=`date +"%H"`
    @ h0 = $h1 - $argv[1]
    set t=`date +"%d-%b-%Y"`
    set t0="'$t $h0"':00:00'"'"
    set t1="'$t $h1"':00:00'"'"
    set time_format='/vms'
    set fname=`date +"%Y%m%d_%H%M"`
else if ($#argv == 2) then
    set t0=$argv[1]
    set t1=$argv[2]
    set time_format='/mission'
    set fname=`date +"%Y%m%d_%H%M"`
endif

set outfile_a="sima_telem_"$fname".png"
set idl_cmd_a="plot_telem,$t0,$t1,/sima,$time_format,/silent,psfile='$campaign_web_page_home/images/$outfile_a'"

set outfile_b="simb_telem_"$fname".png"
set idl_cmd_b="plot_telem,$t0,$t1,/simb,$time_format,/silent,psfile='$campaign_web_page_home/images/$outfile_b'"

echo "    running:"
echo "      $idl_cmd_a"
echo "      $idl_cmd_b"

idl  -quiet -e "device,decompose=0 & $idl_cmd_a &  $idl_cmd_b"
echo " "

/bin/cp -f $campaign_web_page_home/images/$outfile_a $campaign_web_page_home/sima_telem.png
/bin/cp -f $campaign_web_page_home/images/$outfile_b $campaign_web_page_home/simb_telem.png
touch $campaign_web_page_home/index.html

echo " "
echo "Monitor results at: http://dsweb:8080/sorcesds/status/monitor/campaign_watcher/"
echo "====JOB COMPLETE===="

