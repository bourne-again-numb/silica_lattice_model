#/bin/bash

VERSION='v2.8'
RUN='run'
MASTER=$(pwd)/
INPUT="input.dat"
#JOB_QUEUE="pbs.script"
JOB_QUEUE="lsf.script"

#since the dissociation constant for the reaction is very high
#reaction: SN + OH- <=> SI + H2O   Kd = 1.7*10^6
#there would be equal of the amounts of SI and TAA
#in this script we don't distingiush between SI and TAA 

#---------variables-------
i_sisn="-0.80";   	d_sisn="-0.10";		f_sisn=${i_sisn}
i_sitaa="-2.00";  	d_sitaa="-0.10";	f_sitaa=${i_sitaa}
i_teos="48";     	d_teos="5";     	f_teos=${i_teos}
i_taaoh="0";     	d_taaoh="4";   	        f_taaoh=${i_taaoh}
#i_h2o="9500";  	d_h2o="1";   	        f_h2o=${i_h2o}

manual_concen_flag="1"
#here we disregard the concentration of water and just enter TEOS and TAAOH molecule number.
#the amount of ionic silica produces is the same as TAAOH, cause we assume that the dissociation
#constant of the above reaction is too high

i_t="0.15";       	d_t="0.15";       	f_t=${i_t}
i_x="8";     	d_x="1";       	f_x=${i_x}
i_y="16";     	d_y="1";       	f_y=${i_y}
i_z="59";     	d_z="1";       	f_z=${i_z}

RUNS="3"

nsweeps="10000000";nprint="50000";nsnapshot="1000000";neqsweeps="5000000"
pen3="0.00"; pen4="0.00"
tstarf="0.6d0";m0="0.3";m1="0.35";m2="0.5"
start_config="0";time_limit="11.8"

#testing variables
#i_x="8";     	d_x="1";       	f_x=${i_x}
#i_y="8";     	d_y="1";       	f_y=${i_y}
#i_z="8";     	d_z="1";       	f_z=${i_z}

nsweeps="200";nprint="50";nsnapshot="100";neqsweeps="300"
time_limit="0.0014"
#-------------------------

# make the directory structure
for sisn in $(seq ${i_sisn} ${d_sisn} ${f_sisn});do # SISNinteraction

    for sitaa in $(seq ${i_sitaa} ${d_sitaa} ${f_sitaa});do #SITAAinteraction
    
	for ntaaoh in $(seq ${i_taaoh} ${d_taaoh} ${f_taaoh});do #SIconcentration
	    
	    for nteos in $(seq ${i_teos} ${d_teos} ${f_teos});do #SNconcentration
	    		
		for t in $(seq ${i_t} ${d_t} ${f_t});do #temperature
		    
		    for cux in $(seq ${i_x} ${d_x} ${f_x});do #x length
			for cuy in $(seq ${i_y} ${d_y} ${f_y});do #y length
			    for cuz in $(seq ${i_z} ${d_z} ${f_z});do #z length
			
				CURRENT_VERSION=${VERSION}_SISN${sisn}_SITAA${sitaa}_nSN${nteos}_nSI${ntaaoh}_nTAA${ntaaoh}_tstar${t}_${cux}x${cuy}x${cuz}
				mkdir -p ${CURRENT_VERSION}
				
			    #enter the current version directory
				cd ${CURRENT_VERSION}
				for j in $(seq 1 1 ${RUNS});do #run no.
				    
         		    # inside the current version directory create the run number directory
				    CURRENT_RUN=${RUN}_${j}
				    mkdir -p ${CURRENT_RUN}
				    
			    #enter the run directory
				    cd ${CURRENT_RUN}
				    
			     #once inside the current version directiory do the following
			     # copy the xcutable,*.mod,pbs file from the master directory
			     # create the input.dat 
			     # create the pbs.submit file
				    cp -f ${MASTER}xcutable .
				    cp -f ${MASTER}*.mod .
				    cp -f ${MASTER}${JOB_QUEUE} ./${JOB_QUEUE}
				    cp -f ${MASTER}*.gnu .
				    cp -f ${MASTER}analysis.sh .
				    
			    # -------------input.dat file creation---------------
				    echo "${cux},${cuy},${cuz}" > ${INPUT}
				    echo "${t}d0" >> ${INPUT}
				    echo "${nsweeps},${nprint},${nsnapshot},${neqsweeps}" >> ${INPUT}
				    echo "${nteos},${ntaaoh},${ntaaoh},${manual_concen_flag}" >> ${INPUT}
				    echo "${pen3}d0,${pen4}d0" >> ${INPUT}
				    echo "${sitaa}d0,-1.00d0,${sisn}d0" >> ${INPUT}
				    echo "${tstarf},${m0},${m1},${m2}" >> ${INPUT}
				    echo "${start_config}","${time_limit}" >> ${INPUT}
			    # ---------------------------------------------------
			    # -------------edit the job.submit file--------------
				    sed -i s/abcdef/"${CURRENT_VERSION}_${j}"/g ${MASTER}${CURRENT_VERSION}/${CURRENT_RUN}/${JOB_QUEUE}
			    # ---------------------------------------------------
				    
			    #submit the job
			    #qsub ${JOB_QUEUE}
			    #bsub < ${JOB_QUEUE}
			    #./xcutable
				
                            #exit the run directory
		       		    cd ${MASTER}/${CURRENT_VERSION}
				done #run no.			    
				cd ${MASTER}

			    done #z length
			done #y length
		    done #z length
		done #temperature
	    done #TAA number
	done #SN number
    done # SI number
done #SITAA interaction