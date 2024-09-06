#/bin/bash

VERSION='v2.8'
RUN='run'
MASTER=$(pwd)/
INPUT="input.dat"
#JOB_QUEUE="pbs.script"
#JOB_QUEUE="lsf.script"

#since the dissociation constant for the reaction is very high
#reaction: SN + OH- <=> SI + H2O   Kd = 1.7*10^6
#there would be equal of the amounts of SI and TAA
#in this script we don't distingiush between SI and TAA 

#---------variables-------
i_sisn="-0.80";   	d_sisn="-0.10";		f_sisn=${i_sisn}
i_sitaa="-2.00";  	d_sitaa="-0.10";	f_sitaa=${i_sitaa}
i_teos="40";     	d_teos="1";     	f_teos=${i_teos}
i_taaoh="9";     	d_taaoh="1";   	        f_taaoh=${i_taaoh}
#i_h2o="9500";  	d_h2o="1";   	        f_h2o=${i_h2o}

manual_concen_flag="1"
#here we disregard the concentration of water and just enter TEOS and TAAOH molecule number.
#the amount of ionic silica produces is the same as TAAOH, cause we assume that the dissociation
#constant of the above reaction is too high

i_t="0.15";       	d_t="0.15";       	f_t=${i_t}
i_size="100";     	d_size="10";       	f_size=${i_size}

RUNS="3"

nsweeps="10000000";nprint="50000";nsnapshot="1000000";neqsweeps="5000000"
pen3="0.60"; pen4="0.30"
tstarf="0.6d0";m0="0.3";m1="0.35";m2="0.5"
start_config="0";time_limit="115.0"

#testing variables
i_size="8";     	d_size="10";       	f_size=${i_size}
nsweeps="200";nprint="50";nsnapshot="100";neqsweeps="300"
time_limit="0.0014"
#-------------------------

# make the directory structure
for sisn in $(seq ${i_sisn} ${d_sisn} ${f_sisn});do # SISNinteraction

    for sitaa in $(seq ${i_sitaa} ${d_sitaa} ${f_sitaa});do #SITAAinteraction
    
	for nteos in $(seq ${i_teos} ${d_teos} ${f_teos});do #SIconcentration
	    
	    for ntaaoh in $(seq ${i_taaoh} ${d_taaoh} ${f_taaoh});do #SNconcentration
	    		
		for t in $(seq ${i_t} ${d_t} ${f_t});do #temperature
		    
		    for cu in $(seq ${i_size} ${d_size} ${f_size});do #system size
			
			CURRENT_VERSION=${VERSION}_SISN${sisn}_SITAA${sitaa}_nSN${nteos}_nSI${ntaaoh}_nTAA${ntaaoh}_tstar${t}_${cu}cu
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
			    echo "${cu},${cu},${cu}" > ${INPUT}
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

		    done #system size
		done #temperature
	    done #TAA number
	done #SN number
    done # SI number
done #SITAA interaction