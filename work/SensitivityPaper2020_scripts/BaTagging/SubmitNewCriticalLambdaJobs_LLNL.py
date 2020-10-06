#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity/work/SensitivityPaper2020_scripts/BaTagging/"
outputdir = "/p/lustre1/lenardo1/sensitivity_output/October5_2020_critical_lambda_ba_tagging_no_shape_error_histats_final_cuts/"
outputname = ""

base = "BaTagging_Sensitivity_Baseline2019_"


# Number of toy datasets to run for each hypothesis
num_datasets=20000



for num in range(0,1000):

	basename = base + str(num)
	
	scriptfilename = outputdir + "Sub/" +  base + str(num) + ".sub"
	os.system( "rm -f " + scriptfilename )
	outfilename = outputdir + "Out/" + base + str(num) + ".out"
	os.system( "rm -f " + outfilename )

	if num % 5 == 0:
		hyp_val = (float(num) + 0.0001)/8.
		iteration_num = 0
	elif num % 5 == 1:
		hyp_val = (float(num-1) + 0.0001)/8.
		iteration_num = 1
	elif num % 5 == 2:
		hyp_val = (float(num-2) + 0.0001)/8.
		iteration_num = 2
	elif num % 5 == 3:
		hyp_val = (float(num-3) + 0.0001)/8.
		iteration_num = 3
	elif num % 5 == 4:
		hyp_val = (float(num-4) + 0.0001)/8.
		iteration_num = 4
	else:
		hyp_val = (float(num-5) + 0.0001)/8.	
		iteration_num = 5
		
	thescript = "#!/bin/bash\n" + \
		"#SBATCH -t 02:00:00\n" + \
		"#SBATCH -A nuphys\n" + \
		"#SBATCH -e " + outfilename + "\n" + \
		"#SBATCH -o " + outfilename + "\n" + \
		"#SBATCH --mail-type=fail\n" + \
		"#SBATCH -J " + base + "\n" + \
		"#SBATCH --export=ALL \n" + \
		"source /usr/gapps/nexo/setup.sh \n" + \
		"source /g/g20/lenardo1/localpythonpackages/bin/activate \n" + \
		"cd " + execdir + "\n" + \
		"export STARTTIME=`date +%s`\n" + \
		"echo Start time $STARTTIME\n" + \
		"python3 ComputeCriticalLambdaForFixedNumSignal.py {} {} {} {}\n".format(iteration_num,hyp_val,num_datasets,outputdir) + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open( scriptfilename, 'w' )
	scriptfile.write( thescript )
	scriptfile.close()
	
	os.system( "sbatch " + scriptfilename )


