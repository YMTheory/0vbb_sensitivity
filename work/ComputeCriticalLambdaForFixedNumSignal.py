# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os
sys.path.append('../modules')


######################################################################
# Check arguments and load inputs
if len(sys.argv) == 4:
	iteration_num = int(sys.argv[1])
	input_num_signal = float(sys.argv[2])
	output_dir = sys.argv[3]
	if not os.path.exists(output_dir):
		sys.exit('\nERROR: path to output_dir does not exist\n')
else:
	print('\n\nERROR: ComputeCriticalLambdaForNumSignal.py requires 3 arguments')
	print('Usage:')
	print('\tpython ComputeCriticalLambdaForNumSignal.py ' + \
		'<iteration_num> <input_num_signal> </path/to/output/directory/>')
	sys.exit('\n')
######################################################################



# Import useful libraries for analysis
import pandas as pd
import histlite as hl
import numpy as np
from matplotlib import pyplot as plt
import copy

# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood

# Import iMinuit
from iminuit import Minuit


# Create the workspace
workspace = nEXOFitWorkspace.nEXOFitWorkspace()
workspace.LoadComponentsTableFromFile('../tables/ComponentsTable_D-005_v25_2020-01-21.h5')
workspace.CreateGroupedPDFs()

# Create the likelihood object
likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs)

# Get the initial values; set the BB0n num_signal to the user-provided input
initial_values = likelihood.GetVariableValues()
initial_values[ likelihood.GetVariableIndex('Bb0n') ] = input_num_signal

# Initialize the model
likelihood.model_obj.UpdateVariables(initial_values)
likelihood.model_obj.GenerateModelDistribution()
likelihood.AddDataset( likelihood.model_obj.GenerateDataset() )

CONSTRAINTS = True
PAR_LIMITS = True

if PAR_LIMITS:
	# Next, set limits so that none of the PDFs go negative in the fit.
	for var in likelihood.variable_list:
	    if 'Bb0n' in var['Name']:
	        likelihood.SetVariableLimits( var['Name'], \
	                                  lower_limit = -100., \
	                                  upper_limit = 100.)
	    else: 
	        likelihood.SetVariableLimits( var['Name'], \
	                                  lower_limit = 0., \
	                                  upper_limit = var['Value']*10.)


# Increase the step size for the Bb0n variable
likelihood.SetFractionalMinuitInputError('Num_FullLXeBb0n', 0.01/0.0001)

# Set the number of datasets to generate
num_datasets = 2000

output_row = dict()
output_df_list = []

import time

start_time = time.time()
last_time = start_time

for j in range(0,num_datasets):
	#print('Running dataset {}'.format(j))
	best_fit_converged = True
	fixedSig_fit_converged = True
	this_lambda = -1.
	output_row = dict()
	best_fit_parameters = None
	best_fit_errors = None
	fixed_fit_parameters = None
	fixed_fit_errors = None

	likelihood.model_obj.UpdateVariables(initial_values)
	likelihood.model_obj.GenerateModelDistribution()
	likelihood.AddDataset( likelihood.model_obj.GenerateDataset() )

	# Save input values as dict
	input_parameters = dict()
	for var in likelihood.variable_list:
		if 'Bb0n' in var['Name']:
			input_parameters[ var['Name'] ] = input_num_signal
		else:
			input_parameters[ var['Name'] ] = float(var['Value'])

	likelihood.SetAllVariablesFloating()
	likelihood.SetVariableFixStatus('Num_FullTPC_Co60',True)
	
	if CONSTRAINTS:
		rn222_idx = likelihood.GetVariableIndex('Rn222')
		# Fluctuate Rn222 constraint
		rn222_constraint_val = (np.random.randn()*0.1 + 1)*initial_values[rn222_idx]
		# Set Rn222 constraint
		likelihood.SetGaussianConstraintAbsolute(likelihood.variable_list[rn222_idx]['Name'],\
							 rn222_constraint_val, \
	                	                         0.1 * initial_values[rn222_idx])

	print('\n\nRunning dataset {}....\n'.format(j))
	likelihood.PrintVariableList()

	print('\nConstraints:')
	for constraint in likelihood.constraints:
		print('\t{}'.format(constraint))
	print('\n')

	print('\nBest fit:\n')
	best_fit_converged, best_fit_covar_flag, best_fit_iterations = \
			likelihood.CreateAndRunMinuitFit( initial_values, print_level=1 )
	likelihood.PrintVariableList()
	best_fit_parameters = dict( likelihood.fitter.values ) 
	best_fit_errors = dict( likelihood.fitter.errors )

	nll_best = likelihood.fitter.fval
	 
	print('\n\nFit with signal value fixed at {:3.3} cts:\n'.format(input_num_signal))

	likelihood.SetVariableFixStatus('Num_FullLXeBb0n',True)

	fixed_fit_converged, fixed_fit_covar_flag, fixed_fit_iterations = \
			likelihood.CreateAndRunMinuitFit( initial_values, print_level=1 )
	fixed_fit_parameters = dict( likelihood.fitter.values )
	fixed_fit_errors = dict( likelihood.fitter.errors ) 

	output_row['num_signal'] = input_num_signal
	output_row['lambda'] = this_lambda
	output_row['best_fit_converged'] = best_fit_converged
	output_row['best_fit_covar'] = best_fit_covar_flag
	output_row['fixed_fit_converged'] = fixed_fit_converged
	output_row['fixed_fit_covar'] = fixed_fit_covar_flag
	output_row['best_fit_parameters'] = best_fit_parameters
	output_row['best_fit_errors'] = best_fit_errors
	output_row['fixed_fit_parameters'] = fixed_fit_parameters
	output_row['fixed_fit_errors'] = fixed_fit_errors
	output_row['input_parameters'] = input_parameters
	#output_row['dataset'] = likelihood.dataset

	output_df_list.append(output_row)	
	
	print('Dataset {} finished at {:4.4}s'.format(j,time.time()-last_time))
	last_time = time.time()

output_df = pd.DataFrame(output_df_list)
#print(output_df.head())
print('Saving file to output directory: {}'.format(output_dir))
output_df.to_hdf('{}/critical_lambda_calculation_num_sig_{:3.3}_file_{}.h5'.format(output_dir,input_num_signal,iteration_num),key='df')

print('Elapsed: {:4.4}s'.format(time.time()-start_time))
