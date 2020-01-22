import histlite as hl
import pandas as pd
import numpy as np
import sys
import yaml
import time
#import iminuit

import nEXOExcelTableReader


class nEXOFitWorkspace:

   ##########################################################################
   # Constructor just initializes empty objects.
   ##########################################################################
   def __init__( self, config='./config/TUTORIAL_config.yaml'):
      
      # Initialize the class members to null objects.
      #     'df_components' will be a pandas dataframe which
      #              contains the individual simulation
      #              PDFs as histlite histograms and the
      #              relevant info from the materials DB
      #              (component mass, activity, etc.)
      #     'df_group_pdfs' will be a pandas dataframe
      #              which contains the distributions for
      #              each group, weighted appropriately
      #     'neg_log_likelihood' will be a customized
      #              object which computes a likelihood
      #              between the grouped PDFs and some toy
      #              dataset
      #     'minimizer' will be an iMinuit object
      #     'signal_counts' allows the user to manually set
      #              the number of 0nuBB events in the 
      #              dataset
      #     'livetime' is the expected livetime of the 
      #              experiment in units of seconds. 
      #              #TOFIX: this will need to be modified 
      #                      for isotopes that have livetimes
      #                      on the order of years or less    
      #     'histogram_axis_names' contains the names of the
      #              variables we're fitting. These should
      #              be already associated with the PDFs at
      #              an earlier stage in the data processing -
      #              the present code should be totally agnostic
      #              to the variables/binning/etc. 

      self.df_components = pd.DataFrame()       
      self.df_group_pdfs = pd.DataFrame()
      self.neg_log_likelihood = None
      self.minimizer = None
      self.signal_counts=0.0001
      self.livetime = 10. * 365. * 24. * 60. * 60.
      self.histogram_axis_names = None

      config_file = open(config,'r')
      self.config = yaml.load( config_file, Loader=yaml.SafeLoader )
      config_file.close()

   ##########################################################################
   # Loads the input dataframe from the h5 file generated by the
   # ConvertExcel2DataFrame.py script
   ##########################################################################
   def LoadComponentsTableFromFile( self, input_filename ):

      print('\nLoading input data froma previously-generated components table....')

      if not (self.df_components.empty):
         print('\nWARNING: There is already an input dataframe loaded. It ' +\
               'will be overwritten.')

      try:
         self.df_components = pd.read_hdf(input_filename,key='Components')
      except OSError as e:
          print('\nERROR: The input file must be an HDF5 file.\n')
          #sys.exit()
          return
      except KeyError as e:
          print('\n')
          print(e)
          print('\nERROR: The input file should contain the component-wise activity and PDF ' +\
                'information. This can be generated from the Materials Database excel ' +\
                'spreadsheets using the ConvertExcel2DataFrame.py script.\n')
          #sys.exit()
          return

      print('\nLoaded dataframe with {} components.'.format(len(self.df_components)))
      print('Contains the following quantities of interest:')
      for column_name in self.df_components.columns:
          print('\t{}'.format(column_name))

      if 'HistogramAxisNames' in self.df_components.columns:     
         self.histogram_axis_names = self.df_components['HistogramAxisNames'].iloc[0]
         print('\nFit variables:\n' +\
               '\t{}'.format(self.histogram_axis_names))
      else:
         print('WARNING: We don\'t have axis names for the histograms. ' +\
               'Please ensure they are set in the config file, or we might run ' +\
               'into problems.')

      return
   ################# End of LoadComponentsTableFromFile() ###################

   ##########################################################################
   # Creates the components table from the Excel spreadsheet.
   ##########################################################################
   def CreateComponentsTableFromXLS( self, excelFile, histogramsFile ):

      # Check and see if we're overwriting something here.
      if not (self.df_components.empty):
         print('\nWARNING: There is already an input dataframe loaded. It ' +\
               'will be overwritten.')
 
      start_time = time.time()      

      excelTableReader = nEXOExcelTableReader.nEXOExcelTableReader( excelFile, \
                                                              histogramsFile, \
                                                              config = self.config)
      try: 
         excelTableReader.ConvertExcel2DataFrame()
      except KeyError:
         sys.exit()
   
      self.df_components = excelTableReader.components

      # Print out some useful info.
      print('\nLoaded dataframe with {} components.'.format(len(self.df_components)))
      print('Contains the following quantities of interest:')
      for column_name in self.df_components.columns:
          print('\t{}'.format(column_name))
      if 'HistogramAxisNames' in self.df_components.columns:     
         self.histogram_axis_names = self.df_components['HistogramAxisNames'].iloc[0]
         print('\nFit variables:\n' +\
               '\t{}'.format(self.histogram_axis_names))
      else:
         print('WARNING: We don\'t have axis names for the histograms. ' +\
               'Please ensure they are set in the config file, or we might run ' +\
               'into problems.')

      # Store the components table in a file in case you want to
      # go back and look at it later.
      nameTag = excelTableReader.GetExcelNameTag( excelFile )
      outTableName = 'ComponentsTable_' + nameTag + '.h5'
      print('\n\nWriting table to file {}\n'.format(outTableName))
      excelTableReader.components.to_hdf( outTableName, key='Components' )

      end_time = time.time()
      print('Elapsed time = {:3.3} seconds ({:3.3} minutes).'.format( \
                              end_time-start_time, \
                             (end_time-start_time)/60. ) )
      return
 
   ##########################################################################
   # Creates grouped PDFs from the information contained in the input
   # dataframe.
   ##########################################################################
   def CreateGroupedPDFs( self ):

       print('\nCreating grouped PDFs....')
      
       if not (self.df_group_pdfs.empty):
          print('\nWARNING: Group PDFs have already been generated. ' +\
                'They are going to be overwritten.')
 
       self.df_group_pdfs = pd.DataFrame(columns = ['Group',\
                                                    'Histogram',\
                                                    'TotalExpectedCounts'])

       # Loop over rows in df_components, add histograms to the appropriate group.
       for index,row in self.df_components.iterrows():

         if row['Isotope']=='bb0n':
             totalExpectedCounts = self.signal_counts
         else:
             if row['SpecActiv'] > 0.:
                totalExpectedCounts = row['Total Mass or Area'] * \
                                      row['SpecActiv']/1000. * \
                                      row['TotalHitEff_K'] / row['TotalHitEff_N'] * \
                                      self.livetime
             else:
                totalExpectedCounts = 0
         
         
         if not (row['Group'] in self.df_group_pdfs['Group'].values):

             new_group_row = { 'Group' : row['Group'], \
                               'Histogram' : row['Histogram'].normalize((0,1,2),integrate=False) * \
                                             totalExpectedCounts, \
                               'TotalExpectedCounts' : totalExpectedCounts }

             self.df_group_pdfs = self.df_group_pdfs.append(new_group_row,ignore_index=True)
 
         else:

             group_mask = row['Group']==self.df_group_pdfs['Group']            
 
             self.df_group_pdfs['Histogram'].loc[ group_mask ] = \
                  self.df_group_pdfs['Histogram'].loc[ group_mask ] + \
                  row['Histogram'].normalize( (0,1,2), integrate=False) * \
                  totalExpectedCounts
             
             self.df_group_pdfs['TotalExpectedCounts'].loc[ group_mask ] = \
                  self.df_group_pdfs['TotalExpectedCounts'].loc[ group_mask ] + \
                  totalExpectedCounts

       # One more loop accomplishes two things:
       #     1. Normalize the grouped PDFs (now that each component has been added with
       #        the correct weight)
       #     2. Generate a summed PDF (including the expected weights) and append it.
       total_sum_row = {}
       for index,row in self.df_group_pdfs.iterrows():

         self.df_group_pdfs['Histogram'].loc[index] = row['Histogram'].normalize( (0,1,2), integrate=False )

         # TOFIX: Need better handling of negative totals, and need to figure out why 'Far' gives
         # me weird stuff. For now, I'm ignoring these.
         if (row['TotalExpectedCounts']>0.)&(row['Group']!='Off'): #&(row['Group']!='Far'):
             if len(total_sum_row)==0:
                 total_sum_row = {'Group' : 'Total Sum',\
                                  'Histogram' : row['Histogram'],\
                                  'TotalExpectedCounts' : row['TotalExpectedCounts']}
             else:
                 total_sum_row['Histogram'] = total_sum_row['Histogram'] + row['Histogram']
                 total_sum_row['TotalExpectedCounts'] = total_sum_row['TotalExpectedCounts'] + row['TotalExpectedCounts']
      
       self.df_group_pdfs = self.df_group_pdfs.append(total_sum_row,ignore_index=True)

       # Print out all the groups and their expected counts.
       print('\t{:<10} \t{:>15}'.format('Group:','Expected Counts:'))
       for index,row in self.df_group_pdfs.iterrows():
           print('\t{:<10} \t{:>15.2f}'.format( row['Group'], row['TotalExpectedCounts']))

       return 
   ######################## End of CreateGroupPDFs() ########################

   
   ##########################################################################
   # Creates negative log likelihood object
   ##########################################################################
   #def CreateNegLogLikelihood( self ):
       



