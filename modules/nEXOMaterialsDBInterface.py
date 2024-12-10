# -*- coding: utf-8 -*-
# ^this line is required for correctly interpreting the cm^2 in some specific activities

import numpy as np

import re, json, jsonschema
from cloudant.client import CouchDB

from nEXOMaterialsDBSchemas import schemas, validate

class nEXOMaterialsDBInterface:


  ########################################################################################
  def __init__(self, server_choice='default'):

    if server_choice == 'local':
      from getpass import getpass
      password = getpass()
      client = CouchDB('admin', password, url='http://127.0.0.1:5984/', connect=True) #, use_basic_auth=True)
    
    elif server_choice == 'nexo':  
      username = input("Username: ")
      from getpass import getpass
      password = getpass()
      client = CouchDB(username, password, url='https://nexo.ph.ua.edu/', connect=True, use_basic_auth=True)

    elif server_choice == 'default':
      client = CouchDB('xenon','136', url='https://nexo.ph.ua.edu/', connect=True, use_basic_auth=True)
    
    self.db = client['material_database']
    self.keymap = {}   # key-to-id map

    # Masses and halflives are mainly from Google/Wikipedia
    self.halflives_yrs = { 'U-238': 4.468e9,
                       'Th-232': 1.405e10,
                       'K-40': 1.251e9,
                       'Co-60': 5.26,
                       'Cs-137': 30.17,
                       'Rn-222': (3.822 / 365.),
                       'Bi-214': (19.9 / 60. / 24. / 365.),
                       'Xe-137': (3.8 / 60. / 24. / 365. ),
                       'Al-26': 7.17e5,
                       'Ag-110m': (249.95 / 365.),
                       'Co-56': (77.233 / 365.),
                       'bb2n': 2.165e21,
                       'bb0n': 1.e28,
                       'B8nu': 1.e100 }
                     
    self.masses_amu = { 'U-238': 238.0507882,
                       'Th-232': 232.0380553,
                       'K-40': 39.96399817,
                       'Co-60': 59.9338171,
                       'Cs-137': 136.9070895,
                       'Rn-222': 222.0175777,
                       'Bi-214': 213.998712,
                       'Xe-137': 136.911562,
                       'Al-26': 25.98689186,
                       'Ag-110m': 109.906107,
                       'Co-56': 55.9398393,
                       'bb2n': 135.907219,
                       'bb0n': 135.907219,
                       'B8nu': 0. }


  ########################################################################################
  # Get raw CouchDB object 
  def GetDB( self ):
    return self.db

  # Retrieve/Update/Delete document



  ########################################################################################
  # Presort documents by human-readable key
  def PrepareDB( self ):

    for i,doc in enumerate(self.db):
      #if i % 10 == 0: 
         #print(i)
         #print('{}, {}: {}'.format(doc['key'],doc['_id'],doc.keys()))
         #print('key: {}'.format(doc['key']))
         #print('doctype: {}'.format(doc['doctype']))
         #try:
         #    if 'MC' in doc['key']: print(doc['rootfiles'])
         #except KeyError:
         #    print('No rootfiles in keys.')
      if 'deleted' not in doc.keys() and 'hidden' not in doc.keys() and 'key' in doc.keys():
        k = doc['key']
        if k in self.keymap.keys():
          self.keymap[k].append(doc['_id'])
          print('WARNING! Duplicate keys.', k, self.keymap[k])
        else:
          self.keymap[k] = [doc['_id']]

  ########################################################################################
  # Print out all the available tagged geometries and their titles, for bookkeeping. 
  def PrintAllGeometryTags( self ):
      for key, doc_id in self.keymap.items():
          if 'D-' in key:
             doc = self.GetDoc( key )
             print('{}: {}'.format(key,doc['title']))

  
  ########################################################################################
  def GetTaggedGeometry( self, tag ):
      doc = self.GetDoc( tag )
      return doc

  ########################################################################################
  # Get multiple documents by human-readable keys
  def GetDocs( self, key ):
    mainkey = key.split('.')[0]
    if len(self.keymap) == 0:
      self.PrepareDB()
    docids = self.keymap[mainkey]
    return [self.db[docid] for docid in docids]


  ########################################################################################
  # Get a single document by human-readable key
  def GetDoc( self, key ):
    mainkey = key.split('.')[0]
    return self.GetDocs(mainkey)[0]


  ########################################################################################
  def SearchRadioassayDocs( self, query, word=True):
    # query is dict of { field: search_term }
    hits = {}
    for i, doc in enumerate(self.db):
      if 'deleted' in doc.keys() or 'hidden' in doc.keys() or 'key' not in doc.keys(): continue
      if 'doctype' not in doc.keys() or (doc['doctype'] != 'radioassay' and doc['doctype'] != 'retro'): continue
      for search_field in query.keys():
        if search_field not in doc.keys(): continue
        if (word and query[search_field].lower() in re.split('\W+', doc[search_field].lower()))\
            or (not word and query[search_field].lower() in doc[search_field].lower()):
          #hits.append({'key': doc['key'], 'doc': doc})
          key = doc['key']
          hits[key] = {'doc': doc}
    return hits


  ########################################################################################
  def SearchRadioassayData( self, query, word=True ):
    hits = {}
    radioassay_doc_hits = self.SearchRadioassayDocs( query, word )
    for key, doc in radioassay_doc_hits.items():
      doc = doc['doc']  # Dunno why this is a dict...maybe will be useful later on
      if key[0] == 'R':
        for isample, sample in enumerate(doc['samples']):
          for icounting, counting in enumerate(sample['counting']):
            for ianalysis, analysis in enumerate(counting['analysis_result']):
              result_id = '%s.%i.%i.%i' % (doc['key'], isample+1, icounting+1, ianalysis+1)
              hits[result_id] = { 'doc': doc, 'sample': sample, 'counting': counting, \
                                  'measurement': analysis, 'results': analysis['results_new'] }
      elif key[0] == 'P':
        key = '{}'.format(doc['key'])
        hits[key] = {'doc': doc, 'results': doc['measurements']}
    return hits

 
  ########################################################################################
  def GetRadioassayData( self, radioassay_id ):
    # key = 'R-xxx.A.B.C' or 'P-xxx'
    if radioassay_id[0] == 'R': 
      indices = [s.strip() for s in radioassay_id.split('.')]
      sample_id = int(indices[1])-1
      counting_id = int(indices[2])-1
      analysis_id = int(indices[3])-1
      doc = self.GetDoc(indices[0])
      ans = doc['samples'][sample_id]['counting'][counting_id]['analysis_result'][analysis_id]['results_new']
      return ans
    elif radioassay_id[0] == 'P':
      doc = self.GetDoc( radioassay_id )
      return doc['measurements']
    else:
      print('ERROR: The document ID {} is invalid or it does not contain radioassay data.'.format(radioassay_id))
      return None

  # Create document 


  ########################################################################################
  def ConvertSpecActivTo_mBqPerUnitSize( self, specific_activity, unit, isotope ):
      
      kg_per_amu = 1./( 6.02e23 * 1000. )
      seconds_per_year = 60. * 60. * 24. * 365.

      if unit == 'mBq/kg' or unit == 'mBq/cm^2' or unit == u'mBq/cm\u00B2':
         return self.ConvertToFloat( specific_activity )

      mBq_per_kg = 1. / ( self.masses_amu[isotope] * kg_per_amu ) * \
                   1. / ( self.halflives_yrs[isotope]/np.log(2.) * seconds_per_year ) * \
                   1000. # This last factor of 1000 ensures we're in mBq rather than Bq

      # Note: in the materials database, "ppt" and "ppb" are actually per-mass; for
      # example, 1 ppt means something like "1 g U238 per 10^12 g Cu"
      if unit == 'ppt' or unit == 'pg/g':
         new_specific_activity = self.ConvertToFloat( specific_activity ) * 1.e-12 * mBq_per_kg
      elif unit == 'ppb' or unit == 'ng/g':
         new_specific_activity = self.ConvertToFloat( specific_activity ) * 1.e-9  * mBq_per_kg
      elif unit == 'pg/cm^2' or unit == u'pg/cm\u00B2':
         new_specific_activity = self.ConvertToFloat( specific_activity ) * 1.e-12 * (mBq_per_kg/1000.)
      else: 
         print('\nUNITS ERROR: {} not a supported unit - could not do specific activity conversion.'.format(unit))
         raise TypeError     
 
      return new_specific_activity


  ########################################################################################
  def ConvertToFloat(self,s):
      if s.strip() == '': return 0.
      return float(s)






#  ########################################################################################
#  def ConvertToSpecificActivity(self, d, version='R'):  # in mBq/kg
#      
#      kg_per_amu = 1./( 6.02e23 * 1000. )
#      seconds_per_year = 60. * 60. * 24. * 365.
#
#      mBq_per_kg = dict()
#      for isotope in mass_amu.keys():
#          mBq_per_kg[isotope] = 1. / ( self.masses_amu[isotope] * kg_per_amu ) *\
#                                1. / ( self.halflives_yrs[isotope]/np.log(2.) * seconds_per_year ) *\
#                                1000.  # this last factor of 1e3 converts Bq/kg into mBq/kg
#                               
#      concetration_units = { 'ppb': 1.e-9,
#                             'ppt': 1.e-12 }
# 
#      if version[0] == 'R':
#    
#        if d['unit'] in concentration_units.keys() and d['isotope'] in mBq_per_kg.keys():
#
#          scaling_factor = mBq_per_kg[ d['isotope'] ] * concentration_units[ d['unit'] ]
#          d['specific_activity'] = convert_to_float(d['specific_activity']) * scaling_factor
#          d['error'] = convert_to_float(d['error']) * scaling_factor
#          d['asym_lower_error'] = convert_to_float(d['asym_lower_error']) * scaling_factor
#          d['unit'] = 'mBq/kg'
#    
#        elif d['unit'] == 'mBq/kg': 
#
#          d['specific_activity'] = convert_to_float(d['specific_activity'])
#          d['error'] = convert_to_float(d['error'])
#          d['asym_lower_error'] = convert_to_float(d['asym_lower_error'])
#          d['unit'] = 'mBq/kg'
#
#      elif version[0] == 'P':
#
#        if d['unit'] in factor2.keys() and d['isotope'] in factor1.keys():
#
#          factor = factor1[d['isotope']]*factor2[d['unit']]
#          d['value'] = convert_to_float(d['value'])*factor
#          d['error'] = convert_to_float(d['error'])*factor
#          d['unit'] = 'mBq/kg'
#
#        elif d['unit'] == 'mBq/kg': 
#
#          d['value'] = convert_to_float(d['value'])
#          d['error'] = convert_to_float(d['error'])
#          d['unit'] = 'mBq/kg'
#
#      return d
#
#
#
#  ########################################################################################
#def format_data(d,version='R',style='raw',isotope=None):
#
#  if style == 'raw':
#    if version[0] == 'R':
#      if d['error_type'] == 'Upper limit (90% C.L.)':
#        return "<%s %s" % (d['specific_activity'], d['unit'])
#      if d['error_type'] == 'Symmetric error (68% C.L.)':
#        return "%s B1 %s %s" % (d['specific_activity'], d['error'], d['unit'])
#      if d['error_type'] == 'Range':
#        return "[%.3g, %.3g] %s" % (float(d['specific_activity'])-float(d['error']), float(d['specific_activity'])+float(d['error']), d['unit'])
#      if d['error_type'] == 'Asymmetric error':
#        return "%s +%s-%s %s" % (d['specific_activity'], d['error'], d['asym_lower_error'], d['unit'])
#    elif version[0] == 'P':
#      if d['type'] == 'limit':
#        return "<%s %s" % (d['value'], d['unit'])
#      if d['type'] == 'obs':
#        return "%s B1 %s %s" % (d['value'], d['error'], d['unit'])
#
#  elif style == 'specact': # in mBq/kg
#    if isotope == None: return None
#    d = convert_to_specact(d,version)
#    if version[0] == 'R':
#      if d['error_type'] == 'Upper limit (90% C.L.)':
#        return "<%s %s" % (d['specific_activity'], d['unit'])
#      if d['error_type'] == 'Symmetric error (68% C.L.)':
#        return "%s B1 %s %s" % (d['specific_activity'], d['error'], d['unit'])
#      if d['error_type'] == 'Range':
#        return "[%.3g, %.3g] %s" % (float(d['specific_activity'])-float(d['error']), float(d['specific_activity'])+float(d['error']), d['unit'])
#      if d['error_type'] == 'Asymmetric error':
#        return "%s +%s-%s %s" % (d['specific_activity'], d['error'], d['asym_lower_error'], d['unit'])
#    elif version[0] == 'P':
#      if d['type'] == 'limit':
#        return "<%s %s" % (d['value'], d['unit'])
#      if d['type'] == 'obs':
#        return "%s B1 %s %s" % (d['value'], d['error'], d['unit'])
#
#  elif style == 'ul-specact': # Upper limits in mBq/kg
#    if isotope == None: return None
#    d = convert_to_specact(d,version)
#    if version[0] == 'R':
#      if d['error_type'] == 'Upper limit (90% C.L.)':
#        if d['unit'] == 'mBq/kg':
#          return "%s" % (d['specific_activity'])
#        else: 
#          return "<%s %s" % (d['specific_activity'], d['unit'])
#      if d['error_type'] == 'Symmetric error (68% C.L.)':
#        if d['unit'] == 'mBq/kg':
#          return "%s" % (max(0,d['specific_activity']) + 1.64*d['error'])
#        else:
#          return "%s B1 %s %s" % (d['specific_activity'], d['error'], d['unit'])
#      if d['error_type'] == 'Range':
#        return "[%.3g, %.3g] %s" % (float(d['specific_activity'])-float(d['error']), float(d['specific_activity'])+float(d['error']), d['unit'])
#      if d['error_type'] == 'Asymmetric error':
#        return "%s +%s-%s %s" % (d['specific_activity'], d['error'], d['asym_lower_error'], d['unit'])
#    elif version[0] == 'P':
#      if d['type'] == 'limit':
#        if d['unit'] == 'mBq/kg':
#          return "%s" % (d['value'])
#        else: 
#          return "<%s %s" % (d['value'], d['unit'])
#        return "<%s %s" % (d['value'], d['unit'])
#      if d['type'] == 'obs':
#        if d['unit'] == 'mBq/kg':
#          return "%s" % (max(0,d['value']) + 1.64*d['error'])
#        else:
#          return "%s B1 %s %s" % (d['value'], d['error'], d['unit'])
#  return None
#
#
#
#  ########################################################################################
#if __name__ == '__main__':
#
#  #server = 'default'
#  server = 'local'
#  db = nEXOMaterialsDatabase(server)
#
#  #print db.get_next_id('R')
#  
#  #mcdocs = db.get_doc('MC-012')
#  #mcdocs['history'] = ''
#  #print(json.dumps(mcdocs, sort_keys=True, indent=4))
#
#  #docs = db.search_radioassay_data({'title':'teflon'})
#  hits = db.SearchRadioassayData({'title':'copper'})
#  for hit in hits:
#    print(hit['key'])
