########################################################################
##  
## This is a class to map names of components in the spreadsheet
## to names of components in the rootfile names. This should be a 
## temporary thing...
##      - B. Lenardo (7.16.2018)
##
#######################################################################



class NameDict:

  def __init__(self):
    self.data = {}
    self.data["OuterCryoResin"] = 'Outer Cryostat (Resin)'
    self.data["OuterCryoFiber"] = 'Outer Cryostat (Fiber)'
    self.data["OuterCryoSupportResin"] = 'Outer Cryostat Support (Resin)'
    self.data["OuterCryoSupportFiber"] = 'Outer Cryostat Support (Fiber)'
    self.data["InnerCryoResin"] = 'Inner Cryostat (Resin)'
    self.data["InnerCryoFiber"] = 'Inner Cryostat (Fiber)'
    self.data["InnerCryoSupportResin"] = 'Inner Cryostat Support (Resin)'
    self.data["InnerCryoSupportFiber"] = 'Inner Cryostat Support (Fiber)'
    self.data["InnerCryoLiner"] = 'Inner Cryostat Liner'
    self.data["HFE"] = 'HFE'
    self.data["HVTubes"] = 'HV Tubes'
    self.data["HVCables"] = 'HV Cables'
    self.data["HVFeedthru"] = 'HV Feedthrough'
    self.data["HVFeedthruCore"] = 'HV Feedthrough Core'
    self.data["HVPlunger"] = 'HV Plunger'
    self.data["CalibrationGuideTube1"] = 'Calibration Guide Tube 1'
    self.data["CalibrationGuideTube2"] = 'Calibration Guide Tube 2'
    self.data["TPC"] = 'TPC Vessel'
    self.data["TPCSupportCone"] = 'TPC Support Cone'
    self.data["CathodeRn"] = 'Cathode (Radon)'
    self.data["Cathode"] = 'Cathode'
    self.data["Bulge"] = 'Bulge'
    self.data["FieldRings"] = 'Field Rings'
    self.data["SupportSpacers"] = 'Support Rods and Spacers'
    self.data["SiPMStaves"] = 'SiPM Staves'
    self.data["SiPMModule"] = 'SiPM Module (Interposer)'
    self.data["SiPMElectronics"] = 'SiPM Electronics'
    self.data["SiPMCables"] = 'SiPM Cables'
    self.data["SiPMs"] = 'SiPMs'
    self.data["ChargeTilesCables"] = 'Charge Tiles Cables'
    self.data["ChargeTilesElectronics"] = 'Charge Tiles Electronics'
    self.data["ChargeTilesSupport"] = 'Charge Tiles Support'
    self.data["ChargeTilesBacking"] = 'Charge Tiles Backing'
    self.data["LXe"] = 'Full LXe'
    #        self.data["bb2n"] = 'Full LXe'
    #        self.data["bb0n"] = 'Full LXe'
    self.data["ActiveLXe"] = 'Active LXe'
    self.data["InactiveLXe"] = 'Inactive LXe'
    self.data["SolderAnode"] = 'Solder (Anode)'
    self.data["SolderSiPM"] = 'Solder (SiPM)'
