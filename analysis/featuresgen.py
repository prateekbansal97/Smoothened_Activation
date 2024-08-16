import numpy as np
import pickle

distextra1 = ['CA_1', 'CA_9', 'CA_51', 'CA_58', 'CA_60', 'CA_106', 'CA_150', 'CA_233', 'CA_278', 'CA_282', 'CA_324', 'CA_413', 'CA_430', 'CA_445', 'CA_459', 'CA_28', 'CA_49', 'CA_62', 'CA_94', 'CA_149', 'CA_150', 'CA_165', 'CA_174', 'CA_195', 'CA_195', 'CA_199', 'CA_205', 'CA_282', 'CA_282', 'CA_282', 'CA_286', 'CA_323', 'CA_340', 'CA_340', 'CA_362', 'CA_376', 'CA_376', 'CA_386', 'CA_477']
distextra2 = ['CA_2', 'CA_116', 'CA_52', 'CA_154', 'CA_149', 'CA_111', 'CA_152', 'CA_234', 'CA_282', 'CA_478', 'CA_335', 'CA_417', 'CA_434', 'CA_447', 'CA_462', 'CA_76', 'CA_73', 'CA_66', 'CA_106', 'CA_154', 'CA_151', 'CA_458', 'CA_233', 'CA_202', 'CA_211', 'CA_211', 'CA_296', 'CA_365', 'CA_392', 'CA_396', 'CA_392', 'CA_342', 'CA_417', 'CA_420', 'CA_400', 'CA_377', 'CA_378', 'CA_389', 'CA_480']

distfeature_array_1 = ['OH_149','NH1_232','O_277','CZ3_281','CD2_412','OH_429','C_27','CZ3_148','OH_149','ND1_173','CG_194','CD2_204','CG_285','NE2_322','CE2_339','OH_339','OG1_476','NH1_393','SD_474','CB_212','NE1_479','NE1_281','CD2_307','OH_336','CZ_342','NE1_51']
distfeature_array_2 = ['OD2_151','NH1_233','NE1_281','CZ3_477','CZ_416','NE2_433','NZ_75','OE2_153','OE2_150','NH1_232','CG_210','CG_295','SD_391','OH_341','CE1_416','CB_419','NE1_479','NE1_477','NE1_477','NE1_477','CB_189','CG_205','CD2_273','OH_149','CD_460','OH_72']

distfeature_array_1 += distextra1
distfeature_array_2 += distextra2

features_tunnel_extra1 = ['NH2_342','NH2_342','OD1_415','CG_205','CD2_394','NH1_393','CD1_209','CD1_277','CD1_277','CD1_397','CE1_397','CG1_263','CD2_267','CZ_274','CG_281','CG2_208','N_398','NH2_393','CE2_404','ND1_412','CD2_464']
features_tunnel_extra2 = ['OD2_415','OE1_460','OE1_460','CD1_394','CZ3_477','CZ3_477','CE3_477','CD1_397','CD2_209','CD2_209','CE3_477','OH_264','CE_268','CZ2_273','CA_280','CD2_209','CZ_397','CD2_394','CA_405','CD1_413','ND2_463']

features = {'Distances': distfeature_array_1 + distfeature_array_2, 'Distancestunnel': features_tunnel_extra1 + features_tunnel_extra2}

pickle.dump(features, open('./pkl/features.pkl', 'wb'))

parm = {'5L7D': '/home/pdb3/SMO/APO/Analysis/5L7D_Apo_HMass.parm7', '6XBL': '/home/pdb3/SMO/APO/Analysis/6XBL_Apo_HMass.parm7'}
parmstrip = {'5L7D': '/home/pdb3/SMO/APO/Analysis/5L7D_Apo_HMass_stripped.parm7', '6XBL': '/home/pdb3/SMO/APO/Analysis/6XBL_Apo_HMass_stripped.parm7'}
parmpaths = {'strip': parmstrip, 'full': parm}

pickle.dump(parmpaths, open('./pkl/parmpaths.pkl', 'wb'))

