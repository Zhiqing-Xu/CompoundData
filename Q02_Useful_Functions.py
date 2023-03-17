#!/usr/bin/env python
# coding: utf-8
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Microsoft VS header
#--------------------------------------------------#
import os 
import sys
import os.path
from sys import platform
from pathlib import Path
#--------------------------------------------------#
if __name__ == "__main__":
    print("="*80)
    if os.name == 'nt' or platform == 'win32':
        print("Running on Windows")
        if 'ptvsd' in sys.modules:
            print("Running in Visual Studio")
#--------------------------------------------------#
    if os.name != 'nt' and platform != 'win32':
        print("Not Running on Windows")
#--------------------------------------------------#
    if "__file__" in globals().keys():
        try:
            os.chdir(os.path.dirname(__file__))
            print('CurrentDir: ', os.getcwd())
        except:
            print("Problems with navigating to the file dir.")
    else:
        print("Running in python jupyter notebook.")
        try:
            if not 'workbookDir' in globals():
                workbookDir = os.getcwd()
                print('workbookDir: ' + workbookDir)
                os.chdir(workbookDir)
        except:
            print("Problems with navigating to the workbook dir.")
#--------------------------------------------------#
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
#--------------------------------------------------#
import gzip
#--------------------------------------------------#
from AP_funcs import *


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# General Functions for reading SMILES
#====================================================================================================#
class Get_Unique_SMILES:
    # The purpose of this class is to help detect SMILES that represent same chemical compound.
    # All SMILES can be converted to a "unique" string with pre-defined settings.
    # SMILES that represent the same compound are expected to be converted to one unique string.
    # degree of uniqueness is determined by settings (i.e., isomericSmiles, SMARTS_bool, etc.)
    # Uniqueness Ranking:
    # 1. smiles -> mol -> smiles, with isomericSmiles = True  and SMARTS_bool = False
    # 2. smiles -> mol -> smiles, with isomericSmiles = False and SMARTS_bool = False
    # 3. smiles -> mol -> smiles, with isomericSmiles = False and SMARTS_bool = True
    # #1 can be used to parse SMILES string.
    # #2 is normally used to identify same compounds.
    # #3 can be used to deal with SMARTS-like strings output by reaction simulation.
    # #3 strings can be problematic when using RDKIT functions, so USE #3 VERY CAREFULLY.
    #--------------------------------------------------#
    def __init__(self, isomericSmiles = True, kekuleSmiles = False, canonical = True, SMARTS_bool = False):
        self.isomericSmiles = isomericSmiles
        self.kekuleSmiles   = kekuleSmiles 
        self.canonical      = canonical
        self.SMARTS_bool    = SMARTS_bool

    #--------------------------------------------------#
    def GetUniqueSMILES_fromSMILES(self, smiles_x): 

        #------------------------------
        if self.SMARTS_bool == True: # USe MolFromSmarts
            """
            # old_name: UniS(), UniSS()
            # SMILES -> R-mol -> SMARTS -> R-mol -> SMILES 
            """
            SMILES_bool = True
            mol_x=Chem.MolFromSmiles(smiles_x) # MolFromSmiles will NOT return ERROR, no matter what string is input.
            try:
                mol_x_unique=Chem.MolFromSmarts(Chem.MolToSmarts(mol_x))
                unique_smiles=Chem.MolToSmiles(mol_x_unique, 
                                               isomericSmiles = self.isomericSmiles,
                                               kekuleSmiles   = self.kekuleSmiles, 
                                               canonical      = self.canonical)
            except Exception:
                print ("!!!!! Problematic SMILES (Read using MolFromSmarts): ", smiles_x)
                unique_smiles = smiles_x
                SMILES_bool   = False
            return unique_smiles, SMILES_bool
        
        #------------------------------
        if self.SMARTS_bool == False: # USe MolFromSmiles
            """
            old_name: CanS()
            SMILES -> R-mol -> SMILES
            """
            SMILES_bool = True
            #print(smiles_x)
            mol_x = Chem.MolFromSmiles(smiles_x) # MolFromSmiles will NOT return ERROR, no matter what string is input.
            try:
                unique_smiles=Chem.MolToSmiles(mol_x, 
                                               isomericSmiles = self.isomericSmiles,
                                               kekuleSmiles   = self.kekuleSmiles, 
                                               canonical      = self.canonical)
            except Exception:
                print ("!!!!! Problematic SMILES (Read using MolFromSmiles): ", smiles_x)
                unique_smiles = smiles_x
                SMILES_bool   = False
            return unique_smiles, SMILES_bool
        
    #--------------------------------------------------#
    def GetUniqueSMILES_fromInChI(self, InChI_x, sanitize = True, removeHs = True, 
                                                 logLevel = None, treatWarningAsError = False): 
        """
        old_name: - 
        InChI -> R-mol -> SMILES
        """
        InChI_bool = True # Whether this InChI Key is valid.
        # MolFromInchi will return ERROR or NOT ?
        mol_x = Chem.inchi.MolFromInchi(InChI_x, sanitize = sanitize, removeHs = removeHs, 
                                           logLevel = logLevel, treatWarningAsError = treatWarningAsError) 
        
        try:
            unique_smiles=Chem.MolToSmiles(mol_x, 
                                            isomericSmiles = self.isomericSmiles,
                                            kekuleSmiles   = self.kekuleSmiles, 
                                            canonical      = self.canonical)
            
        except Exception:
            print ("!!!!! Problematic InChI (Read using MolFromInchi): ", InChI_x)
            unique_smiles = None
            InChI_bool   = False
        return unique_smiles, InChI_bool

    #--------------------------------------------------#
    def UniS(self, smiles_x): 
        return self.GetUniqueSMILES_fromSMILES(smiles_x)
    
    #--------------------------------------------------#
    def UNQSMI(self, smiles_x):
        # Return Unique SMILES
        return self.UniS(smiles_x)[0]

    #--------------------------------------------------#
    def ValidSMI(self, smiles_x):
        # Return False if input is NOT identified as SMILES
        return self.UniS(smiles_x)[1]
    
    #--------------------------------------------------#
    def UNQSMI_InChI(self, smiles_x):
        # Return Unique SMILES
        smiles_x = self.UniS(smiles_x)[0]
        mol_x = Chem.MolFromSmiles(smiles_x)
        return Chem.MolToInchiKey(mol_x)




#====================================================================================================#


def AP_convert_test1():
    # a, b and c are CoA, d is some random compound.
    smiles_a = "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCS)O"
    smiles_b = "O=C(NCCS)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3OP(=O)(O)O"
    smiles_c = "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O)[C@H](C(=O)NCCC(=O)NCCS)O"
    smiles_d = "[H]N=c1c(c([H])nc(n1[H])C([H])([H])[H])C([H])([H])[n+]1cc(CC(N)C(=O)O)c2ccccc21"

    inchi_a = "1S/C21H36N7O16P3S/c1-21(2,16(31)19(32)24-4-3-12(29)23-5-6-48)8-41-47(38,39)44-46(36,37)40-7-11-15(43-45(33,34)35)14(30)20(42-11)28-10-27-13-17(22)25-9-26-18(13)28/h9-11,14-16,20,30-31,48H,3-8H2,1-2H3,(H,23,29)(H,24,32)(H,36,37)(H,38,39)(H2,22,25,26)(H2,33,34,35)/t11-,14-,15-,16+,20-/m1/s1"
    
    
    inchikey_a = "RGJOEKWQDUBAIZ-IBOSZNHHSA-N"
    #====================================================================================================#
    print ("\nTesting GetUnqSmi.UNQSMI()")
    GetUnqSmi = Get_Unique_SMILES(isomericSmiles = False, kekuleSmiles = False, canonical = True, SMARTS_bool = False)
    print (GetUnqSmi.UNQSMI(smiles_a))
    print (GetUnqSmi.UNQSMI(smiles_b))
    print (GetUnqSmi.UNQSMI(smiles_c))

    print (GetUnqSmi.GetUniqueSMILES_fromInChI(inchi_a, removeHs = False))

    print (GetUnqSmi.UNQSMI_InChI(smiles_a))
    print (GetUnqSmi.UNQSMI_InChI(smiles_b))
    print (GetUnqSmi.UNQSMI_InChI(smiles_c))


    print (GetUnqSmi.UNQSMI("CC1C2C(N(CN2C3=C(N1)N=C(NC3=O)N)C4=CC=C(C=C4)CC(C(C(COC5C(C(C(O5)COP(=O)(O)OC(CCC(=O)O)C(=O)O)O)O)O)O)O)C"))
    print (GetUnqSmi.UNQSMI("CC1Nc2nc(N)[nH]c(=O)c2N2CN(c3ccc(CC(O)C(O)C(O)COC4OC(COP(=O)(O)OC(CCC(=O)O)C(=O)O)C(O)C4O)cc3)C(C)C12"))
    print (GetUnqSmi.UNQSMI("O=C[C@@H](O)[C@@H](O)[C@@H](O)CO"))




AP_convert_test1()

#====================================================================================================#


# #get information about the Molecular
# from bioservices import *
# info = UniChem()
# cid = info.get_src_compound_ids_from_inchikey('KDXKERNSBIXSRK-YFKPBYRVSA-N')[0]

import requests
def get_smiles_from_inchikey(inchikey):
    request = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES/JSON').json()
    return request['PropertyTable']['Properties'][0]['CanonicalSMILES']

print(get_smiles_from_inchikey('RGJOEKWQDUBAIZ-UHFFFAOYSA-N'))











