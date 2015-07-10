#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  protein.py
#  
#  Copyright 2015 farminf <farminf@farminf-3>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


import os

#INPUT
''' 
               PHI      PSI     OMEGA  (ignored)
-------------------------------------- (ignored)
    1  PHE  9999.000  151.378  171.421
    2  ASN   -58.229  155.541 -171.776
    3  MET   -52.501  -43.947  179.557
    4  GLN   -56.571  -60.035 -178.475
    5  CYS   -60.739  -44.202  176.328
    6  GLN   -52.221  -51.881  178.484
    7  ARG   -64.351  -39.167  177.086
    8  ARG   -59.962  -47.298  175.488
    9  PHE   -60.302  -52.855 -178.539
   10  TYR   -51.870  -49.642  176.793
   11  GLU   -59.034  -54.732 -176.974
   12  ALA   -56.160  -42.588  173.851
   13  LEU   -56.987  -53.191 -173.908
   14  HIS   -89.816   55.507 -178.743
   15  ASP  -149.574   75.129 -176.261
   16  PRO   -73.333   -0.355  170.355
   17  ASN   -62.924   14.314  171.277
   18  LEU   -95.567   64.714 -128.574
   19  ASN   -88.405 -174.397 -176.955
   20  GLU   -52.890  -60.711 -176.026
   21  GLU   -58.526  -51.140 -178.737
   22  GLN   -58.176  -47.784 -176.793
   23  ARG   -66.976  -48.690  176.939
   24  ASN   -61.115  -40.110  174.987
   25  ALA   -61.881  -40.932  176.640
   26  LYS   -60.630  -51.193  175.878
   27  ILE   -54.228  -47.882  178.863
   28  LYS   -57.494  -46.893 -179.236
   29  SER   -62.501  -46.749  176.762
   30  ILE   -59.051  -45.000 -178.824
   31  ARG   -66.474  -36.719  175.099
   32  ASP   -62.725  -48.565  178.415
   33  ASP   -60.397  -52.671 -174.902
   34  CYS   -68.000 9999.000 9999.000
'''



# RESIDUES DICTIONARY
'''
 If Only Residue Names are Entered, the Default is to Build an Extended
 Conformation Using L-Amino Acids and Zwitterionic Termini

 Regular Amino Acids:  GLY, ALA, VAL, LEU, ILE, SER, THR, CYS, CYX, PRO,
 PHE, TYR, TRP, HIS, ASP, ASN, GLU, GLN, MET, LYS, ARG, ORN, AIB

 Alternative Protonation States:  CYD, TYD, HID, HIE, ASH, GLH, LYD

 N-Terminal Cap Residues:  H2N=Deprotonated, FOR=Formyl, ACE=Acetyl,
                           PCA=Pyroglutamic Acid
 C-Terminal Cap Residues:  COH=Protonated, NH2=Amide, NME=N-MethylAmide
'''


# OUTPUT  -  'protein.in'
'''
protein                                                  
Protein built with TINKER/Force Field Explorer
ALA  -49.0  -26.0
ALA  -57.0  -47.0
ALA  -57.0  -47.0
ALA  -57.0  -47.0
ALA  -57.0  -47.0
ALA  -57.0  -47.0
ALA  -57.0  -47.0
ALA  -57.0  -47.0


N
'''




class Protein:
    """ Class doc """
    
    def __init__ (self, name = 'teste', 
                       title = 'teste apenas', 
                  forcefield = 'oplsaa.prm',
               phi_psi_table = None, 
                    sequence = None,
                         log = None):

        
        self.name          = name                 #     Name to be Used for Output Files
        self.title         = title                #     Title - first line in the XYZ file
        self.forcefield    = forcefield
        self.phi_psi_table = phi_psi_table        
        self.sequence      = sequence     
        self.log           = log          


        self.one_letter={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
                         'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
                         'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
                         'GLY':'G', 'PRO':'P', 'CYS':'C'}

        self.three_letter={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN',               \
                           'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
                           'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
                           'G':'GLY', 'P':'PRO', 'C':'CYS'}
        self.forcefields  = {
                            'amber94.prm'      : '/home/farminf/Programas/tinker/params/amber94.prm'      ,
                            'amber96.prm'      : '/home/farminf/Programas/tinker/params/amber96.prm'      ,
                            'amber98.prm'      : '/home/farminf/Programas/tinker/params/amber98.prm'      ,
                            'amber99.prm'      : '/home/farminf/Programas/tinker/params/amber99.prm'      ,
                            'amber99sb.prm'    : '/home/farminf/Programas/tinker/params/amber99sb.prm'    ,
                            'amoeba04.prm'     : '/home/farminf/Programas/tinker/params/amoeba04.prm'     ,
                            'amoeba09.prm'     : '/home/farminf/Programas/tinker/params/amoeba09.prm'     ,
                            'amoebabio09.prm'  : '/home/farminf/Programas/tinker/params/amoebabio09.prm'  ,
                            'amoebapro04.prm'  : '/home/farminf/Programas/tinker/params/amoebapro04.prm'  ,
                            'amoebapro13.prm'  : '/home/farminf/Programas/tinker/params/amoebapro13.prm'  ,
                            'charmm19.prm'     : '/home/farminf/Programas/tinker/params/charmm19.prm'     ,
                            'charmm22cmap.prm' : '/home/farminf/Programas/tinker/params/charmm22cmap.prm' ,
                            'charmm22.prm'     : '/home/farminf/Programas/tinker/params/charmm22.prm'     ,
                            'dang.prm'         : '/home/farminf/Programas/tinker/params/dang.prm'         ,
                            'hoch.prm'         : '/home/farminf/Programas/tinker/params/hoch.prm'         ,
                            'iwater.prm'       : '/home/farminf/Programas/tinker/params/iwater.prm'       ,
                            'mm2.prm'          : '/home/farminf/Programas/tinker/params/mm2.prm'          ,
                            'mm3.prm'          : '/home/farminf/Programas/tinker/params/mm3.prm'          ,
                            'mm3pro.prm'       : '/home/farminf/Programas/tinker/params/mm3pro.prm'       ,
                            'mmff.prm'         : '/home/farminf/Programas/tinker/params/mmff.prm'         ,
                            'oplsaal.prm'      : '/home/farminf/Programas/tinker/params/oplsaal.prm'      ,
                            'oplsaa.prm'       : '/home/farminf/Programas/tinker/params/oplsaa.prm'       ,
                            'oplsua.prm'       : '/home/farminf/Programas/tinker/params/oplsua.prm'       ,
                            'smoothaa.prm'     : '/home/farminf/Programas/tinker/params/smoothaa.prm'     ,
                            'smoothua.prm'     : '/home/farminf/Programas/tinker/params/smoothua.prm'     ,
                            'tiny.prm'         : '/home/farminf/Programas/tinker/params/tiny.prm'         ,
                            'water03.prm'      : '/home/farminf/Programas/tinker/params/water03.prm'      ,
                            'water14.prm'      : '/home/farminf/Programas/tinker/params/water14.prm'     
                            }
        self.system  = {
                    'name'      : name      ,         
                    'title'     : title     ,        
                    'forcefield': forcefield,
                    'sequence'  : {
                                       #RES         Phi       Psi     Omega  Chi-Angles
                                   1  :['PHE' , 9999.000,  151.378,  171.421, '' ],
                                   2  :['ASN' ,  -58.229,  155.541, -171.776, '' ],
                                   3  :['MET' ,  -52.501,  -43.947,  179.557, '' ],
                                   4  :['GLN' ,  -56.571,  -60.035, -178.475, '' ],
                                   5  :['CYS' ,  -60.739,  -44.202,  176.328, '' ],
                                   6  :['GLN' ,  -52.221,  -51.881,  178.484, '' ],
                                   7  :['ARG' ,  -64.351,  -39.167,  177.086, '' ],
                                   8  :['ARG' ,  -59.962,  -47.298,  175.488, '' ],
                                   9  :['PHE' ,  -60.302,  -52.855, -178.539, '' ],
                                   10 :['TYR' ,  -51.870,  -49.642,  176.793, '' ],
                                   11 :['GLU' ,  -59.034,  -54.732, -176.974, '' ],
                                   12 :['ALA' ,  -56.160,  -42.588,  173.851, '' ],
                                   13 :['LEU' ,  -56.987,  -53.191, -173.908, '' ],
                                   14 :['HIS' ,  -89.816,   55.507, -178.743, '' ],
                                   15 :['ASP' , -149.574,   75.129, -176.261, '' ],
                                   16 :['PRO' ,  -73.333,   -0.355,  170.355, '' ],
                                   17 :['ASN' ,  -62.924,   14.314,  171.277, '' ],
                                   18 :['LEU' ,  -95.567,   64.714, -128.574, '' ],
                                   19 :['ASN' ,  -88.405, -174.397, -176.955, '' ],
                                   20 :['GLU' ,  -52.890,  -60.711, -176.026, '' ],
                                   21 :['GLU' ,  -58.526,  -51.140, -178.737, '' ],
                                   22 :['GLN' ,  -58.176,  -47.784, -176.793, '' ],
                                   23 :['ARG' ,  -66.976,  -48.690,  176.939, '' ],
                                   24 :['ASN' ,  -61.115,  -40.110,  174.987, '' ],
                                   25 :['ALA' ,  -61.881,  -40.932,  176.640, '' ],
                                   26 :['LYS' ,  -60.630,  -51.193,  175.878, '' ],
                                   27 :['ILE' ,  -54.228,  -47.882,  178.863, '' ],
                                   28 :['LYS' ,  -57.494,  -46.893, -179.236, '' ],
                                   29 :['SER' ,  -62.501,  -46.749,  176.762, '' ],
                                   30 :['ILE' ,  -59.051,  -45.000, -178.824, '' ],
                                   31 :['ARG' ,  -66.474,  -36.719,  175.099, '' ],
                                   32 :['ASP' ,  -62.725,  -48.565,  178.415, '' ],
                                   33 :['ASP' ,  -60.397,  -52.671, -174.902, '' ],
                                   34 :['CYS' ,  -68.000, 9999.000, 9999.000, '' ]
                                   }     
                        }
        
        
        
        self.tinker_protein_in  = None
        self.tinker_protein_out = None
        self.tinker_dynamics_in = None
        
    def ImportPhiPsiFromFile (self, filein):
        """ Function doc """
        print 'ImportPhiPsiFromFile'
        filein = open(filein, 'r')
        for line in filein:
            print line
        
    def ImportSequenceFromFile (self, filein):
        """ Function doc """
        print 'ImportSequenceFromFile'
        filein = open(filein, 'r')
        for line in filein:
            print line
    
    def ImportRestrantsFromFile (self, filein):
        """ Function doc """
        print 'ImportRestrantsFromFile'
        filein = open(filein, 'r')
        for line in filein:
            print line


    def FromTorsionsToSystem (self, filein, name = None, title =  None):
        """ 
        if name  = none, so name   = filename
        if title = none, so title  = filename + title
        
        Import system from a torsions file:
                       PHI      PSI     OMEGA
        --------------------------------------
            1  PHE  9999.000  151.378  171.421
            2  ASN   -58.229  155.541 -171.776
            3  MET   -52.501  -43.947  179.557
            4  GLN   -56.571  -60.035 -178.475
            5  CYS   -60.739  -44.202  176.328
            6  GLN   -52.221  -51.881  178.484
        """
        #NAME
        filename = os.path.basename(filein)
        filename = filename.split('.')
        filename = filename[0]
        self.system['name'] = filename
        #print filename
        
        
        # TORSIONS
        filein = open(filein, 'r')
        for line in filein:
            line2 = line.split()
            
            if len(line2) == 5:
                self.system['sequence'][int(line2[0])] = [line2[1],line2[2],line2[3],line2[4],''] 
        
            if len(line2) == 6:
                self.system['sequence'][int(line2[0])] = [line2[1],line2[2],line2[3],line2[4], line2[5]]         
        
        
    
    
    
    #-------------------------------------#     
    #                                     #
    #      Build Protein Input File       #
    #                                     #
    #-------------------------------------#     
    def BuildProteinInputFile (self, output = None):
        """ Function doc """
        fileout = self.system['name'] + '.in'
        fileout = open(fileout, 'w')
        
        text    = ''
        # 1 Name:
        text += self.system['name'] + '\n'
        
        # 2 title:
        text += self.system['name'] + ' ' + self.system['forcefield'] + '\n'
        
        # 3 Force field
        #text += self.forcefields[self.system['forcefield']] + '\n' 

        # 4 Sequence
        for i in self.system['sequence']:
            #print self.system[i][0], self.system[i][1], self.system[i][2]
            text += '%s %s %s %s %s\n'%(self.system['sequence'][i][0], 
                                    str(self.system['sequence'][i][1]), 
                                    str(self.system['sequence'][i][2]), 
                                    str(self.system['sequence'][i][3]), 
                                    str(self.system['sequence'][i][4]))
        
        # 5 Cyclize the Polypeptide Chain [N] :  
        text += '\n\nN'
        print text        
        
        fileout.write(text)
        fileout.close()
        
        self.tinker_protein_in  = self.system['name'] + '.in'
        self.tinker_protein_out = self.system['name'] + '.out'
    
    def BuildProteinKeyFile (self, output = None):
        """ Function doc """
        fileout = self.system['name'] + '.key'
        fileout = open(fileout, 'w')
        
        text    = ''
        # 1 Name:
        #text += '' + '\n'
        # Output Control
        text += 'ARCHIVE' + '\n'       
        text += 'OVERWRITE' + '\n'         
        text += 'VERBOSE' + '\n'           

        # Force Field Selection
        text += 'PARAMETERS        /home/farminf/MolecularTools/ffe/../tinker/params/oplsaa.prm' + '\n'

        # Nonbonded Cutoff Keywords
        text += 'CHG-CUTOFF        4' + '\n'       
        text += 'VDW-CUTOFF        4' + '\n'       

        
        
        fileout.write(text)
        fileout.close()
        
        
        
        
        
    def BuildXYZModel (self):               #to build the xyz model
        """ 
        if >
        sh: 1: Syntax error: Bad fd number
        
        
        ---------------------------------------------------------------------------------------------
        The problem could be, that in Ubuntu 11.x /bin/sh is linked to /bin/dash and not to bin bash.

        check the link:

        ls -l /bin/sh
        If /bin/sh is a link to /bin/dash, change it to /bin/bash.

        sudo mv /bin/sh /bin/sh.orig
        sudo ln -s /bin/bash /bin/sh
        ---------------------------------------------------------------------------------------------
        """
        #print 'RunProtein'
        
        self.BuildProteinInputFile()
        self.BuildProteinKeyFile()
        
        print self.tinker_protein_in, self.tinker_protein_out
        print 'protein <'+self.tinker_protein_in +'>& '+ self.tinker_protein_out

        os.system('protein <'+self.tinker_protein_in +'>& '+ self.tinker_protein_out)
        os.system('')
def main():
    molecule = Protein()
    molecule.FromTorsionsToSystem ('peptide.ss')
    #molecule.BuildProteinInputFile()
    molecule.BuildXYZModel()
    return 0

if __name__ == '__main__':
    main()

