"""Conversion between pdb and pcd format"""


import os
from time import strftime, localtime
import numpy as np

class _pdb_:

    def _pdb_to_pcd(self, ifile_, color_):
        """ Extracts coordinates from pdb files and convert into pcd format """

        if ifile_[-4:] != '.pdb':
            print("incorrect file extension")
            print("file format may be wrong --> no output")

        atom_list = []
        # get pdb atoms
        try:
            with open(ifile_, "r") as f:
                line = f.readline()            
                while line[0:4]!="ATOM":
                    line = f.readline()
                
                count = 0

                while line:
                    dat_in = line[0:80].split()
                    # delete empty lines
                    if len(dat_in) == 0:
                        line = f.readline()
                        continue
                        
                    if(dat_in[0] == "ATOM"):
                        atom_type = line[12:16].strip()
                        residue_type = line[17:20]
                        chain_id = line[21]
                        # residue sequence number 残基序列号
                        residue_id = int(line[22:26])
                        # 原子类型　C N O S
                        atom_name = line[76:78].strip()
                        # delete H
                        if atom_name != "H":
                            # get the coordinates
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])

                            atom_list.append([x, y, z, atom_type, residue_type, atom_name])
                            count = count + 1

                    line = f.readline()

        except:
            print("Cannot process pdb {}".format(ifile_))
            return -1, None, None


        ofilename = os.path.basename(ifile_).replace("pdb", "pcd")    
        try:
            ofile = open(ofilename, "w")
        except IOError:
            print("pdb_to_pcd: Cannot write to current directory: "
                  "{}. Please check for access rights.".format(os.getcwd()))
            ofile.close()
            return -1, None, None

        properties = []
        colors = []
        ofile.write("VERSION .7\nFIELDS x y z rgb\nSIZE 4 4 4 4\n"
                    "TYPE F F F F\nCOUNT 1 1 1 1\n")
        ofile.write("WIDTH {0}\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\n"
                    "POINTS {0}\nDATA ascii".format(len(atom_list)))
        
        #print(atom_list)
        index = 0

        for atom in atom_list:
            index = index+1
            atom_name = str(atom[5])
            x = float(atom[0])
            y = float(atom[1])
            z = float(atom[2])
            properties.append([index, atom_name])
            colors.append(color_[atom_name])
            ofile.write("\n{} {} {} {}".format(x, y, z, color_[atom_name]))

        ofile.close()                                                    
        
        #print(ofilename)
        return ofilename, properties, colors
    
class _pcd_:

    def _get_coordinates(self, ifile_):

        if ifile_[-4:] != '.pcd':
            print("incorrect file extension")
            print("file format may be wrong --> no output")
            if "DATA ascii" not in open(ifile_, 'r').read():
                return -1

        coordinates = []
        errors = []
        with open(ifile_, 'r') as f:
            data = f.read().split('\n')
        data = [l for l in data[10:] if l != '']

        for line_ in data:
            x, y, z, rgb = line_.split()
            try: 
                x = float(x)
            except ValueError:
                errors.append(-1)
            try:
                y = float(y)
            except ValueError:
                errors.append(-1)
            try:
                z = float(z)
            except ValueError:
                errors.append(-1)
            try:
                rgb = int(rgb)
            except ValueError:
                errors.append(-1)

            if -1 in errors:
                print("Errors while parsing PCD")
                break

            coordinates.append([x, y, z, rgb])

        if -1 not in errors:
            self.coordinates = coordinates
            return coordinates
        else:
            return -1
        
    def _write_pdb_fake(self, ofile_, coordinates_, atom_, atom_type_,
                                                    macromol_="PROTEIN"):

        
        if ofile_[-4:] != '.pdb':
            ofile_ += '.pdb'
        name = os.path.basename(ofile_)
        name = os.path.splitext(name)[0]

        of_string = ""
        of_string += "# Modified by ProCarePDB Version 1.0\n"
        of_string += "# Modification time: {}\n".format(
                                strftime("%a %d %b %Y %H:%M:%S", localtime()))
        of_string += "# Name: {}.pdb\n\n".format(name)

        of_string += "@<TRIPOS>MOLECULE\n"
        of_string += "{}\n".format(name)
        of_string += "{:>5}{:>6}{:>6}{:>6}{:>6}\n".format(
                                                len(coordinates_), 0, 0, 0, 0)
        of_string += "{}\n".format(macromol_)
        of_string += "NO_CHARGES\n"
        of_string += "@<TRIPOS>ATOM"

        # Fake residues in PDB: all RES
        residue = "RES"

        for i, point in enumerate(coordinates_):
            x, y, z, rgb = [*point]
            of_string += ("\n{:>7} {:<8} {:>9.4f} {:>9.4f} {:>9.4f} "
                          "{:<5} {:>5} {:<8} {:>9}".format(i+1,
                                                     atom_[rgb],
                                                     x,
                                                     y,
                                                     z,
                                                     atom_type_[rgb],
                                                     i+1,
                                                     residue+str(i+1),
                                                     0.0000
                                                    ))
        of_string += "\n@<TRIPOS>BOND"
        
        with open(ofile_, 'w') as of:
            of.write(of_string)
        print("written pdb to {}".format(ofile_))
        self.type = "pcd"
        return ofile_

        

class _volsite_cavity_pdb_(_pdb_, _pcd_):
    
    def __init__(self):

        # Version: 1.0 -- Only distinguish C(red), N(green), O(blue), S(yellow); delete H
        __COLOR = {"O":255,
                   "N":65280,
                   "C":16711680,
                   "S":16776960}


        __ATOM = {val:key for key, val in __COLOR.items()}

        # Version: 1.0 -- Only distinguish C, N, O, S; delete H
        __ATOM_TYPE = {"N":"N.am",
                        "O":"O.2",
                        "C":"C.3",
                        "S":"S"}

        # Version: 1.0 -- Ignore residue
        # __RESIDUE = {"OG":"SER",
        #                 "N":"ALA",
        #                 "O":"ALA",
        #                 "NZ":"LYS",
        #                 "CZ":"PHE",
        #                 "CA":"GLY",
        #                 "DU":"CUB",
        #                 "OD1":"ASP",}

        self.COLOR = __COLOR

        self.ATOM = __ATOM

        self.ATOM_TYPE = {key:__ATOM_TYPE[val] 
                                for key, val in __ATOM.items()}
        
        # Version: 1.0 -- Ignore residue
        # self.RESIDUE = {key:__RESIDUE[val] 
        #                         for key, val in __ATOM.items()}

        

    def pdb_to_pcd(self, ifile_):
        return self._pdb_to_pcd(ifile_, self.COLOR)
        


    # donot write pdb
    # def pcd_to_pdb(self, ifile_):
    #     return self._pcd_to_pdb(ifile_, self.ATOM,
    #                                      self.ATOM_TYPE,
    #                                      self.RESIDUE)



    # def write_pdb(self, ofile_, coordinates_):
    #     return self._write_pdb(ofile_, coordinates_, self.ATOM, self.ATOM_TYPE,
    #                 self.RESIDUE)

    def write_pdb_fake(self, ofile_, coordinates_):
        return self._write_pdb_fake(ofile_, coordinates_, self.ATOM, self.ATOM_TYPE)
