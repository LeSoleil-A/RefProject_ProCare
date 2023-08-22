import numpy as np
from sklearn.neighbors import NearestNeighbors
import math

class _similarity_metrics_:

    def tanimoto_similarity(self): # mettre en arg la classe score et utuliser les instances ...   
        similarity = float(self.n_identity)/\
            (self.fitsize + self.refsize - self.n_identity)
        return round(similarity, 4)
   
    def tversky_similarity(self, alpha_=0.95, beta_=0.05):
        similarity = float(self.n_identity)/(alpha_ * (self.fitsize - self.n_identity)\
                                        + beta_ * (self.refsize - self.n_identity)\
                                        + self.n_identity)
        return round(similarity, 4)
   
    def cosine_similarity(self):
        similarity = float(self.n_identity)/(math.sqrt(self.fitsize)*math.sqrt(self.refsize))
        return round(similarity, 4)
  
    def dice_similarity(self):
        similarity = 2 * float(self.n_identity)/(self.fitsize + self.refsize)
        return round(similarity, 4)

    def per_score(self):
        return round(float(self.n_identity)/self.fitsize, 4)
    
    def weighted_score(self):
        # computed from sc-PDB 2016 archive
        """
        TO DO: and the frequency of C N O S?
        """       
        frequency = {"CA": 0.3817,
                     "CZ": 0.1110,
                     "O": 0.0659,
                     "OG": 0.0867,
                     "OD1": 0.0593,
                     "N": 0.1376,
                     "NZ": 0.0579,
                     "DU": 0.0999}

        similarity = float(self.CA)/self.fitsize * 1/frequency["CA"] +\
                         float(self.CZ)/self.fitsize * 1/frequency["CZ"] +\
                         float(self.N)/self.fitsize * 1/frequency["N"] +\
                         float(self.OD1)/self.fitsize * 1/frequency["OD1"] +\
                         float(self.NZ)/self.fitsize * 1/frequency["NZ"] +\
                         float(self.OG)/self.fitsize * 1/frequency["OG"] +\
                         float(self.O)/self.fitsize * 1/frequency["O"] +\
                         float(self.DU)/self.fitsize * 1/frequency["DU"]
        return round(similarity, 4)


class _distances_:

    def hamming_distance(self):
        distance = self.fitsize + self.refsize - (2 * self.n_identity)
        return round(distance, 4)


    def soergel_distance(self):
        distance = float(self.fitsize + self.refsize- (2 * self.n_identity))/\
                                (self.fitsize + self.refsize - self.n_identity)
        return round(distance, 4)


class _ph4_ext_pdb_(_similarity_metrics_, _distances_):
    """
    all neighbors within distance D
    and strict correspondence of properties
    """

    def __init__(self, source_coordinates_, target_coordinates_,
                 source_properties_, target_properties_, distance_threshold_):
        self.source_coordinates = source_coordinates_
        self.target_coordinates = target_coordinates_
        self.source_properties = source_properties_
        self.target_properties = target_properties_
        self.distance_threshold = distance_threshold_
        self.n_identity = 0
        self.n_different = 0
        self.C=0
        self.N=0
        self.O=0
        self.S=0

        if len(self.source_coordinates) > len(self.target_coordinates):
            N = NearestNeighbors(n_neighbors=len(self.source_coordinates),
                                 radius=self.distance_threshold,
                                 algorithm='ball_tree').fit(self.source_coordinates)
            self.distances, self.indices = N.radius_neighbors(self.target_coordinates)
            self.fitProp = self.target_properties
            self.refProp = self.source_properties

        else:
            N = NearestNeighbors(n_neighbors=len(self.target_coordinates),
                                 radius=self.distance_threshold,
                                 algorithm='ball_tree').fit(self.target_coordinates)
            self.distances, self.indices = N.radius_neighbors(self.source_coordinates)
            self.fitProp = self.source_properties
            self.refProp = self.target_properties

        self.fitsize = len(self.fitProp)
        self.refsize = len(self.refProp)
    
    def get_similarity_by_rules(self):

        for i in range(len(self.fitProp)):
            # count the number of atoms with the same atom type in a small area
            # fitProp[i][1]: atom name, such as 'C','O'
            if self.fitProp[i][1] in [self.refProp[e][1] for e in self.indices[i]]:
                self.n_identity += 1
                if self.fitProp[i][1] == "C":
                    self.C += 1
                elif self.fitProp[i][1] == "O":
                    self.O += 1
                elif self.fitProp[i][1] == "N":
                    self.N += 1
                elif self.fitProp[i][1] == "S":
                    self.S += 1
            
        ratio_aligned = round(float(self.n_identity)/self.fitsize, 4)

        if self.n_identity != 0:
            ratio_C_in_aligned = round(float(self.C)/self.n_identity, 4)
            ratio_N_in_aligned = round(float(self.N)/self.n_identity, 4)
            ratio_O_in_aligned = round(float(self.O)/self.n_identity, 4)
            ratio_S_in_aligned = round(float(self.S)/self.n_identity, 4)
        else:
            ratio_C_in_aligned = 0
            ratio_N_in_aligned = 0
            ratio_O_in_aligned = 0
            ratio_S_in_aligned = 0
        
        return ratio_aligned, \
            ratio_C_in_aligned, ratio_N_in_aligned, \
            ratio_O_in_aligned, ratio_S_in_aligned