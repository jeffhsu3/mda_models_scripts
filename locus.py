"""
"""

class Locus(object):
    def __init__(self, idx, theta=4, basepairs=13000,
            mutation_rate=7.6E-8, recombination_rate=0):
        """ Class containing locus metainformation

        Arguments
        ---------
        basepairs : int 
        mutation : float
        recombination : float
            probability of recombination between bp of each locus per
            generation
        perc_cds : float
        """
        self.idx = idx
        self.theta = theta
        self.basepairs = basepairs
        self.mutation_rate = mutation_rate
        self.recombination_rate = recombination_rate



class Loci(object):
    """ Collection of locus
    """
    def __init__(self):
        pass
