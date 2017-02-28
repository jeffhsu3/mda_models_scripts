"""
"""


class Villages(object):
    """ A collection of Villages
    """

    def __init__(self, villages):
        """
        Parameters
        ----------
        villages : a list of figs.Village
        """
        self.villages = villages
        self.current = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.current >= self.villages:
            raise StopIteration
        else:
            self.current += 1
            return self.current - 1


class Village(object):
    """ Basic parameters for the village
    """

    def __init__(self, identifier, hostpopsize = None,
            prevalence = None, distance = None, hours2bite = None, bitesPperson = None,
            bednets = None, bnstart = None, bnstop = None, bncoverage = None, muTrans = None,
            sizeTrans = None):
        """
        """
        self.id =  identifier
        self.hostpopsize = hostpopsize
        self.prev = prevalence
        self.dist = distance #this is distance from village 0
        self.h2b = hours2bite
        self.bpp = bitesPperson
        self.bn = bednets
        self.bnstr = bnstart
        self.bnstp = bnstop
        self.bncov = bncoverage
        self.mu = muTrans
        self.size = sizeTrans

#if __name__ == '__main__':
#    t = [Village(0, hostpopsize=1000, prevalence=0.2, distance=0),
#         Village(1, hostpopsize=200, prevalence=0.1, distance=1000)]

