# TODO
# define default values

class Cell:
    def __init__(self, ID = None, AA = False, PM = 0, CHA = 0, TE = False, type = "A"):
        self.ID = ID
        self.AA = AA
        self.PM = PM
        self.CHA = CHA
        self.TE = TE
        self.type = type # SP NS U A
