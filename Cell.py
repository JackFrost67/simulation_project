class Cell:
    def __init__(self, ID = None, AA = 0, PM = 0, CHA = 0, TE = False, type = "A", index = (-1, -1)):
        self.ID = ID
        self.AA = AA
        self.PM = PM
        self.CHA = CHA
        self.TE = TE
        self.type = type # SP NS U A
        self.index = index
