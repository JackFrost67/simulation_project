from numpy.core.fromnumeric import mean, var
import pygame
import random
import sys
import numpy as np
from operator import attrgetter

# color palette
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
YELLOW = (255, 255, 0)
BLUE = (0, 0, 255)
 
# This sets the margin between each cell
MARGIN = 1

class Cell():
    def __init__(self, ID = None, AA = 1, PM = 0, CHA = 0, TE = False, type = "A", index = None, dir = None):
        self.ID = ID
        self.AA = AA
        self.PM = PM
        self.CHA = CHA
        self.TE = TE
        self.type = type # SP NS U A
        self.index = index
        self.dir = dir #N S W E
        self.neighbors = np.array([])

    def setDirection(self):
        irand = random.randint(0, 3)
        if (irand == 0):
            return "N"
        elif (irand == 1):
            return "W"
        elif (irand == 2):
            return "S"
        elif (irand == 3):
            return "E"
        
    def updateMatrix(self, color):
        # Create a tuple with the new color and assign it.
        self.CHA = 0
        self.PM = 0
        self.AA = 1
        self.dir = None

        if(color == BLUE):
            self.PM = Simulation.defaultPM
            self.type = "SP"
        elif(color == GREEN):
            self.dir = self.setDirection()
            self.CHA = Simulation.defaultCHA
            self.type = "NS"
        elif(color == RED):
            self.AA = 0
            self.type = "U"
        elif(color == WHITE):
            self.type = "A"

class Block(pygame.sprite.Sprite):
    # Constructor. Pass in the color of the block,
    # and its x and y position
    def __init__(self, simulation, color, size = (50, 50), position = (0, 0)):
        # Call the parent class (Sprite) constructor
        super().__init__()

        self._simulation = simulation
        self._color = color

        # Create an image of the block, and fill it with a color.
        # This could also be an image loaded from the disk.
        self.image = pygame.Surface(size)
        self.image.fill(self._color)

        # Fetch the rectangle object that has the dimensions of the image
        # Update the position of this object by setting the values of rect.x and rect.y
        self.rect = self.image.get_rect(topleft = position)

    def update(self):
        self.rect.x = self._simulation.m_position[0]
        self.rect.y = self._simulation.m_position[1]

    def updateColor(self, color, alpha = 255):
        self._color = color
        self.image.fill(color)
        self.image.set_alpha(alpha)

class Simulation():
    ## parameter for the pphysarum simulation
    # parameters for diffusion equation for the cytoplasm
    PMP1 = 0.08
    PMP2 = 0.01

    # parameters for the diffusion of the chemoattractant
    CAP1 = 0.05
    CAP2 = 0.02

    # consumption percentage of the chemoattractant
    CON = 0.95

    # parameter for the attraction to the chemoattractant
    PAP = 0.2

    # threshold of Physarum Mass that encapsulate a NS
    thresholdPM = 0.7

    defaultCHA = 200
    defaultPM = 200

    def __init__(self):
        pygame.display.set_caption('Physarum Polycephalum Simulation') # setting name of the screen
        pygame.mouse.set_visible(0)

        self._running = False
        self._done = False
        self._rows = self._cols = 35
        self._totCHA = 1
        self._size = (800, 800) 
        self.screen = pygame.display.set_mode(self._size)
        self._clock = pygame.time.Clock()
        self.m_position = (pygame.mouse.get_pos()[0], pygame.mouse.get_pos()[1])
        self._all = pygame.sprite.Group()
        self._group = pygame.sprite.Group()
        self._user_group = pygame.sprite.GroupSingle()
        self._block = np.empty((self._rows, self._cols), dtype = Block)
        self._grid = np.empty((self._rows, self._cols), dtype = Cell)
        self._NS = []
        self._SP = []
        self._step = 0

        for row in range(self._rows):
            for col in range(self._cols):
                self._grid[row][col] = Cell(index = (row, col))

    def findNeighbors(self):
        N = self._rows
        M = self._cols

        arr = self._grid
        
        for i in range(N):
            for j in range(M):
                if i == 0:
                    if j == 0:
                        arr[i][j].neighbors = [None, None, arr[i][j + 1], arr[i + 1][j + 1], arr[i + 1][j], None, None, None] 
                    elif j == M - 1:
                        arr[i][j].neighbors = [None, None, None, None, arr[i + 1][j], arr[i + 1][j - 1], arr[i][j - 1], None]
                    else:
                        arr[i][j].neighbors = [None, None, arr[i][j + 1], arr[i + 1][j + 1], arr[i + 1][j], arr[i + 1][j - 1], arr[i][j - 1], None]
                elif i == N - 1:
                    if j == 0:
                        arr[i][j].neighbors = [arr[i - 1][j], arr[i - 1][j + 1], arr[i][j + 1], None, None, None, None, None]
                    elif j == M - 1:
                        arr[i][j].neighbors = [arr[i - 1][j], None, None, None, None, None, arr[i][j - 1], arr[i - 1][j - 1]]
                    else:
                        arr[i][j].neighbors = [arr[i - 1][j], arr[i - 1][j + 1], arr[i][j + 1], None, None, None, arr[i][j - 1], arr[i - 1][j - 1]] 
                elif j == 0:
                    arr[i][j].neighbors = [arr[i - 1][j], arr[i - 1][j + 1], arr[i][j + 1], arr[i + 1][j + 1], arr[i + 1][j], None, None, None]
                elif j == M - 1: 
                    arr[i][j].neighbors = [arr[i - 1][j], None, None, None, arr[i + 1][j], arr[i + 1][j - 1], arr[i][j - 1], arr[i - 1][j - 1]]
                else:
                    arr[i][j].neighbors = [arr[i - 1][j], arr[i - 1][j + 1], arr[i][j + 1], arr[i + 1][j + 1], arr[i + 1][j], arr[i + 1][j - 1], arr[i][j - 1], arr[i - 1][j - 1]]

    def setGrid(self):
        pixelSize = (((self._size[0] - (MARGIN * (self._rows + 1))) / self._rows), ((self._size[1] - (MARGIN * (self._cols + 1))) / self._cols))

        for row in range(self._rows):
            for column in range(self._cols):
                self._block[row][column] = Block(self,
                            WHITE,
                            pixelSize,
                            ((MARGIN + pixelSize[0]) * column + MARGIN, 
                            (MARGIN + pixelSize[1]) * row + MARGIN))
                self._group.add(self._block[row][column])
                self._all.add(self._block[row][column])
    
    def findNS(self):
       for x in self._grid:
           for cell in x:
                if cell.type == "NS":
                    self._NS.append(cell)
    
    def findSP(self):
       for x in self._grid:
           for cell in x:
                if cell.type == "SP":
                    self._SP.append(cell)

    def buildObstacle(self):
        return
        # for cellNS in self._NS:
        #     for neighbor in cellNS.neighbors:
        #         if neighbor != None:
        #             neighbor.type = "U"
        #             self._block[neighbor.index[0]][neighbor.index[1]].updateColor(RED)
                        
        #     if(cellNS.dir == "N" and cellNS.index[1] - 1 >= 0):
        #         self._grid[cellNS.index[0] - 1][cellNS.index[1]].type = self._grid[cellNS.index[0] - 1][cellNS.index[1] - 1].type = self._grid[cellNS.index[0] - 1][cellNS.index[1] + 1].type = "A"
        #         self._block[cellNS.index[0] - 1][cellNS.index[1]].updateColor(WHITE)
        #         self._block[cellNS.index[0] - 1][cellNS.index[1] - 1].updateColor(WHITE)
        #         self._block[cellNS.index[0] - 1][cellNS.index[1] + 1].updateColor(WHITE)
        #     elif(cellNS.dir == "W" and cellNS.index[0] - 1 >= 0):
        #         self._grid[cellNS.index[0]][cellNS.index[1] - 1].type = self._grid[cellNS.index[0] - 1][cellNS.index[1] - 1].type = self._grid[cellNS.index[0] + 1][cellNS.index[1] - 1].type = "A"
        #         self._block[cellNS.index[0]][cellNS.index[1] - 1].updateColor(WHITE)
        #         self._block[cellNS.index[0] - 1][cellNS.index[1] - 1].updateColor(WHITE)
        #         self._block[cellNS.index[0] + 1][cellNS.index[1] - 1].updateColor(WHITE)
        #     elif(cellNS.dir == "E" and cellNS.index[0] + 1 < self._rows):
        #         self._grid[cellNS.index[0]][cellNS.index[1] + 1].type =  self._grid[cellNS.index[0] - 1][cellNS.index[1] + 1].type = self._grid[cellNS.index[0] + 1][cellNS.index[1] + 1].type = "A"
        #         self._block[cellNS.index[0]][cellNS.index[1] + 1].updateColor(WHITE)
        #         self._block[cellNS.index[0] - 1][cellNS.index[1] + 1].updateColor(WHITE)
        #         self._block[cellNS.index[0] + 1][cellNS.index[1] + 1].updateColor(WHITE)
        #     elif(cellNS.dir == "S" and cellNS.index[1] + 1 < self._cols):
        #         self._grid[cellNS.index[0] + 1][cellNS.index[1]].type = self._grid[cellNS.index[0] + 1][cellNS.index[1] - 1].type = self._grid[cellNS.index[0] + 1][cellNS.index[1] + 1].type = "A"
        #         self._block[cellNS.index[0] + 1][cellNS.index[1]].updateColor(WHITE)
        #         self._block[cellNS.index[0] + 1][cellNS.index[1] - 1].updateColor(WHITE)
        #         self._block[cellNS.index[0] + 1][cellNS.index[1] + 1].updateColor(WHITE)

    def diffusionEquation(self):
        for x in self._grid:
            for cell in x:
                if(cell.type != "U" and cell.type != "SP"):
                    maxCHAcell = max((x for x in cell.neighbors if x != None), key = attrgetter("CHA"))

                    i = cell.index[0]
                    j = cell.index[1]

                    N = S = W = E = NW = NE = SW = SE = 0

                    if (maxCHAcell.index == (i - 1, j)):
                        W = self.PAP
                        E = -self.PAP
                    elif (maxCHAcell.index == (i, j - 1)):
                        S = self.PAP
                        N = -self.PAP
                    elif (maxCHAcell.index == (i + 1, j)):
                        E = self.PAP
                        W = -self.PAP
                    elif(maxCHAcell.index == (i, j + 1)):
                        N = self.PAP
                        S = -self.PAP
                    elif (maxCHAcell.index == (i - 1, j - 1)):
                        SW = self.PAP
                        NE = -self.PAP
                    elif (maxCHAcell.index == (i + 1, j - 1)):
                        SE = self.PAP
                        NW = -self.PAP
                    elif (maxCHAcell.index == (i - 1, j + 1)):
                        NW = self.PAP
                        SE = -self.PAP
                    elif (maxCHAcell.index == (i + 1, j + 1)):
                        NE = self.PAP
                        SW = -self.PAP

                    pmVN = sum([(((1 + W) * cell.neighbors[0].PM) - cell.PM) if cell.neighbors[0] is not None and cell.neighbors[0].AA else 0,
                                (((1 + N) * cell.neighbors[2].PM) - cell.PM) if cell.neighbors[2] is not None and cell.neighbors[2].AA else 0,
                                (((1 + E) * cell.neighbors[4].PM) - cell.PM) if cell.neighbors[4] is not None and cell.neighbors[4].AA else 0,
                                (((1 + S) * cell.neighbors[6].PM) - cell.PM) if cell.neighbors[6] is not None and cell.neighbors[6].AA else 0])
                     
                    pmMN = sum([(((1 + NW) * cell.neighbors[1].PM) - cell.PM) if cell.neighbors[1] is not None and cell.neighbors[1].AA else 0,
                                (((1 + NE) * cell.neighbors[3].PM) - cell.PM) if cell.neighbors[3] is not None and cell.neighbors[3].AA else 0,
                                (((1 + SE) * cell.neighbors[5].PM) - cell.PM) if cell.neighbors[5] is not None and cell.neighbors[5].AA else 0,
                                (((1 + SW) * cell.neighbors[7].PM) - cell.PM) if cell.neighbors[7] is not None and cell.neighbors[7].AA else 0])
                    #prop = cell.CHA/self._totCHA 

                    cell.PM = cell.PM + ((self.PMP1 * pmVN) + (self.PMP2 * pmMN))
                    
                    minCHAcell = min((x for x in cell.neighbors if (x != None and x.TE != True)), key = attrgetter("CHA"))

                    """ L'idea iniziale era quella di togliere PM alla cella vicina con il CHA minimo per darlo alla cella
                    corrente. Dato che non funziona ho provato a cercare il minimo della moltiplicazione tra PM e CHA di ogni
                    vicino. Anche questo non funziona.
                    Risultato interessante perchè il PM, nel idea iniziale con un beta superiore a 0.05 , si vede che va proprio
                    nella direzione del NS.  
                     """
                    if (minCHAcell.CHA != 0 and minCHAcell.CHA < cell.CHA and minCHAcell.PM != self.defaultPM and minCHAcell.PM != 0 and minCHAcell is not None):
                        beta = 0.01
                        givePM = minCHAcell.PM * beta
                        minCHAcell.PM = minCHAcell.PM - givePM
                        cell.PM = cell.PM + givePM

                    if cell.PM > self.defaultPM:
                        cell.PM = self.defaultPM
                    elif cell.PM < 0:
                        cell.PM = 0

                if(cell.type != "U" and cell.type != "NS"):
                    chaVN = sum([(cell.neighbors[0].CHA - cell.CHA) if cell.neighbors[0] is not None and cell.neighbors[0].AA else 0,
                                (cell.neighbors[2].CHA - cell.CHA) if cell.neighbors[2] is not None and cell.neighbors[2].AA else 0,
                                (cell.neighbors[4].CHA - cell.CHA) if cell.neighbors[4] is not None and cell.neighbors[4].AA else 0,
                                (cell.neighbors[6].CHA - cell.CHA) if cell.neighbors[6] is not None and cell.neighbors[6].AA else 0])
                    
                    chaMN = sum([(cell.neighbors[1].CHA - cell.CHA) if cell.neighbors[1] is not None and cell.neighbors[1].AA else 0,
                                (cell.neighbors[3].CHA - cell.CHA) if cell.neighbors[3] is not None and cell.neighbors[3].AA else 0,
                                (cell.neighbors[5].CHA - cell.CHA) if cell.neighbors[5] is not None and cell.neighbors[5].AA else 0,
                                (cell.neighbors[7].CHA - cell.CHA) if cell.neighbors[7] is not None and cell.neighbors[7].AA else 0])

                    cell.CHA = (self.CON * cell.CHA) + ((self.CAP1 * chaVN) + (self.CAP2 * chaMN))

                    if cell.CHA > self.defaultCHA:
                        cell.CHA = self.defaultCHA
                    elif cell.CHA < 0:
                        cell.CHA = 0
                
                if cell.PM != 0 and cell.TE != True and cell.type != "SP" and cell.type != "NS" and cell.type != "U":
                    self._block[cell.index[0]][cell.index[1]].updateColor(YELLOW, cell.PM)

    def diffusionCHA(self):
        totCHA = 0
        for x in self._grid:
            for cell in x:
                if(cell.type != "U" and cell.type != "NS"):
                    chaVN = sum([(cell.neighbors[0].CHA - cell.CHA) if cell.neighbors[0] is not None and cell.neighbors[0].AA else 0,
                                (cell.neighbors[2].CHA - cell.CHA) if cell.neighbors[2] is not None and cell.neighbors[2].AA else 0,
                                (cell.neighbors[4].CHA - cell.CHA) if cell.neighbors[4] is not None and cell.neighbors[4].AA else 0,
                                (cell.neighbors[6].CHA - cell.CHA) if cell.neighbors[6] is not None and cell.neighbors[6].AA else 0])
                    
                    chaMN = sum([(cell.neighbors[1].CHA - cell.CHA) if cell.neighbors[1] is not None and cell.neighbors[1].AA else 0,
                                (cell.neighbors[3].CHA - cell.CHA) if cell.neighbors[3] is not None and cell.neighbors[3].AA else 0,
                                (cell.neighbors[5].CHA - cell.CHA) if cell.neighbors[5] is not None and cell.neighbors[5].AA else 0,
                                (cell.neighbors[7].CHA - cell.CHA) if cell.neighbors[7] is not None and cell.neighbors[7].AA else 0])

                    cell.CHA = (self.CON * cell.CHA) + ((self.CAP1 * chaVN) + (self.CAP2 * chaMN))
                    totCHA = totCHA + cell.CHA

                    if cell.CHA > self.defaultCHA:
                        cell.CHA = self.defaultCHA
                    elif cell.CHA < 0:
                        cell.CHA = 0
    
    def cleanCHA(self):
        for x in self._grid:
            for cell in x:
                if (cell.CHA < 1.1e-12):
                    cell.CHA = 0    
    
    def setTE(self):
        for cellNS in self._NS:
            if cellNS.PM >= self.thresholdPM :
                cell = cellNS
                count = 0
                #Aggiunto cellTE per cercare di evitare i loop
                cellTE = []
                while (cell.type != "SP" and count <= 500):
                    cell.TE = True
                    cellTE.append(cell)
                    self._block[cell.index[0]][cell.index[1]].updateColor(BLACK)

                    cell = max((x for x in cell.neighbors if (x != None and x not in cellTE)), key = attrgetter("PM"))
                    count = count + 1
                
                self._NS.remove(cellNS)
                cellNS.CHA = 0
                cellNS.PM = self.defaultPM
                cellNS.type = "SP"
                self._SP.append(cellNS)

                if not self._NS:
                    self._done = True

    def run(self):        
        self._user = Block(self, BLUE, (20, 20))
        self._user_group.add(self._user)
        self._all.add(self._user)

        self._user_group.update()
        self._all.draw(self.screen)

        # main loop
        while True:
            if not self._running: # config phase
                # fill screen 
                self.screen.fill(BLACK)

                # handling events
                for event in pygame.event.get():
                    if event.type == pygame.QUIT:
                        pygame.quit()
                        sys.exit()
                    elif event.type == pygame.MOUSEMOTION:
                        self.m_position = (pygame.mouse.get_pos()[0], pygame.mouse.get_pos()[1])
                    elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                        j = int(np.round(np.interp(self.m_position[0], [0, self._size[0]], [0, self._rows])))
                        i = int(np.round(np.interp(self.m_position[1], [0, self._size[1]], [0, self._cols])))
                        self._block[i][j].updateColor(self._user._color)
                        self._grid[i][j].updateMatrix(self._user._color)
                    elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 3:
                        j = int(np.round(np.interp(self.m_position[0], [0, self._size[0]], [0, self._rows])))
                        i = int(np.round(np.interp(self.m_position[1], [0, self._size[1]], [0, self._cols])))
                        self._block[i][j].updateColor(WHITE)
                        self._grid[i][j].updateMatrix(WHITE)
                    elif event.type == pygame.KEYDOWN: # press enter to start simulation with the configuration 
                        if event.key == pygame.K_RETURN:
                            self.findNeighbors()
                            self.findNS()
                            self.findSP()
                            self.buildObstacle()
                            self._running = True
                            self._user.image.set_alpha(0)
                        elif event.key == pygame.K_RIGHT:
                            if self._user._color == BLUE:
                                self._user.updateColor(GREEN)
                            elif self._user._color == GREEN:
                                self._user.updateColor(RED)
                            elif self._user._color == RED:
                                self._user.updateColor(BLUE)

                # update sprite
                self._user_group.update()
                self._all.draw(self.screen)
                
                # clock cap 60 ticks per seconds
                self._clock.tick(60)
                
                # update
                pygame.display.flip()
            
            elif not self._done: #running the simulation
                pause = True
                for event in pygame.event.get():
                    if event.type == pygame.QUIT:
                        pygame.quit()
                        sys.exit()
                    elif event.type == pygame.KEYDOWN: # press enter to start simulation with the configuration 
                        if event.key == pygame.K_RETURN:
                            self._running = False
                            self._user.image.set_alpha(255)
                            self.__init__()
                            self.setGrid()
                            self.run()
                        elif event.key == pygame.K_s:
                            while pause:
                                for event in pygame.event.get():
                                    if event.key == pygame.K_s:
                                        pause = False
                
                print(self._step)
                # i primi 100 passi facciamo diffondere il CHA
                if(self._step < 100):
                    self.diffusionCHA()

                    if (self._step == 99):
                        self.cleanCHA()
                        tmp = []
                        for i in range(self._rows):
                            tmp.append([])
                            for j in range(self._cols):
                                tmp[i].append(self._grid[i][j].CHA) 
                        print(self._step)

                elif (self._step > 100):
                    self.diffusionEquation()
                elif (self._step == 5000):
                    cellSP = self._SP[len(self._SP) - 2]
                    self._SP.remove(cellSP)

                    for cell in self._SP:
                        cell.type = "NS"
                        cell.PM = 0
                        cell.CHA = self.defaultCHA
                        self._NS.append(cell)

                    self._SP = [cellSP]
                
                if(self._step % 50 == 0):
                    self.setTE()
                        
                if (self._step >= 5000):
                    if (self._step >= 10000): 
                        raise NameError('Too many iterations, giving up')

                # update time step
                self._step = self._step + 1

                # update sprite
                self._user_group.update()
                self._all.draw(self.screen)

                # clock cap n ticks per seconds
                self._clock.tick(120)
                pygame.display.flip()
            
            else:
                for event in pygame.event.get():
                    if event.type == pygame.QUIT:
                        pygame.quit()
                        sys.exit()
                    elif event.type == pygame.KEYDOWN: # press enter to start simulation with the configuration 
                        if event.key == pygame.K_RETURN:
                            self._running = False
                            self._user.image.set_alpha(255)
                            self.__init__()
                            self.setGrid()
                            self.run()
            
if __name__ == "__main__":
    pygame.init()
    simulation = Simulation()
    simulation.setGrid()
    simulation.run()
    pygame.quit()
    sys.exit()
