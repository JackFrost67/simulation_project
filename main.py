from os import set_blocking
from numpy.lib.shape_base import column_stack
import pygame
import random
import sys
import numpy as np
from operator import attrgetter

WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
YELLOW = (255, 255, 0)
 
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

        if(color == YELLOW):
            self.PM = 100
            self.type = "SP"
        elif(color == GREEN):
            self.dir = self.setDirection()
            self.CHA = 100
            self.type = "NS"
        elif(color == RED):
            self.AA = 0
            self.type = "U"
        else:
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
    CAP2 = 0.01

    # consumption percentage of the chemoattractant
    CON = 0.95

    # parameter for the attraction to the chemoattractant
    PAP = 0.7

    # threshold of Physarum Mass that encapsulate a airplaneNS
    thresholdPM = 0.2

    def __init__(self):
        pygame.display.set_caption('Physarum Polycephalum Simulation') # setting name of the screen
        self._running = False
        self._rows = self._cols = 50
        self._size = (1000, 1000) 
        self.screen = pygame.display.set_mode(self._size)
        self._clock = pygame.time.Clock()
        self.m_position = (0, 0)
        self._all = pygame.sprite.Group()
        self._group = pygame.sprite.Group()
        self._user_group = pygame.sprite.GroupSingle()
        self._block = np.empty((self._rows, self._cols), dtype = Block)
        self._grid = np.empty((self._rows, self._cols), dtype = Cell)
        self._NS = []
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
                        arr[i][j].neighbors = [None, None, None, None, arr[i + 1][j], arr[i + 1][j -1], arr[i][j - 1], None]
                    else:
                        arr[i][j].neighbors = [arr[i - 1][j], arr[i - 1][j + 1], arr[i][j + 1], arr[i + 1][j + 1], arr[i + 1][j], None, None, None] 
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
                    arr[i][j].neighbors = [arr[i][j - 1], None, None, None, arr[i + 1][j], arr[i + 1][j - 1], arr[i][j - 1], arr[i - 1][j - 1]]
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
           for y in x:
                if y.type == "NS":
                    self._NS.append(y)

    def buildObstacle(self):
        for cellNS in self._NS:
            for neighbor in cellNS.neighbors:
                if neighbor != None:
                    neighbor.type = "U"
                    self._block[neighbor.index[0]][neighbor.index[1]].updateColor(RED)
                        
            if(cellNS.dir == "N" and cellNS.index[1] - 1 >= 0):
                self._grid[cellNS.index[0] - 1][cellNS.index[1]].type = "A"
                self._block[cellNS.index[0] - 1][cellNS.index[1]].updateColor(WHITE)
            elif(cellNS.dir == "W" and cellNS.index[0] - 1 >= 0):
                self._grid[cellNS.index[0]][cellNS.index[1] - 1].type = "A"
                self._block[cellNS.index[0]][cellNS.index[1] - 1].updateColor(WHITE)
            elif(cellNS.dir == "E" and cellNS.index[0] + 1 < self._rows):
                self._grid[cellNS.index[0]][cellNS.index[1] + 1].type = "A"
                self._block[cellNS.index[0]][cellNS.index[1] + 1].updateColor(WHITE)
            elif(cellNS.dir == "S" and cellNS.index[1] + 1 < self._cols):
                self._grid[cellNS.index[0] + 1][cellNS.index[1]].type = "A"
                self._block[cellNS.index[0] + 1][cellNS.index[1]].updateColor(WHITE)

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

                    pmVN = sum([((1 + N) * cell.neighbors[0].PM) - (cell.neighbors[0].AA * cell.PM) if cell.neighbors[0] is not None else 0,
                                ((1 + E) * cell.neighbors[2].PM) - (cell.neighbors[2].AA * cell.PM) if cell.neighbors[2] is not None else 0,
                                ((1 + S) * cell.neighbors[4].PM) - (cell.neighbors[4].AA * cell.PM) if cell.neighbors[4] is not None else 0,
                                ((1 + W) * cell.neighbors[6].PM) - (cell.neighbors[6].AA * cell.PM) if cell.neighbors[6] is not None else 0])
                    
                    pmMN = sum([((1 + NE) * cell.neighbors[1].PM) - (cell.neighbors[1].AA * cell.PM) if cell.neighbors[1] is not None else 0,
                                ((1 + SE) * cell.neighbors[3].PM) - (cell.neighbors[3].AA * cell.PM) if cell.neighbors[3] is not None else 0,
                                ((1 + SW) * cell.neighbors[5].PM) - (cell.neighbors[5].AA * cell.PM) if cell.neighbors[5] is not None else 0,
                                ((1 + NW) * cell.neighbors[7].PM) - (cell.neighbors[7].AA * cell.PM) if cell.neighbors[7] is not None else 0])
                    
                    cell.PM = cell.PM + self.PMP1 * (pmVN + self.PMP2 * pmMN)

                if(cell.type != "U" and cell.type != "NS"):
                    chaVN = sum([cell.neighbors[0].CHA - (cell.neighbors[0].AA * cell.PM) if cell.neighbors[0] is not None else 0,
                                cell.neighbors[2].CHA - (cell.neighbors[2].AA * cell.PM) if cell.neighbors[2] is not None else 0,
                                cell.neighbors[4].CHA - (cell.neighbors[4].AA * cell.PM) if cell.neighbors[4] is not None else 0,
                                cell.neighbors[6].CHA - (cell.neighbors[6].AA * cell.PM) if cell.neighbors[6] is not None else 0])
                    
                    chaMN = sum([cell.neighbors[1].CHA - (cell.neighbors[1].AA * cell.PM) if cell.neighbors[1] is not None else 0,
                                cell.neighbors[3].CHA - (cell.neighbors[3].AA * cell.PM) if cell.neighbors[3] is not None else 0,
                                cell.neighbors[5].CHA - (cell.neighbors[5].AA * cell.PM) if cell.neighbors[5] is not None else 0,
                                cell.neighbors[7].CHA - (cell.neighbors[7].AA * cell.PM) if cell.neighbors[7] is not None else 0])

                    cell.CHA = self.CON * (cell.CHA + self.CAP1 * (chaVN + self.CAP2 * chaMN))

                    if cell.CHA > 100:
                        cell.CHA = 100
                    elif cell.CHA < 0:
                        cell.CHA = 0
                
                if cell.PM != 0 and (cell.type != "NS" and cell.type != "U"):
                    alpha = np.interp(cell.PM, [0, 1000], [0, 100])
                    self._block[cell.index[0]][cell.index[1]].updateColor(YELLOW, alpha)

    def run(self):
        self._user = Block(self, YELLOW, (20, 20))
        self._user_group.add(self._user)
        self._all.add(self._user)

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
                    elif event.type == pygame.MOUSEBUTTONDOWN:
                        j = int(np.round(np.interp(self.m_position[0], [0, self._size[0]], [0, self._rows - 1])))
                        i = int(np.round(np.interp(self.m_position[1], [0, self._size[1]], [0, self._cols - 1])))
                        self._block[i][j].updateColor(self._user._color)
                        self._grid[i][j].updateMatrix(self._user._color)
                    elif event.type == pygame.KEYDOWN: # press enter to start simulation with the configuration 
                        if event.key == pygame.K_RETURN:
                            self.findNeighbors()
                            self.findNS()
                            self.buildObstacle()
                            self._running = True
                            self._user.image.set_alpha(0)
                        elif event.key == pygame.K_RIGHT:
                            if self._user._color == YELLOW:
                                self._user.updateColor(GREEN)
                            elif self._user._color == GREEN:
                                self._user.updateColor(RED)
                            elif self._user._color == RED:
                                self._user.updateColor(WHITE)
                            else:
                                self._user.updateColor(YELLOW)

                # update sprite
                self._user_group.update()
                self._all.draw(self.screen)
                
                # clock cap 60 ticks per seconds
                self._clock.tick(120)
                
                # update
                pygame.display.flip()
            
            else: #running the simulation
                for event in pygame.event.get():
                    if event.type == pygame.QUIT:
                        pygame.quit()
                        sys.exit()
                    elif event.type == pygame.KEYDOWN: # press enter to start simulation with the configuration 
                        if event.key == pygame.K_RETURN:
                            self._running = False
                            self._user.image.set_alpha(255)
                
                if (self._step % 50 != 0):
                    self.diffusionEquation()
                else:
                    pass

                # update time step
                self._step = self._step + 1

                # update sprite
                self._user_group.update()
                self._all.draw(self.screen)

                # clock cap 60 ticks per seconds
                self._clock.tick(60)
                pygame.display.flip()
            
if __name__ == "__main__":
    pygame.init()
    simulation = Simulation()
    simulation.setGrid()
    simulation.run()
    pygame.quit()
    sys.exit()