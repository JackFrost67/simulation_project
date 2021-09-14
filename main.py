# E = A u U, Entire Area
# A = Available Area, U = Unavailable Area
# S = Place where physarium is introduced
# N = Nutrient Source

# State at time:
# ST = [ AA, PM, CHA, TE]
# where:
#   - AA = Available Area
#   - PM = Physarum Mass
#   - CHA = Chemo Attractant
#   - TE = Tube Existence

# TODO:
# init of parameters
# Apply diffusion equation for certian amount of steps
# Check in any of the N are cover with a certain amount of mould

import sys
import numpy as np
import pygame
import random
from pygame.constants import CONTROLLER_AXIS_INVALID
from pygame.sprite import DirtySprite
from Cell import Cell

size = (width, height) = 300, 300

pygame.init()

screen = pygame.display.set_mode(size, pygame.RESIZABLE) # setting screen size

pygame.display.set_caption('Physarum Polycephalum Simulation') # setting name of the screen

clock = pygame.time.Clock() # setting clock

squareSize = 10
cols, rows = int(screen.get_width() / squareSize), int(screen.get_height() / squareSize)

colorW = (255, 255, 255)
colorY = (255, 255, 0)
colorR = (255, 0, 0)
colorG = (0, 255, 0)
colorB = (0, 0, 0)
grid = [] # empty grid

done = False # false until configuration are done
indexSP = 0

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

# threshold of Physarum Mass that encapsulate a NS
thresholdPM = 0.2

def clip(value, min_, max_):
    value_ = [value]
    value_ = np.clip(value_, min_, max_)
    return value_[0]

def getCHA(i, j):
    if (i < 0 or i >= rows or j < 0 or j >= cols or grid[i * cols + j][2].type == "U"):
        return 0
    else:
        return grid[i * cols + j][2].CHA

def getPM(i, j):
    if (i < 0 or i >= rows or j < 0 or j >= cols or grid[i * cols + j][2].type == "U"):
        return 0
    else:
        return grid[i * cols + j][2].PM

def getAA(i, j):
    if (i < 0 or i >= rows or j < 0 or j >= cols or grid[i * cols + j][2].type == "U"):
        return 0
    else:
        return grid[i * cols + j][2].AA

def setTE(i, j, countLoop = 0):    
    while (grid[i * cols + j][2].type != "SP" and countLoop < 500):  
        grid[i * cols + j][2].TE = True
        grid[i * cols + j] = (grid[i * cols + j][0], colorB, grid[i * cols + j][2])
        
        maxPM = 0
        for x in range(i - 1, i + 2):
            for y in range (j - 1, j + 2):
                if not(x == i and y == j):
                    if (getPM(x, y) > maxPM):
                        l = x
                        m = y 
                        maxPM = getPM(x, y)
        i = l
        j = m 
        countLoop = countLoop + 1

    if(grid[i * cols + j][2].type == "SP"):
        grid[i * cols + j][2].TE = True
        return True

    if (countLoop == 500):
        return False


def diffusion_equation():
    # Now draw the rects. You can unpack the tuples
    # again directly in the head of the for loop.
    for rect, color, cell in grid:
        pygame.draw.rect(screen, color, rect)

    for i in range(rows):
        for j in range(cols):
            if (grid[i * cols + j][2].type != "U" and grid[i * cols + j][2].type != "SP"):
                valuesCHA = [
                    getCHA(i - 1, j - 1), # 0
                    getCHA(i, j - 1), # 1
                    getCHA(i + 1, j - 1), # 2
                    getCHA(i - 1, j), # 3 
                    getCHA(i + 1, j), # 4
                    getCHA(i - 1, j + 1), # 5
                    getCHA(i, j + 1), # 6
                    getCHA(i + 1, j + 1) # 7
                ]

                maxCHA = max(valuesCHA)

                N = S = W = E = NW = NE = SW = SE = 0
                
                if (maxCHA == getCHA(i - 1, j)):
                    W = PAP
                    E = -PAP
                elif (maxCHA == getCHA(i, j - 1)):
                    S = PAP
                    N = -PAP
                elif (maxCHA == getCHA(i + 1, j)):
                    E = PAP
                    W = -PAP
                elif(maxCHA == getCHA(i, j + 1)):
                    N = PAP
                    S = -PAP
                elif (maxCHA == getCHA(i - 1, j - 1)):
                    SW = PAP
                    NE = -PAP
                elif (maxCHA == getCHA(i + 1, j - 1)):
                    SE = PAP
                    NW = -PAP
                elif (maxCHA == getCHA(i - 1, j + 1)):
                    NW = PAP
                    SE = -PAP
                elif (maxCHA == getCHA(i + 1, j + 1)):
                    NE = PAP
                    SW = -PAP

                pmVN = (((1 + W) * getPM(i - 1, j)) - (getAA(i - 1, j) * getPM(i, j))
                    + ((1 + S) * getPM(i, j - 1)) - (getAA(i, j - 1) * getPM(i, j))
                    + ((1 + E) * getPM(i + 1, j)) - (getAA(i + 1, j) * getPM(i, j))
                    + ((1 + N) * getPM(i, j + 1)) - (getAA(i, j + 1) * getPM(i, j))
                )
                pmMN = (((1 + SW) * getPM(i - 1, j - 1)) - (getAA(i - 1, j - 1) * getPM(i, j))
                    + ((1 + SE) * getPM(i + 1, j - 1))- (getAA(i + 1, j - 1) * getPM(i, j))
                    + ((1 + NW) * getPM(i - 1, j + 1)) - (getAA(i - 1, j + 1) * getPM(i, j))
                    + ((1 + NE) * getPM(i + 1, j + 1)) - (getAA(i + 1, j + 1) * getPM(i, j))
                )
                
                grid[i * cols + j][2].PM = getPM(i, j) + PMP1 * (pmVN + PMP2 * pmMN)

            if(grid[i * cols + j][2].type != "U" and grid[i * cols + j][2].type != "NS"):
                
                chaVN = ((getCHA(i - 1, j) - (getAA(i - 1, j) * getCHA(i, j)))
                    + (getCHA(i, j - 1) - (getAA(i, j - 1) * getCHA(i, j)))
                    + (getCHA(i + 1, j) - (getAA(i + 1, j) * getCHA(i, j)))
                    + (getCHA(i, j + 1) - (getAA(i, j + 1) * getCHA(i, j)))
                )
                chaMN = ((getCHA(i - 1, j - 1) - (getAA(i - 1, j - 1) * getCHA(i, j)))
                    + (getCHA(i + 1, j - 1) - (getAA(i + 1, j - 1) * getCHA(i, j)))
                    + (getCHA(i - 1, j + 1) - (getAA(i - 1, j + 1) * getCHA(i, j)))
                    + (getCHA(i + 1, j + 1) - (getAA(i + 1, j + 1) * getCHA(i, j)))
                )          

                grid[i * cols + j][2].CHA = CON * (getCHA(i, j) + CAP1 * (chaVN + CAP2 * chaMN))
                
                if(grid[i * cols + j][2].CHA > 100):
                    grid[i * cols + j][2].CHA = 100
                elif(grid[i * cols + j][2].CHA < 0):
                    grid[i * cols + j][2].CHA = 0
            
            (rect, color, cell) = grid[i * cols + j]

            if(cell.PM != 0 and (cell.type != "NS" and cell.type != "U")):
                alpha = clip((int)(cell.PM), 0, 255)
                color = (255, 255, 255 - alpha)
                grid[i * cols + j] = (rect, color, cell)
            
            #print(grid[indexSP + 1][2].PM)

            pygame.display.flip()

def simulation():

    #TestAA
    testAA=[]
    for i in range(rows):
        testAA.append([])
        for j in range(cols):
            testAA[i].append(grid[i*cols +j][2].AA)

    #Save all NS cell
    cellNS = []
    cellSP = []

    for (_, _ , cell) in grid:
        if (cell.type == "NS"):
            cellNS.append(cell) 

        if (cell.type == "SP"):
            cellSP.append(cell) 

    #Start the simulation
    t = 1
    keepOn = True 
    #lastTwoCellNS[0] last cell
    #lastTwoCellNS[1] second last cell
    lastTwoCellNS = [None, None]
    while (keepOn):
        #for 50 times we compute the equation
        if (t % 50 != 0):
            diffusion_equation()
        else:
        #if list of foods is empty, stop the simulation4
            for cell in cellNS:
                i = cell.index[0]
                j = cell.index[1]
                
                if (cell.type == "NS" and cell.PM >= thresholdPM):
                    keepOn = setTE(i, j)
                    cell.type = "SP"
                    cell.CHA = 0
                    cell.PM = 100
                    cellSP.append(cell)
                    cellNS.remove(cell)
                    lastTwoCellNS[1] = lastTwoCellNS[0]
                    lastTwoCellNS[0] = cell

            TestPM = []
            for a in range(rows):
                TestPM.append([])
                for b in range(cols):
                    TestPM[a].append(grid[a * cols + b][2].PM)

            if (not cellNS):
                return
                    
            if (t >= 5000):
                if (t >= 10000):  
                    return

                if(t == 5000): 
                    for c in cellSP:
                        cell.type = "NS"
                        cell.CHA = 100
                        cellNS.append(cell)
                        cellSP.remove(cell)
                    
                    lastTwoCellNS[1].type = "SP"
                    lastTwoCellNS[1].PM = 100
                    lastTwoCellNS[1].CHA = 0
                    cellSP.add(lastTwoCellNS[1])

                    if (lastTwoCellNS[1] in  cellNS):
                        cellNS.remove(lastTwoCellNS[1])

        print("t", t)
        t = t + 1    
#TODO implemnt the function for build obstacle
#We base this method on the assumption that there is only one SP
def buildObstacle(cellNS, cellSP):
    if(len(cellSP)!= 1):
        return 
    i = cellSP[0].index[0]
    j = cellSP[0].index[1]
    for c in cellNS:
        x = c.index[0]
        y = c.index[1]
        if (c.dir == "N" and x <= i):
            defineObstacle(x + 1, y + 1)
            defineObstacle(x + 1, y)
            defineObstacle(x + 1, y - 1)
            defineObstacle(x    , y + 1)
            defineObstacle(x    , y - 1)
        if (c.dir == "S" and x >= i):
            defineObstacle(x - 1, y + 1)
            defineObstacle(x - 1, y)
            defineObstacle(x - 1, y - 1)
            defineObstacle(x    , y + 1)
            defineObstacle(x    , y - 1)
        if (c.dir == "E" and y >= j):
            defineObstacle(x + 1, y - 1)
            defineObstacle(x    , y - 1)
            defineObstacle(x - 1, y - 1)
            defineObstacle(x + 1, y)
            defineObstacle(x - 1, y)
        if (c.dir == "W" and y <= j):
            defineObstacle(x + 1, y - 1)
            defineObstacle(x    , y - 1)
            defineObstacle(x - 1, y - 1)
            defineObstacle(x + 1, y)
            defineObstacle(x - 1, y)


        
def defineObstacle(i, j):
    index = i * cols + j
    grid[index][2].type = "U"
    grid[index][2].AA = 0
    (rect, _, cell) = grid[index]
    grid[index] = (rect, colorR, cell)


def initializationDirection():
    irand = random.randint(0,3)
    if (irand == 0):
        return "N"
    elif (irand == 1):
        return "W"
    elif (irand == 2):
        return "S"
    elif (irand == 3):
        return "E"

if __name__ == "__main__":
    # generating empty grid
    for y in range(rows):
        for x in range(cols):
            dir = initializationDirection()
            rect = pygame.Rect(x * (squareSize + 1), y * (squareSize + 1), squareSize, squareSize)
            grid.append((rect, colorW, Cell(index=(y,x), dir = dir)))

    simulationOn = True
    while simulationOn:
        for event in pygame.event.get():
            if event.type == pygame.QUIT: # event to exit the simulation
                simulationOn = False
                pygame.quit()
                sys.exit()
            elif event.type == pygame.KEYDOWN: # press enter to start simulation with the configuration 
                if event.key == pygame.K_RETURN:
                    if done:
                        done = False
                    else:
                        done = True

        # Now draw the rects. You can unpack the tuples
        # again directly in the head of the for loop.
        for rect, color, cell in grid:
            pygame.draw.rect(screen, color, rect)

        # draw source, food and other stuff
        if pygame.mouse.get_pressed()[0] and not(done):
            mousePos = pygame.mouse.get_pos()
            for index, (rect, color, cell) in enumerate(grid):
                if rect.collidepoint(mousePos):
                    # Create a tuple with the new color and assign it.
                    cell.CHA = 0
                    cell.PM = 0
                    cell.AA = 1
                    if color == colorY:
                        cell.type = "U" # not available
                        cell.AA = 0
                        grid[index] = (rect, colorR, cell)
                    elif color == colorR:
                        cell.type = "NS" # food
                        cell.CHA = 100
                        grid[index] = (rect, colorG, cell)
                    elif color == colorG:
                        cell.type = "A" # avalaible
                        grid[index] = (rect, colorW, cell)
                    else:
                        cell.type = "SP" # starting point
                        cell.PM = 100
                        indexSP = index
                        grid[index] = (rect, colorY, cell)
        
        if done:
            #TODO do a method for  this type of istruction
            cellNS = []
            cellSP = []
            for (_, _ , c) in grid:
                if (c.type == "NS"):
                    cellNS.append(c) 
                if (c.type == "SP"):
                    cellSP.append(c) 

            buildObstacle(cellNS, cellSP)
            pygame.display.flip()

            simulation() # start simulation
            
            for index, (rect, color, cell) in enumerate(grid):
                if(cell.TE == True):
                    grid[index] = (rect, colorB, cell)
            
        done = False

        pygame.display.flip()
        clock.tick(10)