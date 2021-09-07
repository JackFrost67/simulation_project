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

import pygame
import sys
import numpy as np
from pygame.constants import CONTROLLER_AXIS_INVALID
from Cell import Cell

size = (width, height) = 600, 600

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
PAP = 0.8

# threshold of Physarum Mass that encapsulate a NS
thresholdPM = 500

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

def setTE(i, j):
    alpha0 = alpha1 = alpha2 = alpha3 = alpha4 = alpha5 = alpha6= alpha7 = 1
    while(grid[i * cols + j][2].type != "SP"):
        valuesPM = [
            getPM(i - 1, j - 1) * alpha0, # 0
            getPM(i, j - 1) * alpha1, # 1
            getPM(i + 1, j - 1) * alpha2, # 2
            getPM(i - 1, j) * alpha3, # 3 
            getPM(i + 1, j) * alpha4, # 4
            getPM(i - 1, j + 1) * alpha5, # 5
            getPM(i, j + 1) * alpha6, # 6
            getPM(i + 1, j + 1) * alpha7 # 7
        ]         
        
        grid[i * cols + j][2].TE = True
        grid[i * cols + j] = (grid[i * cols + j][0], colorB, grid[i * cols + j][2])
        maxPMindex = valuesPM.index(max(valuesPM))
        
        alpha0 = alpha1 = alpha2 = alpha3 = alpha4 = alpha5 = alpha6= alpha7 = 1

        if(maxPMindex == 0):
            i = i - 1
            j = j - 1
            alpha7 = 0
        elif(maxPMindex == 1):
            j = j - 1
            alpha6 = 0
        elif(maxPMindex == 2):
            i = i + 1
            j = j - 1
            alpha5 = 0
        elif(maxPMindex == 3):
            i = i - 1
            alpha4 = 0
        elif(maxPMindex == 4):
            i = i + 1
            alpha3 = 1
        elif(maxPMindex == 5):
            i = i - 1
            j = j + 1
            alpha2 = 0
        elif(maxPMindex == 6):
            j = j + 1
            alpha1 = 0
        else:
            i = i + 1
            j = j + 1
            alpha0 = 0

    if(grid[i * cols + j][2].type == "SP"):
        grid[i * cols + j][2].TE = True

def simulation():
    for i in range(rows):
        for j in range(cols):
            if (grid[i * cols + j][2].type != "U"):
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
                
                if (maxCHA == valuesCHA[0]):
                    NW = PAP
                    SE = -PAP
                elif(maxCHA == valuesCHA[1]):
                    N = PAP
                    S = -PAP
                elif (maxCHA == valuesCHA[2]): 
                    NE = PAP
                    SW = -PAP
                elif (maxCHA == valuesCHA[3]):
                    W = PAP
                    E = -PAP
                elif (maxCHA == valuesCHA[4]):
                    E = PAP
                    W = -PAP
                elif (maxCHA == valuesCHA[5]):
                    SW = PAP
                    NE = -PAP
                elif (maxCHA == valuesCHA[6]):
                    S = PAP
                    N = -PAP
                elif (maxCHA == valuesCHA[7]):
                    SE = PAP
                    NW = -PAP

                pmVN = (((1 + W) * getPM(i - 1, j)) - (getAA(i - 1, j) * getPM(i, j))
                    + ((1 + N) * getPM(i, j - 1)) - (getAA(i, j - 1) * getPM(i, j))
                    + ((1 + E) * getPM(i + 1, j)) - (getAA(i + 1, j) * getPM(i, j))
                    + ((1 + S) * getPM(i, j + 1)) - (getAA(i, j + 1) * getPM(i, j))
                )
                pmMN = (((1 + NW) * getPM(i - 1, j - 1)) - (getAA(i - 1, j - 1) * getPM(i, j))
                    + ((1 + NE) * getPM(i + 1, j - 1))- (getAA(i + 1, j - 1) * getPM(i, j))
                    + ((1 + SW) * getPM(i - 1, j + 1)) - (getAA(i - 1, j + 1) * getPM(i, j))
                    + ((1 + SE) * getPM(i + 1, j + 1)) - (getAA(i + 1, j + 1) * getPM(i, j))
                )
                
                grid[i * cols + j][2].PM = getPM(i, j) + PMP1 * (pmVN + PMP2 * pmMN) if getPM(i, j) + PMP1 * (pmVN + PMP2 * pmMN) <= 1000 else 1000

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

                grid[i * cols + j][2].CHA = CON * ((getCHA(i, j) + CAP1 * (chaVN + CAP2 * chaMN)))
                
                if(grid[i * cols + j][2].CHA > 100):
                    grid[i * cols + j][2].CHA = 100
                elif(grid[i * cols + j][2].CHA < 0):
                    grid[i * cols + j][2].CHA = 0

            if(grid[i * cols + j][2].type == "NS" and grid[i * cols + j][2].PM >= thresholdPM):
                setTE(i, j)
                grid[i * cols + j][2].type = "SP"

# generating empty grid
for y in range(rows):
    for x in range(cols):
        rect = pygame.Rect(x * (squareSize + 1), y * (squareSize + 1), squareSize, squareSize)
        grid.append((rect, colorW, Cell()))

while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT: # event to exit the simulation
            pygame.quit()
            sys.exit()
        elif event.type == pygame.KEYDOWN: # press enter to start simulation with the configuration 
            if event.key == pygame.K_RETURN:
                if done:
                    done = False
                else:
                    done = True

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
                    cell.PM = 1000
                    grid[index] = (rect, colorY, cell)
    
    if done:
        simulation() # start simulation
        for index, (rect, color, cell) in enumerate(grid):
            if(cell.PM != 0 and (cell.type != "NS" and cell.type != "U")):
                alpha = clip((int)(cell.PM), 0, 255)
                color = (255, 255, 255 - alpha)
                grid[index] = (rect, color, cell)

            if(cell.TE == True):
                grid[index] = (rect, colorB, cell)
            
    # Now draw the rects. You can unpack the tuples
    # again directly in the head of the for loop.
    for rect, color, cell in grid:
        pygame.draw.rect(screen, color, rect)
    
    pygame.display.flip()
    clock.tick(10)