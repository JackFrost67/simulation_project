# E = A u U, Entire Area
# A = Available Area, U = Unavailable Area
# S = Place where physarium is introduced
# N = Nutrient Source

# State a t time:
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
from Cell import Cell

size = (width, height) = 500, 500

pygame.init()

screen = pygame.display.set_mode(size) # setting screen size

pygame.display.set_caption('Physarum Polycephalum Simulation') # setting name of the screen

clock = pygame.time.Clock() # setting clock

squareSize = 10
cols, rows = int(screen.get_width() / squareSize), int(screen.get_height() / squareSize)

colorW = (255, 255, 255)
colorY = (255, 255, 0)
colorR = (255, 0, 0)
colorG = (0, 255, 0)

# generating empty grid
grid = []
for y in range(rows):
    for x in range(cols):
        rect = pygame.Rect(x * (squareSize + 1), y * (squareSize + 1), squareSize, squareSize)
        grid.append((rect, colorW, Cell()))

while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()

    # draw source, food and other stuff
    if pygame.mouse.get_pressed()[0]:
        mousePos = pygame.mouse.get_pos()
        for index, (rect, color, cell) in enumerate(grid):
            if rect.collidepoint(mousePos):
                # Create a tuple with the new color and assign it.
                if color == colorY:
                    cell.type = "U" # not available
                    grid[index] = (rect, colorR, cell)
                elif color == colorR:
                    cell.type = "NS" # food
                    grid[index] = (rect, colorG, cell)
                elif color == colorG:
                    cell.type = "A" # avalaible
                    grid[index] = (rect, colorW, cell)
                else:
                    cell.type = "SP" # starting point
                    grid[index] = (rect, colorY, cell)

    # Now draw the rects. You can unpack the tuples
    # again directly in the head of the for loop.
    for rect, color, cell in grid:
        pygame.draw.rect(screen, color, rect)

    pygame.display.flip()
    clock.tick(10)