# Grupo:
#  Nome: Abner Eduardo Silveira Santos  NUSP: 10692012
#  Nome: Gyovana Mayara Moriyama        NUSP: 10734387
#  Nome: Henrique Matarazo Camillo      NUSP: 10294943
#  Nome: João Pedro Uchôa Cavalcante    NUSP: 10801169

from __future__ import print_function
import numpy as np
from ortools.linear_solver import pywraplp

namedict = {
    'A': "Andromeda",
    'B': "Black Eye",
    'C': "Cartwheel",
    'D': "Drip Drop",
    'E': "Eye of Sauron",
    'F': "Fireworks",
    'G': "Gusty Garden",
    'H': "Honeyhive",
    'I': "Imelino",
    'J': "Jipperino",
    'K': "Kurzgesagt",
    'L': "Little Sombrero",
    'M': "Medusa Merger",
    'N': "Neigel"
}

class Galaxy:
    count = 0
    def __init__(self, x, y):
        self.position = (x, y)
        self.name = namedict[chr(ord('A') + Galaxy.count)]
        Galaxy.count += 1

class Edge:
    def __init__(self, origin, destination):
        self.origin = origin
        self.destination = destination
        self.inpath = None
        self.distance = distance(origin.position, destination.position)
        self.name = origin.name[0] + destination.name[0].lower()
    
    def __str__(self):
        return self.name + ' = ' + str(self.distance)
    
def distance(a, b):
    return np.sqrt(np.power(a[0]-b[0], 2) + np.power(a[1]-b[1], 2))

def main():
    galaxies = []

    n = int(input())
    for i in range(n):
        point = input().split()
        galaxy = Galaxy(int(point[0]), int(point[1]))
        galaxies.append(galaxy)

    edges = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
           edges[i][j] = Edge(galaxies[i], galaxies[j])

    for i in range(n):
        print(galaxies[i].name)
        for j in range(n):
            print(edges[i][j], end=' ')
        print()

    # Create the linear solver with the GLOP backend.
    solver = pywraplp.Solver.CreateSolver('GLOP')

    # Create the variables.
    for i in range(len(edges)):
        for j in range(len(edges[0])):
            edges[i][j].inpath = solver.BoolVar(edges[i][j].name)
    
    print('Number of variables =', solver.NumVariables())

    # Create a linear constraint.
    for i in range(len(edges)):
        ct = solver.Constraint(1,1)
        for j in range(len(edges[0])):
            if i != j:
                ct.SetCoefficient(edges[i][j].inpath, 1)

    # Create a linear constraint.
    for i in range(len(edges)):
        ct = solver.Constraint(1,1)
        for j in range(len(edges[0])):
            if i != j: 
                ct.SetCoefficient(edges[j][i].inpath, 1)

    print('Number of constraints =', solver.NumConstraints())

    # Create the objective function.
    objective = solver.Objective()
    for i in range(len(edges)):
        for j in range(len(edges[0])):
            if i != j:
                objective.SetCoefficient(edges[j][i].inpath, edges[i][j].distance)
    objective.SetMinimization()

    solver.Solve()
    
    print('Solution:')
    print('Objective value =', objective.Value())
        
    for i in range(len(edges)):
        for j in range(len(edges[0])):
            if (i != j) and (edges[i][j].inpath.solution_value() == 1.0):
                print(edges[i][j].destination.name + ', ', end='')
                break

if __name__ == '__main__':
    main()