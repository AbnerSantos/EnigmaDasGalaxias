# Grupo:
#  Nome: Abner Eduardo Silveira Santos  NUSP: 10692012
#  Nome: Gyovana Mayara Moriyama        NUSP: 10734387
#  Nome: Henrique Matarazo Camillo      NUSP: 10294943
#  Nome: João Pedro Uchôa Cavalcante    NUSP: 10801169

from __future__ import print_function
import numpy as np
from ortools.linear_solver import pywraplp
import matplotlib.pyplot as plt

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
    max_count = len(namedict)

    def __init__(self, x, y):
        self.position = (x, y)
        
        if Galaxy.count < Galaxy.max_count:
            self.name = namedict[chr(ord('A') + Galaxy.count)]
        else:
            self.name = str(Galaxy.count)

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
    # Creates galaxies from input
    galaxies = []

    n = int(input())
    for i in range(n):
        point = input().split()
        galaxy = Galaxy(int(point[0]), int(point[1]))
        galaxies.append(galaxy)

    # Creates edges between all galaxies
    edges = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
           edges[i][j] = Edge(galaxies[i], galaxies[j])

    # Create the linear solver with the GLOP backend.
    solver = pywraplp.Solver.CreateSolver('GLOP')

    u = []

    # Create the variables.
    for i in range(len(edges)):
        u_i = solver.NumVar(1, n-1, 'u' + str(i))
        u.append(u_i)
        for j in range(len(edges[0])):
            if i != j:
                edges[i][j].inpath = solver.BoolVar(edges[i][j].name)
    
    print('Number of variables =', solver.NumVariables())

    # Create constraints
    objective = solver.Objective()
    for i in range(len(edges)):
        ct1 = solver.Constraint(1,1)
        ct2 = solver.Constraint(1,1)
        for j in range(len(edges[0])):
            if i != j:
                
                # From each galaxy, we can only have one exiting path
                ct1.SetCoefficient(edges[i][j].inpath, 1)

               
                # From each galaxy, we can only have one incoming path
                ct2.SetCoefficient(edges[j][i].inpath, 1)

                # # 0 <= Xy + Yx <= 1
                # # Prevents going back 
                # if i > j:
                #     ct3 = solver.Constraint(0,1)
                #     ct3.SetCoefficient(edges[j][i].inpath, 1)
                #     ct3.SetCoefficient(edges[i][j].inpath, 1)

                if(i >= 1 and i < n):
                    ct3 = solver.Constraint(-solver.infinity(), n-1) # <=
                    ct3.SetCoefficient(u[i], 1) # ui
                    
                    ct3.SetCoefficient(u[j], -1) # - uj
                    ct3.SetCoefficient(edges[i][j].inpath, n) # + nCoisa

                # Cost is the sum of the distance of all chosen paths
                objective.SetCoefficient(edges[i][j].inpath, edges[i][j].distance)

    objective.SetMinimization()

    print('Number of constraints =', solver.NumConstraints())

    solver.Solve()
    
    print('Solution:')
    print('Objective value =', objective.Value())

    # Variables used to plot and path recovering
    x_values = []        
    y_values = []        
    names = []
    i = 0
    j = 0

    x_values.append(edges[0][0].origin.position[0])
    y_values.append(edges[0][0].origin.position[1])
    names.append(edges[0][0].origin.name)

    # Recovers path found by the solver
    while j < len(edges):
        if (i != j) and (edges[i][j].inpath.solution_value() == 1.0):
            # print(str(edges[i][j].origin.name[0]) +  str(edges[i][j].destination.name[0]) + ', ', end='')
            # print(edges[i][j].destination.name + ', ', end='')
            print(edges[i][j].origin.name[0], '-> ', end='')

            # Saves point and galaxy name
            x_values.append(edges[i][j].destination.position[0])
            y_values.append(edges[i][j].destination.position[1])
            names.append(edges[i][j].destination.name)

            # Goes to next galaxy
            i = j
            j = 0
            if i == 0:
                break
        # Keeps searching for the found path
        else:
            j += 1

    print(names[0][0])

    # Plots found path    
    plt.plot(x_values, y_values)
    plt.scatter(x_values, y_values, marker='*', c='r', s=130)
    for i in range(len(x_values)):
        plt.annotate(names[i], (x_values[i], y_values[i]))
    plt.show()

if __name__ == '__main__':
    main()