# Grupo:
#  Nome: Abner Eduardo Silveira Santos  NUSP: 10692012
#  Nome: Gyovana Mayara Moriyama        NUSP: 10734387
#  Nome: Henrique Matarazo Camillo      NUSP: 10294943
#  Nome: João Pedro Uchôa Cavalcante    NUSP: 10801169

from __future__ import print_function
import numpy as np
from ortools.linear_solver import pywraplp
import matplotlib.pyplot as plt

# namedict = {
#     'A': "Andromeda",
#     'B': "Black Eye",
#     'C': "Cartwheel",
#     'D': "Drip Drop",
#     'E': "Eye of Sauron",
#     'F': "Fireworks",
#     'G': "Gusty Garden",
#     'H': "Honeyhive",
#     'I': "Imelino",
#     'J': "Jipperino",
#     'K': "Kurzgesagt",
#     'L': "Little Sombrero",
#     'M': "Medusa Merger",
#     'N': "Neigel"
# }

class Galaxy:
    count = 0
    # max_count = len(namedict)

    def __init__(self, x, y):
        self.position = (x, y)
        
        # if Galaxy.count < Galaxy.max_count:
        #     self.name = namedict[chr(ord('A') + Galaxy.count)]
        # else:
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
        galaxy = Galaxy(float(point[0]), float(point[1]))
        galaxies.append(galaxy)

    # Creates edges between all galaxies
    edges = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
           edges[i][j] = Edge(galaxies[i], galaxies[j])

    # Create the linear solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')
    AT_LOWER_BOUND = solver.AT_LOWER_BOUND

    u = []

    # Create the variables.
    for i in range(len(edges)):
        u_i = solver.IntVar(1, n-1, 'u' + str(i))
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

                # ui - uj * n * Xij <= n - 1
                # Prevents subcycles
                if i >= 1 and j >= 1:
                    ct3 = solver.Constraint(-solver.infinity(), n-1) # <=
                    ct3.SetCoefficient(u[i], 1) # ui
                    ct3.SetCoefficient(u[j], -1) # - uj
                    ct3.SetCoefficient(edges[i][j].inpath, n) # + n*X

                # Cost is the sum of the distance of all chosen paths
                objective.SetCoefficient(edges[i][j].inpath, edges[i][j].distance)

    

    objective.SetMinimization()
    objective.BestBound()

    print('Number of constraints =', solver.NumConstraints())

    to_be_visited = []
    hint_edges = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        to_be_visited.append(galaxies[i])

    nearest = 0
    current_galaxy = to_be_visited.pop(nearest)
    
    while len(to_be_visited) > 0:
        # print("Len: " + str(len(to_be_visited)))
        min_distance = float('inf')

        for i in range(len(to_be_visited)):
            cur_distance = distance(current_galaxy.position, to_be_visited[i].position)
            if cur_distance < min_distance:
                min_distance = cur_distance
                nearest = i
                
            
        # print(int(current_galaxy.name))
        # print(int(to_be_visited[nearest].name))
        hint_edges[int(current_galaxy.name)][int(to_be_visited[nearest].name)] = 1.0
        current_galaxy = to_be_visited.pop(nearest)
    hint_edges[int(current_galaxy.name)][0] = 1.0
    
    # Variables used to plot and path recovering
    x_values = []        
    y_values = []        
    names = []
    i = 0
    j = 0
    hint_u = [None for _ in range(n)]

    x_values.append(galaxies[0].position[0])
    y_values.append(galaxies[0].position[1])
    names.append(galaxies[0].name)

    hint_u[0] = 0
    counter = 1
    # Recovers path found by the solver
    while j < len(hint_edges):
        if (i != j) and (hint_edges[i][j] == 1.0):
            # print(str(edges[i][j].origin.name[0]) +  str(edges[i][j].destination.name[0]) + ', ', end='')
            # print(edges[i][j].destination.name + ', ', end='')
            # print(i, '-> ', end='')

            # Saves point and galaxy name
            x_values.append(galaxies[j].position[0])
            y_values.append(galaxies[j].position[1])
            names.append(galaxies[j].name)

            # Goes to next galaxy
            i = j
            counter += 1
            hint_u[i] = counter
            j = 0
            if i == 0:
                break
        # Keeps searching for the found path
        else:
            j += 1

    # print(names[0][0])

    # Plots found path    
    # plt.plot(x_values, y_values)
    # plt.scatter(x_values, y_values, marker='*', c='r', s=130)
    # for i in range(len(x_values)):
    #     plt.annotate(names[i], (x_values[i], y_values[i]))
    # plt.show()

    refs = []
    values = []
    for i in range(n):
        refs.append(u[i])
        values.append(hint_u[i])
        for j in range(n):
            if i != j:
                refs.append(edges[i][j].inpath)
                values.append(hint_edges[i][j])

    solver.SetHint(refs, values)
    # help(solver)
    # help(solver.separating)
    # solver.SetSeparating("SCIP_PARAMSETTING_AGRESSIVE")
    # solver.SCIPsetSeparating(SCIP, "SCIP_PARAMSETTING_AGRESSIVE", quiet)

    solver.set_time_limit(600000)

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
            print(edges[i][j].origin.name, '-> ', end='')

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

    print('Problem solved in %f milliseconds' % solver.wall_time())
    print('Problem solved in %d iterations' % solver.iterations())
    print('Problem solved in %d branch-and-bound nodes' % solver.nodes())

    # Plots found path    
    plt.plot(x_values, y_values)
    plt.scatter(x_values, y_values, marker='.', c='r', s=130)
    for i in range(len(x_values)):
        plt.annotate(names[i], (x_values[i], y_values[i]))
    plt.show()


if __name__ == '__main__':
    main()