# Grupo:
#  Nome: Abner Eduardo Silveira Santos  NUSP: 10692012
#  Nome: Gyovana Mayara Moriyama        NUSP: 10734387
#  Nome: Henrique Matarazo Camillo      NUSP: 10294943
#  Nome: João Pedro Uchôa Cavalcante    NUSP: 10801169

from __future__ import print_function
import numpy as np
from ortools.linear_solver import pywraplp
import matplotlib.pyplot as plt

# A galaxy, with name and position in space
class Galaxy:
    count = 0

    def __init__(self, x, y):
        self.position = (x, y)
        
        self.name = str(Galaxy.count)

        Galaxy.count += 1

# A path between 2 galaxies
class Edge:
    def __init__(self, origin, destination):
        self.origin = origin
        self.destination = destination
        self.inpath = None
        self.distance = distance(origin.position, destination.position)
        self.name = origin.name[0] + destination.name[0].lower()
    
    def __str__(self):
        return self.name + ' = ' + str(self.distance)

# Distance between 2 galaxies
def distance(a, b):
    return np.sqrt(np.power(a[0]-b[0], 2) + np.power(a[1]-b[1], 2))

# Total path distance (objective function)
def route_distance(route):
    dist = 0
    
    for i in range(len(route)-1):
        dist += distance(route[i].position, route[i+1].position)
    dist += distance(route[len(route)-1].position, route[0].position)
    return dist

# ================== 2-OPT Path Improvement Methods ==================

def two_opt_swap(route, i, k):
    new_route = []
    
    for j in range(i):
        new_route.append(route[j])
    for j in range(k, i-1, -1):
        new_route.append(route[j])
    for j in range(k+1, len(route)):
        new_route.append(route[j])

    return new_route     

def two_opt(route):
    best_distance = route_distance(route)

    for i in range(1, len(route)):
        for k in range(i + 1, len(route) - 2):
            new_route = two_opt_swap(route, i, k)
            new_distance = route_distance(new_route)
            
            if new_distance < best_distance:
                best_distance = new_distance
                route = new_route
    
    return route

# ====================== Plot auxiliary methods ==========================

def plot(x_values, y_values, names, fig_size, title, optimal, objective):
    # Plots found path    
    plt.figure(figsize=fig_size)
    plt.title(title)
    plt.plot(x_values, y_values)
    plt.scatter(x_values, y_values, marker='.', c='r', s=130)
    for i in range(len(x_values)):
        plt.annotate(names[i], (x_values[i], y_values[i]))
    plt.gca().invert_xaxis()
    plt.show()

    print('\n')
    print(title + ':')
    print('Objective value = ', objective)
    print('GAP = ', (objective - optimal)/optimal)

def plot_by_galaxies(galaxies, edges, fig_size, title, optimal):
    # Variables used to plot and path recovering
    x_values = []        
    y_values = []        
    names = []
    route = []
    i = 0
    j = 0

    x_values.append(galaxies[0].position[0])
    y_values.append(galaxies[0].position[1])
    names.append(galaxies[0].name)
    route.append(galaxies[0])
    
    while j < len(edges):
        if (i != j) and (edges[i][j] == 1.0):
            # Saves point and galaxy name
            x_values.append(galaxies[j].position[0])
            y_values.append(galaxies[j].position[1])
            names.append(galaxies[j].name)
            route.append(galaxies[j])

            # Goes to next galaxy
            i = j
            j = 0
            if i == 0:
                break
        # Keeps searching for the found path
        else:
            j += 1
            
    plot(x_values, y_values, names, fig_size, title, optimal, route_distance(route))
    
    return route

def plot_by_edges(edges, fig_size, title, optimal, objective):
    # Variables used to plot and path recovering
    x_values = []        
    y_values = []        
    names = []
    route = []
    i = 0
    j = 0

    x_values.append(edges[0][0].origin.position[0])
    y_values.append(edges[0][0].origin.position[1])
    names.append(edges[0][0].origin.name)
    route.append(edges[0][0].origin)
    
    while j < len(edges):
        if (i != j) and (edges[i][j].inpath.solution_value() == 1.0):
            # Saves point and galaxy name
            x_values.append(edges[i][j].destination.position[0])
            y_values.append(edges[i][j].destination.position[1])
            names.append(edges[i][j].destination.name)
            route.append(edges[i][j].destination)

            # Goes to next galaxy
            i = j
            j = 0
            if i == 0:
                break
        # Keeps searching for the found path
        else:
            j += 1

    plot(x_values, y_values, names, fig_size, title, optimal, objective)
    
    return route

def plot_by_route(route, optimal, title, fig_size):
    x_values = []        
    y_values = []        
    names = []
    
    for i in range(len(route)):
        x_values.append(route[i].position[0])
        y_values.append(route[i].position[1])
        names.append(route[i].name)

    plot(x_values, y_values, names, fig_size, title, optimal, route_distance(route))
    
    return route

# =====================================================================

def main():
    # Creates galaxies from input
    galaxies = []

    n = int(input())
    best_solution = int(input())
    fig_size_buffer = input().split()
    fig_size = (float(fig_size_buffer[0]), float(fig_size_buffer[1]))

    for i in range(n):
        point = input().split()
        galaxy = Galaxy(float(point[1]), float(point[0]))
        galaxies.append(galaxy)

    # Creates edges between all galaxies
    edges = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
           edges[i][j] = Edge(galaxies[i], galaxies[j])

# ============================ Solver Setup ==============================

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

# ====================== Heuristic for Warm Start ========================

    to_be_visited = []
    hint_edges = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        to_be_visited.append(galaxies[i])

    nearest = 0
    current_galaxy = to_be_visited.pop(nearest)
    
    while len(to_be_visited) > 0:
        min_distance = float('inf')

        for i in range(len(to_be_visited)):
            cur_distance = distance(current_galaxy.position, to_be_visited[i].position)
            if cur_distance < min_distance:
                min_distance = cur_distance
                nearest = i
        
        hint_edges[int(current_galaxy.name)][int(to_be_visited[nearest].name)] = 1.0
        current_galaxy = to_be_visited.pop(nearest)

    hint_edges[int(current_galaxy.name)][0] = 1.0
    
    route = plot_by_galaxies(galaxies, hint_edges, fig_size, "Nearest Neighbor Heuristic", best_solution)

# ==================== Improve Heuristic with 2-OPT ======================

    # The last connection of Nearest Neighbor heuristic is usually the worst one, and the 2-opt algorithm 
    # can't apply any otimization to the last connection, so we reorder the route to be able to otimize it
    reordered_route = []
    for i in range((len(route)-1)//2, len(route)-1):
        reordered_route.append(route[i])
    for i in range(0, (len(route)-1)//2):
        reordered_route.append(route[i])
    reordered_route.append(reordered_route[0])
            
    route = two_opt(reordered_route)

    plot_by_route(route, best_solution, "Nearest Neighbor With 2-OPT", fig_size)

# ==================== Setup Solver with Warm Start ======================

    hint_edges = [[0.0 for _ in range(n)] for _ in range(n)]

    refs = []
    values = []

    for i in range(len(route)-1):
        hint_edges[int(route[i].name)][int(route[i+1].name)] = 1.0
    hint_edges[int(route[len(route)-1].name)][int(route[0].name)] = 1.0

    for i in range(n):
        for j in range(n):
            if i != j:
                refs.append(edges[i][j].inpath)
                values.append(hint_edges[i][j])

    solver.SetHint(refs, values)

# ==================== SCIP Solution with Warm Start ======================

    solver.set_time_limit(1800000) # 30 mins
    # solver.set_time_limit(10000) # 10s

    solver.Solve()
    
    print('Problem solved in %f milliseconds' % solver.wall_time())
    print('Problem solved in %d iterations' % solver.iterations())
    print('Problem solved in %d branch-and-bound nodes' % solver.nodes())

    route = plot_by_edges(edges, fig_size, "SCIP Solver Solution with Warm Start", best_solution, objective.Value())
    
# ================= Improved SCIP Solution with 2-OPT  ====================

    route.append(route[0])
    route = two_opt(route)
    plot_by_route(route, best_solution, "SCIP Solution With 2-OPT", fig_size)

# Amém
if __name__ == '__main__':
    main()