import re, json, sys, random, copy
import numpy as np
import scipy
import pandas as pd
from collections import defaultdict
from scipy import linalg as LA
'''
ECE1524 – Service Providers Network
Assignment 1
Hatif Sattar – 997063387

This script is to parse file in format below
and save each node as part of a Graph
<nodes>
    <node id="CHINng">
    <description>Chicago_IL</description>
    <location longitude="-87.616700" latitude="-41.833300"/>
    <interfaces>
        <interface id="TO=IPLSng"/>
        <interface id="FROM=NYCMng"/>
    </interfaces>
    </node>
'''
file = open(sys.argv[1], "r+")
out_f = open("topo_graph.txt", "w+")
#file = open ("topo.xml", "r+")

#G = {}
G = defaultdict(list)
#Test Graphs
G2 = {"a" : ["b", "e", "g"],
    "b" : ["a", "e", "h", "c"],
    "c" : ["i", "b", "d", "f"],
    "d" : ["c", "f", "j"],
    "e" : ["a", "b", "g", "h"],
    "f" : ["c", "d", "i", "j"],
    "g" : ["a", "e", "h"],
    "h" : ["g", "e", "b", "i"],
    "i" : ["c", "f", "j", "h"],
    "j" : ["i", "f", "d"]
}

G3 = {"a" : ["b", "c"],
    "b" : ["a", "d", "c"],
    "c" : ["a", "b", "d"],
    "d" : ["c", "b"]
}

#Parse File into Dictionary
for line in file:
    if "node id" in line:
        sp = line.split("\"")
        city = sp[1]
        out_f.write("Node:" + city + "\n")
        myline = "  Nebors:"
        #print ("CITY ", city)
        for line in file:
            if "node>" in line:
                break
            if "FROM" in line and "interface id" in line:
                sp = line.split("\"")
                sp2 = sp[1].split("=")
                #print ("    " , sp2)
                G.setdefault(city, []).append(sp2[1])
                myline += sp2[1]
                myline += ","
        out_f.write(myline + "\n\n")

out_f.close()

#Select the graph to analyze ------------
#G = G3

#print Graph
for key in G:
    print (key)
    print (G[key])

#########################
# Find Node Connectivity
#########################

print ("Number of neboring nodes:\n")

min = 1000000
for key, value in G.items():
    #print value
    num_nebor = len(list(filter(bool, value)))
    print(key, num_nebor)
    if ( num_nebor < min):
        min = num_nebor

print (" ==> NODE CONNECTIVITY =", min)

###########################
# Find Edge Connectivity
###########################

def has_empty_node (g):
    for key, value in g.items():
        num = len(list(filter(bool, value)))
        #print ("Num nebor:",num)
        if (num == 0):
            return True
    return False

def choose_random_key(myG):
    v = random.choice(list(myG.keys()))
    v2 = random.choice(list(myG[v]))
    return v, v2

def karger(myG):
    length = []
    while len(myG) > 2:
        v, v2 = choose_random_key(myG)
        myG[v].extend(myG[v2])
        for x in myG[v2]:
            #print ("remove key:", x, " val:", v2)
            myG[x].remove(v2)
            myG[x].append(v)
        while v in myG[v]:
            myG[v].remove(v)
        del myG[v2]
    for key in myG.keys():
        length.append(len(myG[key]))
    return length[0]

def iterate(n):
    i = 0
    count = 10000   
    while i < n:
        data = copy.deepcopy(G)
        min_cut = karger(data)
        if min_cut < count:
            count = min_cut
        i = i + 1
    return count

print(" ==> EDGE CONNECTIVITY = " + str(iterate(100)), "\n")


##############################
# Find Algebraic Connectivity
##############################

def print_matrix (name, matrix, s):
    print (name, "(size:",s,"): \n")
    print ('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in matrix]))
    print ("\n")

sorted_keys=sorted(G.keys())
size=len(sorted_keys)

## Adjacency Matrix -------------------------------------
adjacent = [ [0]*size for i in range(size) ]
for a,b in [(sorted_keys.index(a), sorted_keys.index(b)) for a, row in G.items() for b in row]:
    if (a == b):
        adjacent[a][b] = 2
    else:
        adjacent[a][b] = 1
print_matrix ("Adjacency Matrix", adjacent, size)

## Diagonal Matrix --------------------------------------
diag = [ [0]*size for i in range(size) ]
for x in range(0,size):
    #print (G.get(sorted_keys[x], "NULL"))
    diag[x][x] = len(G.get(sorted_keys[x], "NULL"))
print_matrix("Diagonal matrix", diag, size)

## Laplacian Matrix ------------------------------------
lap = [ [0]*size for i in range(size) ]
for x in range(0,size):
    for y in range(0,size):
        lap[x][y] = diag[x][y] - adjacent[x][y]

print_matrix ("Laplacian Matrix", lap, size)

##test_matrix = np.matrix( [[1,2,3],[11,12,13],[21,22,23]])
##print (test_matrix)
##print ("\n")

eigen_vals, eigen_vects = LA.eig(lap)

e2 = np.sort(eigen_vals)
print ("\nEigen Values sorted:\n", e2)


## 2D array to store links removed ----------------------------
checked_array = [[0]*2 for i in range(size)]
#print(checked_arry)

##count = sum(len(v) for v in G.values())
#count = (count/2) if ((count%2)==0) else ((count/2) + 1)
##print("Count:", count, "\n")
##assert(count)


# This function to find Node/Link/Algebraic Conn. of G after
# removal of each edge.
def process_graph(key_, value_, theG):
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\nThe Edge - ", key_, "-", value_)
    # Node Connectivity--------
    min = 1000000
    for key, value in theG.items():
        #print value
        num_nebor = len(list(filter(bool, value)))
        #print(key, num_nebor)
        if ( num_nebor < min):
            min = num_nebor
    print (" ==> NODE CONNECTIVITY =", min, "\n")
    # Edge Connectivity---------
    count = 10000
    if (has_empty_node(theG) == True):
        count = 0
    else:
        i,n = 0, 100  
        while i < n:
            data = copy.deepcopy(theG)
            min_cut = karger(data)
            if min_cut < count:
                count = min_cut
            i = i + 1
    print(" ==> EDGE CONNECTIVITY = ", count, "\n")
    sorted_keys=sorted(theG.keys())
    size=len(sorted_keys)
    ## Adjacency Matrix -----------
    adjacent = [ [0]*size for i in range(size) ]
    for a,b in [(sorted_keys.index(a), sorted_keys.index(b)) for a, row in theG.items() for b in row]:
        if (a == b):
            adjacent[a][b] = 2
        else:
            adjacent[a][b] = 1
    #print_matrix ("Adjacency Matrix", adjacent, size) ###----- debug
    ## Diagonal Matrix ---------------
    diag = [ [0]*size for i in range(size) ]
    for x in range(0,size):
        diag[x][x] = len(theG.get(sorted_keys[x], "NULL"))
    #print_matrix("Diagonal matrix", diag, size)   ###-------- debug
    ## Laplacian Matrix ---------------
    lap = [ [0]*size for i in range(size) ]
    for x in range(0,size):
        for y in range(0,size):
            lap[x][y] = diag[x][y] - adjacent[x][y]
    #print_matrix ("Laplacian Matrix", lap, size)  #### -------- debug
    eigen_vals, eigen_vects = LA.eig(lap)
    e2 = np.sort(eigen_vals)
    print ("\nEigen Values sorted:\n", e2)



def mark_array (k, v, array):
    for i in range(len(array)):
        if (array[i][0] == 0 and array[i][1] == 0):
            array[i][0] = key
            array[i][1] = val
                
        

#Remove each edge and find Node/Link/Algebraic Conn. of Graph
Goo = copy.deepcopy(G)
for key, vals in Goo.items():
    #print ("For loop 1 - key = ", key)
    for val in vals:
        #print ("For loop 2 - Val = ", val)
        myG = copy.deepcopy(G)
        #find key/val pair in checked array
        found = False
        for i in range(len(checked_array)):
            if (checked_array[i][0] == key):
                if (checked_array[i][1] == val):
                    found = True
            if (checked_array[i][1] == key):
                if (checked_array[i][0] == val):
                    found = True
        if found == True:
            print ("FOUND Val !!\n")
            continue
        else:
            #Remove edge and save it in array
            #To remove, remove 2 vals from 2 keys
            print ("Key:", key, "\n")
            print ("Val:", val, " in vals: ",vals,"\n")
            myG.pop(key, None)
            vals_new = list(vals)
            vals_new.remove(val)
            ####print ("Now new vals: ",vals_new,"\n")
            ####print ("Check old vals: ",vals,"\n")
            myG[key] = vals_new
            #Remove the edge from other key
            vals_other = myG.get(val)
            print ("other key: ", val, " -> vals: ",vals_other,"\n")
            myG.pop(val, None)
            vals_other.remove(key)
            ####print ("other vals now: ",vals_other,"\n")
            myG[val] = vals_other
            #print(json.dumps(myG, indent = 4))
            process_graph(key, val, myG)
            mark_array(key, val, checked_array)


            
    


input ("press ENTER")
