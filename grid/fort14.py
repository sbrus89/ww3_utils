import numpy as np

class fort14:
  
  def __init__(self,filename='',verbose=False):

    if filename:    
      self.filename = filename
      self.f = open(filename,'r').read().splitlines()

    self.name = None     # name of grid (from file) 

    self.ne = None       # number of elements
    self.nn = None       # number of nodes

    self.xy = None       # x,y coordinates of nodes
    self.depth = None    # depth of nodes
    self.ect = None      # connectivity of nodes
 
    self.nope = None     # number of open boundaries
    self.neta = None     # total number of open boundary nodes
    self.obseg = None    # number of nodes in each open boundary
    self.obnds = None    # nodes in each open boundary

    self.nbou = None     # number of flow boundaries
    self.nvel = None     # total number of flow boundary nodes
    self.fbseg = None    # number of nodes in each flow boundary and boundary type  
    self.fbnds = None    # nodes in each flow boundary

    self.verbose = verbose

########################################################################
########################################################################

  def read_grid(self):
    self.read_header()
    self.read_nodes()
    self.read_connectivity()
    self.read_open_boundaries()
    self.read_flow_boundaries()

########################################################################
########################################################################

  def read_header(self):

    # Read grid name
    self.name = self.f.pop(0)
    line = self.f.pop(0).split()

    # Read number of elements and nodes
    self.ne = int(line[0])
    self.nn = int(line[1])

    if self.verbose:
      print(self.ne)
      print(self.nn)

########################################################################
########################################################################

  def read_nodes(self):

    # Initialize node coordinate and depth arrays
    self.xy = np.zeros((self.nn,2),dtype=np.double)
    self.depth = np.zeros((self.nn,),dtype=np.double)

    # Read coordinates and depths
    for i in range(self.nn):
      line = self.f.pop(0).split()
      j = int(line[0])-1
      self.xy[j,0] = float(line[1])
      self.xy[j,1] = float(line[2])
      self.depth[j] = float(line[3])

    if self.verbose:
      print(self.xy)
      print(self.depth)
      
########################################################################
########################################################################

  def read_connectivity(self):

    # Initialize element connectivity table
    self.ect = np.zeros((self.ne,3),dtype=np.int32)

    # Read element connectivity
    for i in range(self.ne):
      line = self.f.pop(0).split()
      j = int(line[0])-1
      for k in range(int(line[1])):
        self.ect[i,k] = int(line[k+2])

    if self.verbose:
      print(self.ect)

########################################################################
########################################################################

  def read_open_boundaries(self):

    # Read number of open boundaries
    line = self.f.pop(0).split()
    self.nope = int(line[0])

    # Read total number of open boundary nodes
    line = self.f.pop(0).split()
    self.neta = int(line[0])

    # Initialize open boundary arrays
    self.obseg = np.zeros((self.nope,),dtype=np.int32)
    self.obnds = np.zeros((self.nope,self.neta),dtype=np.int32)

    # Read open boundary nodes
    for i in range(self.nope):
      line = self.f.pop(0).split()
      self.obseg[i] = int(line[0])    # number of nodes in boundary
      for j in range(self.obseg):
        line = self.f.pop(0).split()
        self.obnds[i,j] = line[0]     # node numbers in boundary

    if self.verbose:
      print(self.nope)
      print(self.obseg)
      print(self.obnds)

########################################################################
########################################################################

  def read_flow_boundaries(self):

    # Read number of flow boundaries
    line = self.f.pop(0).split()
    self.nbou = int(line[0])

    # Read total number of flow boundary nodes
    line = self.f.pop(0).split()
    self.nvel = int(line[0])

    # Initialize flow boundary arrays
    self.fbseg = np.zeros((self.nbou,2),dtype=np.int32)
    self.fbnds = np.zeros((self.nbou,self.nvel),dtype=np.int32)

    # Read flow boundary nodes
    for i in range(self.nbou):
      line = self.f.pop(0).split()
      self.fbseg[i,0] = int(line[0])   # number of nodes in boundary
      self.fbseg[i,1] = int(line[1])   # boundary type
      for j in range(self.fbseg):
        line = self.f.pop(0).split()
        self.fbnds[i,j] = line[0]      # node numbers in boundary

    if self.verbose:
      print(self.nbou)
      print(self.fbseg)
      print(self.fbnds)

########################################################################
########################################################################

  def write_grid(self,filename):
  
    self.fout = open(filename,'w')
    self.write_header()
    self.write_nodes()
    self.write_connectivity()
    self.write_open_boundaries()
    self.write_flow_boundaries()

########################################################################
########################################################################

  def write_header(self):

    self.fout.write(self.name+'\n')
    self.fout.write('{:8d} {:8d}\n'.format(self.ne,self.nn))

########################################################################
########################################################################

  def write_nodes(self):

    for i in range(self.nn):
      self.fout.write('{:8d} {:24.17e} {:24.17e} {:24.17e}\n'.format(i+1,self.xy[i,0],self.xy[i,1],self.depth[i]))

########################################################################
########################################################################
  
  def write_connectivity(self):

    for i in range(self.ne):
      self.fout.write('{:8d} {:8d} {:8d} {:8d} {:8d}\n'.format(i+1,3,self.ect[i,0],self.ect[i,1],self.ect[i,2]))

########################################################################
########################################################################

  def write_open_boundaries(self):
    
    self.fout.write('{:8d}'.format(self.nope)+19*' '+'! number of open boundaries\n')
    self.fout.write('{:8d}'.format(self.neta)+19*' '+'! number of total open boundary nodes\n')
    for i in range(self.nope):
      self.fout.write('{:8d}'.format(self.obseg[i])+19*' '+'! number of nodes in open boundary {:8d}\n'.format(i+1))
      for j in range(selg.obseg[i]):
        self.fout.write('{:8d}\n'.format(self.obnds[i,j]))

########################################################################
########################################################################

  def write_flow_boundaries(self):
    
    self.fout.write('{:8d}'.format(self.nbou)+19*' '+'! number of flow boundaries\n')
    self.fout.write('{:8d}'.format(self.nvel)+19*' '+'! number of total flow boundary nodes\n')
    for i in range(self.nbou):
      self.fout.write('{:8d} {:8d}'.format(self.fbseg[i,0],self.fbseg[i,1])+19*' '+'! number of nodes in flow boundary {:8d}\n'.format(i+1))
      for j in range(selg.fbseg[i,0]):
        self.fout.write('{:8d}\n'.format(self.fbnds[i,j]))

########################################################################
########################################################################

  def write_gmsh(self,filename):
    
    self.fout = open(filename,'w')
    self.fout.write('$MeshFormat\n')
    self.fout.write('2 0 8\n')
    self.fout.write('$EndMeshFormat\n')
    self.fout.write('$Nodes\n')
    self.fout.write('{:d}\n'.format(self.nn))
    for i in range(self.nn):
      self.fout.write('{:d} {:16.8f} {:16.8f} {:16.8f}\n'.format(i+1,self.xy[i,0],self.xy[i,1],self.depth[i]))
    self.fout.write('$EndNodes\n')
    self.fout.write('$Elements\n')
    self.fout.write('{:d}\n'.format(self.ne))
    for i in range(self.ne):
      self.fout.write('{:d}  2  3  0  {:d}  0  {:d}  {:d}  {:d}\n'.format(i+1,i+1,self.ect[i,0],self.ect[i,1],self.ect[i,2]))
    self.fout.write('$EndElements')
        
