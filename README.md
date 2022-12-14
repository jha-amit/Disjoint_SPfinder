# Disjoint_SPfinder

This is a django based application for computing disjoint pair of least-cost paths separated by a minimum distance.
We realized that in hilly areas where the ground is higly undulating, the available least-cost path routing software
packages may give suboptimal soulution as they usually find two consecutive least-cost paths separated by the minimum
distance which may be more expensive than some combination of two paths. So, we minimize the combined cost of paths 
using a heuristics 'Optimal transmission network topology for resilient power supply' developed by \cite{Zinchenko Y. et al. 2016}.

Input:
![image](https://user-images.githubusercontent.com/57409254/203500167-5553da05-6b52-4c25-a4bd-2165833c4642.png)

The cost is per-unit construction cost in the vicinity of the point located by the latitude and longitudes.

![image](https://user-images.githubusercontent.com/57409254/203500751-d2af72c9-3eaa-471a-8841-82906b1e8786.png)


These costs are scattered data points, we interpolate costs at nodes first, using Inverse distance weight, and then use Bezier surface for cost map.

![image](https://user-images.githubusercontent.com/57409254/203500939-4a6f97b1-903a-48a0-b8a7-08d516c23d0f.png)

Then we compute edge cost:

A flow chart representing the backend of the application is the following:

![image](https://user-images.githubusercontent.com/57409254/203501288-9fe82768-53d1-4cf8-b2f6-66c33885415f.png)


The algorithm uses the following set up and processes:
1. The project area is represented as a diamond shape biconnected digraph and two such orthogonal graphs are combine as a 3D graph using affine transformations.
   This 3D embedding only allows eligible combination of nodes from 2D graph that provide embedded separation constraint and biconnectivity constraint in the 3D graph edges.
   ![2-Paths plan view](https://github.com/jha-amit/Disjoint_SPfinder/blob/main/3Dtransformation%20edited.JPG "3D embedding of two orthogonal 2D graphs")
   
2. A single shortest path is computed in the 3D which overcomes the computational intractibility of finding disjoint paths with constraints on a graph in polynomial time.
3. This single least-cost path in 3D is projectes back to 2D grpahs. Following images may help to clear out the process.
   ![3D to 2D](https://github.com/jha-amit/Disjoint_SPfinder/blob/main/3D%20embedding.JPG)
4. Near the terminals we introduced two complete graphs so that we can avoid the intractability due to path separation constraint in the neighbourhood of the start and end. We relaxed the path separation constraint in the neighbourhood of terminals and computed disjoint paths. Then we computed two least cost paths on the two complete graphs that connect with the disjoint paths computed on the diamond grid.
   
The further info on the parameters to be entered by the user are in the file ![Read me parameters description](https://github.com/jha-amit/Disjoint_SPfinder/blob/main/static/Read_me.txt)
![image](https://user-images.githubusercontent.com/57409254/203498454-8c28d74d-b0db-4de5-ae1a-1b77cd526cc1.png)

Here K is number of layers e.g. in the figure above the number of layers are 7. Number of nodes (N x N) in the 2D graph is proportional to number of layers N =(K-1)/2 + 1

The MIP formulation was taking more than 2hrs for only 200x200 grid!!!

The matlab's Mexfiles and python/cython implementation are comparable.

