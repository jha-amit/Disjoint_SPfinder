1. Diamond graph edge slope: This input is side slopes of the grid with line connecting start to end.

2. Truncation layer gap from outer most circle controls the radius of the radial grid near terminals.
 The tighter the gap the larger the radius.

3. Horizantal spacing : It represents the path spacing measured perpendicular to the line joining start
 and end of diamond grid. Enter this in multiples of 2, e.g. if you want a physical path separation of 51km
 enter d=50.

4. Buffer layer: This is extra layers we want to remove for the diajoint paths to choose from a combination
 of two paths separated by a diatace d. If we provide Buffer layer = 0 the diajoint paths will be pushed to 
 the boundary of the grid in the beginning. We can pick any non negative integer.

5,6,7,8. Latitude, Longitude informaton of start and end points of the diamond grid.

9. Radial increment start/end: This input is for number of rings in the radial graph which controls the density of
 the radial graphs at start and end.

10. API key: Enter your Google Maps API key. Get it from https://developers.google.com/maps/documentation/maps-static

11. Scaling factor for the grid: It scales the 2D grid to the desired edge length.

12. Grid density : This is density of Bezier rectangular patch. E.g, if we want a n x n number of interpolated points
we enter n. This number is also used in integrating the costs along graph edges using simpson's formula. 
Enter a number of format 3k+1 (k= positive integer)

13. Latitude difference for selecting inputs for interpolation: As name signifies, this value determines the box within
within which the known points should be selected for interpolation using inverse distance wight method. e.g. limit of 0.1
will select points within absolute difference of 0.1 in in latitide.

14. Longitude difference for selecting inputs for interpolation: Same as latitude..

15. High value: some high value for replacing infinite costs.

16. Input cost file for terminal graph at start : Upload a CSV file with latitude, longitude, and corresponding costs
in three columns with headers 'Latitude', 'Longitude', and 'Cost'. Make sure that the Latitude Longitudes are in the
viscinity of the radial grid lat-long with suficient density such that there are a few points withon the set limits (within inputs 10 and 11).

17. Input cost file for diamond grid: Input file for diamond grid in the same format as input.

18, 19, 20 are same as the above inputs 16, 14.

21 Input node number to be modified : This is the node ID selected for modifying the edge costs start from the node.
The node ID can be selected from the marker tool provided as a toggle button in the user interface.

22, 23 are modifies costs of edges starting from the selected node. Here the cost Horizantal is the edge which is 
anticlockwise from the line joing start and end. Cost vertical is clockwise from the line joining the start and end.





