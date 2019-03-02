C Code for Simple Estimators for Network Sampling, by Steven
K. Thompson, 2018.

The file fast10.c runs the fast sampling process on a network data set
to calculate the inclusion frequencies f_i for the nodes i = 1, 2,
..., n in the sample of size n.  These values are in the node
variable aadlagav.  The nodes are stored in a simple linked list.  The
links make this an adjacency linked list.  At iteration t the node
inclusion indicator is aad[0].  aad[1] stores the value for iteration
t-1.  These are averaged over iterations to give aadlagav, or f_i in
the paper's terminology.  
