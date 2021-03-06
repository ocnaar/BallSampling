# BallSampling
The function "nBallSampling.m" outputs m random samples within an n-dimensional L-p ball of radius r centered at c.

Samples within L-2 balls are obtained by mapping values from Latin Hypercube Sampling (LHS) into spherical coordinates of n-dimensional hyperspheres. 
Samples of L-1 balls come from applying a rejection method on the points sampled from an L-2 ball of the same radius. 
Samples of L-Inf balls are given by isotropic scaling of vectors from LHS.
