# A-large-covariance-matrix-estimator-under-intermediate-spikiness-regimes
This repository contains data and code relative to the manuscript "A large covariance matrix estimator under intermediate spikiness regimes" by Matteo Farnè and Angela Montanari (https://www.sciencedirect.com/science/article/pii/S0047259X19301216). 

The MATLAB dataset 'supervisory_data.m' contains the covariance matrix and the labels of a selection of Euro Area banking supervisory data. Punctual data are not available due to confidentiality reasons. The dataset contains the covariance matrix 'C', and the relative labels of supervisory indicators, 'Labgood'. The labels may be interpreted exploiting the detailed description at the link https://www.eba.europa.eu/documents/10180/359626/Annex+III+-+FINREP+templates+IFRS.xlsx/049e48a4-e7c2-44c6-89b1-4086447bced9

Two MATLAB functions, called "UNALCE.m" and "POET.m", are provided. The former performs the new procedure for covariance matrix estimation, called UNALCE (UNshrunk ALgebraic Covariance Estimator). The latter performs POET covariance matrix estimation procedure (Fan et al., 2013, http://onlinelibrary.wiley.com/doi/10.1111/rssb.12016/abstract). 
Both functions contain the detailed explanation of input and output arguments. See also https://data.mendeley.com/datasets/nh97vfvhkt/4.
