### Introduction

This code provides tools for infering personality traits from behavior. While personality traits offer considerable insight into the biological basis of individual differences, existing approaches toward understanding personality across species rely on subjective criteria and limited sets of behavioral readouts, resulting in noisy and often inconsistent outcomes. In a recent paper, we introduced a mathematical framework for studying individual differences along dimensions with maximum consistency and discriminative power. This framework was validated in mice, using data from a system for high-throughput longitudinal monitoring of group-housed male mice that yields a variety of readouts from all across an individual’s behavioral repertoire. The code is implemented in Matlab.

### Getting Started

The file `IdentityDemos` is a good starting point as it includes several examples for using the code. The script shows how to:
1. Compute the identity-domains and showing their behavioral correlates
2. Use pre-computed Identity-Domains to produce a map of the personality space

### Prerequisites

The code was tested on Matlab 2017a on both a Windows and MacOS machines, but should in principal work on earlier version of Matlab as well as other operating systems

### Installation

Just put all the code in a directory and cd to this directory within Matlab

### Files

File | Description 
-----|------------
`IdentityDemos.m` |	Includes a set of examples for using the code
`IdentityDomains.m` |	Main class for computing the IDs
`DimReduction.m`	|	A collection of Dimensionality Reduction algorithms
`Auxiliary.m`	|	Auxiliary functions
`Color.m`		|	A collection of colors and tools for colors
`data_table.mat`	|	Sample table
`ID.mat`		|	Precomputed IDs and Archetypes

### Cite

Identity domains in complex behavior: Toward a biology of personality
Oren Forkosh, Stoyo Karamihalev, Simone Roeh, Mareen Engel, Uri Alon, Sergey Anpilov, Markus Nussbaumer, Cornelia Flachskamm, Paul Kaplick, Yair Shemesh, Alon Chen
bioRxiv 395111; doi: https://doi.org/10.1101/395111

### License

The Identity Domains Toolbox for MATLAB is free and open source, distributed under the [MIT license](https://opensource.org/licenses/MIT).

