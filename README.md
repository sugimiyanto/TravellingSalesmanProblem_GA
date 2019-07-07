# TravellingSalesmanProblem (Genetic Algorithm)
- This project is working on Artificial Intelligence area. The main idea is trying to solve Travelling Salesman problem by utilising Genetic Algorithm (GA) using Matlab.
- The added values of GA in this project are the using of Tournament method for parents selection, as well as the involvement of various crossovers.

## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

## Prerequisites
- Matlab version 2015 or later.

## Main program:
1. **TSP_GeneticAlgorithm.m** : uses sub-tour exchange crossover
2. **TSP_GeneticAlgorithm_withEdgeRecombine.m** : uses edge recombination crossover.

## Running the test
How to run:
1. Choose one of the above main programs.
2. If you wish, customize the tournament size, crossover probability ratio, mutation ratio, population size, and iteration size by modifying the variables such as sizeTournament, crosProb, mutRatio, sizePop, and iterSize respectively.
3. Run the program.
4. Select a dataset from folder `dataset`.
5. While running, it shows the progress graph in real-time which displays the optimum shortest path.
6. Done.

## Author
**Sugimiyanto** - *Initial Work* - [Sugimiyanto](https://github.com/sugimiyanto)

## Contributing
This project is far from perfection. I very welcome for any suggestion for the code improvement and knowledge sharing.
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests as appropriate.

## Acknowledgements
- Emad Alamoodi

## Annotation
 This project is generated as one of the tasks of Artificial Intelligence course (CS661) in the Department of Computer Science, Faculty of Computing and Information Technology (FCIT), King Abdulaziz University, Saudi Arabia.