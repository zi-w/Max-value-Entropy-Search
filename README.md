# This repository is under construction.

# Max-value-Entropy-Search
This is the code repository associated with the paper [_Max-value Entropy Search for Efficient Bayesian Optimization_](https://arxiv.org/abs/1703.01968). Please refer to the paper if you need more details on the algorithm.
## How to use the code for parameter tuning
example.m is a simple example using Bayesian optimization to maximize a function. In addition to the proposed Max-value Entropy Search method, we provide several other popular Bayesian optimization methods: Gaussian process upper confidence bound (GP-UCB) by Auer, 2002; Srinivas et al., 2010, probability of improvement (PI) by Kushner, 1964, expected improvement (EI) by Mockus, 1974, and optimization as estimation (EST) by Wang et al., 2016. In the paper, we also compared to [entropy search](http://www.probabilistic-optimization.org/Global.html) (Hennig & Schuler, 2012) and [predictive entropy search](https://bitbucket.org/jmh233/codepesnips2014) (Hernández-Lobato et al., 2014), whose code is open sourced and hence not included in this repository.

Please see the comments in the code for more details about the usage.

## Example functionals for optimization
In test-functions/, we provide some functionals one can use to test Bayesian optimization algorithms. 

0. Approximated functions sampled from Gaussian processes
1. Optimization test functions
2. Tuning hyper-parameters for neural networks 
3. Active learning for robot pushing
4. Tuning the walking speed of a planar bipedal robot

## System requirement
We tested our code with MATLAB R2015b on Ubuntu 14.04 LTS (64-bit) and Mac OS X (64-bit). We developed our code building upon the MATLAB code for basic Gaussian process functionals from predictive entropy search (Hernandez-Lobato et al., 2014), which is developed upon GPstuff (Vanhatalo et al., 2013). To run the example code in example.m, first make sure you installed the GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)). On ubuntu, you can install GSL by 

```
sudo apt-get install libgsl0-dev
```

To run Box2D related code in test-functions/robot-pushing/, please install [Pybox2d](https://github.com/pybox2d/pybox2d).

## References
* Wang, Zi and Jegelka, Stefanie. Max-value Entropy Search for Efficient Bayesian Optimization. arXiv preprint arXiv:1703.01968, 2017.
* Auer, Peter. Using confidence bounds for exploitationexploration tradeoffs. Journal of Machine Learning Research, 3:397–422, 2002.
* Srinivas, Niranjan, Krause, Andreas, Kakade, Sham M, and Seeger, Matthias. Gaussian process optimization in the bandit setting: No regret and experimental design. In International Conference on Machine Learning (ICML), 2010.
* Kushner, Harold J. A new method of locating the maximum point of an arbitrary multipeak curve in the presence of noise. Journal of Fluids Engineering, 86(1):97–106, 1964.
* Mockus, J. On Bayesian methods for seeking the extremum. In Optimization Techniques IFIP Technical Conference, 1974.
* Wang, Zi, Zhou, Bolei, and Jegelka, Stefanie. Optimization as estimation with Gaussian processes in bandit settings. In International Conference on Artificial Intelligence and Statistics (AISTATS), 2016.
* Catto, Erin. Box2d, a 2D physics engine for games. http://box2d.org, 2011.
* Pybox2d, 2D Game Physics for Python. http://github.com/pybox2d/pybox2d.
* Hernández-Lobato, José Miguel, Hoffman, Matthew W, and Ghahramani, Zoubin. Predictive entropy search for efficient global optimization of black-box functions. In Advances in Neural Information Processing Systems (NIPS), 2014. https://bitbucket.org/jmh233/codepesnips2014
* Hennig, Philipp and Schuler, Christian J. Entropy search for information-efficient global optimization. Journal of Machine Learning Research, 13:1809–1837, 2012. http://www.probabilistic-optimization.org/Global.html
* Jarno Vanhatalo, Jaakko Riihimäki, Jouni Hartikainen, Pasi Jylänki, Ville Tolvanen, Aki Vehtari. GPstuff: Bayesian Modeling with Gaussian Processes. Journal of Machine Learning Research, 14(Apr):1175-1179, 2013.
* 