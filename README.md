# Max-value Entropy Search
This is the MATLAB code repository associated with the paper [_Max-value Entropy Search for Efficient Bayesian Optimization_](https://arxiv.org/abs/1703.01968). We propose a new Bayesian optimization technique called **Max-value Entropy Search**, which maximizes the multual information between the selected future observations and the max-value of the function. Please refer to the paper if you need more details on the algorithm.

We developed our code building upon some basic Gaussian process functionals from predictive entropy search (Hernandez-Lobato et al., 2014), which was developed upon GPstuff (Vanhatalo et al., 2013). 

If you have any question, please create an issue [here](https://github.com/zi-w/Max-value-Entropy-Search/issues/new) and remind me via email wangzi@google.com if no reply in a week.

## System Requirement
We tested our code with MATLAB R2015b on Ubuntu 14.04 LTS (64-bit). Our code should work out of box without installing additional packages. However, please make sure you installed the GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)). On Ubuntu, you can install GSL by 

```
sudo apt-get install libgsl0-dev
```

In MATLAB command line, you can mex the c files, for example, by

```
mex chol2invchol.c -lgsl -lblas
```

To run Box2D related code in test-functions/robot-pushing/, please install [Pybox2d](https://github.com/pybox2d/pybox2d).

## Example
example.m is a simple example using Bayesian optimization to maximize a black-box function. Please see the comments in the code for more details about the usage.

gpopt.m is the function for Bayesian optimization with Gaussian processes.

add_gpopt.m is the function for Bayesian optimization with additive Gaussian processes, and is more suitable for high dimensional problems.

## Method Options
In addition to the proposed Max-value Entropy Search method, we also provide several other popular Bayesian optimization methods for the convenience of other researchers. The following is a full list of Bayesian optimization methods that are implemented in this repository.

1. Max-value Entropy Search with Gumbel sampling (MES-G) by Wang & Jegelka, 2017;
2. Max-value Entropy Search with random features (MES-R) by Wang & Jegelka, 2017;
3. Optimization as estimation (EST) by Wang et al., 2016. 
4. Gaussian process upper confidence bound (GP-UCB) by Auer, 2002; Srinivas et al., 2010;
5. Probability of improvement (PI) by Kushner, 1964;
6. Expected improvement (EI) by Mockus, 1974


In the paper, we also used the open sourced code from [entropy search](http://www.probabilistic-optimization.org/Global.html) (Hennig & Schuler, 2012) and [predictive entropy search](https://bitbucket.org/jmh233/codepesnips2014) (Hernández-Lobato et al., 2014) to compare with our method.

For high dimensional Bayesian optimization with additive Gausian processes, we implemented

1. Add-MES-G by Wang & Jegelka, 2017;
2. Add-MES-R by Wang & Jegelka, 2017;
3. Add-EST by Wang & Jegelka, 2017;
4. Add-GP-UCB by Kandasamy et al., 2015.

Our strategy for learning the Gaussian process additive structure is based on Wang et al., 2017.
## Example functionals for optimization
In test_functions/, we provide some functionals one can use to test Bayesian optimization algorithms. 

1. Approximated functions sampled from Gaussian processes: sample_GPprior.m, sample_addGP.m.
2. Optimization test functions: egg.m, michNf.m, shekelMf.m.
3. Tuning hyper-parameters for neural networks: nnet_boston.m, nnet_cancer.m.
4. Active learning for robot pushing: robot_pushing_3.m, robot_pushing_4.m, two_robot_pushing.m.
5. Tuning the walking speed of a planar bipedal robot: walker_speed.m (Westervelt el al., 2007).

## Citation
Please cite our work if you would like to use the code.
```
@inproceedings{wang2017maxvalue,
  title={Max-value Entropy Search for Efficient Bayesian Optimization},
  author={Wang, Zi and Jegelka, Stefanie},
  booktitle={International Conference on Machine Learning (ICML)},
  year={2017}
}
```
## References
* Wang, Zi and Jegelka, Stefanie. Max-value Entropy Search for Efficient Bayesian Optimization. International Conference on Machine Learning (ICML), 2017.
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
* Kandasamy, Kirthevasan, Schneider, Jeff, and Poczos, Barnabas. High dimensional Bayesian optimisation and bandits via additive models. In International Conference on Machine Learning (ICML), 2015.
* Wang, Zi, Li, Chengtao, Jegelka, Stefanie, and Kohli, Pushmeet. Batched High-dimensional Bayesian Optimization via Structural Kernel Learning. International Conference on Machine Learning (ICML), 2017.
* Westervelt, Eric R, Grizzle, Jessy W, Chevallereau, Christine, Choi, Jun Ho, and Morris, Benjamin. Feedback control of dynamic bipedal robot locomotion, volume 28. CRC press, 2007.
