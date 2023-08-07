# DNN-Based-Algorithm
Matlab code to reproduce our DNN-Based Optimization algorithm for the Optimal Released Molecules for Mobile Molecular MIMO Communications. It includes:

- [dataset_generation.m](dataset_generation.m): the file to generate the initial dataset
- [diffusion_rx_runner.m](diffusion_rx_runner.m): the simulator runner
  - [rand_coordinate_generate.m](rand_coordinate_generate.m): the function utilize to find a point that is a certain distance from another point
  - [prepare_vars4_diffusion_runners_PointSrc.m](prepare_vars4_diffusion_runners_PointSrc.m): the simulation parameter of transmitters, receivers, and molecules are setted in this file
  - [CORE_sim_diffusion_3d_P2S_wAbsorption.m](CORE_sim_diffusion_3d_P2S_wAbsorption.m) the main file to simulate the random movements of transmitters, receivers, and molecules
- [Diversification_Generation.m](Diversification_Generation.m): the function is utilized to quantize the output of the DNN to K neighborhood solutions.
- [select_diversification_generation.m](select_diversification_generation.m): the function is utilized to select the best solution among all the neighborhood solutions quantized
  - [error_probability.m](error_probability.m): the formulation of the average BER
  - [probability_mobile.m](probability_mobile.m): the formulation of the channel impulse response
- [Queue.m](Queue.m): define queue structure
  - [Node.m](Node.m): define node of the queue
- [main.m](main.m): the main file to run the DNN-Based Optimization algorithm

# DNN-Based-Algorithm Diagram
![Algorithm](Algorithm(version4).jpg).
