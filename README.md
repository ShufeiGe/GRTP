# GRTP V.0.0.


The package GRTP implements the Bayesian nonparametric methods described in Ge et al., Generalized Random Tessellation Processes. This software constructs a generalized random tessellation processes (forest) for posterior prediction of regression data based on real-valued predictors. For the categorical tasks, we refer readers to the other software tess19 available at https://github.com/ShufeiGe/Tess19. For regression tasks, we apply a two-stage
inference scheme to infer the model to improve computational efficiency, the tree topology of the process is first inferred by an SMC algorithm, and then an MCMC algorithm is implemented to approximate the posterior distribution of parameters associated with each terminal node given the tree topology.  The source, manual and build instructions for the package GRTP are provided in this repository. This software requires the following R packages: optparse, purrr, MCMCPack. This software is released under the open source BSD 2-clause license.


## files
- **LICENSES**    The  *GRTP* software license.
- **demo.R** A simple demonstration to implement the proposed methods.
- **code**  A code folder contains the functions of the implemented methods.
- **code/generative_functins**  The implementation of the generative process of the random tessellation process for the regression models.
- **code/tess_frame_tau.R**   Model inference framework.
- **code/parameter_inference_l.R**   Parameter updates for linear cases.
- **code/parameter_inference_c.R**   Parameter updates for constant cases.
- **data/01**   Data used in the demonstration (demo.R).
 

## Citation
If you use GRTP in your research, please cite the following publication:

S. Ge, S. Wang, L. Wang, Y.W. Teh, L.T. Elliott. Generalized Random Tessellation Processes.

## LICENSES
GRTP v1.0. Copyright (c) 2024. Shufei Ge and Lloyd T. Elliott.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



## Help
Please send bugs, feature requests and comments to geshf@shanghaitech.edu.cn or shufeig@sfu.ca
