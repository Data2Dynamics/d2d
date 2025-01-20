# Examples for how to use model selection for d2d models

All examples are based on the ODE system and data from the test cases 0001 to 0006.
Here, however, we build the model selection problem based on a d2d implementation instead of a Eetab implementation.
All three examples are equivalent formulations of the same model selection problem. This should showcase, how one can generally apply the PEtab Select tools for d2d models.

## 1. Single subspace
In this example, we create a single d2d model that correpsonds to the single model subspace defined in the model_space.tsv table. This model subspace is in the compressed format (similar to test case 0003).

## 2. Two subspaces
In this example, we also create a single d2d model. But here, this model is used as the basis for two separate model subspaces in the model_space.tsv table. Specifically, we create on subspace without the annihilation (k3=0) and one subspace for fitting the annihilation rate (k3='estimate').

## 3. Two distinct subspaces
In this example, we again create two model subspaces with the same meaninga s before. But now, the presence of the annihilation mechanism is encoded within the d2d model definition files (model1.def and model2.def). Therefore, we have two distinct workspaces as the bases for the two model subspaces.