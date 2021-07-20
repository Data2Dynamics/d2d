# PEtabMEIGO

PEtab interface for MEIGO global optimization suite.

## Getting Started

PEtabMEIGO last version is included in MEIGO suite.

### Prerequisites

MATLAB interface for libSBML is the only dependency required and it's up to the user it's installation.

Latest libsSBML MATLAB interface can be downloaded from [this site](https://sourceforge.net/projects/sbml/files/libsbml/MATLAB%20Interface/).

### Installing

As mentioned, PEtabMEIGO it's supplied with MEIGO suite.

## How it works?

Given the directory and model name of a problem described using PEtab format, or the absolute or relative
path pointing to problem's **.yaml** file, PEtabMEIGO function, **getMeigoProblem** returns a struct 
containing the following fields:

* f : Problem's objective function.
* x_L : Optimization parameters lower bounds, unscaled.
* x_U : Optimization parameters upper bounds, unscaled.
* x_0 : Optimization parameters initial values, unscaled.

### Example of use

With problem directory **problem_dir** and model name **model_name** do:

```
    out = getMeigoProblem(problem_dir, model_name);
```

or we can also use problem's **.yaml** file path, **yaml_path**:

```
    out = getMeigoProblem(yaml_path);
```

### Known Issues

* TODO

## Version

* *PEtabMEIGO 20200401*

## Authors

* **Tacio Camba Espí**

## Acknowledgments

* TODO