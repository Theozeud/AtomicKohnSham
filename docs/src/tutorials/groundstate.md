# Computation of the Groundstate

Once you have created your model, discretization and chose the algorithm as explained in the previous section, you can finally compute the groundstate :

```@docs
groundstate(model::KSEModel, discretization::KSEDiscretization, alg::SCFAlgorithm;
            name::String = "", kwargs...)
```

Or similarly, if you have assembly an `AtomProblem`, then you can give it to `groundstate` :

```@docs
groundstate(prob::AtomProblem)
```
