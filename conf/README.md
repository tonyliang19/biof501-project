# Configuration Profiles

This directory stores possible configuration profiles for the pipeline, which allows the flexibility of switching `profile` by simple cli arg, and uses complete different sets of parameters, resource configurations.

I.e switching or even chaining profiles like the following:

```bash

# This is to use test profile, which limits amount of resources available upon running
nextflow run main.nf -profile docker,test
# NOTE: this one below is not implemented
# This could be storing important settings on the compute canada cluster
nextflow run main.nf -profile docker,compute_canada
```

For now, there are two configurations here, the `base` which is sourced by all other profiles to get basic resource settings for the processes. And, there's a `test` config, which represents minimal resource setting to run the pipeline.

In other words, if `test` was not used upon running nextflow, then it uses settings from the `base` config. This gives flexibility to override profiles and have self-contained configurations to repeatedly run the pipeline under certain hyperparameters if any.