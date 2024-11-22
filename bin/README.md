# Executable Binary Scripts

This directory contains R or Python scripts used in [modules](../modules/README.md) processes, such the process meta definition (input and output) and actual analysis/computation is separated.


Each script in this directory should follow this format:

1. Contains a shebang line in the first line of the this dummy script `hello.R` like:

    ```R
    #!usr/bin/env Rscript

    # The line above is shebang, which tells interpreter which command
    # should be used to execute the script, this case it is 'Rscript'
    # And this is ran inside a container that includes the desired binary
    dummy_hello <- paste("This is hello", "world", sep="_")
    ```

2. Make the script executable using `chmod +x hello.R`



This is based on Nextflow's functionality of [managing dependencies on `bin` directory](https://www.nextflow.io/docs/latest/sharing.html#the-bin-directory).