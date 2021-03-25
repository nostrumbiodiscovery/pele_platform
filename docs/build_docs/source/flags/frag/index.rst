FragPELE parameters
===================

These flags are **exclusive to FragPELE** modes.

- **growing_steps**: Number of steps to grow the fragment with.

- **steps_in_gs**: Number of pele steps within each growing step

- **sampling_steps**: Number of pele steps in the final sampling simulation

- **protocol**: Type of protocol. options = [HT, ES]. For more info please refere here.


..  code-block:: yaml

    growing_steps: 6
    steps_in_gs: 6
    sampling_steps: 20
    protocol: HT
    cpus: 24