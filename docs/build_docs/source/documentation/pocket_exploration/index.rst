Allosteric parameters
======================

These flags are exclusive to the ``allosteric: true`` mode.

- **n_components**: Number of clusters after global exploration. In other words, number of inputs for the refinement exploration after the global simulation. Default: 10
- **skip_refinement**: To skip refinement simulation. Default: False

..  code-block:: yaml

    n_components: 10
    skip_refinement: true