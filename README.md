Package to simulate 2D rovibrational spectra of gas-phase samples.

See examples in `examples` directory:

- `methyl_chloride/classify_pathways.py` prints out a mapping between
  polarization-dependence expressions and pathways they correspond to. Modify
  `gen_pathways` call to obtain classifications for different experimental
  schemes.
- `methyl_chloride/peak_picker_three_color.py` and
  `methyl_chloride/peak_picker_three_color.py` allow one to inspect 2D peak
  scatter plots. Clicking on a peak prints out double-sided Feynman diagrams
  that contribute to the selected peak. See comments in the code to modify what
  gets plotted.
- `methyl_chloride/polarization_dependence.py` allows one to inspect
  polarization dependence in the case where two beam polarizations are the same
  and in the case when all polarizations are different. See comments in the code
  to modify the suppressed pathways and excitation scheme.

# Dependencies
- numpy, scipy, anytree, sympy
- [gkpywigxjpf](https://gitlab.com/allisonlab/mdcs/pywigxjpf)
- [knickknacks](https://gitlab.com/allisonlab/mdcs/shed)
- [molspecutils](https://gitlab.com/allisonlab/mdcs/spectroscopy)

# Installation
Install the package from our private GitLab repository by executing:

``` sh
pip install rotsim2d --extra-index-url https://<token_name>:<token>@gitlab.com/api/v4/projects/26140156/packages/pypi
```

where `<token_name>` and `<token>` are obtained by creating a personal token
with `read_api` scope. See [Creating personal access
tokens](https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html#creating-a-personal-access-token)
for details.

[molspecutils](https://gitlab.com/allisonlab/mdcs/spectroscopy) package needs to
be initialized by executing `molspecutils_init` script to download spectroscopic
data and fill the local database.

# Development
Development environment is managed with `pip-tools`. 
