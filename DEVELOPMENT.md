Development environment is managed with `pip-tools`.
Sync current env with runtime, dev and interactive work requirements by calling `pip-refresh.sh`.
This should be done after adding a new dependency.
It won't update any package versions unless necessary.
`pip-refresh force` will remove current reqs and create them from scratch, most likely installing most recent deps.

Local development of rotsim2d, rotsim2d_apps and molspecutils is supported by `local_update.sh` script located at each project's root.
It will build the package, upload it to a local pypi server and install it in the projects that depend on it.
