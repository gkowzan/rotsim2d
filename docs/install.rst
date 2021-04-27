Installation
============

**Dependencies**

- numpy, scipy, anytree, sympy
- `gkpywigxjpf <https://gitlab.com/allisonlab/mdcs/pywigxjpf>`_
- `knickknacks <https://gitlab.com/allisonlab/mdcs/shed>`_
- `molspecutils <https://gitlab.com/allisonlab/mdcs/spectroscopy>`_

**Install**

.. highlight:: sh

Install the package from our private GitLab repository by executing::

  pip install rotsim2d --extra-index-url https://<token_name>:<token>@gitlab.com/api/v4/projects/26140156/packages/pypi

where `<token_name>` and `<token>` are obtained by creating a personal token
with `read_api`` scope. See `Creating personal access tokens
<https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html#creating-a-personal-access-token>`_
for details.

`molspecutils <https://gitlab.com/allisonlab/mdcs/spectroscopy>`_ package needs
to be initialized by executing `molspecutils_init` script to download
spectroscopic data and fill the local database.
