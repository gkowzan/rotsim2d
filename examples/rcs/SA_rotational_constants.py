import psi4
from pathlib import Path

SA_geometry_path = Path("~/notes/.attach/sb/:8c55d411-1a6f-4227-bf33-653ce4043c67/SA_geometry1.txt").expanduser()
with open(SA_geometry_path) as f:
    SA_geometry = psi4.geometry(f.read())
print(SA_geometry.rotational_constants().to_array())

SA_geometry_path = Path("~/notes/.attach/sb/:8c55d411-1a6f-4227-bf33-653ce4043c67/SA_geometry2.txt").expanduser()
with open(SA_geometry_path) as f:
    SA_geometry = psi4.geometry(f.read())
print(SA_geometry.rotational_constants().to_array())
