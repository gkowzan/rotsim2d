"""Profile performance of new classes and database."""
# * Imports
import cProfile
import pstats
from pathlib import Path
from sqlalchemy import create_engine, select
from sqlalchemy.orm import aliased
import matplotlib.pyplot as plt
import molspecutils.happier as h
import molspecutils.molecule
from molspecutils.molecule import CH3ClAlchemyMode, SymTopState
from molspecutils.alchemy import CH3Cl
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.visual as vis
import rotsim2d.propagate as prop
import gkpywigxjpf as wig
plt.ion()

# * Vibrational mode
sql_path = Path(h.hitran_cache) / 'CH3Cl.sqlite3'
engine = create_engine("sqlite:///" + str(sql_path))
ch3cl_mode = CH3ClAlchemyMode(engine)
T = 296.0

# * Pathways
pws = pw.gen_pathways(range(1, 37), [0]*4, meths=[pw.only_SII], rotor='symmetric', klimit=10)
dressed_pws = dl.dress_pws(pws, ch3cl_mode, T)

peaks = dl.peak_list(dressed_pws)

# * Visualize
vis.plot2d_scatter(peaks)

# * Profile
cProfile.run('dl.dress_pws(pws, ch3cl_mode, T)', 'ch3clmode_stats.cprofile')
p = pstats.Stats('ch3clmode_stats.cprofile')
p.strip_dirs().sort_stats(pstats.SortKey.CUMULATIVE).print_stats(30)

# * Test
pair = (SymTopState(nu=0, j=1, k=0), SymTopState(nu=1, j=0, k=0))
rovpp, rovp = aliased(CH3Cl.RovibState), aliased(CH3Cl.RovibState)

subjpp = select(CH3Cl.RotState.id).filter_by(
    j=pair[0].j, k=pair[0].k, l=0, f=0.0
).subquery()
subnupp = select(CH3Cl.VibState.id).filter_by(
    nu1=0, nu2=0, nu3=pair[0].nu, nu4=0, nu5=0, nu6=0
).subquery()
subjp = select(CH3Cl.RotState.id).filter_by(
    j=pair[1].j, k=pair[1].k, l=0, f=-0.0
).subquery()
subnup = select(CH3Cl.VibState.id).filter_by(
    nu1=0, nu2=0, nu3=pair[1].nu, nu4=0, nu5=0, nu6=0
).subquery()
lp = CH3Cl.LineParameters

subq = select(lp.nu, lp.A, lp.gamma_air, lp.delta_air, 
              rovpp.j_id.label("jpp_id"), rovpp.nu_id.label("nupp_id"),
              rovp.j_id.label("jp_id"), rovp.nu_id.label("nup_id"))\
              .execution_options(yield_per=10)\
              .join(lp.transition)\
              .join(rovpp, rovpp.id==CH3Cl.TransitionPair.statepp_id)\
              .join(rovp, rovp.id==CH3Cl.TransitionPair.statep_id)\
              .where(CH3Cl.LineParameters.nu > 500.0, CH3Cl.LineParameters.nu < 800.0).subquery()
result = ch3cl_mode.sess.execute(
    select(subq.c.nu).join(subjpp, subjpp.c.id==subq.c.jpp_id).join(subjp, subjp.c.id==subq.c.jp_id).\
    join(subnupp, subnupp.c.id==subq.c.nupp_id).join(subnup, subnup.c.id==subq.c.nup_id)
)
