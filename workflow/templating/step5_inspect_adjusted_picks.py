import os, logging, warnings, glob
from pathlib import Path

from pyrocko.gui.snuffler.marker import PhaseMarker

from eqcutil.util.logging import basic_logger_config, CriticalExitHandler
from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.viz.eqc_compat import phase_to_pick, pick_to_phase


basic_logger_config(level=logging.INFO)
Logger = logging.getLogger(os.path.split(__file__)[-1])
# ch = logging.StreamHandler()
Logger.addHandler(CriticalExitHandler())

# Ignore FutureWarning
warnings.simplefilter(action='ignore', category=FutureWarning)


# MAP ABSOLUTE PATHS
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Processed Data Directory
PD_DIR = ROOT / "processed_data" / "workflow"
# Template Save Directory
CTSD = PD_DIR / 'templates' / 'single_station' / 'xcc_test'

### USER PARAMETERS ###
cct = 0.5
show_master = True
MIN_MEMBERS = 3

flist = glob.glob(str(CTSD/'*.tgz'))
# Iterate across tribes
for f in flist:
    # Load
    Logger.warning(f'loading {f}')
    ctr = ClusteringTribe().read(f)
    # Iterate across groups
    ser_grp = ctr._c.xcc.value_counts()
    for _gn, _ct in ser_grp.items():
        # Check if group meets minimum membership
        if _ct < MIN_MEMBERS:
            continue
        else:
            pass

        # Subset into group tribe
        _ctr = ctr.get_subset(ctr._c[ctr._c.xcc==_gn].index)
        # Get master template
        Logger.warning(f'fetching master trace based on SNR')
        df_c = _ctr.clusters
        master_row = df_c[df_c.snr == df_c.snr.max()]
        master_name = master_row.index[0]
        Logger.info(f'Master Template is {master_name}')
        master_tmp = _ctr[master_name]
        # Get best-fit shifts from shift_mat from the whole tribe (not the subset)
        shift_line = ctr.shift_mat[master_row.id_no[0], :]
        # Make pick markers for shifted picks
        shifted_markers = []
        altnames = []
        for _e, _tmp_name in enumerate(_ctr._c.index):
            tmp = _ctr[_tmp_name]
            shift_row = _ctr._c.loc[_tmp_name]
            if _tmp_name == master_name:
                Logger.info('making MASTER template marker')
                kind = 1
                altnames.append(f'Grp{shift_row.xcc} - MASTER - {shift_row.etype}')
            else:
                kind = 5
                altnames.append(f'Grp{shift_row.xcc} - {shift_row.etype}')

            pick = tmp.event.picks[0].copy()
            # Use the id_no to pull the specific shift
            shift = ctr.shift_mat[master_row.id_no[0], _ctr._c.loc[_tmp_name,'id_no']]

            pick.time += shift
            pick.evaluation_mode='automatic'
            pm = pick_to_phase(pick, kind=kind)
            shifted_markers.append(pm)

        _ctr.snuffle(altnames=altnames, markers=shifted_markers)        
        breakpoint()
        # # Visualize with snuffler to assess master pick position
        # hdr = [master_name, df_c.loc[master_name, 'etype'], f'Grp{df_c.loc[master_name,"xcc"]}']
        # rtag, markers = master_tmp.snuffle(altname='-'.join(hdr))
        # for marker in markers:
        #     if isinstance(marker, PhaseMarker):

        # breakpoint()
        # Cross check outputs to see if 

