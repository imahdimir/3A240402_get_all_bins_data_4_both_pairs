import numpy as np
from pathlib import Path
from mahdi_env import get_env

PROJ = 'get-all-bins-data-4-both-pairs-240402'

BINS = np.arange(.3 , 1 , .01)

env = get_env()

class Env :
    pgsf = 'pgsf'

e = Env()

class Dir :
    pass

g = Dir()

class SFDir :
    sfn = PROJ + '-SF'
    slf = Path(env[e.pgsf]) / sfn
    inp = slf / 'inp'
    med = slf / 'med'
    out = slf / 'out'
    flt_snps = med / 'flt-snps'

s = SFDir()

class File :
    mrgd_mfi = '/disk/genetics4/ws_dirs/mahdimir/lc/non-proj/ukb/v3/mrgd-mfi.prq'
    # fs and po pairs iids - from "imputed-genotype-corr-240317" proj
    fs_po_ids = s.inp / 'fs-po-ids.txt'

f = File()

class FilePattern :
    flt_snps = s.flt_snps.as_posix() + '/c{chr}.txt'

fp = FilePattern()

class Vars :
    maf_n = '5'
    info_n = '7'
    rsid_n = '1'

    msk = 'msk'
    ms1 = 'ms1'
    chr = 'chr'

v = Vars()
