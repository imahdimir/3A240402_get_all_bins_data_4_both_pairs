"""


    """

import pandas as pd
import sys
from pathlib import Path

sys.path.append(Path.cwd().as_posix())

from prj.lib import BINS
from prj.lib import f
from prj.lib import fp
from prj.lib import v

##
def filter_then_random_draw_snps() :
    """    """

    ##
    df0 = pd.read_parquet(f.mrgd_mfi)
    _df = df0.head()

    ##
    # to not reading again from disk for subsequent runs
    df = df0.copy()

    ##
    # keep snps with maf > 1%
    msk = df[v.maf_n].gt(.01)
    df1 = df[msk]

    ##
    df = df1.copy()

    ##
    df[v.msk] = False

    for b in BINS :
        df[v.ms1] = df[v.info_n].ge(b) & df[v.info_n].lt(b + .01)

        _df = df[df[v.ms1]]

        if len(_df) > 1000 :
            _df = _df.sample(1000)

        df.loc[_df.index , v.msk] = True

    ##
    df = df[df[v.msk]]

    ##
    df = df.drop(columns = [v.msk , v.ms1])
    df.columns = list(range(8)) + [v.chr]

    ##
    print('Number of Total Filtered SNPs' , len(df))

    ##
    def save_subdf_to_txt(df , cn) :
        """ """
        print(cn)
        df = df.iloc[: , :8]
        print(len(df))
        _fn = fp.flt_snps.format(chr = cn)
        df.to_csv(_fn , sep = '\t' , index = False , header = False)

    ##
    gps = df.groupby([v.chr])
    gps.apply(lambda x : save_subdf_to_txt(x , x.name))

##
def main() :
    pass

    ##
    filter_then_random_draw_snps()

    ##

##
if __name__ == '__main__' :
    main()
    print(Path(__file__).relative_to(Path.cwd()) , ' Done.')

def testing_area() :
    """ """

    ##
    df = pd.read_parquet(f.mrgd_mfi)

    ##
    df.columns = list(range(8)) + ['chr']

    ##
    df.to_parquet(f.mrgd_mfi , index = False)

    ##
