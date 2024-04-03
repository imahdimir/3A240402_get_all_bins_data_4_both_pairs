"""


    """

import itertools
import numpy as np
import pandas as pd
from bgen_reader import open_bgen
from pathlib import Path

from lib.lib import BINS , g , s , f , fp , v

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
    gps = df.groupby([v.chr])
    gps.apply(lambda x : save_subdf_to_txt(x , x.name))

    ##

##
def save_subdf_to_txt(df , cn) :
    """ """
    print(cn)
    df = df.iloc[: , :8]
    print(len(df))
    _fp = fpt.snps.format(cn)
    df.to_csv(_fp , sep = '\t' , index = False , header = False)

##
def separate_snps_by_chr() :
    """ """

    ##
    for i in range(1 , 22 + 1) :
        snps_fp = dyr.unpkd_mfi / f'ukb_mfi_chr{i}_v3.txt'

    ##
    df = pd.read_csv(fp.flt_snps_txt , sep = '\t' , header = None)

    ##
    for i in range(1 , 22 + 1) :
        msk = df[0].str.split(':').str[0]

        print(msk)
        # _df = df[msk]
        # _df.to_csv(dyr.med / f'snps_chr{i}.txt' ,
        #            sep = '\t' ,
        #            index = False ,
        #            header = False)
        break

##
def filter_snps_second_approach() :
    """ this approach is for using with bgen files manually, didn't work """

    ##
    df0 = pd.read_parquet(fp.snps)
    _df = df0.head()

    ##
    # to not reading again from disk for subsequent runs
    df = df0.copy()

    ##
    # keep snps with maf > 1%
    msk = df[mc.maf].gt(.01)
    df = df[msk]

    ##
    df[v.msk] = False

    for b in bins.values() :
        df[v.ms1] = df[mc.info].gt(b) & df[mc.info].le(b + .01)
        _df = df[df[v.ms1]]
        if len(_df) > 1000 :
            _df = _df.sample(1000)

        df.loc[_df.index , v.msk] = True

    ##
    df = df[df[v.msk]]

    ##
    for b in bins.values() :
        msk = df[mc.info].gt(b) & df[mc.info].le(b + .01)
        df.loc[msk , v.bin] = b

    ##
    df.to_parquet(fp.flt_snps_1 , index = False)

##
def filter_fs_and_po() :
    """ """

    ##
    df = pd.read_csv(fp.rel , sep = '\s+' , dtype = 'string')

    ##
    types = {
            'full_sibs'        : 'FS' ,
            'parent_offspring' : 'PO'
            }

    msk = df[v.inftype].isin(types.values())

    df = df[msk]

    ##
    df1 = df[['FID1' , 'ID1']]
    df2 = df[['FID2' , 'ID2']]

    df1.columns = ['FID' , 'ID']
    df2.columns = ['FID' , 'ID']

    df = pd.concat([df1 , df2])
    df = df.iloc[: , :2]

    ##
    df.to_csv(fp.fs_po_ids_txt , index = False , header = False , sep = '\t')

##
def filter_fs() :
    """ """

    ##
    df = pd.read_csv(fp.rel , sep = '\s+' , dtype = 'string')

    ##
    types = {
            'full_sibs' : 'FS' ,
            }

    msk = df[v.inftype].isin(types.values())

    df = df[msk]

    ##
    df1 = df[['FID1' , 'ID1']]
    df2 = df[['FID2' , 'ID2']]

    df1.columns = ['FID' , 'ID']
    df2.columns = ['FID' , 'ID']

    df = pd.concat([df1 , df2])
    df = df.iloc[: , :2]

    ##
    df.to_csv(fp.sibs_ids_txt , index = False , header = False , sep = '\t')

##
def make_df_of_iids_from_bgen_open_obj(bg_opn) :
    """ """

    ##
    bg_opn_s = list(bg_opn.samples)
    df = pd.DataFrame({
            'IID' : bg_opn_s
            })

    # get the IID from FID_IID
    df['IID'] = df['IID'].str.split('_').str[1]

    ##
    return df

##
def open_bgen_ret_iid_df_and_prob_arr(bgen_fp) :
    """ """

    ##
    bg_opn = open_bgen(bgen_fp)

    ##
    df_id = make_df_of_iids_from_bgen_open_obj(bg_opn)

    ##
    nd_p = bg_opn.read()

    ##
    return df_id , nd_p

##
def save_dosages_of_all_vars_from_bgen(bgen_fp: Path) :
    """ """

    ##
    df_id , nd_p = open_bgen_ret_iid_df_and_prob_arr(bgen_fp)

    ##
    nd_d = nd_p[: , : , 1] + 2 * nd_p[: , : , 2]

    ##
    df1 = pd.DataFrame(nd_d)

    ##
    df_d = pd.concat([df_id , df1] , axis = 1)

    ##
    _fp = dyr.med / f'dosages_{bgen_fp.stem}.prq'
    df_d.to_parquet(_fp , index = False)

##
def save_hard_calls_of_all_vars_from_bgen(bgen_fp) :
    """ """

    ##
    df_id , nd_p = open_bgen_ret_iid_df_and_prob_arr(bgen_fp)

    ##
    nd_h = np.argmax(nd_p , axis = 2)

    ##
    df1 = pd.DataFrame(nd_h)

    ##
    df_h = pd.concat([df_id , df1] , axis = 1)

    ##
    _fp = dyr.med / f'hard_calls_{bgen_fp.stem}.prq'
    df_h.to_parquet(_fp , index = False)

##
def gat_pairs_ids_in_pairs(identifier) :
    """ """

    df = pd.read_csv(fp.rel , sep = '\s+' , dtype = 'string')

    msk = df['InfType'].eq(identifier)

    df = df[msk]

    df = df[['ID1' , 'ID2']]

    return df

##
def make_prq_fp(gts_type , info_score , pair_suf) :
    """ """
    prq_fp = dyr.med / f'{gts_type}_snps_{info_score}{pair_suf}.prq'
    return prq_fp

##
def make_pairs_gts_dfs(df_gts , df_pairs_ids) :
    """ """

    ##
    dfa = pd.merge(df_pairs_ids[['ID1']] ,
                   df_gts ,
                   left_on = 'ID1' ,
                   right_on = 'IID' ,
                   how = 'left')
    dfb = pd.merge(df_pairs_ids[['ID2']] ,
                   df_gts ,
                   left_on = 'ID2' ,
                   right_on = 'IID' ,
                   how = 'left')

    ##
    dfa = dfa.drop(columns = ['ID1'])
    dfb = dfb.drop(columns = ['ID2'])

    ##
    return dfa , dfb

##
def save_corr_of_sib_pairs_with_gts_df(pair_identifier ,
                                       gts_type ,
                                       info_score ,
                                       pair_type ,
                                       pair_suf
                                       ) :
    """ """

    ##
    df_pairs_ids = gat_pairs_ids_in_pairs(pair_identifier)

    ##
    prq_fp = make_prq_fp(gts_type , info_score , pair_suf)
    df_gts = pd.read_parquet(prq_fp)

    ##
    dfa , dfb = make_pairs_gts_dfs(df_gts , df_pairs_ids)

    ##
    gts1 = dfa.iloc[: , 1 :]
    gts2 = dfb.iloc[: , 1 :]

    ##
    df_cors = gts1.corrwith(gts2 , method = 'pearson')

    ##
    out_fp = dyr.out_dta / f'corr_{pair_type}_{gts_type}_{info_score}.xlsx'
    df_cors.to_excel(out_fp , index = False , header = False)
    print(out_fp)

##
def main() :
    pass

    ##
    info_scores = {
            0 : 30 ,
            1 : 99 ,
            }

    pairs = {
            'sibs'             : ('' , 'FS') ,
            'parent_offspring' : ('_po' , 'PO') ,
            }

    gts_types = {
            0 : 'dosages' ,
            1 : 'hard_calls'
            }

    ##
    prd = itertools.product(info_scores.values() , pairs.values())

    for info , pair in prd :
        print(info , pair)
        fp = dyr.plink_out / f'snps_{info}{pair[0]}.bgen'
        print(fp)

        save_dosages_of_all_vars_from_bgen(fp)

        save_hard_calls_of_all_vars_from_bgen(fp)

    ##
    prd = itertools.product(info_scores.values() ,
                            pairs.keys() ,
                            gts_types.values())

    for info , pair_type , gm in prd :
        pair_iden = pairs[pair_type][1]
        print(info , pair_type , gm)
        pair_suf = pairs[pair_type][0]
        save_corr_of_sib_pairs_with_gts_df(pair_iden ,
                                           gm ,
                                           info ,
                                           pair_type ,
                                           pair_suf)

    ##

    ##

    ##

    ##

    ##

def testing_area() :
    """ """

    ##
    df = pd.read_parquet(f.mrgd_mfi)

    ##
    df.columns = list(range(8)) + ['chr']

    ##
    df.to_parquet(f.mrgd_mfi , index = False)

    ##
    for b in BINS :
        print(b)

    ##
    import sys

    ##
    a = sys

    ##
    from os import environ

    ENV = dict(environ)

    with open(f"{ENV['HOME']}/.export") as f :
        script = f.read()

    script = script.encode()

    ##
    import subprocess

    with subprocess.Popen(['sh'] ,
                          stdin = subprocess.PIPE ,
                          stdout = subprocess.PIPE) as p :
        result = p.communicate(script)

    ##
    result

    ##
    for line in result.splitlines() :
        var , _ , value = line.partition('=')
        print(var , value)
        os.environ[var] = value

    ##
    import subprocess
    import os

    ENV = dict(os.environ)

    filename = f"{ENV['HOME']}/.export"

    ##
    with open(filename) as f :
        script = f.read()
    script = script.encode() + b'\nenv'

    ##
    with subprocess.Popen(['sh'] ,
                          stdin = subprocess.PIPE ,
                          stdout = subprocess.PIPE) as p :
        result = p.communicate(script)

    result

    ##
    for line in result[0].decode().splitlines() :
        var , _ , value = line.partition('=')
        print(var , value)
        os.environ[var] = value

    ##
    os.environ['LC']

    ##


    if e :
        raise StandardError('conf error in {}: {}'.format(filename , e))

    for token in shlex.split(o) :
        parts = token.split('=' , 1)
        if len(parts) == 2 :
            os.environ[parts[0]] = parts[1]


    ##
    o , e = subprocess.Popen(['/bin/bash' , '--restricted' , '--noprofile' ,
                              '--init-file' , filename , '-i' , '-c' ,
                              'declare'] ,
                             env = {
                                     'PATH' : ''
                                     } ,
                             stdout = subprocess.PIPE ,
                             stderr = subprocess.PIPE).communicate()

    print(e)

    ##

    ##
    import importlib

    importlib.reload(os)

    ENV = dict(os.environ)

    ENV['DB']



    ##
    if e :
        raise Exception('conf error in {}: {}'.format(filename , e))

    for token in shlex.split(o) :
        parts = token.split('=' , 1)
        if len(parts) == 2 :
            os.environ[parts[0]] = parts[1]

    ##
    with subprocess.Popen(['sh'] ,
                          stdin = subprocess.PIPE ,
                          stdout = subprocess.PIPE) as p :
        result = p.communicate(script)

    result
    ##
    result[]

    ##
    for line in result[0].decode().splitlines() :
        var , _ , value = line.partition('=')
        print(var , value)
        os.environ[var] = value

    ##
    from pathlib import Path

    x = os.environ['DB']
    Path(x).exists()

    ##
    import os
    x = os.system('zsh -i -c "source ~/.export"')

    ##
    def parse_bash_script(fn) :
        with open(fn) as f :
            for line in f :
                if not line.startswith('#') :  # ignore comments
                    if "export" in line :
                        var , _ , val = line.partition('=')
                        var = var.lstrip()
                        val = val.rstrip()
                        if var == 'DB':
                            print(val)
                        if val.startswith('"') :
                            vals = val.rpartition('"')
                            val = vals[0][1] + vals[2]
                        os.environ[var] = val

    ##
    parse_bash_script(f"{ENV['HOME']}/.export")

    ##
    os.environ['DB']


    ##
    from os import environ

    ENV = dict(environ)
##
ENV['DB']

##
import dotsi
d = dotsi.Dict({"foo": {"bar": "baz"}})

##
d.foo.bar

##
d.extend({"a": 1, "b": 2})

##
from os import environ
import dotsi
from pathlib import Path
import subprocess

def read_env() :
    ENV = dict(environ)

    if 'DB' in ENV :
        print(ENV['DB'], 'there is')
        return dotsi.fy(ENV)

    print('no DB')
    fn = Path(f"{ENV['HOME']}/.export")
    if not fn.exists() :
        fn = Path(f"/homes/nber/mahdimir/.export")

    with open(fn) as f :
        script = f.read()
    script = script.encode() + b'\nenv'

    with subprocess.Popen(['sh'] ,
                          stdin = subprocess.PIPE ,
                          stdout = subprocess.PIPE) as p :
        result = p.communicate(script)

    for line in result[0].decode().splitlines() :
        var , _ , value = line.partition('=')
        print(var , value)
        environ[var] = value

    return dotsi.fy(dict(environ))

##
x = read_env()

##
x.DB

##
