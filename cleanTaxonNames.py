import sys
from datetime import datetime
import os.path
import pandas as pd
import re
import numpy as np
import random
from ast import literal_eval
import csv
#
print("cleanTaxonNames.py: began at", datetime.now().strftime("%H:%M:%S"))
#
# !!! uncomment and use only if debugging !!!
# path = "/media/nlkhd/Nicole/Postdoc/PlantTaxonNameCleaner/"
# full_filename = "/media/nlkhd/Nicole/Postdoc/WorldCuP/Sources/Sweden/Sweden_SKUD_2021.csv"
# old_text_bool = False
# auth_split_bool = True
# cv_bool = True
#
full_filename = sys.argv[1]
filename = full_filename[:-4]
check_file = os.path.exists(full_filename)
if not check_file:
    print(full_filename)
    sys.exit("error: 'FilePath' file does not exist")
# path to Genera and Author aliases
path = sys.argv[2]
check_file = os.path.exists("".join(path + "/GeneraAliases.csv"))
if not check_file:
    print(path)
    sys.exit("error: directory with script files not correctly defined")
print("\nfilename:", full_filename)
def str_to_bool(s):
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
         raise ValueError("can't convert {} to a bool".format(s))
old_text_bool = str_to_bool(sys.argv[3])
if type(old_text_bool) != bool:
    sys.exit("error: 'OldText' argument must be Boolean")
print("\nold text? (replace authority abbrevs and move suspected cultivars):", old_text_bool)
auth_split_bool = str_to_bool(sys.argv[4])
if type(auth_split_bool) != bool:
    sys.exit("error: 'AuthorSplit' argument must be Boolean")
print("authorities listed in 'TaxonName' column?:", auth_split_bool)
cv_bool = str_to_bool(sys.argv[5])
if type(cv_bool) != bool:
    sys.exit("error: 'Cultivar' argument must be Boolean")
print("edit cultivar names?:", cv_bool)
#
# length of vectors shown in console is retricted
# large files will have extremely long, hard-to-read output otherwise
subset_len = 40
print("\nNOTE: console output will show a random subset of", subset_len, "matches to save space!\n\tto show more, edit 'subset_len' parameter in cleanTaxonNames.py")
#
# load functions
#
# replacing values in array based on dictionary, exact matches only
def replace_dict_exact(arr, dic):
    if type(arr) == pd.core.series.Series:
        arr = arr.to_numpy()
    # extract keys and values
    dic_keys = np.array(list(dic.keys()))
    dic_vals = np.array(list(dic.values()))
    # get argsort indices
    s_idx = dic_keys.argsort()
    s_dic_keys = dic_keys[s_idx]
    s_dic_vals = dic_vals[s_idx]
    idx = np.searchsorted(s_dic_keys, arr)
    idx[idx == len(s_dic_vals)] = 0
    mask = s_dic_keys[idx] == arr
    # show output of the values being changed
    rep_idx = [i for i, x in enumerate(mask) if x]
    replaced = arr[rep_idx]
    u_replaced = set(replaced)
    replacement = [dic[x] for x in u_replaced]
    if (len(u_replaced) > 0):
        print("\t", len(u_replaced), "unique values identified to be replaced with", len(replaced), "matches\n\t\treplacing\t", ", ".join(u_replaced), "\n\t\twith\t", ", ".join(replacement))
    else:
        print("\tnothing to be replaced")
    return np.where(mask, s_dic_vals[idx], arr)
# replace author names requires special regexp
def replace_author(arr, dic, alphbtz = False, allow_inits = False):
    if type(arr) == pd.core.series.Series:
        arr = arr.to_numpy()
    # extract keys and values
    dic_keys = np.array(list(dic.keys()))
    dic_vals = np.array(list(dic.values()))
    if alphbtz:
        # sort alphabetically
        s_idx = dic_keys.argsort()
        s_dic_keys = dic_keys[s_idx]
        s_dic_vals = dic_vals[s_idx]
    else:
        s_dic_keys = dic_keys
        s_dic_vals = dic_vals
    for y in range(len(s_dic_keys)):
        this_key_init = re.sub(r'\ ', r' ?', s_dic_keys[y])
        this_key_init = re.sub(r'\.', r'\.?', this_key_init)
        if allow_inits:
            this_key = "".join(['( |^|\(|\.)', this_key_init, '( |$|\)|,)'])
        else:
            this_key = "".join(['( |^|\()', this_key_init, '( |$|\)|,)'])
        this_val = "".join(['\\1', s_dic_vals[y], '\\2'])
        idx = [i for i, v in enumerate(arr) if bool(re.search(this_key, v))]
        if (len(idx) > 0):
            initial = arr[idx]
            initial_samp = set(initial)
            arr = np.array([re.sub(this_key, this_val, v) for v in arr])
            final = arr[idx]
            final_samp = filter(None, final)
            final_samp = set(final_samp)
            print(y, "\treplacing\t", ", ".join(initial_samp), "\n\twith\t", ", ".join(final_samp))
    return arr
# replace substring matches (as regexp) in an array using dictionary
def replace_dict_substr(arr, dic):
    if type(arr) == pd.core.series.Series:
        arr = arr.to_numpy()
    # extract keys and values
    dic_keys = np.array(list(dic.keys()))
    dic_vals = np.array(list(dic.values()))
    # get argsort indices
    s_idx = dic_keys.argsort()
    s_dic_keys = dic_keys[s_idx]
    s_dic_vals = dic_vals[s_idx]
    for y in range(len(s_dic_keys)):
        idx = [i for i, v in enumerate(arr) if s_dic_keys[y] in v]
        if (len(idx) > 0):
            print(y)
            initial = arr[idx]
            arr = np.array([re.sub(re.escape(s_dic_keys[y]), s_dic_vals[y], v) for v in arr])
            final = arr[idx]
            len_initial = len(initial)
            if len_initial > subset_len:
                idx = list(range(len(initial)))
                np.random.shuffle(idx)
                idx_samp = idx[1:subset_len]
                print(y, "\t(", len_initial, " matches) replacing\t", ", ".join(initial[idx_samp]), "\n\twith\t", ", ".join(filter(None, final[idx_samp])), sep = '')
            else:
                print(y, "\t(", len_initial, " matches) replacing\t", ", ".join(initial), "\n\twith\t", ", ".join(filter(None, final)), sep = '')
    return arr
# return equal-length array with only grep matches extracted
def separate_infra(arr, keep_final = True):
    regexp = '(?:^| )(var\.|subvar\.|convar\.|subsp\.|f\.|subf\.)'
    re_sub = '.*(^| )(var\.|subvar\.|convar\.|subsp\.|f\.|subf\.) ?'
    if type(arr) == pd.core.series.Series:
        arr = arr.to_numpy()    
    rank = []
    name = []
    for y in range(len(arr)):
        this_rank = re.findall(regexp, arr[y])
        this_name = re.sub(re_sub, '', arr[y])
        if (len(this_rank) == 0):
            rank.append('')
            name.append(this_name)
        elif (len(this_rank) == 1):
            rank.append(this_rank[0])
            name.append(this_name)
        elif (len(this_rank) > 1):
            if keep_final:
                print("\n\tmultiple matches detected (", ", ".join(this_rank), ") keeping final match (", this_rank[-1],")", sep = '')
                rank.append(this_rank[-1])
                name.append(this_name)
            else:
                print(arr[y])
                raise ValueError('multiple matches in extract output')
    return(rank, name)
# given input array, return array of matched strings ('string'), matched indices ('index'), or equal-length array of booleans indicating regexp match ('bool')
def return_matches(arr, regexp, mode = 'string'):
    check_mode = mode in ['string', 'bool', 'index']
    if not check_mode:
        print(check_mode)
        raise ValueError('mode must be string, bool, or index')
    if type(arr) == pd.core.series.Series:
        arr = arr.to_numpy()    
    check = []
    dummy = True
    if (mode == 'string'):
        for y in range(len(arr)):
            if (len(arr[y]) > 0):
                this_check = re.findall(regexp, arr[y])
                if (len(this_check)):
                    check.append(arr[y])
    elif (mode == 'bool'):
        for y in range(len(arr)):
            if (len(arr[y]) > 0):
                check.append(bool(re.search(regexp, arr[y])))
            else:
                # empty string, assume false
                check.append(False)
                if dummy:
                    print("NOTE: checking empty string and returning false (this error only arises once)\n\tindex\t", y)
                    dummy = False
    elif (mode == 'index'):
        for y in range(len(arr)):
            if (len(arr[y]) > 0):
                this_check = re.findall(regexp, arr[y])
                if (len(this_check)):
                    check.append(y)
    return(check)
# split taxon names in 'TaxonName' column into their parts
def split_taxon(this_taxon_init, authority = auth_split_bool, verbose = False):
    gen_idx = [0]
    spp_idx = []
    infra_rank_idx = []
    infra_idx = []
    auth_idx = []
    infra_auth_idx = []
    temp_auth_idx = []
    cv_idx = []
    this_gen = ''
    this_spp = ''
    this_infra_rank = ''
    this_infra = ''
    this_auth = ''
    this_infra_auth = ''
    temp_auth = ''
    this_cv = ''
    # split taxon name into words for parsing
    this_taxon = literal_eval(f'{this_taxon_init.split()}')
    # test whether species has a nothogenus
    # if there is a space separating the hybrid indicator, remove from list and add to genus
    gen_hyb_test1 = bool(re.search('^[xX\u00D7+]$', this_taxon[0]))
    gen_hyb_test2 = bool(re.search('^[\u00D7+]', this_taxon[0]))
    gen_hyb_test3 = bool(re.search('\+', this_taxon[0]))
    gen_hyb_test = gen_hyb_test1 | gen_hyb_test2
    if gen_hyb_test1:
        if verbose:
            print("\ttaxon is genus hybrid")
        this_taxon = this_taxon[1:]
    this_gen = re.findall('^[xX\u00D7+]?[A-Z]{1}[A-Za-z-]*$', this_taxon[0])
    if (len(this_gen) == 1):
        this_gen = this_gen[0]
        if gen_hyb_test1:
            if gen_hyb_test3:
                this_gen = "".join("+ " + this_gen)
            else:
                this_gen = "".join("\u00D7 " + this_gen)
    else:
        print(this_taxon_init)
        raise ValueError("genus not identified")
    # test whether taxon is a nothospecies
    # first determine whether infraspecific name or cultivar exists to set search range
    infrank_cv_idx = return_matches(arr = this_taxon, regexp = '^(notho)?subsp\.?$|^(notho)?ssp\.?$|^(notho)?(con|sub)?var\.?$|^(notho)?(con|sub)?v\.?$|^(notho)?(sub)?f\.?$|^(notho)?(sub)?f[ao]?\.?$|^(notho)?(sub)?forma\.?$|^(notho)?(sub)?fma\.?$|^[\u03B1\u03B2\u03B3\u03B4\u03B5\u03B6\u03B7\u03B8\u03BB\u03C1]$|^[\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4\'"]|^cv\.?$', mode = 'index')
    if len(infrank_cv_idx) > 0:
        taxon_range = this_taxon[1:infrank_cv_idx[0]]
    else:
        taxon_range = this_taxon[1:]
    spp_hyb_test1 = return_matches(arr = taxon_range, regexp = '^[xX\u00D7]$', mode = "index")
    spp_hyb_test2 = return_matches(arr = taxon_range, regexp = '^[\u00D7]', mode = "bool")
    spp_hyb_test = (len(spp_hyb_test1) > 0) | any(spp_hyb_test2)
    # check whether taxon has been identified to species level
    if len(this_taxon) == 1:
        check_ident_sp = True
        check_ident_cv = False
        this_taxon.append('sp.')
    else:
        check_ident_sp = bool(re.search(r"^spp?.?$", this_taxon[1]))
        check_ident_cv = bool(re.search(r"^[\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4'\"][A-Za-z-0-9À-ÿĀ-ſƀ-ȳ.'\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4]+", this_taxon[1]))
    if (check_ident_sp | check_ident_cv):
        if check_ident_sp:
            if verbose:
                print("\tspecific epithet not defined (sp.)")
            # specific epithet not defined (sp.)
            this_spp = this_taxon[1]
            spp_idx = [1]
        elif check_ident_cv:
            if verbose:
                print("\tspecific epithet not defined (cv only)")
            # specific epithet not defined (cultivar only)
            this_taxon.insert(1, 'sp.')
            this_spp = this_taxon[1]
            spp_idx = [1]
    else:
        # specific epithet is defined
        if (len(spp_hyb_test1) > 0):
            # taxon is a hybrid
            if (len(spp_hyb_test1) == 1):
                if (spp_hyb_test1[0] == 0):
                    if verbose:
                        print("\ttaxon is a hybrid 'Genus x specific epithet'")
                    # check whether there is just a space before standard hybrid name
                    spp_idx.append(spp_hyb_test1[0] + 1)
                    spp_idx.append(spp_hyb_test1[0] + 2)
                else:
                    if verbose:
                        print("\ttaxon is a hybrid 'Genus specific epithet 1 x specific epithet 2'")
                    # or whether hybrid name shows parent epithets
                    spp_idx.append(spp_hyb_test1[0])
                    spp_idx.append(spp_hyb_test1[0] + 1)
                    spp_idx.append(spp_hyb_test1[0] + 2)
            elif (len(spp_hyb_test1) == 2):
                if (spp_hyb_test1[0] == 0):
                    if verbose:
                        print("\ttaxon is a hybrid 'Genus x specific epithet 1 x specific epithet 2'")
                    # check whether there is just a space before standard hybrid name
                    spp_idx.append(spp_hyb_test1[0] + 1)
                    spp_idx.append(spp_hyb_test1[1])
                    spp_idx.append(spp_hyb_test1[1] + 1)
                    spp_idx.append(spp_hyb_test1[1] + 2)
                else:
                    if verbose:
                        print("\ttaxon is a hybrid 'Genus specific epithet 1 x specific epithet 2 x specific epithet 3'")
                    # species has multiple parent epithets
                    spp_idx.append(spp_hyb_test1[0])
                    spp_idx.append(spp_hyb_test1[0] + 1)
                    spp_idx.append(spp_hyb_test1[1])
                    spp_idx.append(spp_hyb_test1[1] + 1)
                    spp_idx.append(spp_hyb_test1[1] + 2)
            else:
                # other hybrid configurations not allowed at the moment
                print(this_taxon_init)
                raise ValueError("hybrid species epithet non-parsable")
            this_spp_init = [e for i, e in enumerate(this_taxon) if i in spp_idx]
            this_spp = ' '.join(this_spp_init)
            if verbose:
                print("\thybrid species epithet:", this_spp)
        else:
            # taxon is not a hybrid
            # species epithet should be in second position, check this
            # epithet should be lowercase and include no special characters
            this_spp = re.findall('^[xX\u00D7]?[a-z-]*$', this_taxon[1])
            if (len(this_spp) == 1):
                this_spp = this_spp[0]
                spp_idx = [1]
                if verbose:
                    print("\tspecies epithet in second position:", this_spp)
            else:
                print(this_taxon_init)
                raise ValueError("specific epithet not identified")
    # first check for cultivars and cultivar groups, which can cause errors recognizing other parts
    # note, a wide range of characters are allowed in cultivar names
    cvgp_idx_init = return_matches(arr = this_taxon, regexp = '^Group$|^Gp$', mode = 'index')
    cvgp_idx_init = [x for x in cvgp_idx_init if bool(re.search('^[A-Z][A-Za-z-]*$', this_taxon[x - 1]))]
    if (len(cvgp_idx_init) == 1):
        # assume group name is the word preceding group
        cv_idx.append(cvgp_idx_init[0] - 1)
        cv_idx.append(cvgp_idx_init[0])
        this_cv = [e for i, e in enumerate(this_taxon) if i in cv_idx]
        this_cv = " ".join(this_cv)
        if verbose:
            print("\ttaxon includes cultivar group:", this_cv)
    elif (len(cvgp_idx_init) > 1):
        print(this_taxon_init)
        raise ValueError("multiple cultivar groups detected")
    cv_idx_init = return_matches(arr = this_taxon, regexp = "^cv\.?$|^cv\.gp\.?$|^[\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4'\"][A-Za-z-0-9À-ÿĀ-ſƀ-ȳ.'\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4]{2,}[\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4'\"]$|^[\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4'\"][A-Za-z-0-9À-ÿĀ-ſƀ-ȳ.'\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4]+$", mode = 'index')
    # account for dutch author 't Hart
    cv_idx_init = [x for x in cv_idx_init if this_taxon[x] != "'t"]
    if (len(cv_idx_init) > 0):
        if (len(cv_idx_init) > 1):
            print("multiple cultivars or cultivar groups detected\n\t", this_taxon_init)
        if (not authority):
            # assume cultivar name is the entire part of name following first cv., cv.gp., or 'Cultivar Name' word
            # to the end of the name, only do this if authorities are not included!
            cv_idx = [i for i, x in enumerate(this_taxon) if i >= cv_idx_init[0]]
        else:
            # if authorities included, only split cv if 'Cultivar Name' method is used
            cv_idx_init2 = return_matches(arr = this_taxon, regexp = "^[\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4'\"][A-Za-z-0-9À-ÿĀ-ſƀ-ȳ.'\u201B\u2019\u2018\u2032\u2033\u275C\u00B4]{2,}[\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4'\"]$|^[\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4'\"][A-Za-z-0-9À-ÿĀ-ſƀ-ȳ.'\u201B\u2019\u2018\u2032\u2033\u275C\u00B4]+$|^[A-Za-z-0-9À-ÿĀ-ſƀ-ȳ.'\u201B\u2019\u2018\u275C\u00B4]+[\u201A\u201B\u201E\u201C\u201F\u201D\u2019\u2018\u2032\u2033\u0022\u275D\u275E\u2E42\u301D\u301E\u301F\uFF02\u275B\u275C\u275F\u00B4'\"]$", mode = 'index')
            if (len(cv_idx_init2) > 0):
                cv_idx = [i for i, x in enumerate(this_taxon) if i >= min(cv_idx_init2) and i <= max(cv_idx_init2)]
        this_cv = [e for i, e in enumerate(this_taxon) if i in cv_idx]
        if (len(cv_idx) > 1):
            this_cv = " ".join(this_cv)
        elif (len(cv_idx) == 1):
            this_cv = this_cv[0]
        elif (len(cv_idx) == 0):
            this_cv = ""
        if verbose:
            print("\ttaxon includes cultivar:", this_cv)
    # check for infraspecific ranks and names
    #
    # ICBN 5.24.1: the name of an infraspecific taxon is a combination of the name of a species and an infraspecific epithet. Saxifraga aizoon subf. surculosa Engl. & Irmsch may also be referred to as Saxifraga aizoon var. aizoon subvar. brevifolia f. multicaulis subf. surculosa Engl. & Irmsch.
    # keep only finest resolution infraspecific part
    # assumes whoever wrote the names didn't put the sections in the wrong order!
    #
    # check for case where identified rank is at the end of the name
    # could be unnamed variant or, for example, authority names like L. f.
    infra_rank_idx = return_matches(arr = this_taxon, regexp = '^(notho)?subsp\.?$|^(notho)?ssp\.?$|^(notho)?(con|sub)?var\.?$|^(notho)?(con|sub)?v\.?$|^(notho)?(sub)?f\.?$|^(notho)?(sub)?f[ao]?\.?$|^(notho)?(sub)?forma\.?$|^(notho)?(sub)?fma\.?$|^[\u03B1\u03B2\u03B3\u03B4\u03B5\u03B6\u03B7\u03B8\u03BB\u03C1]$', mode = 'index')
    infra_rank_idx = [x for x in infra_rank_idx if x != len(this_taxon) - 1]
    # special case where f. as in filius is followed by infraspecific rank
    # note: uppercase added recently to allow capitalized infraspecific names
    # note 2: hybrid symbol added recently to allow infraspecies hybrids
    infra_rank_idx = [x for x in infra_rank_idx if bool(re.search('^(?!^ex$)[A-Za-z-\u00D7]*$', this_taxon[x + 1]))]
    infra_idx = [i + 1 for i in infra_rank_idx]
    if (len(infra_rank_idx) > 0):
        this_infra_rank = [e for i, e in enumerate(this_taxon) if i in infra_rank_idx]
        this_infra = [e for i, e in enumerate(this_taxon) if i in infra_idx]
        if (len(infra_rank_idx) > 1):
            print(this_taxon_init, "\n\tmultiple infraspecific ranks detected (", ", ".join(this_infra_rank), ") keeping finest level (", this_infra_rank[-1],")", sep = '')
        # subset to finest rank
        this_infra_rank = this_infra_rank[-1]
        this_infra = this_infra[-1]
        if verbose:
            print("\tinfraspecific rank (", this_infra_rank, ") and name (", this_infra, ") detected", sep = '')
        # check for species authority before any infraspecific name and extract
        all_idx = [0] + spp_idx + infra_rank_idx + infra_idx + cv_idx
        temp_auth_idx = [i for i in range(len(this_taxon[:infra_idx[0]])) if i not in all_idx]
        if (len(temp_auth_idx) > 0):
            temp_auth_init = [e for i, e in enumerate(this_taxon) if i in temp_auth_idx]
            temp_auth = " ".join(temp_auth_init)
            if verbose:
                print("\tauthority before infraspecific name:", temp_auth)
    elif (len(this_taxon) > (spp_idx[-1] + 1)):
        # if rank is not present, check whether the element after the species name is a lowercase word with no special characters
        # if so, assign as infraspecific name
        # this is prone to error! check manually!
        poten_infra = re.findall("^[a-z-]{3,}$", this_taxon[spp_idx[-1] + 1])
        if (len(poten_infra) > 0):
            print(this_taxon_init, "\n\tcheck manually! potential infraspecific name detected with no rank:", poten_infra[0])
            infra_idx = [spp_idx[-1] + 1]
            this_infra = poten_infra[0]
    # test whether taxon is a nothosubspecies, nothovar, etc.
    if len(infra_rank_idx) > 0:
        infra_hyb_test1 = return_matches(arr = this_taxon[infra_rank_idx[0]:], regexp = '^[xX\u00D7]$', mode = "index")
        infra_hyb_test2 = return_matches(arr = this_taxon[infra_rank_idx[0]:], regexp = '^[\u00D7]', mode = "bool")
    else:
        infra_hyb_test1 = []
        infra_hyb_test2 = [False]
    infra_hyb_test = (len(infra_hyb_test1) > 0) | any(infra_hyb_test2)
    if (len(infra_hyb_test1) > 0):
        infra_hyb_test1 = [x + infra_rank_idx[0] for x in infra_hyb_test1]
        # taxon is a hybrid
        if (len(infra_hyb_test1) == 1):
            if (infra_hyb_test1[0] == 0):
                # check whether there is just a space before standard hybrid name
                infra_idx.append(infra_hyb_test1[0] + 1)
                infra_idx.append(infra_hyb_test1[0] + 2)
            else:
                # or whether hybrid name shows parent epithets
                infra_idx.append(infra_hyb_test1[0])
                infra_idx.append(infra_hyb_test1[0] + 1)
                infra_idx.append(infra_hyb_test1[0] + 2)
        elif (len(infra_hyb_test1) == 2):
            # species has multiple parent epithets
            infra_idx.append(infra_hyb_test1[0])
            infra_idx.append(infra_hyb_test1[0] + 1)
            infra_idx.append(infra_hyb_test1[1])
            infra_idx.append(infra_hyb_test1[1] + 1)
            infra_idx.append(infra_hyb_test1[1] + 2)
        else:
            # other hybrid configurations not allowed at the moment
            print(this_taxon_init)
            raise ValueError("hybrid species epithet non-parsable")
        this_infra_init = [e for i, e in enumerate(this_taxon) if i in infra_idx]
        this_infra = ' '.join(this_infra_init)
    #
    # if authors are included: 
    # authority is assigned to be the text at the end of the name (but before cultivar)
    # ignore any authories in middle levels (e.g., subsp. author for var.)
    # if authors are not included, assign to infra
    name_idx = [0] + spp_idx + temp_auth_idx + infra_rank_idx + infra_idx
    if len(cv_idx) > 0:
        auth_init = this_taxon[max(name_idx) + 1:cv_idx[0]]
    else:
        auth_init = this_taxon[max(name_idx) + 1:]
    if len(auth_init) > 0:
        if temp_auth != '':
            this_auth = temp_auth
            this_infra_auth = " ".join(auth_init)
            if verbose:
                print("\tspecies authority (", this_auth, ") and infraspecies authority (", this_infra_auth,") detected", sep = '')
        elif this_infra != '':
            this_infra_auth = " ".join(auth_init)
            if verbose:
                print("\tonly infraspecies authority (", this_infra_auth,") detected", sep = '')
        else:
            this_auth = " ".join(auth_init)
            if verbose:
                print("\tonly species authority (", this_auth, ") detected", sep = '')
    elif temp_auth != '':
        this_auth = temp_auth
        if verbose:
                print("\tonly species authority (", this_auth, ") detected", sep = '')
    if ((not authority) & (this_auth != '')):
        if verbose:
            print(this_taxon_init, "\n\tauthority detected (", this_auth, "), but auth_split_bool = False\t...reassigning to infraspecific name", sep = '')
        this_infra = this_infra + this_auth
        this_auth = ''
    if this_auth != "":
        test_cap = str.islower(this_auth)
        if test_cap:
            print(this_taxon_init, "\n\tauthority is lowercase!\n\t\t", this_auth)
    output = [this_gen, this_spp, this_infra_rank, this_infra, this_auth, this_infra_auth, this_cv, gen_hyb_test, spp_hyb_test]
    if verbose:
        print("\n", this_taxon_init, "\n\ngenus\t", this_gen, "\nspecies\t", this_spp, "\nauthority\t", this_auth, "\ninfra rank\t", this_infra_rank, "\ninfra\t", this_infra, "\ninfra auth\t", this_infra_auth, "\ncultivar\t", this_cv, "\ngenus hybrid?", gen_hyb_test, "\tspp hybrid?", spp_hyb_test)
    return(output)
#
# typographical ligatures and nonstandard punctuation marks to replace
special_chars = {"\u0132": "IJ",
                "\u0133": "ij",
                "\uFB00": "ff",
                "\uFB03": "ffi",
                "\uFB04": "ffl",
                "\uFB01": "fi",
                "\uFB02": "fl",
                "\u201E": '"',
                "\u201C": '"',
                "\u201F": '"',
                "\u201D": '"',
                "\u2019": "'",
                "\u2018": "'",
                "\u2032": "'",
                "\u2033": '"',
                "\u0022": '"',
                "\u275D": '"',
                "\u275E": '"',
                "\u2E42": '"',
                "\u301D": '"',
                "\u301E": '"',
                "\u301F": '"',
                "\uFF02": '"',
                "\u201A": "'",
                "\u201B": "'",
                "\u275B": "'",
                "\u275C": "'",
                "\u275F": "'",
                "\u00B4": "'",
                "\u0060": "'"}
# ICN 60.6.: diacritical signs are not used in Latin plant names
accent_chars = {"\u00C0": "A",
                "\u00C1": "A",
                "\u00C2": "A",
                "\u00C3": "A",
                "\u00C4": "Ae",
                "\u00C5": "Ao",
                "\u0152": "OE",
                "\u0153": "oe",
                "\u00C6": "AE",
                "\u00E6": "ae",
                "\u00C7": "C",
                "\u00C8": "E",
                "\u00C9": "E",
                "\u00CA": "E",
                "\u00CB": "E",
                "\u00CC": "I",
                "\u00CD": "I",
                "\u00CE": "I",
                "\u00CF": "I",
                "\u00D2": "O",
                "\u00D3": "O",
                "\u00D4": "O",
                "\u00D5": "O",
                "\u00D6": "Oe",
                "\u00D8": "Oe",
                "\u00D9": "U",
                "\u00DA": "U",
                "\u00DB": "U",
                "\u00DC": "Ue",
                "\u00DF": "ss",
                "\u00E0": "a",
                "\u00E1": "a",
                "\u00E2": "a",
                "\u00E3": "a",
                "\u00E4": "ae",
                "\u00E5": "ao",
                "\u00E7": "c",
                "\u00E8": "e",
                "\u00E9": "e",
                "\u00EA": "e",
                "\u00EB": "e",
                "\u00EC": "i",
                "\u00ED": "i",
                "\u00EE": "i",
                "\u00EF": "i",
                "\u00F0": "o",
                "\u00F1": "n",
                "\u00F2": "o",
                "\u00F3": "o",
                "\u00F4": "o",
                "\u00F5": "o",
                "\u00F6": "oe",
                "\u00F8": "o",
                "\u00F9": "u",
                "\u00FA": "u",
                "\u00FB": "u",
                "\u00FC": "ue"}
spacing_chars = {"\u00A0": "\u0020",
                "\u1680": "\u0020",
                "\u2000": "\u0020",
                "\u2001": "\u0020",
                "\u2002": "\u0020",
                "\u2003": "\u0020",
                "\u2004": "\u0020",
                "\u2005": "\u0020",
                "\u2006": "\u0020",
                "\u2007": "\u0020",
                "\u2008": "\u0020",
                "\u2009": "\u0020",
                "\u200A": "\u0020",
                "\u202F": "\u0020",
                "\u205F": "\u0020",
                "\u3000": "\u0020"}
escape_chars = {u'\u000A': u'',
                u'\u000B': u'',
                u'\u000C': u'',
                u'\u000D': u'',
                u'\u0085': u'',
                u'\u0090': u'',
                u'\u2028': u'',
                u'\u2029': u'',
                u'\u003F': u'',
                u'\uFEFF': u''}
#
# load main data file
taxa = pd.read_csv(full_filename, delimiter = ',', encoding = 'utf-8', keep_default_na = False)
print("\nloaded file =", full_filename)
n_row_init = len(taxa)
print("\t", n_row_init, "rows (taxa) in file")
#
# check available columns
full_col_names = ['Genus', 'SpecificEpithet', 'Authority', 'InfraspecificRank', 'InfraspecificName', 'InfraspecificAuthority', 'Cultivar']
full_col_names2 = full_col_names + ['TaxonName']
this_col_names = list(taxa.columns.values)
# taxon name parts included in current file
col_names = [x for x in this_col_names if x in full_col_names]
# all other non-taxon-name columns included in current file
other_cols = [x for x in this_col_names if x not in full_col_names2]
# check for columns that are suspected to be errors
err_col_names = np.array(['Genera', 'Specific Epithet', 'SpecificEpithets', 'Specific Epithets', 'Infraspecific Rank', 'Infraspecific Name', 'IntraspecificRank', 'IntraspecificName', 'Intraspecific Rank', 'Intraspecific Name', 'Specific Epithet', 'Authorities', 'Author', 'Authors', 'Infraspecific Rank', 'InfraspecificRanks', 'Infraspecific Ranks', 'InfraspecificNames','Infraspecific Name', 'Infraspecific Authority', 'InfraspecificAuthorities', 'Infraspecific Authorities', 'InfraspecificAuthor', 'InfraspecificAuthors', 'Cultivars', 'CultivarName', 'CultivarNames'])
check_col_names  = np.array(this_col_names)
check_err_cols = np.isin(err_col_names, check_col_names)
# quit program if error column detected
if any(check_err_cols):
    print(", ".join(err_col_names[check_err_cols]))
    sys.exit("ERROR: suspected error in column names, names must be exact: Genus, SpecificEpithet, Authority, InfraspecificRank, InfraspecificName, InfraspecificAuthority, Cultivar (or TaxonName)")
# check whether TaxonName column and Genus, SpecificEpithet column are included and quit
if all(x in this_col_names for x in ['Genus', 'SpecificEpithet']) and ("TaxonName" in this_col_names):
    sys.exit("ERROR: Genus and SpecificEpithet column detected as well as TaxonName column. choose one option and re-run: TaxonName only or split into parts only")
#
# split into genus, spp, infra columns if necessary
if not all(x in this_col_names for x in ['Genus', 'SpecificEpithet']):
    print("\ntaxon names not split into genus, specific epithet, and infraspecific name")
    if "TaxonName" in this_col_names:
        print("\n'TaxonName' column detected, splitting names\n*this is prone to errors! check output carefully!*")
        genus = []
        species = []
        infra_rank = []
        infra = []
        author = []
        infra_auth = []
        cultivar = []
        for y in range(len(taxa)):
            this_taxon_init = taxa.loc[y, 'TaxonName']
            output = split_taxon(this_taxon_init)
            genus.append(output[0])
            species.append(output[1])
            infra_rank.append(output[2])
            infra.append(output[3])
            author.append(output[4])
            infra_auth.append(output[5])
            cultivar.append(output[6])
        # genus and specific epithet must be stored
        taxa['Genus'] = genus
        taxa['SpecificEpithet'] = species
        # store others only if missing
        if auth_split_bool and "Authority" not in this_col_names:
            taxa['Authority'] = author
        elif auth_split_bool:
            taxa['Authority'] = taxa['Authority'] + author
        if "InfraspecificRank" not in this_col_names:
            taxa['InfraspecificRank'] = infra_rank
        if "InfraspecificName" not in this_col_names:
            taxa['InfraspecificName'] = infra
        if auth_split_bool and "InfraspecificAuthority" not in this_col_names:
            taxa['InfraspecificAuthority'] = infra_auth
        elif auth_split_bool:
            taxa['InfraspecificAuthority'] = taxa['InfraspecificAuthority'] + infra_auth
        if "Cultivar" not in this_col_names:
            taxa['Cultivar'] = cultivar
        this_col_names = list(taxa.columns.values)
        col_names = [x for x in this_col_names if x in full_col_names]
    else:
        sys.exit("ERROR: 'TaxonName' column not identified. check this and re-run")
else:
    print("\nchecking for infraspecific ranks/names in specific epithet column")
    # check whether specific epithet column includes infraspecific rank and name
    infra_regexp = r' (var|v|convar|conv|subvar|subv|subsp|ssp|f[ao]?|fma|forma|subf[oa]?|subfma)[\. ]'
    check_infra = return_matches(arr = taxa['SpecificEpithet'], regexp = infra_regexp)
    if len(check_infra) > 0:
        # if so, place infraspecific name and rank in own columns and remove from specific epithet
        infra_samp = list(set(check_infra))
        if len(infra_samp) > subset_len:
            infra_samp = random.sample(infra_samp, subset_len)
        print("\ninfraspecific ranks/names are present in specific epithet column\t...moving to own column (", len(check_infra), " matches)\n\t\t", ', '.join(infra_samp), sep = '')
        infra_idx = return_matches(arr = taxa['SpecificEpithet'], mode = 'index', regexp = infra_regexp)
        # edit initial epithet to break into three parts
        infedit_regexp = '(.*) (var\.?|v\.|subvar\.?|subv\.?|convar\.?|conv\.?|subsp\.?|ssp\.?|f[ao]?\.|fma\.?|forma |subf[oa]?\.?|subfma\.?) ?(.*)'
        infra_init = check_infra
        infra_name = [re.sub(infedit_regexp, '\\3', v) for v in infra_init]
        infra_rank = [re.sub(infedit_regexp, '\\2', v) for v in infra_init]
        spp_epi = [re.sub(infedit_regexp, '\\1', v) for v in infra_init]
        # add infraspecific rank and name columns if missing
        if 'InfraspecificRank' not in col_names:
            taxa['InfraspecificRank'] = [str() for c in 'c' * taxa.shape[0]]
            col_names.append('InfraspecificRank')
        if 'InfraspecificName' not in col_names:
            taxa['InfraspecificName'] = [str() for c in 'c' * taxa.shape[0]]
            col_names.append('InfraspecificName')
        # if infraspecific name and rank columns already existed, the text will be combined (this case would be unusual)
        taxa.loc[infra_idx, 'InfraspecificName'] = taxa.loc[infra_idx, 'InfraspecificName'] + infra_name
        taxa.loc[infra_idx, 'InfraspecificRank'] = taxa.loc[infra_idx, 'InfraspecificRank'] + infra_rank
        # specific epithet, however, needs to be written over rather than concatenated
        taxa.loc[infra_idx, 'SpecificEpithet'] = spp_epi
    print("\nchecking for cultivar names in specific epithet and infraspecific name")
    col_temp = ['SpecificEpithet', 'InfraspecificName']
    cols = [x for x in col_temp if x in col_names]
    for col in cols:
        print("\tchecking column:", col)
        # check whether specific epithet includes cultivar
        cv_regexp = "(?:^| )(?:cvs?\. ?['\"]?|cvs?\.? ['\"]?|cvgr\.? ?['\"]?|['\"])([,\w .0-9-]*)(?:['\"]| |$)"
        cv_init = taxa[col]
        # remove apostrophes before checking
        cv_init = [re.sub("(\w)'(\w)", r'\1\2', v) for v in cv_init]
        check_cv = return_matches(arr = cv_init, regexp = cv_regexp)   
        if (len(check_cv) > 0):
            cv_samp = list(set(check_cv))
            if len(cv_samp) > subset_len:
                cv_samp = random.sample(cv_samp, subset_len)
            print("\tcultivars identified, moving to cultivar column (", len(cv_samp), " matches)\n\t\t", ', '.join(cv_samp), sep = '')
            cv_idx = return_matches(arr = cv_init, regexp = cv_regexp, mode = 'index')
            # edit initial epithet to break into two parts
            cvedit_regexp = "(.*)(?:^| )((?:cvs?\. ?['\"]?|cvs?\.? ['\"]?|cvg[rp]\.? ?['\"]?|['\"])(?:[,\w .0-9-]*)(?:['\"]| |$))"
            cvedit_init = check_cv
            cv_edit = [re.sub(cvedit_regexp, '\\2', v) for v in cvedit_init]
            name_edit_init = [re.sub(cvedit_regexp, '\\1', v) for v in cvedit_init]
            cleancv_regexp = " ?cvs?\.?( |$)|cvg[rp]\.?( |$)"
            name_edit = [re.sub(cleancv_regexp, '', v) for v in name_edit_init]
            # add cultivar column if missing
            if 'Cultivar' not in col_names:
                taxa['Cultivar'] = [str() for c in 'c' * taxa.shape[0]]
                col_names.append('Cultivar')
            # if cultivar column already existed, the text will be combined
            taxa.loc[cv_idx, 'Cultivar'] = taxa.loc[cv_idx, 'Cultivar'] + cv_edit
            # specific epithet, however, needs to be written over rather than concatenated
            taxa.loc[cv_idx, col] = name_edit
    print("\ncolumn names identified\n\t", ", ".join(col_names), "\ncreating 'TaxonName' column (to become 'OriginalName' column)")
    if 'InfraspecificName' in col_names:
        if 'InfraspecificRank' in col_names:
            if 'Authority' in col_names:
                if 'InfraspecificAuthority' in col_names:
                    if 'Cultivar' in col_names:
                        taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'Authority', 'InfraspecificRank', 'InfraspecificName', 'InfraspecificAuthority', 'Cultivar']].agg(' '.join, axis = 1)
                        print("\tgenus, specific epithet, author, infraspecific rank, infraspecific name, infraspecific author, cultivar")
                    else:
                        taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'Authority', 'InfraspecificRank', 'InfraspecificName', 'InfraspecificAuthority']].agg(' '.join, axis = 1)
                        print("\tgenus, specific epithet, author, infraspecific rank, infraspecific name, infraspecific author")
                else:
                    if 'Cultivar' in col_names:
                        taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'Authority', 'Cultivar']].agg(' '.join, axis = 1)
                        print("\tgenus, specific epithet, infraspecific rank, infraspecific name, author, cultivar")
                    else:
                        taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'Authority']].agg(' '.join, axis = 1)
                        print("\tgenus, specific epithet, infraspecific rank, infraspecific name, author")
            else:
                if 'Cultivar' in col_names:
                    taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'Cultivar']].agg(' '.join, axis = 1)
                    print("\tgenus, specific epithet, infraspecific rank, infraspecific name, cultivar")
                else:
                    taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName']].agg(' '.join, axis = 1)
                    print("\tgenus, specific epithet, infraspecific rank, infraspecific name")
        else:
            if 'Authority' in col_names:
                if 'InfraspecificAuthority' in col_names:
                    if 'Cultivar' in col_names:
                        taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'Authority', 'InfraspecificName', 'InfraspecificAuthority', 'Cultivar']].agg(' '.join, axis = 1)
                        print("\tgenus, specific epithet, author, infraspecific name, infraspecific author, cultivar")
                    else:
                        taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'Authority', 'InfraspecificName', 'InfraspecificAuthority']].agg(' '.join, axis = 1)
                        print("\tgenus, specific epithet, author, infraspecific name, infraspecific author")
                else:
                    if 'Cultivar' in col_names:
                        taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificName', 'Authority', 'Cultivar']].agg(' '.join, axis = 1)
                        print("\tgenus, specific epithet, infraspecific name, author, cultivar")
                    else:
                        taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificName', 'Authority']].agg(' '.join, axis = 1)
                        print("\tgenus, specific epithet, infraspecific name, author")
            else:
                if 'Cultivar' in col_names:
                    taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificName', 'Cultivar']].agg(' '.join, axis = 1)
                    print("\tgenus, specific epithet, infraspecific name, cultivar")
                else:
                    taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificName']].agg(' '.join, axis = 1)
                    print("\tgenus, specific epithet, infraspecific name")
    else:
        if 'Authority' in col_names:
            if 'Cultivar' in col_names:
                print('NOTE: cultivar column identified, but no InfraspecificName column, check manually!')
                taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'Authority', 'Cultivar']].agg(' '.join, axis = 1)
                print("\tgenus, specific epithet, author, cultivar")
            else:
                taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'Authority']].agg(' '.join, axis = 1)
                print("\tgenus, specific epithet, author")
        else:
            if 'Cultivar' in col_names:
                print('NOTE: cultivar column identified, but no InfraspecificName column, check manually!')
                taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'Cultivar']].agg(' '.join, axis = 1)
                print("\tgenus, specific epithet, cultivar")
            else:
                taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet']].agg(' '.join, axis = 1)
                print("\tgenus, specific epithet")
    taxa['TaxonName'] = [re.sub(' +', ' ', v) for v in taxa['TaxonName']]
    taxa['TaxonName'] = [re.sub(' $', '', v) for v in taxa['TaxonName']]
if ((full_col_names[0] not in this_col_names) & (full_col_names[1] not in this_col_names)):
    sys.exit("ERROR: file does not meet minimum requirements, Genus and/or SpecificEpithet columns are missing")
#
# remove rows not defined to species level or below
print("\nchecking for taxa not defined to species level or below")
taxa = taxa[taxa['Genus'].str.strip().astype(bool)]
n_row_gen = len(taxa)
if (n_row_init - n_row_gen > 0):
    print("\tremoved", n_row_init - n_row_gen, "rows missing Genus")
taxa = taxa[taxa['SpecificEpithet'].str.strip().astype(bool)]
n_row_sp = len(taxa)
if (n_row_gen - n_row_sp > 0):
    print("\tremoved", n_row_gen - n_row_sp, "rows missing SpecificEpithet")
# check for taxa above species level, sect., agg., subgen., etc.
asl_regexp = r'\b[Ss]ect\.?\b|\b[Ss]ection\b|\b[Aa]gg\.?\b|\b[Aa]ggregate\b|\b[Ss]ubge?n?\.?\b|\b[Ss]ubgenus\b'
check_asl = return_matches(arr = taxa['SpecificEpithet'], regexp = asl_regexp)
if len(check_asl) > 0:
    print("\t", len(check_asl), "taxa with undefined specific epithets\t...removing\n\t\t", ', '.join(set(check_asl)))
    asl_bool = return_matches(arr = taxa['SpecificEpithet'], mode = 'bool', regexp = asl_regexp)
    taxa = taxa[~np.array(asl_bool)]
# check for undefined specific epithets, spp, sp, etc.
spp_regexp = '^ *[Ss]pecies *$|^ *[Ss]pp?\.? *$|^ *[Ss]pec\.? *$|^ *\?* *$'
check_spp = return_matches(arr = taxa['SpecificEpithet'], regexp = spp_regexp)
if len(check_spp) > 0:
    print("\t", len(check_spp), "taxa with undefined specific epithets\t...removing\n\t\t", ', '.join(set(check_spp)))
    sp_missing_bool = return_matches(arr = taxa['SpecificEpithet'], mode = 'bool', regexp = spp_regexp)
    taxa = taxa[~np.array(sp_missing_bool)]
# check for hybrids not assigned to species level
check_hsp = return_matches(arr = taxa['SpecificEpithet'], regexp = r'^ *[x\u00D7] *$|\b[Hh]ybrids?\b')
if len(check_hsp) > 0:
    print("\t", len(check_hsp), "taxa are hybrids not defined to species level\t...removing\n\t\t", ', '.join(set(check_hsp)))
    hyb_nospp_bool = return_matches(arr = taxa['SpecificEpithet'], mode = 'bool', regexp = r'^ *[x\u00D7] *$|\b[Hh]ybrids?\b')
    taxa = taxa[~np.array(hyb_nospp_bool)]
# check for cultivars not assigned to species
check_csp = return_matches(arr = taxa['SpecificEpithet'], regexp = '^ *cv\.? *$|^ *cultivars? *$')
if len(check_csp) > 0:
    print("\t", len(check_csp), "taxa are cultivars not defined to species level\t...removing\n\t\t", ', '.join(set(check_csp)))
    hyb_nospp_bool = return_matches(arr = taxa['SpecificEpithet'], mode = 'bool', regexp = '^ *cv\.? *$|^ *cultivars? *$')
    taxa = taxa[~np.array(hyb_nospp_bool)]
# taxa names involving open nomenclature not acceptable (cannot identify unique species)
# open names include: sp., spp. (see above), sp. nov., sp. aff., aff., cf., hybrid(s)
opnom_regexp = r'(^| )([Ss]pecies|[Ss]pp?\.? ?[Aa]ff\.?|[Ss]pp?\.? ?[Nn]ov[a.]?|[Ss]pec\.? ?[Nn]ov[a.]?|[Ss]pecies [Nn]ov[a.]|[Aa]ff\.?|[Cc]f\.?)(?:$| )'
for x in range(len(col_names)):
    check_open = return_matches(arr = taxa[col_names[x]], regexp = opnom_regexp)
    if len(check_open) > 0:
        print(col_names[x], "\t", len(check_open), "taxa involving open nomenclature detected\t...removing\n\t\t", ', '.join(set(check_open)))
        open_bool = return_matches(arr = taxa[col_names[x]], mode = 'bool', regexp = opnom_regexp)
        taxa = taxa[~np.array(open_bool)]
taxa = taxa.reset_index()
#
# remove whitespace in all columns of interest
for x in range(len(col_names)):
    for index, row in taxa.iterrows():
        row[col_names[x]] = row[col_names[x]].strip()
print("\nremoved whitespace in all columns")
#
# genus names to sentence case and specific epithets to lowercase
taxa['Genus'] = [x.title() if x else '' for x in taxa['Genus']]
taxa['SpecificEpithet'] = [x.lower() if x else '' for x in taxa['SpecificEpithet']]
if 'InfraspecificName' in col_names:
    taxa['InfraspecificName'] = [x.lower() if x else '' for x in taxa['InfraspecificName']]
print("\ngenus name to sentence case, specific and infraspecific epithets to lowercase")
#
# replace all ligature and alternative quotation characters
print("\nchecking for ligatures and alternative quotation characters")
for x in range(len(col_names)):
    print("\tchecking column:", col_names[x])
    taxa[col_names[x]] = replace_dict_substr(arr = taxa[col_names[x]], dic = special_chars)
#
# change all whitespace characters to unicode space
# remove all linebreak characters
# remove question marks
print("\nchecking for alternative whitespace characters, escape characters, or question marks")
for x in range(len(col_names)):
    print("\tchecking column:", col_names[x])
    taxa[col_names[x]] = replace_dict_substr(arr = taxa[col_names[x]], dic = spacing_chars)
    taxa[col_names[x]] = replace_dict_substr(arr = taxa[col_names[x]], dic = escape_chars)
#
# also remove accented characters for genus, specific epithet, and infraspecific name only
print("\nchecking for accented characters (genus, specific epithet, and infraspecific name only)")
col_temp = ['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName']
cols = [x for x in col_temp if x in col_names]
for col in cols:
    print("\tchecking column:", col)
    taxa[col] = replace_dict_substr(arr = taxa[col], dic = accent_chars)
    test = ', '.join(taxa[col])
    non_ascii = re.findall('[^\x00-\x7F]', test)
    if len(non_ascii) > 0:
        print("\tNOTE: non-ascii characters identified, check manually!\n\t\t", ', '.join(set(non_ascii)))
    else:
        print("\tall remaining characters are ascii")
#
# check for numbers in all parts of names
# these are very unusual and would only be allowed in cultivar names
print("\nchecking for numbers")
for x in range(len(col_names)):
    check_num = return_matches(arr = taxa[col_names[x]], regexp = '[0-9]')
    if len(check_num) > 0:
        if col_names[x] in ['Genus', 'SpecificEpithet', 'Authority', 'InfraspecificAuthority']:
            print(col_names[x], "\n", check_num)
            sys.exit("ERROR: numbers identified in Genus, SpecificEpithet, Authority, or InfraspecificAuthority columns")
        num_samp = list(set(check_num))
        if len(num_samp) > subset_len:
            num_samp = random.sample(num_samp, subset_len)
        print(col_names[x], "\tNOTE: numbers identified, check manually! (", len(check_num)," matches) \n\t\t", ', '.join(num_samp), sep = '')
#
# check for unusual accented characters for authority and cultivar, to be checked manually
print("\nchecking for unusual special characters (authority and cultivar)")
col_temp = ['Authority', 'InfraspecificAuthority', 'Cultivar']
cols = [x for x in col_temp if x in col_names]
for col in cols:
    print("\tchecking column:", col)
    test = ', '.join(taxa[col])
    spec_chars = re.findall("[^A-Za-z-0-9À-ÿĀ-ſƀ-ȳ.,' \n\(\)&]", test)
    if len(spec_chars) > 0:
        print("\tNOTE: unusual special characters identified, check manually!\n\t\t", ', '.join(set(spec_chars)))
    else:
        print("\tno unusual special characters")
#
# clean up extra spaces
print("\nchecking for extra spaces")
for x in range(len(col_names)):
    check_multispace = return_matches(arr = taxa[col_names[x]], regexp = ' {2,}')
    if len(check_multispace) > 0:
        multispace_samp = list(set(check_multispace))
        if len(multispace_samp) > subset_len:
            multispace_samp = random.sample(multispace_samp, subset_len)
        print(col_names[x], "\textra spaces identified\t...standardizing (", len(check_multispace), " matches)\n\t\t", ', '.join(multispace_samp), sep = '')
        taxa[col_names[x]] = [re.sub(' {2,}', r' ', v) for v in taxa[col_names[x]]]
for x in range(len(col_names)):
    check_space = return_matches(arr = taxa[col_names[x]], regexp = '^ | $')
    if len(check_space) > 0:
        space_samp = list(set(check_space))
        if len(space_samp) > subset_len:
            space_samp = random.sample(space_samp, subset_len)
        print(col_names[x], "\tbeginning/ending spaces identified\t...removing (", len(check_space), " matches)\n\t\t", ', '.join(space_samp), sep = '')
        taxa[col_names[x]] = [v.strip() for v in taxa[col_names[x]]]
#
# check for any kind of nonstandard characters in genus, specific epithet, and infraspecific name (anything outside of ASCII alphabet, hyphens, hybrid symbol, space)
print("\nchecking for nonstandard characters in genus, specific epithet")
col_temp = ['Genus', 'SpecificEpithet', 'InfraspecificName']
cols = [x for x in col_temp if x in col_names]
for col in cols:
    print("\tchecking column:", col)
    check_unacc = return_matches(arr = taxa[col], regexp = '[^A-Za-z- \u00d7.]')
    if len(check_unacc) > 0:
        unacc_samp = list(set(check_unacc))
        if len(unacc_samp) > subset_len:
            unacc_samp = random.sample(unacc_samp, subset_len)
        print(col, "\tNOTE: nonstandard characters identified\t...check manually! (", len(check_unacc)," matches) \n\t\t", ', '.join(unacc_samp), sep = '')
#
# clean up spacing around parentheses
# remove synonyms in name
print("\nchecking for parentheses spacing errors and synonyms in name")
for x in range(len(col_names)):
    check_parens1 = return_matches(arr = taxa[col_names[x]], regexp = '([a-zA-Z.]+)\(')
    if len(check_parens1) > 0:
        parens1_samp = list(set(check_parens1))
        if len(parens1_samp) > subset_len:
            parens1_samp = random.sample(parens1_samp, subset_len)
        print(col_names[x], "\tparentheses spacing error identified\t...standardizing (", len(check_parens1), " matches)\n\t\t", ', '.join(parens1_samp), sep = '')
        taxa[col_names[x]] = [re.sub('([a-zA-Z.]+)\(', r'\1 (', v) for v in taxa[col_names[x]]]
    check_parens2 = return_matches(arr = taxa[col_names[x]], regexp = '\)([a-zA-Z.]+)')
    if len(check_parens2) > 0:
        parens2_samp = list(set(check_parens2))
        if len(parens2_samp) > subset_len:
            parens2_samp = random.sample(parens2_samp, subset_len)
        print(col_names[x], "\tparentheses spacing error identified\t...standardizing (", len(check_parens2), " matches)\n\t\t", ', '.join(parens2_samp), sep = '')
        taxa[col_names[x]] = [re.sub('\)([a-zA-Z.]+)', r') \1', v) for v in taxa[col_names[x]]]
    check_parens3 = return_matches(arr = taxa[col_names[x]], regexp = ' \)')
    if len(check_parens3) > 0:
        parens3_samp = list(set(check_parens3))
        if len(parens3_samp) > subset_len:
            parens3_samp = random.sample(parens3_samp, subset_len)
        print(col_names[x], "\tparentheses spacing error identified\t...standardizing (", len(check_parens3), " matches)\n\t\t", ', '.join(parens3_samp), sep = '')
        taxa[col_names[x]] = [re.sub(' \)', r')', v) for v in taxa[col_names[x]]]
    check_syn = return_matches(arr = taxa[col_names[x]], regexp = '([A-Z][^A-Z\(]+)(.*) \(=.*\)')
    if len(check_syn) > 0:
        syn_samp = list(set(check_syn))
        if len(syn_samp) > subset_len:
            syn_samp = random.sample(syn_samp, subset_len)
        print(col_names[x], "\tsynonym in name identified, i.e., (= )\t...removing all content in parens (",len(check_syn), " matches)\n\t\t", ', '.join(syn_samp), sep = '')
        taxa[col_names[x]] = [re.sub('([A-Z][^A-Z\(]+)(.*) \(=.*\)', r'\1\2', v) for v in taxa[col_names[x]]]
    if cv_bool | (col_names[x] != "Cultivar"):   
        check_eq = return_matches(arr = taxa[col_names[x]], regexp = ' ?= ?')
        if len(check_eq) > 0:
            eq_samp = list(set(check_eq))
            if len(eq_samp) > subset_len:
                eq_samp = random.sample(eq_samp, subset_len)
            print(col_names[x], "\tequals sign in name identified\t...removing (", len(check_eq), " matches)\n\t\t", ', '.join(eq_samp), sep = '')
            taxa[col_names[x]] = [re.sub(' ?= ?', '', v) for v in taxa[col_names[x]]]
    check_hybrid = return_matches(arr = taxa[col_names[x]], regexp = '(?: |^)[xX] ')
    if len(check_hybrid) > 0:
        hybrid_samp = list(set(check_hybrid))
        if len(hybrid_samp) > subset_len:
            hybrid_samp = random.sample(hybrid_samp, subset_len)
        print(col_names[x], "\talternative hybrid symbol identified\t...replacing with \u00D7 (", len(check_hybrid), " matches)\n\t\t", ', '.join(hybrid_samp), sep = '')
        taxa[col_names[x]] = [re.sub('^[xX] ', '\u00D7 ', v) for v in taxa[col_names[x]]]
        taxa[col_names[x]] = [re.sub(' [xX] ', ' \u00D7 ', v) for v in taxa[col_names[x]]]
    check_hyb_spacing = return_matches(arr = taxa[col_names[x]], regexp = '(?: |^)[\u00D7+][A-Za-z]')
    if len(check_hyb_spacing) > 0:
        hyb_spacing_samp = list(set(check_hyb_spacing))
        if len(hyb_spacing_samp) > subset_len:
            hyb_spacing_samp = random.sample(hyb_spacing_samp, subset_len)
        print(col_names[x], "\talternative hybrid spacing identified\t...standardizing (", len(check_hyb_spacing), " matches)\n\t\t", ', '.join(hyb_spacing_samp), sep = '')
        taxa[col_names[x]] = [re.sub('^([\u00D7+])([A-Za-z])', '\\1 \\2', v) for v in taxa[col_names[x]]]
        taxa[col_names[x]] = [re.sub(' \u00D7([a-z])', ' \u00D7 \\1', v) for v in taxa[col_names[x]]]
#
# replace genera with corrected genera
# older documents in particular are likely to use alternative spelling
print("\ncorrecting alternative genera spelling")
genus_filename = path + "/GeneraAliases.csv"
genus = pd.read_csv(genus_filename, delimiter = ',', encoding = 'utf-8')
genus = pd.Series(genus.CorrectedGenus.values, index = genus.MisspelledGenus).to_dict()
taxa['Genus'] = replace_dict_exact(arr = taxa['Genus'], dic = genus)
#
# check whether any genera are potential fungal genera using GBIF list
print("\nchecking for fungal genera using GBIF list\n\tnote: some fungi genus taxon names are also vascular plant genus taxon names")
fungi_filename = path + "/GBIF_FungiGenera.csv"
fungi_gen = pd.read_csv(fungi_filename, delimiter = ',', encoding = 'utf-8')
fungi_bool = taxa['Genus'].isin(fungi_gen['taxon_name'])
fungi_idx = [i for i, x in enumerate(fungi_bool) if x]
if len(fungi_idx) > 0:
    match_fungi = taxa['Genus'].iloc[fungi_idx]
    match_fungi = match_fungi.unique()
    match_fungi = np.sort(match_fungi)
    print("\tNOTE: genera found that match fungi genus names! check manually!\n\t", '\n\t'.join(match_fungi), sep = '')
else:
    match_fungi = []
    print("\tno potential fungal genera identified")
#
# check whether any genera are potential algal genera using MOBOT list
print("\nchecking for algae genera using MOBOT list\n\tnote: some algae genus taxon names are also vascular plant genus taxon names")
alga_filename = path + "/GBIF_AlgaeGenera.csv"
alga_gen = pd.read_csv(alga_filename, delimiter = ',', encoding = 'utf-8')
alga_bool = taxa['Genus'].isin(alga_gen['taxon_name'])
alga_idx = [i for i, x in enumerate(alga_bool) if x]
if len(alga_idx) > 0:
    match_alga = taxa['Genus'].iloc[alga_idx]
    match_alga = match_alga.unique()
    match_alga = np.sort(match_alga)
    print("\tNOTE: genera found that match algae genus names! check manually!\n\t", '\n\t'.join(match_alga), sep = '')
else:
    match_alga = []
    print("\tno potential algae genera identified")
#
# check whether any genera are potential bryophyte genera using MOBOT list
print("\nchecking for bryophyte genera using MOBOT list\n\tnote: some bryophyte genus taxon names are also vascular plant genus taxon names")
bryo_filename = path + "/MOBOT_BryophyteGenera.csv"
bryo_gen = pd.read_csv(bryo_filename, delimiter = ',', encoding = 'utf-8')
bryo_bool = taxa['Genus'].isin(bryo_gen['taxon_name'])
bryo_idx = [i for i, x in enumerate(bryo_bool) if x]
if len(bryo_idx) > 0:
    match_bryo = taxa['Genus'].iloc[bryo_idx]
    match_bryo = match_bryo.unique()
    match_bryo = np.sort(match_bryo)
    print("\tNOTE: genera found that match bryophyte genus names! check manually!\n\t", '\n\t'.join(match_bryo), sep = '')
else:
    match_bryo = []
    print("\tno potential bryophyte genera identified")
#
# check whether all genera are present in APG IV and give warnings for mismatches
print("\nchecking for genera in WCVP")
wcvp_filename = path + "/WCVP_Genera_Jun2021.csv"
wcvp_gen = pd.read_csv(wcvp_filename, delimiter = '|', encoding = 'utf-8')
gen_bool = taxa['Genus'].isin(wcvp_gen['taxon_name'])
gen_idx = [i for i, x in enumerate(gen_bool) if not x]
if len(gen_idx) > 0:
    match_gen = taxa['Genus'].iloc[gen_idx]
    match_gen = match_gen.unique()
    match_gen = np.sort(match_gen)
    print("\tNOTE: genera found that are not present in WCVP! check manually!\n\t", '\n\t'.join(match_gen), sep = '')
    # check whether the genera that are not present in WCVP are also likely fungi, algae, or bryophytes
    print("\nchecking whether the genera that are not present in WCVP are likely fungi, algae, or bryophytes")
    match_nonvasc = np.concatenate([match_fungi, match_alga, match_bryo])
    match_vasc = [x for x in match_gen if x not in match_nonvasc]
    if len(match_vasc) > 0:
        match_vasc = np.sort(match_vasc)
        print("\tNOTE: genera found that are not present in WCVP that are also not likely fungi, algae, or bryophytes! check manually!\n\t", '\n\t'.join(match_vasc), sep = '')
    else:
        print("\tall genera that are not present in WCVP may be fungi, algae, or bryophytes")
else:
    print("\tall genera were found in WCVP")
#
if ('Authority' in col_names):
    print("\nstandardizing author names")
    # author initials should not include spaces
    check_author_inits = return_matches(arr = taxa['Authority'], regexp = r'(?:[A-Z]\. )*[A-Z]\. [A-Z]')
    if len(check_author_inits) > 0:
        author_init_samp = list(set(check_author_inits))
        if len(author_init_samp) > subset_len:
            author_init_samp = random.sample(author_init_samp, subset_len)
        print("\tnonstandard spacing of author initials\t...standardizing (", len(check_author_inits), " matches)\n\t\t", ', '.join(author_init_samp), sep = '')
        taxa['Authority'] = [re.sub(r'(?:([A-Z]\.) )?(?:([A-Z]\.) )?([A-Z]\.) ([A-Z])', '\\1\\2\\3\\4', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_author_inits = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r'(?:[A-Z]\. )*[A-Z]\. [A-Z]')
        if len(check_author_inits) > 0:
            author_init_samp = list(set(check_author_inits))
            if len(author_init_samp) > subset_len:
                author_init_samp = random.sample(author_init_samp, subset_len)
            print("\tnonstandard spacing of infraspecific author initials\t...standardizing (", len(check_author_inits), " matches)\n\t\t", ', '.join(author_init_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(r'(?:([A-Z]\.) )?(?:([A-Z]\.) )?([A-Z]\.) ([A-Z])', '\\1\\2\\3\\4', v) for v in taxa['InfraspecificAuthority']]
    # and to ampersand
    check_author_amp = return_matches(arr = taxa['Authority'], regexp = r' and ')
    if len(check_author_amp) > 0:
        author_amp_samp = list(set(check_author_amp))
        if len(author_amp_samp) > subset_len:
            author_amp_samp = random.sample(author_amp_samp, subset_len)
        print("\t'and' instead of ampersand in author\t...standardizing (", len(check_author_amp), " matches)\n\t\t", ', '.join(author_amp_samp), sep = '')
        taxa['Authority'] = [re.sub(r' and ', ' & ', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_author_amp = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r' and ')
        if len(check_author_amp) > 0:
            author_amp_samp = list(set(check_author_amp))
            if len(author_amp_samp) > subset_len:
                author_amp_samp = random.sample(author_amp_samp, subset_len)
            print("\t'and' instead of ampersand in infraspecific author\t...standardizing (", len(check_author_amp), " matches)\n\t\t", ', '.join(author_amp_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(r' and ', ' & ', v) for v in taxa['InfraspecificAuthority']]
    # no period after ex, spaces around ex, lowercase
    check_ex = return_matches(arr = taxa['Authority'], regexp = r'(\.)?(?(1)| )[Ee]x\.? ')
    if len(check_ex) > 0:
        ex_samp = list(set(check_ex))
        if len(ex_samp) > subset_len:
            ex_samp = random.sample(ex_samp, subset_len)
        print("\t'ex' in author citation\t...standardizing (", len(check_ex), " matches)\n\t\t", ', '.join(ex_samp), sep = '')
        taxa['Authority'] = [re.sub(r'(\.)?(?(1)| )[Ee]x\.? ', '\\1 ex ', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_ex = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r'(\.)?(?(1)| )[Ee]x\.? ')
        if len(check_ex) > 0:
            ex_samp = list(set(check_ex))
            if len(ex_samp) > subset_len:
                ex_samp = random.sample(ex_samp, subset_len)
            print("\t'ex' in infraspecific author citation\t...standardizing (", len(check_ex), " matches)\n\t\t", ', '.join(ex_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(r'(\.)?(?(1)| )[Ee]x\.? ', '\\1 ex ', v) for v in taxa['InfraspecificAuthority']]
    # 46.2 When authorship of a name differs from authorship of the publication in which it was validly published, both are sometimes cited, connected by the word “in”. In such a case, “in” and what follows are part of a bibliographic citation and are better omitted unless the place of publication is being cited.
    in_regexp = r' in .*?(\)|$)'
    check_in = return_matches(arr = taxa['Authority'], regexp = in_regexp)
    if len(check_in) > 0:
        in_samp = list(set(check_in))
        if len(in_samp) > subset_len:
            in_samp = random.sample(in_samp, subset_len)
        print("\tpublication included in author citation\t...removing all content including and following 'in' (", len(check_in), " matches)\n\t\t", ', '.join(in_samp), sep = '')
        taxa['Authority'] = [re.sub(in_regexp, '\\1', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_in = return_matches(arr = taxa['InfraspecificAuthority'], regexp = in_regexp)
        if (len(check_in) > 0):
            in_samp = list(set(check_in))
            if len(in_samp) > subset_len:
                in_samp = random.sample(in_samp, subset_len)
            print("\tpublication included in infraspecific author citation\t...removing all content including and following 'in' (", len(check_in), " matches)\n\t\t", ', '.join(in_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(in_regexp, '\\1', v) for v in taxa['InfraspecificAuthority']] 
    # some authorities are listed to unnecessary specificity for the purposes of harmonization, e.g., "Symphyotrichum longifolium sensu G.L. Nesom, non Lam." for these cases, remove all content following "non." this can be particularly misleading, because the standardization method will match to the excluded authority.
    check_anon = return_matches(arr = taxa['Authority'], regexp = r'(^| )auct\.? non .*$')
    if len(check_anon) > 0:
        anon_samp = list(set(check_anon))
        if len(anon_samp) > subset_len:
            anon_samp = random.sample(anon_samp, subset_len)
        print("\t'auct. non' included in author citation\t...removing this and following content (", len(check_anon), " matches\n\t\t", ', '.join(anon_samp), sep = '')
        taxa['Authority'] = [re.sub(r'(^| )auct\.? non .*$', '', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_anon = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r'(^| )auct\.? non .*$')
        if (len(check_anon) > 0):
            anon_samp = list(set(check_anon))
            if len(anon_samp) > subset_len:
                anon_samp = random.sample(anon_samp, subset_len)
            print("\t'auct. non' included in infraspecific author citation\t...removing this and following content (", len(check_anon), " matches)\n\t\t", ', '.join(anon_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(r' auct\.? non .*$', '', v) for v in taxa['InfraspecificAuthority']]
    check_non = return_matches(arr = taxa['Authority'], regexp = r'(^| )non .*$')
    if len(check_non) > 0:
        non_samp = list(set(check_non))
        if len(non_samp) > subset_len:
            non_samp = random.sample(non_samp, subset_len)
        print("\t'non' included in author citation\t...removing this and following content (", len(check_non), " matches)\n\t\t", ', '.join(non_samp), sep = '')
        taxa['Authority'] = [re.sub(r'(^| |\()non .*$', '', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_non = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r'(^| )non .*$')
        if (len(check_non) > 0):
            non_samp = list(set(check_non))
            if len(non_samp) > subset_len:
                non_samp = random.sample(non_samp, subset_len)
            print("\t'non' included in infraspecific author citation\t...removing this and following content (", len(check_non), " matches)\n\t\t", ', '.join(non_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(r'(^| |\()non .*$', '', v) for v in taxa['InfraspecificAuthority']]
        # 47&47A An alteration of the diagnostic characters or of the circumscription of a taxon without the exclusion of the type does not warrant a change of the author citation of the name of the taxon.  When an alteration as mentioned in Art. 47 has been considerable, the nature of the change may be indicated by adding such words, abbreviated where suitable.
    alt_regexp = '(^| | ?\()(p[\. ]{1,2}p\.?|pro hybr\.?|pro sp\.?|pro syn\.?|pro syn\.?|pro nm\.?|nome?n?[\. ]{1,2}inval\.?|nome?n?[\. ]{1,2}ambig\.?|nome?n?[\. ]{1,2}illeg\.?|nome?n?[\. ]{1,2}rej\.?|nome?n?[\. ]{1,2}inq\.?|nome?n?[\. ]{1,2}prov\.?|nome?n?[\. ]{1,2}obsc\.?|nome?n?[\. ]{1,2}nud\.?|nome?n?[\. ]{1,2}dub\.?|nome?n?[\. ]{1,2}cons?[\. ]{1,2}prop\.?|nom[\. ]{1,2}conserv\.?|nome?n?[\. ]{1,2}cons?\.?|s[\. ]{1,2}ampl\.?|s[\. ]{1,2}l\.?|s[\. ]{1,2}s\.|s[\. ]{1,2}lat\.?|sens[\.u] ?lat[\.o]|sens[\.u] ?ampl[\.o]|s[\. ]{1,2}str\.?|sens[\.u] ?str\.?|sens[\.u] ?strict[\.o]|emend\.?|emendavit|mut[\. ]{1,2}char\.?|ined\.?|mutatis characteribus|excl[\. ]{1,2}gen\.?|excl[\. ]{1,2}sp\.?|excl[\. ]{1,2}var\.?).*'
    check_alt = return_matches(arr = taxa['Authority'], regexp = alt_regexp)
    if len(check_alt) > 0:
        alt_samp = list(set(check_alt))
        if len(alt_samp) > subset_len:
            alt_samp = random.sample(alt_samp, subset_len)
        print("\talteration/circumscription included in author citation (e.g., emend.)\t...removing all all content including and following this (", len(check_alt), " matches)\n\t\t", ', '.join(alt_samp), sep = '')
        taxa['Authority'] = [re.sub(alt_regexp, '', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_alt = return_matches(arr = taxa['InfraspecificAuthority'], regexp = alt_regexp)
        if (len(check_alt) > 0):
            alt_samp = list(set(check_alt))
            if len(alt_samp) > subset_len:
                alt_samp = random.sample(alt_samp, subset_len)
            print("\talteration/circumscription included in infraspecific author citation (e.g., emend.)\t...removing all content including and following this (", len(check_alt), " matches)\n\t\t", ', '.join(alt_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(alt_regexp, '', v) for v in taxa['InfraspecificAuthority']]
    # similarly, the "sensu" should be removed from the authority
    check_sens = return_matches(arr = taxa['Authority'], regexp = r'(^| )sensu ')
    if len(check_sens) > 0:
        sens_samp = list(set(check_sens))
        if len(sens_samp) > subset_len:
            sens_samp = random.sample(sens_samp, subset_len)
        print("\t'sensu' included in author citation\t...removing (", len(check_sens), " matches)\n\t\t", ', '.join(sens_samp), sep = '')
        taxa['Authority'] = [re.sub(r'(^| )sensu ', '\\1', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_sens = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r'(^| )sensu ')
        if (len(check_sens) > 0):
            sens_samp = list(set(check_sens))
            if len(sens_samp) > subset_len:
                sens_samp = random.sample(sens_samp, subset_len)
            print("\t'sensu' included in infraspecific author citation\t...removing (", len(check_sens), " matches)\n\t\t", ', '.join(sens_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(r'(^| )sensu ', '\\1', v) for v in taxa['InfraspecificAuthority']]  
    # remove any text in square brackets
    check_sqbr = return_matches(arr = taxa['Authority'], regexp = r' ?\[.*\]')
    if len(check_sqbr) > 0:
        sqbr_samp = list(set(check_sqbr))
        if len(sqbr_samp) > subset_len:
            sqbr_samp = random.sample(sqbr_samp, subset_len)
        print("\tsquare brackets included in author citation\t...removing brackets and contents (", len(check_sqbr), "matches)\n\t\t", ', '.join(sqbr_samp), sep = '')
        taxa['Authority'] = [re.sub(r' ?\[.*\]', '', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_sqbr = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r' ?\[.*\]')
        if (len(check_sqbr) > 0):
            sqbr_samp = list(set(check_sqbr))
            if len(sqbr_samp) > subset_len:
                sqbr_samp = random.sample(sqbr_samp, subset_len)
            print("\tsquare brackets included in infraspecific author citation\t...removing brackets and contents (", len(check_sqbr), " matches)\n\t\t", ', '.join(sqbr_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(r' ?\[.*\]', '', v) for v in taxa['InfraspecificAuthority']]
    # change any et to ampersand
    check_et = return_matches(arr = taxa['Authority'], regexp = r'\b[Ee]t\b(?! al)')
    if len(check_et) > 0:
        et_samp = list(set(check_et))
        if len(et_samp) > subset_len:
            et_samp = random.sample(et_samp, subset_len)
        print("\t'et' in author citation\t...replacing with ampersand (", len(check_et), " matches)\n\t\t", ', '.join(et_samp), sep = '')
        taxa['Authority'] = [re.sub(r'\b[Ee]t\b(?! al)', '&', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_et = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r'\bet\b(?! al)')
        if len(check_et) > 0:
            et_samp = list(set(check_et))
            if len(et_samp) > subset_len:
                et_samp = random.sample(et_samp, subset_len)
            print("\t'et' in infraspecific author citation\t...replacing with ampersand (", len(check_et), " matches)\n\t\t", ', '.join(et_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(r'\bet\b(?! al)', '&', v) for v in taxa['InfraspecificAuthority']]
    # include spaces around ampersand
    check_amp = return_matches(arr = taxa['Authority'], regexp = r' ?& ?')
    if len(check_amp) > 0:
        amp_samp = list(set(check_amp))
        if len(amp_samp) > subset_len:
            amp_samp = random.sample(amp_samp, subset_len)
        print("\tstandardizing spacing around ampersands in author (", len(check_amp), " matches)\n\t\t", ', '.join(amp_samp), sep = '')
        taxa['Authority'] = [re.sub(r' ?& ?', ' & ', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_amp = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r' ?& ?')
        if len(check_amp) > 0:
            amp_samp = list(set(check_amp))
            if len(amp_samp) > subset_len:
                amp_samp = random.sample(amp_samp, subset_len)
            print("\tstandardizing spacing around ampersands in infraspecific author (", len(check_amp), " matches)\n\t\t", ', '.join(amp_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(r' ?& ?', ' & ', v) for v in taxa['InfraspecificAuthority']]
    # remove misplaced commas (at beginning or end of authority)
    com_regexp = '(^ *, *| *, *$)'
    check_com = return_matches(arr = taxa['Authority'], regexp = com_regexp)
    if len(check_com) > 0:
        com_samp = list(set(check_com))
        if len(com_samp) > subset_len:
            com_samp = random.sample(com_samp, subset_len)
        print("\tcommas at beginning and/or end of authority\t...removing (", len(check_com), " matches)\n\t\t", ', '.join(com_samp), sep = '')
        taxa['Authority'] = [re.sub(com_regexp, '', v) for v in taxa['Authority']]
    if 'InfraspecificAuthority' in col_names:
        check_com = return_matches(arr = taxa['InfraspecificAuthority'], regexp = com_regexp)
        if len(check_com) > 0:
            com_samp = list(set(check_com))
            if len(com_samp) > subset_len:
                com_samp = random.sample(com_samp, subset_len)
            print("\tcommas at beginning and/or end of infraspecific authority\t...removing (", len(check_com), " matches)\n\t\t", ', '.join(com_samp), sep = '')
            taxa['InfraspecificAuthority'] = [re.sub(com_regexp, '', v) for v in taxa['InfraspecificAuthority']]
    # check for any words that start with lowercase letters in authority
    lowercase_regexp = r'(?!hort\.)(?!fil\.)\b[a-z][a-zA-Z0-9]{2,}'
    check_lower = return_matches(arr = taxa['Authority'], regexp = lowercase_regexp)
    if len(check_lower) > 0:
        lower_samp = list(set(check_lower))
        if len(lower_samp) > subset_len:
            lower_samp = random.sample(lower_samp, subset_len)
        print("\tNOTE: identified word(s) in authority that start(s) with a lowercase letter! check manually! (", len(check_lower), " matches)\n\t\t", ', '.join(lower_samp), sep = '')
    if 'InfraspecificAuthority' in col_names:
        check_lower = return_matches(arr = taxa['InfraspecificAuthority'], regexp = lowercase_regexp)
        if len(check_lower) > 0:
            lower_samp = list(set(check_lower))
            if len(lower_samp) > subset_len:
                lower_samp = random.sample(lower_samp, subset_len)
            print("\tNOTE: identified word(s) in infraspecific authority that start(s) with a lowercase letter! check manually! (", len(check_lower), " matches)\n\t\t", ', '.join(lower_samp), sep = '')
    # replace authors written with full names to standard forms
    print("\tconverting authority full names to standard forms")
    auth_filename = path + "/AuthorAliases.csv"
    author = pd.read_csv(auth_filename, sep = ',', encoding = 'utf-8', encoding_errors = 'replace')
    # in replace_author() . is converted to \.?
    author['Alternative'] = author['Alternative'] + '.'
    author = pd.Series(author.Standard.values, index = author.Alternative).to_dict()
    taxa['Authority'] = replace_author(arr = taxa['Authority'], dic = author, allow_inits = True)
    if 'InfraspecificAuthority' in col_names:
        taxa['InfraspecificAuthority'] = replace_author(arr = taxa['InfraspecificAuthority'], dic = author, allow_inits = True)
    if old_text_bool:
        # replace authors using list of manually corrected authors
        print("\tcorrecting alternative author abbreviations")
        oldauth_filename = path + "/AuthorAliases_OldText.csv"
        oldauth = pd.read_csv(oldauth_filename, sep = ',', encoding = 'utf-8', encoding_errors = 'replace')
        oldauth = oldauth[oldauth['Ambiguous'] != 'Ambiguous']
        oldauth = pd.Series(oldauth.Standard.values, index = oldauth.Alternative).to_dict()
        taxa['Authority'] = replace_author(arr = taxa['Authority'], dic = oldauth)
        if 'InfraspecificAuthority' in col_names:
            taxa['InfraspecificAuthority'] = replace_author(arr = taxa['InfraspecificAuthority'], dic = oldauth)
        # replace fil with f.
        check_fil = return_matches(arr = taxa['Authority'], regexp = r'\.?(\.| )fil\.?( |$|\)|,)')
        if len(check_fil) > 0:
            fil_samp = list(set(check_fil))
            if len(fil_samp) > subset_len:
                fil_samp = random.sample(fil_samp, subset_len)
            print("\t'fil' identified in authorities\t...standardizing (", len(check_fil), " matches)\n\t\t", ', '.join(fil_samp), sep = '')
            taxa['Authority'] = [re.sub(r'\.?(\.| )fil\.?( |$|\)|,)', '.f.\\2', v) for v in taxa['Authority']]
        if 'InfraspecificAuthority' in col_names:
            check_ifil = return_matches(arr = taxa['InfraspecificAuthority'], regexp = r'\.?(\.| )fil\.?( |$|\)|,)')
            if len(check_ifil) > 0:
                ifil_samp = list(set(check_ifil))
                if len(ifil_samp) > subset_len:
                    ifil_samp = random.sample(ifil_samp, subset_len)
                print("\t'fil' identified in infraspecific authorities\t...standardizing (", len(check_ifil), " matches)\n\t\t", ', '.join(ifil_samp), sep = '')
                taxa['InfraspecificAuthority'] = [re.sub(r'\.?(\.| )fil\.?( |$|\)|,)', '.f.\\2', v) for v in taxa['InfraspecificAuthority']]
#
# check infraspecific ranks for alternative spellings
if 'InfraspecificName' in col_names:
    if 'InfraspecificRank' in col_names:
        rank_col = 'InfraspecificRank'
    else:
        rank_col = 'InfraspecificName'
    print("\nchecking infraspecific ranks")
    check_subsp = return_matches(arr = taxa[rank_col], regexp = '( |^)(subsp|ssp.?)( |$)')
    if len(check_subsp) > 0:
        print("\talternative infraspecific rank spelling identified\t...standardizing (", len(check_subsp), " matches)\n\t\t", ', '.join(set(check_subsp)), sep = '')
        taxa[rank_col] = [re.sub('( |^)(subsp|ssp.?)( |$)', '\\1subsp.\\3', v) for v in taxa[rank_col]]
    check_convar = return_matches(arr = taxa[rank_col], regexp = '( |^)(convar|conv\.?)( |$)')
    if len(check_convar) > 0:
        print("\talternative infraspecific rank spelling identified\t...standardizing (", len(check_convar), " matches)\n\t\t", ', '.join(set(check_convar)), sep = '')
        taxa[rank_col] = [re.sub('( |^)(convar|conv\.?)( |$)', '\\1convar.\\3', v) for v in taxa[rank_col]]
    check_var = return_matches(arr = taxa[rank_col], regexp = '( |^)(var|v\.?)( |$)')
    if len(check_var) > 0:
        print("\talternative infraspecific rank spelling identified\t...standardizing (", len(check_var), " matches)\n\t\t", ', '.join(set(check_var)), sep = '')
        taxa[rank_col] = [re.sub('( |^)(var|v\.?)( |$)', '\\1var.\\3', v) for v in taxa[rank_col]]
    check_subvar = return_matches(arr = taxa[rank_col], regexp = '( |^)(subvar|subv\.?)( |$)')
    if len(check_subvar) > 0:
        print("\talternative infraspecific rank spelling identified\t...standardizing (", len(check_subvar), " matches)\n\t\t", ', '.join(set(check_subvar)), sep = '')
        taxa[rank_col] = [re.sub('( |^)(subvar|subv\.?)( |$)', '\\1subvar.\\3', v) for v in taxa[rank_col]]
    check_f = return_matches(arr = taxa[rank_col], regexp = '( |^)(f[ao]?\.?|forma\.?|fma\.?)( |$)')
    if len(check_f) > 0:
        print("\talternative infraspecific rank spelling identified\t...standardizing (", len(check_f), " matches)\n\t\t", ', '.join(set(check_f)), sep = '')
        taxa[rank_col] = [re.sub('( |^)(f[ao]?\.?|forma\.?|fma\.?)( |$)', '\\1f.\\3', v) for v in taxa[rank_col]]
    check_subf = return_matches(arr = taxa[rank_col], regexp = '( |^)(subf[ao]?\.?|subforma\.?|subfma\.?)( |$)')
    if len(check_subf) > 0:
        print("\talternative infraspecific rank spelling identified\t...standardizing (", len(check_subf), " matches)\n\t\t", ', '.join(set(check_subf)), sep = '')
        taxa[rank_col] = [re.sub('( |^)(subf[ao]?\.?|subforma\.?|subfma\.?)( |$)', '\\1subf.\\3', v) for v in taxa[rank_col]]
    # replace greek letters in infraspecific names (former way to indicate var.)
    check_greek = return_matches(arr = taxa[rank_col], regexp = '[\u03B1\u03B2\u03B3\u03B4\u03B5\u03B6\u03B7\u03B8\u03BB\u03C1]')
    if len(check_greek) > 0:
        print("\talternative intraspecific rank spelling identified\t...standardizing (", len(check_greek), " matches)\n\t\t", ', '.join(set(check_greek)), sep = '')
        taxa[rank_col] = [re.sub('[\u03B1\u03B2\u03B3\u03B4\u03B5\u03B6\u03B7\u03B8\u03BB\u03C1]\.?', 'var.', v) for v in taxa[rank_col]]
    print("\nchecking infraspecific names")
    # remove square brackets
    check_square = return_matches(arr = taxa['InfraspecificName'], regexp = ' ?\[.*?\]')
    if len(check_square) > 0:
        print("\tsquare brackets identified\t...removing brackets and contents (", len(check_square), " matches)\n\t\t", ', '.join(set(check_square)), sep = '')
        taxa['InfraspecificName'] = [re.sub(' ?\[.*?\]', '', v) for v in taxa['InfraspecificName']] 
    # remove parentheses
    check_parens = return_matches(arr = taxa['InfraspecificName'], regexp = ' ?\(.*?\)')
    if len(check_parens) > 0:
        print("\tparentheses identified\t...removing parens and contents (", len(check_parens), " matches)\n\t\t", ', '.join(set(check_parens)), sep = '')
        taxa['InfraspecificName'] = [re.sub(' ?\(.*?\)', '', v) for v in taxa['InfraspecificName']]
    # remove certain elements of names
    # taxonomic circumscription: p.p. (max/min, taxon includes more than one currently recognized entity), sens. lat., sens. strict., nom. cons. (conserved name)
    circ_regexp = r'(^| )(spp[\. ]{1,2}aggr?\.?|aggr?\.?|nom[\. ]{1,2}inval\.?|nom[\. ]{1,2}illeg\.?|nom[\. ]{1,2}rej\.?|nom[\. ]{1,2}inq\.?|nom[\. ]{1,2}obsc\.?|nom[\. ]{1,2}nud\.?|nom[\. ]{1,2}dub\.?|nom[\. ]{1,2}cons?[\. ]{1,2}prop\.?|nom[\. ]{1,2}conserv\.?|nom[\. ]{1,2}cons?\.?|p[\. ]{1,2}p\.?|s[\. ]{1,2}ampl\.?|s[\. ]{1,2}l\.?|s[\. ]{1,2}s\.?|s[\. ]{1,2}lat\.?|sens[\.u] ?lat[\.o]|sens[\.u] ?ampl[\.o]|s[\. ]{1,2}str\.?|sens[\.u] ?str\.?|sens[\.u] ?strict[\.o]|mut[\. ]{1,2}char\.?|excl[\. ]{1,2}gen\.?|excl[\. ]{1,2}sp\.?|excl[\. ]{1,2}var\.?)(?:$| )'
    check_tax = return_matches(arr = taxa['InfraspecificName'], regexp = circ_regexp)
    if len(check_tax) > 0:
        print("\textraneous taxonomic circumscriptions identified\t...extracting (", len(check_tax), " matches)\n\t\t", ', '.join(set(check_tax)), sep = '')
        taxa['InfraspecificName'] = [re.sub(circ_regexp, '', v) for v in taxa['InfraspecificName']]
    if old_text_bool:
        # old texts often list varieties like they are valid species names, e.g. flore-pleno
        # check for these and move to cultivar column if detected
        print("\nchecking for cultivars in infraspecific names")
        cult_regexp = '(?=(?<!subsp)\.)(?=(?<!^convar)\.)(?=(?<!^var)\.)(?=(?<!subvar)\.)(?=(?<!f)\.)(?=(?<!subf)\.)\..|(?<!subsp\. )(?<!convar\. )(?<!^var\. )(?<!subvar\. )(?<!f\. )(?<!subf\. )flore[ -]|(?<!subsp\. )(?<!^convar\. )(?<!^var\. )(?<!subvar\. )(?<!f\. )(?<!subf\. )foliis[ -]|(?<!subsp\. )(?<!^convar\. )(?<!^var\. )(?<!subvar\. )(?<!f\. )(?<!subf\. )fol\.?[ -]|(?<!subsp\. )(?<!^convar\. )(?<!^var\. )(?<!subvar\. )(?<!f\. )(?<!subf\. )fl\.?[ -]|(?<!subsp\. )(?<!^convar\. )(?<!^var\. )(?<!subvar\. )(?<!f\. )(?<!subf\. )fruct[uo][ -]|^cvs?\.? ?[\'"]?|^cvgr\.? ?[\'"]?|^[\'"]'
        check_cult = return_matches(arr = taxa['InfraspecificName'], regexp = cult_regexp)
        cult_idx = return_matches(arr = taxa['InfraspecificName'], mode = 'index', regexp = cult_regexp)
        if 'InfraspecificRank' in col_names:
            rank_idx = [i for i, v in enumerate(taxa['InfraspecificRank']) if v != '']
            cult_idx = [x for x in cult_idx if x not in rank_idx]
        check_cult = taxa.loc[cult_idx, 'InfraspecificName']
        if len(check_cult) > 0:
            cult_samp = list(set(check_cult))
            if len(cult_samp) > subset_len:
                cult_samp = random.sample(cult_samp, subset_len)
            print("\tsuspected cultivar names to cultivar column (", len(check_cult), " matches)\n\t\t", ', '.join(cult_samp), sep = '')
            if 'Cultivar' not in col_names:
                taxa['Cultivar'] = [str() for c in 'c' * taxa.shape[0]]
                col_names.append('Cultivar')
            taxa.loc[cult_idx, 'Cultivar'] = taxa.loc[cult_idx, 'InfraspecificName']
            taxa.loc[cult_idx, 'InfraspecificName'] = ''
            if 'InfraspecificRank' in col_names:
                taxa.loc[cult_idx, 'InfraspecificRank'] = ''
    print("\nchecking for infraspecific ranks in infraspecific names")
    infrank_regexp = '(?:^| )(var\.|subvar\.|convar\.|subsp\.|f\.|subf\.)(?:$| |[a-z])'
    check_infrank = return_matches(arr = taxa['InfraspecificName'], regexp = infrank_regexp)
    if len(check_infrank) > 0:
        rank_infsamp = list(set(check_infrank))
        if len(rank_infsamp) > subset_len:
            rank_infsamp = random.sample(rank_infsamp, subset_len)
        print("\tinfraspecific ranks to rank column (", len(check_infrank), " matches)\n\t\t", ', '.join(rank_infsamp), sep = '')
        sep_output = separate_infra(arr = taxa['InfraspecificName'])
        taxa['InfraspecificName'] = sep_output[1]
        if 'InfraspecificRank' not in col_names:
            taxa['InfraspecificRank'] = sep_output[0]
        elif 'InfraspecificRank' in col_names:
            taxa['InfraspecificRank'] = taxa['InfraspecificRank'] + sep_output[0]
    elif 'InfraspecificRank' not in col_names:
        taxa['InfraspecificRank'] = [str() for c in 'c' * taxa.shape[0]]
    acc_rank = ["", "subsp.", "var.", "convar.", "subvar.", "f.", "subf."]
    rank_bool = ~taxa['InfraspecificRank'].isin(acc_rank)
    acc_rank_idx = [i for i, x in enumerate(rank_bool) if x]
    if len(acc_rank_idx) > 0:
        match_rank = taxa['InfraspecificRank'].iloc[acc_rank_idx]
        match_rank = match_rank.unique()
        match_rank = np.sort(match_rank)
        print("\nNOTE: nonstandard infraspecific ranks detected after cleaning! check manually!\n\t", '\n\t'.join(match_rank), sep = '')
else:
    taxa['InfraspecificRank'] = [str() for c in 'c' * taxa.shape[0]]
    taxa['InfraspecificName'] = [str() for c in 'c' * taxa.shape[0]]
#  
# clean cultivar names
if 'Cultivar' not in col_names:
    print("\ncultivar column not detected")
    taxa['Cultivar'] = [str() for c in 'c' * taxa.shape[0]]
if cv_bool:
    print("\nstandardizing cultivar names, e.g., 'Cultivar Name'")
    initials_regexp = '[A-Z]{2,}'
    check_initials = return_matches(arr = taxa['Cultivar'], regexp = initials_regexp)
    if len(check_initials) > 0:
        initials_samp = list(set(check_initials))
        if len(initials_samp) > subset_len:
            initials_samp = random.sample(initials_samp, subset_len)
        print("\tcultivar names with initials identified (", len(check_initials), " matches).\n\tNOTE: these cultivar names will not be converted to title case!\n\t\t", ', '.join(initials_samp), sep = '')
    # cultivar names standardized 'Cultivar Name'
    for x in range(len(taxa['Cultivar'])):
        this_cv = taxa.loc[x, 'Cultivar']
        if this_cv != '':
            # remove cv. (not standardized)
            this_cv = re.sub("\"|cvs?\.? |cvs?\. ?", "", this_cv)
            # remove apostrophes (except when they represent possession, articles, etc.)
            this_cv = re.sub("(?<![dlomDLOM])(?<!'n)(?<!s )(?<!'[DLOM])'(?!s )(?!n')|^['\"]|['\"]$", "", this_cv)
            # if cultivar name includes initials, do not change to title case
            check_initials = bool(re.search(initials_regexp, this_cv))
            if check_initials:
                this_cv = "".join(["'", this_cv, "'"])
            else:
                this_cv = "".join(["'", this_cv.title(), "'"])
            # some letters/words should not be capitalized in cultivar names, depending on language
            # apostrophe s
            this_cv = re.sub("([A-Za-zÀ-ÿĀ-ſƀ-ȳ])'S( |'$)", "\\1's\\2", this_cv)
            # d' or l' in romance languages
            this_cv = re.sub("D'([A-Za-zÀ-ÿĀ-ſƀ-ȳ])", "d'\\1", this_cv)
            this_cv = re.sub(" L'([A-Za-zÀ-ÿĀ-ſƀ-ȳ])", " l'\\1", this_cv)
            # articles like the, de, des, di, der, etc.
            this_cv = re.sub(" D([eui][sl]?) ", " d\\1 ", this_cv)
            this_cv = re.sub(" D([ei][ermn]) ", " d\\1 ", this_cv)
            this_cv = re.sub(" The ", " the ", this_cv)
            # prepositions and conjuctions
            this_cv = re.sub(" Of ", " of ", this_cv)
            this_cv = re.sub(" A ", " a ", this_cv)
            this_cv = re.sub(" \u00C0 ", " \u00E0 ", this_cv)
            this_cv = re.sub(" And ", " and ", this_cv)
            this_cv = re.sub(" Und ", " und ", this_cv)
            this_cv = re.sub(" V([ao])n ", " v\\1n ", this_cv)
            this_cv = re.sub(" In ", " in ", this_cv)
            this_cv = re.sub(" On ", " on ", this_cv)
            this_cv = re.sub(" En ", " en ", this_cv)
            # mc in irish names
            this_cv = re.sub("(?<=[ ']Mc)([A-Za-zÀ-ÿĀ-ſƀ-ȳ])", lambda pat: pat.group(1).upper(), this_cv)
            taxa.loc[x, 'Cultivar'] = this_cv
else:
    print("\nnot standardizing cultivar names, they will remain as is")
# is taxon a cultivar
cultivar = [bool(row['Cultivar']) for index, row in taxa.iterrows()]
# create final author name
if 'Authority' not in col_names:
    print("\nauthority column not detected")
    taxa['Authority'] = [str() for c in 'c' * taxa.shape[0]]
if 'InfraspecificAuthority' not in col_names:
    print("\ninfraspecific authorities not detected")
    auth_edit = taxa['Authority']
else:
    print("\ninfraspecific and specific authorities detected\n\tfinal authority set to authority of the finest taxonomic level")
    auth_edit = []
    for index, row in taxa.iterrows():
        if row['InfraspecificName']:
            auth_edit.append(row['InfraspecificAuthority'])
        else:
            auth_edit.append(row['Authority'])
taxa['AuthorityCorrected'] = auth_edit
# show autonyms, these can influence alignment across databases
check_autonym = sum(taxa['SpecificEpithet'] == taxa['InfraspecificName'])
if check_autonym > 0:
    print("\nNOTE:", check_autonym, "autonyms in taxon list")
#
# checking genera and epithets for spaces and replacing with hyphens according to IAPT Shenzhen Code 20.3 "The name of a genus may not consist of two words, unless these words are joined by a hyphen" and 23.1 "The name of a species is a binary combination consisting of the name of the genus followed by a single specific epithet in the form of an adjective, a noun in the genitive, or a word in apposition. If an epithet consisted originally of two or more words, these are to be united or hyphenated."
print("\nchecking for genera and epithets with multiple words")
multiword_regexp = '(?<![\u00D7+]) (?![\u00D7+])'
check_mwgen = return_matches(arr = taxa['Genus'], regexp = multiword_regexp)
if len(check_mwgen) > 0:
    print("\tmulti-word genera\t...replacing with hyphen (", len(check_mwgen), " matches)\n\t\t", ', '.join(set(check_mwgen)), sep = '')
    taxa['Genus'] = [re.sub(multiword_regexp, r'-', v) for v in taxa['Genus']]
check_epith = return_matches(arr = taxa['SpecificEpithet'], regexp = multiword_regexp)
if len(check_epith) > 0:
    print("\tmulti-word specific epithets\t...replacing with hyphen (", len(check_epith), " matches)\n\t\t", ', '.join(set(check_epith)), sep = '')
    taxa['SpecificEpithet'] = [re.sub(multiword_regexp, r'-', v) for v in taxa['SpecificEpithet']]
check_infepi = return_matches(arr = taxa['InfraspecificName'], regexp = multiword_regexp)
if len(check_infepi) > 0:
    print("\tmulti-word infraspecific epithets\t...replacing with hyphen (", len(check_infepi), " matches)\n\t\t", ', '.join(set(check_infepi)), sep = '')
    taxa['InfraspecificName'] = [re.sub(multiword_regexp, r'-', v) for v in taxa['InfraspecificName']]
#
# removing periods in epithets
print("\nchecking for periods in species epithets")
check_per = return_matches(arr = taxa['SpecificEpithet'], regexp = '\.')
if len(check_per) > 0:
    print("\tperiods in specific epithets\t...removing (", len(check_per), " matches)\n\t\t", ', '.join(set(check_per)), sep = '')
    taxa['SpecificEpithet'] = [re.sub('\.', r'', v) for v in taxa['SpecificEpithet']]
print("\nchecking for periods in infraspecific epithets")
check_infper = return_matches(arr = taxa['InfraspecificName'], regexp = '\.')
if len(check_infper) > 0:
    print("\tperiods in infraspecific epithets\t...removing (", len(check_infper), " matches)\n\t\t", ', '.join(set(check_infper)), sep = '')
    taxa['InfraspecificName'] = [re.sub('\.', r'', v) for v in taxa['InfraspecificName']]
#
# check for cases where infraspecific rank is present with no name
check_rank = taxa['InfraspecificRank'].str.strip().astype(bool)
check_name = taxa['InfraspecificName'].str.strip().astype(bool)
missing_name = []
for i in range(len(check_rank)):
    # rank present, name missing
    this_name = check_rank[i] & ~check_name[i]
    if this_name:
        taxa.loc[i, 'InfraspecificRank'] = ''
        missing_name.append(i)
print("\nremoved", len(missing_name), "infraspecific ranks because no infraspecific name was present")
#
taxa.rename(columns = {'TaxonName': 'OriginalName'}, inplace = True)
if 'InfraspecificAuthority' not in col_names:
    taxa['InfraspecificAuthority'] = [str() for c in 'c' * taxa.shape[0]]
    # create final combined taxon names
    taxa['TaxonNameCV'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'Cultivar']].agg(' '.join, axis = 1)
    taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName']].agg(' '.join, axis = 1)
    taxa['TaxonNameAuthorCV'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'Authority', 'Cultivar']].agg(' '.join, axis = 1)
    taxa['TaxonNameAuthor'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'Authority']].agg(' '.join, axis = 1) 
else:
    taxa['TaxonNameCV'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'Cultivar']].agg(' '.join, axis = 1)
    taxa['TaxonName'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName']].agg(' '.join, axis = 1)
    taxa['TaxonNameAuthorCV'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'AuthorityCorrected', 'Cultivar']].agg(' '.join, axis = 1)
    taxa['TaxonNameAuthor'] = taxa[['Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'AuthorityCorrected']].agg(' '.join, axis = 1)
# final time removing extra spaces
for col in ['TaxonNameCV', 'TaxonName', 'TaxonNameAuthorCV', 'TaxonNameAuthor']:
    taxa[col] = [re.sub(' +', ' ', v) for v in taxa[col]]
    taxa[col] = [v.strip() for v in taxa[col]]
#
# is taxon an infraspecific hybrid?
hybrid_regexp = '(?: |^)(?:\u00D7|nothosp|grex|gx)(?: |$)'
hybrid_infra = return_matches(arr = taxa['InfraspecificName'], regexp = hybrid_regexp, mode = "bool")
if any(hybrid_infra):
    print('\nWARNING: hybrid symbols in infraspecific name!')
# rearrange columns
cols_save = ['index', 'Genus', 'SpecificEpithet', 'InfraspecificRank', 'InfraspecificName', 'Authority', 'InfraspecificAuthority', 'AuthorityCorrected', 'Cultivar', 'OriginalName', 'TaxonNameAuthorCV', 'TaxonNameAuthor','TaxonNameCV', 'TaxonName']
cols_save = cols_save + other_cols
# remove duplicate rows
dup_cols = cols_save[1:] + other_cols
rownum_init = len(taxa)
taxa_final = taxa.drop_duplicates(subset = dup_cols)
rownum_final = len(taxa_final)
if rownum_init - rownum_final > 0:
    dup_taxa = taxa.loc[taxa.duplicated(subset = dup_cols)]['OriginalName']
    dup_taxa = dup_taxa.tolist()
    print("\n", rownum_init - rownum_final, "duplicate rows removed\n\t", ', '.join(dup_taxa))
# save output
taxa_final = taxa_final[cols_save]
fileoutput = [filename, "_TextCleaned_", datetime.now().strftime("%Y%m%d"), ".csv"]
fileoutput = "".join(fileoutput)
print("\nsaving final cleaned file to\t", fileoutput)
taxa_final.to_csv(fileoutput, lineterminator = '\n', index = False, quoting = csv.QUOTE_NONNUMERIC)