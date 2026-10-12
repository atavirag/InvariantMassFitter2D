"""
Script for the preparation of the samples for the
training of ML models to be used in the HF triggers

\author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
\author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
\author Biao Zhang <biao.zhang@cern.ch>, CCNU
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import uproot
import pandas as pd
from alive_progress import alive_bar
from ROOT import TFile, gRandom  # pylint: disable=no-name-in-module

BITS3P = {"DplusToPiKPi": 0,
          "LcToPKPi": 1,
          "DsToKKPi": 2,
          "XicToPKPi": 3}

CHANNELS3P = {"DplusToPiKPi": 1,
              "DsToKKPi": 2,
              "LcToPKPi": 3,
              "XicToPKPi": 4}

# pylint: disable=too-many-locals
def do_dca_smearing(input_df, n_prongs=None):
    """
    Method to do the DCA smearing for 2 prongs and 3 prongs

    Parameters
    -----------------
    - input_df: pandas dataframe containing all candidates with fFlagOrigin column
    - n_prongs: option to 2 prongs or 3 prongs

    Outputs
    -----------------
    - input_df: New dataframe with the smeared DCA columns
    """

    print(f"Start to do the smearing for {n_prongs} prong")
    # Open the input files
    file_names_data = ["New_DCA_smear/results_dcaXY_LHC22m_pass4_run523308.root",
                       "New_DCA_smear/results_dcaZ_LHC22m_pass4_run523308.root"]
    file_names_mc = ["New_DCA_smear/output_DCAxy_all.root",
                     "New_DCA_smear/output_DCAz_all.root"]
    input_files_data = [TFile.Open(name) for name in file_names_data]
    input_files_mc = [TFile.Open(name) for name in file_names_mc]
    input_files_mean = [TFile.Open(name) for name in file_names_data]

    # Extract DCA resolution histograms
    dca_reso_data, dca_reso_mc, dca_reso_data_meanshift = {}, {}, {}
    for i, par in enumerate(["XY", "Z"]):
        dca_reso_data[par] = input_files_data[i].Get("tge_DCA_res_withoutPVrefit_all")
        dca_reso_mc[par] = input_files_mc[i].Get("tge_DCA_res_withoutPVrefit_all")
        dca_reso_data_meanshift[par] = input_files_mean[i].Get("tge_DCA_res_withoutPVrefit_all")

    # Add smeared DCA columns to the dataframe
    smear_cols = ["XY1", "XY2", "Z1", "Z2"]
    if n_prongs == 3:
        smear_cols.extend(["XY3", "Z3"])
    for col in smear_cols:
        dca_col = f"fDCAPrim{col}"
        pt_col = f"fPT{col[-1]}"

        smear_values = []
        for dca, pt in zip(input_df[dca_col], input_df[pt_col]):
            if pt < 7:
                dca_reso_data_mean_shift = dca_reso_data_meanshift[col[:-1]].Eval(pt) * 1e-4
                dca_reso_data_val = dca_reso_data[col[:-1]].Eval(pt)
                dca_reso_mc_val = dca_reso_mc[col[:-1]].Eval(pt)
            else:
                dca_reso_data_mean_shift = dca_reso_data_meanshift[col[:-1]].Eval(7) * 1e-4
                dca_reso_data_val = dca_reso_data[col[:-1]].Eval(7)
                dca_reso_mc_val = dca_reso_mc[col[:-1]].Eval(7)

            smear_value = gRandom.Gaus(dca + dca_reso_data_mean_shift,
                                       np.sqrt(dca_reso_data_val**2 - dca_reso_mc_val**2) * 1e-4)
            smear_values.append(smear_value)

        input_df[f"{dca_col}_SMEAR"] = smear_values

    # Close the input files
    for file_data, file_mc, file_mean in zip(input_files_data, input_files_mc, input_files_mean):
        file_data.Close()
        file_mc.Close()
        file_mean.Close()

    # Make a figure comparing the DCA variables before and after smearing
    num_cols = len(smear_cols)
    num_rows = num_cols // 2 + num_cols % 2
    fig, axs = plt.subplots(num_rows, 2, figsize=(10, 8))
    axs = axs.flatten()

    for i, col in enumerate(smear_cols):
        dca_col = f"fDCAPrim{col}"
        smear_col = f"{dca_col}_SMEAR"
        axs[i].hist(input_df[dca_col], bins=500, alpha=0.5, label="Before Smearing")
        axs[i].hist(input_df[smear_col], bins=500, alpha=0.3, label="After Smearing")
        axs[i].set_xlabel(col)
        axs[i].set_xlim(-0.05, 0.05)
        axs[i].set_ylim(10e2, None)
        axs[i].set_yscale('log')
        axs[i].legend()

    plt.tight_layout()
    fig.savefig(f"dca_comparison_{n_prongs}prong.png")
    plt.close(fig)

    return input_df

def separate_pt_bins(input_df, pt_bins):
    """
    Method to separate data in pT bins

    Parameters
    -----------------
    - input_df: pandas dataframe containing all candidates
    - pt_bins: list of ints containing pt bins

    Outputs
    -----------------
    - output_df_list: list of dataframes, each df containing a different pT region
    """

    print("Start separating in pT bins")
    output_df_list = []
    for i, pt in enumerate(pt_bins):
        if i-1 >= 0:
            df = input_df.query(f"fPt > {pt_bins[i-1]} and fPt < {pt_bins[i]}")
            output_df_list.append(df)
            print(f"{len(df)} events in range {pt_bins[i-1]} < pt < {pt_bins[i]}")
    return output_df_list

def divide_df_for_origin(input_df, cols_to_remove=None, channel=None):
    """
    Method to divide a dataframe in two (signal, bkg)

    Parameters
    -----------------
    - input_df: pandas dataframe containing all candidates with fFlagOrigin column
    - cols_to_remove: columns to be removed from output dataframes
    - channel: integer corresponding a specific fChannel for signal

    Outputs
    -----------------
    - df_signal: pandas dataframe containing only signal
    - df_bkg: pandas dataframe containing only background candidates
    """

    if cols_to_remove is None:
        cols_to_remove = ['fOriginMcRec']

    df_signal = input_df.query("fOriginMcRec == 1 or fOriginMcRec == 2")
    print(len(df_signal))

    if channel is not None and "fChannel" in df_signal.columns:
        df_signal = df_signal.query(f"fChannel == {channel}")

    df_bkg = input_df.query("fOriginMcRec == 0")
    print(len(df_bkg))

    cols_to_keep = list(df_signal.columns)
    for col in cols_to_remove:
        cols_to_keep.remove(col)
    df_signal = df_signal[cols_to_keep]

    df_bkg = df_bkg[cols_to_keep]

    return df_signal, df_bkg


# pylint: disable=too-many-locals,too-many-branches
# pylint: disable=too-many-nested-blocks,too-many-statements
def main(input_dir, max_files=1000, downscale_bkg=1., force=False, separate_pt=False, is_mc=False, do_smearing=False):
    """
    Main function

    Parameters
    -----------------
    - config: dictionary with configs
    """
    print("Searching for files")
    input_files = []
    for subdir in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, subdir)):
            for subsubdir in os.listdir(os.path.join(input_dir, subdir)):
                if os.path.isdir(os.path.join(input_dir, subdir, subsubdir)):
                    for file in os.listdir(os.path.join(input_dir, subdir, subsubdir)):
                        if "AO2D.root" in file:
                            input_files.append(os.path.join(
                                input_dir, subdir, subsubdir, file))
                        elif os.path.isdir(os.path.join(input_dir, subdir, subsubdir, file)):
                            for file2 in os.listdir(os.path.join(
                                    input_dir, subdir, subsubdir, file)):
                                if "AO2D.root" in file2:
                                    input_files.append(os.path.join(
                                        input_dir, subdir, subsubdir, file, file2))
    print(f"\033[32mFound {len(input_files)} AO2D.root files:\033[0m")
    for f in input_files:
        print(f"\t{f}")

    df_2p = None
    with alive_bar(len(input_files[:max_files])) as bar_alive:
        for file in input_files[:max_files]:
            print(f"\033[32mExtracting dataframes from input "
                  f"{file}\033[0m")

            file_root = uproot.open(file)
            indir = os.path.split(file)[0]

            # 2-prongs --> only D0
            is_d0_filtered = False
            for exfile in os.listdir(indir):
                if "D0ToKPi.parquet.gz" in exfile:
                    is_d0_filtered = True
                    break
            if not is_d0_filtered or force:
                print(file_root.keys())

                # Load all trees separately
                df_base, df_ml, df_par, df_sel, df_mc = None, None, None, None, None

                for tree_name in file_root.keys():
                    if df_base is None and df_ml is None and df_par is None and df_sel is None:
                        if "O2hfd0base" in tree_name:
                            df_base = file_root[tree_name].arrays(library="pd")
                        if "O2hfd0ml" in tree_name:
                            df_ml = file_root[tree_name].arrays(library="pd")
                        if "O2hfd0par" in tree_name:
                            df_par = file_root[tree_name].arrays(library="pd")
                        if "O2hfd0sel" in tree_name:
                            df_sel = file_root[tree_name].arrays(library="pd")
                    else:
                        if "O2hfd0base" in tree_name:
                            df_base = pd.concat([df_base, file_root[tree_name].arrays(library="pd")], ignore_index=True)
                        if "O2hfd0ml" in tree_name:
                            df_ml = pd.concat([df_ml, file_root[tree_name].arrays(library="pd")], ignore_index=True)
                        if "O2hfd0par" in tree_name:
                            df_par = pd.concat([df_par, file_root[tree_name].arrays(library="pd")], ignore_index=True)
                        if "O2hfd0sel" in tree_name:
                            df_sel = pd.concat([df_sel, file_root[tree_name].arrays(library="pd")], ignore_index=True)

                    if is_mc and "O2hfd0mc" in tree_name:
                        if df_mc is None:
                            df_mc = file_root[tree_name].arrays(library="pd")
                        else:
                            df_mc = pd.concat([df_mc, file_root[tree_name].arrays(library="pd")], ignore_index=True)

                print(f"Total candidates in df_base: {len(df_base)}")
                print(f"Total candidates in df_ml: {len(df_ml)}")
                print(f"Total candidates in df_par: {len(df_par)}")
                print(f"Total candidates in df_sel: {len(df_sel)}")
                print(f"Total candidates in df_sel: {len(df_mc)}")


                # Ensure all trees were found
                if is_mc:
                    if df_base is not None and df_ml is not None and df_par is not None and df_sel is not None and df_mc is not None:
                        temp_df = pd.concat([df_base, df_ml, df_par, df_sel, df_mc], axis=1)
                        # Merge the dataframes column-wise (entry-wise)
                        if df_2p is None:
                            df_2p = temp_df
                        else:
                            df_2p = pd.concat([df_2p, temp_df], ignore_index=True)
                    else:
                        print("Error: One or more trees were not found in the file.")
                        continue
                else:
                    if df_base is not None and df_ml is not None and df_par is not None and df_sel is not None:
                        temp_df = pd.concat([df_base, df_ml, df_par, df_sel], axis=1)
                        # Merge the dataframes column-wise (entry-wise)
                        if df_2p is None:
                            df_2p = temp_df
                        else:
                            df_2p = pd.concat([df_2p, temp_df], ignore_index=True)

                        #print(df_2p.head())
                    else:
                        print("Error: One or more trees were not found in the file.")
                        continue
                # Cut on y, -0.5 < y < 0.5
                df_2p = df_2p.query(f"fY > {-0.5} and fY < {0.5}")
                print(f"y cut performed between -0.5 and 0.5")
                # cut on cpa
                df_2p = df_2p.query(f"fCpa > {0.8}")
                print(f"cpa cut performed between to remove events with less than 0.8")
                # cut to select signal (for MC only)
                if is_mc:
                    df_2p = df_2p.query(f"(fCandidateSelFlag == {1} and fFlagMc == {1}) or (fCandidateSelFlag == {2} and fFlagMc == {-1})")
                    print(f"Signal candidates selected")
                if do_smearing:
                    df_2p = do_dca_smearing(df_2p, 2)
                if separate_pt:
                    pt_bins = [1, 2, 3, 5, 7, 10, 24, 50]
                    list_df = separate_pt_bins(df_2p, pt_bins)
                else:
                    list_df = [df_2p]
                for i, df in enumerate(list_df):
                    print("Saving dataframe")
                    if is_mc:
                        df_2p_signal, df_2p_bkg = divide_df_for_origin(df)
                        df_2p_bkg = df_2p_bkg.sample(frac=downscale_bkg, random_state=42)
                        if len(list_df) != 1:
                            #index = list_df.index(df)
                            print(i)
                            df.to_parquet(
                                os.path.join(indir, f"pt{i}/preparation_results/All_D0ToKPi.parquet.gz"),
                                compression="gzip"
                            )
                            df_2p_signal.to_parquet(
                                os.path.join(indir, f"pt{i}/preparation_results/signal/Signal_D0ToKPi.parquet.gz"),
                                compression="gzip"
                            )
                            df_2p_bkg.to_parquet(
                                os.path.join(indir, f"pt{i}/preparation_results/bkg/Bkg_D0ToKPi.parquet.gz"),
                                compression="gzip"
                            )
                        else:
                            df.to_parquet(
                                os.path.join(indir, f"all_pt/preparation_results/All_D0ToKPi.parquet.gz"),
                                compression="gzip"
                            )
                            df_2p_signal.to_parquet(
                                os.path.join(indir, f"all_pt/preparation_results/signal/Signal_D0ToKPi.parquet.gz"),
                                compression="gzip"
                            )
                            df_2p_bkg.to_parquet(
                                os.path.join(indir, f"all_pt/preparation_results/bkg/Bkg_D0ToKPi.parquet.gz"),
                                compression="gzip"
                            )
                    else:
                        if len(list_df) != 1:
                            #index = list_df.index(df)
                            print(i)
                            df.to_parquet(
                                os.path.join(indir, f"pt{i}/preparation_results/bkg/Bkg_D0ToKPi.parquet.gz"),
                                compression="gzip"
                            )
                        else:
                            df.to_parquet(
                                os.path.join(indir, f"all_pt/preparation_results/bkg/Bkg_D0ToKPi.parquet.gz"),
                                compression="gzip"
                            )
                    df = None

            bar_alive()

# Esta función sirve para que "main" se ejecute primero con los argumentos que especifiquemos
# Hace falta porque Python por defecto ejecuta las funciones en el orden en el que estén escritas
if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Arguments")
    PARSER.add_argument("input_dir", metavar="text", default="./",
                        help="input directory with AO2D input files")
    PARSER.add_argument("--max_files", type=int, default=1000,
                        help="max input files to be processed")
    PARSER.add_argument("--downscale_bkg", type=float, default=1.,
                        help="fraction of bkg to be kept")
    PARSER.add_argument("--force", action="store_true", default=False,
                        help="force re-creation of output files")
    PARSER.add_argument("--separatept", action="store_true", default=False,
                        help="separate in pt bins ")
    PARSER.add_argument("--is_mc", action="store_true", default=False,
                        help="is mc")
    PARSER.add_argument("--dosmearing", action="store_true", default=False,
                        help="do smearing on the dca of daughter tracks ")
    ARGS = PARSER.parse_args()

    main(ARGS.input_dir, ARGS.max_files, ARGS.downscale_bkg, ARGS.force, ARGS.separatept, ARGS.is_mc)
