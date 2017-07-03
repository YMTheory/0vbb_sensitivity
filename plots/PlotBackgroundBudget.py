#!/usr/bin/env python

"""PlotBackgroundBudget.py: create figures of background budget

Input is from DB Excel file. Plots are generated grouping by Isotope, Material, or Component."""

import os.path
import sys
from math import floor, log10

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def SqrtSumSq(x):
    x = np.array(x.tolist())
    return np.sqrt(np.sum(np.multiply(x, x)))


def make_plot(df, groupby, filename, xlimit=5.0e-5, ):
    # grouping with custom aggregration
    df2 = df.groupby(groupby + ["CV?"]).agg(
        {'C.V.': np.sum, 'Error': SqrtSumSq, 'Limit 90% C.L.': np.sum, 'CV?': all})

    # add a variable for sorting that considers both CV and 90CL appropriately
    df2['value'] = df2['C.V.'] * df2['CV?'] + df2['Limit 90% C.L.'] * (~df2['CV?'])
    df2.sort_values('value', ascending=True, inplace=True)
    print(df2)

    fig, ax0 = plt.subplots()

    # loop over all the rows to fill the lists used for plotting
    label = []
    xlimit_arrow = 10 ** floor(log10(xlimit)) * 10
    for index, row in df2.iterrows():
        if row['Limit 90% C.L.'] < 1e-5:
            continue

        label.append(' '.join(index[:len(groupby)]))
        label[-1] = label[-1].replace('bb2n', r'$\beta\beta 2\nu$')

        # the counts are in 1 year and 3 tonne so divide by 3 to get the cts/tonne/year
        if row['CV?']:
            xuplims = 0
            value = row['C.V.'] / 3.
            err = row['Error'] / 3.
            color = 'teal'
            marker = '.'
        else:
            xuplims = 1
            value = row['Limit 90% C.L.'] / 3.
            err = value - xlimit_arrow
            color = 'darkorange'
            marker = ''

        # plot one item
        ax0.errorbar(value, len(label) * 2 - 2, xerr=err, xuplims=xuplims,
                     lw=1, capsize=2, capthick=1, color=color, marker=marker)

    # overall plot formatting (labels, etc...)
    nn = len(label)
    y = np.arange(0, 2 * nn, 2)
    ax0.set_yticks(y)
    ax0.set_yticklabels(label)
    ax0.set_xscale('log')
    ax0.set_xlim(xlimit)
    ax0.set_xlabel("cts/ROI/tonne/year")
    plt.subplots_adjust(left=0.32, )
    # plt.show()

    plt.savefig(filename, dpi=200)
    return


def main(argv):
    if len(argv) >= 2:
        inTableName = argv[1]  # '../tables/Summary_v73_2016-09-26_bb2n_0nu.xlsx'
    else:
        sys.exit('Usage: %s xlsx_table_filename' % argv[0])

    if not os.path.exists(inTableName):
        sys.exit('ERROR: File %s was not found!' % inTableName)

    outFolder = os.path.curdir
    if len(argv) == 3:
        outFolder = os.path.expanduser(argv[2])

    table = ''.join(os.path.basename(inTableName).split('.')[:-1])

    df = pd.read_excel(inTableName, sheetname='SS_ExpectedCounts', header=0,
                       skiprows=4, skip_footer=7, parse_cols="A:C,AI:AK", )

    # this is for printing without breaking into multiple lines
    pd.set_option('display.expand_frame_repr', False)

    # remove bb0n rows
    df = df[~df['Isotope'].isin(['bb0n', ])]

    # Add a column that tracks whether the CV>0 or not
    df['CV?'] = pd.Series(df['C.V.'] > 0, index=df.index)
    # print(df)

    # Plot by Material, Isotope, and Component separately
    make_plot(df, ['Material'], os.path.join(outFolder, "BackgroundBudgetByMaterial_" + table + ".png"))
    make_plot(df, ['Isotope'], os.path.join(outFolder, "BackgroundBudgetByIsotope_" + table + ".png"), 5.0e-4)
    make_plot(df, ['Component'], os.path.join(outFolder, "BackgroundBudgetByComponent_" + table + ".png"))

    # This is by Material & Isotope
    make_plot(df, ["Material", "Isotope"], os.path.join(outFolder, "BackgroundBudgetByMaterialIsotope_" + table + '.png'))


if __name__ == "__main__":
    main(sys.argv)
