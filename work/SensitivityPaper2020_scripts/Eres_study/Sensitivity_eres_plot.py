import pandas
import os
from matplotlib import pyplot as plt
import numpy as np
import histlite as hl

dates = ["20_11_19"]
# dates = ["20_10_22", "20_10_21", "20_10_19"]
# root_dir = '/p/lustre2/nexouser/czyz1/output/'
root_dir = '/Users/czyz1/OneDrive - LLNL/Projects/nEXO/testoutput/'

num_dataset = 10
start_it = 0
end_it = 500
eres = ['0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016']

def calc_atoms_136():
    """ Number of Xe136 atoms in nEXO fiducial volume """
    mmass134 = 0.133905395  # kg/mol 134
    mmass136 = 0.135907219  # kg/mol 136
    at_frac = 0.9           # atomic fraction 136 / (136 + 134)
    avog_num = 6.022141E23  # Avogadro's number
    fid_mass = 3281         # mass of fiducial volume [kg]

    atoms136 = (fid_mass * avog_num * at_frac) / ((mmass136 * at_frac) + ((1 - at_frac) * mmass134))

    return atoms136


def sensitivity_calc(atoms136, lt_years, cross_median):
    """Calculate the sensitivty of nEXO in terms of half-life (years)"""
    eff = 0.963  # hit efficiency
    sensitivity = eff * atoms136 * lt_years * np.log(2) / cross_median

    return sensitivity


if __name__ == "__main__":

    lt_years = 10
    atoms136 = calc_atoms_136()
    fig2, ax2 = plt.subplots()

    colors = ['b', 'r', 'g', 'k', 'm']

    for date in dates:
        if date == '20_11_04_DNN1_024':     # special case where the code was ran differently
            end_its = 450
            num_datasets = 50
        else:
            end_its = end_it
            num_datasets = num_dataset

        num_its = end_its - start_it
        num_toys = num_its * len(eres) * num_datasets

        fig, ax = plt.subplots()
        dnn = 'DNN1'
        materialdb = "024"
        failed = 0
        converged = 0
        sensitivity = []

        input_dir = root_dir + "h5/" + date
        output_dir = root_dir + "plots/"# + date

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for plot_index, res in enumerate(eres):
            num_runs = 0
            crossing_masked = []


            for iter in range(start_it, end_its):
                num_runs += 1
                filename = '{}/sens_output_file_90CL_{}_resolution_{}.h5'.format(input_dir, iter, res)
                if not os.path.exists(filename):
                    num_toys -= num_datasets
                    failed += 1
                    print('livetime = {}, iteration = {}, failed = {}'.format(res, iter, failed))
                    print(filename)
                    continue
                df = pandas.read_hdf(filename)
                crossing_masked = crossing_masked + [b for a, b in zip(df['best_fit_converged'],
                                                                       df['90CL_crossing']) if (a and b > 0)]
            converged += len(crossing_masked)
            sensitivity.append(sensitivity_calc(atoms136, lt_years, np.median(crossing_masked)))
            ax.axvline(np.median(crossing_masked), color=colors[plot_index % len(colors)], linestyle='--')
            hteststats = hl.hist(crossing_masked, bins=np.linspace(0, 30, 31))
            hl.plot1d(ax, hteststats,  color=colors[plot_index % len(colors)], label='E_res = {} Y, T_1/2 = {:.3e} Y,'
                                                                                     '\n Median = {:.3}'
                      .format(res, sensitivity[-1], np.median(crossing_masked)))
            ax.set_xlabel('Upper limit on 0nuBB counts at 90% confidence limit', fontsize=16)
            ax.set_ylabel('Counts ({} toys, {} converged)'.format(num_toys, converged), fontsize=16)
            ax.set_xlim([0, 30])
            ax.legend()

        fig.savefig('{}/sens_hist_{}_{}.png'.format(output_dir, dnn, materialdb), dpi=800)

        ax2.plot(eres, sensitivity, '-o', label='MatDB = {}, DNN = {}'.format(materialdb, dnn))
        ax2.set_xlim([0, 12])
        ax2.set_ylim([1E27, 2E28])
        ax2.set_yscale('log')
        ax2.legend()

    fig2.savefig('{}/sens_vs_hl_{}_{}.png'.format(output_dir, dnn, materialdb), dpi=800)

    fig.show()
    fig2.show()
