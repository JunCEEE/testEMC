from lmfit.models import GaussianModel, LinearModel
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
import h5py


def set_plot_style(fig_w, fig_h):
    # Plotting Style Setup
    plt.style.use(['science', 'notebook'])
    # fig_w = 3.5
    # fig_h = 2.625
    # fig_w = 3.3
    # fig_h = 2.7

    size_scale = np.sqrt(fig_w * fig_h / (3.5 * 2.625))
    lablesize = 16 * size_scale
    plt.rcParams.update({
        "axes.labelsize": lablesize,
        "legend.fontsize": 0.8 * lablesize,
        "xtick.labelsize": 0.8 * lablesize,
        "ytick.labelsize": 0.8 * lablesize,
        "legend.title_fontsize": lablesize,
        "axes.titlesize": lablesize,
    })  # specify font size here

    # my_cycle = cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']) * cycler(linestyle=['-', '--', ':', '-.'])
    # my_cycle = cycler('color',
    #                   ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 'gold'
    #                    ]) + cycler(linestyle=['-', ':', '--', '-.', '--'])
    my_cycle = cycler('color',
                      ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 'gold'])

    plt.rcParams["axes.prop_cycle"] = my_cycle


def fitting(x, y):
    offset = LinearModel(prefix='offset_')
    pars = offset.make_params()
    # pars = offset.guess(y, x=x)

    g1 = GaussianModel(prefix='g1_')
    pars.update(g1.make_params())
    pars['g1_center'].set(value=0, min=-5, max=5)
    pars['g1_sigma'].set(value=8)
    pars['g1_amplitude'].set(value=2054)
    # pars['g1_height'].set(value=80, min=75, max=85)

    g2 = GaussianModel(prefix='g2_')
    pars.update(g2.make_params())
    pars['g2_center'].set(value=55, min=25, max=80)
    pars['g2_sigma'].set(value=20)
    pars['g2_amplitude'].set(value=1583)
    # pars['g2_height'].set(value=50, min=40, max=58)

    g3 = GaussianModel(prefix='g3_')
    pars.update(g3.make_params())
    pars['g3_center'].set(value=120, min=100, max=150)
    pars['g3_sigma'].set(value=20)
    pars['g3_amplitude'].set(value=500)
    # pars['g3_height'].set(value=15, min=10)

    g4 = GaussianModel(prefix='g4_')
    pars.update(g4.make_params())
    pars['g4_center'].set(value=160, min=150, max=200)
    pars['g4_sigma'].set(value=20)
    pars['g4_amplitude'].set(value=130)
    # pars['g4_height'].set(value=5, min=0, max=10)

    # mod = offset + g1 + g2 + g3 + g4
    mod = g1 + g2 + g3 + g4
    init = mod.eval(pars, x=x)
    out = mod.fit(y, pars, x=x)
    print(out.fit_report(min_correl=0.5))

    comps = out.eval_components(x=x)

    return init, out.best_fit, comps, out


x = np.loadtxt('./r0112_hist_x')
x = x[:300]
y = np.loadtxt('./r0112_hist_y')

init, best, comps, out = fitting(x, y)

# fitting plot
# fig_w = 16
# fig_h = 8
# # set_plot_style(fig_w, fig_h)

# fig, ax = plt.subplots(1, 2, figsize=(fig_w, fig_h))
# ax[0].plot(x, y, 'o')
# ax[0].plot(x, init, '--', label='init')
# ax[0].plot(x, best, '-', label='best')
# ax[0].legend()

# ax[1].plot(x, y, 'o')
# # ax[1].plot(x, comps['offset_'], '--', label='offset')
# ax[1].plot(x, comps['g1_'], '--', label='g1')
# ax[1].plot(x, comps['g2_'], label='g2')
# ax[1].plot(x, comps['g3_'], label='g3')
# ax[1].plot(x, comps['g4_'], label='g4')
# ax[1].legend()
# # plt.show()
# plt.savefig('photon_peaks_fitting-3.pdf', dpi=300)

# fitting plot
fig_w = 6
fig_h = 4.5
set_plot_style(fig_w, fig_h)

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.plot(x, y, 'o')
# ax.plot(x, best, '--', label='fitting', color='b')
# ax.legend()

ax.plot(x, comps['g1_'], label=r'$\sigma_1$ = {:.1f}'.format(out.params['g1_sigma'].value))
ax.plot(x, comps['g2_'], label=r'$\sigma_2$ = {:.1f}'.format(out.params['g2_sigma'].value))
ax.plot(x, comps['g3_'], label=r'$\sigma_3$ = {:.1f}'.format(out.params['g3_sigma'].value))
ax.plot(x, comps['g4_'], label=r'$\sigma_4$ = {:.1f}'.format(out.params['g4_sigma'].value))
plt.legend()
plt.ylim(bottom=0)
plt.xlabel('ADU')
plt.ylabel('Number of pixels')
plt.savefig('photon_peaks_line.pdf', dpi=300)
plt.show()
