import numpy as np
import matplotlib.pyplot as plt
from functions.benchmark_voigt import voigt_profile as vp
def populate_params(ourfit,bour,dbic,bwh,nwh,vwh,wl,nf,cw,f,g,vres,nstep,stddev, starname,whfit,errors):
    b_WH = bwh
    N_WH = nwh
    v_rad_WH = vwh
    fit = ourfit

    # stating the number of components, this seemed like the easiest way to calculate it
    components = (len(bour) - 1)
    b = np.zeros(components)
    b_err = np.zeros(components)
    N = np.zeros(components)
    N_err = np.zeros(components)
    v_rad = np.zeros(components)
    v_rad_err = np.zeros(components)

    for j in range(components):
        b[j] = fit.params[f'b{j}'].value
        b_err[j] = fit.params[f'b{j}'].stderr
        N[j] = fit.params[f'N{j}'].value
        N_err[j] = fit.params[f'N{j}'].stderr
        v_rad[j] = fit.params[f'v_rad{j}'].value
        v_rad_err[j] = fit.params[f'v_rad{j}'].stderr

    print(b)
    print(N)
    print(v_rad)
    # print results to see if we are satisfied
    print('final no components:', components)
    print('reduced chi-square: ', fit.summary()['redchi'])
    # print('bic', fit.summary()['bic'])
    print('delta bic', dbic)
    # fit.params.pretty_print()

    # plot final fit alongside the raw data
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    # ax = fig.add_subplot(111)
    ax2.errorbar(wl, nf, yerr=np.array([stddev] * len(wl)), mfc='none', c='k', label='Data',
                 alpha=0.2)
    ax2.plot(wl, fit.best_fit, c='b', label='Fit')
    ax1.plot(wl, nf / fit.best_fit, label='residuals')
    # ax1.set_ylim(-0.1, 1.5)
    ax1.set_xlim(4226.8, 4227)
    # ax2.set_ylim(-0.1, 1.5)
    ax2.set_xlim(4226.8, 4227)
    # ax3.set_ylim(-0.1, 1.5)
    ax3.set_xlim(4226.8, 4227)
    fig.suptitle(starname)
    ax3.plot(wl, nf - fit.best_fit, label='actual residuals')
    # plt.legend()
    # plt.show()
    ax2.plot(wl, whfit.best_fit, label='welty and hobbs fit', c = 'orange', alpha = 0.6)
    ax1.legend(loc="upper right")
    ax2.legend(loc="upper right")
    ax3.legend(loc="upper right")
    plt.show()
    #do np.append for b, n and vrad
    #bnew = np.append(b,0.45)
    #Nnew = np.append(N, 2e10)
    #v_radnew = np.append(v_rad,-10)

    #fittedmodelman = vp.fit_multi_voigt_absorptionlines(wavegrid=wl, ydata=nf,
                                                        #restwave=cw, f=f,
                                                        #gamma=g, b=bnew, N=Nnew, v_rad=v_radnew, v_resolution=vres,
                                                        #n_step=nstep, std_dev=errors)
    # plot statement
    #plt.plot(wl, nf, label='data')
    #plt.plot(wl, fittedmodelman.best_fit, label='fitted model')
    #print(fittedmodelman.redchi)  # matches within 3 decimal places of chi squared given in the fit_data txt
    #print(fittedmodelman.bic)
    #print(fittedmodelman.aic)
    #plt.legend()
    #plt.show()

    #print(bnew)
    #print(Nnew)
    #print(v_radnew)
    # get parameters from Welty and Hobbs paper for comparison
    # folder = Path(PYTHONDIR + "/data")
    # filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / (star_name[i] + '.txt')

    # Headers = ["b", "N", "v_rad"]
    # star_parameters = pd.read_csv(filename,
    # delim_whitespace=True,
    # header=None,
    # names=Headers,
    # engine="python",
    # )

    # collect the data from each header in easy to iterate form
    # b_WH = np.asarray(star_parameters["b"])
    # N_WH = np.asarray(star_parameters["N"])
    # v_rad_WH = np.asarray(star_parameters["v_rad"])

    # this makes sure that the arrays are the same size as the fitted array
    # just to allow for comparisons whilst bug fixing
    while (components - len(b_WH)) > 0:
        b_WH = np.append(b_WH, 0)
        N_WH = np.append(N_WH, 0)
        v_rad_WH = np.append(v_rad_WH, 0)
    length = len(b_WH) - components
    while length > 0:
        # print(length)
        new = len(b_WH)
        fit.params.add(f'b{new}', value=0, vary=False)
        fit.params.add(f'N{new}', value=0, vary=False)
        fit.params.add(f'v_rad{new}', value=0, vary=False)
        length -= 1

    # set the v_rad into increasing order and make sure the associated b and N values
    # move with them to the correct positions in their respective arrays
    v_rad_order = v_rad
    v_rad_err_order = v_rad_err
    b_order = b
    b_err_order = b_err
    N_order = N
    N_err_order = N_err
    v_rad.sort()

    # reorder the arrays
    for j in range(components):
        b[np.argwhere(v_rad == v_rad_order[j])[0, 0]] = b_order[j]
        b_err[np.argwhere(v_rad == v_rad_order[j])[0, 0]] = b_err_order[j]
        N[np.argwhere(v_rad == v_rad_order[j])[0, 0]] = N_order[j]
        N_err[np.argwhere(v_rad == v_rad_order[j])[0, 0]] = N_err_order[j]
        v_rad_err[np.argwhere(v_rad == v_rad_order[j])[0, 0]] = v_rad_err_order[j]


    # create and empty array of the right shape for our table that allows strings as an input
    table_data = np.asarray(
        [['                                           ' for x in range(7)] for y in range(len(b_WH))])
    col_labels = ['component', 'bWH', 'NWH', 'v_radWH', 'b', 'N', 'v_rad']
    # row_labels = ['Welty & Hobbs', 'Fitting components']
    n = 0
    # fill table data with all necessary data that we want on the table
    for x in range(7):
        if x == 0:
            for y in range(len(b_WH)):
                table_data[y][x] = f'{n}'
                n += 1
        elif x == 1:
            for y in range(len(b_WH)):
                table_data[y][x] = f'{b_WH[y]:#.3f}'
        elif x == 2:
            for y in range(len(b_WH)):
                table_data[y][x] = f'{N_WH[y]:#.3g}'
        elif x == 3:
            for y in range(len(b_WH)):
                table_data[y][x] = f'{v_rad_WH[y]:#.3f}'
                # col_labels[x] = f'b{x}'
                # col_labels[x + components] = f'N{x}'
                # col_labels[x + components * 2] = f'v_rad{x}'
        elif x == 4:
            for y in range(components):
                table_data[y][x] = '{0:#.3f} \u00B1 {1:#.3f}'.format(b[y], b_err[y])
        elif x == 5:
            for y in range(components):
                table_data[y][x] = '{0:#.3g} \u00B1 {1:#.3g}'.format(N[y], N_err[y])
        else:
            for y in range(components):
                table_data[y][x] = '{0:#.3f} \u00B1 {1:#.3f}'.format(v_rad[y], v_rad_err[y])

    # show table
    table_data2 = np.asarray(
        [['                                           ' for x in range(3)] for y in range(1)])
    col_labels2 = ['Reduced Chi Sqr', 'BIC', 'AIC']
    # row_labels = ['Welty & Hobbs', 'Fitting components']
    # fill table data with all necessary data that we want on the table
    for x in range(3):
        if x == 0:
            for y in range(1):
                table_data2[y][x] = f'{fit.summary()['redchi']:#.3g}'

        elif x == 1:
            for y in range(1):
                table_data2[y][x] = f'{fit.summary()['bic']:#.3g}'
        elif x == 2:
            for y in range(1):
                table_data2[y][x] = f'{fit.summary()['aic']:#.3g}'




    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    table = ax1.table(table_data, colLabels=col_labels, loc='center')
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    table2 = ax2.table(table_data2,colLabels=col_labels2, loc='center')
    # table.scale(1,1)
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table2.auto_set_font_size(False)
    table2.set_fontsize(8)
    ax1.set_frame_on(False)
    ax2.set_frame_on(False)
    fig.suptitle(starname)
    # fig_3.suptitle(star_name[i])
    # save table data as an image so that data can be viewed without running the code everytime
    # plt.savefig('c:/Users/user/edibles/edibles/data/voigt_benchmarkdata/parameter_modelling_data/'+star_name[i]+'_table_data.png')
    plt.show()

    return fit, b, N, v_rad





