import numpy as np

def make_table(components,whb, whn,whv, bs,ns,vs, berr,nerr,verr):
# create and empty array of the right shape for our table that allows strings as an input
    table_data = np.asarray([['                                           ' for x in range(7)] for y in range(components)])
    col_labels = ['component', 'bWH', 'NWH', 'v_radWH', 'b', 'N', 'v_rad']

    n = 0
    # fill table data with all necessary data that we want on the table
    for x in range(7):
        if x == 0:
            for y in range(components):
                table_data[y][x] = f'{n}'
                n += 1
        elif x == 1:
            for y in range(components):
                table_data[y][x] = f'{whb[y]:#.3f}'
        elif x == 2:
            for y in range(components):
                table_data[y][x] = f'{whn[y]:#.3g}'
        elif x == 3:
            for y in range(components):
                table_data[y][x] = f'{whv[y]:#.3f}'
                # col_labels[x] = f'b{x}'
                # col_labels[x + components] = f'N{x}'
                # col_labels[x + components * 2] = f'v_rad{x}'
        elif x == 4:
            for y in range(components):
                table_data[y][x] = '{0:#.3f} \u00B1 {1:#.3f}'.format(bs[y], berr[y])
        elif x == 5:
            for y in range(components):
                table_data[y][x] = '{0:#.3g} \u00B1 {1:#.3g}'.format(ns[y], nerr[y])
        else:
            for y in range(components):
                table_data[y][x] = '{0:#.3f} \u00B1 {1:#.3f}'.format(vs[y], verr[y])



    return table_data, col_labels