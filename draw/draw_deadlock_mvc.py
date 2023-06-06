import matplotlib.pyplot as plt
from parse_deadlock_log import parse_deadlock_2vc3_log
import numpy as np
import matplotlib.ticker as mtick
from itertools import product


def draw_deadlock_mvc(log_file_list: list, figure_file: str, idx):

    cycle_list, apdl_list = parse_deadlock_2vc3_log(log_file_list)
    cycle_list = cycle_list[ :, :idx,:]
    apdl_list = apdl_list[ :, :idx,:]
    
    x = np.linspace(0, idx, idx)
    
    x_ticks = [20, 40, 60, 80, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800]
    x_ticks = [str(i) for i in x_ticks]

    font = {
        'family': 'serif',
        'weight': 'normal',
        'size': 9,
    }
    with plt.style.context(['science', 'ieee']):
        plt.rcParams['font.size'] = 6
        fig, ax = plt.subplots()
        ########################################################
        # ax.plot(x, cycle_list[0, :, 0], label='Same, D=4', marker='s', ms=3, color='k', linestyle='solid') 
        ax.plot(x, cycle_list[1, :, 0], label='Random, D=4', marker='o', ms=3, color='r', linestyle='solid')
        ax.plot(x, cycle_list[2, :, 0], label='Optimal, D=4', marker='o', ms=3, color='b', linestyle='dotted')
        
        # ax.plot(x, cycle_list[0, :, 1], label='Same, D=6', marker='^', ms=3, color='k', linestyle='solid')
        ax.plot(x, cycle_list[1, :, 1], label='Random, D=6', marker='h', ms=3, color='r', linestyle='solid')
        ax.plot(x, cycle_list[2, :, 1], label='Optimal, D=6', marker='h', ms=3, color='b', linestyle='dotted')
        
        # ax.plot(x, cycle_list[0, :, 2], label='Same, D=8', marker='d', ms=3, color='k', linestyle='solid')
        ax.plot(x, cycle_list[1, :, 2], label='Random, D=8', marker='p', ms=3, color='r', linestyle='solid')
        ax.plot(x, cycle_list[2, :, 2], label='Optimal, D=8', marker='p', ms=3, color='b', linestyle='dotted')
        
        ax.legend(prop={'size': 6}, loc=0)
        ax.set_xlabel('Network Size', fontdict=font)
        ax.set_ylabel(ylabel='Latency/cycles', fontdict=font)
        ax.set_xticks(x)
        ax.set_xticklabels(x_ticks[:idx], rotation = 40, fontsize = 'small')
        ax.tick_params(labelsize=6)
        fig.savefig(figure_file + 'multi_root_2vc3_cycle.pdf', dpi=5000)
        
        #########################################################
        plt.rcParams['font.size'] = 6
        fig, ax = plt.subplots()
        # ax.plot(x, apdl_list[0, :, 0], label='Same, D=4', marker='s', ms=3, color='k', linestyle='solid') 
        ax.plot(x, apdl_list[1, :, 0], label='Random, D=4', marker='o', ms=3, color='r', linestyle='solid')
        ax.plot(x, apdl_list[2, :, 0], label='Optimal, D=4', marker='o', ms=3, color='b', linestyle='dotted')
        
        # ax.plot(x, apdl_list[0, :, 1], label='Same, D=6', marker='^', ms=3, color='k', linestyle='solid')
        ax.plot(x, apdl_list[1, :, 1], label='Random, D=6', marker='h', ms=3, color='r', linestyle='solid')
        ax.plot(x, apdl_list[2, :, 1], label='Optimal, D=6', marker='h', ms=3, color='b', linestyle='dotted')
        
        # ax.plot(x, apdl_list[0, :, 2], label='Same, D=8', marker='d', ms=3, color='k', linestyle='solid')
        ax.plot(x, apdl_list[1, :, 2], label='Random, D=8', marker='p', ms=3, color='r', linestyle='solid')
        ax.plot(x, apdl_list[2, :, 2], label='Optimal, D=8', marker='p', ms=3, color='b', linestyle='dotted')
        
        ax.legend(prop={'size': 6}, loc=0)
        ax.set_xlabel('Network Size', fontdict=font)
        ax.set_ylabel(ylabel='Average Packet Delivery Latency/cycles', fontdict=font)
        ax.set_xticks(x)
        ax.set_xticklabels(x_ticks[:idx], rotation = 40, fontsize = 'small')
        ax.tick_params(labelsize=6)
        fig.savefig(figure_file + 'multi_root_2vc3_apdl.pdf', dpi=5000)


if __name__ == '__main__':
    
    idx = 14

    file_list = ['temp/deadlock_free/2vc3_log/{:d}.log'.format(i) for i in range(idx * 3 * 3)]
    draw_deadlock_mvc(file_list, 'temp/deadlock_free/figures/', idx)
    
