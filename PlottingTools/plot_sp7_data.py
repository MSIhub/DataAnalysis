import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PI = 3.14159265359

# region HelperFunctions


def read_data_from_log(log_path):
    headers = ['t', 'x', 'y', 'z', 'roll', 'pitch', 'yaw',
               'vx', 'vy', 'vz', 'vroll', 'vpitch', 'vyaw']
    df = pd.read_csv(log_path, delimiter='\t',
                     names=headers, lineterminator='\r', skiprows=1)
    return df


def plot_sp7_data(df):
    fig = plt.figure(figsize=(16, 8))
    gs = GridSpec(nrows=2, ncols=2)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.plot(df['t'], df['x'], color='blue', label='x')
    ax0.plot(df['t'], df['y'], color='green', label='y')
    ax0.plot(df['t'], df['z'], color='red', label='z')
    ax0.grid()
    ax0.set_title('Position')
    ax0.set_xlabel('timestamp [s]')
    ax0.set_ylabel('Position [m]')
    ax0.legend()

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.plot(df['t'], df['roll']*180/PI, color='blue', label='roll')
    ax1.plot(df['t'], df['pitch']*180/PI, color='green', label='pitch')
    ax1.plot(df['t'], df['yaw']*180/PI, color='red', label='yaw')
    ax1.grid()
    ax1.set_title('Angle')
    ax1.set_xlabel('timestamp [s]')
    ax1.set_ylabel('Angle [deg]')
    ax1.legend()

    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(df['t'], df['vx'], color='blue', label='vx')
    ax2.plot(df['t'], df['vy'], color='green', label='vy')
    ax2.plot(df['t'], df['vz'], color='red', label='vz')
    ax2.grid()
    ax2.set_title('Velocity')
    ax2.set_xlabel('timestamp [s]')
    ax2.set_ylabel('Velocity [m/s]')
    ax2.legend()

    ax3 = fig.add_subplot(gs[1, 1])
    ax3.plot(df['t'], df['vpitch']*180/PI, color='green', label='vpitch')
    ax3.plot(df['t'], df['vroll']*180/PI, color='blue', label='vroll')
    ax3.plot(df['t'], df['vyaw']*180/PI, color='red', label='vyaw')
    ax3.grid()
    ax3.set_title('Angular Velocity')
    ax3.set_xlabel('timestamp [s]')
    ax3.set_ylabel('Angular Velocity [deg/s]')
    ax3.legend()

    fig.align_labels()
    plt.tight_layout()
    plt.show()
# endregion


log_fdr = "../DataAnalysis/log/rt/"  # relative to the current working directory
log_filename = "log_data_1652963710.dat"
log_path = log_fdr+log_filename

df = read_data_from_log(log_path)
plot_sp7_data(df)
