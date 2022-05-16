import matplotlib.pyplot as plt
import pandas as pd

log_fdr = "../DataAnalysis/log/rt/"  # relative to the current working directory
log_filename = "log_data_1652718121.dat"
print(log_fdr+log_filename)
headers = ['t', 'x', 'y', 'z', 'roll', 'pitch', 'yaw',
           'vx', 'vy', 'vz', 'vroll', 'vpitch', 'vyaw']

df = pd.read_csv(log_fdr+log_filename, delimiter='\t',
                 names=headers, lineterminator='\r', skiprows=1)


# region Plotting
# plotting
plt.figure(figsize=(10, 5))
plt.plot(df['t'], df['x'], color='blue', label='x')
plt.plot(df['t'], df['y'], color='green', label='y')
plt.plot(df['t'], df['z'], color='red', label='z')
plt.grid()
plt.title('Position')
plt.xlabel('timestamp [s]')
plt.ylabel('Position [m]')
plt.legend()
plt.tight_layout()


plt.figure(figsize=(10, 5))
plt.plot(df['t'], df['roll'], color='blue', label='roll')
plt.plot(df['t'], df['pitch'], color='green', label='pitch')
plt.plot(df['t'], df['yaw'], color='red', label='yaw')
plt.grid()
plt.title('Angle')
plt.xlabel('timestamp [s]')
plt.ylabel('Angle [rad]')
plt.legend()
plt.tight_layout()


plt.figure(figsize=(10, 5))
plt.plot(df['t'], df['vx'], color='blue', label='vx')
plt.plot(df['t'], df['vy'], color='green', label='vy')
plt.plot(df['t'], df['vz'], color='red', label='vz')
plt.grid()
plt.title('Velocity')
plt.xlabel('timestamp [s]')
plt.ylabel('Velocity [m/s]')
plt.legend()
plt.tight_layout()

plt.figure(figsize=(10, 5))
plt.plot(df['t'], df['vroll'], color='blue', label='vroll')
plt.plot(df['t'], df['vpitch'], color='green', label='vpitch')
plt.plot(df['t'], df['vyaw'], color='red', label='vyaw')
plt.grid()
plt.title('Angular Velocity')
plt.xlabel('timestamp [s]')
plt.ylabel('Angular Velocity [rad/s]')
plt.legend()
plt.tight_layout()
plt.show()
# endregion
