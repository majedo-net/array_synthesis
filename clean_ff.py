import os
import numpy as np

if __name__ == '__main__':
    for f in os.listdir('/data/pso_results/ff'):
        f = f'/data/pso_results/ff/{f}'
        try:
            test = np.genfromtxt(f,delimiter=',',max_rows=181)
        except ValueError:
            print(f'Error reading {f}, deleting far field')
            os.remove(f)


