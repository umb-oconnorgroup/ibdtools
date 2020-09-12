import numpy as np
import matplotlib.pyplot as plt


def qqplot(x, y, label="label"):
    """x and y are two np.arrays of the same size"""
    x = np.sort(x)
    y = np.sort(y)
    max_xy = np.max([x, y])
    min_xy = np.min([x, y])

    plt.figure(figsize=(4, 4))
    plt.scatter(x, y)
    plt.plot([min_xy, max_xy], [min_xy, max_xy], 'r')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(label + '.pdf')


def test_qqplot():
    x = np.random.random(size=1000)
    y = np.random.random(size=1000)
    qqplot(x, y, 'figure')
