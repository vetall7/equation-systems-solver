import matplotlib.pyplot as plt
import numpy as np

def plot_times(data_file, sizes):
    data = np.loadtxt(data_file, skiprows=1)
    times_jacobi = data[:, 1]
    times_gauss_seidel = data[:, 2]
    times_lu = data[:, 3]

    plt.figure(figsize=(10, 6))
    plt.plot(sizes, times_jacobi, label='Jacobi')
    plt.plot(sizes, times_gauss_seidel, label='Gauss-Seidel')
    plt.plot(sizes, times_lu, label='LU Factorization')
    plt.xlabel('Array Size')
    plt.ylabel('Time (ms)')
    plt.title('Execution Time of Different Methods')
    plt.xticks(sizes)
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    sizes = [100, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000]  # Update with your array sizes
    plot_times("data.txt", sizes)
