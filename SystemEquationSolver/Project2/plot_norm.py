import matplotlib.pyplot as plt

# Read data from file
data = []
with open('data.txt', 'r') as file:
    next(file)  # skip header
    for line in file:
        values = line.strip().split('\t')
        if len(values) == 3:
            iteration, jacobi_norm, gauss_seidel_norm = map(float, values)
            data.append((iteration, jacobi_norm, gauss_seidel_norm))
        elif len(values) == 2:
            iteration, jacobi_norm = map(float, values)
            data.append((iteration, jacobi_norm, 0.0))  # Assign default value for Gauss-Seidel norm
        else:
            print("Invalid data format:", line)

# Separate data into lists for plotting
iterations = [d[0] for d in data]
jacobi_norms = [d[1] for d in data]
gauss_seidel_norms = [d[2] for d in data]

# Plot the graph
plt.figure(figsize=(10, 6))
plt.plot(iterations, jacobi_norms, label='Jacobi Norm', marker='o')
plt.plot(iterations, gauss_seidel_norms, label='Gauss-Seidel Norm', marker='x')
plt.yscale('log')  # Set y-axis scale to logarithmic
plt.xlabel('Iteration')
plt.ylabel('Norm')
plt.title('Norm of Residuals vs Iterations')
plt.legend()
plt.grid(True)
plt.show()
