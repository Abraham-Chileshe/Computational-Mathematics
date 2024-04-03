import math
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
from tkinter import Tk, Text, END

def f(x : float):
    return 0.5*x**2 + np.cos(2*x)

def lagrange_interpolation(xi :float, fi:float, x:float):
    n = len(xi)
    result = 0
    for i in range(n):
        term = fi[i]
        for j in range(n):
            if j != i:
                term *= (x - xi[j]) / (xi[i] - xi[j])
        result += term
    return result

# Абсолютная ошибка
def absolute_error(f:float, interpolant:float, x:float):
    return np.abs(interpolant - f(x))

# Относительная ошибка
def relative_error(abs_error, x):
    true_value = f(x)  # Calculate the true function value
    return (abs_error / abs(true_value[0])) * 100

# Оценка остаточного члена

def residual_estimate(f, n, a, b):
    # Step 1: Compute the maximum absolute value of the (n+1)-th derivative of f(x)
    x_vals = np.linspace(a, b, 1000)
    f_n_plus_1 = np.max(np.abs(np.diff(f(x_vals), n+1)))  # Maximum absolute difference of the (n+1)-th derivative

    # Step 2: Compute the factorial of (n+1)
    factorial_n_plus_1 = math.factorial(n+1)

    # Step 3: Compute (b - a) raised to the power of (n+1)
    b_minus_a_power_n_plus_1 = (b - a)**(n+1)

    # Step 4: Compute the final estimate
    r_n = f_n_plus_1 / (factorial_n_plus_1 * b_minus_a_power_n_plus_1)
    return r_n


# Функция для обновления окна с таблицей
def update_table(text_widget, data):
    table = tabulate(data, headers=["n ", "Абсолютная ошибка (Δfn)", "Относительная ошибка (δfn)", "теоретической ошибка (rn)"], tablefmt="grid", floatfmt=".15f")
    text_widget.delete('1.0', END)
    text_widget.insert(END, table)


# Основная функция
def main():
    root = Tk()
    root.title("Interpolation Table and Function Plot")

    # Create a text widget for displaying the table
    text_widget = Text(root, height=30, width=150)
    text_widget.pack()

    a = 0.6  # Start of the interval
    b = 1.1  # End of the interval
    n_values = [3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]  # n from 3 to 100 with step 2

    data = []

    # Plot the target function
    plt.figure(figsize=(6, 6))
    x_vals = np.linspace(a, b, 1000)
    plt.plot(x_vals, f(x_vals), label='Target Function', color='blue')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Target Function Plot')
    plt.grid(True)
    plt.legend()

    for n in n_values:
        xi = np.linspace(a, b, n)
        fi = f(xi)

        # Interpolation
        x_vals = np.linspace(a, b, 1000)
        interpolant = lagrange_interpolation(xi, fi, x_vals)

        # Absolute and relative errors
        abs_error = np.max(absolute_error(f, interpolant, x_vals))
        rel_error = np.max(relative_error(abs_error, x_vals))

        # Residual estimate
        r_n = residual_estimate(f, n, a, b)

        data.append([n, abs_error, rel_error, r_n])

    # Convert data to a numpy array for convenience
    data = np.array(data)

    # Update the table window
    update_table(text_widget, data)

    # Plot the absolute error, relative error, and residual estimate
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 3, 1)
    plt.plot(data[:, 0], data[:, 1], label='Δfn')
    plt.xlabel('n')
    plt.ylabel('Absolute Error')
    plt.title('Absolute Error vs n')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(data[:, 0], data[:, 2], label='δfn')
    plt.xlabel('n')
    plt.ylabel('Relative Error (%)')
    plt.title('Relative Error vs n')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(data[:, 0], data[:, 3], label='rn')
    plt.xlabel('n')
    plt.ylabel('Residual Estimate')
    plt.title('Residual Estimate vs n')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()

    root.mainloop()

if __name__ == "__main__":
    main()
