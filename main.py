import matplotlib.pyplot as plt
from utils import compute_absorption_voigt


if __name__ == '__main__':
    nu, coef = compute_absorption_voigt('CH4', 1e4/12, 1e4/2)

    plt.figure(figsize=(10, 6))
    plt.plot(1e4 / nu, coef, color='gray', linewidth=0.8, label='Voigt Profile')
    plt.title(f'CH4 Absorption (Voigt Profile)')
    plt.xlabel('Wavelength ($\mu m$)')
    plt.ylabel('Absorption Coefficient Cross Section ($cm^2/molecule$)')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    plt.show()

    print("All Job Were Done!")
