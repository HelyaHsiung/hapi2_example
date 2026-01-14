import matplotlib.pyplot as plt
from utils import compute_absorption


if __name__ == '__main__':
    # example of line-by-line data
    Molecule_Name = 'CH4'
    T = 298.6
    P = 1.0
    nu, coef = compute_absorption(Molecule_Name, 1e4/12, 1e4/2, T=T, P=P)

    plt.figure(figsize=(10, 6))
    plt.plot(1e4 / nu, coef, color='gray', linewidth=0.8, label='Voigt Profile')
    plt.title(f'{Molecule_Name} Absorption (Voigt Profile) (T={T}K, T={P}atm)')
    plt.xlabel('Wavelength ($\mu m$)')
    plt.ylabel('Absorption Coefficient ($cm^2/molecule$)')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    plt.savefig(f"Absorbance of {Molecule_Name}.png", dpi=1200)
    plt.show()
    
    # example of cross-section data
    Molecule_Name = 'SF6'
    T = 282.5
    P = 0.96541
    nu, coef = compute_absorption(Molecule_Name, 1e4/12, 1e4/2, T=T, P=P)

    plt.figure(figsize=(10, 6))
    plt.plot(1e4 / nu, coef, color='gray', linewidth=0.8, label='Cross Section')
    plt.title(f'{Molecule_Name} Absorption (T={T}K, T={P}atm)')
    plt.xlabel('Wavelength ($\mu m$)')
    plt.ylabel('Absorption Coefficient ($cm^2/molecule$)')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    plt.savefig(f"Absorbance of {Molecule_Name}.png", dpi=1200)
    plt.show()

    print("All Job Were Done!")
