from hapi2.db.models import Molecule
from hapi2.web import fetch_info, fetch_molecules, fetch_isotopologues, fetch_transitions, fetch_partition_functions
from hapi2.db.sqlalchemy.legacy import storage2cache
from hapi2.opacity.lbl.numba.fast_abscoef import arange_
from hapi2.opacity.lbl.numba import absorptionCoefficient_Voigt


def compute_absorption_voigt(
    molecule_name,
    numin,
    numax,
    T=296.0,
    P=1.0,
    HITRAN_units=True,
    diluent=None,
    step=0.01,
    wing_factor=50
):
    """
    Line-by-line absorption coefficient calculation using Voigt profile.

    Parameters
    ----------
    molecule_name : str
        HITRAN molecule name, e.g. 'CH4'
    numin, numax : float
        Requested spectral range (cm^-1)
    T : float
        Temperature (K)
    P : float
        Pressure (atm)
    diluent : dict or None
        Diluent definition, e.g. {'air':1.0}
        If None, defaults to air-broadening
    step : float
        Wavenumber grid spacing (cm^-1)
    wing_factor : float
        Line-wing extent in units of max Lorentz HWHM
    """

    if diluent is None:
        diluent = {'air': 1.0}

    table_name = f'{molecule_name}_{int(numin)}_{int(numax)}'

    # --- Molecule & isotopologue ---
    try:
        mol = Molecule(molecule_name)
    except:
        fetch_info()
        fetch_molecules()
        mol = Molecule(molecule_name)

    if len(mol.isotopologues) == 0:
        fetch_isotopologues([mol])

    main_iso = max(mol.isotopologues, key=lambda x: x.abundance)

    # Partition function is required for line intensity
    fetch_partition_functions([main_iso])

    # --- Transitions ---
    if main_iso.transitions.count() == 0:
        print(f"Download {molecule_name} spectral line: [{numin}-{numax} cm-1] ...")
        fetch_transitions([main_iso], numin, numax, table_name)

    print(f"Downloaded {table_name}: {main_iso.transitions.count()} n_points.\n")

    # --- Automatic wavenumber grid ---
    nu0 = [tr.nu for tr in main_iso.transitions]
    nu_min = min(nu0)
    nu_max = max(nu0)

    gamma_max = max(tr.gamma_air for tr in main_iso.transitions)
    margin = wing_factor * gamma_max

    wngrid = arange_(nu_min - margin, nu_max + margin, step)

    # --- Cache transitions ---
    storage2cache(table_name)

    # --- Absorption coefficient ---
    wavenumber, alpha = absorptionCoefficient_Voigt(
        SourceTables=[table_name],
        Environment={'T': T, 'p': P},
        HITRAN_units=HITRAN_units,
        Diluent=diluent,
        WavenumberGrid=wngrid
    )

    return wavenumber, alpha


