from hapi2.db.models import Molecule
from hapi2.web import fetch_info, fetch_molecules, fetch_isotopologues, fetch_transitions, \
    fetch_partition_functions, fetch_cross_section_headers, fetch_cross_sections, fetch_cross_section_spectra
from hapi2.db.sqlalchemy.legacy import storage2cache
from hapi2.opacity.lbl.numba.fast_abscoef import arange_
from hapi2.opacity.lbl.numba import absorptionCoefficient_Generic, absorptionCoefficient_Voigt, absorptionCoefficient_Lorentz, absorptionCoefficient_Doppler


def compute_absorption(
    molecule_name,
    numin,
    numax,
    T=296.0,
    P=1.0,
    lbl_profile='Voigt',
    HITRAN_units=True,
    diluent=None,
    step=0.01,
    wing_factor=50
):
    """
    Line-by-line absorption coefficient calculation using Voigt, Lorentz, or Doppler profile.
    Add support for cross-section data.

    Parameters
    ----------
    molecule_name : str
        HITRAN molecule name, e.g. 'CH4', 'SF6'
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

    assert mol.id <= 61, "exceed HITRAN database"

    if mol.id in [30, 35, 42, 55]:  # cross-section: Sulfur Hexafluoride, Chlorine Nitrate, PFC-14, Nitrogen Trifluoride
        headers = fetch_cross_section_headers([mol])
        best_header = None
        min_diff = float('inf')
        for h in headers:
            if h.numin <= numin and h.numax >= numax:
                temp_diff = abs(h.temperature - T)
                if temp_diff < min_diff:
                    min_diff = temp_diff
                    best_header = h
        if best_header is None:
            raise ValueError(f"Can not find {molecule_name} cross-section data in {numin}-{numax} cm-1")
        print(f"Select {molecule_name} cross-section data: T={best_header.temperature}K, range={best_header.numin}-{best_header.numax} cm-1")

        fetch_cross_section_spectra([best_header])
        wavenumber, alpha = best_header.get_data()

        mask = (wavenumber >= numin) & (wavenumber <= numax)
        wavenumber = wavenumber[mask]
        alpha = alpha[mask]

    else:  # line-by-line
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
        if lbl_profile == 'Voigt':
            wavenumber, alpha = absorptionCoefficient_Voigt(
                SourceTables=[table_name],
                Environment={'T': T, 'p': P},
                HITRAN_units=HITRAN_units,
                Diluent=diluent,
                WavenumberGrid=wngrid
            )
        elif lbl_profile == 'Lorentz':
            wavenumber, alpha = absorptionCoefficient_Lorentz(
                SourceTables=[table_name],
                Environment={'T': T, 'p': P},
                HITRAN_units=HITRAN_units,
                Diluent=diluent,
                WavenumberGrid=wngrid
            )
        elif lbl_profile == 'Doppler':
            wavenumber, alpha = absorptionCoefficient_Doppler(
                SourceTables=[table_name],
                Environment={'T': T, 'p': P},
                HITRAN_units=HITRAN_units,
                Diluent=diluent,
                WavenumberGrid=wngrid
            )
        else:
            wavenumber, alpha = absorptionCoefficient_Generic(
                profile='None',
                calcpars={'None': None},
                SourceTables=[table_name],
                Environment={'T': T, 'p': P},
                HITRAN_units=HITRAN_units,
                Diluent=diluent,
                WavenumberGrid=wngrid
            )
    return wavenumber, alpha
