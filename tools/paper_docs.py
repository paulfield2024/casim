"""Per-file / per-routine documentation citing Field et al. (2023, QJRMS,
https://doi.org/10.1002/qj.4414) "Implementation of a double moment cloud
microphysics scheme in the UK Met Office regional numerical weather
prediction model", Appendix A ("CASIM implementation in the Unified
Model"), which explicitly cross-references CASIM source modules.

Each entry is additive documentation only (never touches executable
code). "method" and "reference" are lists of plain text lines (no "!"
prefix - the inserter adds that and wraps them into the standard header
block). "routines" maps a SUBROUTINE/FUNCTION name (as it appears in the
cleaned, upper-cased source) to a short list of comment lines inserted
immediately above its declaration.

Where Field et al. (2023) does not describe a component with a named
equation, that is stated explicitly rather than inventing a citation.
"""

DOCS = {
    "autoconversion.F90": {
        "method": [
            "Autoconversion of cloud liquid to rain: the self-collection of",
            "cloud droplets to form drizzle/rain-sized drops.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.4.1, Eq. (A4)-(A6): mass rate",
            "Praut follows Khairoutdinov & Kogan (2000) (optionally Kogan,",
            "2013); number rate Nwaut = Praut*nw/qw.",
        ],
        "routines": {
            "RAUT": [
                "Praut/Nwaut: Khairoutdinov & Kogan (2000) autoconversion,",
                "Field et al. (2023) Eq. (A4)-(A6).",
            ],
        },
    },
    "accretion.F90": {
        "method": [
            "Accretion (collection) of cloud liquid water by rain drops.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.5.1, Eq. (A10)-(A11):",
            "Pracw = 67(qw*qr)^1.15 following Khairoutdinov & Kogan (2000);",
            "Nracw = Pracw*qw/nw.",
        ],
        "routines": {
            "RACW": [
                "Pracw/Nracw: Khairoutdinov & Kogan (2000) accretion,",
                "Field et al. (2023) Eq. (A10)-(A11).",
            ],
        },
    },
    "aggregation.F90": {
        "method": [
            "Self-collection (aggregation) processes: rain-rain and",
            "snow-snow self collection, which change number concentration",
            "only (no mass change), per hydrometeor type.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.5.3, Eq. (A17)-(A18): rain",
            "self-collection follows Beheng (1994); snow self-collection",
            "follows the double-moment Large Eddy Model formulation of",
            "Gray et al. (2001).",
        ],
        "routines": {
            "RACR": [
                "Rain-rain self-collection number tendency Nracr,",
                "Field et al. (2023) Eq. (A17), after Beheng (1994).",
            ],
            "ICE_AGGREGATION": [
                "Snow-snow self-collection number tendency Nsacs,",
                "Field et al. (2023) Eq. (A18), after Gray et al. (2001).",
            ],
        },
    },
    "breakup.F90": {
        "method": [
            "Mechanical breakup of large snow aggregates once their",
            "mass-weighted mean diameter exceeds a fixed threshold (DSbrk),",
            "increasing number concentration with no change in mass.",
        ],
        "reference": [
            "Field et al. (2023) does not give a named equation for this",
            "breakup process; it is conceptually related to the size",
            "distribution stability limiting discussed qualitatively in",
            "Appendix A.13, but is implemented here as a distinct physical",
            "process rather than the diagnostic re-scaling described there.",
        ],
        "routines": {
            "ICE_BREAKUP": [
                "Snow mechanical breakup above diameter threshold DSbrk;",
                "not given as a numbered equation in Field et al. (2023).",
            ],
        },
    },
    "snow_autoconversion.F90": {
        "method": [
            "Autoconversion of cloud ice crystals to snow aggregates by",
            "diffusional growth past a minimum ice size distribution slope.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.4.2, Eq. (A7)-(A9): Psaut for",
            "lambda_i < lambda_i,min (Di,max = 100 micron), Nsaut assumes",
            "new snow particles have diameter Di,s = 50 micron, and Niaut",
            "is the corresponding cloud-ice number sink.",
        ],
        "routines": {
            "SAUT": [
                "Psaut/Nsaut/Niaut ice-to-snow autoconversion,",
                "Field et al. (2023) Eq. (A7)-(A9).",
            ],
        },
    },
    "evaporation.F90": {
        "method": [
            "Evaporation of rain falling into sub-saturated air.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.8.2, Eq. (A25)-(A27): Prevp",
            "uses the electrostatic-analogue growth-rate framework of",
            "A.8.1 (Eq. A21-A24) with the liquid thermodynamic term ABliq",
            "(Eq. A26); Nrevp assumes homogeneous or inhomogeneous",
            "evaporation depending on l_sigevap.",
        ],
        "routines": {
            "REVP": [
                "Prevp/Nrevp rain evaporation, Field et al. (2023)",
                "Eq. (A25)-(A27).",
            ],
        },
    },
    "ice_deposition.F90": {
        "method": [
            "Vapour deposition growth (or sublimation) of cloud ice, snow",
            "and graupel by diffusion of water vapour.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.8.3, Eq. (A28)-(A29): Pxdep",
            "(or -Pxsub) uses the ice thermodynamic term ABice; number",
            "concentration of snow/graupel decreases at rate -Pxsub*nx/qx",
            "during sublimation, with no number source during deposition.",
        ],
        "routines": {
            "IDEP": [
                "Deposition/sublimation of ice/snow/graupel,",
                "Field et al. (2023) Eq. (A28)-(A29).",
            ],
        },
    },
    "ice_melting.F90": {
        "method": [
            "Melting of cloud ice (instantaneous) and thermal-balance",
            "melting of snow and graupel into rain.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.8.4, Eq. (A30)-(A34): Pimlt",
            "= qi/Delta t; Psmlt/Pgmlt from thermal heat balance (including",
            "the riming heating terms Psacw/Psacr, Pgacw/Pgacr); Nsmlt/Ngmlt",
            "assume mean size is preserved during melting.",
        ],
        "routines": {
            "MELTING": [
                "Ice/snow/graupel melting, Field et al. (2023)",
                "Eq. (A30)-(A34).",
            ],
        },
    },
    "ice_nucleation.F90": {
        "method": [
            "Primary heterogeneous ice nucleation: production of cloud ice",
            "number/mass from water vapour, either by nudging to a",
            "temperature-dependent concentration or via aerosol/dust-based",
            "ice-nucleating-particle parametrisations.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.9.1: default nudges cloud ice",
            "number to Cooper (1987) at water saturation and T < -8C;",
            "aerosol-aware options use DeMott et al. (2010, 2015), Niemand",
            "et al. (2012), Atkinson et al. (2013) or Tobo et al. (2013),",
            "as summarised by Miltenberger et al. (2020). No single",
            "numbered equation is given for these parametrisations.",
        ],
        "routines": {
            "INUC": [
                "Heterogeneous ice nucleation (Cooper 1987 default, or",
                "aerosol-based schemes), Field et al. (2023) Sec. A.9.1.",
            ],
        },
    },
    "homogeneous_freezing.F90": {
        "method": [
            "Heterogeneous immersion freezing of rain drops to graupel",
            "(Bigg, 1953) and homogeneous freezing of cloud droplets to",
            "ice below the homogeneous-freezing temperature threshold.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.9.2, Eq. (A35)-(A37): rain",
            "freezing probability follows Bigg (1953); Appendix A.9.3,",
            "Eq. (A38): homogeneous droplet freezing Nhomc balances the",
            "Squires supersaturation equation for a given updraught,",
            "assuming 50 micron ice spheres are formed.",
        ],
        "routines": {
            "IHOM_RAIN": [
                "Bigg (1953) immersion freezing of rain to graupel,",
                "Field et al. (2023) Eq. (A35)-(A37).",
            ],
            "IHOM_DROPLETS": [
                "Homogeneous freezing of cloud droplets below the -38C",
                "threshold, Field et al. (2023) Eq. (A38).",
            ],
        },
    },
    "ice_multiplication.F90": {
        "method": [
            "Secondary ice production: rime-splintering (Hallett-Mossop)",
            "and (optionally) droplet-shattering / ice-ice collisional",
            "break-up mechanisms.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.10, Eq. (A39): Hallett-Mossop",
            "splinter production is triangular-weighted between -2.5C and",
            "-7.5C (peaking at -5C), assuming 350 splinters per 1e-6 kg of",
            "rimed liquid, each of mass 1e-18 kg (Hallett & Mossop, 1974).",
            "Droplet-shattering and ice-ice collisional break-up are not",
            "described by a numbered equation in Field et al. (2023).",
        ],
        "routines": {
            "HALLET_MOSSOP": [
                "Rime-splintering secondary ice production,",
                "Field et al. (2023) Eq. (A39), after Hallett & Mossop (1974).",
            ],
            "DROPLET_SHATTERING": [
                "Optional secondary ice production by droplet shattering;",
                "not given as a numbered equation in Field et al. (2023).",
            ],
            "ICE_COLLISION": [
                "Optional secondary ice production by ice-ice collisional",
                "break-up; not given as a numbered equation in Field et al.",
                "(2023).",
            ],
        },
    },
    "ice_accretion.F90": {
        "method": [
            "Mixed-phase and cold-phase accretion (collection) between",
            "pairs of different hydrometeor species (e.g. ice-water,",
            "snow-water, graupel-water, snow-ice, graupel-ice/snow).",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.5.2, Eq. (A12)-(A16): uses",
            "either the simple gravitational sweepout (Eq. A12, small",
            "fall-speed collected species) or the full binary collection",
            "equation (Eq. A13-A14, large fall-speed species), both",
            "provided by src/sweepout_rate.F90; collection efficiencies",
            "Exy are listed in their Table A3.",
        ],
        "routines": {
            "IACC": [
                "Mixed-/cold-phase accretion via sweepout or binary",
                "collection, Field et al. (2023) Eq. (A12)-(A16), Table A3.",
            ],
        },
    },
    "graupel_wetgrowth.F90": {
        "method": [
            "Wet/dry growth of graupel: determines whether accreted",
            "liquid freezes fully (dry growth) or is partly shed as rain",
            "because latent heating brings the graupel surface to 0C",
            "(wet growth), following Musil (1970).",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.11, Eq. (A40)-(A43): Pgwet",
            "threshold from Musil (1970), with the (910/rho_g)^0.625",
            "factor adopted from the Met Office Large Eddy Model; shedding",
            "Pgshd computed when Pgacr is negative (Eq. A43).",
        ],
        "routines": {
            "WETGROWTH": [
                "Graupel wet/dry growth and shedding, Field et al. (2023)",
                "Eq. (A40)-(A43), after Musil (1970).",
            ],
        },
    },
    "graupel_embryo.F90": {
        "method": [
            "Formation of new graupel embryos from snow riming past a",
            "critical rate (currently disabled in the UM configuration).",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.11.1, Eq. (A44)-(A45), after",
            "Reisner et al. (1998); the paper notes this process is",
            "currently disabled in this UM implementation.",
        ],
        "routines": {
            "GRAUPEL_EMBRYOS": [
                "Snow-to-graupel embryo formation, Field et al. (2023)",
                "Eq. (A44)-(A45), after Reisner et al. (1998). Disabled by",
                "default in the UM configuration described in the paper.",
            ],
        },
    },
    "sedimentation.F90": {
        "method": [
            "Sedimentation (gravitational settling) of all hydrometeor",
            "mass/number/third moments, including CFL-based substepping",
            "so hydrometeors do not skip whole grid levels in one step.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.12, Eq. (A46)-(A48): terminal",
            "velocity Vx(D) (Eq. A46, Table A1 parameters); analytic flux-",
            "divergence update (Eq. A47-A48) following the exponential-",
            "filter approach of Rotstayn (1997), as already used for",
            "sedimentation in the UM (Wilson & Ballard, 1999).",
        ],
        "routines": {
            "SEDR": [
                "Main sedimentation update, Field et al. (2023)",
                "Eq. (A46)-(A48), after Rotstayn (1997).",
            ],
            "SEDR_1M_2M": [
                "Single-/double-moment variant of the sedimentation",
                "update, Field et al. (2023) Eq. (A46)-(A48).",
            ],
            "TERMINAL_VELOCITY_CFL": [
                "Computes Courant number alpha = Vx*dt/dz (Field et al.,",
                "2023, Eq. A48) to determine sedimentation substepping.",
            ],
        },
    },
    "ventfac.F90": {
        "method": [
            "Ventilation factor for diffusional growth/evaporation of",
            "falling hydrometeors, integrated over the particle size",
            "distribution.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.8.1, Eq. (A23)-(A24): F is",
            "the per-particle ventilation coefficient; the integrated",
            "ventilation factor chi_x combines capacitance, F and the PSD.",
        ],
        "routines": {
            "VENTILATION_3M": [
                "Integrated ventilation factor for a triple-moment PSD,",
                "Field et al. (2023) Eq. (A23)-(A24).",
            ],
            "VENTILATION_1M_2M": [
                "Integrated ventilation factor for single-/double-moment",
                "PSDs, Field et al. (2023) Eq. (A23)-(A24).",
            ],
        },
    },
    "sweepout_rate.F90": {
        "method": [
            "Gravitational sweepout and binary collection kernels used by",
            "accretion/aggregation routines to compute collection rates",
            "between two particle size distributions.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.5.2, Eq. (A12) (simple",
            "sweepout for small fall-speed collected species) and",
            "Eq. (A13)-(A14) (full binary collection equation for large",
            "fall-speed species), referenced there as",
            "src/sweepout_rate.F90 and src/binary_collection.F90 (the",
            "latter's functionality is provided by this module here).",
        ],
        "routines": {
            "SWEEPOUT": [
                "Simple gravitational sweepout kernel, Field et al. (2023)",
                "Eq. (A12).",
            ],
            "SWEEPOUT_1M2M": [
                "Simple gravitational sweepout kernel, single-/double-",
                "moment variant, Field et al. (2023) Eq. (A12).",
            ],
            "BINARY_COLLECTION": [
                "Full binary collection equation, Field et al. (2023)",
                "Eq. (A13)-(A14).",
            ],
            "BINARY_COLLECTION_1M2M": [
                "Full binary collection equation, single-/double-moment",
                "variant, Field et al. (2023) Eq. (A13)-(A14).",
            ],
        },
    },
    "condensation.F90": {
        "method": [
            "Saturation-adjustment condensation/evaporation of cloud",
            "water (used when CASIM's own condensation is active, e.g.",
            "in MONC/KiD, rather than the UM cloud-fraction scheme).",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.6.1 describes the",
            "\"all-or-nothing\" saturation-adjustment scheme used when",
            "CASIM's own condensation step is active; no numbered",
            "equation is given for it there (Eq. A19-A20 instead cover",
            "the UM cloud-fraction coupling, not this module).",
        ],
        "routines": {
            "CONDEVP": [
                "All-or-nothing saturation-adjustment condensation/",
                "evaporation, Field et al. (2023) Sec. A.6.1.",
            ],
        },
    },
    "activation.F90": {
        "method": [
            "Droplet activation from aerosol: generates increments in",
            "cloud droplet number (and mass/aerosol number where coupled",
            "to an aerosol scheme) from a given supersaturation/updraught.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.6.2: droplet number increment",
            "may only increase unless cloud fraction decreases, following",
            "Stevens et al. (1996); the aerosol-based calculation itself",
            "is provided by src/shipway_activation/shipway_activation_mod.F90",
            "(Gordon et al., 2020, adapting Abdul-Razzak & Ghan, 2000).",
        ],
        "routines": {
            "ACTIVATE": [
                "Droplet number/mass activation increment, Field et al.",
                "(2023) Sec. A.6.2, after Stevens et al. (1996).",
            ],
        },
    },
    "distributions.F90": {
        "method": [
            "Maintains the diagnosed particle-size-distribution (PSD)",
            "shape parameters (lambda, mu, N0) per hydrometeor/grid point,",
            "shared by all process-rate routines via module state.",
        ],
        "reference": [
            "Field et al. (2023) Sec. 1/Appendix A: CASIM represents each",
            "hydrometeor's PSD by a generalised gamma function with fixed",
            "shape parameter (one, two, or three prognostic moments); the",
            "size-distribution stability limiting described qualitatively",
            "in Appendix A.13 is applied via this module's lambda/mu state.",
        ],
        "routines": {},
    },
    "lookup.F90": {
        "method": [
            "Gamma-function moment relations (Gfunc/Hfunc) and PSD",
            "parameter inversion (recovering N0, lambda, mu from",
            "prognosed moments) used throughout the process-rate code.",
        ],
        "reference": [
            "Field et al. (2023) does not give numbered equations for",
            "these generic gamma-distribution moment/inversion utilities;",
            "they implement the general PSD framework summarised in",
            "Sec. 1 and used implicitly throughout Appendix A (e.g. the",
            "Gamma(.)/lambda(.) terms in Eq. A12-A18, A36-A38, A44).",
        ],
        "routines": {
            "GFUNC": [
                "Generalised-gamma moment function relating shape",
                "parameter mu to three chosen moment powers p1,p2,p3.",
            ],
            "HFUNC": [
                "Inverse of Gfunc: recovers the third moment given two",
                "known moments and the shape parameter relation.",
            ],
        },
    },
    "m3_incs.F90": {
        "method": [
            "Increments to a hydrometeor's third prognostic moment (e.g.",
            "reflectivity-related moment for triple-moment species) when",
            "mass and/or number are incremented by a process rate.",
        ],
        "reference": [
            "Field et al. (2023) Sec. 1 notes rain/snow/graupel can carry",
            "a third prognostic moment; no numbered equation is given for",
            "these specific moment-increment helper routines.",
        ],
        "routines": {},
    },
    "casim_moments_mod.F90": {
        "method": [
            "Diagnoses higher moments of the ice/snow PSD (e.g. for",
            "coupling to radiation or other schemes needing a specific",
            "moment) from the prognostic mass/number/third moment.",
        ],
        "reference": [
            "Field et al. (2023) Sec. 1 describes the general multi-",
            "moment PSD representation; no numbered equation is given",
            "for this specific moment-diagnosis coupling routine.",
        ],
        "routines": {},
    },
    "lognormal_funcs.F90": {
        "method": [
            "Converts between lognormal aerosol-mode moments (mass M,",
            "number N) and the modal (number-mean) radius Rm.",
        ],
        "reference": [
            "Supports the MURK/ARCL aerosol-to-droplet-number coupling",
            "described qualitatively in Field et al. (2023) Sec. A.7-A.8;",
            "no numbered equation is given there for this conversion.",
        ],
        "routines": {
            "MNTORM": [
                "Lognormal mass/number -> modal radius conversion, used",
                "by the MURK/ARCL aerosol coupling (Field et al., 2023,",
                "Sec. A.7-A.8).",
            ],
        },
    },
    "gauss_4A_func.F90": {
        "method": [
            "Lookup-table evaluation of a Gaussian-integral function used",
            "by the aerosol activation parametrisation.",
        ],
        "reference": [
            "Supports the Gordon et al. (2020) activation scheme",
            "referenced in Field et al. (2023) Sec. A.6.2; no numbered",
            "equation for this specific lookup helper is given there.",
        ],
        "routines": {},
    },
    "special.F90": {
        "method": [
            "Special mathematical functions (Gamma function, error",
            "function/complementary error function and its inverse) used",
            "throughout the PSD and ventilation/deposition calculations.",
        ],
        "reference": [
            "Field et al. (2023) does not describe these generic special-",
            "function implementations directly; they underpin the",
            "Gamma(.) terms used throughout Appendix A (e.g. Eq. A5,",
            "A8-A9, A12-A18, A36-A38, A44).",
        ],
        "routines": {},
    },
    "qsat_casim_func.F90": {
        "method": [
            "Saturation mixing ratio with respect to liquid water",
            "(Qsaturation) and ice (Qisaturation), used by evaporation,",
            "deposition/sublimation and melting rate calculations.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.8.1-A.8.4 use qsat/qisat",
            "(Eq. A21-A22, A25-A26, A28-A29) computed by this module; the",
            "saturation-vapour-pressure formula itself is not spelled out",
            "in the paper.",
        ],
        "routines": {
            "QSATURATION": [
                "Saturation mixing ratio w.r.t. liquid water, used as",
                "qwsat in Field et al. (2023) Eq. (A25)-(A26).",
            ],
            "QISATURATION": [
                "Saturation mixing ratio w.r.t. ice, used as qisat in",
                "Field et al. (2023) Eq. (A28)-(A29), (A38).",
            ],
        },
    },
    "adjust_deposition.F90": {
        "method": [
            "Limits/apportions total vapour deposition onto ice, snow and",
            "graupel so that the combined sink does not exceed the",
            "available supersaturated vapour in a timestep.",
        ],
        "reference": [
            "Supports the deposition/sublimation processes of Field et",
            "al. (2023) Appendix A.8.3 (Eq. A28-A29); the specific",
            "multi-species vapour-budget limiter is a numerical safeguard",
            "not itself described by a numbered equation in the paper.",
        ],
        "routines": {},
    },
    "which_mode_to_use.F90": {
        "method": [
            "Decides, for aerosol processing, whether a droplet/particle",
            "should be treated via the simple mass-threshold method or a",
            "more complete aerosol-mode-tracking method.",
        ],
        "reference": [
            "Supports the aerosol-cloud coupling discussed qualitatively",
            "in Field et al. (2023) Sec. 1 and Appendix A.6.2-A.8; this",
            "internal dispatch logic is not described by a numbered",
            "equation in the paper.",
        ],
        "routines": {},
    },
    "aerosol_routines.F90": {
        "method": [
            "Shared helper routines for transferring aerosol mass/number",
            "between the free-aerosol and in-cloud/in-hydrometeor aerosol",
            "prognostics during activation, nucleation and collection.",
        ],
        "reference": [
            "Field et al. (2023) Sec. 1 describes activated aerosol being",
            "carried as an in-cloud prognostic transported by dynamics",
            "and sedimentation (citing Ghan & Easter, 2006); the specific",
            "bookkeeping routines here are not given numbered equations.",
        ],
        "routines": {},
    },
    "casim_reflec_mod.F90": {
        "method": [
            "Diagnoses radar reflectivity (dBZ) per hydrometeor species",
            "and total, for model verification against observed radar",
            "composites.",
        ],
        "reference": [
            "Field et al. (2023) Sec. 3.3/4.3 compares simulated and",
            "observed radar reflectivity CFADs, but the paper does not",
            "give the reflectivity calculation formula itself, so no",
            "equation citation is given here (do not invent one).",
        ],
        "routines": {},
    },
    "process_routines.F90": {
        "method": [
            "Defines the process_rate derived type and the master list",
            "of process-rate indices (i_praut, i_pracw, ...) shared by",
            "every physics-process module to accumulate tendencies.",
        ],
        "reference": [
            "Provides the software infrastructure for the process-rate",
            "subscript naming convention described in Field et al. (2023)",
            "Table A2; it is CASIM infrastructure, not itself the subject",
            "of a numbered equation.",
        ],
        "routines": {},
    },
    "micro_main.F90": {
        "method": [
            "Top-level microphysics driver (shipway_microphysics /",
            "microphysics_common): loops over grid columns/levels,",
            "preconditions state, dispatches every enabled physics",
            "process in sequence, and gathers diagnostics.",
        ],
        "reference": [
            "Corresponds to the overall CASIM time-step structure shown",
            "schematically in Field et al. (2023) Figure A1 and described",
            "throughout Appendix A; the specific dispatch/loop logic here",
            "is CASIM software infrastructure, not itself a numbered",
            "equation.",
        ],
        "routines": {},
    },
    "mphys_parameters.F90": {
        "method": [
            "Defines the hydro_params derived type (per-hydrometeor",
            "terminal-fall-speed, mass-dimension and shape parameters)",
            "and the module-level cloud/rain/ice/snow/graupel parameter",
            "sets used throughout the process-rate code.",
        ],
        "reference": [
            "Field et al. (2023) Table A1 tabulates exactly these",
            "terminal-fall-speed (a,b,f), mass-dimension (c,d) and shape",
            "(mu) parameters for each of the five hydrometeor species.",
        ],
        "routines": {},
    },
    "mphys_constants.F90": {
        "method": [
            "Physical constants (e.g. specific heats, latent heats, gas",
            "constants, reference density) shared across the process-rate",
            "modules.",
        ],
        "reference": [
            "Provides the thermodynamic constants (cp, Lv, Ls, Lf, Rv,",
            "Ka, rho0, ...) used throughout Field et al. (2023) Appendix",
            "A (e.g. Eq. A22, A26, A29, A31-A32, A36-A38); this module",
            "itself is CASIM infrastructure, not a numbered equation.",
        ],
        "routines": {},
    },
    "mphys_switches.F90": {
        "method": [
            "Module-level logical/integer switches controlling which",
            "hydrometeor species, moments and physics processes are",
            "active for a given CASIM configuration.",
        ],
        "reference": [
            "Implements the configurability described in Field et al.",
            "(2023) Sec. 1 (\"CASIM has been developed in such a way as to",
            "make it configurable\"); the switches themselves are software",
            "infrastructure, not a numbered equation.",
        ],
        "routines": {},
    },
    "initialize.F90": {
        "method": [
            "One-off initialisation entry points: sets microphysics",
            "switches, initialises look-up tables (gamma function, mu",
            "look-up, sedimentation, aerosol) before the first timestep.",
        ],
        "reference": [
            "CASIM software infrastructure supporting the physics",
            "described throughout Field et al. (2023) Appendix A; not",
            "itself the subject of a numbered equation.",
        ],
        "routines": {},
    },
    "mphys_die.F90": {
        "method": [
            "Error-reporting helper (mphys_message) used to print a",
            "diagnostic message and stop when CASIM detects an",
            "unrecoverable configuration or runtime error.",
        ],
        "reference": [
            "CASIM software infrastructure; not discussed in Field et",
            "al. (2023).",
        ],
        "routines": {},
    },
    "mphys_tidy.F90": {
        "method": [
            "Post-process clean-up: removes negligibly small residual",
            "mass/number left in a hydrometeor after the process rates",
            "have been applied, to avoid numerical noise.",
        ],
        "reference": [
            "Related to the threshold/limiting philosophy discussed",
            "qualitatively in Field et al. (2023) Appendix A.13, but the",
            "specific clean-up thresholds here are not given a numbered",
            "equation in the paper.",
        ],
        "routines": {},
    },
    "passive_fields.F90": {
        "method": [
            "Holds passively-diagnosed environmental fields (density,",
            "vertical velocity, exner pressure, temperature, relative",
            "humidity threshold) passed in from the host model each",
            "timestep/column.",
        ],
        "reference": [
            "CASIM software infrastructure providing the environmental",
            "inputs (e.g. rho, w, T) used throughout Field et al. (2023)",
            "Appendix A equations; not itself a numbered equation.",
        ],
        "routines": {},
    },
    "preconditioning.F90": {
        "method": [
            "Decides which grid columns/levels actually need the full",
            "microphysics calculation (i.e. contain hydrometeors or are",
            "sub-saturated), to skip unnecessary work elsewhere.",
        ],
        "reference": [
            "CASIM software/performance infrastructure; not discussed in",
            "Field et al. (2023).",
        ],
        "routines": {},
    },
    "sum_procs.F90": {
        "method": [
            "Sums the individual process-rate contributions (module",
            "sum_process) into net mass/number tendencies per",
            "hydrometeor, applied to update the prognostic fields each",
            "(sub)step.",
        ],
        "reference": [
            "Implements the right-hand-side summation of the tendency",
            "equations Dqv/Dt ... Dng/Dt given in Field et al. (2023)",
            "Eq. (A3); the summation code itself is CASIM infrastructure.",
        ],
        "routines": {},
    },
    "generic_diagnostic_variables.F90": {
        "method": [
            "Defines the casdiags derived type and allocation logic for",
            "all optional diagnostic output fields (process rates, water",
            "paths, surface precipitation rates, radar reflectivity).",
        ],
        "reference": [
            "Diagnostic bookkeeping infrastructure supporting output of",
            "the process rates named in Field et al. (2023) Table A2 and",
            "the water-path/precipitation diagnostics discussed in",
            "Sec. 4; not itself a numbered equation.",
        ],
        "routines": {},
    },
    "derived_constants.F90": {
        "method": [
            "Computes constants derived from the base physical constants",
            "and hydrometeor parameters at initialisation (e.g. combined",
            "shape/gamma-function prefactors reused every timestep).",
        ],
        "reference": [
            "Precomputes constant prefactors used in the Field et al.",
            "(2023) Appendix A rate equations (e.g. the Gamma(.) terms in",
            "Eq. A5, A8-A9, A12-A18); this module itself is CASIM",
            "infrastructure.",
        ],
        "routines": {},
    },
    "type_aerosol.F90": {
        "method": [
            "Defines derived types describing an aerosol species/mode",
            "(mass, number, solubility/hygroscopicity, size distribution",
            "parameters) used by the activation and nucleation code.",
        ],
        "reference": [
            "Supports the aerosol representations (MURK, ARCL, UKCA-mode)",
            "discussed in Field et al. (2023) Sec. A.6.2-A.8; the type",
            "definitions themselves are CASIM infrastructure.",
        ],
        "routines": {},
    },
    "type_process.F90": {
        "method": [
            "Defines the process_name derived type used to identify",
            "which process-rate index a given calculation should",
            "accumulate into.",
        ],
        "reference": [
            "Software infrastructure supporting the process-rate naming",
            "convention in Field et al. (2023) Table A2; not itself a",
            "numbered equation.",
        ],
        "routines": {},
    },
    "variable_precision.F90": {
        "method": [
            "Defines the working-precision kind parameter (wp) used for",
            "all real variables throughout CASIM.",
        ],
        "reference": [
            "CASIM software infrastructure; not discussed in Field et al.",
            "(2023).",
        ],
        "routines": {},
    },
    "precision.F90": {
        "method": [
            "Defines additional numeric kind parameters used where a",
            "specific (not necessarily working) precision is required.",
        ],
        "reference": [
            "CASIM software infrastructure; not discussed in Field et al.",
            "(2023).",
        ],
        "routines": {},
    },
    "thresholds.F90": {
        "method": [
            "Small-value thresholds below which a hydrometeor's mass/",
            "number is treated as negligible/zero, and the size-",
            "distribution min/max slope limits used to keep the PSD",
            "shape stable.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.13 (\"Size distribution",
            "limiting\") describes qualitatively this same stability",
            "safeguard (re-scaling lambda when the mean size exceeds a",
            "large/small threshold); specific numeric threshold values",
            "are implementation choices not tabulated in the paper.",
        ],
        "routines": {},
    },
    "casim_parent.F90": {
        "method": [
            "Records which host model (UM, MONC, KiD, standalone) is",
            "driving CASIM, so that host-specific code paths (e.g. UM",
            "cloud-fraction coupling vs CASIM's own saturation-adjustment",
            "condensation) can be selected.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.6.1 explains that CASIM's own",
            "condensation/activation step is disabled when running inside",
            "the UM (condensation and cloud fraction are instead supplied",
            "by the UM cloud scheme) but active in MONC/KiD; this module",
            "implements that host-selection switch.",
        ],
        "routines": {},
    },
    "casim_runtime.F90": {
        "method": [
            "Runtime bookkeeping (e.g. timestep/substep counters, run-",
            "time flags) shared across CASIM modules.",
        ],
        "reference": [
            "CASIM software infrastructure; not discussed in Field et al.",
            "(2023).",
        ],
        "routines": {},
    },
    "casim_stph.F90": {
        "method": [
            "Optional stochastic physics perturbation hooks for CASIM",
            "process rates (used for ensemble/stochastic-physics",
            "configurations).",
        ],
        "reference": [
            "Stochastic physics is not discussed in Field et al. (2023),",
            "which describes only the deterministic CASIM-2M",
            "configuration; no citation is invented here.",
        ],
        "routines": {},
    },
    "cloud_fraction_dummy.F90": {
        "method": [
            "Stand-in cloud-fraction scheme used when CASIM is run",
            "without a host cloud-fraction scheme (e.g. standalone/",
            "column-model use), providing a trivial 0/1 cloud fraction.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.6.1 describes the",
            "\"all-or-nothing\" (cloud fraction of 1 or 0) saturation-",
            "adjustment assumption used by CASIM's own condensation when",
            "no host cloud-fraction scheme (e.g. the UM bimodal scheme,",
            "Van Weverberg et al., 2021) is coupled in.",
        ],
        "routines": {},
    },
    "dust_hack.F90": {
        "method": [
            "Simplified/placeholder dust-aerosol field handling used",
            "where a full dust-transport scheme is not coupled in.",
        ],
        "reference": [
            "Supports the dust-based ice-nucleating-particle",
            "parametrisations mentioned in Field et al. (2023)",
            "Sec. A.9.1, but this simplified stand-in is not itself",
            "described in the paper.",
        ],
        "routines": {},
    },
    "ship_tracks.F90": {
        "method": [
            "Optional representation of aerosol injection along ship",
            "tracks, for studies of shipping-aerosol effects on cloud.",
        ],
        "reference": [
            "Not discussed in Field et al. (2023), which focuses on the",
            "operational regional NWP configuration; no citation is",
            "invented here.",
        ],
        "routines": {},
    },
    "shipway_activation/shipway_activation_mod.F90": {
        "method": [
            "Aerosol activation scheme: derives the maximum",
            "supersaturation reached in an updraught and the resulting",
            "number of aerosol particles activated to cloud droplets,",
            "accounting for pre-existing droplets competing for water",
            "vapour with un-activated aerosol.",
        ],
        "reference": [
            "Field et al. (2023) Appendix A.6.2: implements the Gordon",
            "et al. (2020) methodology (an adaptation of Abdul-Razzak &",
            "Ghan, 2000) that separately computes the maximum",
            "supersaturation assuming aerosol dominates the vapour sink",
            "and assuming droplets dominate it, taking the lower value;",
            "the updraught used is wact = w + c*sqrt(tke) with c=0 for",
            "the ARCL/MURK configurations described in the paper.",
        ],
        "routines": {
            "CALC_NCCN": [
                "Number of aerosol activated (Nccn) at a given",
                "supersaturation, Field et al. (2023) Sec. A.6.2, after",
                "Gordon et al. (2020) / Abdul-Razzak & Ghan (2000).",
            ],
            "SOLVE_NCCN_HOUSEHOLD": [
                "Householder-iteration solver for the maximum",
                "supersaturation, Field et al. (2023) Sec. A.6.2.",
            ],
            "SOLVE_NCCN_BRENT": [
                "Brent's-method solver for the maximum supersaturation,",
                "Field et al. (2023) Sec. A.6.2.",
            ],
        },
    },
    "shipway_activation/shipway_constants.F90": {
        "method": [
            "Physical constants (gas constant, molecular weight of",
            "water, surface tension, densities) specific to the Shipway",
            "activation scheme.",
        ],
        "reference": [
            "Provides constants used by the Gordon et al. (2020)",
            "activation calculation cited in Field et al. (2023)",
            "Sec. A.6.2; the constants themselves are not tabulated in",
            "the paper.",
        ],
        "routines": {},
    },
    "shipway_activation/shipway_erf.F90": {
        "method": [
            "Error function / complementary error function evaluation",
            "used by the Shipway activation scheme's Kohler-theory",
            "calculations.",
        ],
        "reference": [
            "Generic special-function support for the activation scheme",
            "of Field et al. (2023) Sec. A.6.2; not itself described by",
            "a numbered equation in the paper.",
        ],
        "routines": {},
    },
    "shipway_activation/shipway_lookup.F90": {
        "method": [
            "Lookup-table grid definitions/tolerances used to speed up",
            "the Shipway activation scheme's iterative solve.",
        ],
        "reference": [
            "Numerical-performance infrastructure for the activation",
            "scheme of Field et al. (2023) Sec. A.6.2; not itself",
            "described by a numbered equation in the paper.",
        ],
        "routines": {},
    },
    "shipway_activation/shipway_parameters.F90": {
        "method": [
            "Aerosol-mode parameters (number concentration, dry radius,",
            "geometric standard deviation, solubility/hygroscopicity per",
            "mode) used as input to the Shipway activation scheme.",
        ],
        "reference": [
            "Field et al. (2023) Sec. A.7 describes, for the MURK",
            "coupling, a lognormal ammonium-sulphate mode with mode size",
            "9.5e-8 m, geometric standard deviation 1.4, density",
            "1769 kg/m3 and hygroscopicity B=0.4 supplied via parameters",
            "of this kind.",
        ],
        "routines": {},
    },
}
