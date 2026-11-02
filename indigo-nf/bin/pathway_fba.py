#!/usr/bin/env python3
"""Add indole -> indoxyl -> indigo to an E. coli GEM (BiGG ids, e.g. iML1515) and
compute theoretical yields and a growth/production envelope.
Native tnaA (TRPAS2: trp -> indole) is assumed to be present in the model."""
import argparse, gzip, io, sys
from cobra import Metabolite, Reaction
from cobra.io import read_sbml_model


def load(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return read_sbml_model(io.StringIO(f.read()))
    return read_sbml_model(path)


def add_pathway(m):
    need = ["indole_c", "o2_c", "nadph_c", "nadp_c", "h_c", "h2o_c"]
    miss = [x for x in need if x not in m.metabolites]
    if miss:
        sys.exit(f"model lacks metabolites: {miss}")
    g = m.metabolites.get_by_id
    indoxyl = Metabolite("indoxyl_c", formula="C8H7NO", charge=0, compartment="c")
    indigo = Metabolite("indigo_c", formula="C16H10N2O2", charge=0, compartment="c")

    oxy = Reaction("INDOLE_OXY", name="Indole oxygenase (FMO/NDO type), NADPH-dependent")
    oxy.add_metabolites({g("indole_c"): -1, g("o2_c"): -1, g("nadph_c"): -1, g("h_c"): -1,
                         indoxyl: 1, g("nadp_c"): 1, g("h2o_c"): 1})
    dim = Reaction("INDOXYL_DIMER", name="Indoxyl oxidative dimerisation to indigo (spontaneous)")
    dim.add_metabolites({indoxyl: -2, g("o2_c"): -1, indigo: 1, g("h2o_c"): 2})
    for r in (oxy, dim):
        r.bounds = (0, 1000)
    m.add_reactions([oxy, dim])
    for r in (oxy, dim):
        bal = r.check_mass_balance()
        if bal:
            sys.exit(f"{r.id} unbalanced: {bal}")
    m.add_boundary(indigo, type="demand")  # DM_indigo_c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--biomass", default="BIOMASS_Ec_iML1515_core_75p37M")
    ap.add_argument("--glc-ex", default="EX_glc__D_e")
    ap.add_argument("--trp-ex", default="EX_trp__L_e")
    ap.add_argument("--o2-ex", default="EX_o2_e")
    ap.add_argument("--glc", type=float, default=10.0, help="glucose uptake, mmol/gDW/h")
    ap.add_argument("--trp", type=float, default=2.0, help="max trp uptake in fed scenario")
    ap.add_argument("--o2", type=float, nargs="+", default=[5, 10, 20, 1000])
    ap.add_argument("--out", default=".")
    a = ap.parse_args()

    m = load(a.model)
    add_pathway(m)
    R = m.reactions.get_by_id
    scen = {"glucose_only": {a.glc_ex: -a.glc, a.trp_ex: 0.0},
            "glucose_plus_trp": {a.glc_ex: -a.glc, a.trp_ex: -a.trp}}

    summ = ["scenario\to2_limit\tmax_growth\tmax_indigo\tmol_indigo_per_mol_glc\tcarbon_yield"]
    env = ["scenario\tgrowth_fraction\tmax_indigo"]
    for name, bounds in scen.items():
        if name == "glucose_plus_trp" and a.trp_ex not in m.reactions:
            print(f"WARN: skipping {name}, no {a.trp_ex} in model", file=sys.stderr)
            continue
        bounds = {k: v for k, v in bounds.items() if k in m.reactions}  # glucose_only works without a trp exchange
        if a.glc_ex not in m.reactions:
            sys.exit(f"{a.glc_ex} not in model")
        for o2 in a.o2:
            with m:
                for k, v in bounds.items():
                    R(k).lower_bound = v
                R(a.o2_ex).lower_bound = -o2
                m.objective = R(a.biomass)
                mu = m.slim_optimize()
                m.objective = R("DM_indigo_c")
                sol = m.optimize()
                ind = sol.objective_value if sol.status == "optimal" else float("nan")
                glc = -sol.fluxes[a.glc_ex] if sol.status == "optimal" else float("nan")
                trp = -sol.fluxes[a.trp_ex] if a.trp_ex in m.reactions and sol.status == "optimal" else 0.0
                cin = 6 * glc + 11 * trp
                if cin and 16 * ind / cin > 1.0001:
                    print(f"WARN: carbon yield > 1 in {name} (o2={o2}): model is not carbon-balanced", file=sys.stderr)
                summ.append(f"{name}\t{o2}\t{mu:.4f}\t{ind:.4f}\t{ind / glc if glc else float('nan'):.4f}\t"
                            f"{16 * ind / cin if cin else float('nan'):.4f}")
        # growth-coupled envelope at the most permissive O2 limit
        with m:
            for k, v in bounds.items():
                R(k).lower_bound = v
            R(a.o2_ex).lower_bound = -a.o2[-1]
            m.objective = R(a.biomass)
            mu = m.slim_optimize()
            for i in range(10):
                f = i / 10
                with m:
                    R(a.biomass).lower_bound = f * mu * 0.999
                    m.objective = R("DM_indigo_c")
                    env.append(f"{name}\t{f:.1f}\t{m.slim_optimize():.4f}")

    open(f"{a.out}/fba_summary.tsv", "w").write("\n".join(summ) + "\n")
    open(f"{a.out}/envelope.tsv", "w").write("\n".join(env) + "\n")


if __name__ == "__main__":
    main()
