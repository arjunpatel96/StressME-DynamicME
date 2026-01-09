import yaml
import sys
from pathlib import Path
import pickle

def main(config_path):
    print("→ Loading config file from...", config_path)
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Match run_dynamicme.py behavior
    project_root = Path(config["project_root"])
    sys.path.insert(0, str(project_root))

    print("→ Loading model file from...", config["model_file"])
    with open(config["model_file"], "rb") as f:
        me = pickle.load(f)

    from dynamicme.dynamic import DynamicME
    from cobrame.core.reaction import MetabolicReaction

    # Required sim params (same names as your YAML)
    T  = config["T"]     # not used for step, but we keep it for consistency/logging
    dt = config["dt"]
    V  = config["V"]
    X0 = config["X0"]
    c0_dict = config["c0_dict"]

    print("→ One-step parameters:")
    print(f"   dt    : {dt}")
    print(f"   V     : {V}")
    print(f"   X0    : {X0} g (converted to {X0 / V} g/L)")

    # --- Build extra_rxns_tracked exactly like run_dynamicme.py ---
    extra_rxns_tracked = []

    # Track metabolites -> add metabolic reactions involving _c/_e forms (same logic)
    tracked_met_ids = config.get("tracked_metabolites", [])
    if tracked_met_ids:
        tracked_mets = [me.metabolites.get_by_id(mid) for mid in tracked_met_ids]
        print("→ Tracking metabolites:", tracked_met_ids)

        for met_track in tracked_mets:
            mid_c = met_track.id.replace("_p", "_c")
            mid_e = met_track.id.replace("_p", "_e")
            met_c = me.metabolites.get_by_id(mid_c)
            met_e = me.metabolites.get_by_id(mid_e)

            for rxn in met_track.reactions:
                if (
                    isinstance(rxn, MetabolicReaction)
                    and rxn.keff
                    and (met_c in rxn.metabolites or met_e in rxn.metabolites)
                ):
                    extra_rxns_tracked.append(rxn)

    # Track translation fluxes
    if config.get("tracked_translation_reactions", False) is True:
        rxns_trsl = me.reactions.query("translation_")
        extra_rxns_tracked += list(rxns_trsl)

    # Track biomass reactions
    if config.get("tracked_biomass_to_biomass_reactions", False) is True:
        biomass_rxns = [r for r in me.reactions if r.id.endswith("_biomass_to_biomass")]
        extra_rxns_tracked += biomass_rxns

    # Track complex formation reactions
    if config.get("tracked_complex_formation_reactions", False) is True:
        cplx_rxns = [r for r in me.reactions if r.id.endswith("_formation")]
        extra_rxns_tracked += cplx_rxns

    # Bounds dictionaries (your run_dynamicme.py builds these; if you already have them in config,
    # you can just do lb_dict=config.get('lb_dict', {}) etc.)
    lb_dict = config.get("lb_dict", {})
    ub_dict = config.get("ub_dict", {})

    dyme = DynamicME(me)

    print("→ Running ONE DynamicME timestep (t=0 → t=dt)...")
    step = dyme.simulate_step(
        c0_dict=c0_dict,
        X0=X0 / V,           # g/L (matches run_dynamicme.py)
        dt=dt,
        prec_bs=1e-3,         # matches run_dynamicme.py call
        ZERO_CONC=0.0,        # matches run_dynamicme.py call
        extra_rxns_tracked=extra_rxns_tracked,
        lb_dict=lb_dict,
        ub_dict=ub_dict,
        verbosity=1,          # matches run_dynamicme.py call
    )

    print("=== Done ===")
    print(f"time:   {step['time0']} -> {step['time1']}")
    print(f"biomass:{step['biomass0']} -> {step['biomass1']}")
    print(f"tracked ex_flux keys: {len(step['ex_flux'])}")
    print(f"tracked rxn_flux keys:{len(step['rxn_flux'])}")

if __name__ == "__main__":
    main(sys.argv[1])
