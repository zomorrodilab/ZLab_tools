import cobra

from cobra_utils import load_model


def optimize_model(model_input: str or cobra.Model, add_1ba: bool = False) -> dict:
    """
    Optimizes a multi-species model of metabolism by maximizing the flux through
    all fecal transporter (UFEt) reactions and minimizing the flux through all
    fecal exchange reactions (IEX). The maximized UFEt fluxes are then used to
    set the bounds of the IEX reactions.

    Currently, this function works with multi-species models generated by
    mgPipe.m from the Microbiome Modeling Toolbox by Heinken et al. (2022).

    Parameters
    ----------
    model_input : str or cobra.Model
        Path to a multi-species model in any COBRApy supported format or a
        COBRApy model loaded into memory.
    add_1ba : bool, optional
        If True, will set diet 1ba bounds to the model before optimizing, by
        default False.

    Returns
    -------
    dict
        A dictionary containing the maximized UFEt reaction fluxes.
    dict
        A dictionary containing the minimized IEX reaction fluxes.
    dict
        A dictionary containing the model reaction bounds after minimization of
        the IEX fluxes.

    Raises
    ------
    ValueError
        If the model_input is not a path to a model or a COBRApy model object.
    """
    # Load the model
    if isinstance(model_input, cobra.Model):
        model = model_input
    elif isinstance(model_input, str):
        print("\nLoading the model...")
        model = load_model(model_input)
    else:
        raise ValueError(
            "The model_input must be a path to a model or a COBRApy model object."
        )

    # Add diet 1ba if desired
    if add_1ba:
        for rxn in model.reactions:
            if "Diet_" in rxn.id and rxn.lower_bound != 0:
                if (
                    "dgchol" in rxn.id
                    or "gchola" in rxn.id
                    or "tchola" in rxn.id
                    or "tdchola" in rxn.id
                ):
                    model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 0.0)

    #########################################################
    # Part 1: maximize the flux through all UFEt reactions
    #########################################################

    print(f"\n[STARTED] Part 1: maximizing UFEt fluxes for {model.name}")

    # Fetch all UFEt reactions and store them in a list
    UFEt_rxn_list = []
    for rxn in model.reactions:
        if "UFEt_" in rxn.id:
            UFEt_rxn_list.append(rxn.id)

    # Maximize the flux through all UFEt reactions
    counter = 0
    counter_max = len(UFEt_rxn_list)
    maximized_UFEt_flux_list = []
    for rxn in UFEt_rxn_list:
        counter += 1
        print(
            f"\nMaximizing UFEt reaction {str(counter)} of {str(counter_max)} \
              for {model.name}"
        )
        model.objective = rxn
        solution = model.optimize()
        maximized_UFEt_flux_list.append(solution.objective_value)
        print(f"{rxn}:\t{solution.objective_value}")

    # Create a dictionary of the maximized UFEt fluxes
    maximized_UFEt_flux_dict = dict(
        zip(UFEt_rxn_list, maximized_UFEt_flux_list, strict=False)
    )

    print(f"\n[COMPLETED] Part 1: maximization complete for {model.name}")

    #########################################################
    # Part 2: minimize the flux through all IEX reactions
    #########################################################

    print(f"\n[STARTED] Part 2: minimizing IEX fluxes for {model.name}")

    # Constrain the UFEt reactions by the maximized UFEt fluxes and minimize the
    # flux through all IEX reactions
    counter = 0
    counter_max = len(UFEt_rxn_list)
    minimized_IEX_flux_dict = dict()
    for i in range(len(UFEt_rxn_list)):
        counter += 1
        print(
            f"\nMinimizing IEX reaction {str(counter)} of {str(counter_max)} \
              for {model.name}"
        )
        if maximized_UFEt_flux_list[i] != 0.0:
            # Store the old bounds for the UFEt reaction
            saved_bounds = model.reactions.get_by_id(UFEt_rxn_list[i]).bounds

            # Set the bounds for the UFEt reaction to the calculated maximum
            model.reactions.get_by_id(UFEt_rxn_list[i]).bounds = (
                maximized_UFEt_flux_list[i],
                maximized_UFEt_flux_list[i],
            )

            # Rename the UFEt reaction to match the metabolite name
            metabolite = UFEt_rxn_list[i].replace("UFEt_", "") + "[u]"

            # Iterate over all IEX reactions for each metabolite
            for rxn in model.metabolites.get_by_id(metabolite).reactions:
                # If it is an IEX reaction, minimize the reaction flux
                if "IEX" in rxn.id:
                    model.objective = model.reactions.get_by_id(rxn.id)
                    solution = model.optimize(objective_sense="minimize")
                    minimized_IEX_flux_dict[rxn.id] = solution.objective_value
                    print(f"{rxn.id}:\t{solution.objective_value}")

            # Restore the bounds for the minimized IEX reaction
            model.reactions.get_by_id(UFEt_rxn_list[i]).bounds = saved_bounds

    # Create a dictionary of the minimized IEX fluxes
    model_rxn_bounds_dict = dict()
    for rxn in model.reactions:
        model_rxn_bounds_dict[rxn.id] = rxn.bounds
    print(f"\n[COMPLETED] Part 2: minimization complete for {model.name}")

    return maximized_UFEt_flux_dict, minimized_IEX_flux_dict, model_rxn_bounds_dict
