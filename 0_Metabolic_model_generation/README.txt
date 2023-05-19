Here, the code to contextualize the model Human-GEM 1.14.0 to Ham´s media for the biomass production.

1) Download the Human-GEM 1.14.0 folder from the github repository: https://github.com/SysBioChalmers/Human-GEM/archive/refs/tags/v1.14.0.zip

2) The script 'get_all_models_essential_tasks' contextualizes the model to the essential task (in this case the biomass production).   
The flux through the input exchange reactions of metabolites not involved in the Ham’s media, as provided by Human1 for the biomass production, was set to zero. 

3) Then, the model is simplified by the function simplifyModel of Raven.
