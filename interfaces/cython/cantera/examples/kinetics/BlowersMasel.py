#!/usr/bin/env python
# coding: utf-8

# 
# ```
# conda install -c cantera/label/dev cantera
# ```

# In[15]:


import numpy as np

#%matplotlib inline
#from matplotlib import pyplot as plt

import cantera as ct


# In[29]:


input_file = """
description: |-
  Single OH reaction extracted from GRI-Mech 3.0.

units: {length: cm, time: s, quantity: mol, activation-energy: cal/mol}

phases:
- name: ohmech
  thermo: ideal-gas
  elements: [O, H]
  species: [H2, H, O, OH]
  kinetics: gas
  transport: mixture-averaged
  state: {T: 300.0, P: 1 atm}

species:
- name: H2
  composition: {H: 2}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [2.34433112, 7.98052075e-03, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12,
      -917.935173, 0.683010238]
    - [3.3372792, -4.94024731e-05, 4.99456778e-07, -1.79566394e-10, 2.00255376e-14,
      -950.158922, -3.20502331]
    note: TPIS78
  transport:
    model: gas
    geometry: linear
    well-depth: 38.0
    diameter: 2.92
    polarizability: 0.79
    rotational-relaxation: 280.0
- name: H
  composition: {H: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [2.5, 7.05332819e-13, -1.99591964e-15, 2.30081632e-18, -9.27732332e-22,
      2.54736599e+04, -0.446682853]
    - [2.50000001, -2.30842973e-11, 1.61561948e-14, -4.73515235e-18, 4.98197357e-22,
      2.54736599e+04, -0.446682914]
    note: L7/88
  transport:
    model: gas
    geometry: atom
    well-depth: 145.0
    diameter: 2.05
- name: O
  composition: {O: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [3.1682671, -3.27931884e-03, 6.64306396e-06, -6.12806624e-09, 2.11265971e-12,
      2.91222592e+04, 2.05193346]
    - [2.56942078, -8.59741137e-05, 4.19484589e-08, -1.00177799e-11, 1.22833691e-15,
      2.92175791e+04, 4.78433864]
    note: L1/90
  transport:
    model: gas
    geometry: atom
    well-depth: 80.0
    diameter: 2.75
- name: OH
  composition: {O: 1, H: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [3.99201543, -2.40131752e-03, 4.61793841e-06, -3.88113333e-09, 1.3641147e-12,
      3615.08056, -0.103925458]
    - [3.09288767, 5.48429716e-04, 1.26505228e-07, -8.79461556e-11, 1.17412376e-14,
      3858.657, 4.4766961]
    note: RUS78
  transport:
    model: gas
    geometry: linear
    well-depth: 80.0
    diameter: 2.75

reactions:
- equation: O + H2 <=> H + OH  # Reaction 1
  duplicate: true
  rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
- equation: O + H2 <=> H + OH  # Reaction 2
  duplicate: true
  type: blowers-masel
  rate-constant: {A: 3.87e+04, b: 2.7, Ea0: 6260.0, w0: 1e9}
  note: Example would be clearer if we tweak the Ea0 so the actual Ea matches Reaction 1
"""


# In[30]:


with open('input_file.yaml', 'w') as f:
    f.write(input_file)


# In[31]:


gas = ct.Solution('input_file.yaml')


# In[32]:


gas.species()


# In[37]:


gas.reactions()


# In[38]:


gas.species('OH')


# In[39]:


f"{gas.species('OH').thermo.h(298)/1e6 :.2f} kJ/mol"


# In[46]:


for i, r in enumerate(gas.reactions()):
    print(f"Reaction {i}")
    print(gas.reaction(i))
    print(f"∆Hrxn = {gas.delta_enthalpy[i]/1e6:.2f} kJ/mol")
    print(r.rate)
    print(f"At T = {gas.T} K")
    print(f"k = {gas.forward_rate_constants[i]:.2g} m3/kmol/s")
    print(f"k = {gas.forward_rate_constants[i]*1e6/1e3:.2g} cm3/mol/s")
    print("")


# In[24]:


def change_species_enthalpy(species_name, dH):
    """
    Find the species by name and change it's enthlapy by dH (in J/kmol)
    """
    index = gas.species_index(species_name)

    species = gas.species(index)
    print(f"Initial H(298) = {species.thermo.h(298)/1e6:.1f} kJ/mol")
    dx = dH / ct.gas_constant  # 'dx' is in fact (delta H / R). Note that R in cantera is 8314.462 J/kmol
    assert isinstance(species.thermo, ct.NasaPoly2)
    # print(species.thermo.coeffs)
    perturbed_coeffs = species.thermo.coeffs.copy()
    perturbed_coeffs[6] += dx
    perturbed_coeffs[13] += dx
    
    species.thermo = ct.NasaPoly2(species.thermo.min_temp, species.thermo.max_temp, 
                            species.thermo.reference_pressure, perturbed_coeffs)
    #print(species.thermo.coeffs)
    gas.modify_species(index, species)
    print(f"Modified H(298) = {species.thermo.h(298)/1e6:.1f} kJ/mol")


# In[25]:


change_species_enthalpy('OH', +10e6)


# In[49]:


f"{gas.species('OH').thermo.h(298)/1e6 :.2f} kJ/mol"


# In[55]:


# Change the T just to force Cantera to re-evaluate enthaplies etc.
gas.TP = 310, gas.P
gas.forward_rate_constants
# Then change it back
gas.TP = 300, gas.P
gas.T


# In[47]:


for i, r in enumerate(gas.reactions()):
    print(f"Reaction {i}")
    print(gas.reaction(i))
    print(f"∆Hrxn = {gas.delta_enthalpy[i]/1e6:.2f} kJ/mol")
    print(r.rate)
    print(f"At T = {gas.T} K")
    print(f"k = {gas.forward_rate_constants[i]:.2g} m3/kmol/s")
    print(f"k = {gas.forward_rate_constants[i]*1e6/1e3:.2g} cm3/mol/s")
    print("")


# Notice that Reaction 1, the rate has not changed even though we changed ∆Hrxn.
# But in Reaction 2, increasing ∆Hrxn has increased the barrier height, and the rate has slowed (compared to the cell above)

# In[ ]:




