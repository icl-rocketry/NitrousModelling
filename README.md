# NitrousModelling
Code relating to modelling nitrous oxide, including but not limited to:
- Properties of Nitrous Oxide in different conditions
- Predictions of flow rates for Nitrous Oxide
- Predictions of tank conditions for Nitrous Oxide tanks
A little bit of a mess of scripts, classes and data but useful.

Prefer to use the model built in simscape/simulink using the simscape two phase flow elements (https://github.com/icl-rocketry/PropulsionSystemModel), but this has some useful code snippets

# How to get it working (Requires CoolProp!!)
- Install Python 3.7 https://www.python.org/downloads/release/python-377/
- Install pip for your Python https://pip.pypa.io/en/stable/ (Make sure you install it for the version of python you installed)
- Tell MATLAB where your python is by running "pyversion C:\Users\<Name>\AppData\Local\Programs\Python\Python37\python.exe;"
- In MATLAB run the following command to install CoolProp: "[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']);"
- Verify that it is working and CoolProp can by found by running the command "NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE, FluidProperty.PRESSURE, 40e5, FluidProperty.VAPOR_QUALITY, 0)" and checking that an output is successfully returned.

# How to use it to get Nitrous properties at various conditions
- Make sure NitrousModelling is on the current MATLAB path (or is the folder open as the workspace)
- Use the command: "NitrousFluidCoolProp.getProperty(FluidProperty.<PICK_PROP_TO_GET>, FluidProperty.<PICK_PROP_TO_SPECIFY>, <it's value>, FluidProperty.<PICK_ANOTHER_PROP_TO_SPECIFY>, <it's value>)" to get a thermophysical property given two others. Eg. pressure and vapour quality (This will get a point on the saturation line). This will automatically determine the correct phase when it's possible to do so (or otherwise will assume one of them)
- To get a property for a specific phase you can use "NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.<PICK_PHASE>, FluidProperty.<PICK_PROP_TO_GET>, FluidProperty.<PICK_PROP_TO_SPECIFY>, <it's value>, FluidProperty.<PICK_ANOTHER_PROP_TO_SPECIFY>, <it's value>)".
- You can use tab completion to help pick from the allowed phases, and fluid properties. (Press tab after typing "FluidProperty.").
- Units are the SI units and match those specified in the cool prop documentation (http://www.coolprop.org/coolprop/HighLevelAPI.html#table-of-string-inputs-to-propssi-function). In the case of enthalpy and similar the 'mass' versions are used (not the 'molar' versions).

# Predict injector flow rates via NHNE
- You can use and modify injectorFlowRatesPredictor.m. Be aware that it is using a pretty poor time advancement scheme (A version using ode45 or similar would be better) and so the accuracy of results will depend on having a pretty small timestep.
- The tank model it uses is reasonably simple (does not model phenomenom in the same detail as the two phase receiver accumulator simscape block for example) and whilst it should provide an okay approximation of what happens in the tank it assumes always that the tank contents is entirely nitrous at the same temperature and in equilibrium.
- You can alternative directly 'solve' the value NHNE predicts for given upstream and downstream conditions using:
[X2,T2,v2,h2,G] = SaturatedNitrous.getDownstreamSaturatedNHNEFlowCond(X1,T1,P1,P2,v1,length)
where 1 refers to upstream conditions and 2 refers to downstream conditions and X is vapour quality, T is temperature in K, P is pressure in Pascal, v is velocity in m/s, G is mass flow per unit area in kg/m^2, h is specific enthalpy in J/Kg and length is the characteristic length for NHNE in metres.

# Predict Ball Valve via flow coefficient
- Look at ballValveTest.m

# Predict Linear Valve via flow coefficient
- Look at ballValveTest.m but replace BallValve with LinearValve

# Calculate orifice area exposed to flow for a full bore ball valve at a given opening angle
- Use FullBoreBallValve.getOrificeArea(theta,boreDiameter,ballDiameter) where boreDiameter and ballDiameter are in metres and theta is the angle that the valve is at (with pi/2 corresponding to being fully open). Area returned is in m^2
