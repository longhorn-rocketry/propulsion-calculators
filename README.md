# Propulsion Calculators

These are a series of calculators and simulators for hybrid rocket engines for the Longhorn Rocketry Association.

Most of this code was originally written by Akshay Kulkarni et al. for a senior design project in spring 2019. Josh Collins worked on most of the tank emptying portion of the tank combustion simulator.

## Code

- `regression-simulator` contains code to simulate the regression of a fuel grain cross-section over a burn so that said cross-section can be optimized for, e.g., constant thrust.
- `nozzle-calculator` computes optimal nozzle dimensions given certain input parameters.
- `hybrid-heat-transfer-calculator` is a 1D heat transfer simulation of the casing of a hybrid rocket combustion chamber to help determine how thin said casing can be under given combustion conditions.
- `bolt-calculator` contains a spreadsheet that can compute the safety factor for a given bolt loading.
- `tank-combustion-simulator` is an incomplete collection of ideas and code to completely simulate oxidizer flow and pressures throughout a hybrid engine burn.

## Documentation

Instructions on how to use most of this code are included in `documentation.docx`. Additional documentation, detailing some of the background theory (and the rest of the initial senior design project), is contained in `k_final_report_documentation.pdf`. 