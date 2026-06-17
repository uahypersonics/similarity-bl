# MATLAB Usage

MATLAB can use `simbl` through MATLAB's built-in Python interface. 

## Setup

First, make sure `simbl` is installed in a Python environment:

```bash
pip install similarity-bl
```

In MATLAB, point `pyenv` to that Python executable. For example, if `simbl` is
installed in the `dev` conda environment:

```matlab
pyenv("Version", "/opt/miniconda3/envs/dev/bin/python3")
```

Then import the package:

```matlab
simbl = py.importlib.import_module("simbl");
```

## Examples

=== "Falkner-Skan"

    This example solves an isothermal Mach 5 Falkner-Skan profile, maps the
    similarity coordinate to physical $y$, and plots $u/u_e$ and $T/T_e$.

    ```matlab
    % point MATLAB to the Python environment containing simbl
    pyenv("Version", "/opt/miniconda3/envs/dev/bin/python3")

    % load simbl
    simbl = py.importlib.import_module("simbl");

    % define the local similarity problem
    inputs = simbl.SimilarityInputs( ...
        pyargs( ...
            "mach_edge", 5.0, ...
            "temp_edge", 220.0, ...
            "temp_wall", 300.0, ...
            "wall_bc", "isothermal" ...
        ) ...
    );

    % solve the default Falkner-Skan equations
    result = simbl.solve_similarity(inputs);
    solution = result{1};
    info = result{2};

    % define dimensional local edge scales for this station
    x_station = 0.5;
    dens_edge = 0.02;
    uvel_edge = 500.0;
    visc_edge = 1.5e-5;

    % map eta to physical y
    y_py = simbl.eta2y( ...
        pyargs( ...
            "eta", solution.eta, ...
            "tau", solution.tau, ...
            "x", x_station, ...
            "dens_edge", dens_edge, ...
            "uvel_edge", uvel_edge, ...
            "visc_edge", visc_edge, ...
            "equations", "falkner_skan" ...
        ) ...
    );

    % convert Python/NumPy arrays to MATLAB double arrays
    eta = double(py.array.array('d', solution.eta.tolist()));
    y = double(py.array.array('d', y_py.tolist()));
    y_mm = 1000.0 .* y;
    uvel = double(py.array.array('d', solution.fp.tolist()));
    temp = double(py.array.array('d', solution.tau.tolist()));

    % plot profiles
    subplot(1, 2, 1)
    plot(uvel, y_mm, 'r', 'LineWidth', 2)
    grid on
    xlabel('u/U_e')
    ylabel('y (mm)')
    set(gca, 'FontSize', 16)

    subplot(1, 2, 2)
    plot(temp, y_mm, 'b', 'LineWidth', 2)
    grid on
    xlabel('T/T_e')
    ylabel('y (mm)')
    set(gca, 'FontSize', 16)
    ```

=== "Falkner-Skan-Cooke"

    For swept flows, pass `sweep_angle` in the inputs and select
    `falkner_skan_cooke` in the solver options. The FSC solution includes the
    crossflow profile `solution.g`, which is $w/w_e$. This example also maps
    the profiles to physical $y$ for comparison with CFD or experimental data.

    ```matlab
    % point MATLAB to the Python environment containing simbl
    pyenv("Version", "/opt/miniconda3/envs/dev/bin/python3")

    % load simbl
    simbl = py.importlib.import_module("simbl");

    % define the local swept-flow problem
    inputs = simbl.SimilarityInputs( ...
        pyargs( ...
            "mach_edge", 6.0, ...
            "temp_edge", 55.0, ...
            "temp_wall", 300.0, ...
            "wall_bc", "isothermal", ...
            "sweep_angle", 70.0 ...
        ) ...
    );

    % select the Falkner-Skan-Cooke equation family
    options = simbl.SolverOptions( ...
        pyargs( ...
            "equations", "falkner_skan_cooke" ...
        ) ...
    );

    % solve
    result = simbl.solve_similarity(inputs, options);
    solution = result{1};
    info = result{2};

    % define dimensional local edge scales for this station
    x_station = 0.5;
    dens_edge = 0.02;
    uvel_edge = 500.0;
    visc_edge = 1.5e-5;

    % map eta to physical y
    y_py = simbl.eta2y( ...
        pyargs( ...
            "eta", solution.eta, ...
            "tau", solution.tau, ...
            "x", x_station, ...
            "dens_edge", dens_edge, ...
            "uvel_edge", uvel_edge, ...
            "visc_edge", visc_edge, ...
            "equations", "falkner_skan_cooke" ...
        ) ...
    );

    % convert Python/NumPy arrays to MATLAB double arrays
    eta = double(py.array.array('d', solution.eta.tolist()));
    y = double(py.array.array('d', y_py.tolist()));
    y_mm = 1000.0 .* y;
    uvel = double(py.array.array('d', solution.fp.tolist()));
    wvel = double(py.array.array('d', solution.g.tolist()));
    temp = double(py.array.array('d', solution.tau.tolist()));
    dens = 1.0 ./ temp;

    % plot streamwise and crossflow profiles
    subplot(1, 2, 1)
    plot(uvel, y_mm, 'r', 'LineWidth', 2)
    grid on
    xlabel('u/U_e')
    ylabel('y (mm)')
    set(gca, 'FontSize', 16)

    subplot(1, 2, 2)
    plot(wvel, y_mm, 'b', 'LineWidth', 2)
    grid on
    xlabel('w/W_e')
    ylabel('y (mm)')
    set(gca, 'FontSize', 16)
    ```


## Field Names

The common solver output fields are summarized in the
[Run Solver](../workflow/run_solver.md) workflow page.

For density under the constant-pressure boundary-layer approximation, use
`rho_rhoe = 1.0 ./ temp` in MATLAB, where `temp` is $T/T_e$.
