"""Valid governing equation systems for the similarity solver"""

# single source of truth for the supported equations identifiers
# used by both SolverConfig (schema.py) and SolverOptions (options.py) for validation
VALID_EQUATIONS: frozenset[str] = frozenset({"falkner_skan", "falkner_skan_cooke"})

# short aliases accepted anywhere equations is validated
# normalize these to their canonical form before checking VALID_EQUATIONS
EQUATIONS_ALIASES: dict[str, str] = {
    "fs": "falkner_skan",
    "fsc": "falkner_skan_cooke",
}
