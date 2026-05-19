"""Valid governing equation systems for the similarity solver"""

# single source of truth for the supported equations identifiers
# used by both SolverConfig (schema.py) and SolverOptions (options.py) for validation
VALID_EQUATIONS: frozenset[str] = frozenset({"falkner_skan", "falkner_skan_cooke"})
