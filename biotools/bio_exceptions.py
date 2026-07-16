"""File defining exceptions to be used during sequence creation and validation, reactions,
and more."""

class InvalidSequence(Exception):
    """Raised for an invalid DNA sequence"""
    pass

class InvalidInstantiationException(Exception):
    """Exception for invalid BioMolecules"""
    pass

class InvalidAnnotationException(Exception):
    """Exception for invalid Annotations."""
    pass

class ReactionError(Exception):
    """Class representing reaction errors."""
    pass

class InvalidPaddingException(Exception):
    """Class representing an exception for a lack of padding for enzymes."""
    pass

class SimulationError(Exception):
    """Class representing an error occurring during simulation using
    a BioReactionGraph."""
    pass
