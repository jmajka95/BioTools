# Exceptions file

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
