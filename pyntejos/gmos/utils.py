"""Utilities for reduction of GMOS data"""

def ask_user(question, good_answers):
    """

    Parameters
    ----------
    question : string
    good_answers : list of strings

    Returns
    -------

    """
    answer = raw_input(question)
    if answer in good_answers:
        pass
    else:
        raise ValueError("Stopping process because user said so.")
