import pandas as pd


def build_profile_bridge(df8):
    return df8[["ProfileID", "ModelID", "Datatype"]].copy()
