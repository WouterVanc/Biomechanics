# Replace marker names


# Fill in your specific marker names in the first argument, second argument = required Amsterdam Foot Model label
# Example 
# labels = [s.replace('Right_LatMal', 'ML') for s in labels]

def adjust_marker_names_R(labels): # Rightside
    labels = [s.replace('RTT', 'TT') for s in labels]
    labels = [s.replace('RFH', 'FH') for s in labels]
    labels = [s.replace('RLSHN', 'LSHN') for s in labels]
    labels = [s.replace('RASHN', 'ASHN') for s in labels]
    labels = [s.replace('RLM', 'LM') for s in labels]
    labels = [s.replace('RMM', 'MM') for s in labels]
    labels = [s.replace('RML', 'MM') for s in labels]
    labels = [s.replace('RCALP', 'CALP') for s in labels]
    labels = [s.replace('RCALD', 'CALD') for s in labels]
    labels = [s.replace('RPT', 'PT') for s in labels]
    labels = [s.replace('RST', 'ST') for s in labels]
    labels = [s.replace('RNAV', 'NAV') for s in labels]
    labels = [s.replace('RBM2', 'BM2') for s in labels]
    labels = [s.replace('RBM1', 'BM1') for s in labels]
    labels = [s.replace('RHM1', 'HM1') for s in labels]
    labels = [s.replace('RHM2', 'HM2') for s in labels]
    labels = [s.replace('RHM5', 'HM5') for s in labels]
    labels = [s.replace('RBM5', 'BM5') for s in labels]
    labels = [s.replace('Right HLX1', 'HLX') for s in labels]

    return labels

def adjust_marker_names_L(labels): # Leftside
    labels = [s.replace('LTT', 'TT') for s in labels]
    labels = [s.replace('LFH', 'FH') for s in labels]
    labels = [s.replace('LLSHN', 'LSHN') for s in labels]
    labels = [s.replace('LASHN', 'ASHN') for s in labels]
    labels = [s.replace('LLM', 'LM') for s in labels]
    labels = [s.replace('LMM', 'MM') for s in labels]
    labels = [s.replace('LML', 'MM') for s in labels]
    labels = [s.replace('LCALP', 'CALP') for s in labels]
    labels = [s.replace('LCALD', 'CALD') for s in labels]
    labels = [s.replace('LPT', 'PT') for s in labels]
    labels = [s.replace('LST', 'ST') for s in labels]
    labels = [s.replace('LNAV', 'NAV') for s in labels]
    labels = [s.replace('LBM2', 'BM2') for s in labels]
    labels = [s.replace('LBM1', 'BM1') for s in labels]
    labels = [s.replace('LHM1', 'HM1') for s in labels]
    labels = [s.replace('LHM2', 'HM2') for s in labels]
    labels = [s.replace('LHM5', 'HM5') for s in labels]
    labels = [s.replace('LBM5', 'BM5') for s in labels]
    labels = [s.replace('LHLX', 'HLX') for s in labels]

    return labels