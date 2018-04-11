
import operator

def get_intersection_of_list(list_of_list_features):
    return list(set.intersection(*map(set, list_of_list_features)))

def get_max_of_dict(inp):
    return max(inp.items(), key=operator.itemgetter(1))[0]

def argsort(seq, rev=False):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__, reverse=rev)
